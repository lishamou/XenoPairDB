#!/usr/bin/env python3
"""Read-only OOXML extraction of public Nanno 2021 S1; no hypothesis tests.
Run from project root: python3 revision_20261103/phase3/src/extract_nanno2021_glycan_array.py
Network used only when a named source is absent; existing source bytes never overwritten.
The S1 timepoint is deliberately unresolved because caption and body conflict.
"""
import csv, hashlib, json, math, platform, sys, urllib.request, zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import xml.etree.ElementTree as ET

BASE=Path(__file__).resolve().parents[1]
ROOT=BASE.parent
PUBLIC=BASE/'inputs/public'
RESULTS=BASE/'results'
URLS={
 'nanno2021_article.html':'https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0253029',
 'nanno2021_s001.docx':'https://journals.plos.org/plosone/article/file?type=supplementary&id=10.1371/journal.pone.0253029.s001',
}
W='{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'
MONKEYS=['13GP08','12JP01','13CP03','13CP10','13GP10']
POST_DAYS=[142,46,42,28,35]

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def write_tsv(name,rows):
 p=RESULTS/name
 with p.open('w',newline='') as f:
  wr=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');wr.writeheader();wr.writerows(rows)
 return p

def main():
 PUBLIC.mkdir(parents=True,exist_ok=True);RESULTS.mkdir(parents=True,exist_ok=True)
 sources=[]
 for name,url in URLS.items():
  p=PUBLIC/name
  if not p.exists():
   with urllib.request.urlopen(url,timeout=45) as r: body=r.read()
   with p.open('xb') as f:f.write(body)
  sources.append({'path':str(p.relative_to(ROOT)),'url':url,'bytes':p.stat().st_size,'sha256':sha(p),'accessed_utc':datetime.now(timezone.utc).isoformat(),'reuse_existing_bytes':True})
 for extra in [PUBLIC/'PMC12805290.xml',BASE/'inputs/tier1_frozen.tsv',BASE.parent/'phase2/inputs/legacy_antigen_evidence.md',BASE.parent/'phase2/results/antigen_external_evidence.tsv']:
  if extra.exists(): sources.append({'path':str(extra.relative_to(ROOT)),'bytes':extra.stat().st_size,'sha256':sha(extra),'use':'reference trace or exact-string prior-use audit; not numeric source'})
 doc=PUBLIC/'nanno2021_s001.docx'
 with zipfile.ZipFile(doc) as z:x=ET.fromstring(z.read('word/document.xml'))
 tabs=x.findall('.//'+W+'tbl');assert len(tabs)==1, 'Unexpected table count'
 rows=[[''.join(t.text or '' for t in c.iter(W+'t')) for c in r.findall(W+'tc')] for r in tabs[0].findall(W+'tr')]
 assert rows[0]==['a411 print #','','','IgM signal','','IgG signal'],rows[0]
 assert rows[1][3:8]==MONKEYS and rows[1][9:14]==MONKEYS
 assert len(rows)==410 and all(len(r)==14 for r in rows[1:])
 caption='\n'.join(''.join(t.text or '' for t in p.iter(W+'t')) for p in x.findall(W+'body/'+W+'p')).strip()
 animals=[dict(monkey_id=m,recipient_species='Macaca fascicularis',donor_species='Sus scrofa',graft='adult_porcine_islets',donor_genotype='GGTA1-KO' if m=='13GP08' else 'WT',study_post_sample_day=d,study_post_sample_event='clinical_rejection',S1_signal_timepoint='UNRESOLVED_caption_natural_body_posttransplant',baseline_sample_day='NOT_REPORTED',independent_biological_unit='recipient_monkey',metadata_source='article Methods: Pig islet isolation and transplantation') for m,d in zip(MONKEYS,POST_DAYS)]
 components=[];long=[]
 for rownum,r in enumerate(rows[2:],start=3):
  pid=int(r[0]);assert r[8]=='', (rownum,'nonempty separator')
  components.append(dict(array_version='A411',array_print_id=pid,abbreviation=r[1],description=r[2],explicit_control_by_label='yes' if r[1] in ['Cy3','BSA'] else 'not_classified',source_table_row_1based=rownum))
  for iso,start in [('IgM',3),('IgG',9)]:
   for j,m in enumerate(MONKEYS):
    raw=r[start+j]
    if raw: assert math.isfinite(float(raw))
    long.append(dict(array_version='A411',array_print_id=pid,abbreviation=r[1],description=r[2],monkey_id=m,isotype=iso,timepoint='UNRESOLVED_caption_natural_body_posttransplant',value_reported=raw,missing_in_source=str(raw==''),value_scale='log2_normalized_signal',signal_is_raw='False',technical_replicate='aggregated_not_available',source_docx='nanno2021_s001.docx',source_table_1based=1,source_row_1based=rownum,source_cell_1based=start+j+1))
 ids=[c['array_print_id'] for c in components]
 keys=[(r['array_print_id'],r['monkey_id'],r['isotype'],r['timepoint']) for r in long]
 assert len(components)==408 and len(long)==4080 and len(set(ids))==408 and len(set(keys))==4080
 assert Counter(r['monkey_id'] for r in long)==Counter({m:816 for m in MONKEYS})
 outputs=[write_tsv('nanno2021_glycan_components.tsv',components),write_tsv('nanno2021_monkey_metadata.tsv',animals),write_tsv('nanno2021_glycan_array_long.tsv',long)]
 qc={
  'status':'EXTRACTION_VERIFIED_TIMEPOINT_UNRESOLVED_PREPOST_NOT_COMPUTABLE',
  'S1_caption_verbatim':caption,'n_tables':len(tabs),'n_components':len(components),'n_monkeys':len(MONKEYS),'n_isotypes':2,'n_matrix_cells':len(long),'n_numeric_values':sum(bool(r['value_reported']) for r in long),'n_signal_matrices':1,
  'n_source_missing_numeric_values':sum(not r['value_reported'] for r in long),'missing_cells':[{'print_id':r['array_print_id'],'monkey_id':r['monkey_id'],'isotype':r['isotype']} for r in long if not r['value_reported']],'unique_array_print_ids':len(set(ids)),'array_print_id_range':[min(ids),max(ids)],'missing_print_ids_within_range':sorted(set(range(min(ids),max(ids)+1))-set(ids)),
  'duplicate_abbreviations':{k:v for k,v in Counter(c['abbreviation'] for c in components).items() if v>1},
  'explicit_control_labels':[c['abbreviation'] for c in components if c['explicit_control_by_label']=='yes'],
  'technical_details_from_article':{'serum_dilution':'1:50','wells_per_sample_per_isotype':2,'intensity_floor_before_log2':150,'reported_values':'duplicate-spot averaged, reference normalized, log2 transformed','technical_replicates_in_S1':False},
  'paired_prepost_recomputable_from_S1':False,'separately_labeled_baseline_matrix_present':False,'post_matrix_unambiguously_labeled':False,'peak_timepoint_present':False,
  'analysis_performed':'structural extraction and missingness/uniqueness checks only; no tests or effect estimates',
  'restrictions':['408 array components are not 408 biological replicates','array contains controls and glycopeptides; not 408 unique protein antigens','recipient n=5; genotype WT n=4, GGTA1KO n=1; rejection day differs by recipient','no per-protein glycosylation occupancy','cannot establish ISI/protein ranking validation','never infer baseline or timepoint from signal magnitude']}
 q=RESULTS/'nanno2021_glycan_array_qc.json';q.write_text(json.dumps(qc,ensure_ascii=False,indent=2));outputs.append(q)
 manifest={'time_utc':datetime.now(timezone.utc).isoformat(),'command':' '.join(sys.argv),'python':platform.python_version(),'script_sha256':sha(Path(__file__)),'seed':'not_applicable_no_random_operations','source_license':'PLOS article declares CC0; supplement publicly linked','inputs':sources,'outputs':[{'path':str(p.relative_to(ROOT)),'sha256':sha(p),'bytes':p.stat().st_size} for p in outputs]}
 (RESULTS/'nanno2021_glycan_array_provenance.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
 print(json.dumps({k:qc[k] for k in ['status','n_components','n_monkeys','n_numeric_values','missing_print_ids_within_range','duplicate_abbreviations']}))
if __name__=='__main__':main()
