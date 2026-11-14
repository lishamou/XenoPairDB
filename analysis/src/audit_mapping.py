"""Post hoc identifier sensitivity; frozen primary, labels and scores are not edited."""
from pathlib import Path
import argparse,json,hashlib,sys,ast,re
import numpy as np
import pandas as pd
from scipy import stats
P=argparse.ArgumentParser();P.add_argument('--input',type=Path,required=True);P.add_argument('--output',type=Path,required=True);P.add_argument('--mapping-only',action='store_true');a=P.parse_args();I=a.input;O=a.output;O.mkdir(parents=True,exist_ok=True)
def save(d,n):pd.DataFrame(d).to_csv(O/(n+'.tsv'),sep='\t',index=False,float_format='%.12g',na_rep='NA')
def js(d,n):(O/(n+'.json')).write_text(json.dumps(d,indent=2,ensure_ascii=False)+'\n')
f=pd.read_csv(I/'frozen_antigen.tsv',sep='\t');d=pd.read_csv(I/'primary_mapping.tsv',sep='\t');F=f.set_index('acc');gene2acc=dict(zip(f.gene.dropna(),f.loc[f.gene.notna(),'acc']))
entries={};provenance={}
for file in sorted((I/'uniprot').glob('current_*.json')):
 if file.name.endswith('.meta.json'):continue
 for e in json.loads(file.read_text())['results']:
  acc=e['primaryAccession'];entries[acc]=e;provenance[acc]=str(file.relative_to(I))
for file in sorted((I/'uniprot').glob('supplement_*.json')):
 if file.name.endswith('.meta.json'):continue
 obj=json.loads(file.read_text());rows=obj.get('results',[obj])
 for e in rows:
  if 'primaryAccession' not in e:continue
  acc=e['primaryAccession'];entries[acc]=e;provenance[acc]=str(file.relative_to(I))
byacc={}
for acc,e in entries.items():
 if e.get('organism',{}).get('taxonId')!=9823:continue
 for k in [acc]+e.get('secondaryAccessions',[]):byacc.setdefault(k,set()).add(acc)
# Direct API redirects are explicit retired-accession relationships even when
# the destination record does not repeat the requested accession as secondary.
for file in sorted((I/'uniprot').glob('supplement_*.meta.json')):
 meta=json.loads(file.read_text());payload=file.with_name(file.name.replace('.meta.json','.json'))
 obj=json.loads(payload.read_text())
 current=obj.get('primaryAccession')
 if current and obj.get('organism',{}).get('taxonId')==9823:
  byacc.setdefault(meta['requested'],set()).add(current)
def record(acc):
 hits=byacc.get(acc,set())
 if len(hits)!=1:return None
 e=entries[next(iter(hits))]
 return e if e.get('organism',{}).get('taxonId')==9823 else None
def annot(acc):
 e=record(acc)
 if not e:return dict(accession=acc,current='',geneids=set(),names=set(),primary_names=set(),sequence='',resolved=False)
 gs=e.get('genes',[]);primary={g['geneName']['value'] for g in gs if 'geneName' in g};syn={v['value'] for g in gs for v in g.get('synonyms',[])}
 return dict(accession=acc,current=e['primaryAccession'],geneids={x['id'] for x in e.get('uniProtKBCrossReferences',[]) if x['database']=='GeneID'},names=primary|syn,primary_names=primary,
    sequence=e.get('sequence',{}).get('value',''),resolved=True)
src=sorted(set(';'.join(d.source_accession).split(';')));A={acc:annot(acc) for acc in sorted(set(src)|set(f.acc))}
lookup=[]
for acc,v in A.items():
 e=record(acc) or {};lookup.append(dict(accession=acc,in_source=acc in src,in_frozen=acc in F.index,resolved=v['resolved'],current=v['current'],
  geneids=';'.join(sorted(v['geneids'])),primary_genes=';'.join(sorted(v['primary_names'])),all_gene_aliases=';'.join(sorted(v['names'])),
  taxon=e.get('organism',{}).get('taxonId'),entry_version=e.get('entryAudit',{}).get('entryVersion'),sequence_version=e.get('entryAudit',{}).get('sequenceVersion'),
  last_update=e.get('entryAudit',{}).get('lastAnnotationUpdateDate'),sequence_length=len(v['sequence']),source_file=provenance.get(v['current'],'')))
save(lookup,'identifier_inventory')
A4=dict(A);archive_records=[]
for file in sorted((I/'uniprot').glob('archive_entry_*.txt')):
 text=file.read_text();acc=file.stem.replace('archive_entry_','')
 if A[acc]['resolved']:continue
 assert re.search(r'NCBI_TaxID=9823\b',text),'Archive species mismatch'
 g=set(re.findall(r'^GN\s+Name=([^ ;{]+)',text,re.M));gid=set(re.findall(r'^DR\s+GeneID; ([^;]+);',text,re.M))
 seq=''.join(re.findall(r'[A-Z]+',text.split('\nSQ   ')[1].split('\n',1)[1].split('//')[0]))
 A4[acc]=dict(accession=acc,current='archive:'+acc,geneids=gid,names=g,primary_names=g,sequence=seq,resolved=True)
 archive_records.append(dict(accession=acc,gene=';'.join(g),source=str(file.relative_to(I)),entry_version=json.loads(file.with_suffix('.meta.json').read_text())['entry_version'],sequence_length=len(seq)))
save(archive_records,'archive_resolutions')
def match(acc,level,annotations=A):
 v=annotations[acc];hits=[]
 # An exact stored accession is an identity even when annotation retrieval fails.
 if acc in F.index:return acc,'stored_accession',''
 if not v['resolved']:return None,'unresolved_accession',''
 current=[k for k in F.index if A[k]['current'] and A[k]['current']==v['current']]
 if len(current)==1:return current[0],'current_or_secondary_accession',v['current']
 if len(current)>1:return None,'ambiguous_frozen_current', ';'.join(current)
 if level==1:return None,'no_current_accession_match',''
 if len(v['geneids'])==1:
  same=[k for k in F.index if len(A[k]['geneids'])==1 and A[k]['geneids']==v['geneids']]
  if len(same)==1:return same[0],'shared_unique_GeneID',';'.join(v['geneids'])
  if len(same)>1:return None,'ambiguous_frozen_GeneID',';'.join(same)
 aliases=[]
 for k in F.index:
  t=A[k]
  if not t['resolved'] or not v['names']&t['names']:continue
  if v['geneids'] and t['geneids'] and v['geneids'].isdisjoint(t['geneids']):continue
  aliases.append(k)
 if len(aliases)==1:return aliases[0],'official_gene_alias',';'.join(sorted(v['names']&A[aliases[0]]['names']))
 return None,'ambiguous_alias' if aliases else 'no_frozen_identity',';'.join(aliases)
matches={acc:{level:match(acc,level) for level in [1,2]} for acc in src}
matches4={acc:{2:match(acc,2,A4)} for acc in src}
def group(members,annotations,matches):
 ms=[annotations[m] for m in members];targets=[matches[m][2][0] for m in members]
 if not all(v['resolved'] for v in ms):return None,'group_unresolved_member',''
 gidsets=[v['geneids'] for v in ms];known=[x for x in gidsets if x]
 stable=all(len(x)==1 for x in gidsets) and len(set(next(iter(x)) for x in gidsets))==1
 primary_sets=[v['primary_names'] for v in ms]
 same_name=all(len(x)==1 for x in primary_sets) and len(set(next(iter(x)) for x in primary_sets))==1
 gid_conflict=len(set.union(*known))>1 if known else False
 same_entry=len({v['current'] for v in ms})==1
 if not (stable or same_entry or (same_name and not gid_conflict)):return None,'group_distinct_or_unresolved_gene_identity',''
 if any(t is None for t in targets):return None,'same_gene_group_no_unique_frozen_target',''
 if len(set(targets))!=1:return None,'group_multiple_frozen_targets',''
 return targets[0],('same_gene_group_GeneID' if stable else 'same_gene_group_current' if same_entry else 'same_gene_group_official_name'),(';'.join(sorted(set.union(*(v['geneids'] for v in ms)))) if stable else ';'.join(sorted(ms[0]['primary_names'])))
member_rows=[]
for acc in src:
 for level in [1,2]:
  target,method,evidence=matches[acc][level];member_rows.append(dict(source_accession=acc,level=level,target=target,method=method,evidence=evidence))
save(member_rows,'member_mapping')
allrows=[];baseline_conflicts=[]
for r in d.itertuples():
 members=r.source_accession.split(';');base=gene2acc.get(r.gene) if isinstance(r.gene,str) else None
 v1=base;v2=base;v3=base;reason='baseline' if base else '';evidence=''
 if base and len(members)==1 and A[members[0]]['geneids'] and A[base]['geneids'] and A[members[0]]['geneids'].isdisjoint(A[base]['geneids']):
  baseline_conflicts.append(dict(source_row=r.source_row,source_accession=members[0],baseline=base,reason='disjoint_GeneID'))
 if not base and len(members)==1:
  v1,m1,e1=matches[members[0]][1];v2,m2,e2=matches[members[0]][2];v3=v2;reason=m2;evidence=e2
 elif not base:v3,reason,evidence=group(members,A,matches)
 v4=v3
 if not v4:
  if len(members)>1:v4,r4,e4=group(members,A4,matches4)
  else:v4,r4,e4=matches4[members[0]][2]
  reason=('archive_supported_' if v4 else 'archive_checked_' if any(m in {x['accession'] for x in archive_records} for m in members) else '')+r4;evidence=e4
 row=dict(source_row=r.source_row,source_accession=r.source_accession,source_gene=r.source_gene,peptide=r.peptide,observed_PK=r.observed_PK,observed_E=r.observed_E,log2fc=r.log2fc,
  M0=base,M1=v1,M2=v2,M3=v3,M4=v4,M5=v4 if v4 and r.peptide in A[v4]['sequence'] and A[v4]['sequence'] else None,
  mapping_reason=reason,mapping_evidence=evidence,new_target=v4 if v4!=base else None,
  target_gene=F.loc[v4,'gene'] if v4 else None,target_sequence_contains_peptide=bool(v4 and A[v4]['sequence'] and r.peptide in A[v4]['sequence']))
 allrows.append(row)
md=pd.DataFrame(allrows);save(md,'row_mapping_audit');save(baseline_conflicts,'baseline_identity_conflicts')
summary=dict(source_unique_accessions=len(src),frozen_accessions=len(f),union_accessions=len(A),retrieved_pig_source=sum(A[k]['resolved'] for k in src),
 retrieved_pig_frozen=sum(A[k]['resolved'] for k in F.index),missing_source=[k for k in src if not A[k]['resolved']],missing_frozen=[k for k in F.index if not A[k]['resolved']],
 baseline_conflicts=baseline_conflicts,source_suffix_accessions=[k for k in src if '-' in k or '.' in k],
 source_observed_rows=int((d.observed_E>0).sum()),primary_mapping_observed_counts=d[d.observed_E>0].mapping.value_counts().to_dict(),
 new_observed_rows=md[(md.M4.notna())&(md.M0.isna())&(md.observed_E>0)][['source_row','source_accession','source_gene','new_target','mapping_reason']].to_dict('records'),
 M4_archive_status='Two deleted source entries resolved from stored UniSave flatfiles; archived gene annotations used only in M4/M5',archive_records=archive_records)
js(summary,'mapping_summary')
if a.mapping_only:print('Mapping complete; newly mapped observed rows',len(summary['new_observed_rows']));sys.exit()

# Reuse the previously independently verified statistical implementation as a pure function.
tree=ast.parse((I/'phase4_analysis_code.py').read_text());module=ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['normrank','association']],type_ignores=[])
SEED=20261104;BOOT=10000;PERM=99999
exec(compile(module,'frozen_statistical_helpers','exec'),globals())
metrics=[];flows=[];baseline_cohorts={}
for scheme in ['M0','M1','M2','M3','M4','M5']:
 for endpoint in ['observed_control','source_model_FC']:
  q=md[md[scheme].notna()&(md.observed_PK>0)&np.isfinite(md.log2fc)].copy()
  if endpoint=='observed_control':q=q[q.observed_E>0]
  # Shared peptide exclusion is computed across all available mapped rows before endpoint filtering.
  allmapped=md[md[scheme].notna()];counts=allmapped.groupby('peptide')[scheme].nunique();shared=set(counts[counts>1].index)
  q=q[~q.peptide.isin(shared)]
  ag=q.groupby([scheme,'peptide'],as_index=False).log2fc.median().groupby(scheme,as_index=False).agg(log2fc=('log2fc','median'),peptides=('peptide','size')).rename(columns={scheme:'acc'})
  ag=ag.merge(f[['acc','gene','antigenicity_v2','rank_v2','historical_positive']],on='acc',validate='one_to_one');ag=ag[~ag.historical_positive].sort_values(['gene','acc'],na_position='last').reset_index(drop=True)
  if scheme=='M0':baseline_cohorts[endpoint]=set(ag.acc)
  save(ag,f'cohort_{scheme}_{endpoint}');v=association(ag,f'{scheme}_{endpoint}');v.update(scheme=scheme,endpoint=endpoint,new_genes=';'.join(ag.loc[~ag.acc.isin(pd.read_csv(I/'primary_cohort.tsv',sep='\t').acc),'acc']))
  metrics.append(v);flows.append(dict(scheme=scheme,endpoint=endpoint,mapped_rows=int(md[scheme].notna().sum()),mapped_control_observed_rows=int((md[scheme].notna()&(md.observed_E>0)).sum()),
   shared_peptide_sequences=len(shared),included_peptides=int(ag.peptides.sum()),n=len(ag)))
metrics=pd.DataFrame(metrics);ix=(metrics.scheme!='M0')&metrics.p.notna();metrics.loc[ix,'q_BH_mapping_family']=stats.false_discovery_control(metrics.loc[ix,'p'])
save(metrics,'mapping_sensitivity_statistics');save(flows,'mapping_funnel')
js(dict(command=sys.argv,seed=SEED,bootstrap=BOOT,permutations=PERM,BH_family=int(ix.sum()),python=sys.version,
 input_hashes={str(p.relative_to(I)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(I.rglob('*')) if p.is_file()},
 code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),'run_metadata')
print('Completed',len(metrics),'mapping/outcome scenarios')
