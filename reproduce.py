#!/usr/bin/env python3
"""Offline replay of all reported computational stages into a new output directory."""
from pathlib import Path
import argparse,shutil,subprocess,sys,json,hashlib,platform,time,csv,math,os
P=argparse.ArgumentParser();P.add_argument('--output',type=Path,required=True);p=P.parse_args();ROOT=Path(__file__).resolve().parent;O=p.output.resolve()
if O.exists():raise SystemExit('Output must not exist: '+str(O))
O.mkdir(parents=True);W=O/'analysis';shutil.copytree(ROOT/'analysis',W)
N=W/'inputs/phase4_snapshot';C=N/'inputs/prior_snapshot';steps=[];comparisons=[]
def sha(x):return hashlib.sha256(x.read_bytes()).hexdigest()
# Preserve only explicit curated selection tables as inputs in the results tree.
for d in [W/'results',N/'results',C/'results',C/'phase2/results',C/'phase3/results']:
 for f in d.glob('*'):
  if f.is_dir():shutil.rmtree(f)
  elif f.name!='sequence_input_selection.tsv':f.unlink()
for d in [C/'audit',O/'logs']:d.mkdir(parents=True,exist_ok=True)
commands=[('catalogue',[C/'src/recompute_scores_and_donors.py']),('scoring_benchmarks',[C/'phase2/src/analyze.py','--out',C/'phase2/results']),('structure',[C/'phase2/src/structures.py','--out',C/'phase2/results']),('additional_sensitivities',[C/'phase2/src/supplement.py','--out',C/'phase2/results']),('sequence_organ_antigen',[C/'phase3/src/complete_analyses.py','--output',C/'phase3/results']),('noncircular_score_diagnostic',[C/'phase3/src/isi_leakage_audit.py','--output',C/'phase3/results']),('glycan_array',[C/'phase3/src/extract_nanno2021_glycan_array.py']),('public_proteomics_function',[N/'src/analyze_public_data.py','--input',N/'inputs','--output',N/'results']),('identifier_sensitivity',[W/'src/audit_mapping.py','--input',W/'inputs','--output',W/'results']),('figure_sources',[ROOT/'code/extract_figure_sources.py',O]),('supplementary_displays',[ROOT/'code/extract_review_tables.py',O])]
for label,args in commands:
 print('RUN',label,flush=True);start=time.time();cmd=[sys.executable]+[str(x) for x in args]
 with (O/'logs'/f'{label}.log').open('w') as log:r=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,cwd=O)
 steps.append(dict(stage=label,command=cmd,seconds=time.time()-start,returncode=r.returncode))
 if r.returncode:
  print((O/'logs'/f'{label}.log').read_text()[-5000:]);break
extra=[p for p in W.rglob('*.json') if p.parent.name=='results' and p.name in ['mapping_summary.json','liu_flow.json','liu_top_decile_enrichment.json','isi_leakage_diagnosis.json','nanno2021_glycan_array_qc.json']]+list(W.rglob('results/models/*.pdb'))+[O/'figure_data/checked_review_claims.json']
for f in sorted(list(W.rglob('results/*.tsv'))+list((O/'figures').glob('*.tsv'))+list((O/'figure_data').glob('*.tsv'))+extra):
 if f.name=='sequence_input_selection.tsv':continue
 expected=ROOT/f.relative_to(O)
 if not expected.exists():comparisons.append(dict(file=str(f.relative_to(O)),status='NEW_OUTPUT_NOT_IN_FROZEN_RESULTS'));continue
 exact=sha(f)==sha(expected);status='BYTE_IDENTICAL' if exact else 'DIFFERENT'
 # Numeric serialization may vary slightly across supported scientific-library builds.
 if not exact and f.suffix in ['.json','.pdb']:
  comparisons.append(dict(file=str(f.relative_to(O)),status='DIFFERENT'));continue
 if not exact:
  with f.open(newline='') as q:a=list(csv.reader(q,delimiter='\t'))
  with expected.open(newline='') as q:b=list(csv.reader(q,delimiter='\t'))
  same=len(a)==len(b) and all(len(x)==len(y) for x,y in zip(a,b));maxdiff=0.;diffs=[]
  if same:
   for i,(x,y) in enumerate(zip(a,b)):
    for j,(u,v) in enumerate(zip(x,y)):
     if u==v:continue
     try:
      uu,vv=float(u),float(v);ok=(math.isnan(uu) and math.isnan(vv)) or math.isclose(uu,vv,rel_tol=1e-10,abs_tol=1e-10)
      if not ok:diffs.append([i,j,u,v])
      if math.isfinite(uu) and math.isfinite(vv):maxdiff=max(maxdiff,abs(uu-vv))
     except ValueError:diffs.append([i,j,u[:100],v[:100]])
   same=not diffs
  if same:status='NUMERICALLY_IDENTICAL_1e-10'
  comparisons.append(dict(file=str(f.relative_to(O)),status=status,max_absolute_numeric_difference=maxdiff,first_differences=diffs[:5]))
 else:comparisons.append(dict(file=str(f.relative_to(O)),status=status))
report=dict(python=sys.version,platform=platform.platform(),steps=steps,comparisons=comparisons,passed=all(x['returncode']==0 for x in steps) and len(steps)==len(commands) and all(x['status'] in ['BYTE_IDENTICAL','NUMERICALLY_IDENTICAL_1e-10'] for x in comparisons),note='Frozen curated sequence selection is an input. Numerical result directories were emptied before computation. Fresh runtime/provenance JSON is intentionally not byte-compared to historical timestamps.')
(O/'REPLAY_REPORT.json').write_text(json.dumps(report,indent=2)+'\n');print('RESULT',report['passed'],'tables',len(comparisons),flush=True)
raise SystemExit(0 if report['passed'] else 1)
