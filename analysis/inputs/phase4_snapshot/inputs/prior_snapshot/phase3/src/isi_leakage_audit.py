"""Fixed feature audit, not independent functional validation. No label-derived predictors."""
from pathlib import Path
import argparse,hashlib,json,re,sys,platform,datetime
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score,average_precision_score
BASE=Path(__file__).resolve().parents[1]
SEED=20261103

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def source_aliases():
 aliases={}
 for p in sorted((BASE/'inputs/public').glob('ref_*_europepmc.json')):
  for r in json.loads(p.read_text()).get('resultList',{}).get('result',[]):
   key='doi:'+r['doi'].lower() if r.get('doi') else 'pmid:'+r['pmid'] if r.get('pmid') else None
   if key:
    if r.get('doi'):aliases['doi:'+r['doi'].lower()]=key
    if r.get('pmid'):aliases['pmid:'+r['pmid']]=key
    if r.get('pmcid'):aliases['pmc:'+r['pmcid'].upper()]=key
 return aliases

def source_tokens(s,aliases):
 tokens=[]
 for item in str(s).split(';'):
  x=item.strip().lower().removeprefix('https://doi.org/').removeprefix('doi:').strip()
  if re.fullmatch(r'10\.\d{4,9}/\S+',x):k='doi:'+x
  elif re.fullmatch(r'pmc\d+',x):k='pmc:'+x.upper()
  elif re.fullmatch(r'\d+',x):k='pmid:'+x
  else:k='unresolved:'+x
  tokens.append(aliases.get(k,k))
 return sorted(set(tokens))

def groups(t):
 aliases=source_aliases();parent=list(range(len(t)));seen={};rows=[]
 def find(i):
  while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
  return i
 for i,r in enumerate(t.to_dict('records')):
  genes=sorted(set(re.findall(r'[A-Z][A-Z0-9]+',str(r['gene_key'])+' '+str(r['human_partner']))))
  st=source_tokens(r['source'],aliases);tokens=['g:'+g for g in genes]+['s:'+s for s in st]
  for token in tokens:
   if token in seen:parent[find(i)]=find(seen[token])
   else:seen[token]=i
  rows.append(dict(pair_id=r['pair_id'],source=r['source'],source_tokens=';'.join(st),gene_tokens=';'.join(genes),source_parse_unresolved=any(x.startswith('unresolved:') for x in st)))
 g=[find(i) for i in range(len(t))];canon={v:min(t.iloc[[j for j,x in enumerate(g) if x==v]].pair_id) for v in set(g)}
 for r,x in zip(rows,g):r['dependency_group']=canon[x]
 return pd.DataFrame(rows)

def evaluate(t,B=2000):
 dep=groups(t);y=t.evidence.ne('control').to_numpy(int);s=100-t.identity.to_numpy(float);g=dep.dependency_group.to_numpy();keys=np.unique(g);rng=np.random.default_rng(SEED);stats=[]
 for _ in range(B):
  ix=np.concatenate([np.flatnonzero(g==k) for k in rng.choice(keys,len(keys),replace=True)])
  if len(np.unique(y[ix]))==2:stats.append([roc_auc_score(y[ix],s[ix]),average_precision_score(y[ix],s[ix])])
 a=np.array(stats);return dict(n=len(t),n_noncontrol=int(y.sum()),n_control=int((1-y).sum()),n_dependency_groups=len(keys),auroc=roc_auc_score(y,s),auroc_ci_low=np.quantile(a[:,0],.025),auroc_ci_high=np.quantile(a[:,0],.975),average_precision=average_precision_score(y,s),ap_ci_low=np.quantile(a[:,1],.025),ap_ci_high=np.quantile(a[:,1],.975),prevalence=y.mean(),bootstrap_valid=len(a),bootstrap_requested=B,p_value=np.nan,method='source_normalized_gene_source_dependency_bootstrap',interpretation='fixed_catalogue_feature_association_NOT_functional_validation'),dep

def run(dest):
 dest.mkdir(parents=True,exist_ok=True);p=BASE/'inputs/tier1_frozen.tsv';t=pd.read_csv(p,sep='\t');assert not t.pair_id.duplicated().any()
 for c in ['identity','func_sev','ISI_v2']:t[c]=pd.to_numeric(t[c],errors='coerce')
 rows=[];common=t[t.identity.notna()].copy().reset_index(drop=True)
 for name,z in [('sequence_available',common),('exclude_disputed',common[common.evidence.ne('disputed')].reset_index(drop=True))]:
  m,dep=evaluate(z);rows.append(dict(scope=name,**m));dep['scope']=name
  if name=='sequence_available':dep.to_csv(dest/'isi_dependency_groups.tsv',sep='\t',index=False)
 pd.DataFrame(rows).to_csv(dest/'isi_nonlabel_feature_metrics.tsv',sep='\t',index=False,na_rep='NA',float_format='%.15g')
 t[['pair_id','evidence','identity','func_sev','ISI_v2','incompat_type','source']].assign(sequence_available=t.identity.notna(),control_zero_by_definition=t.evidence.eq('control')&t.func_sev.eq(0)).to_csv(dest/'isi_leakage_axis_audit.tsv',sep='\t',index=False,na_rep='NA')
 y=common.evidence.ne('control').to_numpy(int)
 diagnosis={'all_axes':len(t),'controls':int(t.evidence.eq('control').sum()),'controls_zero_severity':int((t.evidence.eq('control')&t.func_sev.eq(0)).sum()),'noncontrols_positive_severity':int((t.evidence.ne('control')&t.func_sev.gt(0)).sum()),'common_n':len(common),'severity_only_auc':roc_auc_score(y,common.func_sev),'frozen_ISI_auc':roc_auc_score(y,common.ISI_v2),'missing_sequence_controls':t.loc[t.evidence.eq('control')&t.identity.isna(),'pair_id'].tolist(),'independent_validation':False,'severity_auc_permitted_use':'definition_leakage_diagnostic_only','p_value':None}
 (dest/'isi_leakage_diagnosis.json').write_text(json.dumps(diagnosis,indent=2)+'\n')
 return diagnosis
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True,type=Path);a=p.parse_args();d=run(a.output)
 print(json.dumps({'status':'complete','n':d['common_n'],'independent_validation':False}))
