"""Exploratory calculations only. Organ annotations are not calibrated organ weights."""
from pathlib import Path
import argparse,json,hashlib,re,sys,csv
import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner,substitution_matrices
from scipy.stats import mannwhitneyu
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.metrics import roc_auc_score,average_precision_score
BASE=Path(__file__).resolve().parents[1];REV=BASE.parent
SEED=20261103

def out(rows,path):pd.DataFrame(rows).to_csv(path,sep='\t',index=False,na_rep='NA',float_format='%.15g')
def read(path):return pd.read_csv(path,sep='\t')
def range_feature(f):
 try:return int(f['location']['start']['value']),int(f['location']['end']['value'])
 except (KeyError,TypeError,ValueError):return None

def region(record,typ,desc=None):
 return [range_feature(f) for f in record.get('features',[]) if f['type']==typ and (desc is None or desc.lower() in f.get('description','').lower()) and range_feature(f)]
def seq_regions(record,ranges):return ''.join(record['sequence']['value'][a-1:b] for a,b in ranges)
def match_summary(a,b):
 if not a or not b:return dict(status='NO_REGION_ANNOTATION')
 align=PairwiseAligner(mode='global');align.substitution_matrix=substitution_matrices.load('BLOSUM62');align.open_gap_score=-10;align.extend_gap_score=-.5
 if set(a+b)-set('ARNDCQEGHILKMFPSTWYVBZX'):return dict(status='UNSUPPORTED_RESIDUE')
 z=align.align(a,b)[0];paired=ident=0
 for (p,q),(r,s) in zip(z.aligned[0],z.aligned[1]):paired+=q-p;ident+=sum(x==y for x,y in zip(a[p:q],b[r:s]))
 columns=len(a)+len(b)-paired
 return dict(status='OK',human_length=len(a),pig_length=len(b),length_ratio=len(b)/len(a),n_aligned_pairs=paired,alignment_columns=columns,n_identical=ident,identity_paired_pct=100*ident/paired if paired else np.nan,identity_all_columns_pct=100*ident/columns,coverage_human=paired/len(a),coverage_pig=paired/len(b),gap_fraction=(columns-paired)/columns)

def seq_qc(dest):
 directory=BASE/'inputs/public';rows=[];motifs=[]
 selection=read(BASE/'results/sequence_input_selection.tsv')
 for r in selection.to_dict('records'):
  gene=r['gene'];paths={s:directory/f'{gene}_{s}.json' for s in ['pig','human']}
  if not all(p.exists() for p in paths.values()):
   rows.append(dict(gene=gene,region='full_length',status='MISSING_OR_AMBIGUOUS_INPUT'));continue
  records={s:json.loads(p.read_text()) for s,p in paths.items()}
  if records['pig']['organism']['taxonId']!=9823 or records['human']['organism']['taxonId']!=9606:raise ValueError('species mismatch')
  for mode in ['full_length','annotated_extracellular','mature_secreted']:
   ranges={s:region(d,'Topological domain','Extracellular') if mode=='annotated_extracellular' else region(d,'Chain') if mode=='mature_secreted' and not region(d,'Transmembrane') and any('Secreted' in str(c) for c in d.get('comments',[])) else [(1,d['sequence']['length'])] if mode=='full_length' else [] for s,d in records.items()}
   summary=match_summary(seq_regions(records['human'],ranges['human']),seq_regions(records['pig'],ranges['pig']))
   rows.append(dict(gene=gene,region=mode,pig_accession=records['pig']['primaryAccession'],human_accession=records['human']['primaryAccession'],pig_ranges=json.dumps(ranges['pig']),human_ranges=json.dumps(ranges['human']),pig_sequence_sha256=hashlib.sha256(records['pig']['sequence']['value'].encode()).hexdigest(),human_sequence_sha256=hashlib.sha256(records['human']['sequence']['value'].encode()).hexdigest(),**summary))
  # Additional sensitivity transfers human extracellular positions through alignment only;
  # explicitly INFERRED, never treated as pig experimental topology.
  hs=records['human']['sequence']['value'];ps=records['pig']['sequence']['value'];ext=region(records['human'],'Topological domain','Extracellular')
  if ext:
   align=PairwiseAligner(mode='global');align.substitution_matrix=substitution_matrices.load('BLOSUM62');align.open_gap_score=-10;align.extend_gap_score=-.5;z=align.align(hs,ps)[0];mapped=[]
   for (a,b),(c,d) in zip(z.aligned[0],z.aligned[1]):
    mapped.extend((i,j) for i,j in zip(range(a,b),range(c,d)) if any(lo<=i+1<=hi for lo,hi in ext))
   if mapped:rows.append(dict(gene=gene,region='human_ECD_transferred_inference',status='INFERRED_not_pig_annotation',human_length=sum(b-a+1 for a,b in ext),pig_length=len(mapped),n_aligned_pairs=len(mapped),n_identical=sum(hs[i]==ps[j] for i,j in mapped),identity_paired_pct=100*sum(hs[i]==ps[j] for i,j in mapped)/len(mapped),coverage_human=len(mapped)/sum(b-a+1 for a,b in ext),pig_accession=records['pig']['primaryAccession'],human_accession=records['human']['primaryAccession']))
  for species,d in records.items():
   seq=d['sequence']['value'];regions={'signal':region(d,'Signal'),'transmembrane':region(d,'Transmembrane'),'extracellular':region(d,'Topological domain','Extracellular'),'cytoplasmic':region(d,'Topological domain','Cytoplasmic')}
   for i in range(len(seq)-2):
    motif=seq[i:i+3]
    if motif[0]!='N' or motif[1]=='P' or motif[2] not in 'STC':continue
    tags=[kind for kind,rr in regions.items() if any(a<=i+1 and i+3<=b for a,b in rr)]
    motifs.append(dict(gene=gene,species=species,accession=d['primaryAccession'],position_1based=i+1,motif=motif,motif_class='canonical_NXS_T' if motif[2] in 'ST' else 'noncanonical_NXC',topology=';'.join(tags) if tags else 'other_or_unknown',experimentally_occupied='NOT_ESTABLISHED',annotation_release=d.get('entryAudit',{}).get('lastAnnotationUpdateDate')))
 out(rows,dest/'sequence_domain_qc.tsv');out(motifs,dest/'glyco_topology_motifs.tsv')
 m=pd.DataFrame(motifs);out(m.groupby(['gene','species','motif_class','topology']).size().reset_index(name='n_motifs'),dest/'glyco_topology_summary.tsv')

def organ(dest):
 t=read(BASE/'inputs/tier1_frozen.tsv');m=read(BASE/'inputs/axis_donor_masks.tsv');scores=read(REV/'results/axis_scores_long.tsv');rows=[];members=[]
 for organ in ['kidney','heart','liver']:
  tags=t.organs.str.replace(' ','').str.split(',');included=tags.apply(lambda z:bool(set(z)&{organ,'all','immune','vascular'}))
  for i,r in t.iterrows():members.append(dict(organ=organ,pair_id=r.pair_id,organs_recorded=r.organs,included_annotation_subset=bool(included[i]),rule='explicit_organ_or_all_immune_vascular',scope='not_in_subset_is_NOT_biologically_absent'))
  for variant,ss in scores.groupby('variant',sort=True):
   z=t.loc[included,['pair_id','evidence','identity']].merge(ss,on='pair_id',validate='one_to_one');z['score']=pd.to_numeric(z.score,errors='coerce')
   for donor,mask in m.groupby('donor',sort=True):
    q=z.merge(mask[['pair_id','mask']],on='pair_id',validate='one_to_one');valid=q.score.notna();T=q.loc[valid,'score'].sum();A=(q.loc[valid,'score']*q.loc[valid,'mask']).sum();R=T-A
    rows.append(dict(organ=organ,variant=variant,donor=donor,n_axes=len(q),n_noncontrol=(q.evidence!='control').sum(),n_scored=int(valid.sum()),n_missing_noncontrol=int((~valid & q.evidence.ne('control')).sum()),n_original_sequence_missing_noncontrol=int((q.identity.isna() & q.evidence.ne('control')).sum()),n_numeric_despite_missing_sequence_noncontrol=int((q.identity.isna() & valid & q.evidence.ne('control')).sum()),T_scored=T,A_scored=A,R_scored=R,C_scored=A/T if T else np.nan,count_targeted_noncontrol=int(((q.evidence!='control')&(q['mask']==1)).sum()),count_coverage_noncontrol=((q.evidence!='control')&(q['mask']==1)).sum()/(q.evidence!='control').sum(),full_organ_burden=np.nan,interpretation='annotation_stratum_sensitivity_NOT_expression_weighted_risk'))
 out(rows,dest/'organ_annotation_scorecards.tsv');out(members,dest/'organ_annotation_membership.tsv')
 out([dict(source='historical_human_organ_expr_UNVERIFIED.tsv',status='NOT_USED',reason='no_accession_version_sample_identifiers_or_source_matrix_lineage;human_table_contains_GGTA1_label;cannot_treat_as_verified_GTEx_TPM')],dest/'organ_expression_provenance_audit.tsv')

def stats(dest):
 c=read(REV/'phase2/results/isi_common_cohort.tsv');rows=[]
 splits=StratifiedShuffleSplit(n_splits=500,test_size=.2,random_state=42)
 for split,(train,test) in enumerate(splits.split(c,c.label)):
  for model in ['divergence','severity_only','ISI_v2','formula','formula_mixed3','no_expression','no_divergence']:
   y=c.label.iloc[test];s=np.round(c[model].iloc[test],12);rows.append(dict(split=split,model=model,n_test=len(test),test_ids=';'.join(c.pair_id.iloc[test]),auroc=roc_auc_score(y,s),ap=average_precision_score(y,s),method='fixed_score_resampling_NOT_independent_validation_or_exact_legacy62_replay'))
 out(rows,dest/'isi_500_split_diagnostics.tsv');d=pd.DataFrame(rows);out(d.groupby('model').agg(n_splits=('split','count'),mean_auroc=('auroc','mean'),sd_auroc=('auroc','std'),mean_ap=('ap','mean'),sd_ap=('ap','std')).reset_index(),dest/'isi_500_split_summary.tsv')
 a=read(REV/'phase2/results/antigen_priority_all.tsv');known=a[a.historical_positive];background=a[~a.historical_positive];rows=[];rng=np.random.default_rng(SEED)
 for cohort,genes in [('exclude_unverified_CD55',[g for g in known.gene if g!='CD55']),('human_antibody_supported6',['CD46','CD59','PROCR','CD9','CFTR','CDH1'])]:
  p=a[a.gene.isin(genes)];u=mannwhitneyu(p.combined,background.combined,alternative='greater',method='asymptotic');ps=u.statistic/(len(p)*len(background));boot=[]
  for _ in range(2000):
   pp=rng.choice(p.combined,len(p));bb=rng.choice(background.combined,len(background));boot.append(np.mean(pp[:,None]>bb)+.5*np.mean(pp[:,None]==bb))
  null=np.array([rng.choice(background.combined,len(p),replace=False).mean() for _ in range(10000)])
  rows.append(dict(cohort=cohort,n_positive=len(p),n_background=len(background),probability_superiority=ps,bootstrap_ci_low=np.quantile(boot,.025),bootstrap_ci_high=np.quantile(boot,.975),mwu_U=u.statistic,mwu_asymptotic_p=u.pvalue,conditional_background_mean_p=(1+(null>=p.combined.mean()).sum())/10001,interpretation='posthoc_source_sensitivity_mixed_assays_NOT_new_independent_validation'))
 out(rows,dest/'antigen_evidence_sensitivity.tsv')
 # The two extra cohorts are exploratory and two endpoints each, so report one
 # extended family containing all 33 earlier tests and 4 new tests.
 previous=read(REV/'phase2/results/antigen_benchmark_metrics.tsv');family=[]
 for r in previous[previous.p_value.notna()].to_dict('records'):family.append(dict(origin='phase2',comparison=r['scope']+'|'+r['background']+'|'+r['model'],p_value=r['p_value']))
 for r in rows:
  for k in ['mwu_asymptotic_p','conditional_background_mean_p']:family.append(dict(origin='phase3',comparison=r['cohort']+'|'+k,p_value=r[k]))
 f=pd.DataFrame(family);order=np.argsort(f.p_value);p=f.p_value.to_numpy();q=np.empty(len(p));q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);f['q_BH_extended_family']=q;out(f,dest/'antigen_extended_multiplicity.tsv')

def main():
 a=argparse.ArgumentParser();a.add_argument('--output',type=Path,required=True);dest=a.parse_args().output;dest.mkdir(exist_ok=True,parents=True)
 for name,fn in [('organ',organ),('sequence',seq_qc),('stats',stats)]:print(name,flush=True);fn(dest)
 print('complete',flush=True)
if __name__=='__main__':main()
