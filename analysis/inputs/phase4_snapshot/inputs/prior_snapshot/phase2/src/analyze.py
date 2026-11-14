"""Exploratory fixed-catalogue sensitivity; never independent prediction validation."""
from pathlib import Path
import json, csv, re, argparse, itertools
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

BASE=Path(__file__).resolve().parents[1]
SEED=20261103; B=2000; MC=10000

def table(name): return pd.read_csv(BASE/'inputs'/name,sep='\t',keep_default_na=True)
def out(df,path):
    path.parent.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(df).to_csv(path,sep='\t',index=False,na_rep='NA',float_format='%.15g')
def pct(x): return pd.Series(x).rank(method='average',pct=True).to_numpy()
def auc(y,s):
    y=np.asarray(y,dtype=int);s=np.round(np.asarray(s,float),12)
    if not np.isfinite(s).all(): raise ValueError('missing score')
    n=y.sum();m=len(y)-n
    if not n or not m:return np.nan
    return float((rankdata(s)[y==1].sum()-n*(n+1)/2)/(n*m))
def ap(y,s):
    y=np.asarray(y); s=np.round(np.asarray(s,float),12)
    if not y.sum() or y.sum()==len(y): return np.nan
    order=np.argsort(-s,kind='stable');y=y[order];s=s[order]
    ends=np.r_[np.where(s[:-1]!=s[1:])[0],len(s)-1];tp=np.cumsum(y)[ends]
    return float(np.sum(np.diff(np.r_[0,tp])*tp/(ends+1))/y.sum())
def ci(v):
    a=np.array(v,float);a=a[np.isfinite(a)]
    return tuple(np.quantile(a,[.025,.975])) if len(a) else (np.nan,np.nan)
def components(rows):
    # Shared source or gene => connected block; conservative dependence sensitivity.
    n=len(rows);parent=list(range(n));seen={}
    def find(x):
        while parent[x]!=x: parent[x]=parent[parent[x]];x=parent[x]
        return x
    for i,r in enumerate(rows):
        tokens=['g:'+g for g in re.findall(r'[A-Z][A-Z0-9]+',str(r.get('gene_key',''))+' '+str(r.get('human_partner','')))]
        tokens+=['s:'+s for s in re.findall(r'\d{6,9}',str(r.get('source','')))]
        if r.get('pair_key'): tokens+=['p:'+r['pair_key']]
        for t in tokens:
            if t in seen:parent[find(i)]=find(seen[t])
            else:seen[t]=i
    return np.array([find(i) for i in range(n)])
def metric(y,s,groups,rng):
    y=np.array(y,int);s=np.array(s,float);groups=np.array(groups);uniq=np.unique(groups)
    vals=[]
    for _ in range(B):
        ix=np.concatenate([np.flatnonzero(groups==g) for g in rng.choice(uniq,len(uniq),replace=True)])
        if len(np.unique(y[ix]))==2: vals.append([auc(y[ix],s[ix]),ap(y[ix],s[ix])])
    a=np.array(vals).reshape(-1,2)
    return dict(n=len(y),n_positive=int(y.sum()),n_reference=int((1-y).sum()),n_clusters=len(uniq),auroc=auc(y,s),auroc_ci_low=ci(a[:,0])[0],auroc_ci_high=ci(a[:,0])[1],ap=ap(y,s),ap_prevalence_baseline=y.mean(),ap_ci_low=ci(a[:,1])[0],ap_ci_high=ci(a[:,1])[1],bootstrap_valid=len(a),bootstrap_requested=B,p_value=np.nan,method='connected_gene_source_cluster_percentile_bootstrap',interpretation='annotation_consistency_not_independent_validation')
def coverage(scores,masks):
    s=np.asarray(scores,float);m=np.asarray(masks,float)
    if not np.isfinite(s).all() or not np.isfinite(m).all():raise ValueError('missing not zero')
    if (s<0).any() or (m<0).any() or (m>1).any():raise ValueError('range')
    T=float(s.sum());A=float(np.dot(s,m));R=float(np.dot(s,1-m))
    return T,A,R,A/T if T else np.nan

def isi(dest):
    t=table('tier1_frozen.tsv');
    for c in ['identity','func_sev','evid_w','expr_weight','ISI_v2','iface_id']:t[c]=pd.to_numeric(t[c],errors='coerce')
    if t.pair_id.duplicated().any():raise ValueError('duplicate ID')
    t['divergence']=100-t.identity;t['formula']=t.divergence*t.func_sev*t.evid_w*t.expr_weight
    t['no_expression']=t.divergence*t.func_sev*t.evid_w;t['severity_only']=t.func_sev
    t['no_divergence']=t.func_sev*t.evid_w*t.expr_weight
    t['mixed3_divergence']=np.where(t.iface_id.notna(),(200-t.identity-t.iface_id)/2,t.divergence)
    t['formula_mixed3']=t.mixed3_divergence*t.func_sev*t.evid_w*t.expr_weight
    common=t.dropna(subset=['formula','ISI_v2','no_expression','severity_only']).copy()
    common['label']=(common.evidence!='control').astype(int)
    out(common,dest/'isi_common_cohort.tsv');metrics=[];ranks=[]
    for scope,z in [('common',common),('exclude_disputed',common[common.evidence!='disputed']),('exclude_CFH',common[common.pair_id!='XP075'])]:
        groups=components(z.to_dict('records'))
        for name in ['divergence','severity_only','ISI_v2','formula','formula_mixed3','no_expression','no_divergence']:
            metrics.append(dict(scope=scope,model=name,**metric(z.label,z[name],groups,np.random.default_rng(SEED))))
    out(metrics,dest/'isi_ablation_metrics.tsv')
    z=common[common.label==1].copy()
    for name in ['divergence','severity_only','ISI_v2','formula','formula_mixed3','no_expression','no_divergence']:
        z[name+'_rank']=z[name].rank(ascending=False,method='min')
        ranks.append(dict(model=name,n=len(z),spearman_vs_formula=spearmanr(z[name],z.formula).statistic,top10_overlap=len(set(z.nlargest(10,name).pair_id)&set(z.nlargest(10,'formula').pair_id))))
    out(z,dest/'isi_axis_ranks.tsv');out(ranks,dest/'isi_rank_stability.tsv')
    mask=table('axis_donor_masks.tsv');donors=sorted(mask.donor.unique());M=mask.pivot(index='pair_id',columns='donor',values='mask').reindex(index=t.pair_id,columns=donors).to_numpy(float).T
    if not np.isin(M,[0,1]).all():raise ValueError('invalid mask')
    known=t.identity.notna().to_numpy();sev=t.func_sev.to_numpy();ew=t.evid_w.to_numpy();expr=t.expr_weight.to_numpy();div=t.divergence.to_numpy();nonctrl=t.evidence.ne('control').to_numpy()
    rows=[];pairs=[]
    scenarios=[('sequence_available',None),('nominal_missing_0',0),('nominal_missing_25',25),('nominal_missing_50',50),('nominal_missing_75',75),('nominal_missing_100',100)]
    for scope,fill in scenarios:
        ix=known if fill is None else np.ones(len(t),bool);d=div[ix] if fill is None else np.where(known,div,fill)
        for power,expression,eff in itertools.product([.5,1,2],['actual','neutral'],[.5,1]):
            s=ew[ix]*sev[ix]**power*d*(expr[ix] if expression=='actual' else 1)
            current=[]
            for j,donor in enumerate(donors):
                T,A,R,C=coverage(s,M[j,ix]*eff)
                row=dict(scenario=scope,severity_power=power,expression=expression,edit_fraction=eff,donor=donor,n_axes=int(ix.sum()),n_noncontrol=int(nonctrl[ix].sum()),missing_original=int((~known).sum()),T=T,A=A,R=R,C=C,units='weighted_divergence_percentage',full_scope_known=False,functional_rescue='NOT_INFERRED')
                rows.append(row);current.append(row)
            for a,b in itertools.combinations(current,2):pairs.append(dict(scenario=scope,severity_power=power,expression=expression,edit_fraction=eff,donor_a=a['donor'],donor_b=b['donor'],delta_R_a_minus_b=a['R']-b['R'],delta_C_a_minus_b=a['C']-b['C']))
    for name,s,ix,units in [('legacy_stored',t.ISI_v2[known].to_numpy(),known,'legacy_stored'),('without_divergence',ew*sev*expr,np.ones(len(t),bool),'severity_weighted_not_ISI')]:
        for j,donor in enumerate(donors):
            T,A,R,C=coverage(s,M[j,ix]);rows.append(dict(scenario=name,severity_power=1,expression='actual',edit_fraction=1,donor=donor,n_axes=int(ix.sum()),n_noncontrol=int(nonctrl[ix].sum()),missing_original=12,T=T,A=A,R=R,C=C,units=units,full_scope_known=False,functional_rescue='NOT_INFERRED'))
    out(rows,dest/'donor_scenario_scorecards.tsv');out(pairs,dest/'donor_pairwise_differences.tsv')
    baseline=[]
    for name,score in [('phase1_mixed3',t.formula_mixed3),('sequence_only',t.formula),('mixed3_CFH_legacy_exception',t.formula_mixed3.where(t.pair_id!='XP075',t.ISI_v2))]:
        for j,donor in enumerate(donors):
            T,A,R,C=coverage(score[known],M[j,known]);baseline.append(dict(model=name,donor=donor,T=T,A=A,R=R,C=C,n_sequence_available=int(known.sum()),n_noncontrol=int((known&nonctrl).sum()),note='two_missing_sequence_controls_zero_severity_not_ranked;12_missing_noncontrol_unresolved'))
    out(baseline,dest/'isi_phase1_baseline_bridge.tsv')
    rng=np.random.default_rng(SEED);sim=[]
    for k in range(B):
        severity=np.where(nonctrl,np.clip(sev+rng.uniform(-.1,.1,len(t)),0,1),0)
        s=ew[known]*severity[known]*div[known]*expr[known];eff=rng.uniform(.5,1,known.sum())
        for j,donor in enumerate(donors):
            T,A,R,C=coverage(s,M[j,known]*eff);sim.append(dict(iteration=k,donor=donor,T=T,A=A,R=R,C=C))
    sim=pd.DataFrame(sim);out(sim,dest/'donor_perturbation_draws.tsv');summary=[]
    for donor,z in sim.groupby('donor',sort=True):
        for val in ['T','A','R','C']: summary.append(dict(donor=donor,quantity=val,median=z[val].median(),sensitivity_p025=z[val].quantile(.025),sensitivity_p975=z[val].quantile(.975),draws=B,interval_type='parameter_sensitivity_NOT_biological_CI'))
    out(summary,dest/'donor_sensitivity_intervals.tsv')
    pair=[];piv=sim.pivot(index='iteration',columns='donor',values='R')
    for a,b in itertools.combinations(donors,2):
        v=piv[a]-piv[b];pair.append(dict(donor_a=a,donor_b=b,median_delta_R=v.median(),sensitivity_low=v.quantile(.025),sensitivity_high=v.quantile(.975),fraction_a_lower=(v<0).mean(),interpretation='parameter_sensitivity_not_probability_of_clinical_superiority'))
    out(pair,dest/'donor_pairwise_sensitivity.tsv')

def normalize(s):
    s=re.sub(r'\([^)]*\)','',str(s)).strip()
    s=re.sub(r'^human','',s).strip()
    return s if re.fullmatch(r'[A-Z][A-Z0-9]*',s) else None

def tier2(dest):
    atlas=table('tier2_atlas.tsv');atlas['pair_key']=['|'.join(sorted([a,b])) for a,b in zip(atlas.ligand,atlas.receptor)]
    atlas['div_percentile']=pct(atlas.pair_divergence);atlas['glyco_percentile']=pct(atlas.glyco_delta_pig_minus_human)
    atlas['score_recomputed']=100*(.6*atlas.div_percentile+.4*atlas.glyco_percentile)
    atlas['delta_stored']=atlas.score_recomputed-atlas.tier2_score
    atlas['glyco_recomputed']=atlas.lig_pig_NxST+atlas.rec_pig_NxST-atlas.lig_human_NxST-atlas.rec_human_NxST
    atlas['div_recomputed']=100-(atlas.ligand_identity_pct+atlas.receptor_identity_pct)/2
    out(atlas,dest/'tier2_score_audit.tsv');t=table('tier1_frozen.tsv');mapping=[]
    for r in t.to_dict('records'):
        a,b=normalize(r['pig_gene']),normalize(r['human_partner']);key='|'.join(sorted([a,b])) if a and b else None
        hits=atlas[atlas.pair_key==key] if key else atlas.iloc[:0]
        row=dict(pair_id=r['pair_id'],pig_gene=r['pig_gene'],human_partner=r['human_partner'],direction=r['direction'],evidence=r['evidence'],gene_key=r['gene_key'],source=r['source'],pair_key=key,n_atlas_rows=len(hits),eligible=bool(len(hits)),reason='gene_pair_only_direction_unvalidated' if len(hits) else ('compound_or_non_gene_endpoint' if not key else 'no_atlas_match'),label=int(r['evidence']!='control'))
        # averaging reciprocal atlas rows avoids choosing direction on outcome.
        for col in ['div_percentile','glyco_percentile','score_recomputed']:row[col]=hits[col].mean()
        mapping.append(row)
    m=pd.DataFrame(mapping);eligible=m[m.eligible].copy();conflicts=set(eligible.groupby('pair_key').label.nunique().loc[lambda x:x>1].index)
    m['conflicting_catalogue_labels']=m.pair_key.isin(conflicts);out(m,dest/'tier2_axis_crosswalk.tsv')
    metrics=[]
    for scope,z in [('all_mapped',m[m.eligible]),('exclude_disputed',m[m.eligible & (m.evidence!='disputed')]),('exclude_conflicting_gene_pairs',m[m.eligible & ~m.conflicting_catalogue_labels])]:
        if len(z)==0:continue
        for w in [0,.2,.4,.6,.8,1]:
            score=w*z.div_percentile+(1-w)*z.glyco_percentile
            r=metric(z.label,score,components(z.to_dict('records')),np.random.default_rng(SEED));r.update(scope=scope,div_weight=w,unique_pairs=z.pair_key.nunique(),direction_validation='NOT_VALIDATED');metrics.append(r)
    out(metrics,dest/'tier2_retrospective_metrics.tsv')
    main=set(atlas.nlargest(100,'score_recomputed').pair_key);sensitivity=[]
    for w in [0,.2,.4,.6,.8,1]:
        score=w*atlas.div_percentile+(1-w)*atlas.glyco_percentile
        top=set(atlas.loc[score.nlargest(100).index,'pair_key']);sensitivity.append(dict(div_weight=w,n=len(atlas),spearman_vs_06=spearmanr(score,atlas.score_recomputed).statistic,top100_rows_unique_pairs=len(top),jaccard_vs_06=len(top&main)/len(top|main)))
    out(sensitivity,dest/'tier2_weight_sensitivity.tsv')
    out([dict(stage='published_raw_pool',n=np.nan,status='NOT_REPRODUCED_missing_locked_raw_pool'),dict(stage='frozen_atlas_rows',n=len(atlas),status='OBSERVED'),dict(stage='unique_unordered_pairs',n=atlas.pair_key.nunique(),status='OBSERVED_not_directional'),dict(stage='all_axis_candidates',n=len(m),status='OBSERVED'),dict(stage='conservative_gene_pair_overlap',n=m.eligible.sum(),status='OBSERVED'),dict(stage='secreted_branch',n=np.nan,status='NOT_VALIDATED_topology_sequence_inputs_missing')],dest/'tier2_filter_funnel.tsv')

def bh(p):
    p=np.asarray(p,float);order=np.argsort(p);v=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];q=np.empty(len(p));q[order]=np.minimum(v,1);return q

def antigen(dest):
    raw=table('surfaceome.tsv');cols=['divergence','glyco_density','antigenicity_v2','p_len','rank_v2']
    for c in cols:raw[c]=pd.to_numeric(raw[c],errors='coerce')
    x=raw.dropna(subset=cols).copy().reset_index(drop=True)
    if x.gene.duplicated().any() or x.acc.duplicated().any():raise ValueError('surfaceome duplicate')
    x['sequence']=np.round(pct(x.divergence),12);x['glyco']=np.round(pct(x.glyco_density),12);x['combined']=np.round((x.sequence+x.glyco)/2,12);x['recomputed_rank']=x.combined.rank(ascending=False,method='min');x['delta_stored']=x.combined-x.antigenicity_v2
    j=json.loads((BASE/'inputs/antigenicity_benchmark_summary.json').read_text());pos=list(j['covered_positives']);ref=list(j['putative_negatives']);x['historical_positive']=x.gene.isin(pos);x['background_status']=np.where(x.historical_positive,'historical_positive_mixed_species_evidence','untested_background_NOT_true_negative')
    out(x,dest/'antigen_priority_all.tsv');out(raw[raw[cols].isna().any(axis=1)],dest/'antigen_unscoreable.tsv')
    back=x[~x.gene.isin(pos)];pools=[];metrics=[];enrich=[];rng=np.random.default_rng(SEED)
    sets={'historical8':pos,'exclude_complement5':[g for g in pos if g not in ['CD46','CD55','CD59']],'direct_human2':['CDH1','CFTR']}
    for scope,genes in sets.items():
        p=x[x.gene.isin(genes)].copy()
        for cut in [.01,.05,.1]:enrich.append(dict(scope=scope,n_positive=len(p),universe=len(x),top_fraction=cut,n_hits=int((p.recomputed_rank<=np.floor(len(x)*cut)).sum()),fraction_hits=(p.recomputed_rank<=np.floor(len(x)*cut)).mean()))
        for bgname,bg in [('full_untested',back),('legacy12_reference_only',x[x.gene.isin(ref)])]:
            for model in ['combined','sequence','glyco']:
                y=np.r_[np.ones(len(p)),np.zeros(len(bg))];s=np.r_[p[model],bg[model]]
                # Stratified bootstrap per protein; limits and source-cluster sensitivity separately.
                boot=[]
                for _ in range(B):boot.append(auc(y,np.r_[rng.choice(p[model],len(p)),rng.choice(bg[model],len(bg))]))
                null=np.array([rng.choice(bg[model],len(p),replace=False).mean() for _ in range(MC)])
                metrics.append(dict(scope=scope,background=bgname,model=model,n_positive=len(p),n_background=len(bg),probability_superiority=auc(y,s),ci_low=ci(boot)[0],ci_high=ci(boot)[1],mean_score_difference=p[model].mean()-bg[model].mean(),p_value=(1+(null>=p[model].mean()).sum())/(MC+1),test='one_sided_background_mean_Monte_Carlo',ci_method='stratified_protein_percentile_bootstrap',interpretation='enrichment_NOT_classifier_accuracy'))
        for mode,tol,lo,hi in [('matched',5,.8,1.25),('strict_matched',2.5,.9,1.1)]:
            matches=[]
            for r in p.to_dict('records'):
                z=back[(abs(back.divergence-r['divergence'])<=tol)&(back.p_len/r['p_len']>=lo)&(back.p_len/r['p_len']<=hi)]
                matches.append(z)
                for b in z.to_dict('records'):pools.append(dict(scope=scope,matching=mode,positive=r['gene'],positive_accession=r['acc'],background=b['gene'] if pd.notna(b['gene']) else 'UNNAMED:'+b['acc'],background_accession=b['acc'],divergence_delta=b['divergence']-r['divergence'],length_ratio=b['p_len']/r['p_len']))
            if any(len(z)==0 for z in matches):
                metrics.append(dict(scope=scope,background=mode,model='all',n_positive=len(p),n_background=0,p_value=np.nan,interpretation='NO_COMPLETE_MATCH'));continue
            for model in ['combined','sequence','glyco']:
                draws=np.column_stack([rng.choice(z[model],MC,replace=True) for z in matches]);actual=p[model].to_numpy();effects=actual-np.array([z[model].mean() for z in matches]);means=draws.mean(axis=1)
                boot=[rng.choice(effects,len(effects)).mean() for _ in range(B)]
                ps=np.mean([np.mean((a>z[model])+.5*(a==z[model])) for a,z in zip(actual,matches)])
                metrics.append(dict(scope=scope,background=mode,model=model,n_positive=len(p),n_background=len(set().union(*[set(z.gene) for z in matches])),pool_min=min(map(len,matches)),pool_max=max(map(len,matches)),probability_superiority=ps,mean_score_difference=effects.mean(),ci_low=ci(boot)[0],ci_high=ci(boot)[1],p_value=(1+(means>=actual.mean()).sum())/(MC+1),test='one_sided_stratified_conditional_background_Monte_Carlo',ci_method='positive_protein_bootstrap_of_matched_mean_difference',null_p025=np.quantile(means,.025),null_p975=np.quantile(means,.975),interpretation='conditional_enrichment_NOT_true_negative_accuracy'))
    metrics=pd.DataFrame(metrics);metrics['ci_target']=np.where(metrics.background.str.contains('matched'),'mean_score_difference','probability_superiority');ok=metrics.p_value.notna();metrics.loc[ok,'q_BH_all_antigen_tests']=bh(metrics.loc[ok,'p_value']);out(metrics,dest/'antigen_benchmark_metrics.tsv');out(pools,dest/'matched_background_membership.tsv');out(enrich,dest/'antigen_top_fraction_recall.tsv')
    # Same-paper positives are dependent; leave-one-publication-out sensitivity.
    source={'CD9':'21119562','CD46':'21119562','CD59':'21119562','PROCR':'21119562','CDH1':'22497511','CFTR':'22497511','FAM234A':'2018_FAM234A','CD55':'source_unresolved'}
    p=x[x.gene.isin(pos)].copy();p['publication_cluster']=p.gene.map(source);out(p,dest/'antigen_positive_ranks.tsv');loo=[]
    for pub in sorted(set(source.values())):
        z=p[p.publication_cluster!=pub];y=np.r_[np.ones(len(z)),np.zeros(len(back))]
        loo.append(dict(excluded_publication=pub,n_positive=len(z),probability_superiority=auc(y,np.r_[z.combined,back.combined]),mean_score=z.combined.mean(),interpretation='leave_one_publication_out_descriptive'))
    out(loo,dest/'antigen_publication_sensitivity.tsv')
    out([dict(enzyme=g,target=t,protein_score_used=False,evidence='historical_known_glycan_axis_not_new_wetlab',constraint='enzyme_sequence_score_is_not_antibody_binding_to_glycan') for g,t in [('GGTA1','alphaGal'),('CMAH','Neu5Gc'),('B4GALNT2','Sda_like')]],dest/'glycan_enzyme_panel.tsv')

def main():
    a=argparse.ArgumentParser();a.add_argument('--out',type=Path,default=BASE/'results');args=a.parse_args()
    for name,fn in [('ISI',isi),('Tier2',tier2),('Antigen',antigen)]:print(name,flush=True);fn(args.out)
    print('analysis complete',flush=True)
if __name__=='__main__':main()
