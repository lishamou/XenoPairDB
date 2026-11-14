"""Retrospective assay-conditional assessment, not clinical validation.
Frozen scores/labels are read-only. Cell sources, not wells, are Anand units.
"""
from pathlib import Path
import argparse, hashlib, itertools, json, math, re, sys, platform
import numpy as np
import pandas as pd
import scipy
from scipy import stats
import openpyxl

P=argparse.ArgumentParser();P.add_argument('--input',type=Path,required=True);P.add_argument('--output',type=Path,required=True)
args=P.parse_args();I=args.input;O=args.output;O.mkdir(parents=True,exist_ok=True)
SEED=20261104;BOOT=10000;PERM=99999
def save(d,n):
    pd.DataFrame(d).to_csv(O/(n+'.tsv'),sep='\t',index=False,float_format='%.12g',na_rep='NA')
def clean(x):
    if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [clean(v) for v in x]
    if isinstance(x,np.generic):x=x.item()
    if isinstance(x,float) and not np.isfinite(x):return None
    return x
def js(d,n): (O/(n+'.json')).write_text(json.dumps(clean(d),indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def sheet(file,name):return openpyxl.load_workbook(I/'public'/file,data_only=True,read_only=True)[name]
def rows(file,name,header):
    vals=list(sheet(file,name).values);names=vals[header-1]
    return [dict(zip(names,r),source_row=j) for j,r in enumerate(vals[header:],header+1) if r[0] is not None]
def num(x):
    try:v=float(x);return v if np.isfinite(v) else np.nan
    except (ValueError,TypeError):return np.nan
def vals(x):return np.array([num(v) for v in str(x).split('|')])
def positive_count(x):return int(np.sum(np.isfinite(vals(x))&(vals(x)>0)))
def normrank(x,axis=None):
    z=stats.rankdata(x,axis=axis);z=z-z.mean(axis=axis,keepdims=True)
    return z/np.sqrt(np.sum(z*z,axis=axis,keepdims=True))
def association(d,name,partial=False):
    rng=np.random.default_rng(SEED);x=d.antigenicity_v2.to_numpy(float);y=d.log2fc.to_numpy(float);n=len(d)
    out=dict(analysis=name,n=n,rho=np.nan,ci_low=np.nan,ci_high=np.nan,p=np.nan,bootstrap_valid=0,
             ci_method='paired_protein_percentile_bootstrap_conditional_on_assay',test='two_sided_Spearman_permutation_plus_one',seed=SEED)
    if n<4 or len(set(x))<2 or len(set(y))<2:return out
    def rho(a,b,c=None):
        if c is None:return float(np.dot(normrank(a),normrank(b)))
        ar=stats.rankdata(a);br=stats.rankdata(b);cr=np.c_[np.ones(len(c)),stats.rankdata(c)]
        ar=ar-cr@np.linalg.lstsq(cr,ar,rcond=None)[0];br=br-cr@np.linalg.lstsq(cr,br,rcond=None)[0]
        if np.linalg.norm(ar)<1e-10 or np.linalg.norm(br)<1e-10:return np.nan
        return float(np.corrcoef(ar,br)[0,1])
    z=d.global_abundance.to_numpy(float) if partial else None
    obs=rho(x,y,z);out['rho']=obs;boots=[]
    if not np.isfinite(obs):
        out['ci_method']='NOT_ESTIMABLE_constant_rank_residual';return out
    for _ in range(BOOT):
        ix=rng.integers(n,size=n)
        if len(set(x[ix]))<2 or len(set(y[ix]))<2:continue
        v=rho(x[ix],y[ix],z[ix] if partial else None)
        if np.isfinite(v):boots.append(v)
    if not boots:
        out.update(ci_method='NOT_ESTIMABLE_all_bootstraps_degenerate',bootstrap_degenerate=BOOT);return out
    out.update(ci_low=float(np.quantile(boots,.025)),ci_high=float(np.quantile(boots,.975)),bootstrap_valid=len(boots))
    if partial:
        out['test']='partial_Spearman_rank_residuals_no_permutation_p'
        out['bootstrap_degenerate']=BOOT-len(boots)
        if len(boots)<.95*BOOT:
            out.update(ci_low=np.nan,ci_high=np.nan,ci_method='NOT_VALIDATED_more_than_5pct_degenerate_bootstraps')
        return out
    xr=normrank(x);yr=normrank(y);ext=0
    for _ in range(PERM):ext+=abs(float(np.dot(xr,rng.permutation(yr))))>=abs(obs)-1e-12
    out['p']=(ext+1)/(PERM+1);return out

frozen=pd.read_csv(I/'prior_snapshot/phase2/results/antigen_priority_all.tsv',sep='\t')
assert frozen.gene.is_unique and frozen.acc.is_unique
F=frozen.set_index('gene');accmap=dict(zip(frozen.acc,frozen.gene));genes=set(F.index)
def mapping(acc,gene):
    a=str(acc).strip();g='' if gene is None else str(gene).strip()
    if any(x in a for x in [';',',',' ']) or any(x in g for x in [';',',',' ']):return None,'ambiguous_source_group'
    if a in accmap:
        target=accmap[a]
        if g in genes and g!=target:return None,'accession_gene_conflict'
        return target,'exact_accession'
    if g in genes:return g,'unique_explicit_gene'
    return None,'unmapped'

s4='Supplementary Table S4.xlsx';peps=[]
for r in rows(s4,'Peptide_level',1):
    g,m=mapping(r['PG.ProteinAccessions'],r['PG.Genes']);fc=num(r['FoldChange'])
    peps.append(dict(source_file=s4,source_sheet='Peptide_level',source_row=r['source_row'],source_accession=r['PG.ProteinAccessions'],
        source_gene=r['PG.Genes'],gene=g,mapping=m,peptide=r['PEP.StrippedSequence'],fold_change=fc,
        log2fc=math.log2(fc) if fc>0 else np.nan,BFDR=num(r['BFDR']),PK_intensity=r['Intensity'],E_intensity=r['ctrlIntensity'],
        observed_PK=positive_count(r['Intensity']),observed_E=positive_count(r['ctrlIntensity'])))
peps=pd.DataFrame(peps)
shared=peps.dropna(subset=['gene']).groupby('peptide').gene.nunique();shared=set(shared[shared>1].index)
peps['shared_across_mapped_genes']=peps.peptide.isin(shared)
peps['eligible_model']=peps.gene.notna()&~peps.shared_across_mapped_genes&np.isfinite(peps.log2fc)&(peps.observed_PK>0)
peps['eligible_observed']=peps.eligible_model&(peps.observed_E>0)
peps['exclusion']=np.select([peps.gene.isna(),peps.shared_across_mapped_genes,~np.isfinite(peps.log2fc),peps.observed_PK<1,peps.observed_E<1],
    [peps.mapping,'shared_peptide','nonfinite_FC','no_observed_PK','no_observed_E'],default='eligible')
s3rows=list(sheet('Supplementary Table S3.xlsx','Sheet6').values)
enriched=set(r[0] for r in s3rows[2:] if r[0]);peps['author_enriched']=peps.peptide.isin(enriched)
save(peps,'liu_peptide_mapping')
def aggregate(d):
    # One value per stripped peptide, then one per gene: no weighting by duplicate rows.
    q=d.groupby(['gene','peptide'],as_index=False).agg(log2fc=('log2fc','median'),author_enriched=('author_enriched','max'))
    q=q.groupby('gene',as_index=False).agg(log2fc=('log2fc','median'),n_peptides=('peptide','size'),source_enriched=('author_enriched','max'))
    return q.merge(frozen[['gene','acc','antigenicity_v2','rank_v2','historical_positive']],on='gene',validate='one_to_one')
all_observed=aggregate(peps[peps.eligible_observed]);primary=all_observed[~all_observed.historical_positive].copy()
exact=aggregate(peps[peps.eligible_observed&(peps.mapping=='exact_accession')]);exact=exact[~exact.historical_positive]
complete=aggregate(peps[peps.eligible_observed&(peps.observed_E==3)]);complete=complete[~complete.historical_positive]
model=aggregate(peps[peps.eligible_model]);model=model[~model.historical_positive]
save(all_observed,'liu_protein_observed');save(primary,'liu_primary_cohort');save(exact,'liu_exact_accession_cohort');save(complete,'liu_complete_control_cohort');save(model,'liu_source_model_cohort')
metrics=[association(primary,'primary_exclude_legacy8'),association(all_observed,'include_legacy8'),association(exact,'exact_accession_only'),
         association(complete,'all_three_controls_observed'),association(model,'include_source_model_control_missing')]
protein=[]
for r in rows(s4,'Protein_level',2):
    g,m=mapping(r['PG.ProteinAccessions'],r['PG.Genes']);fc=num(r['FoldChange'])
    protein.append(dict(source_row=r['source_row'],source_accession=r['PG.ProteinAccessions'],source_gene=r['PG.Genes'],gene=g,mapping=m,fold_change=fc,
        log2fc=math.log2(fc) if fc>0 else np.nan,observed_PK=positive_count(r['Intensity']),observed_E=positive_count(r['ctrlIntensity'])))
protein=pd.DataFrame(protein);save(protein,'liu_protein_level_mapping')
prot=protein[protein.gene.notna()&(protein.observed_E>0)&(protein.observed_PK>0)&np.isfinite(protein.log2fc)].groupby('gene',as_index=False).log2fc.median().merge(frozen,on='gene',validate='one_to_one')
prot=prot[~prot.historical_positive];save(prot,'liu_protein_level_cohort');metrics.append(association(prot,'protein_level_observed_controls'))

# Abundance: exact accession or a unique source-accession link established above, no fuzzy aliases.
source_links={}
for d in [peps,protein]:
    for r in d.dropna(subset=['gene']).itertuples():source_links.setdefault(r.source_accession,set()).add(r.gene)
ab=[];s5='Supplementary Table S5.xlsx'
for r in rows(s5,'Pig kidney proteins',3):
    a=str(r['PG.ProteinAccessions']);g=accmap.get(a);how='exact_accession'
    if g is None and len(source_links.get(a,set()))==1:g=next(iter(source_links[a]));how='unique_source_accession_gene_link'
    if ';' in a:g=None
    quantities=[num(v) for k,v in r.items() if isinstance(k,str) and k.endswith('.PG.Quantity')]
    q=[v for v in quantities if np.isfinite(v) and v>0]
    ab.append(dict(source_row=r['source_row'],source_accession=a,gene=g,mapping=how if g else 'unmapped_or_group',global_abundance=np.mean(q) if q else np.nan,n_observed=len(q)))
ab=pd.DataFrame(ab);save(ab,'liu_global_abundance_mapping')
aggab=ab.dropna(subset=['gene','global_abundance']).groupby('gene',as_index=False).global_abundance.median()
adj=primary.merge(aggab,on='gene',validate='one_to_one');save(adj,'liu_abundance_complete_cases')
metrics.extend([association(adj,'abundance_complete_cases_unadjusted'),association(adj,'abundance_adjusted_partial',True)])

threshold=int(math.floor(len(frozen)*.10));primary['top_decile']=primary.rank_v2<=threshold
ct=np.array([[sum(primary.top_decile&primary.source_enriched),sum(primary.top_decile&~primary.source_enriched)],
             [sum(~primary.top_decile&primary.source_enriched),sum(~primary.top_decile&~primary.source_enriched)]],dtype=int)
od=stats.contingency.odds_ratio(ct);ci=od.confidence_interval()
enrich=dict(top_rank_threshold=threshold,contingency=ct.tolist(),rows=['frozen_global_top_decile','remaining_measured'],columns=['author_enriched','not_author_selected'],
    odds_ratio=od.statistic,ci_low=ci.low,ci_high=ci.high,p=stats.fisher_exact(ct).pvalue,n=len(primary),method='conditional_exact_OR_CI_two_sided_Fisher',selection='any_eligible_peptide_in_final_S3_Sheet6')
enrich['status']='estimable'
if np.any(ct.sum(axis=0)==0) or np.any(ct.sum(axis=1)==0):
    enrich.update(status='NOT_ESTIMABLE_zero_marginal_no_comparison',odds_ratio=np.nan,ci_low=np.nan,ci_high=np.nan,p=np.nan)
indices=[j for j,m in enumerate(metrics) if j>0 and np.isfinite(m['p'])];has_enrichment=np.isfinite(enrich['p']);pv=[metrics[j]['p'] for j in indices]+([enrich['p']] if has_enrichment else [])
q=stats.false_discovery_control(pv)
for j,v in zip(indices,q):metrics[j]['q_BH_secondary']=v
enrich['q_BH_secondary']=q[-1] if has_enrichment else np.nan;enrich['secondary_family_size']=len(pv)
save(metrics,'liu_associations');js(enrich,'liu_top_decile_enrichment')
gly=[]
for sn in ['PK1','PK2','PK3']:
    ss=list(sheet(s5,sn).values)
    for j,row in enumerate(ss,1):
        if not isinstance(row[0],str) or '|' not in row[0]:continue
        ac=re.findall(r'(?:sp|tr)\|([^|]+)\|',row[0]);gm=re.findall(r'GN=([^\s]+)',row[0]);g,m=mapping(ac[0] if len(ac)==1 else ';'.join(ac),gm[0] if len(gm)==1 else ';'.join(gm))
        gly.append(dict(source_sheet=sn,source_row=j,source_proteins=row[0],peptide=row[1],modified_peptide=row[2],glycan=row[4],gene=g,mapping=m,
            frozen_rank=F.loc[g,'rank_v2'] if g in genes else np.nan,legacy_positive=bool(F.loc[g,'historical_positive']) if g in genes else None))
save(gly,'liu_glycopeptide_concordance')
flow=dict(frozen_rankable=len(frozen),historical_positive=int(frozen.historical_positive.sum()),S4_peptide_rows=len(peps),S4_protein_rows=len(protein),
    S3_author_enriched_unique_peptides=len(enriched),mapping_counts=peps.mapping.value_counts().to_dict(),primary_exclusion_counts=peps.exclusion.value_counts().to_dict(),
    observed_eligible_peptide_rows=int(peps.eligible_observed.sum()),mapped_observed_genes=len(all_observed),legacy_in_observed=all_observed.loc[all_observed.historical_positive,'gene'].tolist(),
    primary_n=len(primary),primary_unique_peptides=int(primary.n_peptides.sum()),exact_accession_n=len(exact),complete_control_n=len(complete),model_n=len(model),abundance_complete_n=len(adj),
    S4_controls_all_missing=int((peps.observed_E==0).sum()),glycopeptide_rows=len(gly),glycopeptide_mapped_genes=sorted(set(r['gene'] for r in gly if r['gene'])),
    scope='conditional_on_final_S4_assay_universe_NOT_4122_ground_truth_negatives',score_column='antigenicity_v2',rank_column='rank_v2')
js(flow,'liu_flow')

# Anand source values and coordinates.
M3='41586_2023_6594_MOESM3_ESM.xlsx';M5='41586_2023_6594_MOESM5_ESM.xlsx';M9='41586_2023_6594_MOESM9_ESM.xlsx'
points=[]
def pt(file,sn,row,col,cid,gen,endpoint,species,time,condition,value,unit):
    if not np.isfinite(num(value)):return
    points.append(dict(file=file,sheet=sn,cell=f'{openpyxl.utils.get_column_letter(col)}{row}',cell_source=str(cid),genotype=gen,
                       endpoint=endpoint,species=species,time=time,condition=condition,value=num(value),unit=unit))
for sn,species in [('3a','human'),('3b','cyno')]:
    for j,r in enumerate(sheet(M5,sn).values,1):
        if j<3 or r[1] is None:continue
        for c,cond in [(3,'serum'),(4,'EDTA')]:pt(M5,sn,j,c,r[0],r[1],'C3b',species,'30min',cond,r[c-1],'MFI_x1000')
for j,r in enumerate(sheet(M5,'3c').values,1):
    if j<3 or r[1] is None:continue
    for c,sp in [(3,'human'),(4,'cyno')]:pt(M5,'3c',j,c,r[0],r[1],'CDC',sp,'45min','serum',r[c-1],'percent_CytotoxRed_positive')
for j,r in enumerate(sheet(M9,'ED Fig. 4b').values,1):
    if j<4 or r[1] is None:continue
    for c,sp,t in [(3,'human','45min'),(4,'cyno','45min'),(5,'human','15h'),(6,'cyno','15h')]:pt(M9,'ED Fig. 4b',j,c,r[0],r[1],'CDC_ED4',sp,t,'serum',r[c-1],'percent_CytotoxRed_positive')
for j,r in enumerate(sheet(M5,'3d').values,1):
    if j<6 or r[1] is None:continue
    for c,cond in [(3,'unblocked'),(4,'CD46_block'),(5,'CD55_block'),(6,'both_block'),(7,'EDTA')]:pt(M5,'3d',j,c,r[0],r[1],'C3b_blocking','human','30min',cond,r[c-1],'MFI_x1000')
points=pd.DataFrame(points);save(points,'anand_source_points')
binding=[]
for j,r in enumerate(sheet(M3,'1b').values,1):
    if j<2:continue
    for antibody,a,b in [('IgG',2,5),('IgM',8,11)]:
        wt=np.array(r[a:a+3],float);ko=np.array(r[b:b+3],float)
        binding.append(dict(species=r[1],antibody=antibody,WT_mean=wt.mean(),KO_mean=ko.mean(),ratio_KO_over_WT=ko.mean()/wt.mean(),
            reduction_percent=100*(1-ko.mean()/wt.mean()),WT_min=wt.min(),WT_max=wt.max(),KO_min=ko.min(),KO_max=ko.max(),
            technical_replicates_per_condition=3,independent_biological_n='not_estimated_from_technical_triplicates',source_file=M3,sheet='1b',
            WT_cells=f'{openpyxl.utils.get_column_letter(a+1)}{j}:{openpyxl.utils.get_column_letter(a+3)}{j}',KO_cells=f'{openpyxl.utils.get_column_letter(b+1)}{j}:{openpyxl.utils.get_column_letter(b+3)}{j}'))
save(binding,'anand_binding_descriptive')
def contrast(a,b,paired=False):
    a=np.array(a,float);b=np.array(b,float);effect=a.mean()-b.mean()
    if paired:
        assert len(a)==len(b);diff=a-b;se=stats.sem(diff);df=len(diff)-1
        ps=[np.mean(diff*np.array(sign)) for sign in itertools.product([-1,1],repeat=len(diff))]
        method='exact_paired_sign_flip';cm='paired_t_model_dependent'
    else:
        va=a.var(ddof=1)/len(a);vb=b.var(ddof=1)/len(b);se=np.sqrt(va+vb);df=(va+vb)**2/(va*va/(len(a)-1)+vb*vb/(len(b)-1))
        z=np.r_[a,b];ps=[]
        for ix in itertools.combinations(range(len(z)),len(a)):
            mask=np.zeros(len(z),bool);mask[list(ix)]=True;ps.append(z[mask].mean()-z[~mask].mean())
        method='exact_label_permutation';cm='Welch_t_model_dependent'
    rad=stats.t.ppf(.975,df)*se
    return dict(n_treatment=len(a),n_reference=len(b),mean_treatment=a.mean(),mean_reference=b.mean(),difference=effect,
                ci_low=effect-rad,ci_high=effect+rad,df=df,p=sum(abs(v)>=abs(effect)-1e-12 for v in ps)/len(ps),
                permutations=len(ps),test=method+'_two_sided',ci_method=cm)
comparisons=[]
for (endpoint,species,time),d in points[points.condition=='serum'].groupby(['endpoint','species','time'],sort=True):
    for tg in ['3KO.7TG','3KO.7TG.RI']:
        a=d[d.genotype==tg];b=d[d.genotype=='3KO']
        if len(a)<2 or len(b)<2:continue
        comparisons.append(dict(comparison=f'{endpoint}_{species}_{time}_{tg}_minus_3KO',endpoint=endpoint,species=species,time=time,treatment=tg,reference='3KO',
            treatment_IDs=';'.join(a.cell_source),reference_IDs=';'.join(b.cell_source),unit=d.unit.iloc[0],**contrast(a.value,b.value)))
for tg in ['3KO.7TG','3KO.7TG.RI']:
    d=points[(points.endpoint=='C3b_blocking')&(points.genotype==tg)].pivot(index='cell_source',columns='condition',values='value')
    for condition in ['CD46_block','CD55_block','both_block']:
        comparisons.append(dict(comparison=f'C3b_blocking_{tg}_{condition}_minus_unblocked',endpoint='C3b_blocking',species='human',time='30min',
            treatment=condition,reference='unblocked',treatment_IDs=';'.join(d.index),reference_IDs=';'.join(d.index),unit='MFI_x1000',**contrast(d[condition],d.unblocked,True)))
cp=pd.DataFrame(comparisons);order=np.argsort(cp.p.to_numpy());adjp=np.minimum(1,np.maximum.accumulate(cp.p.to_numpy()[order]*(len(cp)-np.arange(len(cp)))))
cp['p_Holm']=np.nan;cp.loc[order,'p_Holm']=adjp;save(cp,'anand_contrasts')
js(dict(seed=SEED,bootstrap_draws=BOOT,permutations=PERM,python=sys.version,platform=platform.platform(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__,openpyxl=openpyxl.__version__,
    A1_primary='primary_exclude_legacy8',A1_secondary_family=len(pv),A2_Holm_family=len(cp),A2_ci_caution='Normality_unverifiable_at_n2_3_CI_and_exact_p_use_different_assumptions',
    input_hashes={str(f.relative_to(I)):hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(I.rglob('*')) if f.is_file()},
    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),'run_metadata')
print('Completed public data analysis:',flow['primary_n'],'Liu primary proteins;',len(cp),'Anand contrasts.')
