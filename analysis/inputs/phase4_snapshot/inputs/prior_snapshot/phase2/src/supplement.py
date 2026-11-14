"""Audit supplements; descriptive source triangulation, no label replacement."""
from pathlib import Path
import argparse,json,itertools
import numpy as np
import pandas as pd
from analyze import BASE,out,table,normalize,metric,components,SEED,B,MC,auc,ci

def run(dest):
    t=table('tier1_frozen.tsv');j=json.loads((BASE/'inputs/interface_all_18pairs.json').read_text());records={**j['flagship'],**j['extended']};rows=[]
    for r in t.to_dict('records'):
        a,b=normalize(r['pig_gene']),normalize(r['human_partner']);v=records.get(a,{})
        eligible=bool(v.get('partner')==b and 'iface_identity' in v)
        rows.append(dict(pair_id=r['pair_id'],pig_gene=r['pig_gene'],human_partner=r['human_partner'],direction=r['direction'],gene_key=r['gene_key'],source=r['source'],evidence=r['evidence'],label=int(r['evidence']!='control'),eligible=eligible,interface_identity=v.get('iface_identity') if eligible else np.nan,pdb=v.get('pdb') if eligible else None,legacy_interface_status='historical_summary_coordinates_NOT_revalidated' if eligible else 'no_exact_single_gene_pair_mapping',sequence_identity=r['identity']))
    m=pd.DataFrame(rows);out(m,dest/'isi_interface_crosswalk.tsv');z=m[m.eligible];metrics=[]
    for scope,subset in [('exact_historical_interface_subset',z),('exclude_disputed',z[z.evidence!='disputed'])]:
        for model,col in [('legacy_interface_divergence','interface_identity'),('sequence_divergence','sequence_identity')]:
            values=pd.to_numeric(subset[col],errors='coerce');s=subset[values.notna()];values=values.dropna()
            if len(s): metrics.append(dict(scope=scope,model=model,**metric(s.label,100-values,components(s.to_dict('records')),np.random.default_rng(SEED))))
    out(metrics,dest/'isi_interface_metrics.tsv')
    out([dict(pair_id=r['pair_id'],evidence=r['evidence'],func_sev=r['func_sev'],severity_in_feature=True,control_forces_zero=r['evidence']=='control',labels_are_independent=False,impaired_binding_is_not_label=True,reason='noncontrol_includes_harmful_preserved_recognition;severity_curated_from_related_sources') for r in t.to_dict('records')],dest/'label_feature_leakage.tsv')
    # Actual no-sequence coverage includes known zero controls separately from 12 missing burdens.
    out([dict(pair_id=r['pair_id'],evidence=r['evidence'],missing_sequence=pd.isna(r['identity']),missing_stored_score=pd.isna(r['ISI_v2']),resolution='structural_feature_unknown_zero_severity_control' if r['evidence']=='control' else 'UNRESOLVED_no_imputation_in_main') for r in t[t.identity.isna()].to_dict('records')],dest/'isi_missing_feature_inventory.tsv')
    x=pd.read_csv(dest/'antigen_priority_all.tsv',sep='\t');p=x[x.historical_positive];back=x[~x.historical_positive];pools=[]
    sets={'historical8':list(p.gene),'exclude_complement5':list(p[~p.gene.isin(['CD46','CD55','CD59'])].gene),'direct_human2':['CDH1','CFTR']}
    for scope,genes in sets.items():
        for mode,tol,lo,hi in [('matched',5,.8,1.25),('strict_matched',2.5,.9,1.1)]:
            for r in x[x.gene.isin(genes)].to_dict('records'):
                z=back[(abs(back.divergence-r['divergence'])<=tol)&(back.p_len/r['p_len']>=lo)&(back.p_len/r['p_len']<=hi)]
                pools.append(dict(scope=scope,matching=mode,positive=r['gene'],n_background=len(z),positive_divergence=r['divergence'],background_mean_divergence=z.divergence.mean(),positive_length=r['p_len'],background_mean_length=z.p_len.mean(),mean_abs_divergence_delta=(z.divergence-r['divergence']).abs().mean(),max_abs_divergence_delta=(z.divergence-r['divergence']).abs().max(),mean_length_ratio=(z.p_len/r['p_len']).mean(),status='AVAILABLE' if len(z) else 'NO_MATCH'))
    out(pools,dest/'matched_pool_balance.tsv')
    # Historical rank-percentile fingerprint: never repair upstream stored percentile by hand.
    out([dict(metric='scoreable_universe',value=len(x)),dict(metric='max_absolute_combined_score_difference_vs_frozen',value=x.delta_stored.abs().max()),dict(metric='n_differences_gt_1e_minus8',value=(x.delta_stored.abs()>1e-8).sum()),dict(metric='stored_divergence_percentile_max_error',value=(x.sequence-x.div_pct).abs().max()),dict(metric='stored_glyco_percentile_max_error',value=(x.glyco-x.glyco_pct).abs().max())],dest/'antigen_score_reconciliation.tsv')
    # Publication-stratified positive bootstrap; full background resampled independently.
    pp=pd.read_csv(dest/'antigen_positive_ranks.tsv',sep='\t');groups=sorted(pp.publication_cluster.unique());rng=np.random.default_rng(SEED);boot=[]
    for _ in range(B):
        pos=np.concatenate([pp[pp.publication_cluster==g].combined.to_numpy() for g in rng.choice(groups,len(groups))]);bg=rng.choice(back.combined,len(back));boot.append(auc(np.r_[np.ones(len(pos)),np.zeros(len(bg))],np.r_[pos,bg]))
    out([dict(scope='historical8',n_positive=len(pp),publication_clusters=len(groups),n_background=len(back),probability_superiority=auc(np.r_[np.ones(len(pp)),np.zeros(len(back))],np.r_[pp.combined,back.combined]),ci_low=ci(boot)[0],ci_high=ci(boot)[1],method='positive_publication_cluster_and_background_protein_bootstrap',p_value=np.nan,limitation='only4_clusters_one_unresolved_source_not_independent_validation')],dest/'antigen_source_cluster_ci.tsv')
    # Post-hoc primary-source evidence augmentation: no training or formal relabeling.
    evidence=[]
    for gene in list(p.gene)+['ANXA2','ENO1','SLA-I','SLA-II']:
        historical='21119562' if gene in ['CD9','CD46','CD59','PROCR','ANXA2'] else ('22497511' if gene in ['CDH1','CFTR','ENO1'] else ('10.3969/j.issn.1674-7445.2018.05.006' if gene=='FAM234A' else 'UNRESOLVED'))
        evidence.append(dict(gene=gene,historical_source=historical,historical_source_checked='abstract' if historical in ['21119562','22497511'] else 'not_verified_this_round',historical_assay_species='human_IgG' if historical=='22497511' else ('NHP' if historical=='21119562' else 'unresolved'),new_checked_source='38853364_Table1' if gene in ['CD9','CD46','CD59','PROCR','ANXA2'] else None,new_endpoint='human_preformed_IgM_cell_binding' if gene in ['CD9','CD46','CD59','PROCR','ANXA2'] else None,in_scoreable_universe=gene in set(x.gene),evidence_label_modified=False,independent_heldout_gene=False,limitation='CD55_original_positive_support_unresolved' if gene=='CD55' else 'literature_support_not_new_experiment'))
    out(evidence,dest/'antigen_external_evidence.tsv')
    ext=[]
    for gene,n in [('CD46',24),('CD9',8),('CD59',6),('PROCR',2),('ANXA2',1)]:
        z=x[x.gene==gene];ext.append(dict(gene=gene,source='https://onlinelibrary.wiley.com/doi/10.1111/xen.12865',source_table='Table1_column_sum',serum_population='160_swine_veterinarians',human_IgM_reactive_samples=n,n_serum=160,antigenicity_score=z.combined.iloc[0] if len(z) else np.nan,rank=z.recomputed_rank.iloc[0] if len(z) else np.nan,independent_target_holdout=False,comparison='descriptive_only_expression_assay_and_exposure_confounding',source_internal_caveat='3_antigen_case_Fig2_caption_vs_Table1_disagree_on_CD9_or_PROCR'))
    out(ext,dest/'external_human_IgM_2024_mapping.tsv')
    # Structural quality gates are descriptive, do not discard failed examples.
    st=pd.read_csv(dest/'structure_metrics.tsv',sep='\t');q=[]
    for r in st[st.threshold_A==5].to_dict('records'):
        q.append(dict(gene=r['gene'],mode=r['mode'],retained_fraction=r['n_fit']/r['n_aligned'],severe_fit_selection=r['n_fit']/r['n_aligned']<.5,has_interface_pLDDT_below70=r['interface_min_pLDDT']<70,paralog_template=r['gene'] in ['CD86','ULBP1'],relaxation_performed=False,independent_complex_prediction=False,functional_endpoint_validated=False,role=r['role'],status='NOT_VALIDATED_function;geometry_only'))
    out(q,dest/'structure_validation_limits.tsv')
    old=json.loads((BASE/'inputs/legacy_structure_summary.json').read_text());comp=[]
    for gene in ['ICAM1','ULBP1','CD86']:
        k=next(k for k in old if k.startswith(gene+'_'));z=st[(st.gene==gene)&(st.threshold_A==8)&(st['mode']=='raw_all_CA')].iloc[0]
        comp.append(dict(gene=gene,legacy_rmsd=old[k]['rmsd'],legacy_n_pairs=old[k].get('n_pairs',old[k].get('n_sup_pairs')),legacy_contacts8A=old[k]['n_contacts_8A'],new_raw_rmsd=z.rmsd_fit,new_n_aligned=z.n_aligned,new_contact_residues8A=z.pig_residues,new_residue_pairs8A=z.residue_pairs,model_version='AF_v6',comparison='same_current_protocol_not_exact_historical_model_replay'))
    out(comp,dest/'structure_legacy_comparison.tsv')
    geo=pd.read_csv(dest/'structure_geometry_perturbations.tsv',sep='\t');rows=[]
    for (gene,mode),z in geo.groupby(['gene','mode'],sort=True):rows.append(dict(gene=gene,mode=mode,n_draws=len(z),residue_pairs5A_p025=z.residue_pairs.quantile(.025),median=z.residue_pairs.median(),p975=z.residue_pairs.quantile(.975),interval='rigid_geometry_sensitivity_NOT_statistical_CI'))
    out(rows,dest/'structure_geometry_sensitivity.tsv')
    out([dict(task='secreted_Tier2_branch',status='NOT_VALIDATED',missing='locked_raw_candidates_and_topology_sequences'),dict(task='independent_functional_labels',status='NEEDS_DECISION',missing='per_axis_direction_and_endpoint_adjudication'),dict(task='CD55_antigen_positive',status='NOT_VALIDATED',missing='original_specific_antibody_binding_source'),dict(task='FAM234A_primary_fulltext',status='NOT_VALIDATED',missing='accessible_primary_paper_this_round'),dict(task='external_glycoproteomics',status='TODO',missing='site_level_experimental_occupancy_dataset_and_mapping'),dict(task='relaxation_docking_AFMultimer',status='TODO',missing='configured_server_and_modeling_contract'),dict(task='wetlab',status='TODO',missing='real_samples_assays_and_results')],dest/'validation_gap_inventory.tsv')
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--out',type=Path,default=BASE/'results');run(a.parse_args().out)
