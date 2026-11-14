"""Display-only reductions of frozen outputs; no score refit, relabelling or tests."""
from pathlib import Path
import json,hashlib,sys
import numpy as np
import pandas as pd
B=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parents[1];D=B/'analysis'
F=D/'inputs/phase4_snapshot/inputs/prior_snapshot';P2=F/'phase2';P3=F/'phase3'
O=B/'figure_data';O.mkdir(exist_ok=True);meta={}
def read(p):
    meta[str(p.relative_to(D))]=hashlib.sha256(p.read_bytes()).hexdigest()
    return pd.read_csv(p,sep='\t')
def save(x,name):x.to_csv(O/name,sep='\t',index=False,float_format='%.12g')
t=read(P2/'inputs/tier1_frozen.tsv');m=read(P2/'inputs/axis_donor_masks.tsv');common=read(P2/'results/isi_common_cohort.tsv')
bridge=read(P2/'results/isi_phase1_baseline_bridge.tsv');scenarios=read(P2/'results/donor_scenario_scorecards.tsv')
donors=bridge[bridge.model.eq('phase1_mixed3')].donor.tolist();rows=[]
for code,model in [('M','phase1_mixed3'),('S','sequence_only')]:
    for _,x in bridge[bridge.model.eq(model)].iterrows():
        rows.append(dict(scenario=code,donor=x.donor,n_axes=61,T=x['T'],A=x.A,R=x.R,C=x.C,units='weighted_divergence_percentage',source='isi_phase1_baseline_bridge.tsv'))
# Common-denominator ablation: reduce existing per-axis no_divergence scores using the existing masks.
for donor in donors:
    z=common[['pair_id','no_divergence']].merge(m[m.donor.eq(donor)][['pair_id','mask']],on='pair_id',validate='one_to_one')
    assert len(z)==61
    total=z.no_divergence.sum();addressed=(z.no_divergence*z['mask']).sum();residual=(z.no_divergence*(1-z['mask'])).sum()
    rows.append(dict(scenario='N',donor=donor,n_axes=61,T=total,A=addressed,R=residual,C=addressed/total,units='severity_evidence_expression',source='isi_common_cohort.tsv + axis_donor_masks.tsv; exact row reduction'))
for code,scenario,expression,n in [('E','sequence_available','neutral',61),('L','legacy_stored','actual',63),('F','without_divergence','actual',75)]:
    z=scenarios[scenarios.scenario.eq(scenario)&scenarios.expression.eq(expression)&scenarios.severity_power.eq(1)&scenarios.edit_fraction.eq(1)]
    assert len(z)==6
    for _,x in z.iterrows():
        rows.append(dict(scenario=code,donor=x.donor,n_axes=n,T=x['T'],A=x.A,R=x.R,C=x.C,units=x.units,source='donor_scenario_scorecards.tsv'))
scorecards=pd.DataFrame(rows)
assert np.allclose(scorecards['T'],scorecards.A+scorecards.R,atol=1e-8)
assert np.allclose(scorecards.C,scorecards.A/scorecards['T'])
save(scorecards,'supplementary_scenario_scorecards.tsv')
counts=read(F/'results/donor_scorecard.tsv');counts=counts[counts.variant.eq('legacy_stored')].copy();counts['remaining_noncontrol']=counts.non_control_axes-counts.addressed_non_control_axes
save(counts[['donor','resource_axes','non_control_axes','addressed_resource_axes','addressed_non_control_axes','remaining_noncontrol','axis_count_coverage_fraction','addressed_unscored_non_control_axes']],'supplementary_count_coverage.tsv')
cross=read(P2/'results/tier2_axis_crosswalk.tsv');cross=cross[cross.eligible].copy();assert len(cross)==11
save(cross,'supplementary_tier2_crosswalk.tsv')
qc=read(P3/'results/sequence_domain_qc.tsv');qc=qc[qc.region.eq('full_length')].copy();assert len(qc)==20
save(qc,'supplementary_full_length_qc.tsv')
structure=read(P2/'results/structure_metrics.tsv');contacts=structure[structure.threshold_A.eq(8)].copy();assert len(contacts)==10
save(contacts,'supplementary_8A_contacts.tsv')
org=read(P3/'results/organ_annotation_membership.tsv');cards=read(P3/'results/organ_annotation_scorecards.tsv')
rev=next(d for d in donors if d.startswith('Revivicor'));clon=next(d for d in donors if d.startswith('ClonOrgan-8'))
pivot=m.pivot(index='pair_id',columns='donor',values='mask');diff=pivot.index[pivot[rev]!=pivot[clon]]
expl=t[t.pair_id.isin(diff)][['pair_id','pig_gene','organs','func_sev','identity']].merge(org[org.organ.eq('liver')][['pair_id','included_annotation_subset']],on='pair_id')
assert set(expl.pair_id)=={'XP009','XP052'}
assert not bool(expl.loc[expl.pair_id.eq('XP052'),'included_annotation_subset'].iloc[0])
assert expl.loc[expl.pair_id.eq('XP009'),'func_sev'].iloc[0]==0
save(expl,'liver_mask_explanation.tsv')
facts=dict(total_axes=len(t),unique_pair_ids=t.pair_id.nunique(),evidence=t.evidence.value_counts().to_dict(),sequence_available=int(t.identity.notna().sum()),sequence_noncontrol=int((t.identity.notna()&t.evidence.ne('control')).sum()),stored_numeric=int(t.ISI_v2.notna().sum()),stored_noncontrol=int((t.ISI_v2.notna()&t.evidence.ne('control')).sum()),stored_control=int((t.ISI_v2.notna()&t.evidence.eq('control')).sum()),missing_noncontrol=int((t.identity.isna()&t.evidence.ne('control')).sum()),missing_noncontrol_stored_numeric=int((t.identity.isna()&t.evidence.ne('control')&t.ISI_v2.notna()).sum()),EPO_axes=t.loc[t.gene_key.eq('EPO'),'pair_id'].tolist(),CD99_core_axes=t.loc[t.gene_key.eq('CD99'),'pair_id'].tolist(),XP045=t.loc[t.pair_id.eq('XP045'),'pig_gene'].iloc[0],same_cohort_no_divergence_total=float(common.no_divergence.sum()),THBD_interface_pLDDT_at_5A=float(structure.loc[structure.gene.eq('THBD')&structure.threshold_A.eq(5)&structure['mode'].eq('filtered_CA'),'interface_mean_pLDDT'].iloc[0]))
assert facts['stored_numeric']==63 and facts['stored_noncontrol']==51 and facts['stored_control']==12
assert facts['EPO_axes']==['XP049'] and not facts['CD99_core_axes']
(O/'checked_review_claims.json').write_text(json.dumps(facts,ensure_ascii=False,indent=2)+'\n')
(O/'INPUT_PROVENANCE.json').write_text(json.dumps(meta,indent=2)+'\n')
print('Extracted five supplementary displays and verified disputed review claims; no hypothesis tests rerun')
