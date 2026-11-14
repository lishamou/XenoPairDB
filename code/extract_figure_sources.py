"""Non-graphical extraction of exact frozen rows for the R figure script."""
from pathlib import Path
import json,hashlib,sys
import pandas as pd,numpy as np
B=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parents[1];D=B/'analysis';F=D/'inputs/phase4_snapshot/inputs/prior_snapshot/phase2';O=B/'figures';O.mkdir(exist_ok=True);meta={}
def read(p):meta[str(p.relative_to(B))]=hashlib.sha256(p.read_bytes()).hexdigest();return pd.read_csv(p,sep='\t')
def save(x,n):x.to_csv(O/f'figure{n}_source.tsv',sep='\t',index=False,float_format='%.12g')
t=read(F/'inputs/tier1_frozen.tsv');t['sequence_available']=pd.to_numeric(t.identity,errors='coerce').notna();save(t[['pair_id','system','evidence','sequence_available']],1)
d=read(F/'results/isi_phase1_baseline_bridge.tsv');d=d[d.model.eq('phase1_mixed3')].copy();d['donor_short']=d.donor.str.split(' (',regex=False).str[0];assert np.allclose(d['T'],d.A+d.R);save(d[['donor_short','T','A','R','C']],2)
a=read(F/'results/antigen_priority_all.tsv');c=read(D/'results/cohort_M3_observed_control.tsv');p=read(D/'results/cohort_M0_observed_control.tsv');c['primary']=c.acc.isin(p.acc);save(pd.concat([a[['acc','gene','combined','historical_positive']].assign(panel='a'),c[['acc','gene','antigenicity_v2','log2fc','primary']].assign(panel='b')],ignore_index=True),3)
s=read(F/'results/structure_metrics.tsv');save(s[s.threshold_A.eq(5)],4)
(O/'INPUT_PROVENANCE.json').write_text(json.dumps(meta,indent=2)+'\n')
print('Four figure source tables extracted without numerical refitting')
