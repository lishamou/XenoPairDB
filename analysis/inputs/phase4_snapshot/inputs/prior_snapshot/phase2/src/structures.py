"""Template-superposition geometry audit. No affinity or function prediction."""
from pathlib import Path
import json,copy,argparse
import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import PDBParser,PDBIO,Structure,Model,Chain,Residue
from Bio.SeqUtils import seq1
from analyze import BASE,out,SEED
D=BASE/'inputs/structures'
CASES=[dict(gene='ICAM1',pdb='1MQ8',template='A',partners=['B','D'],role='original_case',caveat='ITGAL_I_domain_only_not_full_heterodimer'),dict(gene='ULBP1',pdb='1KCG',template='C',partners=['AB'],role='original_case',caveat='human_ULBP3_paralog_template'),dict(gene='CD86',pdb='1I8L',template='A',partners=['C','D'],role='original_case',caveat='human_CD80_paralog_template'),dict(gene='CD274',pdb='4ZQK',template='A',partners=['B'],role='preserved_function_control',caveat='PMID17489864_PD1_dependent_inhibition'),dict(gene='THBD',pdb='1DX5',template='I',partners=['MA','NB','OC','PD'],role='impaired_cofactor_control',caveat='PMID18444940_thrombin_binding_preserved_PC_activation_impaired_template_endpoint_mismatch')]
def residues(chain):return [r for r in chain if r.id[0]==' ' and 'CA' in r and seq1(r.resname,custom_map={'MSE':'M'})!='X']
def sequence(rs):return ''.join(seq1(r.resname,custom_map={'MSE':'M'}) for r in rs)
def mm(a,b):
    return np.einsum("ij,jk->ik",a,b) if np.asarray(a).ndim==2 else np.einsum("i,ij->j",a,b)

def fit(p,q):
    p=np.asarray(p,float);q=np.asarray(q,float)
    if len(p)<3 or p.shape!=q.shape:raise ValueError('insufficient matching coordinates')
    pc=p.mean(axis=0);qc=q.mean(axis=0);u,_,vt=np.linalg.svd(mm((p-pc).T,(q-qc)));r=mm(u,vt)
    if np.linalg.det(r)<0:u[:,-1]*=-1;r=mm(u,vt)
    return r,qc-mm(pc,r)

def heavy(rs):
    atoms=[];labels=[]
    for r in rs:
        for a in r:
            if a.element not in ['H','D'] and (a.occupancy is None or a.occupancy>0):atoms.append(a);labels.append(f'{r.parent.id}:{r.id[1]}{r.id[2].strip()}:{r.resname}')
    return np.array([a.coord for a in atoms],float),labels

def contacts(a,b,la,lb,cut):
    if not len(a) or not len(b):return [],dict(atom_pairs=0,residue_pairs=0,pig_residues=0,partner_residues=0,min_distance=np.nan)
    matrix=cKDTree(a).sparse_distance_matrix(cKDTree(b),cut,output_type='coo_matrix')
    pairs={};count=0
    for i,j,d in zip(matrix.row,matrix.col,matrix.data):
        if d>=cut:continue
        count+=1;key=(la[i],lb[j]);pairs[key]=min(pairs.get(key,float('inf')),float(d))
    rows=[dict(pig_residue=k[0],partner_residue=k[1],min_heavy_distance=v) for k,v in sorted(pairs.items())]
    return rows,dict(atom_pairs=count,residue_pairs=len(pairs),pig_residues=len(set(k[0] for k in pairs)),partner_residues=len(set(k[1] for k in pairs)),min_distance=min(pairs.values()) if pairs else np.nan)

def save_model(rs,partners,r,t,path):
    structure=Structure.Structure('template_model');model=Model.Model(0);structure.add(model);chain=Chain.Chain('X');model.add(chain)
    for res in rs:
        nr=Residue.Residue(res.id,res.resname,res.segid)
        for atom in res:nr.add(atom.copy())
        for a in nr:a.set_coord(mm(a.coord,r)+t)
        chain.add(nr)
    for idx,prs in enumerate(partners):
        ch=Chain.Chain(chr(65+idx));model.add(ch)
        for res in prs:
            nr=Residue.Residue(res.id,res.resname,res.segid)
            for atom in res:nr.add(atom.copy())
            ch.add(nr)
    path.parent.mkdir(exist_ok=True,parents=True);io=PDBIO();io.set_structure(structure);io.save(str(path))

def run(dest):
    parser=PDBParser(QUIET=True);aligner=Align.PairwiseAligner(mode='local');aligner.substitution_matrix=substitution_matrices.load('BLOSUM62');aligner.open_gap_score=-10;aligner.extend_gap_score=-.5
    metrics=[];mapping=[];contactrows=[];selection=[];perturb=[];provenance=[];rng=np.random.default_rng(SEED)
    for case in CASES:
        gene,pdb=case['gene'],case['pdb'];h=parser.get_structure(pdb,D/f'{pdb}.pdb')[0];pg=parser.get_structure(gene,D/f'{gene}_pig.pdb')[0];pr=residues(next(pg.get_chains()));hr=residues(h[case['template']]);ps,hs=sequence(pr),sequence(hr)
        uniprot=json.loads((D/f'{gene}_uniprot.json').read_text());af=json.loads((D/'selection.json').read_text());af=next(j['selected_model'] for j in af if j['gene']==gene)
        names=[g.get('geneName',{}).get('value') for g in uniprot.get('genes',[])]
        if uniprot['organism']['taxonId']!=9823 or gene not in names or ps!=uniprot['sequence']['value'] or af['taxId']!=9823:raise ValueError(f'{gene}:species_gene_sequence mismatch')
        provenance.append(dict(gene=gene,uniprot=uniprot['primaryAccession'],taxon=9823,sequence_length=len(ps),model=af['entryId'],version=af['latestVersion'],sequence_exact=True,pdb_template=pdb,template_species=9606,role=case['role'],caveat=case['caveat']))
        al=aligner.align(hs,ps)[0];pairs=[]
        for hb,pb in zip(al.aligned[0],al.aligned[1]):pairs.extend(zip(range(*hb),range(*pb)))
        H=np.array([hr[i]['CA'].coord for i,j in pairs],float);P=np.array([pr[j]['CA'].coord for i,j in pairs],float);r0,t0=fit(P,H)
        keep=np.ones(len(P),bool);rf,tf=r0,t0
        for iteration in range(10):
            dist=np.linalg.norm(mm(P,rf)+tf-H,axis=1);new=keep&(dist<=2.5)
            if new.sum()<10 or np.array_equal(new,keep):break
            keep=new;rf,tf=fit(P[keep],H[keep])
        for k,(i,j) in enumerate(pairs):mapping.append(dict(gene=gene,pdb=pdb,human_chain=case['template'],human_residue=hr[i].id[1],human_aa=hs[i],pig_residue=pr[j].id[1],pig_aa=ps[j],identical=hs[i]==ps[j],pLDDT=pr[j]['CA'].bfactor,kept_for_filtered_fit=bool(keep[k]),raw_distance=np.linalg.norm(mm(P[k],r0)+t0-H[k]),filtered_distance=np.linalg.norm(mm(P[k],rf)+tf-H[k])))
        domain=pr[min(j for i,j in pairs):max(j for i,j in pairs)+1];pc,pl=heavy(domain);hc,hl=heavy(hr)
        # Select physical partner assembly using HUMAN template contacts, never pig score.
        opts=[]
        for candidate in case['partners']:
            partner_lists=[residues(h[ch]) for ch in candidate];partners=sum(partner_lists,[]);bc,bl=heavy(partners);cr,cm=contacts(hc,bc,hl,bl,5)
            opts.append((cm['residue_pairs'],candidate,partner_lists,bc,bl));selection.append(dict(gene=gene,pdb=pdb,template_chain=case['template'],candidate_partner_chains=candidate,human_template_residue_pairs_5A=cm['residue_pairs']))
        _,chosen,partner_lists,bc,bl=max(opts,key=lambda o:(o[0],o[1]))
        for mode,r,t in [('raw_all_CA',r0,t0),('filtered_CA',rf,tf)]:
            xp=mm(pc,r)+t;dist=np.linalg.norm(mm(P,r)+t-H,axis=1);save_model(domain,partner_lists,r,t,dest/'models'/f'{gene}_{mode}.pdb')
            for cutoff in [3,5,8]:
                cr,cm=contacts(xp,bc,pl,bl,cutoff);baseline,hm=contacts(hc,bc,hl,bl,cutoff)
                interface=set(c['pig_residue'] for c in cr);conf=[res['CA'].bfactor for res in domain if f'{res.parent.id}:{res.id[1]}{res.id[2].strip()}:{res.resname}' in interface]
                row=dict(gene=gene,pdb=pdb,role=case['role'],template_chain=case['template'],partner_chains=chosen,mode=mode,threshold_A=cutoff,n_aligned=len(pairs),n_fit=int(keep.sum()) if mode=='filtered_CA' else len(pairs),n_removed=int((~keep).sum()) if mode=='filtered_CA' else 0,rmsd_fit=np.sqrt(np.mean(dist[keep if mode=='filtered_CA' else np.ones(len(keep),bool)]**2)),rmsd_all_aligned=np.sqrt(np.mean(dist**2)),human_resolved_coverage=len(pairs)/len(hr),pig_domain_start=domain[0].id[1],pig_domain_end=domain[-1].id[1],aligned_identity=np.mean([hs[i]==ps[j] for i,j in pairs]),domain_mean_pLDDT=np.mean([z['CA'].bfactor for z in domain]),interface_mean_pLDDT=np.mean(conf) if conf else np.nan,interface_min_pLDDT=min(conf) if conf else np.nan,interface_fraction_pLDDT_lt70=np.mean(np.array(conf)<70) if conf else np.nan,template_residue_pairs=hm['residue_pairs'],template_pig_side_residues=hm['pig_residues'],**cm,caveat=case['caveat'],functional_validation='NOT_ESTABLISHED')
                metrics.append(row)
                for c in cr:contactrows.append(dict(gene=gene,mode=mode,threshold_A=cutoff,**c))
            # Geometry perturbations are a sensitivity distribution, not binding null.
            center=xp.mean(axis=0)
            for it in range(100):
                rr=Rotation.from_rotvec(rng.normal(0,np.deg2rad(5),3)).as_matrix();tt=rng.normal(0,2,3);xq=mm(xp-center,rr)+center+tt
                _,cm=contacts(xq,bc,pl,bl,5);perturb.append(dict(gene=gene,mode=mode,iteration=it,**cm,distribution='rigid_geometry_sensitivity_NOT_functional_null'))
        print(gene,'aligned',len(pairs),'kept',keep.sum(),'chosen',chosen,flush=True)
    out(metrics,dest/'structure_metrics.tsv');out(mapping,dest/'structure_residue_mapping.tsv');out(contactrows,dest/'structure_contacts.tsv');out(selection,dest/'structure_chain_selection.tsv');out(provenance,dest/'structure_input_validation.tsv');out(perturb,dest/'structure_geometry_perturbations.tsv')
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--out',type=Path,default=BASE/'results');args=a.parse_args();run(args.out)
