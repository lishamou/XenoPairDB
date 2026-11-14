"""Download public coordinates only; never uploads project data."""
from pathlib import Path
import urllib.request,json,hashlib,datetime,csv,concurrent.futures
BASE=Path(__file__).resolve().parents[1];D=BASE/'inputs/structures';D.mkdir(exist_ok=True)
CASES=[('ICAM1','F1S3J9','1MQ8'),('ULBP1','I6L5V9','1KCG'),('CD86','K7GMU1','1I8L'),('CD274','A0A8D2BPF9','4ZQK'),('THBD','A0A8W4FEQ7','1DX5')]
def get(url,path):
    if not path.exists():
        req=urllib.request.Request(url,headers={'User-Agent':'XenoPairDB-revision-local-research/1.0'})
        with urllib.request.urlopen(req,timeout=35) as h:data=h.read()
        path.write_bytes(data)
    return dict(file=str(path.relative_to(BASE)),url=url,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),bytes=path.stat().st_size)
def run(case):
    gene,acc,pdb=case;logs=[];state={'gene':gene,'accession':acc,'pdb':pdb}
    tasks=[('af_metadata',f'https://alphafold.ebi.ac.uk/api/prediction/{acc}',D/f'{gene}_af_metadata.json'),('uniprot',f'https://rest.uniprot.org/uniprotkb/{acc}.json',D/f'{gene}_uniprot.json'),('pdb',f'https://files.rcsb.org/download/{pdb}.pdb',D/f'{pdb}.pdb'),('rcsb_entry',f'https://data.rcsb.org/rest/v1/core/entry/{pdb}',D/f'{pdb}_entry.json')]
    for typ,url,path in tasks:
        try:logs.append(get(url,path))
        except Exception as e:logs.append(dict(file=str(path.relative_to(BASE)),url=url,error=str(e)))
    try:
        a=json.loads((D/f'{gene}_af_metadata.json').read_text());state['af_records']=len(a)
        exact=[r for r in a if r.get('uniprotAccession')==acc and r.get('uniprotStart',1)==1]
        if not exact:raise ValueError('no exact full-start model')
        model=sorted(exact,key=lambda r:r.get('uniprotEnd',0),reverse=True)[0]
        state['selected_model']=model;logs.append(get(model['pdbUrl'],D/f'{gene}_pig.pdb'))
    except Exception as e:state['model_error']=str(e)
    try:
        entry=json.loads((D/f'{pdb}_entry.json').read_text())
        for entity in entry['rcsb_entry_container_identifiers']['polymer_entity_ids']:
            logs.append(get(f'https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb}/{entity}',D/f'{pdb}_entity{entity}.json'))
    except Exception as e:state['entity_error']=str(e)
    print(gene,'model_error',state.get('model_error'),'errors',sum('error' in l for l in logs),flush=True)
    return state,logs
if __name__=='__main__':
    results=list(concurrent.futures.ThreadPoolExecutor(max_workers=5).map(run,CASES))
    (D/'selection.json').write_text(json.dumps([r[0] for r in results],indent=2))
    (BASE/'logs/download_manifest.json').write_text(json.dumps(dict(retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),files=[l for r in results for l in r[1]]),indent=2))
