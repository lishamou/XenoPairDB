"""Targeted direct-entry lookup for unresolved source accessions; no biological guesses."""
from pathlib import Path
import json,urllib.request,datetime,hashlib,concurrent.futures
B=Path(__file__).resolve().parents[1];O=B/'inputs/uniprot';ids=json.loads((B/'results/mapping_summary.json').read_text())['missing_source']
def fetch(acc):
 file=O/('supplement_'+acc+'.json');meta=O/('supplement_'+acc+'.meta.json');url='https://rest.uniprot.org/uniprotkb/'+acc+'.json'
 if file.exists() and meta.exists():return
 r=dict(requested=acc,url=url,retrieved_at=datetime.datetime.now(datetime.timezone.utc).isoformat())
 try:
  with urllib.request.urlopen(url,timeout=30) as response:
   raw=response.read();json.loads(raw);file.write_bytes(raw);r.update(status=response.status,final_url=response.url,release=response.headers.get('X-UniProt-Release'),sha256=hashlib.sha256(raw).hexdigest())
 except Exception as e:r.update(status='ERROR',error=str(e),response=e.read().decode() if isinstance(e,urllib.error.HTTPError) else None)
 meta.write_text(json.dumps(r,indent=2)+'\n')
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:list(pool.map(fetch,ids))
print('Targeted direct entries requested',len(ids))
