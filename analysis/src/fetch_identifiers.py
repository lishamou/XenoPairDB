"""Cache authoritative UniProt records without uploading scores or outcomes."""
from pathlib import Path
import argparse,json,hashlib,urllib.request,urllib.parse,datetime,concurrent.futures
import pandas as pd
p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
f=pd.read_csv(a.input/'frozen_antigen.tsv',sep='\t');d=pd.read_csv(a.input/'primary_mapping.tsv',sep='\t')
ids=sorted(set(f.acc)|set(';'.join(d.source_accession).split(';')))
jobs=[ids[i:i+100] for i in range(0,len(ids),100)]
def fetch(job):
 i,accessions=job;name=f'current_{i:03d}';file=a.output/(name+'.json');meta=a.output/(name+'.meta.json')
 if file.exists() and meta.exists():return json.loads(meta.read_text())
 query='('+ ' OR '.join('accession:'+v for v in accessions)+')'
 url='https://rest.uniprot.org/uniprotkb/search?'+urllib.parse.urlencode(dict(query=query,format='json',size=500))
 record=dict(requested=accessions,url=url,retrieved_at=datetime.datetime.now(datetime.timezone.utc).isoformat())
 try:
  with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'XenoPairDB-identifier-audit/1.0'}),timeout=60) as r:
   raw=r.read();data=json.loads(raw);record.update(status=r.status,release=r.headers.get('X-UniProt-Release'),release_date=r.headers.get('X-UniProt-Release-Date'),total=r.headers.get('X-Total-Results'),link=r.headers.get('Link'),received=len(data['results']),sha256=hashlib.sha256(raw).hexdigest())
   assert not record['link'],'Unexpected pagination; do not silently truncate'
   file.write_bytes(raw)
 except Exception as e:record.update(status='ERROR',error=str(e),response=e.read().decode() if isinstance(e,urllib.error.HTTPError) else None)
 meta.write_text(json.dumps(record,indent=2)+'\n');return record
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:
 records=list(pool.map(fetch,enumerate(jobs)))
(a.output/'FETCH_REPORT.json').write_text(json.dumps(dict(requested_unique=len(ids),batches=len(jobs),records=records),indent=2)+'\n')
print('Requested',len(ids),'accessions in',len(jobs),'batches; errors',sum(r['status']=='ERROR' for r in records))
assert all(r['status']!='ERROR' for r in records),'Unresolved fetch failures'
