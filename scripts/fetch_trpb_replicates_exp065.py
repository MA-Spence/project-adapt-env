#!/usr/bin/env python3
"""Retrieve only the declared public TrpB ZIP members with CRC/SHA provenance.
External code is retained for inspection, never executed by this script.
"""
import io,json,requests,zipfile,hashlib
from pathlib import Path
class RemoteZip(io.RawIOBase):
 def __init__(self,url):
  self.url=url;self.pos=0;self.session=requests.Session()
  r=self.session.get(url,headers={'Range':'bytes=0-0'},stream=True,timeout=45)
  print('range response',r.status_code,dict((k,v) for k,v in r.headers.items() if k.lower() in ['content-range','content-length','accept-ranges']),flush=True)
  if r.status_code!=206: r.close();raise RuntimeError('Server does not provide byte ranges')
  self.size=int(r.headers['Content-Range'].split('/')[-1]);r.close()
 def seekable(self):return True
 def readable(self):return True
 def seek(self,n,w=0):
  self.pos=n if w==0 else self.pos+n if w==1 else self.size+n
  return self.pos
 def tell(self):return self.pos
 def read(self,n=-1):
  if n<0:n=self.size-self.pos
  n=min(n,self.size-self.pos)
  if n<=0:return b''
  start=self.pos;end=start+n-1
  # A range-specific query avoids cache servers serving a previous range.
  r=self.session.get(self.url+f'&range={start}-{end}',headers={'Range':f'bytes={start}-{end}'},timeout=90)
  r.raise_for_status()
  assert r.status_code==206 and r.headers['Content-Range'].startswith(f'bytes {start}-{end}/'),(r.status_code,r.headers)
  data=r.content;assert len(data)==n;self.pos+=len(data);return data

if __name__ == "__main__":
 root=Path(__file__).resolve().parents[1]/'data/processed/trpb-replicate-observation-exp065/source'
 root.mkdir(parents=True,exist_ok=True)
 names={'data':['data/figure_data/4-site_merged_replicates/20230827/four-site_summary_AA_data.csv','data/figure_data/4-site_merged_replicates/20230827/four-site_simplified_AA_data.csv'],
 'code':['code/scripts/growth_analysis/growth_analysis/processing_pipeline.py','code/scripts/growth_analysis/growth_analysis/globals.py','code/scripts/growth_analysis/growth_analysis/functions.py','code/notebooks/foursite_seq_processing/4-site_import_combine_data.ipynb','code/notebooks/additional_supplementary/4-site_stop_codon_threshold.ipynb']}
 records=[]
 for archive,paths in names.items():
  url=f'https://data.caltech.edu/records/h5rah-5z170/files/{archive}.zip?download=1';z=zipfile.ZipFile(RemoteZip(url))
  for member in paths:
   info=z.getinfo(member);file=root/Path(member).name
   data=z.read(member);assert len(data)==info.file_size
   file.write_bytes(data)
   records.append({'url':url,'member':member,'file':str(file.relative_to(root.parent)),'size':len(data),'crc32':info.CRC,'sha256':hashlib.sha256(data).hexdigest()})
   print('Retrieved',file.name,len(data),flush=True)
 (root.parent/'source_manifest.json').write_text(json.dumps({'doi':'10.22002/h5rah-5z170','retrieval':'HTTP ranges; ZIP member CRC verified; complete archive MD5 not checked because only selected members retrieved','files':records},indent=2)+'\n')
