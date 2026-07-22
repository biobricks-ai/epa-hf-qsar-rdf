#!/usr/bin/env python3
import argparse,gzip,hashlib,json,re
from pathlib import Path
from urllib.parse import quote
import pyarrow.parquet as pq
B="https://biobricks.ai/epa-hf-qsar-rdf/";RDF="http://www.w3.org/1999/02/22-rdf-syntax-ns#";RDFS="http://www.w3.org/2000/01/rdf-schema#";DCT="http://purl.org/dc/terms/";PROV="http://www.w3.org/ns/prov#";SC="https://schema.org/";SKOS="http://www.w3.org/2004/02/skos/core#";OWL="http://www.w3.org/2002/07/owl#";IAO="http://purl.obolibrary.org/obo/IAO_0000027";DATUM="http://purl.obolibrary.org/obo/IAO_0000109";SRC="https://catalog.data.gov/dataset/toxicity-assessment-for-epas-hydraulic-fracturing-study";CAS_RE=re.compile(r"^\d{2,7}-\d{2}-\d$")
def v(x):s="" if x is None else str(x).strip();return s if s and s.lower() not in {"nan","none","null"} else None
def tr(s,p,o,i=False):return f"<{s}> <{p}> "+(f"<{o}>" if i else json.dumps(str(o),ensure_ascii=False))+" .\n"
def casok(x):
 if not x or not CAS_RE.fullmatch(x):return False
 d=x.replace('-','');return sum(int(c)*i for i,c in enumerate(reversed(d[:-1]),1))%10==int(d[-1])
def kind(table,col):
 x=(table+' '+col).lower()
 if table.startswith('qsar_topkat') or table.startswith('mcda_') and ('score' in x or 'rank' in x) or table=='mcda_episuite_properties' and not 'measured' in x:return 'model-derived prediction or decision-analysis value'
 if any(k in x for k in ('iris','pprtv','atsdr','oehha','hhbp','cicad','measured','rfv','osf')):return 'experimental or authoritative reference value'
 if any(k in x for k in ('disclosures','ingredients','percent')):return 'source occurrence summary'
 return None
def unit(col):
 c=col.lower()
 for s,u in [('mg/kg-day','mg/kg-day'),('per mg/kg-day','(mg/kg-day)^-1'),('mg/kg bw','mg/kg bw'),('g/kg bw','g/kg bw'),('mg/l','mg/L'),('atm-m3/mole','atm-m3/mol'),('hours','h'),('days','d'),('log ratio','log ratio'),('log l/kg','log L/kg')]:
  if s in c:return u
 return 'score-or-count'
def convert(src,out):
 out.mkdir(parents=True,exist_ok=True);z={"tables":{},"source_rows":0,"represented_rows":0,"source_nonempty_cells":0,"mapped_nonempty_cells":0,"semantic_datums":0,"predicted_datums":0,"reference_datums":0}
 with gzip.open(out/'epa-hf-qsar.nt.gz','wt') as sem,gzip.open(out/'epa-hf-qsar.source.nt.gz','wt') as ev:
  for p in sorted(src.glob('*.parquet')):
   pf=pq.ParquetFile(p);t={"source_rows":pf.metadata.num_rows,"represented_rows":0,"source_fields":len(pf.schema.names),"fields_dispositioned":len(pf.schema.names)};z['source_rows']+=pf.metadata.num_rows
   for n,row in enumerate(x for b in pf.iter_batches() for x in b.to_pylist()):
    cas=v(row.get('CASRN'));h=hashlib.sha256((p.stem+str(n)+str(cas)).encode()).hexdigest()[:20];r=B+f"record/{p.stem}/{h}";c=B+(f"substance/cas/{cas}" if casok(cas) else f"substance/record/{h}");ls=[tr(r,RDF+'type',IAO,True),tr(r,RDFS+'label',f"{p.stem} row for {v(row.get('Chemical Name')) or cas}"),tr(r,PROV+'wasDerivedFrom',SRC,True),tr(r,SC+'about',c,True)]
    if v(row.get('Chemical Name')):ls.append(tr(c,RDFS+'label',v(row['Chemical Name'])))
    if casok(cas):ls += [tr(c,SKOS+'exactMatch','https://identifiers.org/cas/'+cas,True),tr(c,OWL+'sameAs','https://identifiers.org/cas/'+cas,True)]
    for col,x in row.items():
     k=kind(p.stem,col)
     if k and isinstance(x,(int,float)) and v(x):
      m=r+'/datum/'+quote(col,safe='');ls += [tr(m,RDF+'type',DATUM,True),tr(m,RDFS+'label',col),tr(m,DCT+'type',k),tr(m,SC+'about',c,True),tr(m,SC+'value',v(x)),tr(m,SC+'unitCode',unit(col)),tr(m,PROV+'wasDerivedFrom',r,True)];z['semantic_datums']+=1;z['predicted_datums']+=int(k.startswith('model'));z['reference_datums']+=int(k.startswith('experimental'))
    evidence=list(ls);non=0
    for col,x in row.items():
     if v(x):non+=1;evidence.append(tr(r,B+'source-column/'+p.stem+'/'+quote(col,safe=''),v(x)))
    sem.write(''.join(ls));ev.write(''.join(evidence));t['represented_rows']+=1;z['represented_rows']+=1;z['source_nonempty_cells']+=non;z['mapped_nonempty_cells']+=non
   t['row_coverage']=t['field_disposition_coverage']=1;z['tables'][p.stem]=t
 z['row_coverage']=z['represented_rows']/z['source_rows'];z['mapped_cell_coverage']=z['mapped_nonempty_cells']/z['source_nonempty_cells'];(out/'coverage.json').write_text(json.dumps(z,indent=2)+'\n');return z
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('output',type=Path);a=p.parse_args();print(json.dumps(convert(a.source,a.output),indent=2))
