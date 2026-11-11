"""Snapshot ChEMBL single-strain, 24 h microdilution MIC records."""
from pathlib import Path
import json, urllib.request, urllib.parse, time, re, hashlib, gzip
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
API = 'https://www.ebi.ac.uk/chembl/api/data/'

def get(url):
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=120) as f: return json.load(f)
        except Exception:
            if attempt == 3: raise
            time.sleep(2**attempt)

def pages(endpoint, params, key):
    url = API + endpoint + '.json?' + urllib.parse.urlencode(dict(params,limit=1000))
    rows, urls = [], []
    while url:
        data = get(url)
        urls.append(url)
        rows.extend(data[key])
        nxt = data['page_meta']['next']
        url = urllib.parse.urljoin('https://www.ebi.ac.uk', nxt) if nxt else None
    return rows, urls

def assay_reason(a):
    desc = a['description'].lower()
    if a.get('assay_organism') != 'Staphylococcus aureus': return 'organism'
    if re.sub(r'\s','',a.get('assay_strain') or '').upper() != 'ATCC29213': return 'strain'
    if a.get('target_chembl_id') != 'CHEMBL352' or a.get('assay_type') != 'F': return 'not_whole_cell_functional'
    if a.get('src_id') != 1: return 'not_literature'
    if 'microdilution' not in desc: return 'method_not_explicit_microdilution'
    if not re.search(r'(?:after|for)\s+24\s*(?:h|hr|hour)',desc): return 'duration_not_explicit_24h'
    if re.search(r'resistant|presence|serum|washout|followed|biofilm|synerg|combination|mutant|irradi|light|liposom|nanopart|encapsul|persister',desc): return 'modified_or_special_conditions'
    if re.search(r'bactericidal|log reduction|90%|alamar|staining',desc): return 'different_readout_or_killing_endpoint'
    if re.search(r'\bat\s+[<>]?\s*\d+(?:\.\d+)?\s*(?:ug|um|mg)',desc): return 'single_concentration_description'
    return ''

def save_gzip(path, data):
    path.write_bytes(gzip.compress(json.dumps(data,ensure_ascii=True).encode(),mtime=0))

def main():
    raw=ROOT/'data/raw'; raw.mkdir(parents=True,exist_ok=True)
    if (raw/'activities.json.gz').exists(): raise FileExistsError('Snapshot exists; do not silently update it')
    assays,urls=pages('assay',{'assay_organism':'Staphylococcus aureus','description__icontains':'29213','order_by':'assay_chembl_id'},'assays')
    audit=[dict(a,selection_reason=assay_reason(a)) for a in assays]
    selected=[a for a in audit if not a['selection_reason']]
    save_gzip(raw/'assays.json.gz',audit)
    records=[]
    ids=[a['assay_chembl_id'] for a in selected]
    for i in range(0,len(ids),20):
        rows,used=pages('activity',{'assay_chembl_id__in':','.join(ids[i:i+20]),'standard_type':'MIC','order_by':'activity_id'},'activities')
        if any(r['assay_chembl_id'] not in ids[i:i+20] for r in rows): raise ValueError('API filter not applied')
        records+=rows; urls+=used
        print(f'Assays {min(i+20,len(ids))}/{len(ids)}; records {len(records)}',flush=True)
    save_gzip(raw/'activities.json.gz',records)
    status=get(API+'status.json')
    meta={'retrieved_utc':datetime.now(timezone.utc).isoformat(),'chembl_status':status,'urls':urls,
          'queried_assays':len(assays),'selected_assays':len(selected),'activity_records':len(records),
          'sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in raw.glob('*.gz')}}
    (raw/'provenance.json').write_text(json.dumps(meta,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in meta.items() if k not in ['urls','chembl_status','sha256']},indent=2))

if __name__=='__main__': main()
