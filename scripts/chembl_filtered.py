#!/usr/bin/env python3
"""ChEMBL bioactivity re-query with the filtering comment 7 requires.

Three distinct-compound counts per target, over the human single-protein target:
  n_current          IC50/Ki/Kd/EC50 <= 10 uM, no other filter (reproduces the
                     published number, for continuity).
  n_pchembl          pChEMBL >= 5. pChEMBL is ChEMBL's curated, standardised
                     activity: it is assigned only to '=' relations, harmonises
                     the assay types, and excludes records with data-validity
                     flags. This is the primary bioactivity-evidence measure.
  n_pchembl_binding  pChEMBL >= 5 AND assay_type = B (direct binding assays only):
                     the high-confidence sensitivity the reviewer asks for.
Compounds counted as distinct molecule_chembl_id. Absence of a human ChEMBL
target entity is recorded as a coverage gap, not biological intractability.
"""
import csv, json, os, time, urllib.parse, urllib.request
BASE="https://www.ebi.ac.uk/chembl/api/data"
S4={r['target_id']:r for r in csv.DictReader(open(os.path.expanduser('~/r2/data/SupplementaryTable_S4_chembl_ligand_evidence.csv')))}
def get(path,**p):
    u=f"{BASE}/{path}.json?"+urllib.parse.urlencode(p)
    for _ in range(4):
        try: return json.load(urllib.request.urlopen(u,timeout=60))
        except Exception: time.sleep(3)
    return None
def distinct(**p):
    """distinct molecule_chembl_id over all pages (only= keeps pages light)."""
    s=set(); off=0
    while off<5000:
        pg=get("activity", only="molecule_chembl_id", limit=1000, offset=off, **p)
        acts=(pg or {}).get("activities",[])
        for a in acts:
            if a.get("molecule_chembl_id"): s.add(a["molecule_chembl_id"])
        if len(acts)<1000: break
        off+=1000
    return len(s)
rows=[]
for tid in sorted(S4):
    r=S4[tid]; cid=r['chembl_id']; gene=r['gene']
    if not cid:
        rows.append(dict(target_id=tid,gene=gene,uniprot=r['uniprot'],chembl_id='',
            status='NO_ENTITY_coverage_gap',n_current=0,n_pchembl=0,n_pchembl_binding=0))
        print(f"  {tid} {gene:9} no human ChEMBL entity  -> coverage gap"); continue
    n_cur=distinct(target_chembl_id=cid, standard_type__in="IC50,Ki,Kd,EC50",
                   standard_units="nM", standard_value__lte=10000)
    n_pc =distinct(target_chembl_id=cid, pchembl_value__gte=5)
    n_pcb=distinct(target_chembl_id=cid, pchembl_value__gte=5, assay_type="B")
    status='OK' if n_pc>0 else 'ENTITY_zero_potent_compounds'
    rows.append(dict(target_id=tid,gene=gene,uniprot=r['uniprot'],chembl_id=cid,status=status,
                     n_current=n_cur,n_pchembl=n_pc,n_pchembl_binding=n_pcb))
    print(f"  {tid} {gene:9} {cid:15} current {n_cur:<4} pchembl {n_pc:<4} pchembl_binding {n_pcb:<4} {status}")
    time.sleep(0.2)
with open(os.path.expanduser('~/r2/out/chembl_filtered.csv'),'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print('\nwrote out/chembl_filtered.csv')
