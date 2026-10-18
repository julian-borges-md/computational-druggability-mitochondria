#!/usr/bin/env python3
"""Quantitative ligand efficiency evidence for the 25 targets, from ChEMBL.

For each target UniProt accession: resolve the Homo sapiens SINGLE PROTEIN
target, then count bioactivity records and distinct compounds with a potency
measurement at or below 10 micromolar. Reproducible against the live API.
"""
import json, time, urllib.parse, urllib.request, csv, sys
BASE = "https://www.ebi.ac.uk/chembl/api/data"
def get(path, **params):
    url = f"{BASE}/{path}.json?" + urllib.parse.urlencode(params)
    for _ in range(4):
        try:
            return json.load(urllib.request.urlopen(url, timeout=60))
        except Exception:
            time.sleep(3)
    return None

TARGETS = [
 ("TGT_001","MT-ND1","P03886"),("TGT_002","MT-ND4","P03905"),("TGT_003","MT-CYB","P00156"),
 ("TGT_004","MT-CO1","P00395"),("TGT_005","MT-ATP6","P00846"),("TGT_006","NDUFV1","P49821"),
 ("TGT_007","SDHA","P31040"),("TGT_008","UQCRC2","P22695"),("TGT_009","COX4I1","P13073"),
 ("TGT_010","ATP5F1A","P25705"),("TGT_011","POLG","P54098"),("TGT_012","TWNK","Q96RR1"),
 ("TGT_013","TFAM","Q00059"),("TGT_014","TIMM23","O14925"),("TGT_015","TOMM40","O96008"),
 ("TGT_016","DNM1L","O00429"),("TGT_017","OPA1","O60313"),("TGT_018","PINK1","Q9BXM7"),
 ("TGT_019","SOD2","P04179"),("TGT_020","PPARGC1A","Q9UBK2"),("TGT_021","SLC25A4","P12235"),
 ("TGT_022","VDAC1","P21796"),("TGT_023","LONP1","P36776"),
 ("TGT_024","KEAP1","Q14145"),   # ligand binding target for the NFE2L2 entry
 ("TGT_025","CLPP","Q16740"),
]
rows=[]
for tid, gene, uni in TARGETS:
    t = get("target", target_components__accession=uni, organism="Homo sapiens", limit=20)
    cands = [x for x in (t or {}).get("targets", []) if x.get("target_type")=="SINGLE PROTEIN"] \
            or (t or {}).get("targets", [])
    if not cands:
        rows.append({"target_id":tid,"gene":gene,"uniprot":uni,"chembl_id":"","n_activities_10uM":0,
                     "n_distinct_compounds":0,"n_activities_any":0,"status":"NO_CHEMBL_TARGET"})
        print(f"{tid} {gene:9s} {uni:8s} no ChEMBL target", flush=True); continue
    cid = cands[0]["target_chembl_id"]
    a10 = get("activity", target_chembl_id=cid, standard_type__in="IC50,Ki,Kd,EC50",
              standard_units="nM", standard_value__lte=10000, limit=1)
    aany = get("activity", target_chembl_id=cid, standard_type__in="IC50,Ki,Kd,EC50", limit=1)
    n10 = (a10 or {}).get("page_meta",{}).get("total_count",0)
    nany= (aany or {}).get("page_meta",{}).get("total_count",0)
    cset=set()
    if n10:
        off=0
        while off < min(n10, 1000):
            page = get("activity", target_chembl_id=cid, standard_type__in="IC50,Ki,Kd,EC50",
                       standard_units="nM", standard_value__lte=10000, limit=1000, offset=off)
            for r in (page or {}).get("activities", []):
                if r.get("molecule_chembl_id"): cset.add(r["molecule_chembl_id"])
            off += 1000
            if not page or len(page.get("activities",[]))<1000: break
    rows.append({"target_id":tid,"gene":gene,"uniprot":uni,"chembl_id":cid,
                 "n_activities_10uM":n10,"n_distinct_compounds":len(cset),
                 "n_activities_any":nany,"status":"OK"})
    print(f"{tid} {gene:9s} {uni:8s} {cid:15s} activities<=10uM={n10:<6d} distinct compounds={len(cset):<5d} any potency records={nany}", flush=True)
    time.sleep(0.3)

with open("out/chembl_ligand_evidence.csv","w",newline="") as f:
    w=csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("\nwrote out/chembl_ligand_evidence.csv")
