import urllib.request, json, csv, time, yaml
SEARCH='https://search.rcsb.org/rcsbsearch/v2/query'
def post(q):
    for _ in range(3):
        try:
            r=urllib.request.Request(SEARCH, data=json.dumps(q).encode(), headers={'Content-Type':'application/json'})
            return json.load(urllib.request.urlopen(r, timeout=40))
        except Exception as e:
            time.sleep(2)
    return None
def api(u):
    for _ in range(3):
        try: return json.load(urllib.request.urlopen(u,timeout=30))
        except Exception: time.sleep(1.5)
    return None
def find(uniprot, limit=25):
    q={"query":{"type":"terminal","service":"text","parameters":{
        "attribute":"rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
        "operator":"exact_match","value":uniprot}},
       "return_type":"entry","request_options":{"paginate":{"start":0,"rows":limit},
       "results_content_type":["experimental"],
       "sort":[{"sort_by":"rcsb_entry_info.resolution_combined","direction":"asc"}]}}
    r=post(q)
    return [x['identifier'] for x in (r or {}).get('result_set',[])]

TARGETS=[('TGT_001','MT-ND1','P03886'),('TGT_002','MT-ND4','P03905'),('TGT_003','MT-CYB','P00156'),
('TGT_004','MT-CO1','P00395'),('TGT_005','MT-ATP6','P00846'),('TGT_006','NDUFV1','P49821'),
('TGT_007','SDHA','P31040'),('TGT_008','UQCRC2','P22695'),('TGT_009','COX4I1','P13073'),
('TGT_010','ATP5F1A','P25705'),('TGT_011','POLG','P54098'),('TGT_012','TWNK','Q96RR1'),
('TGT_013','TFAM','Q00059'),('TGT_014','TIMM23','O14925'),('TGT_015','TOMM40','O96008'),
('TGT_016','DNM1L','O00429'),('TGT_017','OPA1','O60313'),('TGT_018','PINK1','Q9BXM7'),
('TGT_019','SOD2','P04179'),('TGT_020','PPARGC1A','Q9UBK2'),('TGT_021','SLC25A4','P12235'),
('TGT_022','VDAC1','P21796'),('TGT_023','LONP1','P36776'),('TGT_024','NFE2L2','Q16236'),
('TGT_025','CLPP','Q16740')]
EXTRA={'TGT_024':('KEAP1','Q14145')}
rows=[]
for tid,g,uni in TARGETS:
    hits=find(uni)
    note=''
    if tid in EXTRA and not hits:
        pass
    rows.append({'target_id':tid,'gene':g,'uniprot':uni,'n_human_structures':len(hits),'human_pdbs':';'.join(hits[:12])})
    print(f"{tid} {g:10s} {uni:8s} human_structures={len(hits):3d}  {';'.join(hits[:10])}", flush=True)
    if tid in EXTRA:
        g2,u2=EXTRA[tid]; h2=find(u2)
        rows.append({'target_id':tid+'_partner','gene':g2,'uniprot':u2,'n_human_structures':len(h2),'human_pdbs':';'.join(h2[:12])})
        print(f"{tid} {g2:10s} {u2:8s} human_structures={len(h2):3d}  {';'.join(h2[:10])}", flush=True)
with open('work/human_structures.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['target_id','gene','uniprot','n_human_structures','human_pdbs']); w.writeheader(); w.writerows(rows)
print('\nwrote work/human_structures.csv')
