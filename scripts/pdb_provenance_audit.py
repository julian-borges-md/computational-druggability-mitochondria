import urllib.request, json, csv, re, time, yaml
def api(u):
    for _ in range(3):
        try: return json.load(urllib.request.urlopen(u, timeout=30))
        except Exception: time.sleep(1.5)
    return None
dr=yaml.safe_load(open('src/druggability_registry_v1.yaml'))['targets']
pw=yaml.safe_load(open('src/pathway_registry_v1.yaml'))['targets']
ms=open('work/MS_APR7_CSBJ.txt').read()
t1=re.search(r'Target \| Gene \| Tier.*?<<<ENDTABLE>>>',ms,re.S).group(0)
rows=[r.split(' | ') for r in t1.split('\n') if r.startswith('TGT_')]
msmap={r[0].strip():{'gene':r[1].strip(),'tier':r[2].strip(),'pdb_col':r[4].strip()} for r in rows}
# alias sets for orthologue recognition
ALIAS={'MT-ND1':['ND1','NADH-ubiquinone oxidoreductase chain 1','NADH dehydrogenase subunit 1'],
 'MT-ND4':['ND4','NADH-ubiquinone oxidoreductase chain 4'],
 'MT-CYB':['CYTB','cytochrome b','CYTOCHROME B'],
 'MT-CO1':['COX1','COI','cytochrome c oxidase subunit 1','Cytochrome c oxidase subunit 1'],
 'MT-ATP6':['ATP6','ATP synthase subunit a','ATP synthase F0 subunit 6'],
 'NDUFV1':['NADH dehydrogenase [ubiquinone] flavoprotein 1','51 kDa'],
 'SDHA':['succinate dehydrogenase','FAD-binding protein','flavoprotein subunit'],
 'UQCRC2':['Cytochrome b-c1 complex subunit 2','core protein 2'],
 'COX4I1':['Cytochrome c oxidase subunit 4'],'ATP5F1A':['ATP synthase subunit alpha'],
 'POLG':['DNA polymerase subunit gamma'],'TWNK':['Twinkle','TWNK','C10orf2'],
 'TFAM':['Transcription factor A'],'TIMM23':['Mitochondrial import inner membrane translocase subunit Tim23','TIM23'],
 'TOMM40':['Mitochondrial import receptor subunit TOM40','TOM40'],'DNM1L':['Dynamin-1-like','DRP1'],
 'OPA1':['Dynamin-like 120 kDa','OPA1'],'PINK1':['PTEN-induced','PINK1'],
 'SOD2':['Superoxide dismutase [Mn]','superoxide dismutase'],'PPARGC1A':['PGC-1'],
 'SLC25A4':['ADP/ATP translocase','adenine nucleotide translocator'],'VDAC1':['anion-selective channel','VDAC'],
 'LONP1':['Lon protease'],'NFE2L2':['Nuclear factor erythroid','NFE2L2','Kelch-like ECH-associated protein 1','KEAP1'],
 'CLPP':['Clp protease','CLPP']}
out=[]
for tid in sorted(msmap):
    g=msmap[tid]['gene']; uni=pw[tid].get('uniprot')
    ids=re.findall(r'\b([0-9][A-Za-z0-9]{3})\b', msmap[tid]['pdb_col'])
    ids+= [x for x in (dr[tid].get('key_pdb_structures') or []) if re.fullmatch(r'[0-9][A-Za-z0-9]{3}',str(x))]
    seen=[]
    for p in ids:
        if p.upper() in [s.upper() for s in seen]: continue
        seen.append(p)
        e=api(f'https://data.rcsb.org/rest/v1/core/entry/{p}')
        if not e:
            out.append({'target_id':tid,'gene':g,'pdb':p,'verdict':'NOT_FOUND_IN_RCSB','title':'','organisms':'','matched':'','cited_in':'Table1/registry'}); 
            print(f"{tid} {g:10s} {p}  NOT FOUND"); continue
        title=(e.get('struct',{}) or {}).get('title','')
        ents=e['rcsb_entry_container_identifiers']['polymer_entity_ids']
        orgs=set(); accs=[]; descs=[]
        for en in ents:
            pe=api(f'https://data.rcsb.org/rest/v1/core/polymer_entity/{p}/{en}')
            if not pe: continue
            idc=pe.get('rcsb_polymer_entity_container_identifiers',{}) or {}
            accs+= list(idc.get('uniprot_ids') or [])
            descs.append((pe.get('rcsb_polymer_entity',{}) or {}).get('pdbx_description','') or '')
            for s in (pe.get('rcsb_entity_source_organism') or []):
                if s.get('scientific_name'): orgs.add(s['scientific_name'])
        exact = uni in accs
        blob=' | '.join(descs).lower()+' '+title.lower()
        ortho = any(a.lower() in blob for a in ALIAS.get(g,[g]))
        verdict = 'EXACT_HUMAN' if exact else ('ORTHOLOGUE_OR_ALIAS' if ortho else 'PROTEIN_NOT_PRESENT')
        human = 'Homo sapiens' in orgs
        out.append({'target_id':tid,'gene':g,'human_uniprot':uni,'pdb':p,'verdict':verdict,
                    'contains_human':human,'organisms':';'.join(sorted(orgs))[:60],
                    'title':title[:75],'entity_uniprots':';'.join(accs)[:70]})
        print(f"{tid} {g:10s} {p}  {verdict:20s} org={';'.join(sorted(orgs))[:32]:34s} {title[:52]}", flush=True)
with open('out/pdb_provenance_audit.csv','w',newline='') as f:
    k=sorted({x for d in out for x in d}); lead=['target_id','gene','pdb','verdict','contains_human','organisms','title','human_uniprot','entity_uniprots']
    k=lead+[x for x in k if x not in lead]
    w=csv.DictWriter(f,fieldnames=k,extrasaction='ignore'); w.writeheader(); w.writerows(out)
print("\nwrote out/pdb_provenance_audit.csv")
