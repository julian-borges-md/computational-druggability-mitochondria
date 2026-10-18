#!/usr/bin/env python3
"""
MitoCorex Chapter 2B.1, reproducible fpocket v4.0 druggability run.

For each of the 25 targets:
  1. select a HUMAN experimental structure containing the target UniProt accession
     (AlphaFold DB model where no human experimental structure exists);
  2. download mmCIF from RCSB, extract only the chains mapping to that accession;
  3. run fpocket v4.0 with default parameters on the extracted chains;
  4. record every pocket, and report the rank 1 pocket and the maximum
     Druggability Score pocket.

Reproduce with:  python3 fpocket_pipeline.py
Requires: fpocket 4.0, gemmi, network access to RCSB and AlphaFold DB.
"""
import csv, json, os, re, subprocess, sys, time, urllib.request
import gemmi

BASE=os.path.expanduser('~/csbrd'); PD=f'{BASE}/fp/pdb'; OD=f'{BASE}/fp/out'
FPOCKET='/tmp/fpocket/bin/fpocket'
os.makedirs(PD,exist_ok=True); os.makedirs(OD,exist_ok=True)
MAXATOM=90000   # extracted-chain atom ceiling; larger entries are skipped for the next candidate

def api(u):
    for _ in range(3):
        try: return json.load(urllib.request.urlopen(u,timeout=40))
        except Exception: time.sleep(1.5)
    return None

def get_cif(pdb):
    p=f'{PD}/{pdb.upper()}.cif'
    if os.path.exists(p) and os.path.getsize(p)>2000: return p
    try:
        urllib.request.urlretrieve(f'https://files.rcsb.org/download/{pdb.upper()}.cif',p); return p
    except Exception: return None

def get_af(uni):
    p=f'{PD}/AF-{uni}.pdb'
    if os.path.exists(p) and os.path.getsize(p)>2000: return p
    for v in ('v6','v4','v3'):
        try:
            urllib.request.urlretrieve(f'https://alphafold.ebi.ac.uk/files/AF-{uni}-F1-model_{v}.pdb',p); return p
        except Exception: continue
    return None

def chains_for(pdb, uni):
    e=api(f'https://data.rcsb.org/rest/v1/core/entry/{pdb}')
    if not e: return [],'',''
    title=(e.get('struct',{}) or {}).get('title','')
    res=(e.get('rcsb_entry_info',{}) or {}).get('resolution_combined') or []
    method=';'.join((e.get('rcsb_entry_info',{}) or {}).get('experimental_method_names') or [])
    for en in e['rcsb_entry_container_identifiers']['polymer_entity_ids']:
        pe=api(f'https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb}/{en}')
        if not pe: continue
        idc=pe.get('rcsb_polymer_entity_container_identifiers',{}) or {}
        if uni in (idc.get('uniprot_ids') or []):
            return (idc.get('auth_asym_ids') or []), title, f"{res[0] if res else ''}|{method}"
    return [], title, f"{res[0] if res else ''}|{method}"

def extract(src, chains, dest):
    st=gemmi.read_structure(src); st.setup_entities()
    if chains:
        keep=set(chains)
        for model in st:
            for ch in [c.name for c in model]:
                if ch not in keep: model.remove_chain(ch)
    st.remove_ligands_and_waters(); st.remove_hydrogens()
    n=sum(len(r) for m in st for c in m for r in c)
    st.setup_entities()
    doc=st.clone(); doc.write_pdb(dest)
    natoms=sum(1 for m in st for c in m for r in c for a in r)
    return natoms

def parse_info(p):
    out,cur=[],None
    if not os.path.exists(p): return out
    for line in open(p):
        m=re.match(r'Pocket (\d+)\s*:',line.strip())
        if m: cur={'pocket_rank':int(m.group(1))}; out.append(cur); continue
        if cur is not None and ':' in line:
            k,_,v=line.strip().partition(':')
            try: cur[k.strip()]=float(v.strip())
            except ValueError: pass
    return out

hs={r['target_id']:r for r in csv.DictReader(open(f'{BASE}/work/human_structures.csv'))}
man={r['target_id']:r for r in csv.DictReader(open(f'{BASE}/work/manifest.csv'))}
# TGT_024: druggable site is the KEAP1 Kelch domain
SUB={'TGT_024':('KEAP1','Q14145')}

results,dump=[],[]
for tid in sorted(man):
    g=man[tid]['gene']; uni=man[tid]['uniprot']
    gene_used, uni_used = g, uni
    if tid in SUB: gene_used, uni_used = SUB[tid]
    cands=(hs.get(tid+'_partner' if tid in SUB else tid,{}).get('human_pdbs','') or '').split(';')
    cands=[c for c in cands if c]
    chosen=None
    for pdb in cands[:8]:
        ch,title,meta=chains_for(pdb, uni_used)
        if not ch: continue
        cif=get_cif(pdb)
        if not cif: continue
        dest=f'{OD}/{tid}_{gene_used}.pdb'
        try: n=extract(cif, ch[:2], dest)
        except Exception: continue
        if n==0 or n>MAXATOM:
            continue
        chosen=(pdb,ch[:2],title,meta,dest,n,'RCSB_human'); break
    if not chosen:
        af=get_af(uni_used)
        if af:
            dest=f'{OD}/{tid}_{gene_used}.pdb'
            n=extract(af,[],dest)
            chosen=(f'AF-{uni_used}',['A'],'AlphaFold DB predicted model','|AlphaFold2',dest,n,'AlphaFoldDB')
    if not chosen:
        results.append({'target_id':tid,'gene':g,'status':'NO_STRUCTURE'}); print(f'{tid} {g} NO STRUCTURE'); continue
    pdb,ch,title,meta,dest,natoms,srctype=chosen
    tag=os.path.basename(dest)[:-4]; outdir=f'{OD}/{tag}_out'
    if not os.path.isdir(outdir):
        t0=time.time()
        try: subprocess.run([FPOCKET,'-f',dest],capture_output=True,text=True,timeout=1200)
        except subprocess.TimeoutExpired:
            results.append({'target_id':tid,'gene':g,'status':'FPOCKET_TIMEOUT','structure':pdb}); print(f'{tid} {g} TIMEOUT'); continue
        el=time.time()-t0
    else: el=0.0
    pk=parse_info(f'{outdir}/{tag}_info.txt')
    if not pk:
        results.append({'target_id':tid,'gene':g,'status':'NO_POCKETS','structure':pdb}); print(f'{tid} {g} no pockets'); continue
    for p_ in pk: dump.append({'target_id':tid,'gene':gene_used,'structure':pdb,**p_})
    r1=pk[0]; best=max(pk,key=lambda x:x.get('Druggability Score',-1))
    resn,meth=(meta.split('|')+[''])[:2]
    results.append({'target_id':tid,'gene':g,'assessed_protein':gene_used,'uniprot':uni_used,
      'tier':man[tid]['tier'],'structure':pdb,'structure_type':srctype,'chains':';'.join(ch),
      'method':meth,'resolution_A':resn,'n_atoms':natoms,'n_pockets':len(pk),
      'rank1_fpocket_score':r1.get('Score'),'rank1_druggability_score':r1.get('Druggability Score'),
      'rank1_volume_A3':r1.get('Volume'),
      'top_druggable_pocket_rank':best.get('pocket_rank'),
      'max_druggability_score':best.get('Druggability Score'),
      'max_pocket_volume_A3':best.get('Volume'),
      'total_sasa_A2':best.get('Total SASA'),'apolar_sasa_A2':best.get('Apolar SASA'),
      'polar_sasa_A2':best.get('Polar SASA'),'hydrophobicity_score':best.get('Hydrophobicity score'),
      'polarity_score':best.get('Polarity score'),'n_alpha_spheres':best.get('Number of Alpha Spheres'),
      'mean_local_hydrophobic_density':best.get('Mean local hydrophobic density'),
      'runtime_s':round(el,1),'title':title[:80],'status':'OK'})
    print(f"{tid} {gene_used:10s} {pdb:10s} chains={';'.join(ch):6s} atoms={natoms:6d} pockets={len(pk):3d} maxDrug={best.get('Druggability Score')} vol={best.get('Volume')}",flush=True)

cols=['target_id','gene','assessed_protein','uniprot','tier','structure','structure_type','chains','method',
 'resolution_A','n_atoms','n_pockets','rank1_fpocket_score','rank1_druggability_score','rank1_volume_A3',
 'top_druggable_pocket_rank','max_druggability_score','max_pocket_volume_A3','total_sasa_A2','apolar_sasa_A2',
 'polar_sasa_A2','hydrophobicity_score','polarity_score','n_alpha_spheres','mean_local_hydrophobic_density',
 'runtime_s','title','status']
with open(f'{BASE}/out/TableS4_fpocket_druggability.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=cols,extrasaction='ignore'); w.writeheader(); w.writerows(results)
if dump:
    k=sorted({x for d in dump for x in d}); lead=['target_id','gene','structure','pocket_rank']
    k=lead+[x for x in k if x not in lead]
    with open(f'{BASE}/out/TableS5_fpocket_all_pockets.csv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=k,extrasaction='ignore'); w.writeheader(); w.writerows(dump)
print('\nDONE ->  out/TableS4_fpocket_druggability.csv  out/TableS5_fpocket_all_pockets.csv')
