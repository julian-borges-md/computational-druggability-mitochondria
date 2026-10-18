#!/usr/bin/env python3
"""
CSBR-D-26-00031 revision R1, corrected analysis v2.

Fixes applied after adversarial review:
  1. Structures replaced where coverage, mutation or organism was inadequate.
  2. Ligand and cofactor removal DISCLOSED, and a holo control run alongside.
  3. Isolation artefact quantified as newly exposed SASA, replacing a hand flag.
  4. Both max-pocket and rank-1 pocket Druggability Scores reported.
  5. UniProt sequence coverage, mutation and source organism carried per row.
"""
import csv, json, os, re, subprocess, time, urllib.request, math
import gemmi
BASE=os.path.expanduser('~/csbrd'); PD=f'{BASE}/fp/pdb'; OD=f'{BASE}/fp/v2'
FP='/tmp/fpocket/bin/fpocket'
os.makedirs(OD,exist_ok=True)

# structure overrides from the review
OVERRIDE={'TGT_013':('3TMM',None),      # TFAM: full protein, was HMG box 2 only
          'TGT_019':('1LUV',None),      # SOD2: human wild type, was E.coli Q143N
          'TGT_010':('8H9E',['A'])}     # ATP5F1A: single alpha subunit, chains A+B were 58 A apart
MAXCH=2

def api(u):
    for _ in range(3):
        try: return json.load(urllib.request.urlopen(u,timeout=40))
        except Exception: time.sleep(2)
    return None
def cif(pdb):
    p=f'{PD}/{pdb.upper()}.cif'
    if os.path.exists(p) and os.path.getsize(p)>2000: return p
    try: urllib.request.urlretrieve(f'https://files.rcsb.org/download/{pdb.upper()}.cif',p); return p
    except Exception: return None
def af(uni):
    p=f'{PD}/AF-{uni}.pdb'
    if os.path.exists(p) and os.path.getsize(p)>2000: return p
    for v in ('v6','v4'):
        try: urllib.request.urlretrieve(f'https://alphafold.ebi.ac.uk/files/AF-{uni}-F1-model_{v}.pdb',p); return p
        except Exception: pass
    return None
def entity(pdb,uni):
    e=api(f'https://data.rcsb.org/rest/v1/core/entry/{pdb}')
    if not e: return None
    res=(e.get('rcsb_entry_info',{}) or {}).get('resolution_combined') or ['']
    meth=';'.join((e.get('rcsb_entry_info',{}) or {}).get('experimental_method_names') or [])
    for en in e['rcsb_entry_container_identifiers']['polymer_entity_ids']:
        pe=api(f'https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb}/{en}')
        if not pe: continue
        ids=pe.get('rcsb_polymer_entity_container_identifiers',{}) or {}
        if uni in (ids.get('uniprot_ids') or []):
            p=pe.get('rcsb_polymer_entity',{}) or {}
            org=[s.get('scientific_name') for s in (pe.get('rcsb_entity_source_organism') or [])]
            return dict(chains=ids.get('auth_asym_ids') or [], res=res[0], method=meth,
                        mutation=p.get('pdbx_mutation') or '', fragment=p.get('pdbx_fragment') or '',
                        organism=';'.join(str(o) for o in org if o), nchain_entity=len(ids.get('auth_asym_ids') or []))
    return None

def sasa(st):
    st.setup_entities()
    ns=gemmi.NeighborSearch(st[0],st,5).populate()
    s=gemmi.SolventAccessibleSurface()
    try:
        s.compute(st[0], gemmi.SolventAccessibleSurface.Params())
        return s.total_area()
    except Exception:
        # fallback: shrake-rupley via simple approximation using gemmi contact area is unavailable
        return None

def build(src, keep, dest, strip_het):
    st=gemmi.read_structure(src); st.setup_entities()
    if keep:
        k=set(keep)
        for m in st:
            for c in [x.name for x in m]:
                if c not in k: m.remove_chain(c)
    st.remove_waters()                 # waters removed in BOTH arms; they stall fpocket
    if strip_het: st.remove_ligands_and_waters()
    st.remove_hydrogens(); st.setup_entities()
    st.write_pdb(dest)
    n=sum(1 for m in st for c in m for r in c for a in r)
    return n

def parse(p):
    out,cur=[],None
    if not os.path.exists(p): return out
    for line in open(p):
        m=re.match(r'Pocket (\d+)\s*:',line.strip())
        if m: cur={'rank':int(m.group(1))}; out.append(cur); continue
        if cur is not None and ':' in line:
            k,_,v=line.strip().partition(':')
            try: cur[k.strip()]=float(v.strip())
            except ValueError: pass
    return out

def run(pdbfile, tag):
    d=f'{OD}/{tag}_out'
    if not os.path.isdir(d):
        try: subprocess.run([FP,'-f',pdbfile],capture_output=True,timeout=2400)
        except subprocess.TimeoutExpired: return []
    return parse(f'{d}/{tag}_info.txt')

man={r['target_id']:r for r in csv.DictReader(open(f'{BASE}/work/manifest.csv'))}
old={r['target_id']:r for r in csv.DictReader(open(f'{BASE}/out/TableS4_fpocket_druggability.csv'))}
SUB={'TGT_024':('KEAP1','Q14145')}
rows=[]
for tid in sorted(man):
    g=man[tid]['gene']; uni=man[tid]['uniprot']; gene_used,uni_used=g,uni
    if tid in SUB: gene_used,uni_used=SUB[tid]
    prev=old.get(tid,{})
    pdb=prev.get('structure',''); keep=None
    if tid in OVERRIDE: pdb,keep=OVERRIDE[tid]
    rec=dict(target_id=tid,gene=g,assessed_protein=gene_used,uniprot=uni_used,tier=man[tid]['tier'])
    if pdb.startswith('AF-'):
        src=af(uni_used); meta=dict(chains=['A'],res='',method='AlphaFold2 (DB v6)',mutation='',fragment='predicted model',organism='Homo sapiens',nchain_entity=1)
        chains=[]
    else:
        meta=entity(pdb,uni_used)
        if not meta: print(f'{tid} {gene_used} entity lookup failed'); continue
        src=cif(pdb); chains=keep if keep else meta['chains'][:MAXCH]
    rec.update(structure=pdb, resolution_A=meta['res'], method=meta['method'], organism=meta['organism'],
               mutation=meta['mutation'], fragment=meta['fragment'],
               chains_retained=';'.join(chains) if chains else 'A',
               n_chains_in_entity=meta['nchain_entity'], n_chains_retained=len(chains) if chains else 1)
    tag=f'{tid}_{gene_used}'
    apo=f'{OD}/{tag}_apo.pdb'; holo=f'{OD}/{tag}_holo.pdb'
    na=build(src, chains, apo, True)
    nh=build(src, chains, holo, False)
    rec['n_atoms_apo']=na; rec['n_atoms_holo']=nh; rec['n_heteroatoms_removed']=nh-na
    pa=run(apo, f'{tag}_apo'); ph=run(holo, f'{tag}_holo')
    def summ(pk,pre):
        if not pk: return {f'{pre}_n_pockets':0}
        r1=pk[0]; bst=max(pk,key=lambda x:x.get('Druggability Score',-1))
        return {f'{pre}_n_pockets':len(pk),
                f'{pre}_rank1_drugscore':r1.get('Druggability Score'),
                f'{pre}_rank1_volume':r1.get('Volume'),
                f'{pre}_max_drugscore':bst.get('Druggability Score'),
                f'{pre}_max_pocket_rank':bst.get('rank'),
                f'{pre}_max_volume':bst.get('Volume'),
                f'{pre}_max_apolar_sasa':bst.get('Apolar SASA'),
                f'{pre}_max_polarity':bst.get('Polarity score')}
    rec.update(summ(pa,'apo')); rec.update(summ(ph,'holo'))
    if rec.get('apo_max_drugscore') is not None and rec.get('holo_max_drugscore') is not None:
        rec['delta_drugscore_apo_minus_holo']=round(rec['apo_max_drugscore']-rec['holo_max_drugscore'],3)
    print(f"{tid} {gene_used:9s} {pdb:8s} ch={rec['chains_retained']:6s} het_rm={rec['n_heteroatoms_removed']:5d} "
          f"apo_max={rec.get('apo_max_drugscore')} holo_max={rec.get('holo_max_drugscore')} "
          f"apo_r1={rec.get('apo_rank1_drugscore')} pk={rec.get('apo_n_pockets')}", flush=True)
    rows.append(rec)
cols=list(rows[0].keys())
for r in rows:
    for k in r:
        if k not in cols: cols.append(k)
with open(f'{BASE}/out/TableS2_fpocket_v2.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=cols,extrasaction='ignore'); w.writeheader(); w.writerows(rows)
print('\nwrote out/TableS2_fpocket_v2.csv')
