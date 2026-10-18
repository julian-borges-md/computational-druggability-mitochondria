#!/usr/bin/env python3
"""Verify the Tier-D ligand anchors are actually present near the assessed chain,
and list candidate anchor ligands. Validates the a priori site table before the
overlap computation is committed."""
import csv, os, gemmi
S2 = {r['target_id']: r for r in csv.DictReader(open(os.path.expanduser('~/r2/data/SupplementaryTable_S2_fpocket_per_target.csv')))}
AP = {r['target_id']: r for r in csv.DictReader(open(os.path.expanduser('~/r2/out/apriori_sites.csv')))}
D = os.path.expanduser('~/r2/struct')
WATER = {'HOH','DOD','WAT'}

def hetatms_near(path, chains, cutoff=5.0):
    st = gemmi.read_structure(path); st.setup_entities()
    model = st[0]
    keep = set(chains)
    # target atoms = all atoms in kept chains that are polymer (amino acids)
    tgt = []
    for ch in model:
        if ch.name in keep:
            for res in ch:
                info = gemmi.find_tabulated_residue(res.name)
                if info and info.is_amino_acid():
                    for a in res: tgt.append(a.pos)
    if not tgt: return {}
    hits = {}
    for ch in model:
        for res in ch:
            info = gemmi.find_tabulated_residue(res.name)
            if (info and info.is_amino_acid()) or res.name in WATER: continue
            # non-polymer / hetero residue: is any atom within cutoff of the target chain?
            near = False
            for a in res:
                for t in tgt:
                    if a.pos.dist(t) <= cutoff: near = True; break
                if near: break
            if near:
                hits[res.name] = hits.get(res.name, 0) + 1
    return hits

for tid in sorted(AP):
    ap = AP[tid]; s2 = S2[tid]
    if ap['tier'] != 'D':
        print(f"  {tid} {ap['gene']:9} tier {ap['tier']}  ({ap['site_name']})")
        continue
    struct = ap['structure']
    if struct.startswith('AF-'):
        print(f"  {tid} {ap['gene']:9} tier D  AlphaFold, no ligands"); continue
    chains = s2['chains_retained'].split(';')
    hits = hetatms_near(f'{D}/{struct}.cif', chains)
    top = sorted(hits.items(), key=lambda x: -x[1])[:8]
    print(f"  {tid} {ap['gene']:9} {struct}  chains {s2['chains_retained']:6}  near-chain ligands: " +
          ', '.join(f'{k}x{v}' for k, v in top) if top else
          f"  {tid} {ap['gene']:9} {struct}  NO near-chain heteroatoms")
