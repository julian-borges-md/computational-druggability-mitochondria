#!/usr/bin/env python3
"""Download the assessed structures (CIF), AlphaFold models (PDB) and PAE JSON.
files.rcsb.org and alphafold.ebi.ac.uk confirmed reachable."""
import csv, os, urllib.request, time
S2 = list(csv.DictReader(open(os.path.expanduser('~/r2/data/SupplementaryTable_S2_fpocket_per_target.csv'))))
D = os.path.expanduser('~/r2/struct')
def get(url, path):
    if os.path.exists(path) and os.path.getsize(path) > 2000: return 'cached'
    try:
        urllib.request.urlretrieve(url, path); return f'{os.path.getsize(path)//1024} KB'
    except Exception as ex: return f'FAIL {str(ex)[:40]}'
for r in S2:
    s = r['structure']
    if s.startswith('AF-'):
        acc = s.split('-')[1] if '-' in s else r['uniprot']
        m = get(f'https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v6.pdb', f'{D}/{s}.pdb')
        p = get(f'https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-predicted_aligned_error_v6.json', f'{D}/{s}_pae.json')
        print(f"  {r['target_id']} {r['gene']:9} {s:11} model {m:12} pae {p}")
    else:
        c = get(f'https://files.rcsb.org/download/{s}.cif', f'{D}/{s}.cif')
        print(f"  {r['target_id']} {r['gene']:9} {s:11} cif {c}")
    time.sleep(0.2)
print('done')
