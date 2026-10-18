#!/usr/bin/env python3
"""UniProt active-site / binding-site / site features for the assessed protein
of each target: the objective anchor for the a priori therapeutic site (C2)."""
import csv, json, os, time, urllib.request

S2 = list(csv.DictReader(open(os.path.expanduser('~/r2/data/SupplementaryTable_S2_fpocket_per_target.csv'))))
out = {}
for r in S2:
    tid, acc, gene = r['target_id'], r['uniprot'], r['gene']
    url = f'https://rest.uniprot.org/uniprotkb/{acc}.json?fields=ft_act_site,ft_binding,ft_site'
    feats, err = [], ''
    try:
        d = json.load(urllib.request.urlopen(url, timeout=20))
        for f in d.get('features', []):
            loc = f['location']
            feats.append(dict(type=f['type'],
                              start=loc['start'].get('value'), end=loc['end'].get('value'),
                              desc=f.get('description', ''),
                              ligand=(f.get('ligand') or {}).get('name', '')))
    except Exception as ex:
        err = str(ex)[:60]
    out[tid] = dict(gene=gene, uniprot=acc, features=feats)
    res = sorted({f['start'] for f in feats if f['start']})
    print(f"  {tid} {gene:9} {acc:7} {len(feats):3d} feats"
          + (f"  res {res[:12]}" if res else "  (no residue-level site)")
          + (f"  ERR {err}" if err else ""))
    time.sleep(0.2)
json.dump(out, open(os.path.expanduser('~/r2/out/uniprot_sites.json'), 'w'), indent=1)
print('\nwrote out/uniprot_sites.json')
