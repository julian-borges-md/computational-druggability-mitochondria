#!/usr/bin/env python3
"""Per target construct coverage, computed from the deposited structure files
held locally. No network call: RCSB and UniProt are blocked by the session
egress policy, so nothing here is re-fetched.

Reports what the entry itself states:
  uniprot_range_modelled   from _struct_ref_seq db_align_beg/end
  n_residues_observed      counted in the retained chains
  construct_length         from _entity_poly
The full UniProt sequence length is reported only where the local files
establish it (AlphaFold models span the complete sequence by construction).
"""
import csv, os, gemmi

BASE = os.path.expanduser('~/csbrd')
import json
UL = json.load(open(f'{BASE}/work/uniprot_lengths.json'))  # UniProt REST, retrieved 11 Sep 2026
V2   = {r['target_id']: r for r in csv.DictReader(open(f'{BASE}/out/TableS2_fpocket_v2.csv'))}
rows = []

for tid in sorted(V2):
    r = V2[tid]
    pdb, uni = r['structure'], r['uniprot']
    keep = [c for c in r['chains_retained'].split(';') if c]
    rec = dict(target_id=tid, gene=r['gene'], assessed_protein=r['assessed_protein'],
               uniprot=uni, structure=pdb, chains_retained=r['chains_retained'])

    if pdb.startswith('AF-'):
        p = f'{BASE}/fp/pdb/{pdb}.pdb'
        st = gemmi.read_structure(p); st.setup_entities()
        n = sum(1 for c in st[0] for res in c)
        rec.update(source='AlphaFold DB model_v6', uniprot_length=UL[uni],
                   uniprot_range_modelled=f'1-{n}', n_residues_observed=n,
                   construct_length=n, coverage_percent=100.0,
                   n_residues_per_chain=f'A:{n}', n_chains_retained=1,
                   coverage_basis='predicted model spans the complete UniProt sequence by construction')
        rows.append(rec); continue

    p = f'{BASE}/fp/pdb/{pdb}.cif'
    doc = gemmi.cif.read(p); blk = doc.sole_block()

    # entity id whose _struct_ref matches this UniProt accession
    ent, beg, end = None, None, None
    refs = {}
    for rid, acc in blk.find('_struct_ref.', ['id', 'pdbx_db_accession']):
        refs[rid] = acc
    for rid, eid, b, e in blk.find('_struct_ref_seq.',
                                   ['ref_id', '?pdbx_strand_id', 'db_align_beg', 'db_align_end']):
        if refs.get(rid) == uni:
            try:
                b, e = int(b), int(e)
            except ValueError:
                continue
            beg = b if beg is None else min(beg, b)
            end = e if end is None else max(end, e)
    # construct length for the matching entity
    clen = None
    for eid, seq in blk.find('_entity_poly.', ['entity_id', '?pdbx_seq_one_letter_code_can']):
        pass

    st = gemmi.read_structure(p); st.setup_entities()
    per = {}
    for c in st[0]:
        if keep and c.name not in keep:
            continue
        per[c.name] = sum(1 for res in c
                          if (lambda i: i and i.is_amino_acid())(gemmi.find_tabulated_residue(res.name)))
    obs = sum(per.values())

    ul = UL[uni]
    best = max(per.values()) if per else 0
    rec.update(source='RCSB entry _struct_ref_seq alignment; UniProt REST for reference length',
               uniprot_length=ul, 
               uniprot_range_modelled=(f'{beg}-{end}' if beg else 'not stated in the entry'),
               n_residues_observed=obs,
               n_residues_per_chain=';'.join(f'{k}:{v}' for k, v in sorted(per.items())),
               n_chains_retained=len(per),
               construct_length=(end - beg + 1) if beg else '',
               coverage_percent=round(100.0 * best / ul, 1),
               coverage_basis=('amino acid residues resolved in the single best retained chain, divided by '
                               'the UniProt reference sequence length'))
    rows.append(rec)

cols = ['target_id', 'gene', 'assessed_protein', 'uniprot', 'structure', 'chains_retained',
        'source', 'uniprot_length', 'uniprot_range_modelled', 'n_residues_observed',
        'n_residues_per_chain', 'n_chains_retained', 'construct_length', 'coverage_percent',
        'coverage_basis']
out = f'{BASE}/out_v2/coverage_v2.csv'
with open(out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)

for r in rows:
    print(f"{r['target_id']} {r['gene']:10s} {r['structure']:9s} "
          f"range={str(r['uniprot_range_modelled']):12s} obs={r['n_residues_observed']:5d} "
          f"per_chain={str(r['n_residues_per_chain']):18s} "
          f"len={r['uniprot_length']:5d} cov={r['coverage_percent']:5.1f}%"
          + ('   << PARTIAL' if float(r['coverage_percent']) < 60 else ''))
print('\nwrote', out)
