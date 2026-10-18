#!/usr/bin/env python3
"""Final supplementary tables for the second revision.
S2  fpocket per target, apo and holo          (from TableS2_fpocket_v2.csv)
S3  all pockets, both arms                    (already written)
S5  quantitative matrix, plus coverage        (from TableS5_quantitative_matrix_v2.csv)
S12 structural coverage                       (from coverage_v2.csv)
S13 pre-specified sensitivity analyses        (computed here)
"""
import csv, math, os, json
from scipy import stats

BASE = os.path.expanduser('~/csbrd')
OUT  = f'{BASE}/out_v2'
COV  = {r['target_id']: r for r in csv.DictReader(open(f'{OUT}/coverage_v2.csv'))}

# ---------- S2 ----------------------------------------------------------
s2 = list(csv.DictReader(open(f'{BASE}/out/TableS2_fpocket_v2.csv')))
for r in s2:
    c = COV[r['target_id']]
    r['fragment'] = r['fragment'] or 'not stated in the entry'
    r['mutation'] = r['mutation'] or 'none recorded in the entry'
    r['uniprot_length'] = c['uniprot_length']
    r['uniprot_range_declared'] = c['uniprot_range_modelled']
    r['residues_resolved_per_chain'] = c['n_residues_per_chain']
    r['coverage_percent'] = c['coverage_percent']
cols2 = list(s2[0].keys())
with open(f'{OUT}/SupplementaryTable_S2_fpocket_per_target.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=cols2); w.writeheader(); w.writerows(s2)

# ---------- S5 ----------------------------------------------------------
s5 = list(csv.DictReader(open(f'{BASE}/out/TableS5_quantitative_matrix_v2.csv')))
for r in s5:
    c = COV[r['target_id']]
    r['fragment'] = r['fragment'] if r['fragment'] != 'full construct' else (
        'full construct' if float(c['coverage_percent']) >= 95 else 'not stated in the entry')
    r['mutation'] = r['mutation'] if r['mutation'] != 'none' else 'none recorded in the entry'
    r['uniprot_length'] = c['uniprot_length']
    r['coverage_percent'] = c['coverage_percent']
    r['coverage_below_60pc'] = 'yes' if float(c['coverage_percent']) < 60 else 'no'
cols5 = list(s5[0].keys())
with open(f'{OUT}/SupplementaryTable_S5_quantitative_matrix.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=cols5); w.writeheader(); w.writerows(s5)

# ---------- S12 ---------------------------------------------------------
c12 = list(csv.DictReader(open(f'{OUT}/coverage_v2.csv')))
with open(f'{OUT}/SupplementaryTable_S12_structural_coverage.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(c12[0].keys())); w.writeheader(); w.writerows(c12)

# ---------- S13 sensitivity --------------------------------------------
R = {r['target_id']: r for r in s5}
def sp(x, y):
    rho, p = stats.spearmanr(x, y); n = len(x)
    z = math.atanh(rho); se = 1/math.sqrt(n-3)
    return round(rho,3), round(p,4), round(math.tanh(z-1.96*se),3), round(math.tanh(z+1.96*se),3)

SUBSETS = [
 ('all targets', 'no restriction', lambda t: True),
 ('coverage at or above 60 per cent', 'excludes LONP1, OPA1, DNM1L and KEAP1', 
  lambda t: float(COV[t]['coverage_percent']) >= 60),
 ('experimental structures only', 'excludes the three AlphaFold models',
  lambda t: not R[t]['structure'].startswith('AF-')),
 ('single retained chain only', 'excludes every assessment on two chains',
  lambda t: int(R[t]['n_chains_retained']) == 1),
]
RELS = [
 ('tier vs pocket druggability D1', lambda t: int(R[t]['tier_published']), lambda t: float(R[t]['D1_apo_max'])),
 ('tier vs precedent D4',           lambda t: int(R[t]['tier_published']), lambda t: int(R[t]['D4_precedent'])),
 ('tier vs delivery D5',            lambda t: int(R[t]['tier_published']), lambda t: int(R[t]['D5_delivery'])),
 ('D1 vs D5 (withdrawn claim)',     lambda t: float(R[t]['D1_apo_max']),   lambda t: int(R[t]['D5_delivery'])),
]
rows = []
for sname, snote, keep in SUBSETS:
    T = [t for t in sorted(R) if keep(t)]
    for rname, fx, fy in RELS:
        rho, p, lo, hi = sp([fx(t) for t in T], [fy(t) for t in T])
        rows.append(dict(subset=sname, subset_note=snote, n=len(T), relationship=rname,
                         spearman_rho=rho, p_value=p, ci95_low=lo, ci95_high=hi,
                         significant_at_0_05='yes' if p < 0.05 else 'no'))
# coverage as a confound
T = sorted(R)
for nm, fy in (('D1 apo max', lambda t: float(R[t]['D1_apo_max'])),
               ('D2 max pocket volume', lambda t: float(R[t]['D2_apo_max_volume_A3'])),
               ('n pockets detected', lambda t: int(R[t]['n_pockets']))):
    rho, p, lo, hi = sp([float(COV[t]['coverage_percent']) for t in T], [fy(t) for t in T])
    rows.append(dict(subset='confound check', subset_note='does construct coverage drive the structural scores',
                     n=25, relationship=f'coverage percent vs {nm}', spearman_rho=rho, p_value=p,
                     ci95_low=lo, ci95_high=hi, significant_at_0_05='yes' if p < 0.05 else 'no'))
with open(f'{OUT}/SupplementaryTable_S13_sensitivity_analyses.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

json.dump(rows, open(f'{OUT}/sensitivity.json', 'w'), indent=1)
for f_ in ('SupplementaryTable_S2_fpocket_per_target.csv', 'SupplementaryTable_S5_quantitative_matrix.csv',
           'SupplementaryTable_S12_structural_coverage.csv', 'SupplementaryTable_S13_sensitivity_analyses.csv'):
    print(f'{os.path.getsize(f"{OUT}/{f_}"):>8} {f_}')
print()
for r in rows:
    if 'D1' in r['relationship'] or 'coverage' in r['relationship']:
        print(f"{r['subset']:34s} n={r['n']:2d} {r['relationship']:34s} rho={r['spearman_rho']:6.3f} p={r['p_value']:.4f}")
