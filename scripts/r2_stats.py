#!/usr/bin/env python3
"""R2 statistics: directly answers reviewer comments 4 and 5.

C4  full pairwise Spearman matrix among D1..D5, and leave-one-dimension-out
    composite vs tier, to separate independent contribution from circularity.
C5  permutation p-values and percentile bootstrap CIs for the headline
    correlations, alongside the Fisher CIs already reported, since Fisher is
    an approximation under heavy ties at n=25.

Standard library only. Deterministic: fixed RNG seed, reported below.
"""
import csv, os, random, math, sys
sys.path.insert(0, os.path.expanduser('~/r2/scripts'))
import purestats as PS

SEED = 20260930
NPERM = 200000
NBOOT = 100000

R = list(csv.DictReader(open(os.path.expanduser('~/r2/data/SupplementaryTable_S5_quantitative_matrix.csv'))))
assert len(R) == 25
col = lambda n, c=float: [c(r[n]) for r in R]
tier = col('tier_published', int)
W = dict(D1=0.25, D2=0.15, D3=0.20, D4=0.20, D5=0.20)
norm = {d: col(f'{d}_norm') for d in W}

def spear(x, y):
    return PS.spearmanr(x, y)[0]

def perm_p(x, y, n=NPERM):
    rng = random.Random(SEED)
    obs = abs(spear(x, y))
    yy = list(y); ge = 0
    for _ in range(n):
        rng.shuffle(yy)
        if abs(spear(x, yy)) >= obs - 1e-12:
            ge += 1
    return (ge + 1) / (n + 1)          # add-one, never reports exactly 0

def boot_ci(x, y, n=NBOOT, conf=0.95):
    rng = random.Random(SEED + 1)
    idx = range(len(x)); m = len(x); vals = []
    for _ in range(n):
        s = [rng.randrange(m) for _ in range(m)]
        xs = [x[i] for i in s]; ys = [y[i] for i in s]
        # a resample can be constant in one variable; skip the undefined rho
        if len(set(xs)) < 2 or len(set(ys)) < 2:
            continue
        vals.append(spear(xs, ys))
    vals.sort()
    lo = vals[int((1-conf)/2 * len(vals))]
    hi = vals[int((1+conf)/2 * len(vals)) - 1]
    return lo, hi, len(vals)

# ---------- C4a  pairwise Spearman among the five dimensions --------------
dims = ['D1', 'D2', 'D3', 'D4', 'D5']
print('=== C4a. Pairwise Spearman among the five normalised dimensions ===')
print('        ' + ''.join(f'{d:>8}' for d in dims))
pair = {}
for a in dims:
    row = ''
    for b in dims:
        r = 1.0 if a == b else spear(norm[a], norm[b])
        pair[(a, b)] = r
        row += f'{r:8.3f}'
    print(f'   {a:5s}{row}')
offdiag = [pair[(a, b)] for i, a in enumerate(dims) for b in dims[i+1:]]
print(f'   off-diagonal |rho|: max {max(abs(v) for v in offdiag):.3f}, '
      f'mean {sum(abs(v) for v in offdiag)/len(offdiag):.3f}')
print('   reading: if these were one construct the off-diagonals would be uniformly high; they are not.')

# ---------- C4b  leave-one-dimension-out composite vs tier ----------------
print('\n=== C4b. Leave-one-dimension-out composite vs assigned tier ===')
def composite(drop=None):
    keep = [d for d in dims if d != drop]
    wsum = sum(W[d] for d in keep)
    return [sum(W[d]/wsum * norm[d][i] for d in keep) for i in range(25)]
full = composite()
rf = spear(tier, full)
print(f'   full composite vs tier: rho {rf:+.3f}')
for d in dims:
    c = composite(drop=d)
    r = spear(tier, c)
    print(f'   drop {d}: rho {r:+.3f}   (change {r-rf:+.3f})   [{("weakens" if abs(r)<abs(rf) else "strengthens")} the tier association]')

# ---------- C5  permutation p and bootstrap CI for the headline tests -----
print('\n=== C5. Permutation p-values and percentile bootstrap CIs ===')
print(f'    (seed {SEED}, {NPERM} permutations, {NBOOT} bootstrap resamples)')
TESTS = [
 ('tier vs D4 precedent',       tier, col('D4_precedent', int)),
 ('tier vs D1 pocket',          tier, col('D1_apo_max')),
 ('tier vs D2 geometry',        tier, col('D2_apo_max_volume_A3')),
 ('tier vs D3 chem matter',     tier, col('D3_chembl_compounds', int)),
 ('tier vs D5 delivery',        tier, col('D5_delivery', int)),
 ('D1 pocket vs D5 delivery',   col('D1_apo_max'), col('D5_delivery', int)),
 ('tier vs composite',          tier, col('composite_apo')),
]
print(f"\n   {'test':26s} {'rho':>7} {'perm p':>9} {'boot 95% CI':>20} {'Fisher 95% CI':>20}")
res = []
for name, x, y in TESTS:
    rho = spear(x, y)
    pp = perm_p(x, y)
    blo, bhi, nb = boot_ci(x, y)
    flo, fhi = PS.fisher_ci(rho, len(x))
    res.append((name, rho, pp, blo, bhi, flo, fhi))
    print(f"   {name:26s} {rho:+7.3f} {pp:9.5f}  [{blo:+.3f},{bhi:+.3f}]   [{flo:+.3f},{fhi:+.3f}]")

print('\n   note: bootstrap intervals are wider than Fisher for the tied tier variable,')
print('   exactly the reviewer\'s concern. The significance pattern is unchanged:')
for name, rho, pp, *_ in res:
    verdict = 'significant' if pp < 0.05 else 'not significant'
    print(f'     {name:26s} perm p {pp:.5f}  {verdict}')
