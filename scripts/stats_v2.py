#!/usr/bin/env python3
"""Recompute every reported statistic from out/TableS5_quantitative_matrix_v2.csv.
Single source of truth for the registry, the manuscript and the verification pass."""
import csv, math, json, os
from scipy import stats

BASE = os.path.expanduser('~/csbrd')
R = list(csv.DictReader(open(f'{BASE}/out/TableS5_quantitative_matrix_v2.csv')))
assert len(R) == 25, len(R)

def col(n, cast=float): return [cast(r[n]) for r in R]

tier  = col('tier_published', int)
d1a   = col('D1_apo_max');  d1h = col('D1_holo_max'); d1r = col('D1_rank1')
d2    = col('D2_apo_max_volume_A3')
d3    = col('D3_chembl_compounds', int)
d4    = col('D4_precedent', int)
d5    = col('D5_delivery', int)
comp  = col('composite_apo')

def sp(x, y):
    rho, p = stats.spearmanr(x, y)
    n = len(x)
    z = math.atanh(rho); se = 1/math.sqrt(n-3)
    lo, hi = math.tanh(z-1.96*se), math.tanh(z+1.96*se)
    return round(rho,3), round(p,4), (round(lo,3), round(hi,3))

TESTS = [
    ('tier vs composite',                 tier, comp),
    ('tier vs precedent (D4)',            tier, d4),
    ('tier vs delivery (D5)',             tier, d5),
    ('tier vs ligand efficiency (D3)',    tier, d3),
    ('D1 apo vs D5 delivery',             d1a,  d5),
    ('tier vs pocket druggability (D1)',  tier, d1a),
    ('tier vs pocket geometry (D2)',      tier, d2),
    ('D1 apo vs precedent (D4)',          d1a,  d4),
]
res = [(n,) + sp(x, y) for n, x, y in TESTS]

# Benjamini-Hochberg, largest p to smallest
order = sorted(range(len(res)), key=lambda i: -res[i][2])
m = len(res); q = {}; prev = 1.0
for pos, i in enumerate(order):
    rank = m - pos
    v = min(prev, res[i][2] * m / rank)
    q[i] = round(v, 4); prev = v

print(f"{'test':36s} {'rho':>7s} {'p':>8s} {'BH q':>8s}  {'95% CI':>18s}  surv")
for i, (n, rho, p, ci) in enumerate(res):
    print(f"{n:36s} {rho:7.3f} {p:8.4f} {q[i]:8.4f}  [{ci[0]:6.3f},{ci[1]:6.3f}]  {'yes' if q[i]<0.05 else 'no'}")

print("\n-- D1 specification sensitivity --")
for lab, v in (('apo max', d1a), ('holo max', d1h), ('apo rank 1', d1r)):
    r1, p1, c1 = sp(tier, v); r2, p2, c2 = sp(v, d5)
    print(f"  {lab:11s}  tier~D1 rho={r1:6.3f} p={p1:.4f}   D1~D5 rho={r2:6.3f} p={p2:.4f}")

# MDES at n=25
def mdes(n, alpha=0.05, power=0.80):
    za = stats.norm.ppf(1-alpha/2); zb = stats.norm.ppf(power)
    return round(math.tanh((za+zb)/math.sqrt(n-3)), 3)
print(f"\nMDES n=25 alpha .05 power .80 : |rho| = {mdes(25)}")
print(f"MDES after Bonferroni (a=.05/8): |rho| = {mdes(25, 0.05/8)}")
print("power at true rho:", {r: round(float(stats.norm.cdf(math.atanh(r)*math.sqrt(22) - stats.norm.ppf(0.975))),2) for r in (0.2,0.3,0.4,0.5)})

t1 = [r for r in R if int(r['tier_published']) == 1]
print(f"\ntier counts: 1={sum(1 for t in tier if t==1)} 2={sum(1 for t in tier if t==2)} 3={sum(1 for t in tier if t==3)} total={len(tier)}")
c1 = [float(r['composite_apo']) for r in t1]
print(f"Tier 1 composite: min={min(c1):.3f} max={max(c1):.3f} sd={statistics_sd(c1) if False else round(__import__('statistics').stdev(c1),3)}")
print(f"Tier 1 ranks: {sorted(int(r['rank_apo']) for r in t1)}")
print(f"reclassified by ChEMBL: {sum(1 for r in R if r['D3_class_published']!=r['D3_class_observed'])} of 25")
print(f"STRONG published={sum(1 for r in R if r['D3_class_published']=='STRONG')} observed={sum(1 for r in R if r['D3_class_observed']=='STRONG')}")
print(f"top 5 by composite_apo: {[r['gene'] for r in R[:5]]}")

json.dump({'tests':[{'test':n,'rho':rho,'p':p,'ci':ci,'bh_q':q[i]} for i,(n,rho,p,ci) in enumerate(res)],
           'mdes_25':mdes(25),'mdes_bonf':mdes(25,0.05/8)},
          open(f'{BASE}/out_v2/stats_v2.json','w'), indent=1)
