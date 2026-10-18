#!/usr/bin/env python3
import csv, os, sys
sys.path.insert(0, os.path.expanduser('~/r2/scripts')); import purestats as PS
SM = {r['target_id']: r for r in csv.DictReader(open(os.path.expanduser('~/r2/out/site_matched_isolated.csv')))}
M  = {r['target_id']: r for r in csv.DictReader(open(os.path.expanduser('~/r2/data/SupplementaryTable_S5_quantitative_matrix.csv')))}
have = [t for t in SM if SM[t]['tier'] in ('D','U')]           # targets with a defined site
def col(ts, f): return [f(t) for t in ts]
tier = lambda t: int(M[t]['tier_published'])
d4   = lambda t: int(M[t]['D4_precedent'])
smd1 = lambda t: float(SM[t]['site_matched_D1'])
gbd1 = lambda t: float(SM[t]['global_best_D1'])
def show(name, x, y, ts):
    r, p = PS.spearmanr(x, y); lo, hi = PS.fisher_ci(r, len(x))
    print(f"  {name:42s} n={len(ts):2d}  rho {r:+.3f}  p {p:.4f}  CI [{lo:+.3f},{hi:+.3f}]")
print("=== all 25, global-best pocket (current manuscript) ===")
allt = list(M)
show("tier vs global-best D1", col(allt, lambda t: tier(t)), col(allt, lambda t: float(M[t]['D1_apo_max'])), allt)
show("precedent vs global-best D1", col(allt, lambda t: d4(t)), col(allt, lambda t: float(M[t]['D1_apo_max'])), allt)
print(f"\n=== {len(have)} targets with a defined site, isolated chain ===")
show("tier vs site-matched D1", col(have, tier), col(have, smd1), have)
show("precedent vs site-matched D1", col(have, d4), col(have, smd1), have)
show("tier vs global-best D1 (same 17)", col(have, tier), col(have, gbd1), have)
show("precedent vs global-best D1 (same 17)", col(have, d4), col(have, gbd1), have)
print("\n=== site-matched vs global-best agreement (how often the global best IS the site) ===")
same = sum(1 for t in have if abs(smd1(t)-gbd1(t)) < 0.05)
print(f"  global-best pocket coincides with the site (|diff|<0.05) in {same}/{len(have)} targets")
gap = [(SM[t]['gene'], gbd1(t)-smd1(t)) for t in have]
gap.sort(key=lambda x: -x[1])
print("  largest global-minus-site gaps:", ', '.join(f'{g} {d:+.2f}' for g,d in gap[:6]))
