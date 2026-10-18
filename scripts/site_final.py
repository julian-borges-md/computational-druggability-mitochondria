#!/usr/bin/env python3
import csv, os, sys
sys.path.insert(0, os.path.expanduser('~/r2/scripts')); import purestats as PS
ISO={r['target_id']:r for r in csv.DictReader(open(os.path.expanduser('~/r2/out/site_matched_isolated.csv')))}
SUB={r['target_id']:r for r in csv.DictReader(open(os.path.expanduser('~/r2/out/site_matched_subcomplex.csv')))}
M  ={r['target_id']:r for r in csv.DictReader(open(os.path.expanduser('~/r2/data/SupplementaryTable_S5_quantitative_matrix.csv')))}
have=[t for t in ISO if ISO[t]['tier'] in ('D','U')]
def f(x):
    try: return float(x)
    except: return None
tier=lambda t:int(M[t]['tier_published']); d4=lambda t:int(M[t]['D4_precedent'])
iso=lambda t:f(ISO[t]['site_matched_D1']); sub=lambda t:f(SUB[t]['site_matched_D1_sub'])
glob=lambda t:f(ISO[t]['global_best_D1'])
best=lambda t:max(iso(t),sub(t))            # most complete available context
def show(name,xs,ys):
    r,p=PS.spearmanr(xs,ys); lo,hi=PS.fisher_ci(r,len(xs))
    print(f"  {name:46s} n={len(xs):2d}  rho {r:+.3f}  p {p:.4f}  CI [{lo:+.3f},{hi:+.3f}]")
print("Structure vs pharmacological precedent (the central independence claim):")
show("  global-best pocket, all 25", [d4(t) for t in M], [f(M[t]['D1_apo_max']) for t in M])
show("  site-matched, isolated chain", [d4(t) for t in have], [iso(t) for t in have])
show("  site-matched, assembly context", [d4(t) for t in have], [sub(t) for t in have])
show("  site-matched, best available context", [d4(t) for t in have], [best(t) for t in have])
print("\nStructure vs assigned tier (descriptive, partly circular per comment 4):")
show("  site-matched, best available context", [tier(t) for t in have], [best(t) for t in have])
print("\nPer-target site-matched D1 (structure druggability at the therapeutic site):")
print(f"  {'gene':10}{'tier':5}{'iso':>7}{'assembly':>10}{'best':>7}{'global':>8}")
for t in sorted(have, key=lambda t:-best(t)):
    print(f"  {ISO[t]['gene']:10}{ISO[t]['tier']:5}{iso(t):7.3f}{sub(t):10.3f}{best(t):7.3f}{glob(t):8.3f}")
# write combined site-matched table for the matrix rebuild
rows=[]
for t in sorted(ISO):
    tr=ISO[t]['tier']
    rows.append(dict(target_id=t, gene=ISO[t]['gene'], site_tier=tr,
        site_matched_D1_isolated=ISO[t]['site_matched_D1'],
        site_matched_D1_assembly=SUB[t]['site_matched_D1_sub'],
        site_matched_D1_best=('NA' if tr=='N' else round(best(t),3)),
        global_best_D1=ISO[t]['global_best_D1'],
        n_site_res=ISO[t]['n_site_res'], overlap_iso=ISO[t]['site_overlap'],
        interface_pocket=SUB[t].get('interface_pocket','')))
with open(os.path.expanduser('~/r2/out/site_matched_combined.csv'),'w',newline='') as fh:
    w=csv.DictWriter(fh,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print('\nwrote out/site_matched_combined.csv')
