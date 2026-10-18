#!/usr/bin/env python3
"""Derived statistics for the R2 tables and response letter, all from
matrix_v13.csv + site_scores_final.csv, stdlib only. Produces:
  (A) structure vs precedent under the isolated and assembly specifications
  (B) full 5x5 pairwise Spearman matrix among the dimensions (pairwise complete)
  (C) leave one dimension out: composite vs tier, dropping each dimension
Writes out/derived_stats.json. Seed-free (all analytic Spearman)."""
import csv, math, os, sys, json
sys.path.insert(0, os.path.expanduser('~/r2/scripts')); import purestats as PS

V ={r['target_id']:r for r in csv.DictReader(open(os.path.expanduser('~/r2/out/matrix_v13.csv')))}
SS={r['target_id']:r for r in csv.DictReader(open(os.path.expanduser('~/r2/out/site_scores_final.csv')))}
T=sorted(V)
def f(x): return None if x in ('NA','') else float(x)

tier ={t:int(V[t]['tier_published']) for t in T}
d4   ={t:int(V[t]['D4_precedent']) for t in T}
site_iso ={t:f(SS[t]['site_D1_iso'])  for t in T}   # isolated-chain site score
site_best={t:f(SS[t]['site_D1_best']) for t in T}   # best-context (assembly) site score
have=[t for t in T if site_best[t] is not None]      # 17 with a defined site

def sp(xs,ys):
    r,p=PS.spearmanr(xs,ys); lo,hi=PS.fisher_ci(r,len(xs)); return r,p,lo,hi

out={}

# (A) structure vs precedent, isolated vs assembly
print("=== (A) structure vs pharmacological precedent (D4) ===")
specs={}
for name,src in [("isolated", site_iso), ("assembly_best", site_best)]:
    xs=[src[t] for t in have]; ys=[d4[t] for t in have]
    r,p,lo,hi=sp(xs,ys)
    specs[name]=dict(n=len(have), rho=round(r,3), p=round(p,4), ci=[round(lo,3),round(hi,3)])
    print(f"  {name:14s} n={len(have)} rho {r:+.3f} p {p:.4f} CI [{lo:+.3f},{hi:+.3f}]")
out['structure_vs_precedent']=specs

# (B) pairwise 5x5 among dimensions, pairwise-complete on the normalised columns
print("\n=== (B) pairwise Spearman among dimensions (pairwise complete) ===")
dims=['D1','D2','D3','D4','D5']
norm={k:{t:f(V[t][f'{k}_norm']) for t in T} for k in dims}
mat={}; offdiag=[]
for i,a in enumerate(dims):
    for b in dims[i+1:]:
        pairs=[(norm[a][t],norm[b][t]) for t in T if norm[a][t] is not None and norm[b][t] is not None]
        xs=[p[0] for p in pairs]; ys=[p[1] for p in pairs]
        r=PS.spearmanr(xs,ys)[0]
        mat[f"{a}-{b}"]=dict(n=len(pairs), rho=round(r,3))
        offdiag.append(r)
        print(f"  {a}-{b}  n={len(pairs):2d}  rho {r:+.3f}")
out['pairwise']=mat
out['pairwise_summary']=dict(
    min=round(min(offdiag),3), max=round(max(offdiag),3),
    mean_abs=round(sum(abs(x) for x in offdiag)/len(offdiag),3))
print(f"  off-diagonal: min {min(offdiag):+.3f}  max {max(offdiag):+.3f}  mean|.| {sum(abs(x) for x in offdiag)/len(offdiag):.3f}")

# (C) leave-one-dimension-out composite vs tier
print("\n=== (C) leave-one-dimension-out: composite vs tier (n=25) ===")
W={'D1':0.25,'D2':0.15,'D3':0.20,'D4':0.20,'D5':0.20}
def composite(t, drop=None):
    num=den=0.0
    for k in W:
        if k==drop: continue
        v=norm[k][t]
        if v is not None: num+=W[k]*v; den+=W[k]
    return num/den if den else 0.0
base_comp=[composite(t) for t in T]; tiers=[tier[t] for t in T]
r0=PS.spearmanr(tiers, base_comp)[0]
print(f"  full composite            rho {r0:+.3f}")
loo={'full':round(r0,3)}
for d in dims:
    c=[composite(t, drop=d) for t in T]
    r=PS.spearmanr(tiers, c)[0]
    loo[f"drop_{d}"]=round(r,3)
    print(f"  drop {d}                    rho {r:+.3f}")
out['leave_one_out']=loo

json.dump(out, open(os.path.expanduser('~/r2/out/derived_stats.json'),'w'), indent=1)
print("\nwrote out/derived_stats.json")
