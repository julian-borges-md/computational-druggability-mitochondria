#!/usr/bin/env python3
"""Matrix v1.3 for R2. Structural dimensions D1/D2 are the site-matched, best-
available-context scores; NA for the 8 targets with no defined small-molecule
site (comment 8). D3 renamed bioactivity evidence, from pChEMBL counts (comments
1, 7). D4, D5 unchanged. Composite = weight-renormalised mean over available
dimensions. Normalisation equations are explicit (comment 6):
  D1 = Druggability Score, native theoretical range [0,1].
  D2 = (log10 V - min)/(max - min) over targets with a site (observed range).
  D3 = (log10(n+1) - min)/(max - min) over all 25 (observed range).
  D4 = precedent / 3 (theoretical range 0-3).
  D5 = delivery  / 3 (theoretical range 0-3).
Weights: D1 0.25, D2 0.15, D3 0.20, D4 0.20, D5 0.20.
"""
import csv, math, os, sys
sys.path.insert(0, os.path.expanduser('~/r2/scripts')); import purestats as PS
M ={r['target_id']:r for r in csv.DictReader(open(os.path.expanduser('~/r2/data/SupplementaryTable_S5_quantitative_matrix.csv')))}
SS={r['target_id']:r for r in csv.DictReader(open(os.path.expanduser('~/r2/out/site_scores_final.csv')))}
CH={r['target_id']:r for r in csv.DictReader(open(os.path.expanduser('~/r2/out/chembl_filtered.csv')))}
W={'D1':0.25,'D2':0.15,'D3':0.20,'D4':0.20,'D5':0.20}
T=sorted(M)
def num(x):
    try: return float(x)
    except: return None
# raw dimension values
d1={t:(None if SS[t]['site_D1_best']=='NA' else float(SS[t]['site_D1_best'])) for t in T}
d2={t:(None if SS[t]['site_D2_best']=='NA' else float(SS[t]['site_D2_best'])) for t in T}
d3={t:int(CH[t]['n_pchembl']) for t in T}
d4={t:int(M[t]['D4_precedent']) for t in T}
d5={t:int(M[t]['D5_delivery']) for t in T}
tier={t:int(M[t]['tier_published']) for t in T}
def mm(vals):
    v=[x for x in vals if x is not None]; lo,hi=min(v),max(v)
    return lo,hi
# normalisations
n1={t:(None if d1[t] is None else d1[t]) for t in T}                       # native [0,1]
l2={t:(None if d2[t] is None else math.log10(d2[t]+1)) for t in T}; lo2,hi2=mm(l2.values())
n2={t:(None if l2[t] is None else (l2[t]-lo2)/(hi2-lo2)) for t in T}
l3={t:math.log10(d3[t]+1) for t in T}; lo3,hi3=mm(l3.values())
n3={t:(l3[t]-lo3)/(hi3-lo3) for t in T}
n4={t:d4[t]/3 for t in T}; n5={t:d5[t]/3 for t in T}
NORM={'D1':n1,'D2':n2,'D3':n3,'D4':n4,'D5':n5}
def composite(t):
    num_=0.0; den=0.0
    for k in W:
        v=NORM[k][t]
        if v is not None: num_+=W[k]*v; den+=W[k]
    return num_/den
comp={t:composite(t) for t in T}
order=sorted(T, key=lambda t:-comp[t]); rank={t:i+1 for i,t in enumerate(order)}
# write matrix
rows=[]
for t in order:
    rows.append(dict(rank=rank[t], target_id=t, gene=M[t]['gene'], assessed_protein=M[t]['assessed_protein'],
        uniprot=M[t]['uniprot'], tier_published=tier[t], site_tier=SS[t]['tier'],
        composite=round(comp[t],4),
        D1_site_druggability=('NA' if d1[t] is None else round(d1[t],3)),
        D2_site_volume_A3=('NA' if d2[t] is None else round(d2[t],1)),
        D3_bioactivity_pchembl=d3[t], D3_bioactivity_current=int(CH[t]['n_current']),
        D3_bioactivity_binding=int(CH[t]['n_pchembl_binding']),
        D4_precedent=d4[t], D5_delivery=d5[t],
        D1_norm=('NA' if n1[t] is None else round(n1[t],4)),
        D2_norm=('NA' if n2[t] is None else round(n2[t],4)),
        D3_norm=round(n3[t],4), D4_norm=round(n4[t],4), D5_norm=round(n5[t],4),
        global_best_D1=SS[t].get('site_D1_iso',''),  # placeholder; exploratory col added elsewhere
        chembl_status=CH[t]['status']))
with open(os.path.expanduser('~/r2/out/matrix_v13.csv'),'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
# ---- correlations on v1.3 ----
def sp(xs,ys):
    r,p=PS.spearmanr(xs,ys); lo,hi=PS.fisher_ci(r,len(xs)); return r,p,lo,hi
have=[t for t in T if d1[t] is not None]
print("=== v1.3 correlations ===")
for name,xf,yf,ts in [
  ("tier vs composite", tier, comp, T),
  ("tier vs D4 precedent", tier, d4, T),
  ("tier vs D5 delivery", tier, d5, T),
  ("tier vs D3 bioactivity", tier, d3, T),
  ("tier vs D1 site druggability", tier, d1, have),
  ("D1 site vs D4 precedent (central)", d1, d4, have),
  ("D1 site vs D3 bioactivity", d1, d3, have)]:
    xs=[xf[t] for t in ts]; ys=[yf[t] for t in ts]; r,p,lo,hi=sp(xs,ys)
    print(f"  {name:38s} n={len(ts):2d}  rho {r:+.3f}  p {p:.4f}  CI [{lo:+.3f},{hi:+.3f}]")
print("\n=== first wave: composite top 8 ===")
for t in order[:8]:
    print(f"  {rank[t]}. {M[t]['gene']:9} tier{tier[t]} comp {comp[t]:.3f}  "
          f"D1 {str(rows[rank[t]-1]['D1_site_druggability']):>5} D3 {d3[t]:>3} D4 {d4[t]} D5 {d5[t]}  site {SS[t]['tier']}")
print('\nwrote out/matrix_v13.csv')
