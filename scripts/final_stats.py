#!/usr/bin/env python3
"""Authoritative R2 statistics on matrix_v13, with permutation p and percentile
bootstrap CI for every headline test. Seed 20260930. Stdlib only."""
import csv, math, os, random, sys
sys.path.insert(0, os.path.expanduser('~/r2/scripts')); import purestats as PS
V={r['target_id']:r for r in csv.DictReader(open(os.path.expanduser('~/r2/out/matrix_v13.csv')))}
T=sorted(V); SEED=20260930; NP=100000; NB=50000
def f(x): return None if x=='NA' else float(x)
tier=[int(V[t]['tier_published']) for t in T]
comp=[float(V[t]['composite']) for t in T]
d4=[int(V[t]['D4_precedent']) for t in T]; d5=[int(V[t]['D5_delivery']) for t in T]
d3=[int(V[t]['D3_bioactivity_pchembl']) for t in T]
def sub(mask, xs): return [x for x,m in zip(xs,mask) if m]
site_mask=[V[t]['D1_site_druggability']!='NA' for t in T]
d1=[f(V[t]['D1_site_druggability']) for t in T]
d2=[f(V[t]['D2_norm']) for t in T]
def spear(x,y): return PS.spearmanr(x,y)[0]
def perm_p(x,y,n=NP):
    rng=random.Random(SEED); obs=abs(spear(x,y)); yy=list(y); ge=0
    for _ in range(n):
        rng.shuffle(yy)
        if abs(spear(x,yy))>=obs-1e-12: ge+=1
    return (ge+1)/(n+1)
def boot(x,y,n=NB,c=0.95):
    rng=random.Random(SEED+1); m=len(x); v=[]
    for _ in range(n):
        s=[rng.randrange(m) for _ in range(m)]; xs=[x[i] for i in s]; ys=[y[i] for i in s]
        if len(set(xs))<2 or len(set(ys))<2: continue
        v.append(spear(xs,ys))
    v.sort(); return v[int((1-c)/2*len(v))], v[int((1+c)/2*len(v))-1]
TESTS=[
 ("tier vs composite", tier, comp, [True]*len(T)),
 ("tier vs D4 precedent", tier, d4, [True]*len(T)),
 ("tier vs D5 delivery", tier, d5, [True]*len(T)),
 ("tier vs D3 bioactivity", tier, d3, [True]*len(T)),
 ("tier vs D1 site druggability", tier, d1, site_mask),
 ("D1 site vs D4 precedent (central null)", d1, d4, site_mask),
 ("D1 site vs D3 bioactivity", d1, d3, site_mask),
]
print(f"{'test':40s} {'n':>3} {'rho':>7} {'perm p':>9} {'boot 95% CI':>18}")
res=[]
for name,x,y,mask in TESTS:
    xs=sub(mask,x); ys=sub(mask,y); r=spear(xs,ys); pp=perm_p(xs,ys); lo,hi=boot(xs,ys)
    res.append((name,len(xs),r,pp,lo,hi))
    print(f"{name:40s} {len(xs):3d} {r:+7.3f} {pp:9.5f}  [{lo:+.3f},{hi:+.3f}]")
# BH across the family
ps=[p for *_,p,_,_ in [(a,b,c,d,e,g) for (a,b,c,d,e,g) in [(n,nn,r,pp,lo,hi) for (n,nn,r,pp,lo,hi) in res]]]
ps=[r[3] for r in res]; m=len(ps); order=sorted(range(m),key=lambda i:-ps[i]); prev=1.0; q={}
for pos,i in enumerate(order):
    q[i]=min(prev, ps[i]*m/(m-pos)); prev=q[i]
print("\nBenjamini-Hochberg q across the 7-test family:")
for i,(name,n,r,pp,lo,hi) in enumerate(res):
    print(f"  {name:40s} q {q[i]:.4f}  {'survives' if q[i]<0.05 else 'ns'}")
import json
json.dump([dict(test=n,n=nn,rho=round(r,3),perm_p=round(pp,5),ci=[round(lo,3),round(hi,3)],bh_q=round(q[i],4))
           for i,(n,nn,r,pp,lo,hi) in enumerate(res)],
          open(os.path.expanduser('~/r2/out/final_stats.json'),'w'), indent=1)
print("\nwrote out/final_stats.json")
