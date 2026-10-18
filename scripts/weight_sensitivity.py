#!/usr/bin/env python3
"""Weight-robustness / first-wave rank stability (comment 6).
Base weights D1 .25 D2 .15 D3 .20 D4 .20 D5 .20. Sample 20000 weight vectors,
each dimension ~ Uniform(0.5x, 1.5x) base then renormalised, recompute the
available-case composite, and record how often each target is in the top 5 of
the first-wave-eligible set (delivery D5 >= 2, i.e. a viable small-molecule
route). Also equal-weights and precedent-free schemes as named alternatives."""
import csv, math, os, random, sys
sys.path.insert(0, os.path.expanduser('~/r2/scripts')); import purestats as PS
M ={r['target_id']:r for r in csv.DictReader(open(os.path.expanduser('~/r2/data/SupplementaryTable_S5_quantitative_matrix.csv')))}
V ={r['target_id']:r for r in csv.DictReader(open(os.path.expanduser('~/r2/out/matrix_v13.csv')))}
T=sorted(M)
def f(x): 
    return None if x=='NA' else float(x)
NORM={t:{k:f(V[t][f'{k}_norm']) for k in ('D1','D2','D3','D4','D5')} for t in T}
d5={t:int(V[t]['D5_delivery']) for t in T}
gene={t:V[t]['gene'] for t in T}
base={'D1':0.25,'D2':0.15,'D3':0.20,'D4':0.20,'D5':0.20}
def comp(t,w):
    num=den=0.0
    for k in w:
        v=NORM[t][k]
        if v is not None: num+=w[k]*v; den+=w[k]
    return num/den if den else 0.0
def first_wave(w, k=5):
    elig=[t for t in T if d5[t]>=2]
    return set(sorted(elig, key=lambda t:-comp(t,w))[:k])
def base_rank(w):
    return sorted(T, key=lambda t:-comp(t,w))
rng=random.Random(20260930); N=20000
freq={t:0 for t in T}; rankcorr=[]
b_order=base_rank(base); b_idx={t:i for i,t in enumerate(b_order)}
for _ in range(N):
    w={k:base[k]*rng.uniform(0.5,1.5) for k in base}; s=sum(w.values()); w={k:v/s for k,v in w.items()}
    for t in first_wave(w): freq[t]+=1
    o=base_rank(w); rankcorr.append(PS.spearmanr([b_idx[t] for t in T],[o.index(t) for t in T])[0])
print(f"First-wave (top-5 of delivery-eligible) membership over {N} random weightings:")
for t in sorted(T, key=lambda t:-freq[t]):
    if freq[t]==0: continue
    print(f"  {gene[t]:9} {100*freq[t]/N:5.1f}%")
print(f"\nMean Spearman rank correlation vs base weighting: {sum(rankcorr)/len(rankcorr):.3f}")
# named alternative schemes
print("\nFirst wave under named alternative schemes (delivery-eligible top 5):")
schemes={'base':base,
 'equal (0.2 each)':{k:0.2 for k in base},
 'structure-heavy (D1 .40)':{'D1':0.40,'D2':0.20,'D3':0.15,'D4':0.15,'D5':0.10},
 'evidence-heavy (D3+D4 .30 each)':{'D1':0.10,'D2':0.10,'D3':0.30,'D4':0.30,'D5':0.20},
 'precedent-free (D4=0)':{'D1':0.30,'D2':0.20,'D3':0.25,'D4':0.0,'D5':0.25}}
scheme_rows=[]
for nm,w in schemes.items():
    s=sum(w.values()); w={k:v/s for k,v in w.items()}
    fw=sorted([t for t in T if d5[t]>=2], key=lambda t:-comp(t,w))[:5]
    print(f"  {nm:32s}: {', '.join(gene[t] for t in fw)}")
    scheme_rows.append(dict(scheme=nm, first_wave='; '.join(gene[t] for t in fw)))
# dump machine-readable outputs
import json
freq_rows=[dict(gene=gene[t], target_id=t, pct_top5=round(100*freq[t]/N,1))
           for t in sorted(T, key=lambda t:-freq[t]) if freq[t]>0]
json.dump(dict(n_draws=N, mean_rank_corr=round(sum(rankcorr)/len(rankcorr),3),
               membership=freq_rows, schemes=scheme_rows),
          open(os.path.expanduser('~/r2/out/weight_sensitivity.json'),'w'), indent=1)
with open(os.path.expanduser('~/r2/out/weight_sensitivity.csv'),'w',newline='') as fh:
    w2=csv.DictWriter(fh, fieldnames=['gene','target_id','pct_top5']); w2.writeheader(); w2.writerows(freq_rows)
print("\nwrote out/weight_sensitivity.json and .csv")
