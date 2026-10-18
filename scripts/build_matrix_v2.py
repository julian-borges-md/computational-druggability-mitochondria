#!/usr/bin/env python3
"""Corrected 25 x 5 quantitative matrix, v2.
Changes from v1: MT-ND1 contamination removed; TFAM, SOD2, ATP5F1A structures
replaced; MT-CYB and KEAP1 precedent scores corrected against the manuscript's
own rule; apo, holo and rank-1 pocket specifications all carried."""
import csv, math, json
v2 ={r['target_id']:r for r in csv.DictReader(open('out/TableS2_fpocket_v2.csv'))}
old={r['target_id']:r for r in csv.DictReader(open('out/TableS6_quantitative_matrix_25x5.csv'))}
ch ={r['target_id']:r for r in csv.DictReader(open('out/chembl_ligand_evidence.csv'))}
xw ={r['canonical_id']:r for r in csv.DictReader(open('out/id_crosswalk_v1.1.csv'))}

PREC={t:int(old[t]['D4_precedent_score']) for t in old}
NOTE={t:old[t]['D4_precedent_basis'] for t in old}
PREC['TGT_003']=0; NOTE['TGT_003']='none at the assessed site; atovaquone binds Plasmodium cytochrome b, not the human protein'
PREC['TGT_024']=1; NOTE['TGT_024']='validated tool or preclinical; omaveloxolone acts at BTB Cys151, not the assessed Kelch domain'
DEL={'A_no_targeting_required':3,'B_matrix_MTP_or_TPP':2,'C_IMS_or_membrane_embedded':1,
     'D_gene_therapy_primary':0,'E_no_viable_small_molecule':0}
T=sorted(v2)
def nrm(v):
    lo,hi=min(v),max(v); return [0.0 if hi==lo else (x-lo)/(hi-lo) for x in v]
d1=[float(v2[t]['apo_max_drugscore']) for t in T]
d1h=[float(v2[t]['holo_max_drugscore']) for t in T]
d1r=[float(v2[t]['apo_rank1_drugscore']) for t in T]
d2=[math.log10(float(v2[t]['apo_max_volume'])+1) for t in T]
d3=[math.log10(int(ch[t]['n_distinct_compounds'])+1) for t in T]
d4=[PREC[t] for t in T]; d5=[DEL[xw[t]['delivery_category']] for t in T]
n1,n1h,n1r,n2,n3,n4,n5=map(nrm,(d1,d1h,d1r,d2,d3,d4,d5))
W=(0.25,0.15,0.20,0.20,0.20)
comp=[W[0]*a+W[1]*b+W[2]*c+W[3]*d+W[4]*e for a,b,c,d,e in zip(n1,n2,n3,n4,n5)]
comph=[W[0]*a+W[1]*b+W[2]*c+W[3]*d+W[4]*e for a,b,c,d,e in zip(n1h,n2,n3,n4,n5)]
compr=[W[0]*a+W[1]*b+W[2]*c+W[3]*d+W[4]*e for a,b,c,d,e in zip(n1r,n2,n3,n4,n5)]
rows=[]
for i,t in enumerate(T):
    r=v2[t]; ncomp=int(ch[t]['n_distinct_compounds'])
    obs='STRONG' if ncomp>10 else ('MODERATE' if ncomp>=1 else 'ABSENT')
    rows.append(dict(target_id=t, gene=r['gene'], assessed_protein=r['assessed_protein'], uniprot=r['uniprot'],
        tier_published=int(r['tier']), structure=r['structure'], resolution_A=r['resolution_A'],
        organism=r['organism'], mutation=r['mutation'] or 'none', fragment=r['fragment'] or 'full construct',
        n_chains_in_entity=r['n_chains_in_entity'], n_chains_retained=r['n_chains_retained'],
        n_heteroatoms_removed=r['n_heteroatoms_removed'],
        D1_apo_max=round(d1[i],3), D1_holo_max=round(d1h[i],3), D1_rank1=round(d1r[i],3),
        D1_apo_minus_holo=round(d1[i]-d1h[i],3),
        D2_apo_max_volume_A3=round(float(r['apo_max_volume']),1), n_pockets=int(r['apo_n_pockets']),
        D3_chembl_compounds=ncomp, D3_chembl_target=ch[t]['chembl_id'] or 'no human entity',
        D3_class_observed=obs, D3_class_published=old[t]['D3_class_published'],
        D4_precedent=PREC[t], D4_basis=NOTE[t], D4_precedent_v1=int(old[t]['D4_precedent_score']),
        D5_delivery=DEL[xw[t]['delivery_category']], D5_category=xw[t]['delivery_category'],
        D1_norm=round(n1[i],4), D2_norm=round(n2[i],4), D3_norm=round(n3[i],4),
        D4_norm=round(n4[i],4), D5_norm=round(n5[i],4),
        composite_apo=round(comp[i],4), composite_holo=round(comph[i],4), composite_rank1=round(compr[i],4)))
rows.sort(key=lambda x:-x['composite_apo'])
for i,r in enumerate(rows,1): r['rank_apo']=i
for i,r in enumerate(sorted(rows,key=lambda x:-x['composite_holo']),1): r['rank_holo']=i
for i,r in enumerate(sorted(rows,key=lambda x:-x['composite_rank1']),1): r['rank_rank1']=i
rows.sort(key=lambda x:x['rank_apo'])
cols=['rank_apo','rank_holo','rank_rank1','target_id','gene','assessed_protein','uniprot','tier_published',
 'composite_apo','composite_holo','composite_rank1','D1_apo_max','D1_holo_max','D1_rank1','D1_apo_minus_holo',
 'D1_norm','D2_apo_max_volume_A3','D2_norm','n_pockets','D3_chembl_compounds','D3_chembl_target',
 'D3_class_published','D3_class_observed','D3_norm','D4_precedent','D4_precedent_v1','D4_basis','D4_norm',
 'D5_delivery','D5_category','D5_norm','structure','resolution_A','organism','mutation','fragment',
 'n_chains_in_entity','n_chains_retained','n_heteroatoms_removed']
with open('out/TableS5_quantitative_matrix_v2.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=cols,extrasaction='ignore'); w.writeheader(); w.writerows(rows)
print(f"{'rk':3s} {'gene':10s} {'T':2s} {'comp':6s} {'D1apo':6s} {'D1holo':6s} {'D1r1':6s} {'cmpd':5s} {'D4':2s} {'D5':2s} {'v1 rk':6s}")
oldrank={r['gene']:int(r['composite_rank']) for r in old.values()}
for r in rows:
    ch_=f"{oldrank.get(r['gene'],'-'):>5}"
    print(f"{r['rank_apo']:<3d} {r['gene']:10s} {r['tier_published']:<2d} {r['composite_apo']:<6.3f} "
          f"{r['D1_apo_max']:<6.2f} {r['D1_holo_max']:<6.2f} {r['D1_rank1']:<6.2f} {r['D3_chembl_compounds']:<5d} "
          f"{r['D4_precedent']:<2d} {r['D5_delivery']:<2d} {ch_}")
print("\nwrote out/TableS5_quantitative_matrix_v2.csv")
