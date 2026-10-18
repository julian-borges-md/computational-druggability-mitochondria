#!/usr/bin/env python3
"""AlphaFold local confidence for the 3 modelled targets (comment 8).
pLDDT = B-factor column of the model. Report global mean, and mean at the
a priori site residues, plus mean PAE. Justifies treating structure as NA where
the model or the site region is unreliable."""
import csv, json, os, gemmi, statistics
D=os.path.expanduser('~/r2/struct')
S2={r['target_id']:r for r in csv.DictReader(open(os.path.expanduser('~/r2/data/SupplementaryTable_S2_fpocket_per_target.csv')))}
UP=json.load(open(os.path.expanduser('~/r2/out/uniprot_sites.json')))
AF={'TGT_014':('TIMM23','AF-O14925','O14925'),'TGT_020':('PPARGC1A','AF-Q9UBK2','Q9UBK2'),
    'TGT_021':('SLC25A4','AF-P12235','P12235')}
def plddt(struct):
    st=gemmi.read_structure(f'{D}/{struct}.pdb'); m=st[0]
    per={}
    for ch in m:
        for res in ch:
            bs=[a.b_iso for a in res]
            if bs: per[res.seqid.num]=sum(bs)/len(bs)
    return per
rows=[]
for tid,(gene,struct,acc) in AF.items():
    per=plddt(struct)
    vals=list(per.values()); gm=statistics.mean(vals)
    upres=sorted({f['start'] for f in UP[tid]['features'] if f['start']})
    site=[per[u] for u in upres if u in per]
    # PAE mean
    pae=json.load(open(f'{D}/{struct}_pae.json'))
    pae_arr = pae[0]['predicted_aligned_error'] if isinstance(pae,list) else pae.get('predicted_aligned_error',[])
    pae_mean = statistics.mean(statistics.mean(row) for row in pae_arr) if pae_arr else float('nan')
    frac_conf = sum(1 for v in vals if v>=70)/len(vals)
    site_mean = statistics.mean(site) if site else None
    rows.append(dict(target_id=tid, gene=gene, model=struct, uniprot=acc,
        global_pLDDT=round(gm,1), frac_pLDDT_ge70=round(frac_conf,2),
        site_pLDDT=('NA' if site_mean is None else round(site_mean,1)),
        n_site_res=len(site), mean_PAE_A=round(pae_mean,1),
        interpretation=('retained: high confidence at site' if gm>=70 else 'treated as NA: low confidence / disordered')))
    print(f"  {tid} {gene:9} global pLDDT {gm:5.1f}  frac>=70 {frac_conf:4.0%}  "
          f"site pLDDT {statistics.mean(site):5.1f} (n={len(site)})  mean PAE {pae_mean:4.1f} A" if site
          else f"  {tid} {gene:9} global pLDDT {gm:5.1f}  frac>=70 {frac_conf:4.0%}  no UniProt site res  mean PAE {pae_mean:4.1f} A")
import csv as _csv
with open(os.path.expanduser('~/r2/out/af_confidence.csv'),'w',newline='') as fh:
    w=_csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("wrote out/af_confidence.csv")
