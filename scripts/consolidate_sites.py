#!/usr/bin/env python3
"""Consolidate site-matched structural scores: for each target, the site-matched
pocket's Druggability Score AND volume, isolated + assembly arms, reading the
cached fpocket outputs. Adds volume (D2) to the combined table."""
import csv, json, os, re, gemmi
D=os.path.expanduser('~/r2/struct'); FP=os.path.expanduser('~/r2/fp')
S2={r['target_id']:r for r in csv.DictReader(open(os.path.expanduser('~/r2/data/SupplementaryTable_S2_fpocket_per_target.csv')))}
UP=json.load(open(os.path.expanduser('~/r2/out/uniprot_sites.json')))
ANCHOR={'TGT_001':['PLX'],'TGT_003':['HEM','U10'],'TGT_004':['HEA','CU'],'TGT_006':['FMN'],
        'TGT_007':['FAD','MLI'],'TGT_010':['ATP'],'TGT_019':['MN'],'TGT_024':['A1LVB']}
TIER_N={'TGT_002','TGT_005','TGT_008','TGT_009','TGT_013','TGT_014','TGT_015','TGT_020'}
AA=lambda rn:(lambda i:bool(i) and i.is_amino_acid())(gemmi.find_tabulated_residue(rn))
_CACHE={}
def load(struct):
    if struct in _CACHE: return _CACHE[struct].clone()
    p=f'{D}/{struct}.cif' if not struct.startswith('AF-') else f'{D}/{struct}.pdb'
    st=gemmi.read_structure(p); st.setup_entities(); _CACHE[struct]=st; return st.clone()
def uoff(struct,acc,chain):
    if struct.startswith('AF-'): return 0
    blk=gemmi.cif.read(f'{D}/{struct}.cif').sole_block()
    refacc={r[0]:r[1] for r in blk.find('_struct_ref.',['id','pdbx_db_accession'])}
    for r in blk.find('_struct_ref_seq.',['ref_id','?pdbx_strand_id','?pdbx_auth_seq_align_beg','?db_align_beg']):
        if refacc.get(r[0])==acc and r[2] not in ('?','.',None) and r[3] not in ('?','.',None):
            strands=(r[1] or '').replace(' ','').split(',')
            if chain in strands or r[1] in ('?','.',None):
                try: return int(r[2])-int(r[3])
                except: pass
    return 0
def site_res(tid, model, chains):
    if tid in ANCHOR:
        an=set(ANCHOR[tid]); ap=[a.pos for ch in model for res in ch if res.name in an for a in res]
    else:
        acc=S2[tid]['uniprot']; upr=sorted({f['start'] for f in UP[tid]['features'] if f['start']}); ap=[]
        for ch in chains:
            off=uoff(S2[tid]['structure'],acc,ch); want={u+off for u in upr}
            for c in model:
                if c.name!=ch: continue
                for res in c:
                    if AA(res.name) and res.seqid.num in want: ap+=[a.pos for a in res]
    site=set()
    for ch in model:
        for res in ch:
            if AA(res.name) and any(a.pos.dist(p)<=5.0 for a in res for p in ap): site.add((ch.name,res.seqid.num))
    return site
def parse(outdir, tag):
    info=f'{outdir}/{tag}_info.txt'; sc={}; vol={}
    if os.path.exists(info):
        cur=None
        for ln in open(info):
            m=re.match(r'Pocket (\d+)',ln)
            if m: cur=int(m.group(1))
            elif cur and 'Druggability Score' in ln: sc[cur]=float(ln.split(':')[1])
            elif cur and re.match(r'\s*Volume\s*:', ln): vol[cur]=float(ln.split(':')[1])
    pk=[]; pdir=f'{outdir}/pockets'
    if os.path.isdir(pdir):
        for fn in os.listdir(pdir):
            m=re.match(r'pocket(\d+)_atm\.pdb',fn)
            if not m: continue
            rk=int(m.group(1)); res=set()
            for ln in open(f'{pdir}/{fn}'):
                if ln.startswith(('ATOM','HETATM')):
                    try: res.add((ln[21],int(ln[22:26])))
                    except: pass
            pk.append(dict(rank=rk,ds=sc.get(rk,0.0),vol=vol.get(rk,0.0),res=res))
    return pk
def match(pk, site):
    best=None; bo=0
    for p in pk:
        o=len(p['res'] & site)
        if o>bo or (o==bo and best and p['ds']>best['ds']): best,bo=p,o
    return (best,bo) if (best and bo>=1) else (None,0)
rows=[]
for tid in sorted(S2):
    s2=S2[tid]; struct=s2['structure']; chains=s2['chains_retained'].split(';'); gene=s2['gene']
    tier='N' if tid in TIER_N else ('D' if tid in ANCHOR else 'U')
    r=dict(target_id=tid,gene=gene,tier=tier)
    if tier=='N':
        r.update(site_D1_iso='NA',site_D2_iso='NA',site_D1_sub='NA',site_D2_sub='NA',
                 site_D1_best='NA',site_D2_best='NA'); rows.append(r); continue
    # isolated
    st=load(struct); site=site_res(tid, st[0], chains)
    pk_iso=parse(f'{FP}/{tid}_{gene}_iso_out', f'{tid}_{gene}_iso')
    b,_=match(pk_iso, site)
    d1i,d2i=(b['ds'],b['vol']) if b else (0.0,0.0)
    # subcomplex
    subtag=f'{tid}_{gene}_sub'; subdir=f'{FP}/{subtag}_out'
    if os.path.isdir(subdir) and not struct.startswith('AF-'):
        stx=load(struct); m=stx[0]; keep=set(chains)
        tgt=[a.pos for ch in m if ch.name in keep for res in ch if AA(res.name) for a in res]
        nb=[]
        for ch in m:
            if ch.name in keep: continue
            n=sum(1 for res in ch if AA(res.name) and any(a.pos.dist(t)<=5.0 for a in res for t in tgt) )
            if n>=20: nb.append(ch.name)
        for ch in [c.name for c in m]:
            if ch not in keep|set(nb): m.remove_chain(ch)
        site_sub=site_res(tid, stx[0], chains)
        b2,_=match(parse(subdir, subtag), site_sub)
        d1s,d2s=(b2['ds'],b2['vol']) if b2 else (0.0,0.0)
    else:
        d1s,d2s=d1i,d2i
    d1b,d2b=(max(d1i,d1s), d2i if d1i>=d1s else d2s)
    r.update(site_D1_iso=round(d1i,3),site_D2_iso=round(d2i,1),site_D1_sub=round(d1s,3),
             site_D2_sub=round(d2s,1),site_D1_best=round(d1b,3),site_D2_best=round(d2b,1))
    rows.append(r)
    print(f"  {tid} {gene:9} {tier}  D1_best {d1b:.3f}  D2_best {d2b:8.1f}")
with open(os.path.expanduser('~/r2/out/site_scores_final.csv'),'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print('\nwrote out/site_scores_final.csv')
