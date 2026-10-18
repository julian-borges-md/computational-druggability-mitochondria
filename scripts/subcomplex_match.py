#!/usr/bin/env python3
"""Biological-subcomplex site-matched arm (reviewer comment 3).

For each target with a defined site, build the target-containing neighbourhood
= assessed chain(s) + every chain with >=20 atoms within 5 A of them (never the
whole megacomplex), run fpocket on it, and match the a priori site (anchor +/-
5 A, now allowed to include interface residues) to the best-overlapping pocket.
Monomeric targets reduce to the isolated chain.
"""
import csv, json, os, re, subprocess, gemmi
D=os.path.expanduser('~/r2/struct'); OUT=os.path.expanduser('~/r2/fp'); FP='/tmp/fpocket_src/bin/fpocket'
S2={r['target_id']:r for r in csv.DictReader(open(os.path.expanduser('~/r2/data/SupplementaryTable_S2_fpocket_per_target.csv')))}
UP=json.load(open(os.path.expanduser('~/r2/out/uniprot_sites.json')))
ANCHOR={'TGT_001':['PLX'],'TGT_003':['HEM','U10'],'TGT_004':['HEA','CU'],'TGT_006':['FMN'],
        'TGT_007':['FAD','MLI'],'TGT_010':['ATP'],'TGT_019':['MN'],'TGT_024':['A1LVB']}
TIER_N={'TGT_002','TGT_005','TGT_008','TGT_009','TGT_013','TGT_014','TGT_015','TGT_020'}
AA=lambda rn:(lambda i:bool(i) and i.is_amino_acid())(gemmi.find_tabulated_residue(rn))
def load(struct):
    p=f'{D}/{struct}.cif' if not struct.startswith('AF-') else f'{D}/{struct}.pdb'
    st=gemmi.read_structure(p); st.setup_entities(); return st
def uniprot_offset(struct,acc,chain):
    if struct.startswith('AF-'): return 0
    doc=gemmi.cif.read(f'{D}/{struct}.cif'); blk=doc.sole_block()
    refacc={r[0]:r[1] for r in blk.find('_struct_ref.',['id','pdbx_db_accession'])}
    for r in blk.find('_struct_ref_seq.',['ref_id','?pdbx_strand_id','?pdbx_auth_seq_align_beg','?db_align_beg']):
        if refacc.get(r[0])==acc and r[2] not in ('?','.',None) and r[3] not in ('?','.',None):
            strands=(r[1] or '').replace(' ','').split(',')
            if chain in strands or r[1] in ('?','.',None):
                try: return int(r[2])-int(r[3])
                except ValueError: pass
    return 0
def anchor_positions(tid, model, chains):
    if tid in ANCHOR:
        anchors=set(ANCHOR[tid])
        return [a.pos for ch in model for res in ch if res.name in anchors for a in res]
    s2=S2[tid]; acc=s2['uniprot']; upres=sorted({f['start'] for f in UP[tid]['features'] if f['start']})
    pos=[]
    for ch in chains:
        off=uniprot_offset(s2['structure'],acc,ch); want={u+off for u in upres}
        for c in model:
            if c.name!=ch: continue
            for res in c:
                if AA(res.name) and res.seqid.num in want: pos+=[a.pos for a in res]
    return pos
def neighbour_chains(model, chains, cut=5.0, minatoms=20):
    keep=set(chains)
    tgt=[a.pos for ch in model if ch.name in keep for res in ch if AA(res.name) for a in res]
    counts={}
    for ch in model:
        if ch.name in keep: continue
        n=0
        for res in ch:
            if not AA(res.name): continue
            for a in res:
                if any(a.pos.dist(t)<=cut for t in tgt): n+=1; break
        if n>=minatoms: counts[ch.name]=n
    return sorted(counts, key=lambda c:-counts[c])
def build_subcomplex(tid):
    s2=S2[tid]; struct=s2['structure']; chains=s2['chains_retained'].split(';')
    st=load(struct); model=st[0]
    if struct.startswith('AF-'):
        nb=[]
    else:
        nb=neighbour_chains(model,chains)
    keep=set(chains)|set(nb)
    for ch in [c.name for c in model]:
        if ch not in keep: model.remove_chain(ch)
    # site over the retained coords (anchor +/- 5 A, any retained chain)
    apos=anchor_positions(tid, st[0], chains)
    site=set()
    for ch in st[0]:
        for res in ch:
            if AA(res.name) and any(a.pos.dist(ap)<=5.0 for a in res for ap in apos):
                site.add((ch.name,res.seqid.num))
    st.remove_ligands_and_waters(); st.remove_hydrogens(); st.setup_entities()
    dest=f'{OUT}/{tid}_{s2["gene"]}_sub.pdb'; st.write_pdb(dest)
    return dest, site, chains, nb
def run_fpocket(pdb,tag):
    outdir=f'{OUT}/{tag}_out'
    if not os.path.isdir(outdir): subprocess.run([FP,'-f',pdb],capture_output=True,timeout=240)
    return outdir
def parse_pockets(outdir,tag):
    info=f'{outdir}/{tag}_info.txt'; scores={}
    if os.path.exists(info):
        cur=None
        for ln in open(info):
            m=re.match(r'Pocket (\d+)',ln)
            if m: cur=int(m.group(1))
            elif cur and 'Druggability Score' in ln: scores[cur]=float(ln.split(':')[1])
    pockets=[]; pdir=f'{outdir}/pockets'
    if os.path.isdir(pdir):
        for fn in os.listdir(pdir):
            m=re.match(r'pocket(\d+)_atm\.pdb',fn)
            if not m: continue
            rank=int(m.group(1)); res=set(); chset=set()
            for ln in open(f'{pdir}/{fn}'):
                if ln.startswith(('ATOM','HETATM')):
                    ch=ln[21]
                    try: rn=int(ln[22:26])
                    except ValueError: continue
                    res.add((ch,rn)); chset.add(ch)
            pockets.append(dict(rank=rank,drugscore=scores.get(rank,0.0),residues=res,nchains=len(chset)))
    return pockets
rows=[]
for tid in sorted(S2):
    s2=S2[tid]; gene=s2['gene']
    if tid in TIER_N:
        rows.append(dict(target_id=tid,gene=gene,tier='N',n_neighbours=0,site_matched_D1_sub='NA',
                         interface_pocket='',overlap='')); 
        print(f"  {tid} {gene:9} tier N  NA"); continue
    tag=f'{tid}_{gene}_sub'
    dest,site,chains,nb=build_subcomplex(tid)
    pk=parse_pockets(run_fpocket(dest,tag),tag)
    best=None; bestov=0
    for p in pk:
        o=len(p['residues'] & site)
        if o>bestov or (o==bestov and best and p['drugscore']>best['drugscore']): best,bestov=p,o
    smd=round(best['drugscore'],3) if (best and bestov>=1) else 0.0
    iface = (best['nchains']>1) if best else False
    tier='D' if tid in ANCHOR else 'U'
    rows.append(dict(target_id=tid,gene=gene,tier=tier,n_neighbours=len(nb),
                     site_matched_D1_sub=smd,interface_pocket='yes' if iface else 'no',overlap=bestov))
    print(f"  {tid} {gene:9} tier {tier}  +{len(nb):2d} neighbour chains  site_D1_sub {smd:<6}  "
          f"{'INTERFACE' if iface else 'intra'}  overlap {bestov}")
with open(os.path.expanduser('~/r2/out/site_matched_subcomplex.csv'),'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print('\nwrote out/site_matched_subcomplex.csv')
