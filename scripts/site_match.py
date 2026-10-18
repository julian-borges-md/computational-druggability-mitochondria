#!/usr/bin/env python3
"""Site-matched pocket analysis (reviewer comments 2, 3, 8).

Tier D  anchor = the structure's bound cofactor/inhibitor; site = residues of the
        assessed chain within 5 A of it.
Tier U  anchor = UniProt active/binding residues (author numbering via struct_ref
        offset); site = residues within 5 A of those.
Tier N  no defined small-molecule site -> structural score NA.
Isolated-chain arm.
"""
import csv, json, os, re, subprocess, gemmi
D = os.path.expanduser('~/r2/struct'); OUT = os.path.expanduser('~/r2/fp'); os.makedirs(OUT, exist_ok=True)
FP = '/tmp/fpocket_src/bin/fpocket'
S2 = {r['target_id']: r for r in csv.DictReader(open(os.path.expanduser('~/r2/data/SupplementaryTable_S2_fpocket_per_target.csv')))}
UP = json.load(open(os.path.expanduser('~/r2/out/uniprot_sites.json')))
ANCHOR = {'TGT_001':['PLX'],'TGT_003':['HEM','U10'],'TGT_004':['HEA','CU'],'TGT_006':['FMN'],
          'TGT_007':['FAD','MLI'],'TGT_010':['ATP'],'TGT_019':['MN'],'TGT_024':['A1LVB']}
TIER_N = {'TGT_002','TGT_005','TGT_008','TGT_009','TGT_013','TGT_014','TGT_015','TGT_020'}
AA = lambda rn: (lambda i: bool(i) and i.is_amino_acid())(gemmi.find_tabulated_residue(rn))
def load(struct):
    p = f'{D}/{struct}.cif' if not struct.startswith('AF-') else f'{D}/{struct}.pdb'
    st = gemmi.read_structure(p); st.setup_entities(); return st
def uniprot_offset(struct, acc, chain):
    if struct.startswith('AF-'): return 0
    doc = gemmi.cif.read(f'{D}/{struct}.cif'); blk = doc.sole_block()
    refacc = {r[0]: r[1] for r in blk.find('_struct_ref.', ['id','pdbx_db_accession'])}
    for r in blk.find('_struct_ref_seq.', ['ref_id','?pdbx_strand_id','?pdbx_auth_seq_align_beg','?db_align_beg']):
        if refacc.get(r[0])==acc and r[2] not in ('?','.',None) and r[3] not in ('?','.',None):
            strands=(r[1] or '').replace(' ','').split(',')
            if chain in strands or r[1] in ('?','.',None):
                try: return int(r[2])-int(r[3])
                except ValueError: pass
    return 0
def neighbourhood(model, chains, anchor_pos, cut=5.0):
    site=set()
    for ch in model:
        if ch.name not in chains: continue
        for res in ch:
            if AA(res.name) and any(a.pos.dist(ap)<=cut for a in res for ap in anchor_pos):
                site.add((ch.name,res.seqid.num))
    return site
def site_residues(tid):
    s2=S2[tid]; struct=s2['structure']; chains=s2['chains_retained'].split(';')
    if tid in TIER_N: return None,'N'
    st=load(struct); model=st[0]
    if tid in ANCHOR:
        anchors=set(ANCHOR[tid])
        apos=[a.pos for ch in model for res in ch if res.name in anchors for a in res]
        return neighbourhood(model,chains,apos),'D'
    acc=s2['uniprot']; upres=sorted({f['start'] for f in UP[tid]['features'] if f['start']})
    anchor_pos=[]
    for ch in chains:
        off=uniprot_offset(struct,acc,ch); want={u+off for u in upres}
        for c in model:
            if c.name!=ch: continue
            for res in c:
                if AA(res.name) and res.seqid.num in want: anchor_pos+=[a.pos for a in res]
    return neighbourhood(model,chains,anchor_pos),'U'
def build_apo(struct, chains, dest):
    st=load(struct); model=st[0]; keep=set(chains)
    for ch in [c.name for c in model]:
        if ch not in keep: model.remove_chain(ch)
    st.remove_ligands_and_waters(); st.remove_hydrogens(); st.setup_entities(); st.write_pdb(dest)
def run_fpocket(pdb, tag):
    outdir=f'{OUT}/{tag}_out'
    if not os.path.isdir(outdir): subprocess.run([FP,'-f',pdb],capture_output=True,timeout=2400)
    return outdir
def parse_pockets(outdir, tag):
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
            rank=int(m.group(1)); res=set()
            for ln in open(f'{pdir}/{fn}'):
                if ln.startswith(('ATOM','HETATM')):
                    ch=ln[21]
                    try: rn=int(ln[22:26])
                    except ValueError: continue
                    res.add((ch,rn))
            pockets.append(dict(rank=rank,drugscore=scores.get(rank,0.0),residues=res))
    return pockets
rows=[]
for tid in sorted(S2):
    s2=S2[tid]; struct=s2['structure']; chains=s2['chains_retained'].split(';'); gene=s2['gene']
    tag=f'{tid}_{gene}_iso'; site,tier=site_residues(tid)
    apo=f'{OUT}/{tag}.pdb'
    if not os.path.exists(apo): build_apo(struct,chains,apo)
    pk=parse_pockets(run_fpocket(apo,tag),tag)
    global_d1=max((p['drugscore'] for p in pk),default=0.0)
    if tier=='N':
        site_d1,ov,nsite='NA','',0
    else:
        nsite=len(site); best=None; bestov=0
        for p in pk:
            o=len(p['residues'] & site)
            if o>bestov or (o==bestov and best and p['drugscore']>best['drugscore']): best,bestov=p,o
        site_d1=round(best['drugscore'],3) if (best and bestov>=1) else 0.0; ov=bestov
    rows.append(dict(target_id=tid,gene=gene,tier=tier,n_site_res=nsite,n_pockets=len(pk),
                     site_overlap=ov,site_matched_D1=site_d1,global_best_D1=round(global_d1,3)))
    print(f"  {tid} {gene:9} tier {tier}  site_res {nsite:3}  pockets {len(pk):3}  "
          f"overlap {str(ov):>3}  site_D1 {str(site_d1):>6}  global_D1 {global_d1:.3f}")
with open(os.path.expanduser('~/r2/out/site_matched_isolated.csv'),'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print('\nwrote out/site_matched_isolated.csv')
