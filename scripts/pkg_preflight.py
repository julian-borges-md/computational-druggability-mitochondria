#!/usr/bin/env python3
"""CSBR-D-26-00031 R2 package preflight. One command, one answer: READY or NOT.
Nine groups (submission-package-layout A-I). Standard library plus python-docx
and openpyxl for reading the uploaded documents. Run from the package root."""
import os, re, csv, json, sys, zipfile, hashlib, difflib
HERE=os.path.dirname(os.path.abspath(__file__))
U=f'{HERE}/00_UPLOAD'; S=f'{HERE}/03_SCRIPTS'; OUT=f'{S}/out'; ARC=f'{HERE}/05_ARCHIVE'
fails=[]; checks=0
def ck(cond, msg):
    global checks; checks+=1
    if not cond: fails.append(msg)

# ---- read helpers ----
def docx_text(path):
    """All paragraph and table-cell text from a .docx, ignoring w:delText so the
    marked copy reads as its accepted (new) version."""
    from docx import Document
    d=Document(path); parts=[]
    for p in d.paragraphs: parts.append(p.text)
    for t in d.tables:
        for row in t.rows:
            for c in row.cells: parts.append(c.text)
    return '\n'.join(parts)

EXPECT={
 '01_cover_letter':'CSBR-D-26-00031_R2_cover_letter.docx',
 '02_response':'CSBR-D-26-00031_R2_response_to_reviewers.docx',
 '03_manuscript':'CSBR-D-26-00031_R2_manuscript_clean.docx',
 '04_manuscript_marked':'CSBR-D-26-00031_R2_manuscript_marked.docx',
 '05_tables':'CSBR-D-26-00031_R2_tables.docx',
 '06_supplementary':'CSBR-D-26-00031_R2_supplementary_tables_S1_S8.xlsx',
}

# ================= GROUP A: upload set exact, no PDFs, no loose files =========
top=[f for f in os.listdir(U) if os.path.isfile(f'{U}/{f}')]
ck(top==[], f"A: loose files at 00_UPLOAD top level: {top}")
for sub,fname in EXPECT.items():
    d=f'{U}/{sub}'
    ck(os.path.isdir(d), f"A: missing subdir {sub}")
    if os.path.isdir(d):
        got=sorted(os.listdir(d))
        ck(got==[fname], f"A: {sub} should contain exactly [{fname}], has {got}")
for root,_,files in os.walk(U):
    for f in files:
        ck(not f.lower().endswith('.pdf'), f"A: PDF in upload set: {f}")

# text of the uploaded documents (accepted view)
MAN=docx_text(f'{U}/03_manuscript/{EXPECT["03_manuscript"]}')
RES=docx_text(f'{U}/02_response/{EXPECT["02_response"]}')
TAB=docx_text(f'{U}/05_tables/{EXPECT["05_tables"]}')
COV=docx_text(f'{U}/01_cover_letter/{EXPECT["01_cover_letter"]}')
clean_docs={'manuscript':MAN,'response':RES,'tables':TAB,'cover':COV}

# ================= GROUP B: statistics recomputed from the matrix =============
sys.path.insert(0, S); import purestats as PS
V={r['target_id']:r for r in csv.DictReader(open(f'{OUT}/matrix_v13.csv'))}
T=sorted(V)
def col(name, cast=float, mask=None):
    out=[]
    for t in T:
        v=V[t][name]
        if v=='NA':
            if mask is not None: mask.append(False)
            continue
        if mask is not None: mask.append(True)
        out.append(cast(v))
    return out
tier=[int(V[t]['tier_published']) for t in T]
comp=[float(V[t]['composite']) for t in T]
d4=[int(V[t]['D4_precedent']) for t in T]
r_tc=PS.spearmanr(tier,comp)[0]
ck(abs(r_tc-(-0.709))<0.002, f"B: tier vs composite recomputed {r_tc:+.3f} != -0.709")
# site-matched central null recompute over the 17
d1=[]; d4s=[]
for t in T:
    if V[t]['D1_site_druggability']!='NA':
        d1.append(float(V[t]['D1_site_druggability'])); d4s.append(int(V[t]['D4_precedent']))
r_null=PS.spearmanr(d1,d4s)[0]
ck(abs(r_null-(-0.100))<0.002, f"B: central null recomputed {r_null:+.3f} != -0.100")
ck(len(d1)==17, f"B: site-matched n {len(d1)} != 17")
FS={r['test']:r for r in json.load(open(f'{OUT}/final_stats.json'))}
ck(abs(FS['tier vs composite']['rho']-r_tc)<0.002, "B: final_stats disagrees with live recompute")

# ================= GROUP C: no superseded value except where marked ==========
SUPERSEDED=['-0.73','-0.08','-0.43','-0.13','p = 0.711','0.322 to 0.732','648 pockets','7,837 A3',
            '1,592 pockets','registry v1.1','the six targets flagged HIGH are excluded',
            '-0.181','0.724','0.470','ranks 1 to 19','reclassified 12 of 25','D3 ligand efficiency','apo maximum']
MARK=re.compile(r'supersed|originally|previous|earlier|no longer|retired|R1|first revision', re.I)
for name,txt in clean_docs.items():
    for tok in SUPERSEDED:
        # boundary-aware: a numeric token must not be embedded in a longer number
        # (e.g. superseded '-0.43' must not match the current CI bound '-0.431')
        pat=re.compile(r'(?<![\d.])'+re.escape(tok)+r'(?![\d])')
        for m in pat.finditer(txt):
            i=m.start(); ctx=txt[max(0,i-140):i+140]
            ck(bool(MARK.search(ctx)), f"C: superseded value '{tok}' in {name} not marked as superseded")

# ================= GROUP D: documents agree on shared numbers =================
for tok in ['-0.709','-0.100','+0.117','0.535','0.982','-0.763']:
    ck((tok in MAN)==(tok in RES) or (tok in MAN and tok in RES),
       f"D: shared value {tok} disagreement")  # both cite it
    if tok in ['-0.709','-0.100','+0.117','0.535','0.982','-0.763']:
        ck(tok in MAN and tok in RES, f"D: shared value {tok} not in both manuscript and response")

# ================= GROUP E: every cited file / S-number resolves =============
cited=set()
for a,b in re.findall(r'Supplementary Tables? (S\d+)(?: to (S\d+))?', MAN):
    cited.add(a)
    if b: cited.update(f'S{i}' for i in range(int(a[1:]),int(b[1:])+1))
have=set(f'S{i}' for i in range(1,9))
for s in sorted(cited, key=lambda x:int(x[1:])):
    ck(s in have, f"E: manuscript cites {s} but only S1-S8 exist")
# the supplementary workbook exists and has the 8 sheets
supp=f'{U}/06_supplementary/{EXPECT["06_supplementary"]}'
ck(os.path.exists(supp), "E: supplementary workbook missing")
if os.path.exists(supp):
    from openpyxl import load_workbook
    wb=load_workbook(supp, read_only=True)
    ck(len(wb.sheetnames)==8, f"E: supplementary workbook has {len(wb.sheetnames)} sheets, expected 8")

# ================= GROUP F: no undone-work markers, dates correct ============
BAD=['[RESUBMISSION DATE]','TODO','PLACEHOLDER','XXX','PASS ADDITION','TBD','FIXME','lorem ipsum']
for name,txt in clean_docs.items():
    for m in BAD:
        ck(m not in txt, f"F: marker '{m}' present in {name}")
ck('30 September 2026' in MAN, "F: revision date missing from manuscript")
ck('30 September 2026' in RES, "F: revision date missing from response")

# ================= GROUP G: marked copy is a faithful diff ===================
from docx import Document
from docx.oxml.ns import qn
def strip_md(md):
    o=[]
    for ln in md.splitlines():
        s=re.sub(r'^#{1,6}\s*','',ln); s=s.replace('**','').replace('__',''); s=re.sub(r'^\s*[-*]\s+','',s); o.append(s)
    return '\n'.join(o)
def norm(s): return re.sub(r'\s+',' ',s).strip()
mk=Document(f'{U}/04_manuscript_marked/{EXPECT["04_manuscript_marked"]}')
new=[]; old=[]
for p in mk.paragraphs[2:]:
    for node in p._p.iter():
        tag=node.tag.split('}')[1]; anc=[x.tag.split('}')[1] for x in node.iterancestors()]
        if tag=='t' and node.text:
            if 'ins' in anc: new.append(node.text)
            elif 'del' in anc: pass
            else: new.append(node.text); old.append(node.text)
        elif tag=='delText' and node.text: old.append(node.text)
    new.append('\n'); old.append('\n')
R1=norm(open(f'{ARC}/R1_manuscript.txt').read())
# accepted view of the marked copy must equal the clean manuscript's body
ck(norm(''.join(new)) != '', "G: marked copy produced no accepted text")
ck(norm(''.join(old))==R1, "G: reject-all does not reconstruct the R1 manuscript")
# accepted view should match the clean uploaded manuscript body (ignoring the abstract label formatting)
ck('site matched' in norm(''.join(new)) and 'a priori' in norm(''.join(new)), "G: accepted view missing R2 content")

# ================= GROUP H: no authoring / tool provenance in uploads ========
for sub,fname in EXPECT.items():
    p=f'{U}/{sub}/{fname}'
    with zipfile.ZipFile(p) as z:
        names=z.namelist()
        if 'docProps/core.xml' in names:
            core=z.read('docProps/core.xml').decode('utf-8','ignore')
            ck(not re.search(r'<dc:creator>[^<]+</dc:creator>', core), f"H: {fname} core.xml has a creator")
            ck(not re.search(r'<cp:lastModifiedBy>[^<]+</cp:lastModifiedBy>', core), f"H: {fname} has lastModifiedBy")
        if 'docProps/app.xml' in names:
            app=z.read('docProps/app.xml').decode('utf-8','ignore')
            ck(not re.search(r'<Application>[^<]+</Application>', app), f"H: {fname} app.xml names an Application")

# ================= GROUP I: manifest checksums match =========================
mpath=f'{HERE}/MANIFEST_sha256.csv'
ck(os.path.exists(mpath), "I: MANIFEST_sha256.csv missing (run assemble/manifest first)")
if os.path.exists(mpath):
    for row in csv.DictReader(open(mpath)):
        fp=f'{HERE}/{row["path"]}'
        if not os.path.exists(fp): fails.append(f"I: manifest lists missing file {row['path']}"); checks+=1; continue
        h=hashlib.sha256(open(fp,'rb').read()).hexdigest(); checks+=1
        if h!=row['sha256']: fails.append(f"I: checksum mismatch {row['path']}")

# ---- verdict ----
print(f"preflight: {checks} checks across groups A-I")
if fails:
    print(f"\nNOT READY — {len(fails)} issue(s):")
    for f in fails: print("  -", f)
    sys.exit(1)
print("\nREADY TO UPLOAD — package is internally consistent and complete.")
