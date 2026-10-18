#!/usr/bin/env python3
"""Marked-changes manuscript as TRUE Word tracked changes (w:ins / w:del), so
that in Word 'Accept all' reconstructs the R2 manuscript and 'Reject all'
reconstructs the R1 manuscript (preflight group G). Word-level diff of the R1
manuscript against the R2 manuscript. Produces build/R2_manuscript_marked.docx.
"""
import re, os, difflib
from docx import Document
from docx.shared import Pt
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
R=os.path.expanduser('~/r2')
AUTHOR='J Y V Borges'; DATE='2026-09-30T00:00:00Z'

def strip_md(md):
    out=[]
    for ln in md.splitlines():
        s=re.sub(r'^#{1,6}\s*','',ln)
        s=s.replace('**','').replace('__','')
        s=re.sub(r'^\s*[-*]\s+','',s)
        out.append(s)
    return '\n'.join(out)

r1=open(f'{R}/docs/R1_manuscript.txt',encoding='utf-8').read()
r2=strip_md(open(f'{R}/build/R2_manuscript.md',encoding='utf-8').read())
def tok(t): return re.findall(r'\n|\S+|[ ]+', t)
a=tok(r1); b=tok(r2)
sm=difflib.SequenceMatcher(a=a,b=b,autojunk=False)

doc=Document()
doc.styles['Normal'].font.name='Calibri'; doc.styles['Normal'].font.size=Pt(10.5)
h=doc.add_paragraph().add_run('CSBR-D-26-00031 — Revised Manuscript with Marked Changes (R1 to R2)')
h.bold=True; h.font.size=Pt(12)
lg=doc.add_paragraph()
lg.add_run('Word tracked changes. Accept All reconstructs the revised (R2) manuscript; '
           'Reject All reconstructs the previous (R1) manuscript.').italic=True

_id=[0]
def rid():
    _id[0]+=1; return str(_id[0])
def run_el(text, deleted=False):
    r=OxmlElement('w:r')
    t=OxmlElement('w:delText' if deleted else 'w:t')
    t.set(qn('xml:space'),'preserve'); t.text=text
    r.append(t); return r
def ins_el(text):
    ins=OxmlElement('w:ins'); ins.set(qn('w:id'),rid()); ins.set(qn('w:author'),AUTHOR); ins.set(qn('w:date'),DATE)
    ins.append(run_el(text,deleted=False)); return ins
def del_el(text):
    d=OxmlElement('w:del'); d.set(qn('w:id'),rid()); d.set(qn('w:author'),AUTHOR); d.set(qn('w:date'),DATE)
    d.append(run_el(text,deleted=True)); return d

cur=doc.add_paragraph(); pel=cur._p
def emit(tokens, kind):
    global cur, pel
    buf=''
    def flush():
        nonlocal buf
        if not buf: return
        if kind=='eq': pel.append(run_el(buf))
        elif kind=='ins': pel.append(ins_el(buf))
        elif kind=='del': pel.append(del_el(buf))
        buf=''
    for t in tokens:
        if t=='\n':
            flush(); cur=doc.add_paragraph(); pel=cur._p
        else:
            buf+=t
    flush()

for tag,i1,i2,j1,j2 in sm.get_opcodes():
    if tag=='equal': emit(a[i1:i2],'eq')
    elif tag=='delete': emit(a[i1:i2],'del')
    elif tag=='insert': emit(b[j1:j2],'ins')
    elif tag=='replace': emit(a[i1:i2],'del'); emit(b[j1:j2],'ins')

doc.save(f'{R}/build/R2_manuscript_marked.docx')
print("wrote build/R2_manuscript_marked.docx (true w:ins/w:del tracked changes)")
print(f"similarity {sm.ratio():.2f}")
