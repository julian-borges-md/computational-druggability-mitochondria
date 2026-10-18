#!/usr/bin/env python3
"""Assemble the CSBR-D-26-00031 R2 upload package to the standing governance
layout (submission-package-layout). Generates the directory tree, strips
document provenance from every uploaded file, and writes the generated manifest
and upload order. Idempotent: rebuilds the package from scratch each run."""
import os, re, csv, shutil, zipfile, hashlib, subprocess, datetime
R=os.path.expanduser('~/r2'); B=f'{R}/build'; OUT=f'{R}/out'; TB=f'{OUT}/tables'
PKG=f'{R}/CSBR-D-26-00031_R2_PACKAGE'
if os.path.exists(PKG): shutil.rmtree(PKG)

DIRS=['','00_UPLOAD/01_cover_letter','00_UPLOAD/02_response','00_UPLOAD/03_manuscript',
      '00_UPLOAD/04_manuscript_marked','00_UPLOAD/05_tables','00_UPLOAD/06_supplementary',
      '01_ZENODO','02_PORTAL_FIELDS','03_SCRIPTS','04_AUDIT','05_ARCHIVE']
for d in DIRS: os.makedirs(f'{PKG}/{d}', exist_ok=True)

# ---------- provenance scrub for uploaded OOXML (group H) ----------
def scrub_docx(path):
    tmp=path+'.tmp'
    with zipfile.ZipFile(path) as zin, zipfile.ZipFile(tmp,'w',zipfile.ZIP_DEFLATED) as zout:
        for it in zin.namelist():
            data=zin.read(it)
            if it=='docProps/core.xml':
                data=('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                      '<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" '
                      'xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" '
                      'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'
                      '<dc:creator></dc:creator><cp:lastModifiedBy></cp:lastModifiedBy>'
                      '<dc:title></dc:title></cp:coreProperties>').encode()
            elif it=='docProps/app.xml':
                data=('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                      '<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties">'
                      '<Application></Application></Properties>').encode()
            zout.writestr(it,data)
    os.replace(tmp,path)

def scrub_xlsx(path):
    tmp=path+'.tmp'
    with zipfile.ZipFile(path) as zin, zipfile.ZipFile(tmp,'w',zipfile.ZIP_DEFLATED) as zout:
        for it in zin.namelist():
            data=zin.read(it)
            if it=='docProps/core.xml':
                data=('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                      '<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" '
                      'xmlns:dc="http://purl.org/dc/elements/1.1/"><dc:creator></dc:creator>'
                      '<cp:lastModifiedBy></cp:lastModifiedBy></cp:coreProperties>').encode()
            elif it=='docProps/app.xml':
                data=('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                      '<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties">'
                      '<Application></Application></Properties>').encode()
            zout.writestr(it,data)
    os.replace(tmp,path)

# ---------- cover letter -> docx ----------
subprocess.run(['pandoc',f'{B}/R2_cover_letter.md','-o',f'{B}/R2_cover_letter.docx'],check=True)

# ---------- place uploaded deliverables ----------
UPLOAD={
 '00_UPLOAD/01_cover_letter/CSBR-D-26-00031_R2_cover_letter.docx': f'{B}/R2_cover_letter.docx',
 '00_UPLOAD/02_response/CSBR-D-26-00031_R2_response_to_reviewers.docx': f'{B}/R2_response.docx',
 '00_UPLOAD/03_manuscript/CSBR-D-26-00031_R2_manuscript_clean.docx': f'{B}/R2_manuscript_clean.docx',
 '00_UPLOAD/04_manuscript_marked/CSBR-D-26-00031_R2_manuscript_marked.docx': f'{B}/R2_manuscript_marked.docx',
 '00_UPLOAD/05_tables/CSBR-D-26-00031_R2_tables.docx': f'{B}/R2_tables.docx',
 '00_UPLOAD/06_supplementary/CSBR-D-26-00031_R2_supplementary_tables_S1_S8.xlsx': f'{TB}/Supplementary_Tables_R2.xlsx',
}
for dst,src in UPLOAD.items():
    shutil.copy2(src, f'{PKG}/{dst}')
    if dst.endswith('.docx'): scrub_docx(f'{PKG}/{dst}')
    elif dst.endswith('.xlsx'): scrub_xlsx(f'{PKG}/{dst}')

# ---------- 01_ZENODO ----------
shutil.copy2(f'{OUT}/druggability_registry_v1.3.yaml', f'{PKG}/01_ZENODO/druggability_registry_v1.3.yaml')
shutil.copy2(f'{TB}/Supplementary_Tables_R2.xlsx', f'{PKG}/01_ZENODO/Supplementary_Tables_S1_S8.xlsx')
for f in os.listdir(TB):
    if f.endswith('.csv'): shutil.copy2(f'{TB}/{f}', f'{PKG}/01_ZENODO/{f}')
open(f'{PKG}/01_ZENODO/HOW_TO_DEPOSIT.md','w').write(
 "# Zenodo deposit (new version)\n\n"
 "Create a NEW VERSION of the existing record https://doi.org/10.5281/zenodo.19389356.\n\n"
 "Contents of this new version:\n"
 "- druggability_registry_v1.3.yaml (supersedes v1.0-v1.2)\n"
 "- Supplementary_Tables_S1_S8.xlsx and the per-table CSVs\n"
 "- the full 03_SCRIPTS reproduction and verification set\n\n"
 "Keep the concept DOI in the manuscript; the version DOI updates automatically.\n")

# ---------- 02_PORTAL_FIELDS ----------
man=open(f'{B}/R2_manuscript.md').read()
abs=re.search(r'## Abstract\n(.+?)\nKeywords:', man, re.S).group(1)
abs_plain=re.sub(r'\*\*','',abs).strip()
open(f'{PKG}/02_PORTAL_FIELDS/abstract.txt','w').write(abs_plain+'\n')
kw=re.search(r'Keywords:(.+)', man).group(1).strip()
open(f'{PKG}/02_PORTAL_FIELDS/keywords.txt','w').write(kw+'\n')
open(f'{PKG}/02_PORTAL_FIELDS/correspondence.txt','w').write(
 "Corresponding author: Julian Yin Vieira Borges, MD\n"
 "ORCID: 0009-0001-9929-3135\n"
 "Email: jyborges@bu.edu\n"
 "Affiliation: Department of Computer Science, Boston University, Boston, MA, USA\n"
 "Manuscript: CSBR-D-26-00031 (second revision, R2)\n")

# ---------- 03_SCRIPTS ----------
for f in os.listdir(f'{R}/scripts'):
    if f.endswith('.py'): shutil.copy2(f'{R}/scripts/{f}', f'{PKG}/03_SCRIPTS/{f}')
# key data the scripts read
os.makedirs(f'{PKG}/03_SCRIPTS/out', exist_ok=True)
for f in ['matrix_v13.csv','site_scores_final.csv','apriori_sites.csv','chembl_filtered.csv',
          'final_stats.json','derived_stats.json','weight_sensitivity.json','af_confidence.csv','provenance.csv']:
    shutil.copy2(f'{OUT}/{f}', f'{PKG}/03_SCRIPTS/out/{f}')

# ---------- 04_AUDIT ----------
shutil.copy2(f'{B}/CHANGE_SUMMARY.md', f'{PKG}/04_AUDIT/CHANGE_SUMMARY.md')

# ---------- 05_ARCHIVE (verbatim prior versions, read-only) ----------
for f in os.listdir(f'{R}/docs'):
    if f.startswith('R1_') : shutil.copy2(f'{R}/docs/{f}', f'{PKG}/05_ARCHIVE/{f}')
open(f'{PKG}/05_ARCHIVE/WHY_THIS_IS_HERE.md','w').write(
 "# Archive (read only)\n\nVerbatim R1 submission documents, retained for the marked-changes diff and "
 "for audit. Never edited or typographically normalised.\n")

print("package tree assembled at", PKG)
