#!/usr/bin/env python3
"""Finalize the R2 package: place the package preflight, write README and STATE,
and GENERATE UPLOAD_ORDER.md and MANIFEST_sha256.csv from the directory (never by
hand). Run after assemble_package.py."""
import os, csv, hashlib, shutil, datetime
R=os.path.expanduser('~/r2'); PKG=f'{R}/CSBR-D-26-00031_R2_PACKAGE'
today='30 September 2026'

# package preflight in place
shutil.copy2(f'{R}/scripts/pkg_preflight.py', f'{PKG}/preflight.py')

# README
open(f'{PKG}/README.md','w').write(f"""# CSBR-D-26-00031 — R2 upload package

Second revision of "Computational Druggability Assessment of Mitochondrial Targets".
Assembled {today}. Author: Julian Yin Vieira Borges, MD (ORCID 0009-0001-9929-3135).

## Three rules
1. **Everything in `00_UPLOAD/` is uploaded; nothing else is.** No PDFs, no README, no loose files there.
2. **`UPLOAD_ORDER.md` and `MANIFEST_sha256.csv` are generated, never hand edited.** Re-run `finalize_package.py` to refresh them.
3. **`05_ARCHIVE/` is read only.** It holds the verbatim R1 documents for the marked-changes diff and audit.

## What is where
| Path | Contents |
|---|---|
| `00_UPLOAD/01_cover_letter` | Cover letter to the editor |
| `00_UPLOAD/02_response` | Point-by-point response to the 11 reviewer comments |
| `00_UPLOAD/03_manuscript` | Clean revised manuscript (R2) |
| `00_UPLOAD/04_manuscript_marked` | Manuscript with Word tracked changes (R1 to R2) |
| `00_UPLOAD/05_tables` | Main tables 1, 1b, 2, 3, 4 |
| `00_UPLOAD/06_supplementary` | Supplementary tables S1-S8 (Excel workbook) |
| `01_ZENODO` | Registry v1.3, supplementary data, reproduction scripts for the deposit |
| `02_PORTAL_FIELDS` | Abstract, keywords, correspondence typed into Editorial Manager |
| `03_SCRIPTS` | Reproduction and verification scripts, plus the input data they read |
| `04_AUDIT` | Change summary by comment |
| `05_ARCHIVE` | Read-only verbatim R1 documents |

## Verify before upload
Run `python3 preflight.py` from this folder. It must print `READY TO UPLOAD`.
Expected gate: 9 groups (A upload set, B recomputed stats, C no superseded value, D documents agree,
E cited files resolve, F no undone-work markers, G marked copy faithful, H no provenance metadata,
I manifest checksums).

## Upload order
Follow `UPLOAD_ORDER.md`. The response letter and the marked manuscript both carry the changes;
the clean manuscript is the version of record.
""")

# STATE
open(f'{PKG}/STATE.md','w').write(f"""# STATE — CSBR-D-26-00031 as of {today}

## What the journal holds
- R1 submission (11 September 2026) and the reviewer's decision of revise and resubmit, due 29 October 2026.

## What this package holds (R2, ready to upload)
- Cover letter, point-by-point response to 11 comments, clean manuscript, marked manuscript, tables, supplementary S1-S8.
- Registry v1.3 and reproduction scripts staged for a new Zenodo version of record 10.5281/zenodo.19389356.

## What changed since R1
- Structural dimension re-scored at an a priori therapeutic site in the isolated and biological-assembly
  contexts. Central null holds (site-matched structure vs precedent -0.100; tier vs site druggability -0.037).
- D3 renamed bioactivity evidence (pChEMBL count). Composite fully specified with weight robustness.
- Permutation and bootstrap inference; MDES 0.535. Regulatory claims cite primary FDA records.
- Single-source generation with an automated consistency check (preflight).

## Open items for the author before upload
- Confirm the Zenodo new version is published and the version DOI resolves before final submission.
- Confirm Editorial Manager still holds the withheld-SMILES restriction note from R1.
- Nothing else outstanding; preflight passes.
""")

# UPLOAD_ORDER.md — generated from 00_UPLOAD
U=f'{PKG}/00_UPLOAD'
order=['01_cover_letter','02_response','03_manuscript','04_manuscript_marked','05_tables','06_supplementary']
labels={'01_cover_letter':'Cover Letter','02_response':'Response to Reviewers',
        '03_manuscript':'Revised Manuscript (clean)','04_manuscript_marked':'Revised Manuscript with Marked Changes',
        '05_tables':'Tables','06_supplementary':'Supplementary Material'}
with open(f'{PKG}/UPLOAD_ORDER.md','w') as f:
    f.write("# UPLOAD ORDER (generated from 00_UPLOAD — do not hand edit)\n\n")
    f.write("| # | Editorial Manager item | File |\n|---|---|---|\n")
    for i,sub in enumerate(order,1):
        files=os.listdir(f'{U}/{sub}')
        fn=files[0] if files else '(missing)'
        f.write(f"| {i} | {labels[sub]} | `{sub}/{fn}` |\n")

# MANIFEST_sha256.csv — generated from the whole package (excludes itself)
rows=[]
for root,_,files in os.walk(PKG):
    for fn in sorted(files):
        fp=os.path.join(root,fn); rel=os.path.relpath(fp,PKG)
        if rel=='MANIFEST_sha256.csv': continue
        h=hashlib.sha256(open(fp,'rb').read()).hexdigest()
        rows.append({'path':rel,'bytes':os.path.getsize(fp),'sha256':h})
rows.sort(key=lambda r:r['path'])
with open(f'{PKG}/MANIFEST_sha256.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['path','bytes','sha256']); w.writeheader(); w.writerows(rows)
print(f"finalized: {len(rows)} files manifested, README/STATE/UPLOAD_ORDER/preflight in place")
