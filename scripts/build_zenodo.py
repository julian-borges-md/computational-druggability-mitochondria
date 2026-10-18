#!/usr/bin/env python3
"""Assemble the Zenodo v1.3 new-version deposit bundle for CSBR-D-26-00031,
mirroring the R1 deposit layout (README, HOW_TO_DEPOSIT, version notes, data/,
scripts/). Everything is drawn from the verified R2 outputs."""
import os, shutil, json
R=os.path.expanduser('~/r2'); OUT=f'{R}/out'; TB=f'{OUT}/tables'
Z=f'{R}/ZENODO_v1.3'
if os.path.exists(Z): shutil.rmtree(Z)
os.makedirs(f'{Z}/data', exist_ok=True); os.makedirs(f'{Z}/data/carried_forward_from_R1', exist_ok=True)
os.makedirs(f'{Z}/scripts', exist_ok=True)

# ---- data: registry, supplementary tables, underlying outputs ----
shutil.copy2(f'{OUT}/druggability_registry_v1.3.yaml', f'{Z}/data/druggability_registry_v1.3.yaml')
shutil.copy2(f'{TB}/Supplementary_Tables_R2.xlsx', f'{Z}/data/Supplementary_Tables_S1_S8.xlsx')
for f in os.listdir(TB):
    if f.startswith('S') and f.endswith('.csv'): shutil.copy2(f'{TB}/{f}', f'{Z}/data/{f}')
for f in ['matrix_v13.csv','site_scores_final.csv','apriori_sites.csv','chembl_filtered.csv','provenance.csv',
          'final_stats.json','derived_stats.json','weight_sensitivity.json','af_confidence.csv']:
    shutil.copy2(f'{OUT}/{f}', f'{Z}/data/{f}')

# ---- scripts: reproduction + verification ----
for f in os.listdir(f'{R}/scripts'):
    if f.endswith('.py'): shutil.copy2(f'{R}/scripts/{f}', f'{Z}/scripts/{f}')

# ---- Zenodo metadata (API schema, for reference/automation) ----
meta={
 "metadata":{
  "upload_type":"dataset",
  "title":"MitoCorex mitochondrial druggability assessment: registry v1.3, supplementary tables and reproduction scripts (CSBR-D-26-00031, R2)",
  "version":"1.3",
  "language":"eng",
  "license":"cc-by-4.0",
  "access_right":"open",
  "creators":[{"name":"Borges, Julian Yin Vieira","orcid":"0009-0001-9929-3135",
               "affiliation":"Department of Computer Science, Boston University"}],
  "keywords":["druggability","mitochondrial targets","binding site","fpocket","ChEMBL",
              "computational pharmacology","drug delivery","MitoCorex"],
  "related_identifiers":[
    {"identifier":"10.5281/zenodo.19389356","relation":"isVersionOf","scheme":"doi"}],
  "description":("Version of record for the second revision (R2) of CSBR-D-26-00031. Adds the "
    "site-matched structural re-analysis: the structural dimension is scored at an a priori "
    "therapeutic site in the isolated chain and the biological assembly context; the D3 dimension "
    "is redefined as a curated pChEMBL bioactivity count; inference uses permutation and bootstrap "
    "methods; and the whole package is generated from a single matrix and verified by an automated "
    "consistency check. Supersedes registry versions 1.0 to 1.2, which must not be cited.")
 }
}
json.dump(meta, open(f'{Z}/metadata.json','w'), indent=2)

# ---- human-readable metadata to type into the web form ----
open(f'{Z}/ZENODO_METADATA.txt','w').write(
"""ZENODO NEW VERSION — fields to enter
====================================
Record (concept DOI): 10.5281/zenodo.19389356   [use "New version" on this record]
Upload type: Dataset
Title:  MitoCorex mitochondrial druggability assessment: registry v1.3, supplementary
        tables and reproduction scripts (CSBR-D-26-00031, R2)
Version: 1.3
License: Creative Commons Attribution 4.0 International (CC BY 4.0)
Access:  Open
Author:  Borges, Julian Yin Vieira  |  ORCID 0009-0001-9929-3135
         Department of Computer Science, Boston University
Keywords: druggability; mitochondrial targets; binding site; fpocket; ChEMBL;
          computational pharmacology; drug delivery; MitoCorex

Description (paste):
Version of record for the second revision (R2) of CSBR-D-26-00031. This version adds the
site-matched structural re-analysis requested at review: the structural dimension is scored at
an a priori therapeutic site, defined before pocket comparison, in the isolated chain and the
biological assembly context. The bioactivity dimension is a curated pChEMBL compound count.
Inference uses permutation tests and percentile bootstrap intervals. The registry, supplementary
tables and scripts are generated from a single quantitative matrix and verified by an automated
consistency check (preflight) before release. Supersedes registry versions 1.0 to 1.2, which
must not be cited.
""")

# ---- README ----
open(f'{Z}/README.md','w').write(
"""# Zenodo deposit v1.3 — CSBR-D-26-00031 (R2)

Version of record for the second revision. Concept DOI 10.5281/zenodo.19389356 (stable across versions);
a new version DOI is minted when this bundle is published.

## Contents
- `data/druggability_registry_v1.3.yaml` — full per-target registry, site-matched R2 layer, supersedes v1.0-v1.2.
- `data/Supplementary_Tables_S1_S8.xlsx` and `data/S1..S8_*.csv` — the supplementary tables cited in the manuscript.
- `data/matrix_v13.csv`, `site_scores_final.csv`, `apriori_sites.csv`, `chembl_filtered.csv`, `provenance.csv`,
  `final_stats.json`, `derived_stats.json`, `weight_sensitivity.json`, `af_confidence.csv` — the single-source
  outputs every reported number derives from.
- `data/carried_forward_from_R1/` — the restricted SMILES compound library (identifiers withheld pending patent)
  and the Crossref reference verification log, unchanged from R1.
- `scripts/` — reproduction and verification: `build_matrix_v13.py`, `derived_stats.py`, `final_stats.py`,
  `build_tables.py`, `build_registry_v13.py`, `build_marked.py`, `preflight.py`, `pkg_preflight.py`,
  `purestats.py` (stdlib statistics), and the retrieval scripts.

## Reproduce
`python3 scripts/final_stats.py && python3 scripts/derived_stats.py && python3 scripts/build_tables.py`
then `python3 scripts/preflight.py` must print PASSED. `purestats.py` uses the standard library only.

## Changelog v1.2 -> v1.3
- Structural dimension scored at an a priori therapeutic site (isolated + biological assembly). Central null holds.
- D3 renamed bioactivity evidence (pChEMBL count). Composite fully specified with a weight Monte Carlo.
- Permutation + bootstrap inference; minimum detectable correlation 0.535.
- Single-source generation with an automated consistency check.
""")

# ---- HOW_TO_DEPOSIT ----
open(f'{Z}/HOW_TO_DEPOSIT.md','w').write(
"""# How to publish this as a new Zenodo version

Do this AFTER the R2 manuscript is submitted in Editorial Manager.

1. Open the existing record: https://doi.org/10.5281/zenodo.19389356
2. Click **New version** (top right). This preserves the concept DOI and mints a new version DOI.
3. Remove the previous version's files if you want a clean version, or keep them; then upload the
   contents of this bundle: the whole `data/` and `scripts/` folders (drag the files in).
4. Set the metadata from `ZENODO_METADATA.txt` (title, version 1.3, license CC BY 4.0, keywords, description).
5. Under **Related/alternate identifiers**, confirm `isVersionOf 10.5281/zenodo.19389356` is present.
6. Save, then **Publish**. Publishing mints the version DOI and cannot be undone (you can only add further versions).
7. The manuscript cites the concept DOI, which does not change, so no manuscript edit is needed after publishing.

Note: publishing a Zenodo version is public and permanent. Review the file list and metadata before you click Publish.
""")

# count
nfiles=sum(len(f) for _,_,f in os.walk(Z))
print(f"Zenodo bundle assembled at {Z}")
print(f"data files: {len(os.listdir(f'{Z}/data'))-1} + carried_forward; scripts: {len(os.listdir(f'{Z}/scripts'))}")
