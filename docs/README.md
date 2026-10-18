# Zenodo deposit v1.3 — CSBR-D-26-00031 (R2)

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
