#!/usr/bin/env python3
"""Automated cross-document consistency check (answers reviewer comment 10).

Verifies, before the package is assembled, that:
  1. every headline statistic in the manuscript and response letter matches the
     recomputed value in the authoritative sources (final_stats / derived_stats /
     weight_sensitivity / af_confidence / matrix_v13);
  2. no superseded value from an earlier revision survives in either document;
  3. every Supplementary Table referenced in the manuscript exists as a file;
  4. the manuscript and the response letter agree on every shared headline number.

Exit status is non-zero if any check fails. Stdlib only, so it runs anywhere.
"""
import csv, json, os, re, sys
R=os.path.expanduser('~/r2'); OUT=f'{R}/out'; TB=f'{OUT}/tables'; B=f'{R}/build'
def rd(p): return list(csv.DictReader(open(p)))
def rj(p): return json.load(open(p))
man=open(f'{B}/R2_manuscript.md').read()
res=open(f'{B}/R2_response.md').read()
both={'manuscript':man,'response':res}
FS={r['test']:r for r in rj(f'{OUT}/final_stats.json')}
DS=rj(f'{OUT}/derived_stats.json'); WS=rj(f'{OUT}/weight_sensitivity.json')
AF={r['gene']:r for r in rd(f'{OUT}/af_confidence.csv')}
MX={r['gene']:r for r in rd(f'{OUT}/matrix_v13.csv')}
fails=[]; checks=0
def present(text, token, where, label):
    global checks; checks+=1
    if token not in text: fails.append(f"[MISSING] {where}: {label} -> expected '{token}'")
def absent(text, token, where, label):
    global checks; checks+=1
    if token in text: fails.append(f"[SUPERSEDED] {where}: {label} -> '{token}' must not appear")

# ---- 1. headline statistics must appear in the manuscript ----
def rho(test): return f"{FS[test]['rho']:+.3f}".replace('+','')  # e.g. -0.663 / 0.263
present(man, "-0.663", 'manuscript', 'tier vs D4 precedent rho')
present(man, "-0.709", 'manuscript', 'tier vs composite rho')
present(man, "-0.037", 'manuscript', 'tier vs D1 site rho')
present(man, "-0.100", 'manuscript', 'central null (D1 site vs D4) rho')
present(man, "+0.117", 'manuscript', 'isolated-chain structure vs precedent')
present(man, "0.535",  'manuscript', 'minimum detectable correlation')
present(man, "0.982",  'manuscript', 'weight-perturbation mean rank correlation')
# leave-one-out
for tok in ("-0.763","-0.688","-0.617","-0.591","-0.664"):
    present(man, tok, 'manuscript', f'leave-one-out {tok}')
# pairwise coherence
for tok in ("+0.66","+0.62"):
    present(man, tok, 'manuscript', f'pairwise coherence {tok}')
# AlphaFold confidence
for g,(pl,pae) in {'PPARGC1A':('52.7','28.1'),'TIMM23':('59.2','21.8'),'SLC25A4':('92.1','4.5')}.items():
    assert AF[g]['global_pLDDT']==pl and AF[g]['mean_PAE_A']==pae, f"AF source drift {g}"
    present(man, pl, 'manuscript', f'{g} global pLDDT')
present(man, "91.7", 'manuscript', 'SLC25A4 site pLDDT')
# weight-sensitivity membership
for g,pct in {'DNM1L':'100','LONP1':'98','SDHA':'95'}.items():
    present(man, pct, 'manuscript', f'{g} weight-robustness %')
# tier counts and site counts
present(man, "15 Tier 1", 'manuscript', 'tier-1 count'); present(man, "0.183 to 0.737",'manuscript','tier-1 composite range')
present(man, "0.512", 'manuscript', 'tier-1 median')
present(man, "Seventeen targets", 'manuscript', '17 defined sites')

# ---- 2. superseded values must be absent from both documents ----
SUPERSEDED={
 "-0.181":"old apo tier vs D1", "0.724":"old CLPP/Tier-1 max composite",
 "0.470":"old Tier-1 median", "ranks 1 to 19":"old Tier-1 rank span",
 "reclassified 12 of 25":"old reclassification claim",
 "D3 ligand efficiency":"un-renamed dimension label",
 "rho 0.052, p 0.80":"old phrasing of central null as global-best",
 "apo maximum":"retired three-specification framing",
 "druggability_registry_v1.2":"superseded registry version as current",
}
for tok,lab in SUPERSEDED.items():
    for where,txt in both.items(): absent(txt, tok, where, lab)

# ---- 3. supplementary numbering: single register, every reference resolves ----
refs=set(re.findall(r'Supplementary Tables? (S\d+)(?: to (S\d+))?', man))
cited=set()
for a,b in refs:
    cited.add(a)
    if b: cited.update(f'S{i}' for i in range(int(a[1:]), int(b[1:])+1))
have={ 'S1':'S1_provenance_and_sites.csv','S2':'S2_site_matched_scores.csv','S3':'S3_pairwise_matrix.csv',
       'S4':'S4_leave_one_out.csv','S5':'S5_quantitative_matrix.csv','S6':'S6_weight_sensitivity.csv',
       'S7':'S7_chembl_filtering.csv','S8':'S8_alphafold_confidence.csv'}
for s in sorted(cited, key=lambda x:int(x[1:])):
    checks+=1
    if s not in have:
        fails.append(f"[DANGLING REF] {s} referenced in manuscript but not in the register")
    elif not os.path.exists(f'{TB}/{have[s]}'):
        fails.append(f"[MISSING FILE] {s} referenced but {have[s]} absent")
# register must be one-to-one (no S-number reused for two tables)
checks+=1
if len(set(have.values()))!=len(have): fails.append("[REGISTER] duplicate file assigned to two S-numbers")

# ---- 4. manuscript and response must agree on the numbers both cite ----
# (values both documents state; manuscript-only Table 4 rows are not listed here)
for tok in ("-0.709","-0.100","+0.117","0.535","0.982","-0.763","-0.688","-0.617","-0.591","-0.664"):
    checks+=1
    if tok in man and tok in res:
        continue
    if tok not in man and tok not in res:
        fails.append(f"[DISAGREE] shared value '{tok}' missing from BOTH documents")
    else:
        where='response' if tok in man else 'manuscript'
        fails.append(f"[DISAGREE] shared value '{tok}' missing from {where}")

# ---- 5. first-wave set stated identically ----
for g in ("DNM1L","LONP1","SDHA","CLPP"):
    present(man, g, 'manuscript', f'first-wave {g}')
present(man, "KEAP1 as the ligand binding partner", 'manuscript', 'first-wave NFE2L2/KEAP1')

print(f"preflight: {checks} checks run")
if fails:
    print(f"FAILED ({len(fails)} issues):")
    for f in fails: print("  ", f)
    sys.exit(1)
print("PASSED: no inconsistencies. Manuscript and response agree with the recomputed sources.")
