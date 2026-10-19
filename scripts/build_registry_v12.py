#!/usr/bin/env python3
"""Druggability registry v1.2.

Supersedes v1.1, which was assembled before the adversarial second pass and
therefore carries:
  - MT-ND1 fpocket output computed on the entire respiratory Complex I
    (648 pockets, 7837.219 A3) rather than on the MT-ND1 chain (27, 1773.902)
  - pre-replacement structures for TFAM, SOD2 and ATP5F1A
  - precedent scores for MT-CYB and KEAP1 that breach the manuscript's own rule
  - a hand assigned isolation_artefact_risk flag that no script produced
  - a tier concordance statistic (rho -0.078, p 0.711) computed on the
    contaminated matrix
"""
import csv, json, os, yaml, collections

BASE = os.path.expanduser('~/csbrd')
reg  = yaml.safe_load(open(f'{BASE}/out/druggability_registry_v1.1.yaml'))
v2   = {r['target_id']: r for r in csv.DictReader(open(f'{BASE}/out/TableS2_fpocket_v2.csv'))}
mx   = {r['target_id']: r for r in csv.DictReader(open(f'{BASE}/out/TableS5_quantitative_matrix_v2.csv'))}
cov  = {r['target_id']: r for r in csv.DictReader(open(f'{BASE}/out_v2/coverage_v2.csv'))}
st   = json.load(open(f'{BASE}/out_v2/stats_v2.json'))
fw   = json.load(open(f'{BASE}/work/first_wave_v2.json'))
S2W  = {r['stage2'] if False else r['gene']: r for r in fw['tier1']}

def f(x):
    try: return float(x)
    except (TypeError, ValueError): return None
def i(x):
    try: return int(x)
    except (TypeError, ValueError): return None

# ---- structure replacements made after the second pass -------------------
REPLACED = {
 'TGT_013': ('3FGH', '3TMM', 'HMG box 2 only, 27 per cent UniProt coverage; it also anchored the D1 '
                             'normalisation floor. Replaced with the full length protein.'),
 'TGT_019': ('9NSJ', '1LUV', 'Q143N active site mutant expressed in Escherichia coli. Replaced with the '
                             'human wild type enzyme.'),
 'TGT_010': ('8H9E chains A;B', '8H9E chain A',
             'The two retained chains lie 58 A apart; the beta subunits carrying every nucleotide site had '
             'been deleted from between them. Reassessed on a single alpha subunit.'),
}
PREC_FIX = {
 'TGT_003': dict(was=3, now=0,
   reason='Atovaquone binds Plasmodium falciparum cytochrome b. Methods 2.4 states that a compound acting '
          'on a pathogen orthologue is not precedent for the human protein. Section 3.5 of the submitted '
          'manuscript said so while the table scored it 3.'),
 'TGT_024': dict(was=3, now=1,
   reason='Omaveloxolone acts at Cys151 in the KEAP1 BTB domain. The assessed structure 8XGK spans residues '
          '322 to 612, the Kelch domain. The FDA label states the mechanism is not fully established. No '
          'Kelch domain protein protein interaction inhibitor is approved in any jurisdiction.'),
}

for tid, rec in reg['targets'].items():
    r, m = v2[tid], mx[tid]

    # ---- fpocket block, both arms, no hand assigned flag -----------------
    rec['fpocket_v4_0'] = collections.OrderedDict([
        ('executed', 'fpocket 4.0 (Le Guilloux et al. 2009), default parameters, re-executed 1 September 2026'),
        ('determinism_check', 'three independent runs on an identical input file return identical output'),
        ('structure', r['structure']),
        ('chains_retained', r['chains_retained']),
        ('n_chains_in_entity', i(r['n_chains_in_entity'])),
        ('n_chains_retained', i(r['n_chains_retained'])),
        ('chain_cap', 'at most two chains of the matching polymer entity are retained; where the entity '
                      'carries more, this is a deliberate and disclosed truncation'),
        ('apo', collections.OrderedDict([
            ('preparation', 'ligands, cofactors, metals, lipids and waters removed'),
            ('n_heteroatoms_removed', i(r['n_heteroatoms_removed'])),
            ('n_pockets_detected', i(r['apo_n_pockets'])),
            ('rank1_druggability_score', f(r['apo_rank1_drugscore'])),
            ('rank1_pocket_volume_A3', f(r['apo_rank1_volume'])),
            ('max_druggability_score', f(r['apo_max_drugscore'])),
            ('max_druggability_pocket_rank', i(f(r['apo_max_pocket_rank']))),
            ('max_pocket_volume_A3', f(r['apo_max_volume'])),
            ('max_pocket_apolar_sasa_A2', f(r['apo_max_apolar_sasa'])),
            ('max_pocket_polarity', f(r['apo_max_polarity'])),
        ])),
        ('holo', collections.OrderedDict([
            ('preparation', 'waters removed only; ligands, cofactors, metals and lipids retained'),
            ('n_pockets_detected', i(r['holo_n_pockets'])),
            ('rank1_druggability_score', f(r['holo_rank1_drugscore'])),
            ('rank1_pocket_volume_A3', f(r['holo_rank1_volume'])),
            ('max_druggability_score', f(r['holo_max_drugscore'])),
            ('max_druggability_pocket_rank', i(f(r['holo_max_pocket_rank']))),
            ('max_pocket_volume_A3', f(r['holo_max_volume'])),
            ('max_pocket_apolar_sasa_A2', f(r['holo_max_apolar_sasa'])),
            ('max_pocket_polarity', f(r['holo_max_polarity'])),
        ])),
        ('delta_drugscore_apo_minus_holo', f(r['delta_drugscore_apo_minus_holo'])),
        ('analysis_set', 'all 25 targets enter every reported comparison; no target is excluded on any '
                         'structural ground'),
        ('full_pocket_listing', 'SupplementaryTable_S3_fpocket_all_pockets.csv'),
    ])

    # ---- structural evidence -------------------------------------------
    se = rec['structural_evidence']
    se['assessed_structure']   = r['structure']
    se['chains_assessed']      = r['chains_retained']
    se['resolution_angstrom']  = f(r['resolution_A'])
    se['organism']             = r['organism'] or 'Homo sapiens'
    se['mutation']             = r['mutation'] or 'none recorded in the entry'
    se['fragment']             = r['fragment'] or 'not stated in the entry'
    cv = cov[tid]
    se['uniprot_coverage'] = collections.OrderedDict([
        ('uniprot_sequence_length', i(cv['uniprot_length'])),
        ('uniprot_range_declared_by_entry', cv['uniprot_range_modelled']),
        ('residues_resolved_per_retained_chain', cv['n_residues_per_chain']),
        ('coverage_percent', f(cv['coverage_percent'])),
        ('basis', cv['coverage_basis']),
        ('partial', f(cv['coverage_percent']) < 60.0),
    ])
    if tid in REPLACED:
        was, now, why = REPLACED[tid]
        se['replaced_after_review'] = {'was': was, 'now': now, 'reason': why}

    # ---- precedent corrections -----------------------------------------
    if tid in PREC_FIX:
        fx = PREC_FIX[tid]
        rec['precedent_score_correction'] = {
            'score_in_v1_1': fx['was'], 'score_in_v1_2': fx['now'], 'reason': fx['reason']}

    # ---- quantitative scores -------------------------------------------
    rec['quantitative_scores'] = collections.OrderedDict([
        ('D1_pocket_druggability', {'apo_max': f(m['D1_apo_max']), 'holo_max': f(m['D1_holo_max']),
                                    'apo_rank1': f(m['D1_rank1']), 'normalised_apo_max': f(m['D1_norm'])}),
        ('D2_pocket_geometry',     {'apo_max_volume_A3': f(m['D2_apo_max_volume_A3']),
                                    'normalised_log10': f(m['D2_norm'])}),
        ('D3_ligand_efficiency',   {'chembl_target': m['D3_chembl_target'],
                                    'n_compounds_at_or_below_10uM': i(m['D3_chembl_compounds']),
                                    'class_published': m['D3_class_published'],
                                    'class_observed': m['D3_class_observed'],
                                    'reclassified': m['D3_class_published'] != m['D3_class_observed'],
                                    'normalised_log10': f(m['D3_norm'])}),
        ('D4_pharmacological_precedent', {'score': i(m['D4_precedent']), 'score_in_v1_1': i(m['D4_precedent_v1']),
                                    'basis': m['D4_basis'], 'normalised': f(m['D4_norm'])}),
        ('D5_delivery_feasibility',{'score': i(m['D5_delivery']), 'category': m['D5_category'],
                                    'normalised': f(m['D5_norm'])}),
        ('weights', {'D1': 0.25, 'D2': 0.15, 'D3': 0.20, 'D4': 0.20, 'D5': 0.20}),
        ('composite', {'apo_max': f(m['composite_apo']), 'holo_max': f(m['composite_holo']),
                       'apo_rank1': f(m['composite_rank1'])}),
        ('rank_of_25', {'apo_max': i(m['rank_apo']), 'holo_max': i(m['rank_holo']),
                        'apo_rank1': i(m['rank_rank1'])}),
    ])
    g = rec['gene_symbol']
    if g in S2W:
        rec['first_wave_selection'] = {'stage_2_verdict': S2W[g]['stage2'], 'basis': S2W[g]['basis']}

# ---- metadata -----------------------------------------------------------
md = reg['metadata']
md['document_id'] = 'DRUGGABILITY_REGISTRY_MITOCHONDRIAL_v1.2'
md['version'] = '1.2'
md['supersedes'] = 'druggability_registry_v1.1.yaml (version 1.1, 2026-08-23), which must not be cited'
md['date_created'] = '2026-09-11'
md['assessment_framework'][0] = ('fpocket 4.0 (Le Guilloux et al. 2009), default parameters, executed '
                                 '23 August 2026 and re-executed on corrected inputs 1 September 2026')
md['changes_from_v1_1'] = [
 'MT-ND1 fpocket output in v1.1 described the entire respiratory Complex I, not the MT-ND1 chain. An '
 'abandoned script wrote the full complex to the pipeline path and left fpocket running on it; the pipeline '
 'later overwrote the file with the correctly extracted chain but skipped re-execution because an output '
 'directory already existed. Reported values were 648 pockets and a 7837.219 A3 cavity. The correct values '
 'are 27 pockets and a 1773.902 A3 maximum pocket. Scope of the defect is MT-ND1 only; all 25 targets were '
 're-executed and every other target reproduces exactly.',
 'Ligand, cofactor, metal and lipid removal before pocket detection is now disclosed, and a holo arm '
 'retaining them is reported alongside the apo arm for every target. The effect is bidirectional: removal '
 'inflates the score for SDHA by 0.244 and deflates it for TOMM40 by 0.380.',
 'TFAM reassessed on 3TMM (full protein) rather than 3FGH (HMG box 2 only). SOD2 reassessed on 1LUV (human '
 'wild type) rather than 9NSJ (Q143N mutant expressed in Escherichia coli). ATP5F1A reassessed on 8H9E '
 'chain A rather than chains A and B, which lie 58 A apart.',
 'MT-CYB precedent corrected from 3 to 0 and KEAP1 from 3 to 1; both breached the rule stated in Methods 2.4.',
 'The isolation_artefact_risk field is removed. It was assigned by hand, no script produced it, and it did '
 'not correlate with newly exposed surface area. Three pocket selection specifications are reported instead.',
 'Tier concordance recomputed on the corrected matrix: rho -0.181, p 0.387 (v1.1 reported -0.078, p 0.711).',
 'Every correlation now carries a 95 per cent confidence interval and the panel carries a minimum '
 'detectable effect size.',
 'Per target UniProt sequence coverage computed and carried: reference length from UniProt, the residue '
 'range the PDB entry declares, residues actually resolved in each retained chain, and the resulting '
 'coverage percentage. Four assessments fall below 60 per cent: LONP1 19.8 (proteolytic domain only), '
 'OPA1 34.4, DNM1L 46.5 and KEAP1 46.6. Three of these four carry first wave targets, so the coverage '
 'limitation is stated wherever the first wave is reported.',
 'The construct fragment field no longer defaults to "full construct" where the PDB entry leaves '
 'pdbx_fragment empty. An absent annotation is now reported as "not stated in the entry"; the coverage '
 'block carries the substantive answer.',
]
md['analysis_set'] = ('All 25 targets enter every reported comparison. No exclusion rule is applied. '
                      'Methods 2.2 of the first revision stated that six targets were excluded; that '
                      'exclusion was never implemented and the statement has been withdrawn.')

# ---- summary ------------------------------------------------------------
def get(name):
    return next(t for t in st['tests'] if t['test'] == name)
reg['druggability_summary']['fpocket_vs_tier_concordance'] = collections.OrderedDict([
 ('specification', 'apo maximum Druggability Score pocket, all 25 targets'),
 ('spearman_rho_tier_vs_max_druggability', get('tier vs pocket druggability (D1)')['rho']),
 ('p_value', get('tier vs pocket druggability (D1)')['p']),
 ('ci_95', get('tier vs pocket druggability (D1)')['ci']),
 ('benjamini_hochberg_q', get('tier vs pocket druggability (D1)')['bh_q']),
 ('sensitivity', {'holo_max': -0.159, 'apo_rank1': -0.209}),
 ('minimum_detectable_rho_at_n25_power_80', st['mdes_25']),
 ('interpretation',
  'No association between pocket druggability and the assigned tier is detectable at this sample size. '
  'The confidence interval spans -0.538 to +0.231, so the data remain compatible with a moderate '
  'association in either direction; this is a null result, not evidence of absence. Tier assignment is '
  'dominated by pharmacological precedent (rho -0.663, q 0.0012), which survives multiplicity control.'),
])
reg['druggability_summary']['withdrawn_finding'] = collections.OrderedDict([
 ('claim', 'Pocket druggability and delivery feasibility are inversely related among mitochondrial targets.'),
 ('as_reported_in_first_revision', 'rho -0.43, p 0.032'),
 ('corrected', {'apo_max': -0.343, 'p': 0.0931, 'benjamini_hochberg_q': 0.149,
                'holo_max': -0.289, 'apo_rank1': 0.172}),
 ('status', 'WITHDRAWN. Not significant under any specification, does not survive multiplicity control, '
            'and the sign reverses under the rank 1 pocket rule. The first revision value depended on the '
            'contaminated MT-ND1 record and on selecting the maximum scoring pocket.'),
])
reg['druggability_summary']['multiplicity_control'] = {
 'method': 'Benjamini-Hochberg, 8 tests, q < 0.05',
 'tests': st['tests'],
}
reg['druggability_summary']['statistical_power'] = {
 'n': 25, 'alpha': 0.05, 'target_power': 0.80,
 'minimum_detectable_rho': st['mdes_25'],
 'minimum_detectable_rho_bonferroni': st['mdes_bonf'],
 'power_at_true_rho': {0.2: 0.16, 0.3: 0.31, 0.4: 0.51, 0.5: 0.73},
}
reg['druggability_summary']['first_wave'] = {
 'selected': fw['first_wave'],
 'basis': 'two stage selection: composite rank, then a direction of effect and chemical matter check',
 'structural_coverage_caveat': (
   'Three of the five first wave targets are assessed on partial constructs: LONP1 at 19.8 per cent of the '
   'UniProt sequence (proteolytic domain only), DNM1L at 46.5 per cent and NFE2L2 via KEAP1 at 46.6 per '
   'cent. Their D1 and D2 values describe the resolved portion of the protein, not the whole protein.'),
}
reg['druggability_summary']['structural_coverage'] = {
 'method': 'residues resolved in the single best retained chain, divided by the UniProt reference length',
 'uniprot_reference_retrieved': '11 September 2026',
 'below_60_percent': {cov[t]['gene']: float(cov[t]['coverage_percent'])
                      for t in sorted(cov) if float(cov[t]['coverage_percent']) < 60.0},
 'full_table': 'SupplementaryTable_S12_structural_coverage.csv',
}

out = f'{BASE}/out_v2/SupplementaryData_S1_druggability_registry_v1.2.yaml'
with open(out, 'w') as fh:
    fh.write('# Mitochondrial druggability registry, version 1.2\n')
    fh.write('# Supersedes druggability_registry_v1.1.yaml. Manuscript CSBR-D-26-00031.\n')
    yaml.safe_dump(json.loads(json.dumps(reg, default=str)), fh,
                   sort_keys=False, width=110, default_flow_style=False, allow_unicode=True)
print('wrote', out, os.path.getsize(out), 'bytes')
chk = yaml.safe_load(open(out))
print('targets:', len(chk['targets']))
t1 = chk['targets']['TGT_001']['fpocket_v4_0']
print('TGT_001 apo pockets:', t1['apo']['n_pockets_detected'], 'max vol:', t1['apo']['max_pocket_volume_A3'])
print('648 present anywhere:', '648' in open(out).read().replace('0.648',''))
print('7837 present anywhere:', '7837.219' in open(out).read() and 'in changes_from_v1_1 only')
