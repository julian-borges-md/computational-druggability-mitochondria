#!/usr/bin/env python3
"""Build every R2 table (main + supplementary) from the authoritative sources,
so no number is transcribed by hand. Emits CSV per table into out/tables/ and a
single supplementary workbook out/tables/Supplementary_Tables_R2.xlsx.

Sources (single source of truth):
  matrix_v13.csv            composite + per-dimension values
  site_scores_final.csv     isolated / assembly / best site scores
  apriori_sites.csv         a priori therapeutic site definitions
  chembl_filtered.csv       three-count ChEMBL filtering
  final_stats.json          permutation p / bootstrap CI / BH q (7-test family)
  derived_stats.json        pairwise matrix + leave-one-out
  weight_sensitivity.json   weight-robustness membership + schemes
  af_confidence.csv         AlphaFold pLDDT / PAE
"""
import csv, json, os
OUT=os.path.expanduser('~/r2/out'); TB=f'{OUT}/tables'; os.makedirs(TB, exist_ok=True)
def rd(p): return list(csv.DictReader(open(f'{OUT}/{p}')))
def rdj(p): return json.load(open(f'{OUT}/{p}'))
M ={r['target_id']:r for r in rd('matrix_v13.csv')}
SS={r['target_id']:r for r in rd('site_scores_final.csv')}
AP={r['target_id']:r for r in rd('apriori_sites.csv')}
CH={r['target_id']:r for r in rd('chembl_filtered.csv')}
FS=rdj('final_stats.json'); DS=rdj('derived_stats.json'); WS=rdj('weight_sensitivity.json')
AF=rd('af_confidence.csv')
order=sorted(M, key=lambda t:int(M[t]['rank']))
def wcsv(name, header, rows):
    with open(f'{TB}/{name}','w',newline='') as f:
        w=csv.writer(f); w.writerow(header); w.writerows(rows)
    return header, rows

TABLES={}

# ---------- Table 1: five-dimension assessment, all 25, site-matched ----------
h=['Rank','Target','Protein','Pub tier','Site tier','D1 site (isolated)','D1 site (assembly)',
   'D2 site volume A3','D3 bioactivity (pChEMBL)','D4 precedent','D5 delivery','Composite']
rows=[]
for t in order:
    m=M[t]; s=SS[t]
    rows.append([m['rank'], m['gene'], m['assessed_protein'], m['tier_published'], s['tier'],
                 s['site_D1_iso'], s['site_D1_best'], s['site_D2_best'],
                 m['D3_bioactivity_pchembl'], m['D4_precedent'], m['D5_delivery'], m['composite']])
TABLES['Table1_five_dimension']=wcsv('Table1_five_dimension.csv', h, rows)

# ---------- Table 3: Tier-1 ranked, two-stage selection ----------
STAGE2={
 'DNM1L':('PASS','Both dominant-negative missense and loss-of-function variants; fission modulation is therapeutically relevant in both directions, so no single direction is imposed'),
 'NFE2L2':('PASS','KEAP1 antagonism stabilises Nrf2 (gain of antioxidant signalling); largest chemical-matter set in the panel; direction does not depend on NFE2L2 genotype'),
 'LONP1':('PASS','CODAS-syndrome variants are largely missense with reduced proteolytic activity; chaperone/protease modulation, with activator and inhibitor chemotypes noted'),
 'SDHA':('PASS','Leigh syndrome / complex II deficiency variants are loss-of-function; direction is stabilisation or functional rescue; a purely inhibitory precedent would not qualify'),
 'CLPP':('PASS','Perrault syndrome is loss-of-function, so activation is required; an acyldepsipeptide hyperactivator exists, but rescue of a near-null genotype is a hypothesis to test, not an established direction'),
 'MT-CO1':('FAIL','mtDNA-encoded, gene-therapy delivery category, no rescue-directed chemical matter'),
 'POLG':('FAIL','Nucleoside-analogue precedent is inhibition with mitochondrial toxicity, whereas POLG disease requires polymerase rescue'),
 'MT-CYB':('FAIL','Atovaquone binds Plasmodium cytochrome b, not the human protein; human inhibition is an adverse effect. Precedent scored 0'),
 'ATP5F1A':('HOLD','Only 2 potent compounds; the F1 nucleotide site is recovered in the assembly context but chemical matter is thin'),
 'OPA1':('HOLD','GTPase-site selectivity unresolved'),
 'VDAC1':('HOLD','Large channel with multiple sub-sites; no single defined pharmacophore'),
 'PINK1':('HOLD','No PINK1 activator has entered a registered trial'),
 'NDUFV1':('HOLD','Zero potent ChEMBL compounds; prior IACS-010759 attribution withdrawn'),
 'SLC25A4':('HOLD','CATR/BKA cavity is tight; assessed on a high-confidence AlphaFold model (site pLDDT 92)'),
 'SOD2':('HOLD','Zero potent ChEMBL compounds; Mn active site with no drug-like precedent'),
}
h=['Rank','Gene','Composite','D1 site (assembly)','D3 bioactivity (pChEMBL)','D4','D5','Stage 2','Basis']
rows=[]
for t in order:
    if M[t]['tier_published']!='1': continue
    m=M[t]; g=m['gene']; st,basis=STAGE2.get(g,('',''))
    rows.append([m['rank'], g, m['composite'], SS[t]['site_D1_best'],
                 m['D3_bioactivity_pchembl'], m['D4_precedent'], m['D5_delivery'], st, basis])
TABLES['Table3_tier1_ranked']=wcsv('Table3_tier1_ranked.csv', h, rows)

# ---------- Table 4: what the tier label encodes (7-test family) ----------
LABEL={
 'tier vs composite':'Tier vs composite score',
 'tier vs D4 precedent':'Tier vs D4 pharmacological precedent',
 'tier vs D5 delivery':'Tier vs D5 delivery feasibility',
 'tier vs D3 bioactivity':'Tier vs D3 bioactivity evidence',
 'tier vs D1 site druggability':'Tier vs D1 site-matched pocket druggability',
 'D1 site vs D4 precedent (central null)':'D1 site-matched druggability vs D4 precedent (central null)',
 'D1 site vs D3 bioactivity':'D1 site-matched druggability vs D3 bioactivity',
}
h=['Test','n','Spearman rho','permutation p','bootstrap 95% CI','BH q','Survives BH (q<0.05)']
rows=[]
for r in FS:
    ci=r['ci']
    rows.append([LABEL.get(r['test'],r['test']), r['n'], f"{r['rho']:+.3f}", f"{r['perm_p']:.4f}",
                 f"[{ci[0]:+.3f}, {ci[1]:+.3f}]", f"{r['bh_q']:.4f}",
                 'yes' if r['bh_q']<0.05 else 'no'])
TABLES['Table4_tier_encoding']=wcsv('Table4_tier_encoding.csv', h, rows)

# ---------- Supplementary S1: structure provenance + a priori site ----------
PV={r['target_id']:r for r in rd('provenance.csv')}
h=['Target','Gene','Assessed protein','UniProt','PDB / model','Resolution A','Coverage %','Chains kept / entity',
   'Site tier','A priori site','Anchor','Basis']
rows=[[AP[t]['target_id'],AP[t]['gene'],AP[t]['assessed_protein'],AP[t]['uniprot'],
       PV[t]['pdb'],PV[t]['resolution_A'],PV[t]['coverage_pct'],PV[t]['chains_kept_in_entity'],
       AP[t]['tier'],AP[t]['site_name'],AP[t]['anchor'] or 'UniProt active/binding residues',AP[t]['basis']]
      for t in sorted(AP)]
TABLES['S1_provenance_and_sites']=wcsv('S1_provenance_and_sites.csv', h, rows)

# ---------- Supplementary S2: site-matched scores (iso / assembly / best) ----------
h=['Target','Gene','Site tier','D1 isolated','D2 isolated A3','D1 assembly','D2 assembly A3','D1 best','D2 best A3']
rows=[[SS[t]['target_id'],SS[t]['gene'],SS[t]['tier'],SS[t]['site_D1_iso'],SS[t]['site_D2_iso'],
       SS[t]['site_D1_sub'],SS[t]['site_D2_sub'],SS[t]['site_D1_best'],SS[t]['site_D2_best']]
      for t in sorted(SS)]
TABLES['S2_site_matched_scores']=wcsv('S2_site_matched_scores.csv', h, rows)

# ---------- Supplementary S3: pairwise Spearman matrix ----------
dims=['D1','D2','D3','D4','D5']
pm=DS['pairwise']
def cell(a,b):
    if a==b: return '1.000'
    key=f"{a}-{b}" if f"{a}-{b}" in pm else f"{b}-{a}"
    return f"{pm[key]['rho']:+.3f}"
h=['Dimension']+dims
rows=[[a]+[cell(a,b) for b in dims] for a in dims]
sm=DS['pairwise_summary']
rows.append(['off-diagonal range', f"min {sm['min']:+.3f}", f"max {sm['max']:+.3f}", f"mean|.| {sm['mean_abs']:.3f}", '', ''])
TABLES['S3_pairwise_matrix']=wcsv('S3_pairwise_matrix.csv', h, rows)

# ---------- Supplementary S4: leave-one-dimension-out ----------
loo=DS['leave_one_out']
DNAME={'D1':'pocket druggability','D2':'pocket geometry','D3':'bioactivity evidence',
       'D4':'pharmacological precedent','D5':'delivery feasibility'}
h=['Composite variant','Spearman rho vs tier (n=25)','Change vs full']
rows=[['Full composite (all 5 dimensions)', f"{loo['full']:+.3f}", 'reference']]
for d in dims:
    delta=loo[f'drop_{d}']-loo['full']
    direction='strengthens' if delta<0 else ('weakens' if delta>0 else 'unchanged')
    rows.append([f'Drop {d} ({DNAME[d]})', f"{loo['drop_'+d]:+.3f}", f"{delta:+.3f} ({direction})"])
TABLES['S4_leave_one_out']=wcsv('S4_leave_one_out.csv', h, rows)

# ---------- Supplementary S5: full quantitative matrix v1.3 (normalised values) ----------
h=['Rank','Target','Gene','Pub tier','Site tier','D1 norm','D2 norm','D3 norm','D4 norm','D5 norm','Composite','ChEMBL status']
rows=[[M[t]['rank'],M[t]['target_id'],M[t]['gene'],M[t]['tier_published'],M[t]['site_tier'],
       M[t]['D1_norm'],M[t]['D2_norm'],M[t]['D3_norm'],M[t]['D4_norm'],M[t]['D5_norm'],
       M[t]['composite'],M[t]['chembl_status']] for t in order]
TABLES['S5_quantitative_matrix']=wcsv('S5_quantitative_matrix.csv', h, rows)

# ---------- Supplementary S6: weight sensitivity ----------
h=['Gene','Target','Top-5 membership over 20,000 random weightings (%)']
rows=[[m['gene'],m['target_id'],f"{m['pct_top5']:.1f}"] for m in WS['membership']]
rows.append(['','',''])
rows.append([f"Mean Spearman rank correlation vs base weighting: {WS['mean_rank_corr']:.3f}",'',''])
for s in WS['schemes']:
    rows.append([f"First wave, {s['scheme']}", s['first_wave'],''])
TABLES['S6_weight_sensitivity']=wcsv('S6_weight_sensitivity.csv', h, rows)

# ---------- Supplementary S7: ChEMBL 3-count filtering ----------
h=['Target','Gene','UniProt','ChEMBL target','Status','Count <=10 uM','pChEMBL >=5','pChEMBL >=5 binding']
rows=[[CH[t]['target_id'],CH[t]['gene'],CH[t]['uniprot'],CH[t]['chembl_id'],CH[t]['status'],
       CH[t]['n_current'],CH[t]['n_pchembl'],CH[t]['n_pchembl_binding']] for t in sorted(CH)]
TABLES['S7_chembl_filtering']=wcsv('S7_chembl_filtering.csv', h, rows)

# ---------- Supplementary S8: AlphaFold confidence ----------
h=['Target','Gene','Model','UniProt','Global pLDDT','Frac pLDDT>=70','Site pLDDT','n site res','Mean PAE (A)','Interpretation']
rows=[[r['target_id'],r['gene'],r['model'],r['uniprot'],r['global_pLDDT'],r['frac_pLDDT_ge70'],
       r['site_pLDDT'],r['n_site_res'],r['mean_PAE_A'],r['interpretation']] for r in AF]
TABLES['S8_alphafold_confidence']=wcsv('S8_alphafold_confidence.csv', h, rows)

# ---------- assemble supplementary workbook ----------
try:
    from openpyxl import Workbook
    from openpyxl.styles import Font, Alignment
    wb=Workbook(); wb.remove(wb.active)
    SHEETS=[('S1 provenance+sites','S1_provenance_and_sites'),('S2 site-matched','S2_site_matched_scores'),
            ('S3 pairwise','S3_pairwise_matrix'),('S4 leave-one-out','S4_leave_one_out'),
            ('S5 matrix v1.3','S5_quantitative_matrix'),('S6 weight sensitivity','S6_weight_sensitivity'),
            ('S7 ChEMBL filtering','S7_chembl_filtering'),('S8 AlphaFold','S8_alphafold_confidence')]
    for sheet,key in SHEETS:
        header,rows=TABLES[key]; ws=wb.create_sheet(sheet[:31])
        ws.append(header)
        for c in ws[1]: c.font=Font(bold=True)
        for row in rows: ws.append(row)
        for col in ws.columns:
            width=max((len(str(c.value)) for c in col if c.value is not None), default=10)
            ws.column_dimensions[col[0].column_letter].width=min(width+2, 60)
    wb.save(f'{TB}/Supplementary_Tables_R2.xlsx')
    print("wrote Supplementary_Tables_R2.xlsx")
except Exception as e:
    print("xlsx skipped:", e)

print(f"wrote {len(TABLES)} table CSVs to {TB}")
for k in TABLES: print("  ", k)
