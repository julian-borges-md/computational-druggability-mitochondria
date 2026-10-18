#!/usr/bin/env python3
"""A priori therapeutic-site definition per target (comment 2).

Method, stated so it is reviewer-defensible and objective:
  Tier D (defined, ligand-anchored): the structure contains the biologically or
     pharmacologically relevant cofactor/ligand. Site = residues within 5 A of it,
     computed from the structure at overlap time. Anchor ligands listed below.
  Tier U (defined, UniProt-anchored): no such ligand in the structure, but UniProt
     annotates active/binding-site residues. Site = those residues (+/- neighbours).
  Tier N (no defined small-molecule site): structural subunit or intrinsically
     disordered. Site-matched structural score = NA (not a low number). This is the
     honest answer the reviewer's comment 8 also points at.
"""
import csv, json, os
S2 = {r['target_id']: r for r in csv.DictReader(open(os.path.expanduser('~/r2/data/SupplementaryTable_S2_fpocket_per_target.csv')))}
UP = json.load(open(os.path.expanduser('~/r2/out/uniprot_sites.json')))

# per-target: tier, anchor ligand (3-letter HET or peptide), site name, basis
SITE = {
 'TGT_001': ('D', 'quinone/UQ region', 'Complex I quinone cavity (ND1)',
             'rotenone and piericidin bind the ND1 Q-site; anchor by the quinone tunnel residues'),
 'TGT_002': ('N', '', 'proton-pumping TMH bundle', 'no defined small-molecule site; structural subunit'),
 'TGT_003': ('D', 'HEM/UQ', 'Qi and Qo ubiquinone sites', 'antimycin/stigmatellin site; anchor by b-hemes and quinone'),
 'TGT_004': ('D', 'HEA/CU', 'binuclear CuB-heme a3 centre', 'catalytic centre; anchor by heme a3 and Cu'),
 'TGT_005': ('D', 'oligomycin region', 'a-subunit proton channel', 'oligomycin binds the c-ring/a interface; anchor by the proton hemichannel residues'),
 'TGT_006': ('D', 'FMN/SF4', 'FMN / NADH site', 'catalytic flavin; anchor by FMN'),
 'TGT_007': ('D', 'FAD', 'FAD / dicarboxylate site', 'succinate/malonate site; anchor by FAD'),
 'TGT_008': ('N', '', 'structural core subunit', 'no defined small-molecule site'),
 'TGT_009': ('N', '', 'allosteric regulatory subunit', 'no validated small-molecule site'),
 'TGT_010': ('D', 'ANP/ADP', 'F1 nucleotide site', 'catalytic/non-catalytic nucleotide site; anchor by bound nucleotide'),
 'TGT_011': ('U', '', 'pol gamma active site', 'polymerase and 3-5 exonuclease catalytic residues (UniProt)'),
 'TGT_012': ('U', '', 'helicase ATP site', 'Walker A/B ATP site (UniProt)'),
 'TGT_013': ('N', '', 'HMG-box DNA interface', 'DNA-binding protein; no small-molecule site'),
 'TGT_014': ('N', '', 'TM import channel', 'membrane channel; no validated small-molecule site'),
 'TGT_015': ('N', '', 'beta-barrel import channel', 'no validated small-molecule site'),
 'TGT_016': ('U', '', 'GTPase site', 'G1-G4 GTP-binding motifs (UniProt)'),
 'TGT_017': ('U', '', 'GTPase site', 'GTP-binding motifs (UniProt)'),
 'TGT_018': ('U', '', 'kinase ATP cleft', 'ATP-binding and catalytic residues (UniProt)'),
 'TGT_019': ('D', 'MN', 'Mn active site', 'superoxide dismutase Mn centre; anchor by Mn'),
 'TGT_020': ('N', '', 'intrinsically disordered coactivator', 'AlphaFold low confidence, disordered; treat structure as NA (also comment 8)'),
 'TGT_021': ('U', '', 'ADP/ATP carrier cavity', 'substrate cavity residues (UniProt); CATR/BKA site'),
 'TGT_022': ('U', '', 'channel pore / NTD', 'pore and N-terminal helix residues (UniProt)'),
 'TGT_023': ('U', '', 'protease + AAA+ ATPase sites', 'Ser-Lys protease dyad and Walker motifs (UniProt)'),
 'TGT_024': ('D', 'Nrf2 peptide / Kelch pocket', 'KEAP1 Kelch central pocket',
             'Nrf2 ETGE/DLG binding pocket in the Kelch beta-propeller; anchor by bound peptide/ligand or the canonical Kelch pocket residues R415 R483 S363 N382 S508 S602 Y525 Y572 Q530'),
 'TGT_025': ('U', '', 'ClpP protease active site', 'Ser153-His178-Asp227 catalytic triad (UniProt)'),
}

rows = []
for tid in sorted(SITE):
    tier, anchor, name, basis = SITE[tid]
    up = UP[tid]
    upres = sorted({f['start'] for f in up['features'] if f['start']})
    s2 = S2[tid]
    rows.append(dict(target_id=tid, gene=up['gene'], assessed_protein=s2['assessed_protein'],
                     uniprot=up['uniprot'], structure=s2['structure'], tier=tier,
                     site_name=name, anchor=anchor,
                     uniprot_site_residues=';'.join(map(str, upres)) if upres else '',
                     basis=basis))
with open(os.path.expanduser('~/r2/out/apriori_sites.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

TIERNAME = {'D': 'defined, ligand-anchored', 'U': 'defined, UniProt-anchored', 'N': 'no small-molecule site (NA)'}
from collections import Counter
c = Counter(r['tier'] for r in rows)
print(f"{'tid':8}{'gene':10}{'str':10}{'tier':5} site")
for r in rows:
    print(f"{r['target_id']:8}{r['gene']:10}{r['structure']:10}{r['tier']:5} {r['site_name']}")
print()
for t in ('D', 'U', 'N'):
    genes = [r['gene'] for r in rows if r['tier'] == t]
    print(f"  Tier {t} ({TIERNAME[t]}): {c[t]:2d}  {', '.join(genes)}")
print('\nwrote out/apriori_sites.csv')
