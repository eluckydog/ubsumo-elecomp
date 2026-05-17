"""Extract validation sites from DeepPCT + Build enzyme specificity filter.

Task 1: Extract experimentally confirmed Ub/SUMO K sites from DeepPCT training data.
Task 2: Build sequence-based enzyme specificity filter (SUMO ΨKxE motif + Ub E3 patterns).
"""
import csv
import sys
import os
from collections import Counter, defaultdict
import re
import math
import numpy as np

# =========================================================================
# TASK 1: DeepPCT Validation Data Extraction
# =========================================================================

DEEPPCT_TRAIN = "archives/subagent_artifacts/temp_deeppct/datasets/training-set.csv"
DEEPPCT_TEST = "archives/subagent_artifacts/temp_deeppct/datasets/indepentent-test-set.csv"

# Normalize PTM names (case inconsistencies)
PTM_NORMALIZE = {
    "ubiquitination": "Ub",
    "Ubiquitination": "Ub",
    "sumoylation": "SUMO",
    "Sumoylation": "SUMO",
    "SUMOylation": "SUMO",
    "phosphorylation": "Phospho",
    "Phosphorylation": "Phospho",
    "acetylation": "Ac",
    "Acetylation": "Ac",
    "methylation": "Me",
    "Methylation": "Me",
    "O-GlcNAcylation": "OGlcNAc",
    "O-GlcNAc": "OGlcNAc",
}


def parse_site(site_str):
    """Parse site string like 'K163' -> ('K', 163)."""
    aa = site_str[0]
    pos = int(site_str[1:])
    return aa, pos


def extract_ub_sumo_sites(csv_path):
    """Extract all K sites annotated as Ub or SUMO targets.
    
    Returns list of dicts with:
        protein_name, uniprot_id, site_pos, ptm_type, 
        context (what other PTM is involved), relationship,
        evidence (PMID)
    """
    sites = []
    
    with open(csv_path, "r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            ptm1_norm = PTM_NORMALIZE.get(row["PTM1"].strip(), row["PTM1"].strip())
            ptm2_norm = PTM_NORMALIZE.get(row["PTM2"].strip(), row["PTM2"].strip())
            
            # Site1 as Ub/SUMO target
            if ptm1_norm in ("Ub", "SUMO"):
                aa, pos = parse_site(row["Site1"].strip())
                if aa == "K":
                    sites.append({
                        "protein": row["Protein_name"].strip(),
                        "uniprot": row["UniProt ID"].strip(),
                        "k_site": pos,
                        "ptm": ptm1_norm,
                        "context_ptm": ptm2_norm,
                        "context_site": row["Site2"].strip(),
                        "relationship": row["Relationship (PTM1 to PTM2)"].strip(),
                        "evidence": row["PMID"].strip(),
                    })
            
            # Site2 as Ub/SUMO target
            if ptm2_norm in ("Ub", "SUMO"):
                aa, pos = parse_site(row["Site2"].strip())
                if aa == "K":
                    sites.append({
                        "protein": row["Protein_name"].strip(),
                        "uniprot": row["UniProt ID"].strip(),
                        "k_site": pos,
                        "ptm": ptm2_norm,
                        "context_ptm": ptm1_norm,
                        "context_site": row["Site1"].strip(),
                        "relationship": row["Relationship (PTM1 to PTM2)"].strip(),
                        "evidence": row["PMID"].strip(),
                    })
    
    return sites


# Extract sites
train_sites = extract_ub_sumo_sites(DEEPPCT_TRAIN)
test_sites = extract_ub_sumo_sites(DEEPPCT_TEST)
all_sites = train_sites + test_sites

print("=" * 60)
print("TASK 1: DeepPCT Ub/SUMO K-site extraction")
print("=" * 60)
print("Training set Ub/SUMO K sites: %d" % len(train_sites))
print("Test set Ub/SUMO K sites: %d" % len(test_sites))
print("Total: %d" % len(all_sites))

# Breakdown by PTM type
ptm_counts = Counter(s["ptm"] for s in all_sites)
print("  Ub: %d, SUMO: %d" % (ptm_counts.get("Ub", 0), ptm_counts.get("SUMO", 0)))

# Unique proteins
uniprots = set(s["uniprot"] for s in all_sites)
proteins = set(s["protein"] for s in all_sites)
print("Unique proteins (UniProt): %d" % len(uniprots))
print("Unique proteins (name): %d" % len(proteins))

# Top proteins by number of annotated sites
protein_site_counts = Counter(s["protein"] for s in all_sites)
print("\nTop 20 proteins by annotated Ub/SUMO K sites:")
for name, count in protein_site_counts.most_common(20):
    # Get UniProt IDs for this protein
    uids = set(s["uniprot"] for s in all_sites if s["protein"] == name)
    ub_count = sum(1 for s in all_sites if s["protein"] == name and s["ptm"] == "Ub")
    sumo_count = sum(1 for s in all_sites if s["protein"] == name and s["ptm"] == "SUMO")
    print("  %-12s (%s): %d sites (Ub=%d SUMO=%d)" % (name, ",".join(uids), count, ub_count, sumo_count))

# Sites with crosstalk between Ub and SUMO (most relevant for our model)
ub_sumo_crosstalk = [s for s in all_sites 
                     if s["ptm"] == "Ub" and s["context_ptm"] == "SUMO"
                     or s["ptm"] == "SUMO" and s["context_ptm"] == "Ub"]
print("\nUb-SUMO crosstalk sites: %d" % len(ub_sumo_crosstalk))
for s in ub_sumo_crosstalk[:10]:
    print("  %s (%s) K%d: %s <- %s at %s (%s)" % (
        s["protein"], s["uniprot"], s["k_site"],
        s["ptm"], s["context_ptm"], s["context_site"], s["relationship"]))


# =========================================================================
# TASK 2: Enzyme Specificity Filter
# =========================================================================

print("\n" + "=" * 60)
print("TASK 2: Enzyme specificity filter")
print("=" * 60)

# SUMO consensus: ΨKxE/D (Ψ = I, L, V, F, M)
# This is recognized by Ubc9 (SUMO E2)
SUMO_MOTIF_PATTERN = re.compile(r"[ILVFM]K.[DE]")

# Ub E3 motifs (various degrons)
# KEN box: KENxxx(P/S)
# D-box: RxxLxxxxN
# N-end rule: depends on N-terminal residue
# PEST sequence: rich in P, E, S, T
UB_MOTIFS = {
    "KEN_box": re.compile(r"KEN..."),
    "D_box": re.compile(r"R..L....N"),
    "PEST_like": re.compile(r"[PEST]{4,}"),  # Simplified
}


class EnzymeSpecificityFilter:
    """Sequence-based filter for Ub vs SUMO preference.
    
    SUMO E2 (Ubc9) recognizes ΨKxE/D consensus motif.
    Ub E3 ligases recognize various degrons (KEN box, D-box, etc.)
    
    This filter takes a protein sequence window around a target K
    and returns a log-odds score for Ub vs SUMO preference.
    """
    
    def __init__(self):
        self.sumo_consensus = SUMO_MOTIF_PATTERN
    
    def has_sumo_motif(self, seq_window, k_position_in_window):
        """Check if K is in a SUMO consensus ΨKxE/D.
        
        SUMO consensus: hydrophobic(I/L/V/F/M) - K - any - acidic(D/E)
        K must be at position 1 in the 4-residue window (positions 0,1,2,3).
        """
        if len(seq_window) < 4:
            return False
        
        # The K must be at position 1 of the 4-mer
        # Need to check all possible 4-mer windows containing this K
        for i in range(len(seq_window) - 3):
            if i + 1 != k_position_in_window:
                continue
            tetramer = seq_window[i:i+4]
            if self.sumo_consensus.match(tetramer):
                return True
        return False
    
    def score(self, seq_window, k_position_in_window):
        """Return SUMO score (positive = SUMO favored, negative = Ub favored).
        
        Based on:
        - SUMO consensus motif presence
        - Charge context around K (acidic patches favor SUMO)
        - Hydrophobic patch proximity (favors Ub E3 recognition)
        """
        score = 0.0
        
        # 1. SUMO motif: strong positive signal
        if self.has_sumo_motif(seq_window, k_position_in_window):
            score += 2.0
        
        # 2. Local charge: acidic context favors SUMO (mimics ΨKxE/D broader context)
        k_pos = k_position_in_window
        neg_charge_count = sum(1 for i, aa in enumerate(seq_window) 
                              if aa in "DE" and abs(i - k_pos) <= 5)
        pos_charge_count = sum(1 for i, aa in enumerate(seq_window) 
                              if aa in "KRH" and abs(i - k_pos) <= 5 and i != k_pos)
        local_charge = neg_charge_count - pos_charge_count
        score += 0.3 * local_charge
        
        # 3. Hydrophobic context: Ub E3 ligases often recognize hydrophobic patches
        hydro_count = sum(1 for i, aa in enumerate(seq_window) 
                         if aa in "ILVFM" and abs(i - k_pos) <= 5)
        score -= 0.2 * hydro_count  # Ub favored
        
        return score
    
    def get_probability(self, seq_window, k_position_in_window, 
                       electrostatic_delta_e, temperature=2.5):
        """Combine electrostatic + sequence specificity into probability.
        
        Modified softmax:
            P(SUMO) = sigmoid(electrostatic_DE/T + sequence_score)
            P(Ub) = 1 - P(SUMO)
        """
        raw_score = electrostatic_delta_e / temperature + self.score(seq_window, k_position_in_window)
        p_sumo = 1.0 / (1.0 + math.exp(-raw_score))
        return np.clip(p_sumo, 0.01, 0.99)


# Test the filter with known examples
print("\nTest: Enzyme specificity on known Ub/SUMO sites")

filter_ = EnzymeSpecificityFilter()

# PCNA K164: known competition site, SUMO in unstressed
# PCNA sequence around K164: ...LMDLDVEQLG IP EQE...
# Actually let me extract from PDB
def read_sequence_from_pdb(pdb_path):
    """Extract CA residue sequence from PDB."""
    residues = []
    seen = set()
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ENDMDL"):
                break
            if line.startswith("ATOM") and " CA " in line[12:16]:
                res_name = line[17:20].strip()
                res_id = line[22:26].strip()
                key = res_id
                if key not in seen:
                    seen.add(key)
                    # Convert 3-letter to 1-letter
                    aa_map = {
                        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
                        "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
                        "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
                        "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
                        "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
                    }
                    residues.append((res_id, aa_map.get(res_name, "X")))
    return residues

# Analyze PCNA
pcna_seq = read_sequence_from_pdb("data/pdb/pcna_1VYM.pdb")
print("\nPCNA (1VYM) sequence analysis:")
print("  Total residues (chain A): %d" % len([r for r in pcna_seq if r[0].isdigit() and int(r[0]) <= 260]))

# Find K164 flanking
for i, (rid, aa) in enumerate(pcna_seq):
    if rid.strip() == "164" and aa == "K":
        # Get window ±7 residues
        start = max(0, i - 7)
        end = min(len(pcna_seq), i + 8)
        window = "".join(r[1] for r in pcna_seq[start:end])
        k_pos = i - start
        print("  K164 window (±7): %s" % window)
        print("  K164 position in window: %d" % k_pos)
        
        has_s = filter_.has_sumo_motif(window, k_pos)
        s_score = filter_.score(window, k_pos)
        print("  SUMO motif (VKxE/D)? %s" % has_s)
        print("  Sequence score: %+.2f (positive=SUMO)" % s_score)
        
        # Combined with electrostatic
        combined_p = filter_.get_probability(window, k_pos, 3.48, 2.5)
        print("  E+Seq combined P(SUMO): %.3f" % combined_p)
        break

# Check PCNA K117, K138 too (other surface sites)
for target_pos in ["117", "138", "254", "77"]:
    for i, (rid, aa) in enumerate(pcna_seq):
        if rid.strip() == target_pos and aa == "K":
            start = max(0, i - 7)
            end = min(len(pcna_seq), i + 8)
            window = "".join(r[1] for r in pcna_seq[start:end])
            k_pos = i - start
            s_score = filter_.score(window, k_pos)
            has_s = filter_.has_sumo_motif(window, k_pos)
            print("  K%s: motif=%s score=%+.2f  window=%s" % (target_pos, has_s, s_score, window))
            break

print("\nDone!")