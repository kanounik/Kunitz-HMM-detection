## Kunitz Domain Profile HMM

## Table of Contents

- [Overview](#overview)
- [Background](#background)
- [Repository Structure](#repository-structure)
- [Pipeline](#pipeline)
- [Methods](#methods)
- [Results](#results)
- [Dependencies](#dependencies)
- [Data Sources](#data-sources)
- [Project Report](#project-report)
- [Authors](#authors)


## Overview
This project builds and validates a **Profile Hidden Markov Model (pHMM)** for the detection of the **Kunitz-type protease inhibitor domain** (PFAM: PF00014) in protein sequences. The trained model is applied to annotate Kunitz domains across the entire Swiss-Prot database.

### Key goals:
- Build a custom pHMM for the Kunitz domain from structural and sequence data
- Validate the model against manually curated Swiss-Prot annotations using 2-fold cross-validation
- Scan full Swiss-Prot and analyze the distribution of Kunitz domains across organisms

## Background
The **Kunitz domain** is a short (~58 residues), disulfide-rich alpha+beta fold found in serine protease inhibitors. The canonical example is **BPTI** (Bovine Pancreatic Trypsin Inhibitor, PDB: 3TGI, chain I), which tightly binds trypsin via its Lys15 residue. The fold is stabilized by three conserved disulfide bonds (Cys5–Cys55, Cys14–Cys38, Cys30–Cys51).

### Kunitz domains appear in proteins such as:
* Aprotinin (BPTI) — antifibrinolytic drug (Trasylol)
* Alzheimer's amyloid precursor protein (APP)
* Tissue factor pathway inhibitor (TFPI)

## Repository Structure

```text
kunitz-hmm-profile/
│
├── pdb_files/                   # Downloaded PDB structures (26 files)
├── kunitz_chains/               # Extracted Kunitz domain chains
│
├── 3TGI.pdb                     # Canonical BPTI reference structure
│
├── kunitz_with_pdb.fasta        # Sequences with PDB structures
├── kunitz_all.fasta             # All collected Kunitz sequences
├── kunitz_sequences.fasta       # Filtered Kunitz sequences for alignment
├── kunitz_aligned.fasta         # MAFFT multiple sequence alignment
│
├── positives.fasta              # Full Swiss-Prot positive set (368 sequences)
├── negatives.fasta              # Full Swiss-Prot negative set (564,993 sequences)
├── test_positives.fasta         # Balanced test positive set (368 sequences)
├── test_negatives.fasta         # Balanced test negative set (368 sequences)
│
├── kunitz.hmm                   # Trained HMM profile (102 nodes)
│
├── hits_pos.tbl                 # hmmsearch results on positive test set
├── hits_neg.tbl                 # hmmsearch results on negative test set
├── swissprot_all.fasta          # Full Swiss-Prot sequences
├── swissprot_all.tbl            # hmmsearch results on full Swiss-Prot (379 hits)
│
├── download_pdb.py              # Download 26 PDB structures
├── extract_chains.py            # Extract Kunitz chains from PDB files
├── validate_hmm.py              # Run hmmsearch & compute confusion matrix
├── evalue_analysis.py           # E-value threshold analysis & metrics
│
├── environment.yml
├── requirements.txt
└── README.md
```

## Methods
### 1. Structure Selection
26 PDB structures containing Kunitz domains were retrieved using download_pdb.py, starting from the canonical BPTI structure (3TGI, chain I). Kunitz chains were then extracted with extract_chains.py. Sequences with PF00014 annotation were additionally retrieved from UniProt to supplement the training set.
   
### 2. Sequence Alignment & Seed Preparation
The selected Kunitz sequences were aligned using MAFFT, producing kunitz_aligned.fasta. The resulting multiple sequence alignment was used directly as input for HMMER.

### 3. HMM Training
The profile HMM was built using HMMER hmmbuild, yielding a model with 102 nodes:

