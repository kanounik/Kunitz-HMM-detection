## Kunitz Domain Profile HMM

## Table of Contents

- [Overview](#overview)
- [Background](#background)
- [Repository Structure](#repository-structure)
- [Quick Start](#quick-start)
- [Methods](#methods)
- [Results](#results)
- [Dependencies](#dependencies)
- [Data Sources](#data-sources)
- [Project Report](#project-report)
- [Authors](#authors)


## Overview
This project builds and validates a **Profile Hidden Markov Model (pHMM)** for the detection of the Kunitz-type protease inhibitor domain (PFAM: PF00014) in protein sequences. The trained model is then applied to annotate Kunitz domains across the entire Swiss-Prot database.

### Key goals:
- Build a custom pHMM for the Kunitz domain from structural data
- Validate the model against manually curated Swiss-Prot annotations
- Annotate Kunitz domains across Swiss-Prot and analyze their distribution

## Background
The Kunitz domain is a short (~58 residues), disulfide-rich alpha+beta fold found in serine protease inhibitors. The canonical example is BPTI (Bovine Pancreatic Trypsin Inhibitor, PDB: 3TGI, chain I), which tightly binds trypsin via its Lys15 residue. The fold is stabilized by three conserved disulfide bonds (Cys5–Cys55, Cys14–Cys38, Cys30–Cys51).

### Kunitz domains appear in proteins such as:
* Aprotinin (BPTI) — antifibrinolytic drug (Trasylol)
* Alzheimer's amyloid precursor protein (APP)
* Tissue factor pathway inhibitor (TFPI)

## Repository Structure

```text
kunitz-hmm-profile/
│
├── data/
│   ├── structures/              # 26 downloaded PDB files
│   ├── sequences/               # FASTA files (positive/negative sets)
│   ├── alignments/
│   │   └── kunitz_aligned.fasta # MAFFT multiple sequence alignment
│   └── swissprot/               # Swiss-Prot subsets used for testing
│
├── models/
│   └── kunitz.hmm               # Trained HMM profile (102 nodes)
│
├── results/
│   ├── hmmsearch_train.out      # Consistency test: 26/26 recovered
│   ├── hmmsearch_test.out       # Validation results (368 pos + 368 neg)
│   ├── hmmsearch_swissprot.out  # Full Swiss-Prot scan (379 hits)
│   ├── confusion_matrix.tsv     # TP=363, FP=0, TN=368, FN=5
│   ├── performance.tsv          # Metrics per E-value threshold
│   └── figures/                 # ROC curve, score distributions
│
├── scripts/
│   ├── download_pdb.py          # Download 26 PDB structures
│   ├── extract_chains.py        # Extract Kunitz chains from PDB files
│   ├── 03_build_hmm.sh          # Run hmmbuild
│   ├── 04_search_hmm.sh         # Run hmmsearch on datasets
│   ├── 05_evaluate.py           # Confusion matrix & metrics
│   └── 06_annotate_swissprot.py # Apply model to full Swiss-Prot
│
├── report/
│   └── kunitz_hmm_report.pdf    # Final project report
│
├── environment.yml
├── requirements.txt
└── README.md
