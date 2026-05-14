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
The **Kunitz domain** is a short (~58 residues), disulfide-rich alpha+beta fold found in serine protease inhibitors. The canonical example is **BPTI** (Bovine Pancreatic Trypsin Inhibitor, PDB: [`3TGI.pdb`](data/3TGI.pdb), chain I), which tightly binds trypsin via its Lys15 residue. The fold is stabilized by three conserved disulfide bonds (Cys5–Cys55, Cys14–Cys38, Cys30–Cys51).

### Kunitz domains appear in proteins such as:
* Aprotinin (BPTI) — antifibrinolytic drug (Trasylol)
* Alzheimer's amyloid precursor protein (APP)
* Tissue factor pathway inhibitor (TFPI)

## Repository Structure

```text
kunitz-hmm-profile/
│
├── pdb_files/                         # 26 downloaded PDB structures
├── kunitz_chains/                     # Extracted Kunitz domain chains
│
├── [kunitz_with_pdb.fasta](kunitz_with_pdb.fasta)   # Kunitz sequences with PDB structures
├── [kunitz_all.fasta](kunitz_all.fasta)             # All collected Kunitz sequences
├── [kunitz_sequences.fasta](kunitz_sequences.fasta) # Filtered sequences used for alignment
├── [kunitz_aligned.fasta](kunitz_aligned.fasta)     # MAFFT multiple sequence alignment
├── [kunitz.hmm](kunitz.hmm)                         # Trained HMM profile (102 nodes)
├── [3TGI.pdb](3TGI.pdb)                             # Canonical BPTI reference structure
│
├── [positives.fasta](positives.fasta)               # Swiss-Prot positive set (368, PF00014)
├── negatives.fasta                                  # Swiss-Prot negative set (not included)
├── [test_positives.fasta](test_positives.fasta)     # Balanced test positive set
├── [test_negatives.fasta](test_negatives.fasta)     # Balanced test negative set
│
├── results/
│   ├── [hits_pos.tsv](results/hits_pos.tsv)
│   ├── [hits_neg.tsv](results/hits_neg.tsv)
│   ├── [confusion_matrix.tsv](results/confusion_matrix.tsv)
│   ├── [performance_metrics.tsv](results/performance_metrics.tsv)
│   ├── [evalue_analysis.tsv](results/evalue_analysis.tsv)
│   ├── [swissprot_hits.tsv](results/swissprot_hits.tsv)
│   └── figures/
│       ├── confusion_matrices.png
│       ├── roc_curve.png
│       ├── mcc_thresholds.png
│       └── hmm_logo.png
│
├── scripts/
│   ├── [download_pdb.py](scripts/download_pdb.py)
│   ├── [extract_chains.py](scripts/extract_chains.py)
│   ├── [validate_hmm.py](scripts/validate_hmm.py)
│   ├── [evalue_analysis.py](scripts/evalue_analysis.py)
│   ├── [cross_validation.py](scripts/cross_validation.py)
│   ├── [roc_curve.py](scripts/roc_curve.py)
│   └── [scan_swissprot.py](scripts/scan_swissprot.py)
│
├── [README.md](README.md)
├── [requirements.txt](requirements.txt)
└── [environment.yml](environment.yml)
```

## Methods
### 1. Structure Selection
26 PDB structures containing Kunitz domains were retrieved using download_pdb.py, starting from the canonical BPTI structure (3TGI, chain I). Kunitz chains were then extracted with extract_chains.py. Sequences with PF00014 annotation were additionally retrieved from UniProt to supplement the training set.
   
### 2. Sequence Alignment & Seed Preparation
The selected Kunitz sequences were aligned using MAFFT, producing kunitz_aligned.fasta. The resulting multiple sequence alignment was used directly as input for HMMER.

### 3. HMM Training
The profile HMM was built using HMMER hmmbuild, yielding a model with 102 nodes:

