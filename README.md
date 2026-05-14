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
├── data/
│   ├── pdb_files/                    # 26 downloaded PDB structures
│   ├── kunitz_chains/                # Extracted Kunitz domain chains
│   ├── 3TGI.pdb                      # Canonical BPTI reference structure
│   ├── kunitz_with_pdb.fasta         # Kunitz sequences with PDB structures
│   ├── kunitz_all.fasta              # All collected Kunitz sequences
│   ├── kunitz_sequences.fasta        # Filtered sequences used for alignment
│   └── kunitz_aligned.fasta          # MAFFT multiple sequence alignment
│
├── models/
│   └── kunitz.hmm                    # Trained HMM profile (102 nodes)
│
├── results/
│   ├── positives_ids.tsv             # IDs + lengths of 368 positive sequences (PF00014)
│   ├── negatives_ids.tsv             # IDs + lengths of 368 sampled negative sequences
│   ├── hits_pos.tsv                  # hmmsearch results on positive test set (363 hits)
│   ├── hits_neg.tsv                  # hmmsearch results on negative test set (no hits)
│   ├── confusion_matrix.tsv          # TP=363, FP=0, TN=368, FN=5
│   ├── performance_metrics.tsv       # Sensitivity, Specificity, Precision, F1, MCC
│   ├── evalue_analysis.tsv           # Metrics across 9 E-value thresholds
│   ├── swissprot_hits.tsv            # 379 Kunitz hits across full Swiss-Prot
│   └── figures/
│       ├── confusion_matrices.png    # 2-fold cross-validation confusion matrices
│       ├── roc_curve.png             # ROC curve (AUC = 0.995, both folds)
│       ├── mcc_thresholds.png        # MCC vs E-value threshold
│       └── hmm_logo.png              # HMM sequence logo (102 positions)
│
├── scripts/
│   ├── download_pdb.py               # Download 26 PDB structures
│   ├── extract_chains.py             # Extract Kunitz chains from PDB files
│   ├── validate_hmm.py               # hmmsearch + confusion matrix + metrics → TSV
│   ├── evalue_analysis.py            # E-value threshold analysis → TSV
│   ├── cross_validation.py           # 2-fold cross-validation
│   ├── roc_curve.py                  # ROC curve generation (AUC = 0.995)
│   └── swissprot.py                  # Scan full Swiss-Prot → TSV
│
├── README.md
├── requirements.txt
└── environment.yml
```

## Methods
