# Kunitz Domain Profile HMM

## Table of Contents

- [Overview](#overview)
- [Background](#background)
- [Repository Structure](#repository-structure)
- [Workflow](#workflow)
- [Results](#results)
- [Swiss-Prot Annotation](#swiss-prot-annotation)
- [Technologies](#technologies)
- [Authors](#authors)


## Overview
A structurally informed Profile Hidden Markov Model (pHMM) for detecting the Kunitz-type protease inhibitor domain (PF00014) in protein sequences.
The model was trained using experimentally resolved protein structures, validated against curated Swiss-Prot annotations, compared with BLAST, and applied to scan the entire Swiss-Prot database.

### Key goals:
- Build a custom pHMM for the Kunitz domain from structural and sequence data
- Validate the model against manually curated Swiss-Prot annotations using 2-fold cross-validation
- Scan the complete Swiss-Prot database and identify Kunitz-domain proteins

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
│   ├── pdb_files/                       # 26 downloaded PDB structures
│   ├── kunitz_chains/                   # Extracted Kunitz domain chains
│   ├── 3TGI.pdb                         # Canonical BPTI reference structure
│   ├── kunitz_with_pdb.fasta            # Kunitz sequences with PDB structures
│   ├── kunitz_all.fasta                 # All collected Kunitz sequences (398)
│   ├── kunitz_sequences.fasta           # Filtered sequences used for alignment
│   ├── kunitz_aligned.fasta             # PDBe-fold structural alignment
│   └── human_kunitz_db.fasta            # 18 human Kunitz sequences (BLAST reference DB)
│
├── models/
│   └── kunitz.hmm                       # Trained HMM profile (102 nodes)
│
├── results/
│   │
│   ├── hmm/
│   │   ├── positives_ids.tsv            # Accessions + lengths of 362 positive sequences
│   │   ├── negatives_ids.tsv            # Accessions + lengths of 362 sampled negatives
│   │   ├── hits_pos.tsv                 # hmmsearch results on positive test set
│   │   ├── hits_neg.tsv                 # hmmsearch results on negative test set
│   │   ├── confusion_matrix.tsv         # TP=357, FP=0, TN=362, FN=5
│   │   ├── performance_metrics.tsv      # Sensitivity, Specificity, Precision, F1, MCC
│   │   ├── evalue_analysis.tsv          # HMM metrics across 9 E-value thresholds
│   │   ├── cross_validation_results.tsv # 2-fold cross-validation results
│   │   └── swissprot_hits.tsv           # 379 Kunitz hits across full Swiss-Prot
│   │
│   ├── blast/
│   │   ├── blast_positives_ids.tsv      # Accessions of 350 positive sequences
│   │   ├── blast_negatives_ids.tsv      # Accessions of 350 negative sequences
│   │   ├── blast_pos.tsv                # BLAST results on positive test set
│   │   ├── blast_neg.tsv                # BLAST results on negative test set
│   │   ├── blast_confusion_matrix.tsv   # TP=349, FP=2, TN=348, FN=1
│   │   ├── blast_performance.tsv        # BLAST metrics at optimal threshold
│   │   └── blast_evalue_analysis.tsv    # BLAST metrics across 9 E-value thresholds
│   │
│   └── figures/
│       ├── confusion_matrices.png       # 2-fold cross-validation confusion matrices
│       ├── roc_curve.png                # ROC curve (AUC = 0.9945, both folds)
│       ├── mcc_thresholds.png           # MCC vs E-value threshold
│       └── hmm_logo.png                 # HMM sequence logo (102 positions)
│
├── scripts/
│   ├── download_pdb.py                  # Download 26 PDB structures
│   ├── extract_chains.py                # Extract Kunitz chains from PDB files
│   ├── validate_hmm.py                  # HMM: hmmsearch + confusion matrix 
│   ├── evalue_analysis.py               # HMM: E-value threshold analysis 
│   ├── cross_validation.py              # HMM: 2-fold cross-validation
│   ├── roc_curve.py                     # HMM: ROC curve (AUC = 0.9945)
│   ├── swissprot.py                     # HMM: Scan full Swiss-Prot 
│   └── blast.py                         # BLAST: prediction + metrics
│
├── Report/
│ └── REPORT.pdf                         # Final MSc project report
│
└── README.md
```

## Workflow

```text
PDB Structures (26)
        │
        ▼
PDBe-fold Structural Alignment
        │
        ▼
Multiple Sequence Alignment
        │
        ▼
HMMER hmmbuild
        │
        ▼
Profile HMM (102 Match States)
        │
        ├── Validation Dataset
        │       ├── Positive Set
        │       └── Negative Set
        │
        ▼
Performance Evaluation
        │
        ▼
Swiss-Prot Full Scan
```

Using :
- hmmbuild [`kunitz.hmm`](models/kunitz.hmm) , [`kunitz_aligned.fasta`](data/kunitz_aligned.fasta) ➔ **Build HMM**
- python [`validate_hmm.py`](scripts/validate_hmm.py) ➔ **Validate Model**
- python [`cross_validation.py`](scripts/cross_validation.py) ➔ **Cross Validation**
- python [`swissprot.py`](scripts/swissprot.py)  ➔ **Swiss-Prot Scan**

## Results

**HMM Performance**

| Metric | Value |
|--------|------|
| Sensitivity | 0.9862 |
| Specificity | 1.0000 |
| Precision | 1.0000 |
| F1 | 0.9930 |
| MCC | 0.9863 |

**Confusion Matrix**

| | Predicted + | Predicted - |
| :--- | :--- | :--- |
| **Actual +** | TP = 357 | FN = 5 |
| **Actual -** | FP = 0 | TN = 362 |

**Cross-Validation**

| Metric | Fold 1 | Fold 2|
|--------|--------|-------|
| **MCC** | 0.9892 | 0.9838 |
| **AUC** | 0.9945 | 0.9945 |

**BLAST Comparison**

| Metric | HMM | BLAST |
|--------|------|------|
| Sensitivity | 0.9862 | **0.9971** |
| Specificity | **1.0000** | 0.9943 |
| Precision | **1.0000** | 0.9943 |
| F1 | 0.9930 | **0.9957** |
| MCC | 0.9863 | **0.9914** |

The HMM achieved perfect specificity with zero false positives, while BLAST achieved slightly higher sensitivity.

## Swiss-Prot Annotation

Scanning the complete Swiss-Prot database identified:
- 379 Kunitz-domain proteins
- 362 previously annotated proteins
- 17 potential unannotated candidates
These results demonstrate the usefulness of structurally informed pHMMs for large-scale protein annotation.

## Technologies

Python 3.12\
HMMER 3.4\
BLAST+\
Biopython\
NumPy\
Matplotlib\
PDBe-fold

## Author

**Kimia Kanouni**

MSc Bioinformatics, University of Bologna

Contact: [kimia.kanouni@studio.unibo.it]

Supervisor: Prof. Emidio Capriotti — http://biofold.org/

## License
MIT License



