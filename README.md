# Kunitz Domain Profile HMM

## Table of Contents

- [Overview](#overview)
- [Background](#background)
- [Repository Structure](#repository-structure)
- [Methods](#methods)
-  [Results](#results)- [Dependencies](#dependencies)
- [Data Sources](#data-sources)
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
├── README.md
├── requirements.txt
└── environment.yml
```

## Methods
**1. Structure & Sequence Collection**

26 PDB structures containing Kunitz domains were downloaded using [`download_pdb.py`](scripts/download_pdb.py), starting from the canonical BPTI structure ([`3TGI.pdb`](data/3TGI.pdb), chain I). Kunitz chains were extracted with [`extract_chains.py`](scripts/extract_chains.py). Sequences with a PF00014 annotation were additionally retrieved from UniProt to enrich the training set, stored in [`kunitz_all.fasta`](data/kunitz_all.fasta) (398 sequences, including 18 human proteins).

**2. Structural Alignment**

The selected Kunitz domain structures were structurally aligned using **PDBe-fold**. The resulting alignment was saved as [`kunitz_aligned.fasta`](data/kunitz_aligned.fasta) and used directly as input for HMMER.

**3. HMM Training**

The profile HMM was built from  [`kunitz_aligned.fasta`](data/kunitz_aligned.fasta) using HMMER hmmbuild, producing  [`kunitz.hmm`](models/kunitz.hmm) with 102 nodes:
```bash
hmmbuild kunitz.hmm kunitz_aligned.fasta
```
A ***consistency test*** confirmed all **26/26** training sequences are recovered by the model.

**4.  BLAST Reference Database**

18 human Kunitz proteins from [`kunitz_all.fasta`](data/kunitz_all.fasta) were extracted to build [`human_kunitz_db.fasta`](data/human_kunitz_db.fasta), used as the BLAST reference database. A sequence is predicted positive if it produces a hit below the E-value threshold.

**5. Validation Dataset**

The validation set was extracted from **Swiss-Prot** (reviewed, manually curated, no fragments). To prevent data leakage, human proteins present in the training set were excluded from the test set.

| Dataset Component | HMM | BLAST |
|-------------------|-----|-------|
| Positive set | 362 non-human Kunitz proteins | 350 non-human Kunitz proteins |
| Negative set | 362 random non-Kunitz proteins | 350 random non-Kunitz proteins |
| Human proteins excluded | 6 proteins | 18 proteins |

Sequence accessions available in:

[`positives_ids.tsv`](results/hmm/positives_ids.tsv), [`negatives_ids.tsv`](results/hmm/negatives_ids.tsv), [`blast_positives_ids.tsv`](results/blast/blast_positives_ids.tsv), [`blast_negatives_ids.tsv`](results/blast/blast_negatives_ids.tsv).



**6. Performance Evaluation**

Both methods were evaluated across 9 E-value thresholds using:

| Metric | Formula |
| :--- | :--- |
| **Sensitivity (TPR)** | $TPR = \frac{TP}{TP + FN}$ |
| **Specificity (TNR)** | $TNR = \frac{TN}{TN + FP}$ |
| **Precision (PPV)** | $PPV = \frac{TP}{TP + FP}$ |
| **F1 Score** | $F1 = 2 \times \frac{Precision \times Sensitivity}{Precision + Sensitivity}$ |
| **MCC** | $MCC = \frac{(TP \times TN) - (FP \times FN)}{\sqrt{(TP+FP)(TP+FN)(TN+FP)(TN+FN)}}$ |


**7. Swiss-Prot Full Scan**

All 565,361 Swiss-Prot sequences were scanned using the pHMM via  [`swissprot.py`](scripts/swissprot.py), saving results to  [`swissprot_hits.tsv`](results/swissprot_hits.tsv).

## Results
**Consistency Test**

All **26/26** training sequences were recovered by the model — HMM built correctly ✅

**E-value Threshold Analysis**

Full results in [` evalue_analysis.tsv`](results/evalue_analysis.tsv).

| E-value | TP | FN | FP | TN | Sens | Spec | Prec | F1 | MCC |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1e-50 | 3 | 359 | 0 | 362 | 0.0083 | 1.0000 | 1.0000 | 0.0164 | 0.0645 |
| 1e-30 | 18 | 344 | 0 | 362 | 0.0497 | 1.0000 | 1.0000 | 0.0947 | 0.1602 |
| 1e-20 | 36 | 327 | 0 | 362 | 0.0967 | 1.0000 | 1.0000 | 0.1763 | 0.2251 |
| 1e-10 | 276 | 86 | 0 | 362 | 0.7624 | 1.0000 | 1.0000 | 0.8651 | 0.7849 |
| 1e-05 | 352 | 10 | 0 | 362 | 0.9724 | 1.0000 | 1.0000 | 0.9860 | 0.9727 |
| **1e-03** | **357** | **5** | **0** | **362** | **0.9862** | **1.0000** | **1.0000** | **0.9930** | **0.9863** |
| 1e-02 | 357 | 5 | 0 | 362 | 0.9862 | 1.0000 | 1.0000 | 0.9930 | 0.9863 |
| 0.1 | 357 | 5 | 0 | 362 | 0.9862 | 1.0000 | 1.0000 | 0.9930 | 0.9863 |
| 1.0 | 358 | 4 | 0 | 362 | 0.9890 | 1.0000 | 1.0000 | 0.9944 | 0.9891 |

MCC reaches a stable plateau from **1e-5** onwards. The optimal threshold is **E-value ≤ 1e-3**.

<img width="600" height="300" alt="mcc_thresholds" src="https://github.com/user-attachments/assets/29e4cd5f-0414-4ff8-8968-c165686163ef" />

**HMM — Overall Validation (E-value ≤ 1e-3)**

Full confusion matrix → [` confusion_matrix.tsv`](results/confusion_matrix.tsv)
Full metrics → [` performance_metrics.tsv`](results/performance_metrics.tsv)

| | Predicted + | Predicted - |
| :--- | :--- | :--- |
| **Actual +** | TP = 357 | FN = 5 |
| **Actual -** | FP = 0 | TN = 362 |


**HMM —2-Fold Cross-Validation (E-value ≤ 1e-3)**

| | Fold 1 | Fold 2 |
| :--- | :--- | :--- |
| **TP** | 182 | 181 |
| **FN** | 2 | 3 |
| **FP** | 0 | 0 |
| **TN** | 184 | 184 |
| **Sensitivity** | 0.9891 | 0.9837 |
| **Specificity** | 1.0000 | 1.0000 |
| **Precision** | 1.0000 | 1.0000 |
| **MCC** | 0.9892 | 0.9838 |
| **AUC** | 0.9945 | 0.9945 |

<img width="683" height="317" alt="confusion_matrices" src="https://github.com/user-attachments/assets/4c7a29de-a8d3-492b-8e40-5627847f9276" />



<img width="400" height="400" alt="roc_curve" src="https://github.com/user-attachments/assets/53e98e3a-da13-47f4-9e59-2046dc944920" />

**BLAST — E-value Threshold Analysis**

Full results in [`blast_evalue_analysis.tsv`](results/blast/blast_evalue_analysis.tsv).

| E-value | TP | FN | FP | TN | Sens | Spec | Prec | F1 | MCC |
|--------|----|----|----|----|------|------|------|------|------|
| 1e-50 | 39 | 311 | 0 | 350 | 0.1114 | 1.0000 | 1.0000 | 0.2005 | 0.2429 |
| 1e-30 | 49 | 301 | 0 | 350 | 0.1400 | 1.0000 | 1.0000 | 0.2456 | 0.2744 |
| 1e-20 | 62 | 288 | 0 | 350 | 0.1771 | 1.0000 | 1.0000 | 0.3010 | 0.3117 |
| 1e-10 | 339 | 11 | 1 | 349 | 0.9686 | 0.9971 | 0.9971 | 0.9826 | 0.9661 |
| 1e-05 | 346 | 4 | 1 | 349 | 0.9886 | 0.9971 | 0.9971 | 0.9928 | 0.9858 |
| **1e-03** | **349** | **1** | **2** | **348** | **0.9971** | **0.9943** | **0.9943** | **0.9957** | **0.9914** |
| 1e-02 | 349 | 1 | 3 | 347 | 0.9971 | 0.9914 | 0.9915 | 0.9943 | 0.9886 |
| 0.1 | 350 | 0 | 30 | 320 | 1.0000 | 0.9143 | 0.9211 | 0.9589 | 0.9177 |
| 1.0 | 350 | 0 | 155 | 195 | 1.0000 | 0.5571 | 0.6931 | 0.8187 | 0.6214 |

MCC reaches a stable plateau from 1e-5 onwards. Optimal threshold: **E-value ≤ 1e-3.**

**BLAST — Overall Validation (E-value ≤ 1e-3)**

Full confusion matrix → [`blast_confusion_matrix.tsv`](results/blast/blast_confusion_matrix.tsv).
Full metrics → [`blast_performance.tsv`](results/blast/blast_performance.tsv).

|              | Predicted + | Predicted − |
|--------------|-------------|-------------|
| **Actual +** | TP = 349    | FN = 1      |
| **Actual −** | FP = 2      | TN = 348    |

**HMM vs BLAST — Comparison**

| Metric | HMM | BLAST |
|--------|------|------|
| Sensitivity | 0.9862 | **0.9971** |
| Specificity | **1.0000** | 0.9943 |
| Precision | **1.0000** | 0.9943 |
| F1 | 0.9930 | **0.9957** |
| MCC | 0.9863 | **0.9914** |
| False Positives | **0** | 2 |
| False Negatives | 5 | **1** |


**Key observations:**

BLAST achieves **higher sensitivity** (0.9971 vs 0.9862) — misses only 1 Kunitz protein
HMM achieves **perfect specificity** (1.0000 vs 0.9943) — zero false positives
Both methods perform excellently with MCC > 0.98
HMM is more suitable when precision is critical; BLAST when recall is critical



**False Negatives**

**HMM** — 5 proteins not detected: highly divergent sequences with low similarity to the training set, likely evolutionary outliers of the Kunitz fold.

**BLAST** — 1 protein not detected: insufficient sequence similarity to any human Kunitz reference.


**Swiss-Prot Full Scan (HMM)**

Full results →  [` swissprot_hits.tsv`](results/hmm/swissprot_hits.tsv)

Scanning all **565,361 Swiss-Prot sequences** at E-value ≤ 1e-3:

- **379 proteins** detected as containing a Kunitz domain
- **17 novel hits** beyond the 362 known annotated proteins — potentially unannotated Kunitz domains
- Distribution spans spiders, snakes, humans, and other metazoa, consistent with the known restriction of this fold to the animal kingdom

**HMM Sequence Logo**

The trained model's position conservation visualized across 102 nodes:

<img width="604" height="205" alt="hmm_logo" src="https://github.com/user-attachments/assets/bef82c69-c780-48d7-8a08-dcd948b51149" />

## Dependencies

**Python ≥ 3.8**

(Also biopython, numpy, matplotlib, requests)

**External Tools**
| Tool | Version | Purpose | Install |
| :--- | :--- | :--- | :--- |
| **HMMER** | $\ge$ 3.3 | Build & search HMM profiles | `conda install -c bioconda hmmer` |
| **BLAST+** | $\ge$ 2.12 | BLAST-based prediction | `conda install -c bioconda blast` |
| **PDBe-fold** | web | Structural alignment | [ebi.ac.uk/msd-srv/ssm/](https://www.ebi.ac.uk/msd-srv/ssm/) |

## Data Sources

| Resource | URL | Usage |
| :--- | :--- | :--- |
| **Protein Data Bank** | [rcsb.org](https://www.rcsb.org) | 26 Kunitz domain structures |
| **UniProt / Swiss-Prot** | [uniprot.org](https://www.uniprot.org) | Positive set (PF00014) & full scan |
| **PFAM via InterPro** | [PF00014](https://www.ebi.ac.uk/interpro/entry/pfam/PF00014/) | Domain annotation reference |

## Authors

**Kimia Kanouni**

MSc Bioinformatics, University of Bologna

Contact: [kimia.kanouni@studio.unibo.it]

Prof. Emidio Capriotti — http://biofold.org/

## License
This project is licensed under the MIT License — see [`LICENSE`](LICENSE) for details.
Data from PDB, UniProt, and PFAM are subject to their respective terms of use.

