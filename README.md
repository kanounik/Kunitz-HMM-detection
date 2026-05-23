# Kunitz Domain Profile HMM

## Table of Contents

- [Overview](#overview)
- [Background](#background)
- [Repository Structure](#repository-structure)
- [Methods](#methods)
-  [Results](#results)
- [Dependencies](#dependencies)
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
│   ├── positives_ids.tsv             # Accessions + lengths of 368 positive sequences
│   ├── negatives_ids.tsv             # Accessions + lengths of 368 sampled negatives
│   ├── hits_pos.tsv                  # hmmsearch results on positive test set (363/368 hits)
│   ├── hits_neg.tsv                  # hmmsearch results on negative test set (0 hits)
│   ├── confusion_matrix.tsv          # TP=363, FP=0, TN=368, FN=5
│   ├── performance_metrics.tsv       # Sensitivity, Specificity, Precision, F1, MCC
│   ├── evalue_analysis.tsv           # Full metrics across 9 E-value thresholds
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
│   ├── validate_hmm.py               # hmmsearch + IDs + confusion matrix + metrics → TSV
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
**1. Structure & Sequence Collection**

26 PDB structures containing Kunitz domains were downloaded using [`download_pdb.py`](scripts/download_pdb.py), starting from the canonical BPTI structure ([`3TGI.pdb`](data/3TGI.pdb), chain I). Kunitz chains were extracted with [`extract_chains.py`](scripts/extract_chains.py). Sequences with a PF00014 annotation were additionally retrieved from UniProt to enrich the training set, stored in [`kunitz_all.fasta`](data/kunitz_all.fasta).

**2. Multiple Sequence Alignment**

The selected Kunitz domain structures were structurally aligned using **PDBe-fold**. The resulting alignment was saved as [`kunitz_aligned.fasta`](data/kunitz_aligned.fasta) and used directly as input for HMMER.

**3. HMM Training**

The profile HMM was built from  [`kunitz_aligned.fasta`](data/kunitz_aligned.fasta) using HMMER hmmbuild, producing  [`kunitz.hmm`](model/kunitz.hmm) with 102 nodes:
```bash
hmmbuild kunitz.hmm kunitz_aligned.fasta
```
A ***consistency test*** confirmed all **26/26** training sequences are recovered by the model.

**4. Validation Dataset**

The validation set was extracted from **Swiss-Prot** (reviewed, manually curated entries only, fragments excluded). Sequence accessions are available in [`positives.tsv`](results/positives.tsv) and [`negatives.tsv`](results/negatives.tsv) for full reproducibility.

- Positive set — 368 proteins with an annotated PF00014 domain in UniProt, To prevent data leakage, 6 human proteins present in the training set were excluded, leaving 362
- Negative set — 368 proteins randomly sampled (random seed = 42) from 564,993 non-Kunitz Swiss-Prot sequences,  balanced to match the positive set.

**5. Performance Evaluation**

[`validate_hmm.py`](scripts/validate_hmm.py) ran hmmsearch on both sets and saved the results to [`hits_pos.tsv`](results/hits_pos.tsv) and [`hits_neg.tsv`](results/hits_neg.tsv). Metrics were computed across 9 E-value thresholds via [`evalue_analysis.py`](scripts/evalue_analysis.py) and saved to [`evalue_analysis.tsv`](results/evalue_analysis.tsv). 2-fold cross-validation was performed with cross_validation.py and ROC curves generated with [`roc_curve.py`](scripts/roc_curve.py).

| Metric | Formula |
| :--- | :--- |
| **Sensitivity (TPR)** | $TPR = \frac{TP}{TP + FN}$ |
| **Specificity (TNR)** | $TNR = \frac{TN}{TN + FP}$ |
| **Precision (PPV)** | $PPV = \frac{TP}{TP + FP}$ |
| **F1 Score** | $F1 = 2 \times \frac{Precision \times Sensitivity}{Precision + Sensitivity}$ |
| **MCC** | $MCC = \frac{(TP \times TN) - (FP \times FN)}{\sqrt{(TP+FP)(TP+FN)(TN+FP)(TN+FN)}}$ |


**6. Swiss-Prot Full Scan**

All 565,361 Swiss-Prot sequences were scanned using [`swissprot.py`](scripts/swissprot.py). Results were saved directly to [`swissprot_hits.tsv`](results/swissprot_hits.tsv) with no intermediate files produced.

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

**Overall Validation (E-value ≤ 1e-3)**

Full confusion matrix → [` confusion_matrix.tsv`](results/confusion_matrix.tsv)

| | Predicted + | Predicted - |
| :--- | :--- | :--- |
| **Actual +** | TP = 357 | FN = 5 |
| **Actual -** | FP = 0 | TN = 362 |

Full metrics → [` performance_metrics.tsv`](results/performance_metrics.tsv)

| Metric | Value |
| :--- | :--- |
| Sensitivity | 0.9862 |
| Specificity | 1.0000 |
| Precision | 1.0000 |
| F1 | 0.9930 |
| MCC | 0.9863 |



**2-Fold Cross-Validation (E-value ≤ 1e-3)**

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
| **AUC** | 0.995 | 0.995 |

<img width="683" height="317" alt="confusion_matrices" src="https://github.com/user-attachments/assets/4c7a29de-a8d3-492b-8e40-5627847f9276" />



<img width="400" height="400" alt="roc_curve" src="https://github.com/user-attachments/assets/53e98e3a-da13-47f4-9e59-2046dc944920" />

---

**False Negatives**

5 Kunitz-domain proteins were not detected at the optimal threshold. All correspond to highly divergent sequences with low similarity to the training set, likely representing evolutionary outliers of the Kunitz fold.

**Swiss-Prot Full Scan**

Full results →  [` swissprot_hits.tsv`](results/swissprot_hits.tsv)

Scanning all **565,361 Swiss-Prot sequences** at E-value ≤ 1e-3:

- **379 proteins** detected as containing a Kunitz domain
- **17 novel hits** beyond the 3622 known annotated proteins — potentially unannotated Kunitz domains
- Distribution spans spiders, snakes, humans, and other metazoa, consistent with the known restriction of this fold to the animal kingdom

**HMM Sequence Logo**

The trained model's position conservation visualized across 102 nodes:

<img width="604" height="205" alt="hmm_logo" src="https://github.com/user-attachments/assets/bef82c69-c780-48d7-8a08-dcd948b51149" />

## Dependencies

**Python ≥ 3.8**

(Also biopython, numpy, matplotlib, seaborn, requests)

**External Tools**
| Tool | Version | Purpose | Install |
| :--- | :--- | :--- | :--- |
| **HMMER** | $\ge$ 3.3 | Build & search HMM profiles | `conda install -c bioconda hmmer` |
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


