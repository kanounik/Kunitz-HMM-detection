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
This project builds and validates a Profile Hidden Markov Model (pHMM) for the detection of the Kunitz-type protease inhibitor domain (PFAM: PF00014) in protein sequences. The trained model is then applied to annotate Kunitz domains across the entire Swiss-Prot database.

Key goals:
1- Build a custom pHMM for the Kunitz domain from structural data
2- Validate the model against manually curated Swiss-Prot annotations
3- Annotate Kunitz domains across Swiss-Prot and analyze their distribution


## Background
The Kunitz domain is a short (~58 residues), disulfide-rich alpha+beta fold found in serine protease inhibitors. The canonical example is BPTI (Bovine Pancreatic Trypsin Inhibitor, PDB: 3TGI, chain I), which tightly binds trypsin via its Lys15 residue. The fold is stabilized by three conserved disulfide bonds (Cys5–Cys55, Cys14–Cys38, Cys30–Cys51).

Kunitz domains appear in proteins such as:
* Aprotinin (BPTI) — antifibrinolytic drug (Trasylol)
* Alzheimer's amyloid precursor protein (APP)
* Tissue factor pathway inhibitor (TFPI)

