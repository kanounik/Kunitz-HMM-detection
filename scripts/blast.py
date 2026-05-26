#!/usr/bin/env python3
"""
blast.py
----------------
BLAST-based method to predict the presence of the BPTI/Kunitz domain.

Strategy:
  - Reference DB : human Kunitz proteins extracted from kunitz_all.fasta
  - Positive set : non-human Swiss-Prot Kunitz proteins (positives.fasta)
                   with human training proteins excluded to prevent data leakage
  - Negative set : random non-Kunitz Swiss-Prot proteins (negatives.fasta),
                   balanced to match the positive set
  - Prediction   : a sequence is predicted positive if it has a BLAST hit
                   against the human Kunitz DB below the E-value threshold

Usage:
    python blast.py

Outputs:
    human_kunitz_db.fasta         Human Kunitz reference sequences
    blast_pos.tsv                 BLAST results on positive test set
    blast_neg.tsv                 BLAST results on negative test set
    blast_confusion_matrix.tsv    Confusion matrix
    blast_performance.tsv         Sensitivity, Specificity, Precision, F1, MCC
    blast_evalue_analysis.tsv     Metrics across multiple E-value thresholds
"""

import subprocess
import csv
import os
import random
import tempfile

RANDOM_SEED      = 42
EVALUE_THRESHOLD = 1e-3

# Human Kunitz accessions — used as BLAST reference database
# These are excluded from the test set to prevent data leakage
HUMAN_ACCESSIONS = {
    "O43278", "O43291", "O95428", "O95925", "P02760",
    "P05067", "P10646", "P12111", "P48307", "Q02388",
    "Q06481", "Q2UY09", "Q8TEU8", "Q96NZ8", "Q8IUA0",
    "Q9BQY6", "P49223", "Q6UDR6"
}

# ─── Fasta helpers ────────────────────────────────────────────────────────────

def parse_fasta(path):
    records = {}
    header = None
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    records[header] = "".join(seq)
                parts = line[1:].split("|")
                header = parts[1] if len(parts) >= 2 else line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if header:
        records[header] = "".join(seq)
    return records

def parse_fasta_full(path):
    """Parse FASTA keeping full header line for makeblastdb."""
    records = {}
    header = None
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    records[header] = "".join(seq)
                header = line[1:]
                seq = []
            else:
                seq.append(line)
    if header:
        records[header] = "".join(seq)
    return records

def write_fasta(records, path):
    with open(path, "w") as f:
        for acc, seq in records.items():
            f.write(f">{acc}\n{seq}\n")

def save_ids_tsv(records, label, path):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["accession", "length", "label"])
        for acc, seq in records.items():
            writer.writerow([acc, len(seq), label])
    print(f"  Saved: {path}  ({len(records)} sequences)")

# ─── BLAST helpers ────────────────────────────────────────────────────────────

def make_blast_db(fasta_path, db_path):
    """Build a BLAST protein database."""
    subprocess.run(
        ["makeblastdb", "-in", fasta_path, "-dbtype", "prot",
         "-out", db_path, "-title", "human_kunitz"],
        capture_output=True
    )
    print(f"  BLAST DB built: {db_path}")

def run_blastp(query_fasta, db_path, tsv_path, evalue=10):
    """
    Run blastp and save results as TSV.
    Output format: qseqid sseqid evalue bitscore pident length
    """
    subprocess.run(
        ["blastp",
         "-query",   query_fasta,
         "-db",      db_path,
         "-out",     tsv_path,
         "-outfmt",  "6 qseqid sseqid evalue bitscore pident length",
         "-evalue",  str(evalue),
         "-num_threads", "4"],
        capture_output=True
    )

    # Parse best hit per query
    hits = {}
    rows = []
    with open(tsv_path) as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 6:
                continue
            qseqid  = cols[0]
            sseqid  = cols[1]
            evalue  = float(cols[2])
            bitscore = float(cols[3])
            pident  = float(cols[4])
            length  = int(cols[5])
            if qseqid not in hits or evalue < hits[qseqid]:
                hits[qseqid] = evalue
            rows.append({
                "query":    qseqid,
                "subject":  sseqid,
                "e_value":  evalue,
                "bitscore": bitscore,
                "pident":   pident,
                "length":   length
            })

    # Overwrite with clean TSV
    with open(tsv_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["query", "subject", "e_value", "bitscore", "pident", "length"],
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"  Saved: {tsv_path}  ({len(hits)} queries with hits)")
    return hits

# ─── Metrics ──────────────────────────────────────────────────────────────────

def metrics(tp, fp, tn, fn):
    sens  = tp / (tp + fn) if (tp + fn) else 0
    spec  = tn / (tn + fp) if (tn + fp) else 0
    prec  = tp / (tp + fp) if (tp + fp) else 0
    f1    = 2 * prec * sens / (prec + sens) if (prec + sens) else 0
    denom = ((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)) ** 0.5
    mcc   = (tp*tn - fp*fn) / denom if denom else 0
    return sens, spec, prec, f1, mcc

def evaluate(pos_ids, neg_ids, pos_hits, neg_hits, threshold):
    TP = sum(1 for acc in pos_ids if pos_hits.get(acc, 999) <= threshold)
    FN = len(pos_ids) - TP
    FP = sum(1 for acc in neg_ids if neg_hits.get(acc, 999) <= threshold)
    TN = len(neg_ids) - FP
    return TP, FN, FP, TN

def save_confusion_matrix(tp, fp, tn, fn, path):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["", "Predicted_Positive", "Predicted_Negative"])
        writer.writerow(["Actual_Positive", tp, fn])
        writer.writerow(["Actual_Negative", fp, tn])
    print(f"  Saved: {path}")

def save_performance(sens, spec, prec, f1, mcc, threshold, path):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Threshold", "Sensitivity", "Specificity",
                         "Precision", "F1", "MCC"])
        writer.writerow([threshold, f"{sens:.4f}", f"{spec:.4f}",
                         f"{prec:.4f}", f"{f1:.4f}", f"{mcc:.4f}"])
    print(f"  Saved: {path}")

def save_evalue_analysis(pos_ids, neg_ids, pos_hits, neg_hits, path):
    thresholds = [1e-50, 1e-30, 1e-20, 1e-10, 1e-5, 1e-3, 1e-2, 0.1, 1.0]
    fieldnames = ["e_value_threshold", "TP", "FN", "FP", "TN",
                  "Sensitivity", "Specificity", "Precision", "F1", "MCC"]

    print(f"\n{'E-value':<12} {'TP':>5} {'FN':>5} {'FP':>5} {'TN':>5} "
          f"{'Sens':>7} {'Spec':>7} {'Prec':>7} {'F1':>7} {'MCC':>7}")
    print("-" * 80)

    rows = []
    for t in thresholds:
        TP, FN, FP, TN = evaluate(pos_ids, neg_ids, pos_hits, neg_hits, t)
        sens, spec, prec, f1, mcc = metrics(TP, FP, TN, FN)
        print(f"{t:<12.0e} {TP:>5} {FN:>5} {FP:>5} {TN:>5} "
              f"{sens:>7.4f} {spec:>7.4f} {prec:>7.4f} {f1:>7.4f} {mcc:>7.4f}")
        rows.append({
            "e_value_threshold": t,
            "TP": TP, "FN": FN, "FP": FP, "TN": TN,
            "Sensitivity": round(sens, 4),
            "Specificity": round(spec, 4),
            "Precision":   round(prec, 4),
            "F1":          round(f1,   4),
            "MCC":         round(mcc,  4)
        })

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"\n  Saved: {path}")
    return rows

# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    random.seed(RANDOM_SEED)

    # 1. Build human Kunitz reference database
    print("=== Step 1: Building human Kunitz BLAST database ===")
    all_kunitz = parse_fasta_full("kunitz_all.fasta")

    human_records = {
        h: s for h, s in all_kunitz.items()
        if any(acc in h for acc in HUMAN_ACCESSIONS)
    }
    print(f"  Human Kunitz sequences found: {len(human_records)}")
    write_fasta(human_records, "human_kunitz_db.fasta")
    make_blast_db("human_kunitz_db.fasta", "human_kunitz_db")

    # 2. Load and filter positive set (exclude human proteins)
    print("\n=== Step 2: Preparing test sets ===")
    positives_all = parse_fasta("positives.fasta")
    positives = {acc: seq for acc, seq in positives_all.items()
                 if acc not in HUMAN_ACCESSIONS}
    print(f"  Positives: {len(positives_all)} total, "
          f"{len(positives_all)-len(positives)} human removed, "
          f"{len(positives)} for testing")

    # 3. Load and sample negative set
    print("  Loading negatives (large file, please wait)...")
    negatives_all = parse_fasta("negatives.fasta")
    n = len(positives)
    sampled_neg = dict(random.sample(list(negatives_all.items()), n))
    print(f"  Negatives: {n} sampled to match positives")

    # Save test FASTA files
    write_fasta(positives,   "blast_test_positives.fasta")
    write_fasta(sampled_neg, "blast_test_negatives.fasta")

    # Save ID lists
    save_ids_tsv(positives,   "positive", "blast_positives_ids.tsv")
    save_ids_tsv(sampled_neg, "negative", "blast_negatives_ids.tsv")

    # 4. Run BLAST
    print("\n=== Step 3: Running BLAST ===")
    print("  Searching positives against human Kunitz DB...")
    pos_hits = run_blastp("blast_test_positives.fasta",
                          "human_kunitz_db", "blast_pos.tsv")

    print("  Searching negatives against human Kunitz DB...")
    neg_hits = run_blastp("blast_test_negatives.fasta",
                          "human_kunitz_db", "blast_neg.tsv")

    # 5. Evaluate at optimal threshold
    print(f"\n=== Step 4: Evaluation at E-value ≤ {EVALUE_THRESHOLD} ===")
    pos_ids = list(positives.keys())
    neg_ids = list(sampled_neg.keys())

    TP, FN, FP, TN = evaluate(pos_ids, neg_ids, pos_hits, neg_hits, EVALUE_THRESHOLD)
    sens, spec, prec, f1, mcc = metrics(TP, FP, TN, FN)

    print(f"  TP={TP}  FN={FN}  FP={FP}  TN={TN}")
    print(f"\n  Sensitivity : {sens:.4f}")
    print(f"  Specificity : {spec:.4f}")
    print(f"  Precision   : {prec:.4f}")
    print(f"  F1 Score    : {f1:.4f}")
    print(f"  MCC         : {mcc:.4f}")

    save_confusion_matrix(TP, FP, TN, FN, "blast_confusion_matrix.tsv")
    save_performance(sens, spec, prec, f1, mcc,
                     EVALUE_THRESHOLD, "blast_performance.tsv")

    # 6. E-value threshold analysis
    print("\n=== Step 5: E-value threshold analysis ===")
    save_evalue_analysis(pos_ids, neg_ids, pos_hits, neg_hits,
                         "blast_evalue_analysis.tsv")

    # 7. False negatives
    print("\n=== False Negatives (Kunitz missed by BLAST) ===")
    count = 0
    for acc in pos_ids:
        if pos_hits.get(acc, 999) > EVALUE_THRESHOLD:
            print(f"  {acc}  evalue={pos_hits.get(acc, 'not found')}")
            count += 1
            if count >= 10:
                print("  ...")
                break

    print("\nDone. All results saved.")
    print("\n=== Summary ===")
    print(f"  Reference DB  : {len(human_records)} human Kunitz sequences")
    print(f"  Positive set  : {len(pos_ids)} non-human Kunitz sequences")
    print(f"  Negative set  : {len(neg_ids)} random non-Kunitz sequences")
    print(f"  Threshold     : E-value ≤ {EVALUE_THRESHOLD}")
    print(f"  TP={TP}  FN={FN}  FP={FP}  TN={TN}")
    print(f"  MCC={mcc:.4f}  F1={f1:.4f}")

if __name__ == "__main__":
    main()