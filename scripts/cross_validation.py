import random
import subprocess
import csv
import os
import tempfile

RANDOM_SEED      = 42
EVALUE_THRESHOLD = 1e-3

# Human Kunitz accessions used in training — excluded to prevent data leakage
# Must match the same set used in validate_hmm.py
TRAINING_HUMAN_ACCESSIONS = {
    "O43278", "O43291", "O95428", "O95925", "P02760",
    "P05067", "P10646", "P12111", "P48307", "Q02388",
    "Q06481", "Q2UY09", "Q8TEU8", "Q96NZ8", "Q8IUA0",
    "Q9BQY6", "P49223", "Q6UDR6"
}

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

def write_fasta(records, path):
    with open(path, "w") as f:
        for acc, seq in records.items():
            f.write(f">{acc}\n{seq}\n")

def run_hmmsearch(hmm, fasta):
    """Run hmmsearch using a temp file. Returns dict of {acc: best_evalue}."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tmp", delete=False) as tmp:
        tmp_path = tmp.name

    subprocess.run(
        ["hmmsearch", "--tblout", tmp_path, "-E", "10", hmm, fasta],
        capture_output=True
    )

    hits = {}
    with open(tmp_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.split()
            if len(cols) < 5:
                continue
            acc    = cols[0]
            evalue = float(cols[4])
            if acc not in hits or evalue < hits[acc]:
                hits[acc] = evalue

    os.remove(tmp_path)
    return hits

def metrics(tp, fp, tn, fn):
    sens  = tp / (tp + fn) if (tp + fn) else 0
    spec  = tn / (tn + fp) if (tn + fp) else 0
    prec  = tp / (tp + fp) if (tp + fp) else 0
    f1    = 2 * prec * sens / (prec + sens) if (prec + sens) else 0
    denom = ((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)) ** 0.5
    mcc   = (tp*tn - fp*fn) / denom if denom else 0
    return sens, spec, prec, f1, mcc

def main():
    random.seed(RANDOM_SEED)

    # Load and filter positives
    print("Loading positives...")
    positives_all = parse_fasta("positives.fasta")
    positives = {acc: seq for acc, seq in positives_all.items()
                 if acc not in TRAINING_HUMAN_ACCESSIONS}
    removed = len(positives_all) - len(positives)
    print(f"  {len(positives_all)} total → {removed} human removed → {len(positives)} for testing")

    # Load negatives
    print("Loading negatives (large file, please wait)...")
    negatives = parse_fasta("negatives.fasta")
    print(f"  {len(negatives)} sequences")

    n = len(positives)
    sampled_neg = dict(random.sample(list(negatives.items()), n))
    print(f"\nSampled {n} negatives to match positives")

    # Split into 2 folds
    pos_items = list(positives.items())
    neg_items = list(sampled_neg.items())
    random.shuffle(pos_items)
    random.shuffle(neg_items)

    half   = n // 2
    folds  = [
        (dict(pos_items[:half]), dict(neg_items[:half])),
        (dict(pos_items[half:]), dict(neg_items[half:])),
    ]

    # Run both folds and collect results
    print(f"\n{'Fold':<6} {'TP':>5} {'FN':>5} {'FP':>5} {'TN':>5} "
          f"{'Sens':>7} {'Spec':>7} {'Prec':>7} {'F1':>7} {'MCC':>7}")
    print("-" * 75)

    all_rows = []
    for i, (pos_fold, neg_fold) in enumerate(folds, 1):

        # Write temp FASTA files
        write_fasta(pos_fold, f"_fold{i}_pos.fasta")
        write_fasta(neg_fold, f"_fold{i}_neg.fasta")

        # Run hmmsearch
        pos_hits = run_hmmsearch("kunitz.hmm", f"_fold{i}_pos.fasta")
        neg_hits = run_hmmsearch("kunitz.hmm", f"_fold{i}_neg.fasta")

        # Compute metrics
        nf = len(pos_fold)
        TP = sum(1 for ev in pos_hits.values() if ev <= EVALUE_THRESHOLD)
        FN = nf - TP
        FP = sum(1 for ev in neg_hits.values() if ev <= EVALUE_THRESHOLD)
        TN = nf - FP

        sens, spec, prec, f1, mcc = metrics(TP, FP, TN, FN)

        print(f"Fold {i}  {TP:>5} {FN:>5} {FP:>5} {TN:>5} "
              f"{sens:>7.4f} {spec:>7.4f} {prec:>7.4f} {f1:>7.4f} {mcc:>7.4f}")

        all_rows.append({
            "Fold":        i,
            "Threshold":   EVALUE_THRESHOLD,
            "TP":          TP,
            "FN":          FN,
            "FP":          FP,
            "TN":          TN,
            "Sensitivity": round(sens, 4),
            "Specificity": round(spec, 4),
            "Precision":   round(prec, 4),
            "F1":          round(f1,   4),
            "MCC":         round(mcc,  4),
        })

        # Clean up temp FASTA files
        os.remove(f"_fold{i}_pos.fasta")
        os.remove(f"_fold{i}_neg.fasta")

    # Save single combined TSV
    output_path = "cross_validation_results.tsv"
    fieldnames  = ["Fold", "Threshold", "TP", "FN", "FP", "TN",
                   "Sensitivity", "Specificity", "Precision", "F1", "MCC"]

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(all_rows)

    print(f"\n  Saved: {output_path}")
    print("\nDone.")

if __name__ == "__main__":
    main()