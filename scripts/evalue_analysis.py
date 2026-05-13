import csv
 
def parse_fasta_ids(path):
    ids = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                parts = line[1:].split("|")
                acc = parts[1] if len(parts) >= 2 else line[1:].split()[0]
                ids.append(acc.strip())
    return ids
 
def parse_tsv(path):
    """Read hits_pos.tsv or hits_neg.tsv — produced by validate_hmm.py"""
    hits = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            acc    = row["target_name"]
            evalue = float(row["e_value"])
            if acc not in hits or evalue < hits[acc]:
                hits[acc] = evalue
    return hits
 
def metrics(tp, fp, tn, fn):
    sens = tp / (tp + fn) if (tp + fn) else 0
    spec = tn / (tn + fp) if (tn + fp) else 0
    prec = tp / (tp + fp) if (tp + fp) else 0
    f1   = 2*prec*sens / (prec+sens) if (prec+sens) else 0
    denom = ((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))**0.5
    mcc  = (tp*tn - fp*fn) / denom if denom else 0
    return sens, spec, prec, f1, mcc
 
# Load
pos_ids  = parse_fasta_ids("test_positives.fasta")
neg_ids  = parse_fasta_ids("test_negatives.fasta")
pos_hits = parse_tsv("hits_pos.tsv")
neg_hits = parse_tsv("hits_neg.tsv")
n = len(pos_ids)
 
# Test multiple thresholds
thresholds = [1e-50, 1e-30, 1e-20, 1e-10, 1e-5, 1e-3, 1e-2, 0.1, 1.0]
 
# Print to terminal (unchanged)
print(f"{'E-value':<12} {'TP':>5} {'FN':>5} {'FP':>5} {'TN':>5} {'Sens':>7} {'Spec':>7} {'Prec':>7} {'F1':>7} {'MCC':>7}")
print("-" * 80)
 
rows = []
for t in thresholds:
    TP = sum(1 for acc in pos_ids if pos_hits.get(acc, 10) <= t)
    FN = n - TP
    FP = sum(1 for acc in neg_ids if neg_hits.get(acc, 10) <= t)
    TN = n - FP
    sens, spec, prec, f1, mcc = metrics(TP, FP, TN, FN)
    print(f"{t:<12.0e} {TP:>5} {FN:>5} {FP:>5} {TN:>5} {sens:>7.4f} {spec:>7.4f} {prec:>7.4f} {f1:>7.4f} {mcc:>7.4f}")
    rows.append({
        "e_value_threshold": t,
        "TP": TP, "FN": FN, "FP": FP, "TN": TN,
        "Sensitivity": round(sens, 4),
        "Specificity": round(spec, 4),
        "Precision":   round(prec, 4),
        "F1":          round(f1,   4),
        "MCC":         round(mcc,  4)
    })
 
# Save to TSV
tsv_path = "evalue_analysis.tsv"
fieldnames = ["e_value_threshold", "TP", "FN", "FP", "TN",
              "Sensitivity", "Specificity", "Precision", "F1", "MCC"]
with open(tsv_path, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
 
print(f"\nSaved: {tsv_path}")