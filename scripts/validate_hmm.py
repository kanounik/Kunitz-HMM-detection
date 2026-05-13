import random
import subprocess
 
RANDOM_SEED = 42
EVALUE_THRESHOLD = 1e-3
 
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
 
def run_hmmsearch(hmm, fasta, tbl):
    subprocess.run(["hmmsearch", "--tblout", tbl, "-E", "10", hmm, fasta],
                   capture_output=True)
    hits = {}
    with open(tbl) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.split()
            if len(cols) < 5:
                continue
            acc, evalue = cols[0], float(cols[4])
            if acc not in hits or evalue < hits[acc]:
                hits[acc] = evalue
    return hits
 
def metrics(tp, fp, tn, fn):
    sens = tp / (tp + fn) if (tp + fn) else 0
    spec = tn / (tn + fp) if (tn + fp) else 0
    prec = tp / (tp + fp) if (tp + fp) else 0
    f1   = 2 * prec * sens / (prec + sens) if (prec + sens) else 0
    denom = ((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)) ** 0.5
    mcc  = (tp*tn - fp*fn) / denom if denom else 0
    return sens, spec, prec, f1, mcc
 
def main():
    random.seed(RANDOM_SEED)
 
    print("Loading positives...")
    positives = parse_fasta("positives.fasta")
    print(f"  {len(positives)} sequences")
 
    print("Loading negatives (large file, please wait)...")
    negatives = parse_fasta("negatives.fasta")
    print(f"  {len(negatives)} sequences")
 
    n = len(positives)
    sampled_neg = dict(random.sample(list(negatives.items()), n))
    print(f"\nSampled {n} negatives to match positives")
 
    write_fasta(positives,   "test_positives.fasta")
    write_fasta(sampled_neg, "test_negatives.fasta")
    print("Written: test_positives.fasta, test_negatives.fasta")
 
    print("\nRunning hmmsearch on positives...")
    pos_hits = run_hmmsearch("kunitz.hmm", "test_positives.fasta", "hits_pos.tbl")
    print("Running hmmsearch on negatives...")
    neg_hits = run_hmmsearch("kunitz.hmm", "test_negatives.fasta", "hits_neg.tbl")
 
    TP = sum(1 for ev in pos_hits.values() if ev <= EVALUE_THRESHOLD)
    FN = n - TP
    FP = sum(1 for ev in neg_hits.values() if ev <= EVALUE_THRESHOLD)
    TN = n - FP
 
    print(f"\nE-value threshold: {EVALUE_THRESHOLD}")
    print(f"  TP={TP}  FN={FN}  TN={TN}  FP={FP}")
 
    sens, spec, prec, f1, mcc = metrics(TP, FP, TN, FN)
    print(f"\n  Sensitivity : {sens:.4f}")
    print(f"  Specificity : {spec:.4f}")
    print(f"  Precision   : {prec:.4f}")
    print(f"  F1 Score    : {f1:.4f}")
    print(f"  MCC         : {mcc:.4f}")
 
    print("\nFalse Negatives (Kunitz missed):")
    count = 0
    for acc in positives:
        if pos_hits.get(acc, 10) > EVALUE_THRESHOLD:
            print(f"  {acc}  evalue={pos_hits.get(acc, 'not found')}")
            count += 1
            if count >= 5:
                print("  ...")
                break
 
    print("\nFalse Positives (non-Kunitz detected):")
    count = 0
    for acc, ev in neg_hits.items():
        if ev <= EVALUE_THRESHOLD:
            print(f"  {acc}  evalue={ev:.2e}")
            count += 1
            if count >= 5:
                print("  ...")
                break
 
if __name__ == "__main__":
    main()