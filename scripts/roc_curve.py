import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
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
    """Read hits TSV — produced by validate_hmm.py"""
    hits = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            acc    = row["target_name"]
            evalue = float(row["e_value"])
            if acc not in hits or evalue < hits[acc]:
                hits[acc] = evalue
    return hits

# Load data
pos_ids  = parse_fasta_ids("test_positives.fasta")
neg_ids  = parse_fasta_ids("test_negatives.fasta")
pos_hits = parse_tsv("hits_pos.tsv")
neg_hits = parse_tsv("hits_neg.tsv")
n = len(pos_ids)

# Build scores: use -log10(evalue) as score
# Higher score = more likely Kunitz
def get_score(acc, hits):
    ev = hits.get(acc, 10.0)  # default high evalue if not found
    ev = max(ev, 1e-100)      # avoid log(0)
    return -np.log10(ev)

pos_scores = [get_score(acc, pos_hits) for acc in pos_ids]
neg_scores = [get_score(acc, neg_hits) for acc in neg_ids]

# Build ROC curve by varying threshold
all_scores = np.linspace(
    min(pos_scores + neg_scores) - 0.1,
    max(pos_scores + neg_scores) + 0.1,
    500
)

tprs = []
fprs = []

for threshold in all_scores:
    TP = sum(1 for s in pos_scores if s >= threshold)
    FN = n - TP
    FP = sum(1 for s in neg_scores if s >= threshold)
    TN = n - FP
    tpr = TP / (TP + FN) if (TP + FN) else 0
    fpr = FP / (FP + TN) if (FP + TN) else 0
    tprs.append(tpr)
    fprs.append(fpr)

# Compute AUC using trapezoidal rule
auc = abs(np.trapz(tprs, fprs))

# Plot
fig, ax = plt.subplots(figsize=(6, 6))
ax.plot(fprs, tprs, color='steelblue', lw=2, label=f'HMM (AUC = {auc:.3f})')
ax.plot([0, 1], [0, 1], color='gray', lw=1, linestyle='--', label='Random')
ax.set_xlabel('False Positive Rate', fontsize=12)
ax.set_ylabel('True Positive Rate', fontsize=12)
ax.set_title('ROC Curve - Kunitz HMM Classifier', fontsize=13)
ax.legend(loc='lower right', fontsize=11)
ax.set_xlim([0, 1])
ax.set_ylim([0, 1.02])
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('roc_curve.png', dpi=150)
print(f"ROC curve saved: roc_curve.png")
print(f"AUC = {auc:.4f}")