#!/usr/bin/env python3
"""
swissprot.py
-----------------
Scans all Swiss-Prot sequences against the Kunitz HMM profile.
Saves results directly as TSV.

Usage:
    python swissprot.py
"""

import subprocess
import csv
import os

HMM        = "kunitz.hmm"
FASTA      = "swissprot_all.fasta"
OUTPUT_TSV = "swissprot_hits.tsv"
TMP_TBL    = "_tmp_swissprot.tbl"
THRESHOLD  = 1e-3

def main():
    print("=== Scanning Swiss-Prot with Kunitz HMM ===")
    print(f"  HMM    : {HMM}")
    print(f"  Input  : {FASTA}")
    print(f"  Output : {OUTPUT_TSV}")
    print()

    # Run hmmsearch
    print("Running hmmsearch (this may take a few minutes)...")
    subprocess.run(
        ["hmmsearch", "--tblout", TMP_TBL, "-E", "1000",
         "--cpu", "4", HMM, FASTA],
        capture_output=True
    )
    print("  Done.")

    # Parse and save as TSV
    rows = []
    with open(TMP_TBL) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.split()
            if len(cols) < 6:
                continue
            rows.append({
                "target_name": cols[0],
                "accession":   cols[1],
                "e_value":     float(cols[4]),
                "score":       float(cols[5]),
                "description": " ".join(cols[18:]) if len(cols) > 18 else ""
            })

    # Save ALL hits (any e-value) to TSV
    with open(OUTPUT_TSV, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["target_name", "accession", "e_value", "score", "description"],
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(rows)

    # Delete temporary tbl immediately
    os.remove(TMP_TBL)

    # Summary
    hits_at_threshold = [r for r in rows if r["e_value"] <= THRESHOLD]
    print(f"\n=== Results ===")
    print(f"  Total sequences scanned : see {FASTA}")
    print(f"  Total hits (any e-value): {len(rows)}")
    print(f"  Hits at E-value <= {THRESHOLD} : {len(hits_at_threshold)}")
    print(f"\n  Saved: {OUTPUT_TSV}")

    print(f"\n=== Top 10 hits (E-value <= {THRESHOLD}) ===")
    print(f"{'target_name':<20} {'e_value':>12} {'score':>8}  description")
    print("-" * 70)
    for r in sorted(hits_at_threshold, key=lambda x: x["e_value"])[:10]:
        print(f"{r['target_name']:<20} {r['e_value']:>12.2e} {r['score']:>8.1f}  {r['description'][:30]}")

    print("\nDone.")

if __name__ == "__main__":
    main()