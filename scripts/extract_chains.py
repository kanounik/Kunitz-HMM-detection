from Bio import PDB
from Bio.PDB import PDBParser, PDBIO, Select, PPBuilder
import os

KUNITZ_CHAINS = {
    "1BPI": ("A", None, None,  "BPTI bovine pancreatic trypsin inhibitor"),
    "1AAP": ("A", 1,    56,    "APP human amyloid precursor Kunitz domain"),
    "1OWT": ("A", None, None,  "APP mouse amyloid precursor Kunitz domain"),
    "1CA0": ("I", None, None,  "APP rat amyloid precursor Kunitz domain"),
    "1ADZ": ("A", None, None,  "TFPI1 human Kunitz domain"),      # fixed: take whole chain A
    "1ZR0": ("A", 1,    58,    "TFPI2 human Kunitz domain"),
    "1QNJ": ("A", 1,    64,    "AMBP human bikunin Kunitz domain"),
    "2WBM": ("A", 1,    60,    "Collagen alpha-3(VI) Kunitz domain"),
    "1RW6": ("A", 1,    57,    "APLP2 human Kunitz domain"),
    "1DTX": ("A", None, None,  "Alpha-dendrotoxin"),
    "1DX4": ("A", 1,    60,    "Dendrotoxin K"),
    "1DEM": ("A", None, None,  "Dendrotoxin I"),
    "2BS1": ("B", 1,    60,    "Beta-bungarotoxin B2 chain"),     # fixed: trim to 60
    "1SHL": ("A", 1,    60,    "PI-stichotoxin She2a"),
    "1XMT": ("A", 1,    58,    "KappaPI-theraphotoxin-Hs1a"),
    "2WNL": ("A", 1,    63,    "Carboxypeptidase inhibitor SmCI"),
    "2CNU": ("A", 1,    60,    "Conkunitzin-S1"),
    "2UUX": ("A", None, None,  "Textilinin-1"),
    "2ODY": ("B", 1,    60,    "Boophilin Kunitz domain"),
    "4HYN": ("A", 1,    60,    "MitTx-alpha"),
    "1CL1": ("A", 1,    60,    "Calcicludine"),
    "1YEL": ("A", 1,    58,    "SPINT1 human Kunitz domain 1"),   # fixed: trim to 58
    "3E8N": ("A", 1,    58,    "SPINT2 human Kunitz domain 1"),   # fixed: trim to 58
    "4BDS": ("A", 1,    60,    "Mambaquaretin-1"),
    "3WMP": ("A", 1,    60,    "PPTI Kunitz domain"),             # fixed: trim to 60
    "1YNT": ("A", 1,    63,    "Kunitz-type serine protease inhibitor IX"),
}

class ChainSelect(Select):
    def __init__(self, chain_id):
        self.chain_id = chain_id
    def accept_chain(self, chain):
        return chain.id == self.chain_id

def extract_kunitz_chains(pdb_dir="pdb_files",
                          out_pdb_dir="kunitz_chains",
                          out_fasta="kunitz_sequences.fasta"):

    os.makedirs(out_pdb_dir, exist_ok=True)
    parser  = PDBParser(QUIET=True)
    io      = PDBIO()
    builder = PPBuilder()

    fasta_records = []
    success = []
    failed  = []

    print("Extracting Kunitz domain chains...")
    print("=" * 65)

    for pdb_id, (chain_id, start_res, end_res, description) in KUNITZ_CHAINS.items():
        pdb_file = os.path.join(pdb_dir, f"{pdb_id}.pdb")

        if not os.path.exists(pdb_file):
            print(f"❌ {pdb_id}: file not found")
            failed.append(pdb_id)
            continue

        try:
            structure = parser.get_structure(pdb_id, pdb_file)
            model     = structure[0]
            chain_ids = [c.id for c in model]

            if chain_id not in chain_ids:
                print(f"⚠️  {pdb_id}: chain {chain_id} not found! "
                      f"Available: {chain_ids}")
                failed.append(pdb_id)
                continue

            # Save chain as PDB
            io.set_structure(structure)
            out_pdb = os.path.join(out_pdb_dir, f"{pdb_id}_{chain_id}.pdb")
            io.save(out_pdb, ChainSelect(chain_id))

            # Extract sequence
            chain    = model[chain_id]
            polypeps = builder.build_peptides(chain)
            full_seq = "".join(str(pp.get_sequence()) for pp in polypeps)

            # Trim to Kunitz domain region if specified
            if start_res is not None and end_res is not None:
                sequence = full_seq[start_res-1 : end_res]
            else:
                sequence = full_seq

            if len(sequence) == 0:
                print(f"⚠️  {pdb_id}: empty sequence after trimming")
                failed.append(pdb_id)
                continue

            length_ok = "✅" if len(sequence) <= 80 else "⚠️ "
            fasta_records.append(f">{pdb_id}_{chain_id} {description}")
            fasta_records.append(sequence)

            print(f"✅ {pdb_id} chain {chain_id}: {len(sequence):>3} residues "
                  f"{length_ok} | {description[:40]}")
            success.append(pdb_id)

        except Exception as e:
            print(f"❌ {pdb_id}: ERROR — {e}")
            failed.append(pdb_id)

    # Write FASTA
    with open(out_fasta, "w") as f:
        f.write("\n".join(fasta_records) + "\n")

    print("\n" + "=" * 65)
    print(f"✅ Successfully extracted: {len(success)}/{len(KUNITZ_CHAINS)}")
    print(f"❌ Failed:                 {len(failed)}")
    if failed:
        print(f"   Failed IDs: {failed}")
    print(f"\nOutput:")
    print(f"  → {out_pdb_dir}/          (individual Kunitz chain PDB files)")
    print(f"  → {out_fasta}  (all sequences in FASTA format)")

if __name__ == "__main__":
    extract_kunitz_chains(
        pdb_dir     = "pdb_files",
        out_pdb_dir = "kunitz_chains",
        out_fasta   = "kunitz_sequences.fasta"
    )
