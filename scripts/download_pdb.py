import urllib.request
import os

# List of PDB IDs for proteins containing Kunitz domain
# Selected from UniProt Swiss-Prot entries with PF00014 annotation
# These will be used to build the structural alignment for HMM training
pdb_ids = [
    "1BPI",  # BPTI - bovine pancreatic trypsin inhibitor (reference structure)
    "1AAP",  # APP human - amyloid precursor protein
    "1OWT",  # APP mouse
    "1CA0",  # APP rat
    "1ADZ",  # TFPI1 human - tissue factor pathway inhibitor
    "1ZR0",  # TFPI2 human
    "1QNJ",  # AMBP human
    "2WBM",  # Collagen alpha-3(VI) human
    "1RW6",  # APLP2 human
    "1DTX",  # Alpha-dendrotoxin
    "1DX4",  # Dendrotoxin K
    "1DEM",  # Dendrotoxin I
    "2BS1",  # Beta-bungarotoxin
    "1SHL",  # PI-stichotoxin
    "1XMT",  # Theraphotoxin
    "2WNL",  # Carboxypeptidase inhibitor SmCI
    "2CNU",  # Conkunitzin-S1
    "2UUX",  # Textilinin-1
    "2ODY",  # Boophilin
    "4HYN",  # MitTx-alpha
    "1CL1",  # Calcicludine
    "1YEL",  # SPINT1 human
    "3E8N",  # SPINT2 human
    "4BDS",  # Mambaquaretin-1
    "3WMP",  # PPTI
    "1YNT",  # Kunitz inhibitor IX
]

os.makedirs("pdb_files", exist_ok=True)

failed = []
for pdb_id in pdb_ids:
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    out = f"pdb_files/{pdb_id}.pdb"
    if os.path.exists(out):
        print(f"⏭️  Skipping {pdb_id} (already exists)")
        continue
    try:
        urllib.request.urlretrieve(url, out)
        print(f"✅ Downloaded {pdb_id}")
    except Exception as e:
        print(f"❌ Failed {pdb_id}: {e}")
        failed.append(pdb_id)

print(f"\nDownloaded successfully: {len(pdb_ids) - len(failed)}/{len(pdb_ids)}")
if failed:
    print(f"Failed: {failed}")