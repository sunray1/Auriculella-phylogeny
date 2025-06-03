#!/usr/bin/env python3
import sys
import json
import os
import csv
from collections import defaultdict

def sanitize(s):
    return s.replace(" ", "_").replace("/", "_") if s else "NA"

def main(json_path, out_dir):
    os.makedirs(out_dir, exist_ok=True)

    # Load JSON records
    records = []
    with open(json_path) as f:
        for line in f:
            line = line.strip()
            if line:
                records.append(json.loads(line))

    # Organize sequences by locus
    sequences_by_locus = defaultdict(list)
    summary_counts = defaultdict(lambda: defaultdict(int))

    for rec in records:
        seq = rec.get("nuc", "").replace("-", "").upper()
        if not seq:
            continue
        
        # Get genus and species
        genus = rec.get("genus", "")
        species = rec.get("species")
        
        # If species is missing, substitute with genus_only label
        if not species or species.lower() in ["na", "null"]:
            species_raw = f"{genus}_sp"
        else:
            species_raw = species

        species = sanitize(species_raw)
        locus = sanitize(rec.get("marker_code"))
        
        processid_raw = rec.get("processid", "")
        if processid_raw.startswith("GB"):
            continue  # Skip GenBank records
        
        processid = sanitize(processid_raw)

        sampleid = sanitize(rec.get("sampleid"))

        header = f">{processid}_{sampleid}_{species}"
        sequences_by_locus[locus].append((header, seq))
        summary_counts[locus][species] += 1

    # Write one FASTA file per locus
    for locus, seqs in sequences_by_locus.items():
        fasta_path = os.path.join(out_dir, f"{locus}.fas")
        with open(fasta_path, "w") as fout:
            for header, seq in seqs:
                fout.write(f"{header}\n{seq}\n")

    # Write summary matrix
    all_loci = sorted(summary_counts.keys())
    all_species = sorted({species for locus in summary_counts for species in summary_counts[locus]})

    matrix_path = os.path.join(out_dir, "summary_by_locus.csv")
    with open(matrix_path, "w", newline="") as fout:
        writer = csv.writer(fout)
        writer.writerow(["species"] + all_loci)
        for species in all_species:
            row = [species]
            for locus in all_loci:
                row.append(summary_counts[locus].get(species, 0))
            writer.writerow(row)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: ./json_to_fasta.py input.json output_dir")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
