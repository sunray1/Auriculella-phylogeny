#!/usr/bin/env python3
import argparse
import re

MODEL_MAP = {
    "GTR": "GTRCAT",
    "GTR+G": "GTRGAMMA",
    "GTR+I+G": "GTRGAMMAI",
}

def parse_args():
    parser = argparse.ArgumentParser(description="Pick RAxML model and extract partition file from best_scheme.txt")
    parser.add_argument("-i", "--input", required=True, help="Path to best_scheme.txt")
    parser.add_argument("-m", "--model-out", required=True, help="Path to write RAxML model (e.g. GTRGAMMA)")
    parser.add_argument("-p", "--partition-out", required=True, help="Path to write RAxML partition file")
    return parser.parse_args()

def extract_subset_table(text):
    match = re.search(r"Subset\s+\|.*?Scheme Description", text, re.DOTALL)
    if not match:
        raise RuntimeError("Could not find subset model section.")
    return match.group(0).strip().splitlines()

def extract_partition_block(text):
    match = re.search(r"RaxML-style partition definitions.*?(DNA,.*?)\n\n", text, re.DOTALL)
    if not match:
        raise RuntimeError("Could not extract RAxML partition block.")
    return [line.strip() for line in match.group(1).splitlines() if line.strip().startswith("DNA,")]

def main():
    args = parse_args()

    with open(args.input) as f:
        text = f.read()

    # Print subset table
    subset_table = extract_subset_table(text)
    print("\nSubset table:\n")
    for line in subset_table:
        print(line)

    # Present model options
    print("\nSelect a rate heterogeneity model for RAxML:")
    for i, (pf_model, raxml_model) in enumerate(MODEL_MAP.items(), 1):
        print(f"{i}: {pf_model} → {raxml_model}")

    try:
        choice = int(input("\nPick model (enter number): "))
        pf_model = list(MODEL_MAP.keys())[choice - 1]
        raxml_model = MODEL_MAP[pf_model]
    except (IndexError, ValueError):
        raise SystemExit("Invalid selection.")

    # Write selected model
    with open(args.model_out, "w") as f:
        f.write(raxml_model + "\n")

    # Write partition file
    partition_lines = extract_partition_block(text)
    with open(args.partition_out, "w") as f:
        for line in partition_lines:
            f.write(line + "\n")

if __name__ == "__main__":
    main()
