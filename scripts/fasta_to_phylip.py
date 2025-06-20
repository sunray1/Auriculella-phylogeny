import argparse
from Bio import AlignIO

def fasta_to_strict_phylip(input_fasta, output_phylip):
    alignment = AlignIO.read(input_fasta, "fasta")
    nseq = len(alignment)
    length = alignment.get_alignment_length()

    with open(output_phylip, "w") as out:
        out.write(f"{nseq} {length}\n")
        for record in alignment:
            name = record.id
            sequence = str(record.seq).replace("-", "N")  # Optional: standardize gaps
            out.write(f"{name} {sequence}\n")

def main():
    parser = argparse.ArgumentParser(description="Convert FASTA to PHYLIP (strict format).")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output PHYLIP file")
    args = parser.parse_args()

    # Perform the conversion
    fasta_to_strict_phylip(args.input, args.output)

if __name__ == "__main__":
    main()
