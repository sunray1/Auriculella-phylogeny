# Auriculella Phylogenomics Pipeline

This repository implements a reproducible bioinformatics workflow for generating a species-level phylogeny of the Hawaiian land snail genus *Auriculella* using publicly available sequence data from the Barcode of Life Data Systems (BOLD). The pipeline includes automated data retrieval, sequence processing, alignment, supermatrix construction, phylogenetic inference, rooting, and tip collapse.

The workflow is executed via GNU Make and is fully containerized within a Conda environment to ensure reproducibility across platforms.

## Overview

The pipeline performs the following steps:

1. Query BOLD for sequences from *Auriculella* and relevant outgroup taxa (*Tornatellidinae*).
2. Parse and filter BOLD records, excluding ambiguous taxa and GenBank imports.
3. Write filtered sequences to FASTA format, grouped by marker (e.g., 16S, COI).
4. Perform multiple sequence alignment with MAFFT for each locus.
5. Concatenate alignments into a supermatrix using FASconCAT-G.
6. Infer a maximum likelihood phylogeny using IQ-TREE.
7. Root the tree using an explicit outgroup (MRCA of selected *Tornatellidinae* tips).
8. Rename tip labels to species names.
9. Collapse redundant tips by species to produce a simplified tree.

## Requirements

- GNU Make
- Conda (miniconda or anaconda)
- Bash-compatible shell

The following software tools are managed within the Conda environment:

- mafft
- iqtree
- perl (for FASconCAT-G)
- python (>=3.7) with `biopython`
- jq
- curl
- nw_utils (for `nw_reroot`, `nw_condense`)

## Setup

### Clone the repository

```bash
git clone https://github.com/yourusername/auriculella-phylogeny.git
cd auriculella-phylogeny```

### Create and activate the Conda environment
```bash
conda env create -f environment.yml
conda activate auriculella```

### Usage
To execute the full pipeline:

```bash
make```

This will produce a collapsed, species-labeled tree in:

```bash
outputs/08_tree_collapsed/auriculella_species_collapsed.tre```

Intermediate results, logs, and alignments are organized under the outputs/ directory by processing stage.

### Directory Structure
```graphql
outputs/
├── 01_download/               # Raw JSON data from BOLD
├── 02_fasta/                  # Filtered FASTA sequences by marker
├── 03_alignment/              # MAFFT alignments
├── 04_supermatrix/            # Concatenated alignments
├── 05_tree/                   # IQ-TREE inference
├── 06_tree_rooted/            # Rooted tree (MRCA of outgroup)
├── 07_tree_species_renamed/   # Tree with species-only tip labels
└── 08_tree_collapsed/         # Final collapsed phylogeny```

### Notes
* Only high-confidence records are retained (e.g., species-level IDs excluding Auriculella_sp, and excluding GenBank records with processid beginning GB).
* The rooting step explicitly uses three Tornatellidinae records defined in the Makefile.
* The tree collapsing step assumes monophyly of species and simplifies downstream visualization and interpretation.