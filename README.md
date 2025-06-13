# Auriculella Phylogenomics Pipeline

This repository implements a reproducible bioinformatics workflow for generating a species-level phylogeny of the Hawaiian land snail genus *Auriculella* using publicly available sequence data from the Barcode of Life Data Systems (BOLD). The pipeline includes automated data retrieval, sequence processing, alignment, supermatrix construction, phylogenetic inference, rooting, and tip collapse.

The workflow is executed via GNU Make and is fully containerized within a Conda environment to ensure reproducibility across platforms.

## Reference

This pipeline conceptually replicates the phylogenetic analysis presented in:

> Yeung, N. W., & Hayes, K. A. (2020).  
> Overlooked but not forgotten: the first new extant species of Hawaiian land snail described in 60 years, *Auriculella gagneorum* sp. nov. (Achatinellidae, Auriculellinae).  
> *ZooKeys*, 915: 137–162.  
> [https://doi.org/10.3897/zookeys.915.50669](https://doi.org/10.3897/zookeys.915.50669)

## Overview

The pipeline performs the following steps:

1. Query BOLD for sequences from *Auriculella* and relevant outgroup taxa (*Tornatellidinae*).
2. Write filtered sequences to FASTA format, grouped by marker (e.g., 16S, COI).
3. Perform multiple sequence alignment with MAFFT for each locus.
    - Includes a manual verification pause where the user can inspect and adjust alignments (e.g., check reading frames) before concatenation.
4. Concatenate alignments into a supermatrix using FASconCAT-G.
5. Run maximum likelihood inference using IQ-TREE across five different model schemes:

    * No partitioning
    * No partitioning with site filtering (removing gappy (>.9) and uninformative sites)
    * Partitioning by locus
    * Partitioning by locus with site filtering
    * Partitioning by locus with COI codon positions

    Log-likelihoods for all trees are computed, and the best-supported topology is selected automatically.

6. Root the tree using an explicit outgroup (MRCA of selected *Tornatellidinae* tips).
7. Rename tip labels to species names.
8. Collapse redundant tips by species to produce a simplified tree.

## Requirements

- GNU Make
- Conda (miniconda or anaconda)
- Bash-compatible shell

The following software tools are managed within the Conda environment:

- mafft
- iqtree
- clipkit
- perl (for FASconCAT-G)
- python (>=3.7) with `biopython`
- jq
- curl
- nw_utils (for `nw_reroot`, `nw_condense`)

## Setup

### Clone the repository

```bash
git clone https://github.com/yourusername/auriculella-phylogeny.git
cd auriculella-phylogeny
```

### Create and activate the Conda environment
```bash
conda env create -f environment.yml
conda activate auriculella
```

### Usage
To execute the full pipeline:

```bash
make
```

This will produce a collapsed, species-labeled tree in:

```bash
outputs/08_tree_collapsed/auriculella_species_collapsed.tre
```

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
└── 08_tree_collapsed/         # Final collapsed phylogeny
```

### Notes
* Only high-confidence records are retained (e.g., species-level IDs excluding Auriculella_sp, and excluding GenBank records with processid beginning GB).
* The rooting step explicitly uses three Tornatellidinae records defined in the Makefile.
* The tree collapsing step assumes monophyly of species and simplifies downstream visualization and interpretation.
