# Main target
all: outputs/08_tree_collapsed/auriculella_species_collapsed.tre

# Step 01A: Download Auriculella
# Genbank has the same sequences
outputs/01_download_files/auriculella.json:
	mkdir -p outputs/01_download_files
	curl -s "https://portal.boldsystems.org/api/query/preprocessor?query=tax:Auriculella" \
	| jq -r '.successful_terms[0].matched' \
	| xargs -I{} curl -s "https://portal.boldsystems.org/api/query?query={}" \
	| jq -r '.query_id' \
	| xargs -I{} curl -s -o outputs/01_download_files/auriculella.json \
	"https://portal.boldsystems.org/api/documents/{}/download?format=json"

# Step 01B: Download Tornatellidinae
# Use this for rooting
outputs/01_download_files/tornatellidinae.json:
	mkdir -p outputs/01_download_files
	curl -s "https://portal.boldsystems.org/api/query/preprocessor?query=tax:Tornatellidinae" \
	| jq -r '.successful_terms[0].matched' \
	| xargs -I{} curl -s "https://portal.boldsystems.org/api/query?query={}" \
	| jq -r '.query_id' \
	| xargs -I{} curl -s -o outputs/01_download_files/tornatellidinae.json \
	"https://portal.boldsystems.org/api/documents/{}/download?format=json"

# Step 01C: Concatenate both JSON files
outputs/01_download_files/auriculella_tornatellidinae.json: outputs/01_download_files/auriculella.json outputs/01_download_files/tornatellidinae.json
	cat outputs/01_download_files/auriculella.json outputs/01_download_files/tornatellidinae.json > outputs/01_download_files/auriculella_tornatellidinae.json

# Step 02: Convert to FASTA (filtered and summarized)
outputs/02_fasta/done: outputs/01_download_files/auriculella_tornatellidinae.json
	mkdir -p outputs/02_fasta
	./scripts/json_to_fasta.py outputs/01_download_files/auriculella_tornatellidinae.json outputs/02_fasta
	touch outputs/02_fasta/done

# Step 03: Align each locus
outputs/03_alignment/done: outputs/02_fasta/done
	mkdir -p outputs/03_alignment
	for file in outputs/02_fasta/*.fas; do \
		base=$$(basename $$file .fas); \
		mafft --auto $$file > outputs/03_alignment/$$base.aln.fas; \
	done
	touch outputs/03_alignment/done

# Step 04: Build supermatrix from aligned loci
outputs/04_supermatrix/FcC_supermatrix.fas: outputs/03_alignment/done
	mkdir -p outputs/04_supermatrix
	cp outputs/03_alignment/*.aln.fas outputs/04_supermatrix/
	cd outputs/04_supermatrix && perl ../../scripts/FASconCAT-G_v1.02.pl -s

# Step 05: Build tree from the supermatrix
outputs/05_tree/auriculella_phylogeny.tre: outputs/04_supermatrix/FcC_supermatrix.fas
	mkdir -p outputs/05_tree
	iqtree -s outputs/04_supermatrix/FcC_supermatrix.fas -nt AUTO -bb 1000 -pre outputs/05_tree/auriculella_phylogeny
	mv outputs/05_tree/auriculella_phylogeny.treefile outputs/05_tree/auriculella_phylogeny.tre

# Step 06: Re-root tree using MRCA of three outgroup tips
outputs/06_tree_rooted/auriculella_rooted.tre: outputs/05_tree/auriculella_phylogeny.tre
	mkdir -p outputs/06_tree_rooted
	nw_reroot $< ACHA002-20_HLS20002_Tornatellides_sp ACHA216-22_PCMB54681_Tornatellides_subperforatus ACHA001-20_HLS20486_Tornatellaria_sp > $@

# Step 07: Rename tips in tree to species only
outputs/07_tree_species_renamed/auriculella_species.tre: outputs/06_tree_rooted/auriculella_rooted.tre
	mkdir -p outputs/07_tree_species_renamed
	./scripts/rename_tree_tips_by_species.py $< $@

# Step 08: Collapse identical species labels
outputs/08_tree_collapsed/auriculella_species_collapsed.tre: outputs/07_tree_species_renamed/auriculella_species.tre
	mkdir -p outputs/08_tree_collapsed
	nw_condense $< > $@