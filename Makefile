# Main target
all: outputs/08_tree_collapsed/auriculella_species_collapsed.tre

# Tree Models
TREE_LABEL_NO_PART = no_partitioning_without_site_removal
TREE_LABEL_SITE_REMOVED = no_partitioning_with_site_removal
TREE_LABEL_PART = locus_partitioning_without_site_removal
TREE_LABEL_PART_REMOVED = locus_partitioning_with_site_removal
TREE_LABEL_CODON = locus_codon_partition_without_site_removal

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
	
# Human verification of alignments
verify_coi:
	@bash -c ' \
		echo "Please verify outputs/03_alignment/COI-5P.aln.fas reading frame and other locus alignments"; \
		read -p "Press Enter to continue..."; \
	'

# Step 04: Build supermatrix from aligned loci
outputs/04_supermatrix/filtered/FcC_supermatrix.fas \
outputs/04_supermatrix/original/FcC_supermatrix.fas \
: outputs/03_alignment/done verify_coi
	mkdir -p outputs/04_supermatrix/original
	mkdir -p outputs/04_supermatrix/filtered

	cp outputs/03_alignment/*.aln.fas outputs/04_supermatrix/original/

	for aln in outputs/04_supermatrix/original/*.aln.fas; do \
		base=$$(basename $$aln); \
		clipkit $$aln -o outputs/04_supermatrix/filtered/$$base -m kpi-gappy; \
	done

	cd outputs/04_supermatrix/original && perl ../../../scripts/FASconCAT-G_v1.02.pl -s
	cd outputs/04_supermatrix/filtered && perl ../../../scripts/FASconCAT-G_v1.02.pl -s

# Step 05: Build trees from the supermatrix

# Build tree without partitioning and without site removal
outputs/05_tree/$(TREE_LABEL_NO_PART)/auriculella_phylogeny_$(TREE_LABEL_NO_PART).iqtree: outputs/04_supermatrix/original/FcC_supermatrix.fas
	mkdir -p outputs/05_tree/$(TREE_LABEL_NO_PART)
	iqtree -s outputs/04_supermatrix/original/FcC_supermatrix.fas \
	   -redo -nt AUTO -bb 1000 \
	   -pre outputs/05_tree/$(TREE_LABEL_NO_PART)/auriculella_phylogeny_$(TREE_LABEL_NO_PART)
	mv outputs/05_tree/$(TREE_LABEL_NO_PART)/auriculella_phylogeny_$(TREE_LABEL_NO_PART).treefile outputs/05_tree/$(TREE_LABEL_NO_PART)/auriculella_phylogeny_$(TREE_LABEL_NO_PART).tre

# Build tree without partitioning and with site removal
outputs/05_tree/$(TREE_LABEL_SITE_REMOVED)/auriculella_phylogeny_$(TREE_LABEL_SITE_REMOVED).iqtree: outputs/04_supermatrix/filtered/FcC_supermatrix.fas
	mkdir -p outputs/05_tree/$(TREE_LABEL_SITE_REMOVED)
	iqtree -s outputs/04_supermatrix/filtered/FcC_supermatrix.fas \
	   -redo -nt AUTO -bb 1000 \
	   -pre outputs/05_tree/$(TREE_LABEL_SITE_REMOVED)/auriculella_phylogeny_$(TREE_LABEL_SITE_REMOVED)
	mv outputs/05_tree/$(TREE_LABEL_SITE_REMOVED)/auriculella_phylogeny_$(TREE_LABEL_SITE_REMOVED).treefile outputs/05_tree/$(TREE_LABEL_SITE_REMOVED)/auriculella_phylogeny_$(TREE_LABEL_SITE_REMOVED).tre

# Build tree partitioning by locus and without site removal
outputs/05_tree/$(TREE_LABEL_PART)/auriculella_phylogeny_$(TREE_LABEL_PART).iqtree: outputs/04_supermatrix/original/FcC_supermatrix.fas outputs/04_supermatrix/original/FcC_info.xls
	mkdir -p outputs/05_tree/$(TREE_LABEL_PART)
	awk -F '\t' ' \
		$$1 ~ /\.aln$$/ { \
			sub(/\.aln/, "", $$1); \
			print "DNA, " $$1 " = " $$2 "-" $$3; \
		}' outputs/04_supermatrix/original/FcC_info.xls > outputs/05_tree/$(TREE_LABEL_PART)/auriculella_partitions.txt
	iqtree -s outputs/04_supermatrix/original/FcC_supermatrix.fas \
	   -p outputs/05_tree/$(TREE_LABEL_PART)/auriculella_partitions.txt \
	   -redo -nt AUTO -bb 1000 \
	   -pre outputs/05_tree/$(TREE_LABEL_PART)/auriculella_phylogeny_$(TREE_LABEL_PART)
	mv outputs/05_tree/$(TREE_LABEL_PART)/auriculella_phylogeny_$(TREE_LABEL_PART).treefile outputs/05_tree/$(TREE_LABEL_PART)/auriculella_phylogeny_$(TREE_LABEL_PART).tre

# Build tree partitioning by locus and with site removal
outputs/05_tree/$(TREE_LABEL_PART_REMOVED)/auriculella_phylogeny_$(TREE_LABEL_PART_REMOVED).iqtree: outputs/04_supermatrix/filtered/FcC_supermatrix.fas outputs/04_supermatrix/filtered/FcC_info.xls
	mkdir -p outputs/05_tree/$(TREE_LABEL_PART_REMOVED)
	awk -F '\t' ' \
		$$1 ~ /\.aln$$/ { \
			sub(/\.aln/, "", $$1); \
			print "DNA, " $$1 " = " $$2 "-" $$3; \
		}' outputs/04_supermatrix/filtered/FcC_info.xls > outputs/05_tree/$(TREE_LABEL_PART_REMOVED)/auriculella_partitions.txt
	iqtree -s outputs/04_supermatrix/filtered/FcC_supermatrix.fas \
	   -p outputs/05_tree/$(TREE_LABEL_PART_REMOVED)/auriculella_partitions.txt \
	   -redo -nt AUTO -bb 1000 \
	   -pre outputs/05_tree/$(TREE_LABEL_PART_REMOVED)/auriculella_phylogeny_$(TREE_LABEL_PART_REMOVED)
	mv outputs/05_tree/$(TREE_LABEL_PART_REMOVED)/auriculella_phylogeny_$(TREE_LABEL_PART_REMOVED).treefile outputs/05_tree/$(TREE_LABEL_PART_REMOVED)/auriculella_phylogeny_$(TREE_LABEL_PART_REMOVED).tre

# Build tree partitioning by locus and COI codon and without site removal
outputs/05_tree/$(TREE_LABEL_CODON)/auriculella_phylogeny_$(TREE_LABEL_CODON).iqtree: outputs/04_supermatrix/original/FcC_supermatrix.fas outputs/04_supermatrix/original/FcC_info.xls
	mkdir -p outputs/05_tree/$(TREE_LABEL_CODON)
	awk -F '\t' ' \
		$$1 ~ /\.aln$$/ { \
			sub(/\.aln/, "", $$1); \
			start = $$2; end = $$3; \
			if ($$1 == "COI-5P") { \
				print "DNA, COI-5P_pos1 = " start "-" end "\\3"; \
				print "DNA, COI-5P_pos2 = " start+1 "-" end "\\3"; \
				print "DNA, COI-5P_pos3 = " start+2 "-" end "\\3"; \
			} else { \
				print "DNA, " $$1 " = " start "-" end; \
			} \
		}' outputs/04_supermatrix/original/FcC_info.xls > outputs/05_tree/$(TREE_LABEL_CODON)/auriculella_partitions.txt
	iqtree -s outputs/04_supermatrix/original/FcC_supermatrix.fas \
	   -p outputs/05_tree/$(TREE_LABEL_CODON)/auriculella_partitions.txt \
	   -redo -nt AUTO -bb 1000 \
	   -pre outputs/05_tree/$(TREE_LABEL_CODON)/auriculella_phylogeny_$(TREE_LABEL_CODON)
	mv outputs/05_tree/$(TREE_LABEL_CODON)/auriculella_phylogeny_$(TREE_LABEL_CODON).treefile outputs/05_tree/$(TREE_LABEL_CODON)/auriculella_phylogeny_$(TREE_LABEL_CODON).tre
	
outputs/05_tree/log_likelihood_summary.tsv: \
	outputs/05_tree/$(TREE_LABEL_NO_PART)/auriculella_phylogeny_$(TREE_LABEL_NO_PART).iqtree \
	outputs/05_tree/$(TREE_LABEL_SITE_REMOVED)/auriculella_phylogeny_$(TREE_LABEL_SITE_REMOVED).iqtree \
	outputs/05_tree/$(TREE_LABEL_PART)/auriculella_phylogeny_$(TREE_LABEL_PART).iqtree \
	outputs/05_tree/$(TREE_LABEL_PART_REMOVED)/auriculella_phylogeny_$(TREE_LABEL_PART_REMOVED).iqtree \
	outputs/05_tree/$(TREE_LABEL_CODON)/auriculella_phylogeny_$(TREE_LABEL_CODON).iqtree \

	printf "label\tlog_likelihood\n" > $@
	for label in \
		$(TREE_LABEL_NO_PART) \
		$(TREE_LABEL_SITE_REMOVED) \
		$(TREE_LABEL_PART) \
		$(TREE_LABEL_PART_REMOVED) \
		$(TREE_LABEL_CODON) \
		$(TREE_LABEL_CODON_REMOVED); do \
		ll=$$(grep "Log-likelihood of consensus tree" outputs/05_tree/$$label/auriculella_phylogeny_$$label.iqtree | awk '{print $$NF}'); \
		printf "%s\t%s\n" "$$label" "$$ll" >> $@; \
	done

# Step 06: Re-root tree using MRCA of three outgroup tips
outputs/06_tree_rooted/done: outputs/05_tree/log_likelihood_summary.tsv
	mkdir -p outputs/06_tree_rooted

	# Extract best label (highest log-likelihood = closest to 0)
	best=$$(tail -n +2 $< | sort -k2,2n | head -n1 | cut -f1); \
	tree_path=outputs/05_tree/$$best/auriculella_phylogeny_$$best.tre; \
	base=$$(basename $$tree_path .tre); \
	nw_reroot $$tree_path ACHA002-20_HLS20002_Tornatellides_sp ACHA216-22_PCMB54681_Tornatellides_subperforatus ACHA001-20_HLS20486_Tornatellaria_sp \
	> outputs/06_tree_rooted/$$base.rooted.tre

	touch $@


# Step 07: Rename tips in tree to species only
outputs/07_tree_species_renamed/auriculella_species.tre: outputs/06_tree_rooted/done
	mkdir -p outputs/07_tree_species_renamed

	best=$$(tail -n +2 outputs/05_tree/log_likelihood_summary.tsv | sort -k2,2n | head -n1 | cut -f1); \
	rooted_tree=outputs/06_tree_rooted/auriculella_phylogeny_$$best.rooted.tre; \
	./scripts/rename_tree_tips_by_species.py $$rooted_tree $@

# Step 08: Collapse identical species labels
outputs/08_tree_collapsed/auriculella_species_collapsed.tre: outputs/07_tree_species_renamed/auriculella_species.tre
	mkdir -p outputs/08_tree_collapsed
	nw_condense $< > $@