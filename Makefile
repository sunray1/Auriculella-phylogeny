# Main target
all: outputs/09_tree_collapsed/auriculella_species_collapsed.tre

# Tree Models
IQTREE_LABEL_NO_PART = iqtree_no_partitioning_without_site_removal
IQTREE_LABEL_SITE_REMOVED = iqtree_no_partitioning_with_site_removal
IQTREE_LABEL_PART = iqtree_locus_partitioning_without_site_removal
IQTREE_LABEL_PART_REMOVED = iqtree_locus_partitioning_with_site_removal
IQTREE_LABEL_CODON = iqtree_locus_codon_partition_without_site_removal
RAXML_LABEL_PART = raxml_partitioning_without_site_removal
RAXML_LABEL_PART_REMOVED = raxml_partitioning_with_site_removal

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
: outputs/03_alignment/done \
verify_coi
	mkdir -p outputs/04_supermatrix/original
	mkdir -p outputs/04_supermatrix/filtered

	cp outputs/03_alignment/*.aln.fas outputs/04_supermatrix/original/

	for aln in outputs/04_supermatrix/original/*.aln.fas; do \
		base=$$(basename $$aln); \
		clipkit $$aln -o outputs/04_supermatrix/filtered/$$base -m kpi-gappy; \
	done

	@if ls outputs/04_supermatrix/original/FcC* 1> /dev/null 2>&1; then \
		rm outputs/04_supermatrix/original/FcC*; \
	fi
	@if ls outputs/04_supermatrix/filtered/FcC* 1> /dev/null 2>&1; then \
		rm outputs/04_supermatrix/filtered/FcC*; \
	fi
	cd outputs/04_supermatrix/original && perl ../../../scripts/FASconCAT-G_v1.02.pl -s
	cd outputs/04_supermatrix/filtered && perl ../../../scripts/FASconCAT-G_v1.02.pl -s

# Step 05: PartitionFinder2 for RAxML

outputs/05_partitionfinder/original/analysis/best_scheme.txt \
outputs/05_partitionfinder/filtered/analysis/best_scheme.txt \
: outputs/04_supermatrix/original/FcC_supermatrix.fas \
outputs/04_supermatrix/filtered/FcC_supermatrix.fas
    # original alignment
	@if [ -d outputs/05_partitionfinder/original ]; then \
		rm -rf outputs/05_partitionfinder/original; \
	fi
	mkdir -p outputs/05_partitionfinder/original
	python scripts/fasta_to_phylip.py -i outputs/04_supermatrix/original/FcC_supermatrix.fas -o outputs/05_partitionfinder/original/FcC_supermatrix.phy
	
	# Make config file
	awk -F '\t' '\
	  $$1 ~ /\.aln$$/ {\
	    gene = $$1; sub(/\.aln$$/, "", gene);\
	    start = $$2; end = $$3;\
	    if (gene == "COI-5P") {\
	      print gene "_pos1 = " start "-" end "\\3;";\
	      print gene "_pos2 = " start+1 "-" end "\\3;";\
	      print gene "_pos3 = " start+2 "-" end "\\3;";\
	    } else {\
	      print gene " = " start "-" end ";";\
	    }\
	  }' outputs/04_supermatrix/original/FcC_info.xls > outputs/05_partitionfinder/original/data_blocks.txt
	printf "%s\n" \
"# ALIGNMENT FILE #" \
"alignment = FcC_supermatrix.phy;" \
"# BRANCHLENGTHS #" \
"branchlengths = unlinked;" \
"# MODELS OF EVOLUTION #" \
"models = all;" \
"model_selection = aicc;" \
"# DATA BLOCKS #" \
"[data_blocks]" \
>> outputs/05_partitionfinder/original/partition_finder.cfg
	cat outputs/05_partitionfinder/original/data_blocks.txt >> outputs/05_partitionfinder/original/partition_finder.cfg
	printf "%s\n" \
"# SCHEMES #" \
"[schemes]" \
"search = greedy;" \
>> outputs/05_partitionfinder/original/partition_finder.cfg
	rm outputs/05_partitionfinder/original/data_blocks.txt
	python scripts/partitionfinder-2.1.1/PartitionFinder.py outputs/05_partitionfinder/original --raxml --force-restart
	
	# Filtered alignment
	@if [ -d outputs/05_partitionfinder/filtered ]; then \
		rm -rf outputs/05_partitionfinder/filtered; \
	fi
	mkdir -p outputs/05_partitionfinder/filtered
	python scripts/fasta_to_phylip.py -i outputs/04_supermatrix/filtered/FcC_supermatrix.fas -o outputs/05_partitionfinder/filtered/FcC_supermatrix.phy
	# Make config file
	awk -F '\t' '\
	  $$1 ~ /\.aln$$/ {\
	    gene = $$1; sub(/\.aln$$/, "", gene);\
	    start = $$2; end = $$3;\
	    if (gene == "COI-5P") {\
	      print gene "_pos1 = " start "-" end "\\3;";\
	      print gene "_pos2 = " start+1 "-" end "\\3;";\
	      print gene "_pos3 = " start+2 "-" end "\\3;";\
	    } else {\
	      print gene " = " start "-" end ";";\
	    }\
	  }' outputs/04_supermatrix/filtered/FcC_info.xls > outputs/05_partitionfinder/filtered/data_blocks.txt
	printf "%s\n" \
"# ALIGNMENT FILE #" \
"alignment = FcC_supermatrix.phy;" \
"# BRANCHLENGTHS #" \
"branchlengths = unlinked;" \
"# MODELS OF EVOLUTION #" \
"models = all;" \
"model_selection = aicc;" \
"# DATA BLOCKS #" \
"[data_blocks]" \
>> outputs/05_partitionfinder/filtered/partition_finder.cfg
	cat outputs/05_partitionfinder/filtered/data_blocks.txt >> outputs/05_partitionfinder/filtered/partition_finder.cfg
	printf "%s\n" \
"# SCHEMES #" \
"[schemes]" \
"search = greedy;" \
>> outputs/05_partitionfinder/filtered/partition_finder.cfg
	rm outputs/05_partitionfinder/filtered/data_blocks.txt
	python scripts/partitionfinder-2.1.1/PartitionFinder.py outputs/05_partitionfinder/filtered --raxml --force-restart

# Step 06: Build trees from the supermatrix

# RAxML: Build tree with partitioning (PartitionFinder2) and without site removal
outputs/06_tree/$(RAXML_LABEL_PART)/RAxML_info.$(RAXML_LABEL_PART) \
: outputs/04_supermatrix/original/FcC_supermatrix.fas \
outputs/05_partitionfinder/original/analysis/best_scheme.txt

	@if [ -d outputs/06_tree/$(RAXML_LABEL_PART) ]; then \
		rm -rf outputs/06_tree/$(RAXML_LABEL_PART); \
	fi
	mkdir -p outputs/06_tree/$(RAXML_LABEL_PART)
	
	# Get best model and make the partition file
	python scripts/extract_model_and_partitions.py \
		-i outputs/05_partitionfinder/original/analysis/best_scheme.txt \
		-m outputs/06_tree/$(RAXML_LABEL_PART)/selected_model.txt \
		-p outputs/06_tree/$(RAXML_LABEL_PART)/partition_file.txt

	cd outputs/06_tree/$(RAXML_LABEL_PART) && \
	MODEL=$$(cat selected_model.txt) && \
	SEED=$$(date +%s) && \
	raxmlHPC -s ../../04_supermatrix/original/FcC_supermatrix.fas \
	   -f a -m $$MODEL -q partition_file.txt \
	   -p $$SEED -x $$SEED -# 100 \
	   -n $(RAXML_LABEL_PART)
	mv outputs/06_tree/$(RAXML_LABEL_PART)/RAxML_bipartitions.$(RAXML_LABEL_PART) outputs/06_tree/$(RAXML_LABEL_PART)/RAxML_bipartitions.$(RAXML_LABEL_PART).tre

# RAxML: Build tree with partitioning (PartitionFinder2) and with site removal
outputs/06_tree/$(RAXML_LABEL_PART_REMOVED)/RAxML_info.$(RAXML_LABEL_PART_REMOVED) \
: outputs/04_supermatrix/filtered/FcC_supermatrix.fas \
outputs/05_partitionfinder/filtered/analysis/best_scheme.txt

	@if [ -d outputs/06_tree/$(RAXML_LABEL_PART_REMOVED) ]; then \
		rm -rf outputs/06_tree/$(RAXML_LABEL_PART_REMOVED); \
	fi
	mkdir -p outputs/06_tree/$(RAXML_LABEL_PART_REMOVED)
	
	# Get best model and make the partition file
	python scripts/extract_model_and_partitions.py \
		-i outputs/05_partitionfinder/filtered/analysis/best_scheme.txt \
		-m outputs/06_tree/$(RAXML_LABEL_PART_REMOVED)/selected_model.txt \
		-p outputs/06_tree/$(RAXML_LABEL_PART_REMOVED)/partition_file.txt

	cd outputs/06_tree/$(RAXML_LABEL_PART_REMOVED) && \
	MODEL=$$(cat selected_model.txt) && \
	SEED=$$(date +%s) && \
	raxmlHPC -s ../../04_supermatrix/filtered/FcC_supermatrix.fas \
	   -f a -m $$MODEL -q partition_file.txt \
	   -p $$SEED -x $$SEED -# 100 \
	   -n $(RAXML_LABEL_PART_REMOVED)
	mv outputs/06_tree/$(RAXML_LABEL_PART_REMOVED)/RAxML_bipartitions.$(RAXML_LABEL_PART_REMOVED) outputs/06_tree/$(RAXML_LABEL_PART_REMOVED)/RAxML_bipartitions.$(RAXML_LABEL_PART_REMOVED).tre

# IQ-Tree: Build tree without partitioning and without site removal
outputs/06_tree/$(IQTREE_LABEL_NO_PART)/auriculella_phylogeny_$(IQTREE_LABEL_NO_PART).iqtree: outputs/04_supermatrix/original/FcC_supermatrix.fas
	mkdir -p outputs/06_tree/$(IQTREE_LABEL_NO_PART)
	iqtree -s outputs/04_supermatrix/original/FcC_supermatrix.fas \
	   -redo -nt AUTO -bb 1000 \
	   -pre outputs/06_tree/$(IQTREE_LABEL_NO_PART)/auriculella_phylogeny_$(IQTREE_LABEL_NO_PART)
	mv outputs/06_tree/$(IQTREE_LABEL_NO_PART)/auriculella_phylogeny_$(IQTREE_LABEL_NO_PART).treefile outputs/06_tree/$(IQTREE_LABEL_NO_PART)/auriculella_phylogeny_$(IQTREE_LABEL_NO_PART).tre

# IQ-Tree: Build tree without partitioning and with site removal
outputs/06_tree/$(IQTREE_LABEL_SITE_REMOVED)/auriculella_phylogeny_$(IQTREE_LABEL_SITE_REMOVED).iqtree: outputs/04_supermatrix/filtered/FcC_supermatrix.fas
	mkdir -p outputs/06_tree/$(IQTREE_LABEL_SITE_REMOVED)
	iqtree -s outputs/04_supermatrix/filtered/FcC_supermatrix.fas \
	   -redo -nt AUTO -bb 1000 \
	   -pre outputs/06_tree/$(IQTREE_LABEL_SITE_REMOVED)/auriculella_phylogeny_$(IQTREE_LABEL_SITE_REMOVED)
	mv outputs/06_tree/$(IQTREE_LABEL_SITE_REMOVED)/auriculella_phylogeny_$(IQTREE_LABEL_SITE_REMOVED).treefile outputs/06_tree/$(IQTREE_LABEL_SITE_REMOVED)/auriculella_phylogeny_$(IQTREE_LABEL_SITE_REMOVED).tre

# IQ-Tree: Build tree partitioning by locus and without site removal
outputs/06_tree/$(IQTREE_LABEL_PART)/auriculella_phylogeny_$(IQTREE_LABEL_PART).iqtree: outputs/04_supermatrix/original/FcC_supermatrix.fas outputs/04_supermatrix/original/FcC_info.xls
	mkdir -p outputs/06_tree/$(IQTREE_LABEL_PART)
	awk -F '\t' ' \
		$$1 ~ /\.aln$$/ { \
			sub(/\.aln/, "", $$1); \
			print "DNA, " $$1 " = " $$2 "-" $$3; \
		}' outputs/04_supermatrix/original/FcC_info.xls > outputs/06_tree/$(IQTREE_LABEL_PART)/auriculella_partitions.txt
	iqtree -s outputs/04_supermatrix/original/FcC_supermatrix.fas \
	   -p outputs/06_tree/$(IQTREE_LABEL_PART)/auriculella_partitions.txt \
	   -redo -nt AUTO -bb 1000 \
	   -pre outputs/06_tree/$(IQTREE_LABEL_PART)/auriculella_phylogeny_$(IQTREE_LABEL_PART)
	mv outputs/06_tree/$(IQTREE_LABEL_PART)/auriculella_phylogeny_$(IQTREE_LABEL_PART).treefile outputs/06_tree/$(IQTREE_LABEL_PART)/auriculella_phylogeny_$(IQTREE_LABEL_PART).tre

# IQ-Tree: Build tree partitioning by locus and with site removal
outputs/06_tree/$(IQTREE_LABEL_PART_REMOVED)/auriculella_phylogeny_$(IQTREE_LABEL_PART_REMOVED).iqtree: outputs/04_supermatrix/filtered/FcC_supermatrix.fas outputs/04_supermatrix/filtered/FcC_info.xls
	mkdir -p outputs/06_tree/$(IQTREE_LABEL_PART_REMOVED)
	awk -F '\t' ' \
		$$1 ~ /\.aln$$/ { \
			sub(/\.aln/, "", $$1); \
			print "DNA, " $$1 " = " $$2 "-" $$3; \
		}' outputs/04_supermatrix/filtered/FcC_info.xls > outputs/06_tree/$(IQTREE_LABEL_PART_REMOVED)/auriculella_partitions.txt
	iqtree -s outputs/04_supermatrix/filtered/FcC_supermatrix.fas \
	   -p outputs/06_tree/$(IQTREE_LABEL_PART_REMOVED)/auriculella_partitions.txt \
	   -redo -nt AUTO -bb 1000 \
	   -pre outputs/06_tree/$(IQTREE_LABEL_PART_REMOVED)/auriculella_phylogeny_$(IQTREE_LABEL_PART_REMOVED)
	mv outputs/06_tree/$(IQTREE_LABEL_PART_REMOVED)/auriculella_phylogeny_$(IQTREE_LABEL_PART_REMOVED).treefile outputs/06_tree/$(IQTREE_LABEL_PART_REMOVED)/auriculella_phylogeny_$(IQTREE_LABEL_PART_REMOVED).tre

# IQ-Tree: Build tree partitioning by locus and COI codon and without site removal
outputs/06_tree/$(IQTREE_LABEL_CODON)/auriculella_phylogeny_$(IQTREE_LABEL_CODON).iqtree: outputs/04_supermatrix/original/FcC_supermatrix.fas outputs/04_supermatrix/original/FcC_info.xls
	mkdir -p outputs/06_tree/$(IQTREE_LABEL_CODON)
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
		}' outputs/04_supermatrix/original/FcC_info.xls > outputs/06_tree/$(IQTREE_LABEL_CODON)/auriculella_partitions.txt
	iqtree -s outputs/04_supermatrix/original/FcC_supermatrix.fas \
	   -p outputs/06_tree/$(IQTREE_LABEL_CODON)/auriculella_partitions.txt \
	   -redo -nt AUTO -bb 1000 \
	   -pre outputs/06_tree/$(IQTREE_LABEL_CODON)/auriculella_phylogeny_$(IQTREE_LABEL_CODON)
	mv outputs/06_tree/$(IQTREE_LABEL_CODON)/auriculella_phylogeny_$(IQTREE_LABEL_CODON).treefile outputs/06_tree/$(IQTREE_LABEL_CODON)/auriculella_phylogeny_$(IQTREE_LABEL_CODON).tre
	
outputs/06_tree/iqtree_log_likelihoods.tsv \
outputs/06_tree/raxml_log_likelihoods.tsv: \
	outputs/06_tree/$(RAXML_LABEL_PART)/RAxML_info.$(RAXML_LABEL_PART) \
	outputs/06_tree/$(RAXML_LABEL_PART_REMOVED)/RAxML_info.$(RAXML_LABEL_PART_REMOVED) \
	outputs/06_tree/$(IQTREE_LABEL_NO_PART)/auriculella_phylogeny_$(IQTREE_LABEL_NO_PART).iqtree \
	outputs/06_tree/$(IQTREE_LABEL_SITE_REMOVED)/auriculella_phylogeny_$(IQTREE_LABEL_SITE_REMOVED).iqtree \
	outputs/06_tree/$(IQTREE_LABEL_PART)/auriculella_phylogeny_$(IQTREE_LABEL_PART).iqtree \
	outputs/06_tree/$(IQTREE_LABEL_PART_REMOVED)/auriculella_phylogeny_$(IQTREE_LABEL_PART_REMOVED).iqtree \
	outputs/06_tree/$(IQTREE_LABEL_CODON)/auriculella_phylogeny_$(IQTREE_LABEL_CODON).iqtree

	# IQ-TREE log-likelihoods
	printf "label\tlog_likelihood\n" > outputs/06_tree/iqtree_log_likelihoods.tsv
	for label in \
		$(IQTREE_LABEL_NO_PART) \
		$(IQTREE_LABEL_SITE_REMOVED) \
		$(IQTREE_LABEL_PART) \
		$(IQTREE_LABEL_PART_REMOVED) \
		$(IQTREE_LABEL_CODON) \
		$(IQTREE_LABEL_CODON_REMOVED); do \
		ll=$$(grep "Log-likelihood of consensus tree" outputs/06_tree/$$label/auriculella_phylogeny_$$label.iqtree | awk '{print $$NF}'); \
		printf "%s\t%s\n" "$$label" "$$ll" >> outputs/06_tree/iqtree_log_likelihoods.tsv; \
	done

	# RAxML log-likelihoods
	printf "label\tlog_likelihood\n" > outputs/06_tree/raxml_log_likelihoods.tsv
	for label in \
		$(RAXML_LABEL_PART) \
		$(RAXML_LABEL_PART_REMOVED); do \
		ll=$$(grep "Final ML Optimization Likelihood" outputs/06_tree/$$label/RAxML_info.$$label | awk '{print $$NF}'); \
		printf "%s\t%s\n" "$$label" "$$ll" >> outputs/06_tree/raxml_log_likelihoods.tsv; \
	done

#Human input to choose the tree to move forward with	
outputs/06_tree/selected_tree_label.txt: \
outputs/06_tree/iqtree_log_likelihoods.tsv \
outputs/06_tree/raxml_log_likelihoods.tsv
	@bash -c ' \
		echo "Please choose which tree to use for downstream steps:"; \
		echo "1) RAxML   - $(RAXML_LABEL_PART)"; \
		echo "2) RAxML   - $(RAXML_LABEL_PART_REMOVED)"; \
		echo "3) IQTREE - $(IQTREE_LABEL_NO_PART)"; \
		echo "4) IQTREE - $(IQTREE_LABEL_SITE_REMOVED)"; \
		echo "5) IQTREE - $(IQTREE_LABEL_PART)"; \
		echo "6) IQTREE - $(IQTREE_LABEL_PART_REMOVED)"; \
		echo "7) IQTREE - $(IQTREE_LABEL_CODON)"; \
		read -p "Enter the number of your choice: " choice; \
		case $$choice in \
			1) echo "$(RAXML_LABEL_PART)" > outputs/06_tree/selected_tree_label.txt ;; \
			2) echo "$(RAXML_LABEL_PART_REMOVED)" > outputs/06_tree/selected_tree_label.txt ;; \
			3) echo "$(IQTREE_LABEL_NO_PART)" > outputs/06_tree/selected_tree_label.txt ;; \
			4) echo "$(IQTREE_LABEL_SITE_REMOVED)" > outputs/06_tree/selected_tree_label.txt ;; \
			5) echo "$(IQTREE_LABEL_PART)" > outputs/06_tree/selected_tree_label.txt ;; \
			6) echo "$(IQTREE_LABEL_PART_REMOVED)" > outputs/06_tree/selected_tree_label.txt ;; \
			7) echo "$(IQTREE_LABEL_CODON)" > outputs/06_tree/selected_tree_label.txt ;; \
			*) echo "Invalid choice"; exit 1 ;; \
		esac; \
		echo "Saved selected label to outputs/06_tree/selected_tree_label.txt"; \
	'

# Step 07: Re-root tree using MRCA of three outgroup tips
outputs/07_tree_rooted/done: outputs/06_tree/selected_tree_label.txt
	mkdir -p outputs/07_tree_rooted

	# Use user-selected tree label
	label=$$(cat $<); \
	tree_path=outputs/06_tree/$$label/auriculella_phylogeny_$$label.tre; \
	base=$$(basename $$tree_path .tre); \
	nw_reroot $$tree_path ACHA002-20_HLS20002_Tornatellides_sp ACHA216-22_PCMB54681_Tornatellides_subperforatus ACHA001-20_HLS20486_Tornatellaria_sp \
	> outputs/07_tree_rooted/$$base.rooted.tre

	touch $@


# Step 08: Rename tips in tree to species only
outputs/08_tree_species_renamed/auriculella_species.tre: outputs/07_tree_rooted/done
	mkdir -p outputs/07_tree_species_renamed

	best=$$(tail -n +2 outputs/06_tree/log_likelihood_summary.tsv | sort -k2,2n | head -n1 | cut -f1); \
	rooted_tree=outputs/07_tree_rooted/auriculella_phylogeny_$$best.rooted.tre; \
	./scripts/rename_tree_tips_by_species.py $$rooted_tree $@

# Step 09: Collapse identical species labels
outputs/09_tree_collapsed/auriculella_species_collapsed.tre: outputs/08_tree_species_renamed/auriculella_species.tre
	mkdir -p outputs/09_tree_collapsed
	nw_condense $< > $@