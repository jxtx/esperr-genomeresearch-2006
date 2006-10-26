SHELL=bash

# ---- Locations of various software ----------------------------------------

# PHAST binaries from Adam Siepel, needed to build initial substitution model
PH=~/phast/bin

# The ESPERR scripts (installed somewhere, in this case under your home directory)
RP=~/$(MACHTYPE)/bin

# The bx-python scripts (installed somewhere, in this case under your home directory)
BX=~/$(MACHTYPE)/bin

# ---- Other variables you might need to change -----------------------------

# The chromosomes of your reference species
CHROMOSOMES=1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y

# The substitution model for the tree and species you are working with (this
# is in PHAST format and contains the tree and rate matrix). If this is set 
# to 'random_1mb_windows.mod' the Makefile will attempt to build the model
# using PHAST. Otherwise it will expect that file to exist.
MOD=random_1mb_windows.mod

# The species you want to use
SPECIES=hg17,panTro1,rheMac2,mm7,rn3,canFam2,bosTau2

# Large phylogenetic tree from which we will select a subset
BIG_TREE=bigtree.nh

# Names of those species in the tree file, we will rename them
TREE_NAMES=human,chimp,macaque,mouse,rat,dog,cow

# The alignments training data should be extracted from, these are in maf 
# format named $(CHROMOSOME).maf and have been indexed with maf_build_index.py
ALIGN=/depot/data1/cache/human/hg17/align/17way/maf

# This file contains the length of each chromosome in the reference sequence
# and location of a nib file with the sequence. It is optional for everything
# except generating the 1mb windows for inferring substitution model.
CHROMLEN=chromlen.txt

# Intervals (in the reference species) of your positive and negative training
# datasets (the extension .bed is assumed)
POS=positive
NEG=negative

# This name will be used as the prefix for all generated files
BASENAME=example

# Maximal order of VOMM
ORDER=2

# Model spec (this default is to use a VOMM)
MODEL=tree_pruned:D=.01,N=10

# ---- Targets --------------------------------------------------------------

# The default target does everything up to the search
atoms : $(BASENAME).ent_agglom.75.mapping

# This target actually runs the search
search: $(BASENAME).ent_agglom.75.adapt
	
# Get the best mapping produced by the search
best: $(BASENAME).ent_agglom.75.adapt.best.mapping

# Train a scoring matrix for the best mapping
sm: $(BASENAME).ent_agglom.75.adapt.best.sm

# Removes all generated data files
clean :
	rm -rf *.maf *.col_counts *.adists *.clusters *.mapping *.adapt *.sm random_1mb_windows* tree.nh

# ---- Building substitution models -----------------------------------------

# Mapping from tree names to build names
RENAME=$(shell python -c "print '; '.join( [ '%s -> %s' % ( a, b ) for a, b in zip( '$(TREE_NAMES)'.split(','), '$(SPECIES)'.split(',') ) ] )")

# Extract subset of big tree and rename nodes
tree.nh : $(BIG_TREE)
	$(PH)/tree_doctor --prune-all-but "$(TREE_NAMES)" \
	                  --rename "$(RENAME)" \
					  --tree-only $< > $@

# Creae random sample of genomic windows
random_1mb_windows.bed : 
	$(RP)/genome_windows.py hg17 1000000 < $(CHROMLEN) | randomLines 50 > $@

# Extract MAF blocks for random windows
random_1mb_windows_split_maf: random_1mb_windows.bed
	mkdir -p $@
	rm -f $@/*.maf
	cat $< \
		| prefix_lines "hg17." \
		| maf_extract_ranges_indexed.py $(ALIGN)/chr*.maf \
		| $(BX)/maf_limit_to_species.py $(SPECIES) \
		| $(BX)/maf_split_by_src.py -c 0 -o $@/

random_1mb_windows.ss : random_1mb_windows_split_maf
	$(PH)/msa_view --aggregate $(SPECIES) -i MAF --unordered-ss -o SS $</*.maf > $@

# Estimate tree from sample
%.mod : %.ss tree.nh
	$(PH)/phyloFit --tree tree.nh --scale-only --gaps-as-bases --subst-mod HKY85+Gap -i SS $< --out-root $*					  

# ---- Building training data -----------------------------------------------
	
%.maf : %.bed
	cat $< \
	| $(BX)/prefix_lines "hg17." \
	| $(BX)/maf_tile_2.py --strand $(SPECIES) /dev/null $(foreach c, $(CHROMOSOMES), $(ALIGN)/chr$c.maf) \
	| $(BX)/maf_chop.py \
	| $(BX)/maf_translate_chars.py \
	| $(BX)/maf_filter_max_wc.py 50 3 \
	> $@

$(NEG).maf : $(NEG).bed $(POS).maf
	cat $< \
	| $(BX)/prefix_lines "hg17." \
	| $(BX)/maf_tile_2.py --strand $(SPECIES) /dev/null $(foreach c, $(CHROMOSOMES), $(ALIGN)/chr$c.maf) \
	| $(BX)/maf_chop.py \
	| $(BX)/maf_translate_chars.py \
	| $(BX)/maf_filter_max_wc.py 50 3 \
	| $(BX)/maf_randomize.py `$(BX)/maf_count.py < $(POS).maf` \
	> $@	

# For the first stages of inferrence we just merge the training sets
$(BASENAME).maf : $(POS).maf $(NEG).maf
	cat $^ > $@

# Create column counts for all columns with at most 3 missing species
%.col_counts : %.maf
	cat $^ | $(BX)/maf_col_counts_all.py -m 3 > $@

# Infer ancestral distributions for selected columns using substitution model
%.adists : %.col_counts $(MOD)
	$(RP)/ancestral_dists.py $(MOD) $(SPECIES) $< > $@

# Run 'entropy agglomeration' to group ancestral distributions, using an 
# initial k-nearest-neighbors agglomeration for columns that occur less than
# 10 times
%.clusters : %.adists
	$(RP)/entropy_agglomeration.py $< $@ 10

# Extract from the 'entropy agglomeration' a clustering for 75 clusters
%.ent_agglom.75.mapping : %.adists %.clusters
	$(RP)/make_mapping_from_clusters.py $^ 75 > $@

# Run the search starting from the 'entropy agglomeration' based clustering
%.adapt : %.mapping
	mkdir -p $@
	rm -f $@/*
	$(RP)/rp_adapt_mc_mpi.py $(POS).maf $(NEG).maf -d $@ -f maf -a $< -o $(ORDER) -M $(MODEL)
	
# Get the best mapping found by the search
%.adapt.best.mapping : %.adapt
	$(RP)/atom_mapping_to_complete.py $*.mapping `ls $</*.mapping | tail -1` > $@
	
# Train a scoring matrix using that mapping
%.sm : %.mapping $(POS).maf $(NEG).maf 
	$(RP)/rp_train.py $(POS).maf $(NEG).maf $@ -m $< -o $(ORDER) -f maf -M $(MODEL)
	
# Special rule to prevent 'make' from deleting output
.SECONDARY:
