AUTOSOMES=1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
CHROMOSOMES=$(AUTOSOMES) X Y
ALIGN=/depot/data1/cache/human/hg17/align/17way/maf

CHROMLEN=chromlen.txt

REF=hg17
SPECIES=hg17,panTro1,rheMac2,mm7,rn3,canFam2,bosTau2
REGIONS=$(CHROMOSOMES)

REORDER=0,1,2,3,4,5,6

MOD=example.ent_agglom.75.adapt.best
MODEL=tree_pruned:N=10,D=0.01

all_scores : $(foreach r, $(REGIONS), chr$r.scores.bz2)
all_truncated : $(foreach r, $(REGIONS), chr$r.scores.truncated.bz2)
all_covered : $(foreach r, $(REGIONS), chr$r.covered)

chr%.covered : $(ALIGN)/chr%.maf
	./maf_covered_ranges.py $(SPECIES) 4 < $< > $@

chr%.scores.bz2 : chr%.covered $(ALIGN)/chr%.maf 
	$(BX)/prefix_lines.py "$(REF)." < $< \
	| $(BX)/maf_tile_2.py $(SPECIES) chromlen.txt $(ALIGN)/chr$*.maf \
	| $(BX)/maf_translate_chars.py \
	| $(RP)/rp_score_maf.py /dev/stdin $(MOD).sm /dev/stdout -m $(MOD).mapping -w 100 -s 1 -M $(MODEL) \
	| bzip2 \
	> $@

chr%.scores.truncated.bz2 : chr%.scores.bz2
	bzcat $< | ../truncate.py | bzip2 > $@

.SECONDARY:
