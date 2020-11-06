VARIANTSDIR = ../gatk4/
REF = ../e_mel_3/e_mel_3.fa
REFDICT = $(basename $(REF)).dict
RAWALIGNMENT = ../alignment.dedup.bam
KBBQALIGNMENT = recalibrated.bam
REPEATBED = ../liftover/e_mel_3_repeatmask.bed
RAWVARIANTS = var-calls.vcf.gz
SCAFFOLDS = $(addprefix scaffold_,$(shell seq 1 11))
SCAFFOLDEDVARIANTS = $(addsuffix .var-calls.vcf.gz, $(SCAFFOLDS))

SCRIPTDIR = ~/bin/zypy/scripts
FLTWITHREP = $(SCRIPTDIR)/filt_with_replicates.pl
DIPLOIDIFY = $(SCRIPTDIR)/diploidify.py
VCFTOFA = $(SCRIPTDIR)/vcf2fa.sh
PLOTTREE = $(SCRIPTDIR)/plot_tree.R

ALL = 01-first-11-scaffolds.vcf.gz 02-depth.vcf.gz 03-excesshet.vcf.gz \
        04-near-indel.vcf.gz 05-biallelic.vcf.gz 06-repeat.vcf.gz \
        07-replicate.vcf.gz 08-variable.vcf.gz

all: $(ALL) num_variants.txt

.SHELL: bash

.PHONY: all

%.vcf.gz: %.vcf
	bgzip -c $< >$@

%.vcf.gz.csi: %.vcf.gz
	bcftools index $<

%.bam.bai: %.bam
	samtools index $<

.SECONDARY:

.PRECIOUS:

$(KBBQALIGNMENT): $(RAWALIGNMENT)
	~/bin/cbbq/build/cbbq -t12 -k32 -g 605951511 -c 240 $< > $@ 2> cbbq.log

$(SCAFFOLDEDVARIANTS): %.var-calls.vcf.gz: $(KBBQALIGNMENT) $(REF) $(KBBQALIGNMENT).bai $(REFDICT)
	~/.local/bin/gatk HaplotypeCaller -I $< -R $(word 2, $^) --heterozygosity 0.025 -O $@ --QUIET -L $*

$(RAWVARIANTS): $(REFDICT) $(SCAFFOLDEDVARIANTS)
	picard SortVcf O=$@ SEQUENCE_DICTIONARY=$< $(addprefix I=,$(wordlist 2,$(words $^),$^))

01-first-11-scaffolds.vcf.gz: $(RAWVARIANTS) $(RAWVARIANTS).csi
	bcftools view -r scaffold_1,scaffold_2,scaffold_3,scaffold_4,scaffold_5,scaffold_6,scaffold_7,scaffold_8,scaffold_9,scaffold_10,scaffold_11 -Oz -o $@ $<

02-depth.vcf.gz: 01-first-11-scaffolds.vcf.gz 01-first-11-scaffolds.vcf.gz.csi
	bcftools filter -i 'INFO/DP <= 500' -O z -o $@ $<

03-excesshet.vcf.gz: 02-depth.vcf.gz 02-depth.vcf.gz.csi
	bcftools filter -i 'ExcessHet <= 40' -O z -o $@ $<

#testing shows -g is applied last
04-near-indel.vcf.gz: 03-excesshet.vcf.gz 03-excesshet.vcf.gz.csi
	bcftools filter -g 50 -O z -o $@ $<

05-biallelic.vcf.gz: 04-near-indel.vcf.gz 04-near-indel.vcf.gz.csi
	bcftools view -m2 -M2 -v snps -O z -o $@ $<

06-repeat.vcf.gz: 05-biallelic.vcf.gz $(REPEATBED) 05-biallelic.vcf.gz.csi
	vcftools --gzvcf $< --exclude-bed $(word 2, $^) --remove-filtered-all --recode --recode-INFO-all --stdout | bgzip -c > $@

07-replicate.vcf.gz: 06-repeat.vcf.gz 06-repeat.vcf.gz.csi
	zcat $< | $(FLTWITHREP) -s -g 3 | bgzip -c > $@

08-variable.vcf.gz: 07-replicate.vcf.gz
	zcat $< | $(DIPLOIDIFY) -t vcf -v | bgzip -c > $@

num_variants.txt: $(ALL)
	parallel --tag -I{} 'bcftools view -H {} | wc -l' ::: $(ALL) > $@


