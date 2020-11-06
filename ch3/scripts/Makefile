VARIANTSDIR = ../gatk4/
REF = ../e_mel_3/e_mel_3.fa
REFDICT = $(basename $(REF)).dict
RAWALIGNMENT = ../alignment.dedup.bam
REPEATBED = ../liftover/e_mel_3_repeatmask.bed
RAWVARIANTS = discoRes_k_31_c_3_D_100_P_3_b_0_coherent.vcf
SCAFFOLDS = $(addprefix scaffold_,$(shell seq 1 11))
SCAFFOLDEDVARIANTS = $(addsuffix .var-calls.vcf.gz, $(SCAFFOLDS))
DATADIR = ../data/

SCRIPTDIR = ~/bin/zypy/scripts
FLTWITHREP = $(SCRIPTDIR)/filt_with_replicates.pl
DIPLOIDIFY = $(SCRIPTDIR)/diploidify.py
VCFTOFA = $(SCRIPTDIR)/vcf2fa.sh
PLOTTREE = $(SCRIPTDIR)/plot_tree.R

SAMPLENUMS := 1 2 3 4 5 6 7 8
SAMPLEREPS := a b c
SAMPLES := $(foreach NUM,$(SAMPLENUMS),$(addprefix M$(NUM),$(SAMPLEREPS)))
SAMPLECMD := $(foreach SAM,$(SAMPLES),echo $(wildcard $(DATADIR)/*$(SAM)*/$(SAM)*.fastq | tr ' ' '\t' > $(SAM)_fof.txt))
SAMPLEFOFS := $(addsuffix _fof.txt,$(SAMPLES))

ALL = 01-first-11-scaffolds.vcf.gz 02-depth.vcf.gz \
        04-near-indel.vcf.gz 05-biallelic.vcf.gz 06-repeat.vcf.gz \
        07-replicate.vcf.gz 08-variable.vcf.gz

all: $(ALL) num_variants.txt

.SHELL: bash

.PHONY: all disco

%.vcf.gz: %.vcf
	bgzip -c $< >$@

%.vcf.gz.csi: %.vcf.gz
	bcftools index $<

%.vcf.gz.tbi: %.vcf.gz
	bcftools index -t $<

%.bam.bai: %.bam
	samtools index $<

.SECONDARY:

.PRECIOUS:

%_fof.txt:
	echo $(wildcard $(DATADIR)/*$**/$**.fastq) $% |\
	xargs -n1 > $@

fof.txt: $(SAMPLEFOFS)
	echo $^ | xargs -n1 > $@

disco: fof.txt $(REF)
	~/bin/DiscoSnp/run_discoSnp++.sh -r $< -T -G $(word 2, $^) -R

disco_sample_names.txt: discoRes_read_files_correspondance.txt
	sed 's/C_/G/; s/[[:space:]].*\(M[[:digit:]][abc]\).*/ \1/g' $< >$@

#sorted.vcf.gz: $(RAWVARIANTS).gz disco_sample_names.txt
#	bcftools reheader --samples $(word 2, $^) $< | \
#	bcftools sort -m64G -Oz -o $@ -

01-first-11-scaffolds.vcf.gz: $(RAWVARIANTS).gz disco_sample_names.txt $(REF).fai
	bcftools reheader --samples $(word 2, $^) -f $(word 3, $^) $< | \
	bcftools view -t scaffold_1,scaffold_2,scaffold_3,scaffold_4,scaffold_5,scaffold_6,scaffold_7,scaffold_8,scaffold_9,scaffold_10,scaffold_11 -Ou | \
	bcftools sort -m 64G -Oz -o $@ -

02-depth.vcf.gz: 01-first-11-scaffolds.vcf.gz 01-first-11-scaffolds.vcf.gz.csi
	bcftools filter -i 'SUM(FORMAT/DP) <= 500' -O z -o $@ $<

#03-excesshet.vcf.gz: 02-depth.vcf.gz 02-depth.vcf.gz.csi
#	bcftools filter -i 'ExcessHet <= 40' -O z -o $@ $<

#testing shows -g is applied last
04-near-indel.vcf.gz: 02-depth.vcf.gz 02-depth.vcf.gz.csi
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
	parallel --tag -I{} "bcftools view -e 'GT=\".\"' -H {} | wc -l" ::: $(ALL) > $@


