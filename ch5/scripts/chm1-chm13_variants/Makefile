DATADIR = ../../../
VCFIN = $(DATADIR)/chr1.vcf.gz
BEDIN = $(DATADIR)/chr1_confident.bed.gz
REFIN = $(DATADIR)/chr1.renamed.fa
RAWBAM =  ../chr1_only_confident.fixmate.sorted.bam
PROCESSEDBAM = raw.bam
HAPDIPJS = $(HOME)/bin/hapdip/hapdip.js

INPUTINCREMENT = 20
INPUTFIRST = 0
INPUTLAST = 100

SHELL=bash
FPRS := $(shell seq $(INPUTFIRST) $(INPUTINCREMENT) $(INPUTLAST))
#FNR_FPR
ADD_FPRS = $(addsuffix _$(FPR),$(FPRS))
ERRATES := $(foreach FPR,$(FPRS),$(ADD_FPRS))
ERRATES := $(ERRATES) cbbq
ERBAMS := $(addprefix ../,$(addsuffix .recal.bam,$(ERRATES)))
ERVCFS := $(addsuffix .vcf.gz,$(ERRATES))
RTGCALLS := $(addsuffix /calls.vcf.gz,$(ERRATES))
SNPSUMS := $(addsuffix /hapdip_summary.txt,$(ERRATES))
OUTSUMS := $(addsuffix .summary.txt,$(ERRATES))

#index helpers
%.bed.gz.tbi: %.bed.gz
	tabix -p bed $<

%.vcf.gz.tbi: %.vcf.gz
	bcftools index -t $<

%.bam.bai: %.bam
	samtools index $<

%.fa.sdf: %.fa
	rtg format -o $@ $<

.PHONY: all

.DEFAULT: all

.PRECIOUS:

.SECONDARY:

all: $(OUTSUMS)

$(ERVCFS): %.vcf.gz: ../%.recal.bam $(REFIN) $(BEDIN)
	gatk HaplotypeCaller -R $(word 2, $^) -I $< -L $(word 3, $^) --native-pair-hmm-threads 16 -O $@ --QUIET true

$(PROCESSEDBAM): $(RAWBAM)
	gatk RevertBaseQualityScores -I $< -O $@

raw.vcf.gz: $(PROCESSEDBAM) $(REFIN) $(BEDIN)
	gatk HaplotypeCaller -R $(word 2, $^) -I $< -L $(word 3, $^) --native-pair-hmm-threads 16 -O $@ --QUIET true

#benchmarking result directories
$(RTGCALLS): %/calls.vcf.gz: %.vcf.gz $(VCFIN) $(BEDIN) $(REFIN).sdf
	rm -rf $(dir $@) && rtg vcfeval -b $(word 2, $^) -c $< -e $(word 3, $^) -o $(dir $@) -t $(word 4, $^) --output-mode=annotate -T 1

raw/calls.vcf.gz: raw.vcf.gz $(VCFIN) $(BEDIN) $(REFIN).sdf
	rm -rf $(dir $@) && rtg vcfeval -b $(word 2, $^) -c $< -e $(word 3, $^) -o $(dir $@) -t $(word 4, $^) --output-mode=annotate -T 8

$(SNPSUMS): %/hapdip_summary.txt: %/calls.vcf.gz
	k8 $(HAPDIPJS) rtgeval $(BEDIN) $(dir $@) > $@

raw/hapdip_summary.txt: raw/calls.vcf.gz
	k8 $(HAPDIPJS) rtgeval $(BEDIN) $(dir $@) > $@

$(OUTSUMS): %.summary.txt: %/hapdip_summary.txt
	cp $< $@

raw.summary.txt: raw/hapdip_summary.txt
	cp $< $@
