DATADIR = ../../../../
VCFIN = $(DATADIR)/chr1.vcf.gz
BEDIN = $(DATADIR)/chr1_confident.bed.gz
REFIN = $(DATADIR)/chr1.renamed.fa
HAPDIPJS = $(HOME)/bin/hapdip/hapdip.js
TAG = QUAL

INPUTINCREMENT = 20
INPUTFIRST = 0
INPUTLAST = 100

SHELL=bash
FPRS := $(shell seq $(INPUTFIRST) $(INPUTINCREMENT) $(INPUTLAST))
#FNR_FPR
ADD_FPRS = $(addsuffix _$(FPR),$(FPRS))
ERRATES := $(foreach FPR,$(FPRS),$(ADD_FPRS))
ERRATES := $(ERRATES) raw cbbq
RTGROCS := $(addsuffix _$(TAG)/snp_roc.tsv.gz,$(ERRATES))
OUTROCS := $(addsuffix _$(TAG).snp_roc.tsv.gz,$(ERRATES))

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

all: $(OUTROCS)

$(RTGROCS): %_$(TAG)/snp_roc.tsv.gz: ../%.vcf.gz $(VCFIN) $(BEDIN) $(REFIN).sdf
	rm -rf $(dir $@) && rtg vcfeval -b $(word 2, $^) -c $< -e $(word 3, $^) -t $(word 4, $^) --output-mode=roc-only --vcf-score-field=$(TAG) -T 4 -o $(dir $@)

$(OUTROCS): %_$(TAG).snp_roc.tsv.gz: %_$(TAG)/snp_roc.tsv.gz
	cp $< $@

