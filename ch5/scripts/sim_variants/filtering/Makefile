DATADIR = ../../../../
VCFIN = $(DATADIR)/chr1.vcf.gz
BEDIN = $(DATADIR)/chr1_confident.bed.gz
REFIN = $(DATADIR)/chr1.renamed.fa
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
ERVCFS := $(addprefix ../,$(addsuffix .vcf.gz,$(ERRATES)))
FLTVCFS := $(addsuffix _flt.vcf.gz,$(ERRATES))
RTGCALLS := $(addsuffix _flt/calls.vcf.gz,$(ERRATES))
SNPSUMS := $(addsuffix _flt/hapdip_summary.txt,$(ERRATES))
OUTSUMS := $(addsuffix _flt.summary.txt,$(ERRATES))

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

#filter: only snps
#filter: DP > 35 & DP < 65 (+-2 SDEV from mean of ~50)
#filter: QUAL > 75; this is ~bottom 1% and the maximum F
#   value from the previous analysis is between 75 and 150.
$(FLTVCFS): %_flt.vcf.gz: ../%.vcf.gz $(BEDIN)
	bcftools view -v snps -i'DP>35 & DP < 65 & QUAL > 75' -R $(word 2, $^) -Oz -o $@ $<

raw_flt.vcf.gz: ../raw.vcf.gz $(BEDIN)
	bcftools view -v snps -i'DP>35 & DP < 65 & QUAL > 75' -R $(word 2, $^) -Oz -o $@ $<

#benchmarking result directories
$(RTGCALLS): %/calls.vcf.gz: %.vcf.gz $(VCFIN) $(BEDIN) $(REFIN).sdf %.vcf.gz.tbi
	rm -rf $(dir $@) && rtg vcfeval -b $(word 2, $^) -c $< -e $(word 3, $^) -o $(dir $@) -t $(word 4, $^) --output-mode=annotate -T 4

raw_flt/calls.vcf.gz: raw_flt.vcf.gz $(VCFIN) $(BEDIN) $(REFIN).sdf raw_flt.vcf.gz.tbi
	rm -rf $(dir $@) && rtg vcfeval -b $(word 2, $^) -c $< -e $(word 3, $^) -o $(dir $@) -t $(word 4, $^) --output-mode=annotate -T 4

$(SNPSUMS): %/hapdip_summary.txt: %/calls.vcf.gz
	k8 $(HAPDIPJS) rtgeval $(BEDIN) $(dir $@) > $@

raw_flt/hapdip_summary.txt: raw_flt/calls.vcf.gz
	k8 $(HAPDIPJS) rtgeval $(BEDIN) $(dir $@) > $@

$(OUTSUMS): %.summary.txt: %/hapdip_summary.txt
	cp $< $@

raw_flt.summary.txt: raw_flt/hapdip_summary.txt
	cp $< $@
