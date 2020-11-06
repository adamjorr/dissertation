DATADIR = ../../

VCFIN = $(DATADIR)/chr1.vcf.gz
BAMIN = $(DATADIR)/chr1.bam
BEDIN = $(DATADIR)/chr1_confident.bed.gz
REFIN = $(DATADIR)/chr1.renamed.fa
PROCESSEDBAM = chr1_only_confident.fixmate.sorted.bam

VCFLINES := $(shell bcftools index -n $(VCFIN))
INPUTINCREMENT = 20
INPUTFIRST = 0
INPUTLAST = 100

SHELL=bash
FPRS := $(shell seq $(INPUTFIRST) $(INPUTINCREMENT) $(INPUTLAST))
#FNR_FPR
ADD_FPRS = $(addsuffix _$(FPR),$(FPRS))
ERRATES := $(foreach FPR,$(FPRS),$(ADD_FPRS))
ERBENCHMARKS := $(addsuffix .benchmark,$(ERRATES))
ERPLOTS := $(addsuffix .recal.after.pdf,$(ERRATES))
ERCSVS := $(addsuffix .recal.after.csv,$(ERRATES))

#functions to get FNR and FPR of a given file named FNR_FPR.*
getfnr = $(word 1, $(subst _, ,$(1)))
getfpr = $(word 2, $(subst ., ,$(subst _, ,$(1))))

#actual negative and actual positive files
ANBED = ansites.bed.gz
APBED = apsites.bed.gz
ANLINES = $(shell zcat $(ANBED) | wc -l)
APLINES = $(shell zcat $(APBED) | wc -l)

.PHONY: plots csvs all

plots: $(ERPLOTS)

csvs: $(ERCSVS)

all: plots csvs

.DEFAULT: all

.PRECIOUS:

.SECONDARY:

#index helpers
%.bed.gz.tbi: %.bed.gz
	tabix -p bed $<

%.vcf.gz.tbi: %.vcf.gz
	bcftools index -t $<

%.bam.bai: %.bam
	samtools index $<

#process raw BAM
chr1_only_confident.namesorted.bam: $(BAMIN) $(BEDIN)
	samtools view -h $< -F 3844 -L $(word 2, $^) | samtools sort -@8 -n -o $@ -O bam -

chr1_only_confident.fixmate.bam: chr1_only_confident.namesorted.bam
	samtools fixmate -r --threads 8 $< $@

chr1_only_confident.fixmate.onlypaired.bam: chr1_only_confident.fixmate.bam
	samtools view -h -b -o $@ -f 1 $<

$(PROCESSEDBAM): chr1_only_confident.fixmate.onlypaired.bam
	samtools sort -@ 8 -O bam -o $@ $<

#--- BED files containing the appropraite FNR_FPR ---
#false negative rate = FN / AP
$(APBED): $(VCFIN) $(BEDIN)
	bcftools query -R $(word 2, $^) -f '%CHROM\t%POS0\t%END\n' $< | \
	sort -k1,1 -k2,2n | bedtools merge -i - | \
	bedtools intersect -sorted -a - -b $(word 2, $^) | \
	bedtools makewindows -b - -w 1 | \
	sort -k1,1 -k2,2n | bgzip >$@

%.fn.bed.gz: $(APBED)
	zcat $< | \
	shuf -n $(shell echo $$(($(APLINES) * (100 - $(call getfnr,$*))/100))) | \
	sort -k1,1 -k2,2n | bgzip > $@

#false positive rate =  FP / AN
$(ANBED): $(BEDIN) $(APBED)
	bedtools subtract -a $< -b $(word 2, $^) | \
	bedtools makewindows -b - -w 1 | \
	sort -k1,1 -k2,2n | bgzip > $@

%.fp.bed.gz: $(ANBED)
	zcat $< | shuf -n $(shell echo $$(($(ANLINES) * $(call getfpr,$*) / 100))) | \
	sort -k1,1 -k2,2n | bgzip > $@

#combine the false positives with the appropriate number
#of false negatives
%.toskip.bed.gz: %.fp.bed.gz %.fn.bed.gz
	zcat $(word 2, $^) | cat <(zcat $<) - | \
	sort -k1,1 -k2,2n | bedtools merge -i - | bgzip >$@

# --- Do calibration ---
%.recal.txt: $(PROCESSEDBAM) $(REFIN) $(BEDIN) %.toskip.bed.gz $(PROCESSEDBAM).bai %.toskip.bed.gz.tbi
	gatk BaseRecalibrator -I $< -R $(word 2, $^) -L $(word 3, $^) --known-sites $(word 4, $^) --use-original-qualities -O $@

%.recal.bam: $(PROCESSEDBAM) $(REFIN) %.recal.txt
	gatk ApplyBQSR -I $< -R $(word 2, $^) --bqsr-recal-file $(word 3, $^) --use-original-qualities -O $@

cbbq.recal.bam: $(PROCESSEDBAM)
	~/bin/cbbq/build/cbbq -t6 -k32 -g 214206308 -a .15 --use-oq $< > $@ 2> cbbqlog.txt


# --- Benchmark with GATK ---
# Use actual positives for this test.
%.recal.after.txt: %.recal.bam $(REFIN) $(BEDIN) $(APBED) $(APBED).tbi
	gatk BaseRecalibrator -I $< -R $(word 2, $^) -L $(word 3, $^) --known-sites $(word 4, $^) -O $@

%.recal.after.pdf: %.recal.txt %.recal.after.txt
	gatk AnalyzeCovariates -before $< -after $(word 2, $^) -plots $@

raw.recal.after.pdf: raw.recal.after.txt
	gatk AnalyzeCovariates -bqsr $< -plots $@

%.recal.after.csv: %.recal.txt %.recal.after.txt
	gatk AnalyzeCovariates -before $< -after $(word 2, $^) -csv $@

raw.recal.after.csv: raw.recal.after.txt
	gatk AnalyzeCovariates -bqsr $< -csv $@

allplots.pdf: $(ERPLOTS)
	gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$@ $^
