DATADIR = ../../

VCFIN = $(DATADIR)/chr1.vcf.gz
BAMIN = ../2020-09-13/chr1_only_confident.fixmate.sorted.bam
BEDIN = $(DATADIR)/chr1_confident.bed.gz
REFIN = $(DATADIR)/chr1.renamed.fa

#this bam has the quality scores reverted since
#NGM can't handle OQ tags
PROCESSEDBAM = ../2020-09-13/variants/raw.bam
CHIMPBAM = chimpaln.bam
CHIMPSKIP = chimpaln.vcf.gz
CHIMPREF = GCF_002880755.1_Clint_PTRv2_genomic.fa

SHELL=bash

.PHONY: all

all: chimpaln.recal.after.pdf chimpaln.recal.after.csv chimpaln_30.recal.after.pdf chimpaln_30.recal.after.csv

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

%.fa.fai: %.fa
	samtools faidx $<

%.dict: %.fa
	picard CreateSequenceDictionary R=$< O=$@

# --- Realign reads to chimp ref and get skips ---
collated_reads.bam: $(PROCESSEDBAM)
	samtools collate -f -u -@4 -o $@ $<

rglines.txt: collated_reads.bam
	samtools view -H collated_reads.bam | grep ^@RG > $@

rgargs.txt: rglines.txt
	cat $< | sed 's/^@RG\t// ; y/\t/ / ; s/\([[:upper:]]\{2\}\):/--rg-\L\1=/g' > $@

rgids.txt: rglines.txt
	cat $< | sed 's/^@RG\t// ; s/ID:\([^[:space:]]*\)\t.*/\1/g' > $@

$(CHIMPBAM): $(CHIMPREF) rgids.txt rgargs.txt collated_reads.bam
	samtools split -f "%!.%." -@8 $(word 4, $^) && \
	parallel --link 'ngm --no-progress -r $< -p -q {1}.bam {2} | samtools sort -m 2G -O bam -o {1}.chimpaln.bam - ' :::: $(word 2, $^) :::: $(word 3, $^) && \
	cat $(word 2, $^) | sed 's/$$/.chimpaln.bam/' | xargs samtools merge -@8 -c -p $@

$(CHIMPSKIP): $(CHIMPBAM) $(CHIMPREF) $(CHIMPBAM).bai $(patsubst %.fa,%.dict,$(CHIMPREF)) $(CHIMPREF).fai
	gatk HaplotypeCaller -R $(word 2, $^) -I $< -stand-call-conf 50 --native-pair-hmm-threads 8 -O $@

#use -stand-call-conf 30 instead of 50 or 40 (default)
chimpaln_30.vcf.gz: $(CHIMPBAM) $(CHIMPREF) $(CHIMPBAM).bai $(patsubst %.fa,%.dict,$(CHIMPREF)) $(CHIMPREF).fai
	gatk HaplotypeCaller -R $(word 2, $^) -I $< -stand-call-conf 30 --native-pair-hmm-threads 8 -O $@

# --- Do calibration ---
chimpaln.recal.txt: $(CHIMPBAM) $(CHIMPREF) $(CHIMPSKIP) $(CHIMPBAM).bai
	gatk BaseRecalibrator -I $< -R $(word 2, $^) --known-sites $(word 3, $^) -O $@

chimpaln_30.recal.txt: $(CHIMPBAM) $(CHIMPREF) chimpaln_30.vcf.gz $(CHIMPBAM).bai
	gatk BaseRecalibrator -I $< -R $(word 2, $^) --known-sites $(word 3, $^) -O $@

# Recalibrate the alignment to human using
# the model from the Chimp alignment
%.recal.bam: $(BAMIN) $(REFIN) %.recal.txt
	gatk ApplyBQSR -I $< -R $(word 2, $^) --bqsr-recal-file $(word 3, $^) --use-original-qualities -O $@

# --- Benchmark with GATK ---
# Use actual positives for this test.
%.recal.after.txt: %.recal.bam $(REFIN) $(BEDIN) $(VCFIN) $(VCFIN).tbi
	gatk BaseRecalibrator -I $< -R $(word 2, $^) -L $(word 3, $^) --known-sites $(word 4, $^) -O $@

%.recal.after.pdf: %.recal.txt %.recal.after.txt
	gatk AnalyzeCovariates -before $< -after $(word 2, $^) -plots $@

%.recal.after.csv: %.recal.txt %.recal.after.txt
	gatk AnalyzeCovariates -before $< -after $(word 2, $^) -csv $@
