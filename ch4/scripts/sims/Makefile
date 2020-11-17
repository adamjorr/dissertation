SHELL=bash
EGRANDISURL1 = https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/AU/SX/AUSX01/AUSX01.1.fsa_nt.gz
EGRANDISURL2 = https://sra-download.ncbi.nlm.nih.gov/traces/wgs03/wgs_aux/AU/SX/AUSX01/AUSX01.2.fsa_nt.gz
EGRANDISPT1 = $(notdir $(EGRANDISURL1))
EGRANDISPT2 = $(notdir $(EGRANDISURL2))
SAMPLEDREFSIZE = 5000000
ARTDIR = ./art_bin_MountRainier/

.PRECIOUS:

.SECONDARY:

#index helpers
%.fa.gz.fai: %.fa.gz
	samtools faidx $<

%.bam.bai: %.bam
	samtools index $<

%.vcf.gz: %.vcf
	bgzip -c $< > $@

%.vcf.gz.tbi: %.vcf.gz
	bcftools index -t $<

%.dict: %.fa
	gatk CreateSequenceDictionary -R $<

%.dict: %.fa.gz
	gatk CreateSequenceDictionary -R $<

#download AUSX00000000.1, the e. grandis 1.0 contigs.
#this would represent an assembly where scaffolds haven't been joined into contigs
$(EGRANDISPT1):
	wget $(EGRANDISURL1) -O $@

$(EGRANDISPT2):
	wget $(EGRANDISURL2) -O $@

full_grandis_contigs.fa.gz: $(EGRANDISPT1) $(EGRANDISPT2)
	cat <(zcat $<) <(zcat $(word 2, $^)) | bgzip -c > $@

#subsample contigs so we have a reasonable amount of data
#GATK recommends at least 100M bp per read group
#to have 100M bp with 150bp reads, we need 100000000/150 = 666,667 reads
#at 20X coverage, which is high but reasonable for a non-model organism, this requires
# 20 = sequencedbp / gsize; gsize = sequencedbp/20
#So we need a genome size of 5M with 20x to read 100M bases.
sampled_contigs.txt: full_grandis_contigs.fa.gz.fai
	cat $< | shuf | awk '{if (x < $(SAMPLEDREFSIZE)) {x += $$2; print $$1} else exit}' > $@

sampled_contigs.fa.gz: sampled_contigs.txt full_grandis_contigs.fa.gz
	cat $< | xargs samtools faidx $(word 2, $^) | bgzip -c > $@

#theta 0.025
#5,000,000 bp * 0.025 = 125,000
#titv ratio of 2 is an appropriate estimate for eucs
simug.simseq.genome.fa simug.refseq2simseq.SNP.vcf &: sampled_contigs.fa.gz
	simuG.pl -seed 20201112 -titv_ratio 2.0 -snp_count $$(( $(SAMPLEDREFSIZE) * 25 / 1000 )) -refseq $< -prefix simug

#concatenate both references to generate a file with both
#haplotypes in it
#adds a .s suffix to the sampled contigs so they can be differentiated
combined_ref.fa: sampled_contigs.fa.gz simug.simseq.genome.fa
	zcat $< | sed '/^>/s/$$/.s/' | cat - $(word 2, $^) > $@

#simulate reads from the reference with ART
art1.fq art2.fq art.sam &: combined_ref.fa
	$(ARTDIR)/art_illumina -sam -na -ss HS25 -p -l 101 -m 300 -s 20 -rs 20201112 -i $< -f 20 -o art

#align with ngm
ngm.bam: sampled_contigs.fa.gz art1.fq art2.fq
	ngm -r $< -1 $(word 2, $^) -2 $(word 3, $^) -t 4 --rg-id sim --rg-sm sim --rg-lb sim --rg-pl ILLUMINA --rg-pu sim | samtools sort -m 2G -O bam -o $@ -

#recalibrate with true variants
%.recal.txt: %.bam sampled_contigs.fa.gz simug.refseq2simseq.SNP.vcf.gz %.bam.bai sampled_contigs.dict simug.refseq2simseq.SNP.vcf.gz.tbi
	gatk BaseRecalibrator -I $< -R $(word 2, $^) --known-sites $(word 3, $^) -O $@

#recalibrate with true variants
%.recal.bam: %.bam sampled_contigs.fa.gz %.recal.txt
	gatk ApplyBQSR -I $< -R $(word 2, $^) --bqsr-recal-file $(word 3, $^) -O $@

#benchmark recalibration with true variants
%.recal.after.txt: %.recal.bam sampled_contigs.fa.gz simug.refseq2simseq.SNP.vcf.gz %.recal.bam.bai sampled_contigs.dict simug.refseq2simseq.SNP.vcf.gz.tbi
	gatk BaseRecalibrator -I $< -R $(word 2, $^) --known-sites $(word 3, $^) -O $@

%.recal.after.csv: %.recal.txt %.recal.after.txt
	gatk AnalyzeCovariates -before $< -after $(word 2, $^) -csv $@

%.recal.after.pdf: %.recal.txt %.recal.after.txt
	gatk AnalyzeCovariates -before $< -after $(word 2, $^) -plots $@

kbbq-ngm.recal.bam: ngm.bam
	~/bin/cbbq/build/cbbq -t6 -k32 -g $(SAMPLEDREFSIZE) $< > $@ 2> cbbqlog.txt

kbbq-ngm.recal.after.csv: kbbq-ngm.recal.after.txt
	gatk AnalyzeCovariates -after $< -csv $@

#get initial calls with HC and use them to recalibrate
initial-calls.vcf.gz initial-calls.vcf.gz.tbi &: ngm.bam sampled_contigs.fa.gz ngm.bam.bai sampled_contigs.dict
	gatk HaplotypeCaller -R $(word 2, $^) -I $< -stand-call-conf 50 --native-pair-hmm-threads 8 -O $@ --heterozygosity 0.025

initial-calls.recal.txt: ngm.bam sampled_contigs.fa.gz initial-calls.vcf.gz ngm.bam.bai sampled_contigs.dict initial-calls.vcf.gz.tbi
	gatk BaseRecalibrator -I $< -R $(word 2, $^) --known-sites $(word 3, $^) -O $@

initial-calls.recal.bam: ngm.bam sampled_contigs.fa.gz initial-calls.recal.txt
	gatk ApplyBQSR -I $< -R $(word 2, $^) --bqsr-recal-file $(word 3, $^) -O $@
