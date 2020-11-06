
#VARIANTSDIR = ../../2018-07-24_callability/
#REF = ../../e_mel_3/e_mel_3.fa
#ALIGNMENT = $(VARIANTSDIR)/alignment.bam
#REPEATBED = ../../liftover/e_mel_3_repeatmask.bed
#RAWVARIANTS = $(VARIANTSDIR)/gatk/var-calls.vcf.gz
#INREPEATS = $(VARIANTSDIR)/gatk/inrepeats.txt
DNGCALLS = ../../dng/filter/goodsites.vcf

SCRIPTDIR = ~/bin/zypy/scripts
FLTWITHREP = $(SCRIPTDIR)/filt_with_replicates.pl
DIPLOIDIFY = $(SCRIPTDIR)/diploidify.py
VCFTOFA = $(SCRIPTDIR)/vcf2fa.sh
PLOTTREE = $(SCRIPTDIR)/plot_tree.R
CALLABILITYRES = $(SCRIPTDIR)/spiked_mutation_results.py

ALL = 01-first-11-scaffolds.vcf.gz 02-depth.vcf.gz \
	04-near-indel.vcf.gz 05-biallelic.vcf.gz 06-repeat.vcf.gz \
	07-replicate.vcf.gz 08-variable.vcf.gz

RESULTS = $(addsuffix .result, $(basename $(basename $(ALL))))
FNRESULTS = $(addsuffix .fn.result, $(basename $(basename $(ALL))))
TPRESULTS = $(addsuffix .tp.result, $(basename $(basename $(ALL))))

results: fn.txt tp.txt

.PHONY: results

.SECONDARY:

%.vcf.gz: %.vcf
	bgzip -c $< >$@

%.vcf.gz.csi: %.vcf.gz
	bcftools index $<

#Calls not in the callset but in DNGCALLS are FN
%.fn.result: ../%.vcf.gz
	bcftools view -H -T ^$< $(DNGCALLS) | wc -l > $@

#Calls that match DNGCALLS are TPs
%.tp.result: ../%.vcf.gz
	bcftools view -H -T $(DNGCALLS) $< | wc -l > $@

fn.txt: $(FNRESULTS)
	parallel -I{} --tag 'cat {}' ::: $^ > $@

tp.txt: $(TPRESULTS)
	parallel -I{} --tag 'cat {}' ::: $^ > $@
