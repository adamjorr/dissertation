GATKVARIANTS = ../gatk4/repeat-filter.vcf
SCRIPTDIR = ~/bin/zypy/scripts
FLTWITHREP = $(SCRIPTDIR)/filt_with_replicates.pl
DIPLOIDIFY = $(SCRIPTDIR)/diploidify.py
VCFTOFA = $(SCRIPTDIR)/vcf2fa.sh
PLOTTREE = $(SCRIPTDIR)/plot_tree.R
GATKFPR = $(SCRIPTDIR)/gatkfpr.py

ALL = 01-first-11-scaffolds.vcf.gz 02-depth.vcf.gz \
        04-near-indel.vcf.gz 05-biallelic.vcf.gz 06-repeat.vcf.gz \
        07-replicate.vcf.gz 08-variable.vcf.gz

TXTS = $(addsuffix .txt, $(basename $(basename $(ALL))))

.DEFAULT: n_incomplete.txt

.PHONY: pdfs

%.vcf.gz: %.vcf
	bgzip -c $< >$@

%.vcf.gz.csi: %.vcf.gz
	bcftools index $<

%.txt: ../%.vcf.gz
	bcftools view -e 'GT="."' -Ou $< | bcftools view -H -i'GT="./."' | wc -l > $@

.SECONDARY:

n_sites_missing.txt: $(TXTS)
	parallel --tag -I{} 'cat {}' ::: $^ > $@
