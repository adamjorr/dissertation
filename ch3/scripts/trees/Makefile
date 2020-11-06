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

PDFS = $(addsuffix .pdf, $(basename $(basename $(ALL))))

pdfs: $(PDFS)

.PHONY: pdfs

%.vcf.gz: %.vcf
	bgzip -c $< >$@

%.vcf.gz.csi: %.vcf.gz
	bcftools index $<

%.fa: ../%.vcf.gz
	bcftools view -e 'GT="."' $< | $(VCFTOFA) | $(DIPLOIDIFY) -v >$@

.SECONDARY:

%.nwk: %.fa
	export RXMLTMPDIR=$(shell mktemp -d --tmpdir=$(CURDIR) raxml_XXXXXX); \
	raxml -T 8 -f a -s $< -n nwk -m ASC_GTRGAMMA -w $$RXMLTMPDIR --asc-corr=lewis -p 12345 -x 12345 -# 100 && \
	cat $$RXMLTMPDIR/RAxML_bestTree.nwk >$@ &&\
	rm -rf $$RXMLTMPDIR

%.pdf: %.nwk
	$(PLOTTREE) $< $@

%.png: %.pdf
	gs -sDEVICE=pngalpha -o $@ $<
