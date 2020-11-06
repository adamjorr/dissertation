GATKVARIANTS = ../01-first-11-scaffolds.vcf.gz
CALLEDVARIANTS = ../08-variable.vcf.gz
JOBNUMS=$(shell seq 1 100)
VCFNUMS = $(addsuffix .vcf.gz, $(JOBNUMS))
JOBCALLS=$(addprefix %/,$(VCFNUMS))

SCRIPTDIR = ~/bin/zypy/scripts
DIPLOIDIFY = $(SCRIPTDIR)/diploidify.py
FPRFILT = $(SCRIPTDIR)/gatkfpr.py
SHUF_LABELS = $(SCRIPTDIR)/label_permutation.R
SHELL=bash
DNGCALLS = ../../dng/filter/goodsites.vcf
##########################################################

ALL = 01-first-11-scaffolds.vcf.gz 02-depth.vcf.gz 03-excesshet.vcf.gz \
        04-near-indel.vcf.gz 05-biallelic.vcf.gz 06-repeat.vcf.gz \
        07-replicate.vcf.gz 08-variable.vcf.gz

RESULTS = $(addsuffix .result, $(basename $(basename $(ALL))))
RESFOLDERS = $(addsuffix /, $(basename $(basename $(ALL))))
ALLVCFS = $(foreach folder, $(RESFOLDERS), $(addprefix $(folder), $(VCFNUMS)))
ALLTREES = $(addprefix trees/, $(addsuffix .ped,$(JOBNUMS)))

default: all

all: $(RESULTS) results.txt

.PHONY: all default

.DELETE_ON_ERROR:

.INTERMEDIATE:

.SECONDARY:

#########################################################
# Helper Rules

%.vcf.gz.csi: %.vcf.gz
	bcftools index $<

#

trees/:
	mkdir -p $@

trees/%.ped: trees/
	echo "##PEDNG v1.0" > $@; echo -e "M\t.\t.\t0" | paste - <( $(SHUF_LABELS) ) >> $@

$(RESFOLDERS):
	mkdir -p $@

vpath %.ped trees/

$(ALLVCFS): | $(RESFOLDERS)
	$(FPRFILT) -n 3 -v ../$(subst /,,$(dir $@)).vcf.gz -p trees/$(notdir $(basename $(basename $@))).ped | \
	$(DIPLOIDIFY) -t vcf -v | bgzip > $@

$(RESULTS): %.result : $(DNGCALLS) $(JOBCALLS)
	rm -f $@; for f in $(wordlist 2,$(words $^),$^) ; do bcftools view -H -T ^$< $$f | wc -l >> $@; done

#FP, FP/P (FDR)
results.txt: $(RESULTS)
	parallel -I{} --tag 'cat {} | awk "{x+=\$$1}END{print x,x/NR}"' ::: $^ > $@

