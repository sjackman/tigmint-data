# Correct misassemblies using Tigmint

# Number of threads
t=16

# Options for biomake
QsubArgs=-c$t --mem-per-cpu=7900M

# gzip compression program
gzip=pigz -p$t

# Parameters of ARCS
c=5
e=30000
r=0.05

# Parameters of ABySS-Scaffold
s=5000
n=20

# Parameters of LINKS
a=0.1
l=10
z=500

# Report run time and memory usage
time=command time -v -o $@.time
export SHELL=zsh -e -o pipefail
export REPORTTIME=1
export TIMEFMT=time user=%U system=%S elapsed=%E cpu=%P memory=%M job=%J

.DELETE_ON_ERROR:
.SECONDARY:

all: reports tables

abyss2_bam:
	$(MAKE) draft=$@ \
		abyss2.abyss-fac.tsv \
		abyss2.scaftigs.abyss-fac.tsv \
		abyss2.scaftigs.GRCh38.samtobreak.tsv \
		abyss2.hg004.c$c_e$e_r$r.arcs.a$a_l$l.links.abyss-fac.tsv \
		abyss2.hg004.c$c_e$e_r$r.arcs.a$a_l$l.links.scaftigs.abyss-fac.tsv \
		abyss2.hg004.c$c_e$e_r$r.arcs.a$a_l$l.links.scaftigs.GRCh38.samtobreak.tsv \
		abyss2.hg004.bx.as100.nm5.bam.bai \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.tsv \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.bed.bam.bai \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.summary.html \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.genomecov.stats.tsv \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.fa \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.abyss-fac.tsv \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.scaftigs.abyss-fac.tsv \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.scaftigs.GRCh38.samtobreak.tsv \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.hg004.c$c_e$e_r$r.arcs.a$a_l$l.links.abyss-fac.tsv \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.hg004.c$c_e$e_r$r.arcs.a$a_l$l.links.scaftigs.abyss-fac.tsv \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.hg004.c$c_e$e_r$r.arcs.a$a_l$l.links.scaftigs.GRCh38.samtobreak.tsv

# Aggregate the results.
abyss2_aggregate:
	$(MAKE) draft=abyss2 \
		abyss2.depth.100.starts.1-5.samtobreak.gscore.tsv \
		abyss2.depth.80-120.starts.2-4.arcs.abyss-fac.tsv \
		abyss2.depth.80-120.starts.2-4.arcs.samtobreak.tsv

abyss2_bionano_arcs:
	$(MAKE) draft=$@ \
		abyss2_bionano_arcs.abyss-fac.tsv \
		abyss2_bionano_arcs.scaftigs.abyss-fac.tsv \
		abyss2_bionano_arcs.scaftigs.GRCh38.samtobreak.tsv \
		abyss2_bionano_arcs.hg004.bx.as100.nm5.bam.bai \
		abyss2_bionano_arcs.hg004.bx.as100.nm5.bam.mi.bx.molecule.tsv \
		abyss2_bionano_arcs.hg004.bx.as100.nm5.bam.mi.bx.molecule.bed.bam.bai \
		abyss2_bionano_arcs.hg004.bx.as100.nm5.bam.mi.bx.molecule.summary.html \
		abyss2_bionano_arcs.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.genomecov.stats.tsv \
		abyss2_bionano_arcs.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.fa \
		abyss2_bionano_arcs.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.abyss-fac.tsv \
		abyss2_bionano_arcs.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.scaftigs.abyss-fac.tsv \
		abyss2_bionano_arcs.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.scaftigs.GRCh38.samtobreak.tsv \
		abyss2_bionano_arcs.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.hg004.c$c_e$e_r$r.arcs.a$a_l$l.links.abyss-fac.tsv \
		abyss2_bionano_arcs.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.hg004.c$c_e$e_r$r.arcs.a$a_l$l.links.scaftigs.abyss-fac.tsv \
		abyss2_bionano_arcs.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.hg004.c$c_e$e_r$r.arcs.a$a_l$l.links.scaftigs.GRCh38.samtobreak.tsv

supernova_bam:
	$(MAKE) draft=$@ \
		supernova.abyss-fac.tsv \
		supernova.scaftigs.abyss-fac.tsv \
		supernova.scaftigs.GRCh38.samtobreak.tsv \
		supernova.hg004.bx.as100.nm5.bam.bai \
		supernova.hg004.bx.as100.nm5.bam.mi.bx.molecule.tsv \
		supernova.hg004.bx.as100.nm5.bam.mi.bx.molecule.bed.bam.bai \
		supernova.hg004.bx.as100.nm5.bam.mi.bx.molecule.summary.html \
		supernova.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.genomecov.stats.tsv \
		supernova.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.fa \
		supernova.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.abyss-fac.tsv \
		supernova.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.scaftigs.abyss-fac.tsv \
		supernova.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.scaftigs.GRCh38.samtobreak.tsv \
		supernova.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.hg004.c$c_e$e_r$r.arcs.a$a_l$l.links.abyss-fac.tsv \
		supernova.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.hg004.c$c_e$e_r$r.arcs.a$a_l$l.links.scaftigs.abyss-fac.tsv \
		supernova.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.hg004.c$c_e$e_r$r.arcs.a$a_l$l.links.scaftigs.GRCh38.samtobreak.tsv

abyss2 discovardenovo discovardenovo-besst sga soapdenovo supernova:
	$(MAKE) draft=$@ \
		$@.abyss-fac.tsv \
		$@.scaftigs.GRCh38.samtobreak.tsv \
		$@.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.genomecov.stats.tsv \
		$@.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.abyss-fac.tsv \
		$@.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.scaftigs.GRCh38.samtobreak.tsv \
		$@.hg004.c$c_e$e_r$r.arcs.a$a_l$l.links.abyss-fac.tsv \
		$@.hg004.c$c_e$e_r$r.arcs.a$a_l$l.links.scaftigs.GRCh38.samtobreak.tsv \
		$@.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.hg004.c$c_e$e_r$r.arcs.a$a_l$l.links.abyss-fac.tsv \
		$@.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.hg004.c$c_e$e_r$r.arcs.a$a_l$l.links.scaftigs.GRCh38.samtobreak.tsv

sim:
	$(MAKE) draft=sim.abyss sample=sim.lr \
		sim.abyss.abyss-fac.tsv \
		sim.abyss.scaftigs.GRCh38.samtobreak.tsv \
		sim.abyss.sim.lr.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.genomecov.stats.tsv \
		sim.abyss.sim.lr.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.abyss-fac.tsv \
		sim.abyss.sim.lr.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.scaftigs.GRCh38.samtobreak.tsv \
		sim.abyss.sim.lr.c$c_e$e_r$r.arcs.a$a_l$l.links.abyss-fac.tsv \
		sim.abyss.sim.lr.c$c_e$e_r$r.arcs.a$a_l$l.links.scaftigs.GRCh38.samtobreak.tsv \
		sim.abyss.sim.lr.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.sim.lr.c$c_e$e_r$r.arcs.a$a_l$l.links.abyss-fac.tsv \
		sim.abyss.sim.lr.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.sim.lr.c$c_e$e_r$r.arcs.a$a_l$l.links.scaftigs.GRCh38.samtobreak.tsv

reports: \
	hg004.window2000.span20.s5000_n20.quast.metrics.html

tables: \
	hg004.window2000.span20.s5000_n20.quast.metrics.tsv.md

nxrepair: \
	abyss2.hg004.nxrepair.fa \
	abyss2.hg004.nxrepair.abyss-fac.tsv \
	abyss2.hg004.nxrepair.scaftigs.GRCh38.samtobreak.tsv \
	abyss2.hg004.nxrepair.hg004.c$c_e$e_r$r.arcs.a$a_l$l.links.abyss-fac.tsv \
	abyss2.hg004.nxrepair.hg004.c$c_e$e_r$r.arcs.a$a_l$l.links.scaftigs.GRCh38.samtobreak.tsv

# Download assemblies from NCBI GIAB.

# ABySS 2.0
abyss2.fa:
	curl -o $@ ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/BCGSC_HG004_ABySS2.0_assemblies_12082016/abyss-2.0/scaffolds.fa

# DISCOVARdenovo
discovardenovo.fa:
	curl -o $@ ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/BCGSC_HG004_ABySS2.0_assemblies_12082016/discovar/contigs.fa

# DISCOVARdenovo + ABySS-scaffold
discovardenovo-abyss.fa:
	curl -o $@ https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/BCGSC_HG004_ABySS2.0_assemblies_12082016/discovar/abyss-scaffolds.fa

# DISCOVARdenovo + BESST
discovardenovo-besst.fa:
	curl -o $@ ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/BCGSC_HG004_ABySS2.0_assemblies_12082016/discovar/besst-scaffolds.fa

# SGA
sga.fa:
	curl -o $@ ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/BCGSC_HG004_ABySS2.0_assemblies_12082016/sga/contigs.fa

# SOAPdenovo
soapdenovo.fa:
	curl -o $@ ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/BCGSC_HG004_ABySS2.0_assemblies_12082016/soapdenovo/scaffolds.fa

# Falcon
falcon.fa:
	curl -o $@ ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/MtSinai_PacBio_Assembly_falcon_03282016/hg004_p_and_a_ctg.fa

# PBcR
pbcr.fa:
	curl https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/549/595/GCA_001549595.1_GIAB_Ashkenazim_Mother_HG004_NA24143_hu8E87A9_PacBio_Assembly_with_PBcR/GCA_001549595.1_GIAB_Ashkenazim_Mother_HG004_NA24143_hu8E87A9_PacBio_Assembly_with_PBcR_genomic.fna.gz | gunzip -c >$@

# Rename the LINKS scaffolds.
abyss2_bionano_arcs.fa: abyss2_bionano_arcs.orig.fa
	gsed 's/scaffold//;s/,[^\t]*//' $< >$@

# Supernova 1.1

# Assemble the linked reads with Supernova
supernova.stamp: data/hg004lr/stamp
	supernova run --id=supernova --fastqs=$(<D) --localcores=64 --localmem=500
	touch $@

# Supernova 2

# Assemble the linked reads with Supernova
supernova2.stamp: data/hg004lr/stamp
	supernova run --id=supernova2 --fastqs=$(<D) --localcores=$t --localmem=256
	touch $@

# Generate the assembled pseudohap FASTA file.
supernova2.fa: supernova2.stamp
	supernova mkoutput --style=pseudohap --asmdir=supernova2/outs/assembly --outprefix=supernova2
	gunzip -c supernova2.fasta.gz >$@

# NA12878
# See https://github.com/nanopore-wgs-consortium/NA12878

# Downloads the 10x Genomics Chromium linked reads.
# See https://support.10xgenomics.com/de-novo-assembly/datasets/2.0.0/wfu/
data/wfu_fastqs.tar:
	curl -o $@ http://s3-us-west-2.amazonaws.com/10x.files/samples/assembly/2.0.0/wfu/wfu_fastqs.tar

# Extract the 10x Genomics Chromium linked reads.
data/wfu.stamp: data/wfu_fastqs.tar
	tar -C data -xf $<
	mv data/wfu/HNJJKCCXX data/HNJJKCCXXlr
	mv data/wfu/HNKKVCCXX data/HNKKVCCXXlr
	rmdir data/wfu
	touch $@

# Create stamp files for each flowcell.
data/HNJJKCCXXlr.stamp data/HNKKVCCXXlr.stamp: data/wfu.stamp
	touch $@

# Concatenate the 10x Genomics Chromium linked reads.
na12878.lrbasic.fq.gz: HNJJKCCXX.lrbasic.fq.gz HNKKVCCXX.lrbasic.fq.gz
	cat $^ >$@

# Download the Canu assembly.
# See https://genomeinformatics.github.io/NA12878-nanopore-assembly/
na12878.canu.fa:
	curl -o $@ http://s3.amazonaws.com/nanopore-human-wgs/canu.35x.contigs.polished2.fasta

# Download the Supernova 2 assembly.
# See https://support.10xgenomics.com/de-novo-assembly/datasets/2.0.0/wfu
na12878.supernova2.fa:
	curl http://cf.10xgenomics.com/samples/assembly/2.0.0/wfu/wfu_pseudohap.fasta.gz | gunzip -c >$@

# Longranger

# Extract 10x Chromium barcodes using longranger basic.
%.lrbasic.fq.gz: data/%lr.stamp
	longranger basic --id=$*_lrbasic --fastqs=data/$*lr
	ln -sf $*_lrbasic/outs/barcoded.fastq.gz $@

# Symlink the longranger basic FASTQ file.
sim.lr.lrbasic.fq.gz: simlr_lrbasic/outs/barcoded.fastq.gz
	ln -sf $< $@

# Barcodes

# Add the barcode to the read ID, and skip reads without barcodes.
%.bx.fq.gz: %.lrbasic.fq.gz
	gunzip -c $< | gawk ' \
		{ bx = "NA" } \
		match($$0, "BX:Z:([ACGT]*)-1", x) { bx = x[1] } \
		bx == "NA" { getline; getline; getline; next } \
		{ print $$1 "_" bx " " $$2; getline; print; getline; print; getline; print }' \
		| $(gzip) >$@

# BWA

# Index the genome.
%.fa.bwt: %.fa
	bwa index $<

# Align mate-pair reads to the draft genome and sort.
$(draft).%.mp.sort.bam: %.mp.fq.gz $(draft).fa.bwt
	rm -f $(TMPDIR)/$@.* || true
	bwa mem -t$t -p $(draft).fa $< | samtools sort -@$t -T$(TMPDIR)/$@ -o $@

# Align linked reads to the draft genome and sort.
$(draft).%.bam: %.fq.gz $(draft).fa.bwt
	rm -f $(TMPDIR)/$@.* || true
	bwa mem -t$t -pC $(draft).fa $< | samtools view -u -F4 | samtools sort -@$t -T$(TMPDIR)/$@ -o $@

# samtools

# Index a FASTA file.
%.fa.fai: %.fa
	samtools faidx $<

# Index a BAM file.
%.bam.bai: %.bam
	samtools index $<

# Remove alignments with an alignment score below a threshold.
%.as100.bam: %.bam
	samtools view -h -F4 $< | gawk -F'\t' ' \
			/^@/ { print; next } \
			{ as = 0 } \
			match($$0, "AS:.:([^\t]*)", x) { as = x[1] } \
			as >= 100' \
		| samtools view -@$t -b -o $@

# Select alignments with number of mismatches below a threshold.
nm=5
%.nm$(nm).bam: %.bam
	samtools view -h -F4 $< \
		| gawk -F'\t' ' \
			/^@/ { print; next } \
			{ nm = 999999999 } \
			match($$0, "NM:i:([^\t]*)", x) { nm = x[1] } \
			nm < $(nm)' \
		| samtools view -@$t -b -o $@

# Extract the alignment score, barcode and molecule identifier.
%.bam.bx.tsv: %.bam
	samtools view -F4 $< | gawk -F'\t' ' \
		BEGIN { print "Flags\tRname\tPos\tMapq\tAS\tNM\tBX\tMI" } \
		{ as = bx = mi = nm = "NA" } \
		match($$0, "AS:.:([^\t]*)", x) { as = x[1] } \
		match($$0, "NM:.:([^\t]*)", x) { nm = x[1] } \
		match($$0, "BX:Z:([^\t]*)", x) { bx = x[1] } \
		match($$0, "MI:i:([^\t]*)", x) { mi = x[1] } \
		{ print $$2 "\t" $$3 "\t" $$4 "\t" $$5 "\t" as "\t" nm "\t" bx "\t" mi }' >$@

# Group reads into molecules and add molecule identifiers.
%.bam.mi.bx.tsv: %.bam.bx.tsv
	./mi.r $< $@

# Create a TSV file of molecule extents.
%.bx.molecule.tsv: %.bx.tsv
	mlr --tsvlite \
		then stats1 -g BX,MI,Rname -a count,min,p50,max -f Pos,Mapq,AS,NM \
		then rename Pos_min,Start,Pos_max,End,Mapq_p50,Mapq_median,AS_p50,AS_median,NM_p50,NM_median,Pos_count,Reads \
		then put '$$Size = $$End - $$Start' \
		then cut -o -f Rname,Start,End,Size,BX,MI,Reads,Mapq_median,AS_median,NM_median \
		then filter '$$Reads >= 4' \
		$< >$@

# Create a BED file of molecule extents.
%.bx.molecule.bed: %.bx.molecule.tsv
	mlr --tsvlite --headerless-csv-output \
		put '$$Start = $$Start - 1; $$End = $$End - 1' \
		then put '$$Name = "Reads=" . $$Reads . ",Size=" . $$Size . ",Mapq=" . $$Mapq_median . ",AS=" . $$AS_median . ",NM=" . $$NM_median . ",BX=" . $$BX . ",MI=" . $$MI' \
		then cut -o -f Rname,Start,End,Name,Reads $< >$@

# Report summary statistics of a Chromium library
%.bx.molecule.summary.html: %.bx.molecule.tsv
	Rscript -e 'rmarkdown::render("summary.rmd", "html_document", "$@", params = list(input_tsv="$<", output_tsv="$*.summary.tsv"))'

# bedtools

# Convert BED to BAM.
%.bed.bam: %.bed $(draft).fa.fai
	awk '$$2 != $$3' $< | bedtools bedtobam -i - -g $(draft).fa.fai | samtools sort -@$t -Obam -o $@

# Compute the depth of coverage of a BAM file.
%.bam.depth.tsv: %.bam
	(printf "Rname\tPos\tDepth\n"; bedtools genomecov -d -ibam $<) >$@

# Compute the depth of coverage of a BED file.
%.bed.depth.tsv: %.bed $(draft).fa.fai
	(printf "Rname\tPos\tDepth\n"; awk '$$2 != $$3' $< | bedtools genomecov -d -g $(draft).fa.fai -i -) >$@

# Calculate depth of coverage statistics.
%.depth.stats.tsv: %.depth.tsv
	mlr --tsvlite stats1 -a count,p25,p50,p75,mean,stddev -f Depth $< >$@

# Calculate depth of coverage statistics per sequence.
%.depth.seqstats.tsv: %.depth.tsv
	mlr --tsvlite stats1 -a count,p25,p50,p75,mean,stddev -f Depth -g Rname $< >$@

# Compute statistics on the depth of coverage of a BED file.
%.bed.genomecov.tsv: %.bed $(draft).fa.fai
	(printf "Rname\tDepth\tCount\tRsize\tFraction\n"; awk '$$2 != $$3' $< | bedtools genomecov -g $(draft).fa.fai -i -) >$@

# Calculate depth of coverage statistics from bedtools genomecov.
%.genomecov.stats.tsv: %.genomecov.tsv
	mlr --tsvlite \
		then filter '$$Rname == "genome" && $$Depth > 0' \
		then step -a rsum -f Fraction \
		then put -q '@Depth_count += $$Count; if (is_null(@p25) && $$Fraction_rsum >= 0.25) { @p25 = $$Depth }; if (is_null(@p50) && $$Fraction_rsum >= 0.50) { @p50 = $$Depth }; if (is_null(@p75) && $$Fraction_rsum >= 0.75) { @p75 = $$Depth } end { emitf @Depth_count, @p25, @p50, @p75 }' \
		then rename p25,Depth_p25,p50,Depth_p50,p75,Depth_p75 \
		then put '$$Depth_IQR = $$Depth_p75 - $$Depth_p25' \
		$< >$@

# Identify breakpoints

# Select BED records of a given size or larger.
size_threshold=2000
%.size$(size_threshold).bed: %.bed
	awk '$$3 - $$2 >= $(size_threshold)' $< >$@

# Count start positions of molecules larger than a threshold size.
starts_size_threshold=2000
%.molecule.starts.tsv: %.molecule.tsv
	mlr --tsvlite \
		then filter '$$Size >= $(starts_size_threshold)' \
		then count-distinct -f Rname,Start \
		then rename Start,Pos,count,Starts \
		then sort -f Rname -n Pos \
		$< >$@

# Select position below the coverage threshold.
depth_threshold=100
%.depth.$(depth_threshold).tsv: %.depth.tsv
	mlr --tsvlite filter '$$Depth < $(depth_threshold)' $< >$@

# Join the tables of depth of coverage and number of molecule starts.
%.molecule.size$(size_threshold).depth.$(depth_threshold).starts.tsv: %.molecule.size$(size_threshold).bed.depth.$(depth_threshold).tsv %.molecule.starts.tsv
	mlr --tsvlite join -u -j Rname,Pos -f $^ >$@

# Select positions with low depth of coverage and high numer of molecule starts.
starts_threshold=2
pos_threshold=1000
%.depth.$(depth_threshold).starts.$(starts_threshold).breakpoints.tsv: %.depth.$(depth_threshold).starts.tsv
	mlr --tsvlite filter '$$Depth < $(depth_threshold) && $$Starts >= $(starts_threshold) && $$Pos >= $(pos_threshold)' $< >$@

# Identify breakpoints with low depth of coverage and high number of molecule starts.
%.size$(size_threshold).depth.starts.breakpoints.tsv: %.size$(size_threshold).bed.depth.tsv %.starts.tsv
	Rscript -e 'rmarkdown::render("breakpoints.rmd", "html_notebook", "$*.depth.starts.breakpoints.nb.html", params = list(depth_tsv="$<", starts_tsv="$*.starts.tsv", depth_starts_tsv="$*.depth.starts.tsv", breakpoints_tsv="$@"))'

# Group breakpoints by position.
group_threshold=1000
%.breakpoints.grouped.tsv: %.breakpoints.tsv
	mlr --tsvlite \
		then step -a delta -g Rname -f Pos \
		then put '$$NewGroup = $$Pos_delta > 0 && $$Pos_delta < $(group_threshold) ? 0 : 1' \
		then step -a rsum -f NewGroup \
		then cut -x -f NewGroup \
		then rename NewGroup_rsum,BreakpointID \
		$< >$@

# Count the number of breakpoints.
%.breakpoints.grouped.count.tsv: %.breakpoints.grouped.tsv
	mlr --tsvlite then stats1 -a max -f BreakpointID then rename BreakpointID_max,PP then put '$$File = FILENAME' $< >$@

# Determine coordinates of the breaktig subsequences.
%.breakpoints.tigs.bed: %.breakpoints.tsv $(draft).fa.fai
	Rscript -e 'rmarkdown::render("breaktigs.rmd", "html_notebook", "$*.breakpoints.tigs.nb.html", params = list(input_tsv="$<", input_fai="$(draft).fa.fai", output_bed="$@"))'

# Break scaffolds at the breakpoints.
%.breakpoints.tigs.fa: %.breakpoints.tigs.bed $(draft).fa
	bedtools getfasta -name -fi $(draft).fa -bed $< \
		| sed 's/::/ /;s/^NN*//;s/NN*$$//;s/^$$/N/' >$@

################################################################################
# Tigmint-span

# Parameters of Tigmint-span
sample=hg004
reads=$(sample).bx
dist=50000
span=20
window=2000

# Correct misassemblies using Tigmint.
%.$(reads).as0.65.nm5.molecule.size2000.trim0.window$(window).span$(span).breaktigs.fa: %.fa.bwt $(reads).fq.gz
	command time -v -o $@.all.time /home/sjackman/work/tigmint/tigmint-make t=$t dist=$(dist) window=$(window) span=$(span) draft=$* reads=$(reads) $@

# Correct misassemblies using Tigmint, setting the dist parameter.
%.$(reads).as0.65.nm5.dist$(dist).molecule.size2000.trim0.window$(window).span$(span).breaktigs.fa: %.fa.bwt $(reads).fq.gz
	command time -v -o $@.all.time /home/sjackman/work/tigmint/tigmint-make t=$t dist=$(dist) window=$(window) span=$(span) draft=$* reads=$(reads) $@

# Symlink the results of Tigmint-span and ARCS.
%.tigmint-span.arcs.fa: %.$(reads).as0.65.nm5.molecule.size2000.trim0.window$(window).span$(span).breaktigs.$(sample).c$c_e$e_r$r.arcs.a$a_l$l.links.fa
	ln -sf $< $@

# Symlink the result of Tigmint-span.
%.tigmint-span.fa: %.$(reads).as0.65.nm5.molecule.size2000.trim0.window$(window).span$(span).breaktigs.fa
	ln -sf $< $@

# Calculate assembly contiguity and correctness metrics using QUAST.
%.tigmint-span.quast.tsv: %.fa %.tigmint-span.fa %.arcs.fa %.tigmint-span.arcs.fa
	$(time) ~/.linuxbrew/bin/quast.py -t$t -se --fast --large --scaffold-gap-max-size 100000 --min-identity 90 -R $(ref_fa) -o $*.tigmint-span.quast $^
	ln -sf $*.tigmint-span.quast/transposed_report.tsv $@

# Calculate assembly contiguity and correctness metrics using QUAST.
%.tigmint-span.window$(window).span$(span).quast.tsv: %.$(reads).as0.65.nm5.molecule.size2000.trim0.window$(window).span$(span).breaktigs.fa %.$(reads).as0.65.nm5.molecule.size2000.trim0.window$(window).span$(span).breaktigs.$(sample).c$c_e$e_r$r.arcs.a$a_l$l.links.fa
	$(time) ~/.linuxbrew/bin/quast.py -t$t -se --fast --large --scaffold-gap-max-size 100000 --min-identity 90 -R $(ref_fa) -o $*.tigmint-span.window$(window).span$(span).quast $^
	ln -sf $*.tigmint-span.window$(window).span$(span).quast/transposed_report.tsv $@

# Aggregate the QUAST results.
assemblies.tigmint-span.quast.tsv: abyss2.tigmint-span.quast.tsv discovardenovo-besst.tigmint-span.quast.tsv supernova.tigmint-span.quast.tsv
	mlr --tsvlite cat $^ >$@

# Count the number of cuts made by Tigmint in scaffolds larger than the QUAST's sequence length threshold.
%.$(reads).as0.65.nm5.molecule.size2000.trim0.window$(window).span$(span).breaktigs.cuts.tsv: %.$(reads).as0.65.nm5.molecule.size2000.trim0.window$(window).span$(span).breaktigs.bed %.fa.fai
	(printf "Assembly\tCuts\n$*.$(reads).as0.65.nm5.molecule.size2000.trim0.window$(window).span$(span).breaktigs\t"; \
	mlr --nidx --fs tab --no-mmap --from $< \
		then cut -f 1,4 \
		then join -f <(cut -f 1,2 $*.fa.fai) -j 1 \
		then filter '$$2 >= 3000' \
		| grep -c -- '-[0-9]*[24680]$$') >$@

################################################################################
# Calculate assembly contiguity and correctness metrics.

# Reference genome
ref=GRCh38
ref_fa=$(ref).fa

# Size of the reference genome with Ns
GwithN=3088269832

# Path to ABySS executables
abyss_bin=/gsc/btl/linuxbrew/Cellar/abyss/2.0.2/bin

# BWA

# Align an assembly to the reference using BWA-MEM.
%.$(ref).sam.gz: %.fa
	time bwa mem -xintractg -t$t $(ref_fa) $< | $(gzip) >$@

# Align paired-end reads to the draft genome and do not sort.
%.$(sample).bx.sortn.bam: %.fa.bwt $(sample).bx.fq.gz
	$(time) bwa mem -t$t -pC $*.fa $(sample).bx.fq.gz | samtools view -@$t -F4 -b -o $@

# ABySS

# Convert scaffolds to scaftigs.
%.scaftigs.fa: %.fa
	seqtk seq $< | tr _ '~' | $(abyss_bin)/abyss-fatoagp -f $@ >$@.agp

# Calculate assembly contiguity metrics with abyss-fac.
%.abyss-fac.tsv: %.fa
	$(abyss_bin)/abyss-fac -G$(GwithN) -t500 $< >$@

# Calculate assembly contiguity and correctness metrics with abyss-samtobreak.
%.samtobreak.tsv: %.sam.gz
	gunzip -c $< | bin/abyss-samtobreak -G$(GwithN) -q10 -l500 >$@

# Calculate assembly contiguity and correctness metrics with abyss-samtobreak.
%.samtobreak.txt: %.sam.gz
	bin/abyss-samtobreak --text -G$(GwithN) -q10 -l500 $< >$@

# Convert samtobreak.txt to TSV using Miller.
%.samtobreak.txt.tsv: %.samtobreak.txt
	mlr --ixtab --ips ': ' --otsvlite --from $< \
		then rename 'Number of unmapped contigs,Unmapped_contigs' \
		then rename 'Total length of unmapped contigs,Unmapped_bases' \
		then rename 'Mapped contig bases,Mapped_bases' \
		then rename 'Contig NGA50,Contig_NGA50' \
		then rename 'Number of Q10 break points longer than 500 bp,Contig_breakpoints' \
		then rename 'Scaffold NG50,Scaffold_NG50' \
		then rename 'Scaffold NGA50,Scaffold_NGA50' \
		then rename 'Number of Q10 scaffold breakpoints longer than 500 bp,Scaffold_breakpoints' \
		then cut -r -x -f ' ' \
		then put '$$Total_breakpoints = $$Contig_breakpoints + $$Scaffold_breakpoints' \
		>$@

# ARCS

# Create a graph of linked contigs using ARCS.
%.$(sample).c$c_e$e_r$r.arcs_original.gv %.$(sample).c$c_e$e_r$r.arcs.dist.gv %.$(sample).c$c_e$e_r$r.arcs.dist.tsv: %.$(sample).bx.sortn.bam %.fa
	bin/arcs -s98 -c$c -l0 -z500 -m4-20000 -d0 -e$e -r$r -v \
		-f $*.fa \
		-b $*.$(sample).c$c_e$e_r$r.arcs \
		-g $*.$(sample).c$c_e$e_r$r.arcs.dist.gv \
		--tsv=$*.$(sample).c$c_e$e_r$r.arcs.dist.tsv \
		--barcode-counts=$<.barcode-counts.tsv \
		$<

# Convert the ARCS graph to LINKS TSV format.
%.$(sample).c$c_e$e_r$r.arcs.links.tsv: %.$(sample).c$c_e$e_r$r.arcs_original.gv %.fa
	bin/arcs-makeTSVfile $< $@ $*.fa

# Scaffold the assembly using the ARCS graph and LINKS.
%.$(sample).c$c_e$e_r$r.arcs.a$a_l$l_z$z.links.scaffolds.fa %.$(sample).c$c_e$e_r$r.arcs.a$a_l$l_z$z.links.assembly_correspondence.tsv: %.$(sample).c$c_e$e_r$r.arcs.links.tsv %.fa
	cp $< $*.$(sample).c$c_e$e_r$r.arcs.a$a_l$l_z$z.links.tigpair_checkpoint.tsv
	bin/LINKS -k20 -l$l -t2 -a$a -x1 -z$z -s /dev/null -f $*.fa -b $*.$(sample).c$c_e$e_r$r.arcs.a$a_l$l_z$z.links

# Rename the scaffolds.
%.links.fa: %.links.scaffolds.fa
	gsed -r 's/^>scaffold([^,]*),(.*)/>\1 scaffold\1,\2/' $< >$@

# Scaffold the assembly using ARCS and ABySS-scaffold.
%.$(sample).c$c_e$e_r$r.arcs.n$n.abyss-scaffold.path: %.$(sample).c$c_e$e_r$r.arcs.dist.gv
	abyss-scaffold -k144 -s1000-100000 -n$n -G$(GwithN) -o $@ $<

%.$(sample).c$c_e$e_r$r.arcs.n$n.abyss-scaffold.fa: %.fa %.fa.fai %.$(sample).c$c_e$e_r$r.arcs.n$n.abyss-scaffold.path
	MergeContigs -v -k144 -o $@ $^

%.$(sample).c$c_e$e_r$r.arcs.s$s_n$n.abyss-scaffold.path: %.$(sample).c$c_e$e_r$r.arcs.dist.gv
	abyss-scaffold -k144 -s$s -n$n -G$(GwithN) -o $@ $<

%.$(sample).c$c_e$e_r$r.arcs.s$s_n$n.abyss-scaffold.fa: %.fa %.fa.fai %.$(sample).c$c_e$e_r$r.arcs.s$s_n$n.abyss-scaffold.path
	MergeContigs -v -k144 -o $@ $^

# Aggregate the results.

# Symlink the assemblies.
%.tigmint.arcs.fa: %.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.hg004.c$c_e$e_r$r.arcs.a$a_l$l.links.fa
	ln -sf $< $@

%.tigmint.fa: %.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.fa
	ln -sf $< $@

%.arcs.fa: %.$(sample).c$c_e$e_r$r.arcs.a$a_l$l.links.fa
	ln -sf $< $@

%.depth.100.starts.2.abyss-fac.tsv: \
		%.abyss-fac.tsv \
		%.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.abyss-fac.tsv \
		%.hg004.c$c_e$e_r$r.arcs.a$a_l$l.links.abyss-fac.tsv \
		%.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.hg004.c$c_e$e_r$r.arcs.a$a_l$l.links.abyss-fac.tsv
	mlr --tsvlite cat $^ >$@

%.depth.100.starts.2.scaftigs.abyss-fac.tsv: \
		%.scaftigs.abyss-fac.tsv \
		%.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.scaftigs.abyss-fac.tsv \
		%.hg004.c$c_e$e_r$r.arcs.a$a_l$l.links.scaftigs.abyss-fac.tsv \
		%.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.hg004.c$c_e$e_r$r.arcs.a$a_l$l.links.scaftigs.abyss-fac.tsv
	mlr --tsvlite cat $^ >$@

%.depth.100.starts.2.samtobreak.tsv: \
		%.scaftigs.GRCh38.samtobreak.tsv \
		%.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.scaftigs.GRCh38.samtobreak.tsv \
		%.hg004.c$c_e$e_r$r.arcs.a$a_l$l.links.scaftigs.GRCh38.samtobreak.tsv \
		%.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.hg004.c$c_e$e_r$r.arcs.a$a_l$l.links.scaftigs.GRCh38.samtobreak.tsv
	mlr --tsvlite cat $^ >$@

abyss2.depth.100.starts.1-5.breakpoints.count.tsv: \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.1.breakpoints.grouped.count.tsv \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.grouped.count.tsv \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.3.breakpoints.grouped.count.tsv \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.4.breakpoints.grouped.count.tsv \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.5.breakpoints.grouped.count.tsv
	mlr --tsvlite put 'FILENAME =~ "[.]depth[.]([0-9]*)[.]starts[.]([0-9]*)[.]"; $$Depth = "\1"; $$Starts = "\2"' $^ >$@

abyss2.depth.100.starts.1-5.samtobreak.tsv: \
		abyss2.scaftigs.GRCh38.samtobreak.tsv \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.1.breakpoints.tigs.scaftigs.GRCh38.samtobreak.tsv \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.scaftigs.GRCh38.samtobreak.tsv \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.3.breakpoints.tigs.scaftigs.GRCh38.samtobreak.tsv \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.4.breakpoints.tigs.scaftigs.GRCh38.samtobreak.tsv \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.5.breakpoints.tigs.scaftigs.GRCh38.samtobreak.tsv
	mlr --tsvlite put 'FILENAME =~ "[.]depth[.]([0-9]*)[.]starts[.]([0-9]*)[.]"; $$Depth = "\1"; $$Starts = "\2"' $^ >$@

abyss2.depth.80-120.starts.2-4.arcs.abyss-fac.tsv: \
		abyss2.hg004.c5_e30000_r0.05.arcs.a0.1_l10.links.abyss-fac.tsv \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.80.starts.2.breakpoints.tigs.hg004.c5_e30000_r0.05.arcs.a0.1_l10.links.abyss-fac.tsv \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.hg004.c5_e30000_r0.05.arcs.a0.1_l10.links.abyss-fac.tsv \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.120.starts.2.breakpoints.tigs.hg004.c5_e30000_r0.05.arcs.a0.1_l10.links.abyss-fac.tsv \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.hg004.c5_e30000_r0.05.arcs.a0.1_l10.links.abyss-fac.tsv \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.3.breakpoints.tigs.hg004.c5_e30000_r0.05.arcs.a0.1_l10.links.abyss-fac.tsv \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.4.breakpoints.tigs.hg004.c5_e30000_r0.05.arcs.a0.1_l10.links.abyss-fac.tsv
	mlr --tsvlite put 'FILENAME =~ "[.]depth[.]([0-9]*)[.]starts[.]([0-9]*)[.]"; $$Depth = "\1"; $$Starts = "\2"' $^ >$@

abyss2.depth.80-120.starts.2-4.arcs.samtobreak.tsv: \
		abyss2.hg004.c5_e30000_r0.05.arcs.a0.1_l10.links.scaftigs.GRCh38.samtobreak.tsv \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.80.starts.2.breakpoints.tigs.hg004.c5_e30000_r0.05.arcs.a0.1_l10.links.scaftigs.GRCh38.samtobreak.tsv \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.hg004.c5_e30000_r0.05.arcs.a0.1_l10.links.scaftigs.GRCh38.samtobreak.tsv \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.120.starts.2.breakpoints.tigs.hg004.c5_e30000_r0.05.arcs.a0.1_l10.links.scaftigs.GRCh38.samtobreak.tsv \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.2.breakpoints.tigs.hg004.c5_e30000_r0.05.arcs.a0.1_l10.links.scaftigs.GRCh38.samtobreak.tsv \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.3.breakpoints.tigs.hg004.c5_e30000_r0.05.arcs.a0.1_l10.links.scaftigs.GRCh38.samtobreak.tsv \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.4.breakpoints.tigs.hg004.c5_e30000_r0.05.arcs.a0.1_l10.links.scaftigs.GRCh38.samtobreak.tsv
	mlr --tsvlite put 'FILENAME =~ "[.]depth[.]([0-9]*)[.]starts[.]([0-9]*)[.]"; $$Depth = "\1"; $$Starts = "\2"' $^ >$@

assemblies.depth.100.starts.2.abyss-fac.tsv: \
		abyss2.depth.100.starts.2.abyss-fac.tsv \
		discovardenovo-besst.depth.100.starts.2.abyss-fac.tsv \
		supernova.depth.100.starts.2.abyss-fac.tsv
	mlr --tsvlite cat $^ >$@

assemblies.depth.100.starts.2.samtobreak.tsv: \
		abyss2.depth.100.starts.2.samtobreak.tsv \
		discovardenovo-besst.depth.100.starts.2.samtobreak.tsv \
		soapdenovo.depth.100.starts.2.samtobreak.tsv \
		supernova.depth.100.starts.2.samtobreak.tsv
	mlr --tsvlite cat $^ >$@

# NxRepair

# Concatenate the mate-pair data.
hg004.mp.fq.gz: MPHG004-23100079.mp.fq.gz MPHG004-23110109.mp.fq.gz
	cat $^ >$@

# Correct assembly errors with NxRepair.
%.hg004.nxrepair.fa: %.hg004.mp.sort.bam %.fa %.hg004.mp.sort.bam.bai
	$(time) ~/.linuxbrew/bin/nxrepair $< $*.fa $*.hg004.nxrepair.tsv $@

# Correct assembly errors with NxRepair and set parameter T.
nxrepair_T=2.6
%.hg004.T$(nxrepair_T).nxrepair.fa: %.hg004.mp.sort.bam %.fa %.hg004.mp.sort.bam.bai
	$(time) ~/.linuxbrew/bin/nxrepair -T-$(nxrepair_T) $< $*.fa $*.hg004.T$(nxrepair_T).nxrepair.tsv $@

# wgsim

# Simulate paired-end reads using wgsim.
sim.pe.1.fq:
	wgsim -e 0.001 -r 0 -d 400 -s 90 -N 434296528 -1 250 -2 250 -S 1 $(ref_fa) sim.pe.1.fq sim.pe.2.fq

# Simulate mate-pair reads using wgsim.
sim.mp.1.fq:
	wgsim -e 0.001 -r 0 -d 6000 -s 1400 -N 350564272 -1 125 -2 125 -S 1 $(ref_fa) sim.mp.2.fq sim.mp.1.fq

# Interleave and compress paired-end reads.
%.fq.gz: %.1.fq %.2.fq
	seqtk mergepe $^ | $(gzip) >$@

# LRSIM

# Simulate linked reads using LRSIM.
sim.lr: $(ref).fa
	simulateLinkedReads -z $t -x 524 -d 1 -1 0 -4 0 -7 0 -0 0 -g $< -p $@
	touch $@

# ABySS
abyss=/gsc/btl/linuxbrew/Cellar/abyss/2.0.1-k192/bin

# Assemble the paired-end and mate-pair reads using ABySS.
abyss/%-scaffolds.fa: %.pe.fq.gz %.mp.fq.gz
	mkdir -p $(@D)
	$(time) $(abyss)/abyss-pe -C $(@D) \
		name=sim np=$t v=-v \
		k=144 q=15 B=26G H=4 kc=3 \
		l=40 s=1000 n=10 \
		S=1000-10000 N=20 mp6k_de=--mean mp6k_n=1 \
		lib=pe400 pe400=$(PWD)/$*.pe.fq.gz \
		mp=mp6k mp6k=$(PWD)/$*.mp.fq.gz

# Symlink the ABySS assembly.
sim.abyss.fa: abyss/sim-scaffolds.fa
	ln -sf $< $@

# QUAST

# Calculate assembly contiguity and correctness metrics using QUAST.
%.quast.tsv: %.fa
	~/.linuxbrew/bin/quast.py -t$t -se --fast --large --scaffold-gap-max-size 100000 --min-identity 90 -R $(ref_fa) -o $*.quast $<
	cp $*.quast/transposed_report.tsv $@

# ABySS-Scaffold with optimized s

# Aggregate the QUAST results of one assembler.
%.window$(window).span$(span).n$n.quast.tsv: \
		%.quast.tsv \
		%.$(sample).bx.as0.65.nm5.molecule.size2000.trim0.window$(window).span$(span).breaktigs.quast.tsv \
		%.$(sample).c$c_e$e_r$r.arcs.n$n.abyss-scaffold.quast.tsv \
		%.$(sample).bx.as0.65.nm5.molecule.size2000.trim0.window$(window).span$(span).breaktigs.$(sample).c$c_e$e_r$r.arcs.n$n.abyss-scaffold.quast.tsv
	mlr --tsvlite cat $^ >$@

# Aggregate the QUAST results of the HG004 assemblies.
hg004.window$(window).span$(span).n$n.quast.tsv: \
		abyss2.window$(window).span$(span).n$n.quast.tsv \
		discovardenovo-besst.window$(window).span$(span).n$n.quast.tsv \
		supernova2.window$(window).span$(span).n$n.quast.tsv
	mlr --tsvlite cat $^ >$@

# Aggregate the QUAST results of the NA12878 assemblies.
na12878.window$(window).span$(span).n$n.quast.tsv: \
		na12878.canu.window$(window).span$(span).n$n.quast.tsv \
		na12878.supernova2.window$(window).span$(span).n$n.quast.tsv
	mlr --tsvlite cat $^ >$@

# ABySS-Scaffold with fixed s

# Aggregate the QUAST results of one assembler.
%.window$(window).span$(span).s$s_n$n.quast.tsv: \
		%.quast.tsv \
		%.$(sample).bx.as0.65.nm5.molecule.size2000.trim0.window$(window).span$(span).breaktigs.quast.tsv \
		%.$(sample).c$c_e$e_r$r.arcs.s$s_n$n.abyss-scaffold.quast.tsv \
		%.$(sample).bx.as0.65.nm5.molecule.size2000.trim0.window$(window).span$(span).breaktigs.$(sample).c$c_e$e_r$r.arcs.s$s_n$n.abyss-scaffold.quast.tsv
	mlr --tsvlite cat $^ >$@

# Aggregate the QUAST results of the HG004 assemblies.
hg004.window$(window).span$(span).s$s_n$n.quast.tsv: \
		abyss2.window$(window).span$(span).s$s_n$n.quast.tsv \
		discovardenovo-abyss.window$(window).span$(span).s$s_n$n.quast.tsv \
		discovardenovo-besst.window$(window).span$(span).s$s_n$n.quast.tsv \
		supernova2.window$(window).span$(span).s$s_n$n.quast.tsv
	mlr --tsvlite cat $^ >$@

# Aggregate the QUAST results of the NA12878 assemblies.
na12878.window$(window).span$(span).s$s_n$n.quast.tsv: \
		na12878.canu.window$(window).span$(span).s$s_n$n.quast.tsv \
		na12878.supernova2.window$(window).span$(span).s$s_n$n.quast.tsv
	mlr --tsvlite cat $^ >$@

# Aggregate the QUAST results of the long reads assemblies.
sms.window$(window).span$(span).s$s_n$n.quast.tsv: \
		falcon.window$(window).span$(span).s$s_n$n.quast.tsv \
		na12878.canu.window$(window).span$(span).s$s_n$n.quast.tsv
	mlr --tsvlite cat $^ >$@

# Aggregate the QUAST results of some assemblies.
assemblies.window$(window).span$(span).s$s_n$n.quast.tsv: \
		abyss2.window$(window).span$(span).s$s_n$n.quast.tsv \
		falcon.window$(window).span$(span).s$s_n$n.quast.tsv \
		na12878.canu.window$(window).span$(span).s$s_n$n.quast.tsv \
		discovardenovo-besst.window$(window).span$(span).s$s_n$n.quast.tsv \
		supernova2.window$(window).span$(span).s$s_n$n.quast.tsv
	mlr --tsvlite cat $^ >$@

# Aggregate the abyss-samtobreak results of one assembler.
%.window$(window).span$(span).s$s_n$n.scaftigs.GRCh38.samtobreak.tsv: \
		%.scaftigs.GRCh38.samtobreak.tsv \
		%.$(sample).bx.as0.65.nm5.molecule.size2000.trim0.window$(window).span$(span).breaktigs.scaftigs.GRCh38.samtobreak.tsv \
		%.$(sample).c$c_e$e_r$r.arcs.s$s_n$n.abyss-scaffold.scaftigs.GRCh38.samtobreak.tsv \
		%.$(sample).bx.as0.65.nm5.molecule.size2000.trim0.window$(window).span$(span).breaktigs.$(sample).c$c_e$e_r$r.arcs.s$s_n$n.abyss-scaffold.scaftigs.GRCh38.samtobreak.tsv
	mlr --tsvlite cat $^ >$@

# Aggregate the abyss-samtobreak results of the HG004 assemblies.
hg004.window$(window).span$(span).s$s_n$n.scaftigs.GRCh38.samtobreak.tsv: \
		abyss2.window$(window).span$(span).s$s_n$n.scaftigs.GRCh38.samtobreak.tsv \
		falcon.window$(window).span$(span).s$s_n$n.scaftigs.GRCh38.samtobreak.tsv \
		discovardenovo-abyss.window$(window).span$(span).s$s_n$n.scaftigs.GRCh38.samtobreak.tsv \
		discovardenovo-besst.window$(window).span$(span).s$s_n$n.scaftigs.GRCh38.samtobreak.tsv \
		supernova2.window$(window).span$(span).s$s_n$n.scaftigs.GRCh38.samtobreak.tsv
	mlr --tsvlite cat $^ >$@

# Aggregate the abyss-samtobreak results of the NA12878 assemblies.
na12878.window$(window).span$(span).s$s_n$n.scaftigs.GRCh38.samtobreak.tsv: \
		na12878.canu.window$(window).span$(span).s$s_n$n.scaftigs.GRCh38.samtobreak.tsv \
		na12878.supernova2.window$(window).span$(span).s$s_n$n.scaftigs.GRCh38.samtobreak.tsv
	mlr --tsvlite cat $^ >$@

# LINKS

# Aggregate the QUAST results of one assembler.
%.window$(window).span$(span).a$a_l$l_z$z.quast.tsv: \
		%.quast.tsv \
		%.$(sample).bx.as0.65.nm5.molecule.size2000.trim0.window$(window).span$(span).breaktigs.quast.tsv \
		%.$(sample).c$c_e$e_r$r.arcs.a$a_l$l_z$z.links.quast.tsv \
		%.$(sample).bx.as0.65.nm5.molecule.size2000.trim0.window$(window).span$(span).breaktigs.$(sample).c$c_e$e_r$r.arcs.a$a_l$l_z$z.links.quast.tsv
	mlr --tsvlite cat $^ >$@

# Aggregate the QUAST results of the HG004 assemblies.
hg004.window$(window).span$(span).a$a_l$l_z$z.quast.tsv: \
		abyss2.window$(window).span$(span).a$a_l$l_z$z.quast.tsv \
		discovardenovo-besst.window$(window).span$(span).a$a_l$l_z$z.quast.tsv \
		supernova.window$(window).span$(span).a$a_l$l_z$z.quast.tsv
	mlr --tsvlite cat $^ >$@

# Aggregate the QUAST results of the NA12878 assemblies.
na12878.window$(window).span$(span).a$a_l$l_z$z.quast.tsv: \
		na12878.canu.window$(window).span$(span).a$a_l$l_z$z.quast.tsv \
		na12878.supernova2.window$(window).span$(span).a$a_l$l_z$z.quast.tsv
	mlr --tsvlite cat $^ >$@

# Aggregate the QUAST results of the HG004 parameter sweep.
hg004.parameters.s5000_n20.quast.tsv: \
		abyss2.window200.span20.s5000_n20.quast.tsv \
		abyss2.window500.span20.s5000_n20.quast.tsv \
		abyss2.window1000.span20.s5000_n20.quast.tsv \
		abyss2.window2000.span20.s5000_n20.quast.tsv \
		abyss2.window5000.span20.s5000_n20.quast.tsv \
		abyss2.window10000.span20.s5000_n20.quast.tsv \
		abyss2.window2000.span2.s5000_n20.quast.tsv \
		abyss2.window2000.span5.s5000_n20.quast.tsv \
		abyss2.window2000.span10.s5000_n20.quast.tsv \
		abyss2.window2000.span20.s5000_n20.quast.tsv \
		abyss2.window2000.span50.s5000_n20.quast.tsv \
		abyss2.window2000.span100.s5000_n20.quast.tsv \
		discovardenovo-abyss.hg004.bx.as0.65.nm5.molecule.size2000.trim0.window200.span20.breaktigs.quast.tsv \
		discovardenovo-abyss.hg004.bx.as0.65.nm5.molecule.size2000.trim0.window500.span20.breaktigs.quast.tsv \
		discovardenovo-abyss.hg004.bx.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.quast.tsv \
		discovardenovo-abyss.hg004.bx.as0.65.nm5.molecule.size2000.trim0.window2000.span20.breaktigs.quast.tsv \
		discovardenovo-abyss.hg004.bx.as0.65.nm5.molecule.size2000.trim0.window5000.span20.breaktigs.quast.tsv \
		discovardenovo-abyss.hg004.bx.as0.65.nm5.molecule.size2000.trim0.window10000.span20.breaktigs.quast.tsv \
		discovardenovo-abyss.hg004.bx.as0.65.nm5.molecule.size2000.trim0.window2000.span2.breaktigs.quast.tsv \
		discovardenovo-abyss.hg004.bx.as0.65.nm5.molecule.size2000.trim0.window2000.span5.breaktigs.quast.tsv \
		discovardenovo-abyss.hg004.bx.as0.65.nm5.molecule.size2000.trim0.window2000.span10.breaktigs.quast.tsv \
		discovardenovo-abyss.hg004.bx.as0.65.nm5.molecule.size2000.trim0.window2000.span20.breaktigs.quast.tsv \
		discovardenovo-abyss.hg004.bx.as0.65.nm5.molecule.size2000.trim0.window2000.span50.breaktigs.quast.tsv \
		discovardenovo-abyss.hg004.bx.as0.65.nm5.molecule.size2000.trim0.window2000.span100.breaktigs.quast.tsv \
		discovardenovo-besst.window200.span20.s5000_n20.quast.tsv \
		discovardenovo-besst.window500.span20.s5000_n20.quast.tsv \
		discovardenovo-besst.window1000.span20.s5000_n20.quast.tsv \
		discovardenovo-besst.window2000.span20.s5000_n20.quast.tsv \
		discovardenovo-besst.window5000.span20.s5000_n20.quast.tsv \
		discovardenovo-besst.window10000.span20.s5000_n20.quast.tsv \
		discovardenovo-besst.window2000.span2.s5000_n20.quast.tsv \
		discovardenovo-besst.window2000.span5.s5000_n20.quast.tsv \
		discovardenovo-besst.window2000.span10.s5000_n20.quast.tsv \
		discovardenovo-besst.window2000.span20.s5000_n20.quast.tsv \
		discovardenovo-besst.window2000.span50.s5000_n20.quast.tsv \
		discovardenovo-besst.window2000.span100.s5000_n20.quast.tsv
	mlr --tsvlite cat $^ | sort -u >$@

# Aggregate the QUAST results of the simulation parameter sweep.
sim.abyss.parameters.s5000_n20.quast.tsv: \
		sim.abyss.sim.lr.bx.as0.65.nm5.molecule.size2000.trim0.window1000.span2.breaktigs.quast.tsv \
		sim.abyss.sim.lr.bx.as0.65.nm5.molecule.size2000.trim0.window1000.span5.breaktigs.quast.tsv \
		sim.abyss.sim.lr.bx.as0.65.nm5.molecule.size2000.trim0.window1000.span10.breaktigs.quast.tsv \
		sim.abyss.sim.lr.bx.as0.65.nm5.molecule.size2000.trim0.window1000.span20.breaktigs.quast.tsv \
		sim.abyss.sim.lr.bx.as0.65.nm5.molecule.size2000.trim0.window1000.span50.breaktigs.quast.tsv \
		sim.abyss.sim.lr.bx.as0.65.nm5.molecule.size2000.trim0.window1000.span100.breaktigs.quast.tsv \
		sim.abyss.sim.lr.bx.as0.65.nm5.molecule.size2000.trim0.window2000.span20.breaktigs.quast.tsv
	mlr --tsvlite cat $^ >$@

# RMarkdown reports

# Compute the assembly metrics for a set of assemblies.
%.quast.metrics.html %.quast.metrics.tsv: %.quast.tsv
	Rscript -e 'rmarkdown::render("quast.rmd", "html_document", "$*.quast.metrics.html", params = list(input_tsv="$<", output_tsv="$*.quast.metrics.tsv"))'

# Compute the assembly metrics for a set of assemblies.
%.metrics.html %.metrics.tsv: %.abyss-fac.tsv %.samtobreak.tsv
	Rscript -e 'rmarkdown::render("metrics.rmd", "html_document", "$*.metrics.html", params = list(abyss_fac_tsv="$<", samtobreak_tsv="$*.samtobreak.tsv", output_tsv="$*.metrics.tsv"))'

# Compute the precision, recall, and G-score.
%.samtobreak.gscore.html %.samtobreak.gscore.tsv: %.breakpoints.count.tsv %.samtobreak.tsv
	Rscript -e 'rmarkdown::render("precision-recall.rmd", "html_document", "$*.samtobreak.gscore.html", params = list(breakpoints_count_tsv="$<", samtobreak_tsv="$*.samtobreak.tsv", output_tsv="$*.samtobreak.gscore.tsv"))'

# Compute the assembly metrics for parameter sensitivty of Tigmint-span.
%.quast.parameters.html %.quast.parameters.tsv: %.quast.tsv
	Rscript -e 'rmarkdown::render("tigmint-parameters.rmd", "html_document", "$*.quast.parameters.html", params = list(input_tsv="$<", output_tsv="$*.quast.parameters.tsv"))'

# Compute the assembly metrics for parameter sensitivty.
%.parameters.html %.parameters.tsv: %.abyss-fac.tsv %.samtobreak.tsv
	Rscript -e 'rmarkdown::render("parameters.rmd", "html_document", "$*.parameters.html", params = list(abyss_fac_tsv="$<", samtobreak_tsv="$*.samtobreak.tsv", output_tsv="$*.parameters.tsv"))'

# Convert TSV to Markdown with Miller.
%.tsv.md: %.tsv
	mlr --itsvlite --omd cat $< >$@
