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

# Parameters of LINKS
a=0.1
l=10

# Report run time and memory usage
time=env time -v -o $@.time
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

assemblies: abyss2 discovardenovo discovardenovo-besst sga soapdenovo supernova

quast: assemblies.quast.tsv

reports: \
	abyss2.depth.80-120.starts.2-4.arcs.parameters.html \
	assemblies.depth.100.starts.2.metrics.html

tables: \
	abyss2.depth.80-120.starts.2-4.arcs.parameters.tsv.md \
	assemblies.depth.100.starts.2.metrics.tsv.md

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

# DISCOVARdenovo + BESST
discovardenovo-besst.fa:
	curl -o $@ ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/BCGSC_HG004_ABySS2.0_assemblies_12082016/discovar/besst-scaffolds.fa

# SGA
sga.fa:
	curl -o $@ ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/BCGSC_HG004_ABySS2.0_assemblies_12082016/sga/contigs.fa

# SOAPdenovo
soapdenovo.fa:
	curl -o $@ ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/BCGSC_HG004_ABySS2.0_assemblies_12082016/soapdenovo/scaffolds.fa

# Rename the LINKS scaffolds.
abyss2_bionano_arcs.fa: abyss2_bionano_arcs.orig.fa
	gsed 's/scaffold//;s/,[^\t]*//' $< >$@

# Supernova

# Assemble the linked reads with Supernova
NA24143-giab/outs/summary.csv:
	supernova run --id=NA24143-giab --fastqs=/projects/btl/datasets/hsapiens/giab/HG004/10XGenomics/NA24143-chromium-raw --localcores=64 --localmem=500

# Longranger

# Extract 10x Chromium barcodes using longranger basic.
%_lrbasic/outs/barcoded.fastq.gz: data/%/stamp
	longranger basic --id=$*_lrbasic --fastqs=$(<D)

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
	bwa mem -t$t -p $(draft).fa $< | samtools sort -@$t -o $@

# Align linked reads to the draft genome and sort.
$(draft).%.bam: %.fq.gz $(draft).fa.bwt
	bwa mem -t$t -pC $(draft).fa $< | samtools view -h -F4 | samtools sort -@$t -o $@

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
# Calculate assembly contiguity and correctness metrics.

# Reference genome
ref=GRCh38
ref_fa=$(ref).fa

# Size of the reference genome with Ns
GwithN=3088269832

# Size of the reference genome without Ns
GwithoutN=2937639113

# Path to ABySS executables
abyss_bin=/gsc/btl/linuxbrew/Cellar/abyss/2.0.2/bin

# BWA

# Align an assembly to the reference using BWA-MEM.
%.$(ref).sam: %.fa
	time bwa mem -xintractg -t$t $(ref_fa) $< >$@

# Align paired-end reads to the draft genome and do not sort.
sample=hg004
%.$(sample).bx.sortn.bam: %.fa.bwt $(sample).bx.fq.gz
	bwa mem -t$t -pC $*.fa $(sample).bx.fq.gz | samtools view -@$t -h -F4 -o $@

# ABySS

# Convert scaffolds to scaftigs.
%.scaftigs.fa: %.fa
	seqtk seq $< | tr _ '~' | $(abyss_bin)/abyss-fatoagp -f $@ >$@.agp

# Calculate assembly contiguity metrics with abyss-fac.
%.abyss-fac.tsv: %.fa
	$(abyss_bin)/abyss-fac -G$(GwithN) -t500 $< >$@

# Calculate assembly contiguity and correctness metrics with abyss-samtobreak.
%.samtobreak.txt: %.sam
	(echo "File: $<"; abyss-samtobreak --text -G$(GwithN) -q10 -l500 $<) >$@

# Convert samtobreak.txt to TSV using Miller.
%.samtobreak.tsv: %.samtobreak.txt
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
%.$(sample).c$c_e$e_r$r.arcs.a$a_l$l.links.scaffolds.fa %.$(sample).c$c_e$e_r$r.arcs.a$a_l$l.links.assembly_correspondence.tsv: %.$(sample).c$c_e$e_r$r.arcs.links.tsv %.fa
	cp $< $*.$(sample).c$c_e$e_r$r.arcs.a$a_l$l.links.tigpair_checkpoint.tsv
	LINKS -k20 -l$l -t2 -a$a -x1 -s /dev/null -f $*.fa -b $*.$(sample).c$c_e$e_r$r.arcs.a$a_l$l.links

# Rename the scaffolds.
%.links.fa: %.links.scaffolds.fa
	gsed -r 's/^>scaffold([^,]*),(.*)/>\1 scaffold\1,\2/' $< >$@

# Aggregate the results.
s=$(starts_threshold)

# Symlink the assemblies.
%.tigmint.arcs.fa: %.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.$s.breakpoints.tigs.hg004.c$c_e$e_r$r.arcs.a$a_l$l.links.fa
	ln -sf $< $@

%.tigmint.fa: %.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.$s.breakpoints.tigs.fa
	ln -sf $< $@

%.arcs.fa: %.hg004.c$c_e$e_r$r.arcs.a$a_l$l.links.fa
	ln -sf $< $@

%.depth.100.starts.$s.abyss-fac.tsv: \
		%.abyss-fac.tsv \
		%.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.$s.breakpoints.tigs.abyss-fac.tsv \
		%.hg004.c$c_e$e_r$r.arcs.a$a_l$l.links.abyss-fac.tsv \
		%.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.$s.breakpoints.tigs.hg004.c$c_e$e_r$r.arcs.a$a_l$l.links.abyss-fac.tsv
	mlr --tsvlite cat $^ >$@

%.depth.100.starts.$s.scaftigs.abyss-fac.tsv: \
		%.scaftigs.abyss-fac.tsv \
		%.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.$s.breakpoints.tigs.scaftigs.abyss-fac.tsv \
		%.hg004.c$c_e$e_r$r.arcs.a$a_l$l.links.scaftigs.abyss-fac.tsv \
		%.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.$s.breakpoints.tigs.hg004.c$c_e$e_r$r.arcs.a$a_l$l.links.scaftigs.abyss-fac.tsv
	mlr --tsvlite cat $^ >$@

%.depth.100.starts.$s.samtobreak.tsv: \
		%.scaftigs.GRCh38.samtobreak.tsv \
		%.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.$s.breakpoints.tigs.scaftigs.GRCh38.samtobreak.tsv \
		%.hg004.c$c_e$e_r$r.arcs.a$a_l$l.links.scaftigs.GRCh38.samtobreak.tsv \
		%.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.$s.breakpoints.tigs.hg004.c$c_e$e_r$r.arcs.a$a_l$l.links.scaftigs.GRCh38.samtobreak.tsv
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
	$(time) nxrepair $< $*.fa $*.nxrepair $@

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
		S=1000-10000 N=7 mp6k_de=--mean mp6k_n=1 \
		lib=pe400 pe400=$(PWD)/$*.pe.fq.gz \
		mp=mp6k mp6k=$(PWD)/$*.mp.fq.gz

# Symlink the ABySS assembly.
sim.abyss.fa: abyss/sim-scaffolds.fa
	ln -sf $< $@

# QUAST

# Calculate assembly contiguity and correctness metrics using QUAST.
%.quast/transposed_report.tsv: %.fa %.tigmint.fa %.arcs.fa %.tigmint.arcs.fa
	~/.linuxbrew/bin/quast.py -t$t -se --fast --large --min-identity 90 --scaffold-gap-max-size 1000000 -R $(ref_fa) -o $(@D) $^

# Symlink the QUAST results.
%.quast.tsv: %.quast/transposed_report.tsv
	ln -sf $< $@

# Aggregate the QUAST results.
assemblies.quast.tsv: abyss2.quast.tsv discovardenovo-besst.quast.tsv supernova.quast.tsv
	mlr --tsvlite cat $^ >$@

# RMarkdown reports

# Compute the assembly metrics for a set of assemblies.
%.metrics.html %.metrics.tsv: %.abyss-fac.tsv %.samtobreak.tsv
	Rscript -e 'rmarkdown::render("metrics.rmd", "html_document", "$*.metrics.html", params = list(abyss_fac_tsv="$<", samtobreak_tsv="$*.samtobreak.tsv", output_tsv="$*.metrics.tsv"))'

# Compute the assembly metrics for a set of assemblies.
%.quast.nb.html %.quast.out.tsv: %.quast.tsv
	Rscript -e 'rmarkdown::render("quast.rmd", "html_notebook", "$*.quast.html", params = list(input_tsv="$<", output_tsv="$*.quast.out.tsv"))'

# Compute the precision, recall, and G-score.
%.samtobreak.gscore.html %.samtobreak.gscore.tsv: %.breakpoints.count.tsv %.samtobreak.tsv
	Rscript -e 'rmarkdown::render("precision-recall.rmd", "html_document", "$*.samtobreak.gscore.html", params = list(breakpoints_count_tsv="$<", samtobreak_tsv="$*.samtobreak.tsv", output_tsv="$*.samtobreak.gscore.tsv"))'

# Compute the assembly metrics for parameter sensitivty.
%.parameters.html %.parameters.tsv: %.abyss-fac.tsv %.samtobreak.tsv
	Rscript -e 'rmarkdown::render("parameters.rmd", "html_document", "$*.parameters.html", params = list(abyss_fac_tsv="$<", samtobreak_tsv="$*.samtobreak.tsv", output_tsv="$*.parameters.tsv"))'

# Convert TSV to Markdown with Miller.
%.tsv.md: %.tsv
	mlr --itsvlite --omd cat $< >$@
