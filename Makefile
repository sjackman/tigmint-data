# Correct misassemblies using Tigmint

# Number of threads
t=64

# gzip compression program
gzip=pigz -p$t

.DELETE_ON_ERROR:
.SECONDARY:

all: abyss2 abyss2_bionano_arcs

abyss2:
	$(MAKE) draft=$@ \
		abyss2.abyss-fac.tsv \
		abyss2.scaftigs.abyss-fac.tsv \
		abyss2.scaftigs.GRCh38.samtobreak.tsv \
		abyss2.hg004.bx.as100.nm5.bam.bai \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.tsv \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.bed.bam.bai \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.summary.html \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.depth.stats.tsv \
		abyss2.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.breakpoints.tsv

abyss2_bionano_arcs:
	$(MAKE) draft=$@ \
		abyss2_bionano_arcs.abyss-fac.tsv \
		abyss2_bionano_arcs.scaftigs.abyss-fac.tsv \
		abyss2_bionano_arcs.scaftigs.GRCh38.samtobreak.tsv \
		abyss2_bionano_arcs.hg004.bx.as100.nm5.bam.bai \
		abyss2_bionano_arcs.hg004.bx.as100.nm5.bam.mi.bx.molecule.tsv \
		abyss2_bionano_arcs.hg004.bx.as100.nm5.bam.mi.bx.molecule.bed.bam.bai \
		abyss2_bionano_arcs.hg004.bx.as100.nm5.bam.mi.bx.molecule.summary.html \
		abyss2_bionano_arcs.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.bed.depth.stats.tsv \
		abyss2_bionano_arcs.hg004.bx.as100.nm5.bam.mi.bx.molecule.size2000.depth.100.starts.breakpoints.tsv

# Rename the LINKS scaffolds.
abyss2_bionano_arcs.fa: %.fa: %.orig.fa
	gsed 's/scaffold//;s/,[^\t]*//' $< >$@

# Add the barcode to the read ID, and skip reads without barcodes.
%.bx.fq.gz: %.longranger.basic.fq.gz
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

# Align paired-end reads to the draft genome and sort.
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
starts_threshold=4
pos_threshold=1000
%.depth.$(depth_threshold).starts.breakpoints.tsv: %.depth.$(depth_threshold).starts.tsv
	mlr --tsvlite filter '$$Depth < $(depth_threshold) && $$Starts >= $(starts_threshold) && $$Pos >= $(pos_threshold)' $< >$@

# Identify breakpoints with low depth of coverage and high number of molecule starts.
%.size$(size_threshold).depth.starts.breakpoints.tsv: %.size$(size_threshold).bed.depth.tsv %.starts.tsv
	Rscript -e 'rmarkdown::render("breakpoints.rmd", "html_notebook", "$*.depth.starts.breakpoints.nb.html", params = list(depth_tsv="$<", starts_tsv="$*.starts.tsv", depth_starts_tsv="$*.depth.starts.tsv", breakpoints_tsv="$@"))'

################################################################################
# Calculate assembly contiguity and correctness metrics.

# Reference genome
ref=GRCh38
ref_fa=$(ref).fa
ref_gff=/projects/btl/reference_genomes/H_sapiens/GRCh38/Homo_sapiens.GRCh38.86.chr.gff3

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

# ABySS

# Convert scaffolds to scaftigs.
%.scaftigs.fa: %.fa
	seqtk seq $< | tr _ '~' | $(abyss_bin)/abyss-fatoagp -f $@ >$@.agp

# Calculate assembly contiguity metrics with abyss-fac.
%.abyss-fac.tsv: %.fa
	$(abyss_bin)/abyss-fac -G$(GwithN) -t500 $< >$@

# Calculate assembly contiguity and correctness metrics.
%.samtobreak.txt: %.sam
	(echo '==> $< <=='; bin/abyss-samtobreak -G$(GwithN) -l500 $<) >$@

# Convert samtobreak.txt to TSV.
%.samtobreak.tsv: %.samtobreak.txt
	bin/abyss-samtobreak-to-tsv $< >$@
