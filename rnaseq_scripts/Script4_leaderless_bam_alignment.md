# Install packages locally (terminal)

Using the public transcriptome:
[SRR13362821](https://www.ncbi.nlm.nih.gov/sra/SRX9786900%5baccn)

And the matching genome for A. phalloides: 88mAP (10721)

    cd /Users/songs005/Documents/Programs
    curl --output sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-mac64.tar.gz
    tar -vxzf sratoolkit.tar.gz
    export PATH=$PWD/sratoolkit.3.2.1-mac-x86_64/bin:$PATH
    which fastq-dump
    # OUTPUT
    /Users/songs005/Documents/Programs/sratoolkit.3.2.1-mac-x86_64/bin/fastq-dump
    # test:
    fastq-dump --stdout -X 2 SRR390728
    # it produced the output!

    # set up sratoolkit config:
    vdb-config --interactive
    # local cache will save to:
    # /Users/songs005/Documents/Programs/sra_cache
    cd /Users/songs005/Documents/Programs/sra_cache
    prefetch SRR13362821 --max-size 45g # this took ~8 minutes to run on my laptop
    fasterq-dump --split-files SRR13362821
    vdb-validate SRR13362821

    # install hisat2
    cd /Users/songs005/Documents/Programs
    # conda
    conda
    conda install -c bioconda hisat2
    # check it
    hisat2 --version
    # /Users/songs005/miniconda3/bin/hisat2-align-s version 2.2.1
    # 64-bit
    # Built on VM-8b0b7afd3b
    # Wed Nov 27 22:00:39 GMT 2024
    # Compiler: clang-16: error: linker command failed with exit code 1 (use -v to see invocation)
    # Options: -O3   -funroll-loops -g3 -std=c++11
    # Sizeof {int, long, long long, void*, size_t, off_t}: {4, 8, 8, 8, 8, 8}

    # install fastqc
    conda install bioconda::fastqc
    fastqc --version
    # FastQC v0.12.1

    # install fastp
    conda install bioconda::fastp
    # got an error:
    # Channels:
    #  - bioconda
    #  - defaults
    # Platform: osx-arm64
    # Collecting package metadata (repodata.json): done
    # Solving environment: failed
    # 
    # LibMambaUnsatisfiableError: Encountered problems while solving:
    #   - nothing provides libdeflate >=1.20,<1.21.0a0 needed by fastp-0.23.4-h5510893_5
    # 
    # Could not solve for environment specs
    # The following packages are incompatible
    # └─ fastp =* * is not installable because there are no viable options
    #    ├─ fastp [0.23.4|0.24.0] would require
    #    │  └─ libdeflate >=1.20,<1.21.0a0 *, which does not exist (perhaps a missing channel);
    #    └─ fastp [0.24.0|0.24.1|...|1.0.1] would require
    #       └─ libcxx >=18 *, which does not exist (perhaps a missing channel).
    # tried this instead:
    conda config --add channels conda-forge
    conda config --set channel_priority strict
    conda install -c bioconda fastp
    fastp --version
    # fastp 0.24.0

    # install multiqc
    conda install bioconda::multiqc
    multiqc --version
    # multiqc, version 1.30

    # install bwa
    conda install bioconda::bwa
    bwa
    # Program: bwa (alignment via Burrows-Wheeler transformation)
    # Version: 0.7.18-r1243-dirty
    # Contact: Heng Li <hli@ds.dfci.harvard.edu>

    # install samtools using homebrew https://formulae.brew.sh/formula/samtools
    brew install samtools
    samtools --version
    # samtools 1.22.1
    # Using htslib 1.22.1
    # Copyright (C) 2025 Genome Research Ltd.

# Align public RNAseq data to genome

    cd /Users/songs005/Documents/Programs/sra_cache
    prefetch SRR13362821 --max-size 45g
    fasterq-dump --split-files SRR13362821
    # spots read      : 170,345,047
    # reads read      : 340,690,094
    # reads written   : 340,690,094
    vdb-validate SRR13362821
    # 2025-08-06T16:45:33 vdb-validate.3.2.1 info: Validating '/Users/songs005/Documents/Programs/sra_cache/sra/SRR13362821.sra'...
    # 2025-08-06T16:45:33 vdb-validate.3.2.1 info: Database 'SRR13362821.sra' metadata: md5 ok
    # 2025-08-06T16:45:33 vdb-validate.3.2.1 info: Table 'SEQUENCE' metadata: md5 ok
    # 2025-08-06T16:45:33 vdb-validate.3.2.1 info: Column 'ALTREAD': md5 ok
    # 2025-08-06T16:45:39 vdb-validate.3.2.1 info: Column 'QUALITY': md5 ok
    # 2025-08-06T16:45:58 vdb-validate.3.2.1 info: Column 'READ': md5 ok
    # 2025-08-06T16:45:58 vdb-validate.3.2.1 info: Column 'SPOT_GROUP': md5 ok
    # 2025-08-06T16:45:58 vdb-validate.3.2.1 info: Database '/Users/songs005/Documents/Programs/sra_cache/sra/SRR13362821.sra' contains only unaligned reads
    # 2025-08-06T16:45:58 vdb-validate.3.2.1 info: Database 'SRR13362821.sra' is consistent

    # now move files to aphal directory on my laptop
    cd /Volumes/Liv/Lab_projects/Amanita_phalloides/88mAP_rnaseq
    ls /Users/songs005/Documents/Programs/sra_cache
    mv /Users/songs005/Documents/Programs/sra_cache/SRR13362821_1.fastq .
    mv /Users/songs005/Documents/Programs/sra_cache/SRR13362821_2.fastq .

    # now run fastp - this took ~678 seconds or  ~11 min
    fastp -c -q 10 -l 50 --detect_adapter_for_pe \
      -i SRR13362821_1.fastq \
      -I SRR13362821_2.fastq \
      -o SRR13362821_1_trim.fastq\
      -O SRR13362821_2_trim.fastq \
      --html $dir/${samplename}_trim.fastp.html

    # now fastqc - this takes approx. ~ 10 min per file
    mkdir fastqc_output
    fastqc -o fastqc_output SRR13362821_1.fastq
    fastqc -o fastqc_output SRR13362821_2.fastq
    fastqc -o fastqc_output SRR13362821_1_trim.fastq
    fastqc -o fastqc_output SRR13362821_2_trim.fastq

    # now multiqc
    multiqc fastqc_output

    # next - hisat2 to align to genome
    # build hisat index - 35 seconds
    gunzip /Volumes/Liv/Lab_projects/Amanita_phalloides/aphalloides_reff/88mAP/88mAP_Trimmed_Scaffolds.fna.gz

    hisat2-build /Volumes/Liv/Lab_projects/Amanita_phalloides/aphalloides_reff/88mAP/88mAP_Trimmed_Scaffolds.fna 88mAP_aphal_hisat_index

    # align to 88 mAP reference - takes longer
    hisat2 -x 88mAP_aphal_hisat_index -1 SRR13362821_1_trim.fastq -2 SRR13362821_2_trim.fastq -S 88mAP.sam --phred33 --dta -t --rna-strandness RF

    # Time loading forward index: 00:00:01
    # Time loading reference: 00:00:00
    # Multiseed full-index search: 00:55:26
    # 161044471 reads; of these:
    #   161044471 (100.00%) were paired; of these:
    #     6405723 (3.98%) aligned concordantly 0 times
    #     150860000 (93.68%) aligned concordantly exactly 1 time
    #     3778748 (2.35%) aligned concordantly >1 times
    #     ----
    #     6405723 pairs aligned concordantly 0 times; of these:
    #       647638 (10.11%) aligned discordantly 1 time
    #     ----
    #     5758085 pairs aligned 0 times concordantly or discordantly; of these:
    #       11516170 mates make up the pairs; of these:
    #         9587000 (83.25%) aligned 0 times
    #         1553798 (13.49%) aligned exactly 1 time
    #         375372 (3.26%) aligned >1 times
    # 97.02% overall alignment rate
    # Time searching: 00:55:26
    # Overall time: 00:55:27

    # convert to bam
    samtools view -uS 88mAP.sam | samtools sort -o 88mAP.bam

    samtools index 88mAP.bam

    # also try using bwa mem - using the same code from Mickey's notes
    bwa index /Volumes/Liv/Lab_projects/Amanita_phalloides/aphalloides_reff/88mAP/88mAP_Trimmed_Scaffolds.fna

    bwa mem -t 4 -R "@RG\tID:mAF29\tSM:29\tPL:ILLUMINA" /Volumes/Liv/Lab_projects/Amanita_phalloides/aphalloides_reff/88mAP/88mAP_Trimmed_Scaffolds.fna SRR13362821_1_trim.fastq SRR13362821_2_trim.fastq > 88mAP_bwa.sam
    # Real time: 2398.883 sec; CPU: 9723.128 sec

    samtools view -uS 88mAP_bwa.sam | samtools sort -o 88mAP_bwa.bam

    samtools index 88mAP_bwa.bam

    # when I ran bwa mem, I noticed many lines said something like:
    # candidate unique pairs for (FF, FR, RF, RR): (2, 146301, 5, 0)

    # I therefore reran my hisat2 using FR orientation
    # align to 88 mAP reference - takes longer
    hisat2 -x 88mAP_aphal_hisat_index -1 SRR13362821_1_trim.fastq -2 SRR13362821_2_trim.fastq -S 88mAP_FR.sam --phred33 --dta -t --rna-strandness FR
    # 161044471 reads; of these:
    #   161044471 (100.00%) were paired; of these:
    #     6405723 (3.98%) aligned concordantly 0 times
    #     150860000 (93.68%) aligned concordantly exactly 1 time
    #     3778748 (2.35%) aligned concordantly >1 times
    #     ----
    #     6405723 pairs aligned concordantly 0 times; of these:
    #       647638 (10.11%) aligned discordantly 1 time
    #     ----
    #     5758085 pairs aligned 0 times concordantly or discordantly; of these:
    #       11516170 mates make up the pairs; of these:
    #         9587000 (83.25%) aligned 0 times
    #         1553798 (13.49%) aligned exactly 1 time
    #         375372 (3.26%) aligned >1 times
    # 97.02% overall alignment rate
    # Time searching: 00:53:27
    # Overall time: 00:53:28

    # convert to bam
    samtools view -uS 88mAP_FR.sam | samtools sort -o 88mAP_FR.bam

    samtools index 88mAP_FR.bam


    # remove all the sam files as they are very large.
    rm 88mAP_FR.sam
    rm 88mAP.sam
    rm 88mAP_bwa.sam
    rm 88mAP_aphal_hisat_index*

# Prepare 88mAP / 10721 mushroom genome in R

    # get chromosome
    chr <- "NODE_399_length_17301_cov_50.2263_ID_797"
    chr

    ## [1] "NODE_399_length_17301_cov_50.2263_ID_797"

    # define start and stop range
    start <- 9849
    stop <- 9982

    # get range
    range <- paste0(chr,":",start-350,"-",stop+150)
    range

    ## [1] "NODE_399_length_17301_cov_50.2263_ID_797:9499-10132"

    # load genome data
    library(Biostrings)

    ## Loading required package: BiocGenerics

    ## Loading required package: generics

    ## 
    ## Attaching package: 'generics'

    ## The following objects are masked from 'package:base':
    ## 
    ##     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
    ##     setequal, union

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
    ##     unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## Loading required package: XVector

    ## Loading required package: GenomeInfoDb

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

    # read the FASTA file
    fasta <- readDNAStringSet("/Volumes/Liv/Lab_projects/Amanita_phalloides/aphalloides_reff/88mAP/88mAP_Trimmed_Scaffolds.fna")

    # get chromosome/contig names
    names <- names(fasta)

    # get sequence lengths
    lengths <- width(fasta)

    # combine into a data frame that can be used in karyoploteR
    # see tutorial - https://bernatgel.github.io/karyoploter_tutorial//Tutorial/CustomGenomes/CustomGenomes.html
    library(karyoploteR)

    ## Loading required package: regioneR

    ## Loading required package: GenomicRanges

    custom.genome <- data.frame(chr = names, start = 1, end = lengths)

    # view length for chromosome of interest
    length <- custom.genome[custom.genome$chr == chr,"end"]
    length

    ## [1] 17301

# Panels A and B: RNA read density & chromosome diagram

    # load alignment file
    library(GenomicAlignments)

    ## Loading required package: SummarizedExperiment

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'MatrixGenerics'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## 
    ## Attaching package: 'Biobase'

    ## The following object is masked from 'package:MatrixGenerics':
    ## 
    ##     rowMedians

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     anyMissing, rowMedians

    ## Loading required package: Rsamtools

    library(Rsamtools)

    # custom BAM file (must be indexed, i.e., .bai file present)
    # first - bwa mem bam file
    bam <- "/Volumes/Liv/Lab_projects/Amanita_phalloides/88mAP_rnaseq/88mAP_bwa.bam"
    # then - hisat2 FR file
    bam2 <- "/Volumes/Liv/Lab_projects/Amanita_phalloides/88mAP_rnaseq/88mAP_FR.bam"
    # finally - hisat2 RF file
    bam3 <- "/Volumes/Liv/Lab_projects/Amanita_phalloides/88mAP_rnaseq/88mAP.bam"

    # plot forward
    pdf("/Volumes/Liv/Lab_projects/Amanita_phalloides/88mAP_rnaseq/bam_plot.pdf", width = 7.08661, height = 3, family="ArialMT")  # Specify desired width and height

    kp <- plotKaryotype(genome = custom.genome, chromosomes = chr, zoom=toGRanges(range))
    kpAddBaseNumbers(kp, tick.dist = 200, add.units = TRUE)
    kp <- kpPlotBAMCoverage(kp, data=bam3)
    kpAxis(kp, ymax=kp$latest.plot$computed.values$max.coverage)
    kpRect(kp,
           chr = chr,
           x0 = start, x1 = stop,
           y0 = 0, y1 = 1,
           col = NA, border = "red", lwd = 2)
    dev.off()

    ## quartz_off_screen 
    ##                 2

# Notes on the panel A - introns

When I inspect this bam in IGV, I can see two regions where there are
virtually no aligned reads and are likely therefore introns: ![IGV
screenshot](/Users/songs005/Documents/igv.png) For some reason unknown
to me, the plotKaryotype() and kpPlotBAMCoverage() functions are not
showing the low read density at the introns: ![kpPlotBAMCoverage
output](/Volumes/Liv/Lab_projects/Amanita_phalloides/88mAP_rnaseq/bam_plot.pdf)
I therefore opened this plot in illustrator and manually indicated the
introns for the figure in the manuscript. I also flipped the plot
horizontally, since the reading frames we are interested in are in
reverse, which is harder to read.

# Panel C: Amino acid translations

    # prepare amino acid translation
    # MLGFLVLP is on the reverse strand; the potential previous mutated leader sequence starts 16 aa upstream # so add 48 nt downstream of the stop codon and translate entire region
    # also want to include upstream region until it reaches the nearest stop codon
    region <- GRanges(seqnames = chr, ranges = IRanges(start = start, end = stop+156))
    seq <- reverseComplement(getSeq(fasta, region))
    seq

    ## DNAStringSet object of length 1:
    ##     width seq
    ## [1]   290 GTCAGTAGATTGTGCTATGGTACGATTTCCCTCC...GACTGATGTGTTATGCGTTAGCCTTTGTTAAAT

    # translate to amino acids (3 forward frames)
    aa1 <- translate(seq, if.fuzzy.codon="solve")

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : last 2 bases were ignored

    aa2 <- translate(subseq(seq, start=2), if.fuzzy.codon="solve")

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : last base was ignored

    aa3 <- translate(subseq(seq, start=3), if.fuzzy.codon="solve")
    # Give each frame a name
    names(aa1) <- "+1"
    names(aa2) <- "+2"
    names(aa3) <- "+3"

    # check the seqs
    aa1

    ## AAStringSet object of length 1:
    ##     width seq                                               names               
    ## [1]    96 VSRLCYGTISLRNQVISHEWDG*...RGEPTTHPTTYPWTDVLCVSLC* +1

    aa2

    ## AAStringSet object of length 1:
    ##     width seq                                               names               
    ## [1]    96 SVDCAMVRFPSAIRS*ATSGTDK...EVSPPPIQRRIHGLMCYALAFVK +2

    aa3

    ## AAStringSet object of length 1:
    ##     width seq                                               names               
    ## [1]    96 Q*IVLWYDFPPQSGHKPRVGRIK...R*AHHPSNDVSMD*CVMR*PLLN +3

    library(ggmsa)

    ## Registered S3 methods overwritten by 'ggalt':
    ##   method                  from   
    ##   grid.draw.absoluteGrob  ggplot2
    ##   grobHeight.absoluteGrob ggplot2
    ##   grobWidth.absoluteGrob  ggplot2
    ##   grobX.absoluteGrob      ggplot2
    ##   grobY.absoluteGrob      ggplot2

    ## ggmsa v1.14.0  Document: http://yulab-smu.top/ggmsa/
    ## 
    ## If you use ggmsa in published research, please cite:
    ## L Zhou, T Feng, S Xu, F Gao, TT Lam, Q Wang, T Wu, H Huang, L Zhan, L Li, Y Guan, Z Dai*, G Yu* ggmsa: a visual exploration tool for multiple sequence alignment and associated data. Briefings in Bioinformatics. DOI:10.1093/bib/bbac222

    # combine into one AAStringSet
    # exon 1 is in +1 and exon 2 is in +2
    aaset <- AAStringSet(c(aa1, aa2, aa3))

    # make custom colors
    my_custom <- data.frame(names = c(LETTERS[1:26],"-"), 
                             color = "gray80", 
                             stringsAsFactors = FALSE)
    head(my_custom)

    ##   names  color
    ## 1     A gray80
    ## 2     B gray80
    ## 3     C gray80
    ## 4     D gray80
    ## 5     E gray80
    ## 6     F gray80

    # Override M color
    my_custom[my_custom$names == "M", "color"] <- "#03FF00"

    #ggmsa(aaset, char_width = 0.5, seq_name = T)
    library(ggplot2)
    ggmsa(aaset, char_width = 0.5, seq_name = T, start = 1, custom_color = my_custom) +
         geom_rect(mapping = aes(xmin = 51.5, xmax = 73.5,ymin = 0.5, ymax = 1.5), # LABEL EXON 1
                  alpha=0.2, 
                  color = "black", 
                  size = 1, 
                  fill = NA) +
       geom_rect(mapping = aes(xmin = 92.5, xmax = 96.5,ymin = 2.5, ymax = 3.5), # LABEL EXON 2
                  alpha=0.2, 
                  color = "black", 
                  size = 1, 
                  fill = NA) +
         geom_rect(mapping = aes(xmin = 35.5, xmax = 51.5,ymin = 0.5, ymax = 1.5), # LABEL MUTATED LEADER
                  alpha=0.2, 
                  color = "#005AB5", 
                  size = 1, 
                  fill = NA)  +
         geom_rect(mapping = aes(xmin = 51.5, xmax = 59.5,ymin = 0.5, ymax = 1.5), # LABEL CORE SEQ
                  alpha=0.2, 
                  color = "#D81B60", 
                  size = 1, 
                  fill = NA)

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(dd$x): no non-missing arguments to min; returning Inf

    ## Warning in min(dd$y): no non-missing arguments to min; returning Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(dd$x): no non-missing arguments to min; returning Inf

    ## Warning in min(dd$y): no non-missing arguments to min; returning Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(dd$x): no non-missing arguments to min; returning Inf

    ## Warning in min(dd$y): no non-missing arguments to min; returning Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(dd$x): no non-missing arguments to min; returning Inf

    ## Warning in min(dd$y): no non-missing arguments to min; returning Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(dd$x): no non-missing arguments to min; returning Inf

    ## Warning in min(dd$y): no non-missing arguments to min; returning Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(dd$x): no non-missing arguments to min; returning Inf

    ## Warning in min(dd$y): no non-missing arguments to min; returning Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(dd$x): no non-missing arguments to min; returning Inf

    ## Warning in min(dd$y): no non-missing arguments to min; returning Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(dd$x): no non-missing arguments to min; returning Inf

    ## Warning in min(dd$y): no non-missing arguments to min; returning Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(dd$x): no non-missing arguments to min; returning Inf

    ## Warning in min(dd$y): no non-missing arguments to min; returning Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(dd$x): no non-missing arguments to min; returning Inf

    ## Warning in min(dd$y): no non-missing arguments to min; returning Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(dd$x): no non-missing arguments to min; returning Inf

    ## Warning in min(dd$y): no non-missing arguments to min; returning Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(x, na.rm = na.rm): no non-missing arguments to min; returning
    ## Inf

    ## Warning in max(x, na.rm = na.rm): no non-missing arguments to max; returning
    ## -Inf

    ## Warning in min(dd$x): no non-missing arguments to min; returning Inf

    ## Warning in min(dd$y): no non-missing arguments to min; returning Inf

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## Warning: No shared levels found between `names(values)` of the manual scale and the
    ## data's fill values.

![](Script4_leaderless_bam_alignment_files/figure-markdown_strict/amino_acids-1.png)

    ggsave("/Volumes/Liv/Lab_projects/Amanita_phalloides/88mAP_rnaseq/aa_translation.pdf",width=7.08661,height=1) 

    ## Warning: No shared levels found between `names(values)` of the manual scale and the
    ## data's fill values.
