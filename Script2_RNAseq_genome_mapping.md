### Data QC

First I uploaded the data to MSI from my laptop to filezilla. I am using
the reference genome 10511_Aphal_PT_AllpathsLG_LINKS_jelly_pilon.fna,
which is from the 72mAP mushroom and assembled from long-read
sequencing.

## Rough pipeline

On Minnesota supercomputer (MSI):

-   fastqc
-   fastp trimming
-   multiqc
-   hisat
-   featurecounts

Then transfer to laptop for DESeq2 in R.

I will use slurm arrays for fast processing on MSI. First I will
generate a config file that specifies my samples, based on the fastq
reads:

``` r
# get filenames
filenames <- list.files("/Users/songs005/Library/CloudStorage/Box-Box/Drott_lab_shared_backup/Data/Aphalloides_rnaseq_for_leaderless/RAWTRANS/raw_data/201016_AHL2WCDSXY",pattern=".gz")

# extract sample ID column without underscores
Sample_ID <- unique(paste0(sapply(strsplit(filenames,split="_"), `[`, 1)))

SampleName <- unique(paste0(sapply(strsplit(filenames,split="_"), `[`, 1),
                            "_",
                            sapply(strsplit(filenames,split="_"), `[`, 2),
                            "_",
                            sapply(strsplit(filenames,split="_"), `[`, 3)))


# make empty config file
config <- as.data.frame(cbind(rep(1:length(SampleName)), SampleName,Sample_ID))

# add dir column to specify output dirs
colnames(config)[1] <- "ArrayTaskID"

# save output
write.table(config,"slurm_scripts/aphal_rnaseq_config.txt",row.names=FALSE)

# remember to transfer this config file to MSI and specify its filepath in the following slurm scripts.
```

### set up qc environment

Ensure config and data files are transferred to MSI. Log onto MSI. Set
up environments & packages for analysis.

``` bash
module load conda
conda --version
conda install -c bioconda fastp
conda install -c bioconda fastqc
conda install -c bioconda multiqc

#Next create an environment and install tools the tools
conda create --yes -n qc fastp fastqc multiqc

#To activate the env, in an interactive terminal you would type:
conda activate qc

# if you want to activate the environment in a bash script for slurm stuff, type:
source activate qc

# check versions!
fastp -version

# if -version fails, try -v or --version or -help or --help -h
fastqc -version
multiqc -version
```

#### run_fastqc.sh

``` bash
#!/bin/bash -l
#SBATCH --job-name=fastqc_raw_trim    # Optional job name
#SBATCH --ntasks=1                      # Total number of tasks
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g                     # max memory limit
#SBATCH --time=1:00:00                # Walltime limit
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=songs005@umn.edu
#SBATCH --array=1-39                  # this will run one job for each sample

## This pulls from a tab-delimited file (array_config.txt) which has two columns. 
## First column is labeled ArrayTaskID, second is Sample as the header and first row
## example:
## head(config.txt)
##  ArrayTaskID     SampleName Sample_ID
## 1           1 10745_S87_L001     10745

# define config file & variables
config=/scratch.global/songs005/aphal_rnaseq/aphal_rnaseq_config.txt

# sample names / prefixes for reads
samplename=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

# location of raw reads
data=/scratch.global/songs005/aphal_rnaseq/raw_reads

# output dir
dir=/scratch.global/songs005/aphal_rnaseq/raw_reads/fastqc

# make output directory if it doesn't exist already
mkdir -p $dir

# then run QC on the raw reads
fastqc -o $dir ${data}/${samplename}_R1_001.fastq.gz
fastqc -o $dir ${data}/${samplename}_R2_002.fastq.gz
```

#### run_fastp.sh

This runs very fast - 3-5 min per read file typically

``` bash
#!/bin/bash -l
#SBATCH --job-name=fastp_trim         # Optional job name
#SBATCH --ntasks=1                      # Total number of tasks
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g                     # max memory limit
#SBATCH --time=0:20:00                # Walltime limit
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=songs005@umn.edu
#SBATCH --array=1-39                  # this will run one job for each sample

# define config file & variables
config=/scratch.global/songs005/aphal_rnaseq/aphal_rnaseq_config.txt

# sample names / prefixes for reads
samplename=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

# location of raw reads
data=/scratch.global/songs005/aphal_rnaseq/raw_reads

# output dir
dir=/scratch.global/songs005/aphal_rnaseq/trimmed_data

# make output directory
mkdir -p $dir

# run fastp; this will autodetect adapters 
# -c is for base correction based on R1/R2 overlap
# -q is for qualified_quality_phred; minimal acceptable quality value
# -l is for length filtering options; reads shorter than this length will be discarded
# --detect_adapter_for_pe will specifically look for adapters; can result in a slightly more accurate & cleaner trim; slows down script a bit
  
fastp -c -q 10 -l 50 --detect_adapter_for_pe \
  -i $data/${samplename}_R1_001.fastq.gz \
  -I $data/${samplename}_R2_001.fastq.gz \
  -o $dir/${samplename}_R1_trim.fastq.gz \
  -O $dir/${samplename}_R2_trim.fastq.gz \
  --html $dir/${samplename}_trim.fastp.html \
  --json $dir/${samplename}_trim.fastp.json 
  
```

#### run_trim_qc.sh

``` bash
#!/bin/bash -l
#SBATCH --job-name=fastqc_trim        # Optional job name
#SBATCH --ntasks=1                      # Total number of tasks
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g                     # max memory limit
#SBATCH --time=0:40:00                # Walltime limit
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=songs005@umn.edu
#SBATCH --array=1-39                  # this will run one job for each sample

# define config file & variables
config=/scratch.global/songs005/aphal_rnaseq/aphal_rnaseq_config.txt

# sample names / prefixes for reads
samplename=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

# location of reads
data=/scratch.global/songs005/aphal_rnaseq/trimmed_data

# output dir
dir=/scratch.global/songs005/aphal_rnaseq/trimmed_data/fastqc

# make output directory
mkdir -p $dir

# then run QC on the raw reads
fastqc -o $dir ${data}/${samplename}_R1_trim.fastq.gz
fastqc -o $dir ${data}/${samplename}_R2_trim.fastq.gz
```

After running fastq on the raw and trimmed reads, I ran multiqc on each
directory to summarize fastQC results. The command for this is simply:

``` bash
# load conda
module load conda
conda activate qc
# cd to dir with fastqc
cd /scratch.global/songs005/aphal_rnaseq/trimmed_data/fastqc
multiqc .

cd /scratch.global/songs005/aphal_rnaseq/raw_data/fastqc
multiqc .
```

I transferred the QC results to my laptop using filezilla.

### Hisat2

This took about 40 minutes to run hisat2 for one sample, and another 40
minutes or so to convert the sam to bam and index it.

``` bash
# start interactive session
# this is to troubleshoot the code and see what works
srun -N 1 --ntasks-per-node=4 --mem=60gb --gres=gpu:a40:1 -t 2:00:00 -p interactive-gpu --tmp 50gb --pty bash

# load hisat2 and samtools
module load hisat2
hisat2 --version
# OUTPUT:
# /common/software/install/migrated/hisat2/2.1.0/bin/hisat2-align-s version 2.1.0
# 64-bit
# Built on swan.msi.umn.edu
# Tue Sep 19 14:45:21 CDT 2017
# Compiler: gcc version 4.4.7 20120313 (Red Hat 4.4.7-18) (GCC)

module load samtools
samtools --version
# OUTPUT:
# samtools 1.16.1
# Using htslib 1.16

# generate index - takes like 30 seconds
hisat2-build reff/10511_Aphal_PT_AllpathsLG_LINKS_jelly_pilon.fna reff/aphal_hisat_index

# test run on one sample - took about 40 min
hisat2 -x reff/aphal_hisat_index -1 trimmed_data/10708_S250_L002_trim_R1.fastq.gz -2 trimmed_data/10708_S250_L002_trim_R2.fastq.gz -S test.sam --phred33 --dta -t --rna-strandness RF

# convert to bam - took about 40 min
samtools view -uS test.sam | samtools sort -o test.bam
samtools index test.bam
```

#### hisat2.sh

``` bash
#!/bin/bash -l
#SBATCH --job-name=hisat2         # Optional job name
#SBATCH --ntasks=1                # Total number of tasks
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00            # Walltime limit
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=songs005@umn.edu
#SBATCH --array=1-39              # this will run one job for each sample

# config file
config=/users/2/songs005/fastq_species_detector/aphal_config.txt
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

# set dirs
genomedir="/scratch.global/songs005/aphal_rnaseq/reff"
fastqdir="/home/mdrott/shared/liv/trimmed_data"
outdir="/scratch.global/songs005/aphal_rnaseq/17june_hisat2" # with hisat2 output
outdir2="/scratch.global/songs005/aphal_rnaseq/17june_hisat2_unmapped" # with hisat2 output

# make a specific dir for this job
#mkdir -p $outdir
mkdir -p $outdir/$sample

# load modules
module load hisat2
module load samtools

# build index - only needs to be done done once, can comment out
#hisat2-build $genomedir/10511_Aphal_PT_AllpathsLG_LINKS_jelly_pilon.fna $genomedir/aphal_hisat_index

# run hisat2
hisat2 -x $genomedir/aphal_hisat_index \
-1 $fastqdir/${sample}_trim_R1.fastq.gz \
-2 $fastqdir/${sample}_trim_R2.fastq.gz \
-S $outdir/$sample/${sample}.sam \
--phred33 \
--rna-strandness FR \
--dta \
--max-intronlen 3000 \
--un-conc $outdir/$sample \
-t # print walltime

# convert sam to bam format
samtools view -uS $outdir/$sample/${sample}.sam | samtools sort -o $outdir/$sample/${sample}.bam
samtools index $outdir/$sample/${sample}.bam
```

### Featurecounts

#### Subread installation

Next we run featurecounts from the subread package. To use
featurecounts, I must install the subread package on MSI.

``` bash
# exit qc environment
conda deactivate

# make new environment
conda create --name subread_env
conda activate subread_env

# install subread
conda install bioconda::subread
conda list
# OUTPUT BELOW:
# packages in environment at /users/2/songs005/.conda/envs/subread_env:
#
# Name                    Version                   Build  Channel
# _libgcc_mutex             0.1                        main
# _openmp_mutex             5.1                       1_gnu
# libgcc-ng                 11.2.0               h1234567_1
# libgomp                   11.2.0               h1234567_1
# subread                   2.0.1                h5bf99c6_1    bioconda
# zlib                      1.2.13               h5eee18b_1
```

#### featurecounts.sh

This script will only run after running the following on the console:

``` bash
module load conda
conda activate subread_env
```

bash script below:

``` bash
#!/bin/bash -l
#SBATCH --job-name=featurecounts            # Optional job name
#SBATCH --ntasks=1                          # Total number of tasks
#SBATCH --cpus-per-task=8
#SBATCH --mem=64g                           # max memory limit
#SBATCH --time=2:00:00                      # Walltime limit
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=songs005@umn.edu
#SBATCH --array=1-39

config=/scratch.global/songs005/aphal_rnaseq/reff/config.txt
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

# set dirs
genomedir="/scratch.global/songs005/aphal_rnaseq/reff"
fastqdir="/scratch.global/songs005/aphal_rnaseq/trimmed_data"
outdir="/scratch.global/songs005/aphal_rnaseq/hisat2_output_msdins_array" # with star output
countdir="/scratch.global/songs005/aphal_rnaseq/hisat2_featurecounts_output_msdins" # for featureCounts output

# mkdir
mkdir -p $countdir
mkdir -p $countdir/$sample

featureCounts -T 20 \
        -a $genomedir/aphal_reff_with_msdins.gff3 \
        -o $countdir/$sample/${sample}_featurecounts.txt \
        --minOverlap 20 \
        -g Parent \
        -B \
        -C \
        -p \
        $outdir/$sample/${sample}.bam

# arguments for featureCounts
# -a genome annotation file
# -o output filename
# -g parses the attribute 'Parent' in the gff3 file
# -b include only fragments with both ends successfully aligned
# -C do not count chimeric fragments where two ends align to diff chr
# -p specify paired-end reads
# input file is listed last with no flag
```

### Fastqc on unmapped reads

Since Hisat2 output a file containing the unmapped reads for each
sample, I ran fastqc on these as well to see if there were any sequence
features that popped out.

``` bash
#!/bin/bash -l
#SBATCH --job-name=1unmapped_fastqc         # Optional job name
#SBATCH --ntasks=1                          # Total number of tasks
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g                           # max memory limit
#SBATCH --time=1:00:00                      # Walltime limit
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --mail-user=songs005@umn.edu
#SBATCH --array=1-39
config=/scratch.global/songs005/aphal_rnaseq/reff/unmapped_aphal_config.txt
sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

# source activate qc

# location of unmapped reads

data=/scratch.global/songs005/aphal_rnaseq/hisat2_output_msdins_array_unmapp

# output dir
dir=/scratch.global/songs005/aphal_rnaseq/unmapped_fastqc

# make output directory
mkdir -p $dir

# then run QC on the raw reads
fastqc -o $dir/${sample} $data/*.1.fastq
fastqc -o $dir/${sample} $data/*.2.fastq
```

Next I switched to R to extract the overrepresented sequences from each
unmapped file.

Fastqc automatically identifies duplicate sequences within each fasta,
so I can extract these from the fastqc_data.txt files and then blast
them or align them to a reference sequence.

``` r
unmapp_path = "/Volumes/Liv/Lab_projects/Amanita_phalloides/MSI_outputs/unmapped_fastqc"

list_unmapped <- list.files(path=unmapp_path,recursive=TRUE,pattern="fastqc_data.txt",full.names=TRUE)

read_fastqc <- function(filename) {
  #filename=list_unmapped[1]
  filename_simp <- strsplit(filename, "/")[[1]][9]
  temp <- readLines(filename)
  
  # find row with overrepresented seqs
  
  temp <- as.data.frame(Map(function(i,j) read.table(text=temp[(i+1):(j-2)], sep='', header=FALSE), grep('>>Overrepresented', temp), grep('>>Adapter', temp)))
  
  # rename columns
  colnames(temp) <- c("Sequence","Count","Percentage","Possible_source")
  temp$header <- paste0(">",rownames(temp))
  
  # extract seqs as a fasta
  sequence_list <- as.vector(rbind(temp$header, temp$Sequence))
  
  # save temp as output in a folder
  outdir = "/Volumes/Liv/Lab_projects/Amanita_phalloides/MSI_outputs/unmapped_fastqc/overrep_stats/"

  outdir2 = "/Volumes/Liv/Lab_projects/Amanita_phalloides/MSI_outputs/unmapped_fastqc/overrep_seqs/"
  
  write.csv(temp, paste0(outdir,filename_simp,"_fastqc_overrep_stats.csv"),row.names=FALSE)
  write.table(sequence_list, 
            file = paste0(outdir2, filename_simp, "_fastqc_overrep_seqs.fasta"), 
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE)
  print(paste0("Extracted overrepresented sequences for: ",filename_simp))
  
}

# apply function to all entries on list
lapply(list_unmapped, read_fastqc)
```

    ## [1] "Extracted overrepresented sequences for: 10708_S250_L002_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10708_S250_L002_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10709_S251_L002_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10709_S251_L002_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10710_S252_L002_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10710_S252_L002_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10711_S253_L002_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10711_S253_L002_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10713_S254_L002_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10713_S254_L002_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10716_S255_L002_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10716_S255_L002_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10717_S256_L002_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10717_S256_L002_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10720_S257_L002_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10720_S257_L002_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10721_S279_L004_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10721_S279_L004_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10745_S280_L004_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10745_S280_L004_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10745_S87_L001_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10745_S87_L001_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10749_S281_L004_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10749_S281_L004_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10751_S282_L004_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10751_S282_L004_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10754_S283_L004_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10754_S283_L004_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10755_S284_L004_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10755_S284_L004_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10763_S285_L004_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10763_S285_L004_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10764_S286_L004_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10764_S286_L004_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10765_S258_L002_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10765_S258_L002_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10767_S259_L002_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10767_S259_L002_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10769_S260_L002_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10769_S260_L002_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10770_S261_L002_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: 10770_S261_L002_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: Aso1_S262_L002_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: Aso1_S262_L002_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: Champ2_S263_L002_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: Champ2_S263_L002_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: Doc1_S264_L002_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: Doc1_S264_L002_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: Doc2_S265_L002_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: Doc2_S265_L002_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: Doc3_S266_L002_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: Doc3_S266_L002_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: Doc4_S287_L004_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: Doc4_S287_L004_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: Gron7_S267_L002_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: Gron7_S267_L002_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: Gron8_S268_L002_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: Gron8_S268_L002_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: Gron9_S269_L002_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: Gron9_S269_L002_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: Lamp1_S270_L002_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: Lamp1_S270_L002_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: Lamp2_S271_L002_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: Lamp2_S271_L002_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: Lamp4_S272_L002_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: Lamp4_S272_L002_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: Lamp5_S273_L002_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: Lamp5_S273_L002_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: Lamp6_S274_L002_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: Lamp6_S274_L002_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: Lamp8_S275_L002_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: Lamp8_S275_L002_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: Quill1_S276_L002_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: Quill1_S276_L002_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: Tisza2_S277_L002_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: Tisza2_S277_L002_unmapped.2_fastqc"
    ## [1] "Extracted overrepresented sequences for: Tisza3_S278_L002_unmapped.1_fastqc"
    ## [1] "Extracted overrepresented sequences for: Tisza3_S278_L002_unmapped.2_fastqc"

    ## [[1]]
    ## [1] "Extracted overrepresented sequences for: 10708_S250_L002_unmapped.1_fastqc"
    ## 
    ## [[2]]
    ## [1] "Extracted overrepresented sequences for: 10708_S250_L002_unmapped.2_fastqc"
    ## 
    ## [[3]]
    ## [1] "Extracted overrepresented sequences for: 10709_S251_L002_unmapped.1_fastqc"
    ## 
    ## [[4]]
    ## [1] "Extracted overrepresented sequences for: 10709_S251_L002_unmapped.2_fastqc"
    ## 
    ## [[5]]
    ## [1] "Extracted overrepresented sequences for: 10710_S252_L002_unmapped.1_fastqc"
    ## 
    ## [[6]]
    ## [1] "Extracted overrepresented sequences for: 10710_S252_L002_unmapped.2_fastqc"
    ## 
    ## [[7]]
    ## [1] "Extracted overrepresented sequences for: 10711_S253_L002_unmapped.1_fastqc"
    ## 
    ## [[8]]
    ## [1] "Extracted overrepresented sequences for: 10711_S253_L002_unmapped.2_fastqc"
    ## 
    ## [[9]]
    ## [1] "Extracted overrepresented sequences for: 10713_S254_L002_unmapped.1_fastqc"
    ## 
    ## [[10]]
    ## [1] "Extracted overrepresented sequences for: 10713_S254_L002_unmapped.2_fastqc"
    ## 
    ## [[11]]
    ## [1] "Extracted overrepresented sequences for: 10716_S255_L002_unmapped.1_fastqc"
    ## 
    ## [[12]]
    ## [1] "Extracted overrepresented sequences for: 10716_S255_L002_unmapped.2_fastqc"
    ## 
    ## [[13]]
    ## [1] "Extracted overrepresented sequences for: 10717_S256_L002_unmapped.1_fastqc"
    ## 
    ## [[14]]
    ## [1] "Extracted overrepresented sequences for: 10717_S256_L002_unmapped.2_fastqc"
    ## 
    ## [[15]]
    ## [1] "Extracted overrepresented sequences for: 10720_S257_L002_unmapped.1_fastqc"
    ## 
    ## [[16]]
    ## [1] "Extracted overrepresented sequences for: 10720_S257_L002_unmapped.2_fastqc"
    ## 
    ## [[17]]
    ## [1] "Extracted overrepresented sequences for: 10721_S279_L004_unmapped.1_fastqc"
    ## 
    ## [[18]]
    ## [1] "Extracted overrepresented sequences for: 10721_S279_L004_unmapped.2_fastqc"
    ## 
    ## [[19]]
    ## [1] "Extracted overrepresented sequences for: 10745_S280_L004_unmapped.1_fastqc"
    ## 
    ## [[20]]
    ## [1] "Extracted overrepresented sequences for: 10745_S280_L004_unmapped.2_fastqc"
    ## 
    ## [[21]]
    ## [1] "Extracted overrepresented sequences for: 10745_S87_L001_unmapped.1_fastqc"
    ## 
    ## [[22]]
    ## [1] "Extracted overrepresented sequences for: 10745_S87_L001_unmapped.2_fastqc"
    ## 
    ## [[23]]
    ## [1] "Extracted overrepresented sequences for: 10749_S281_L004_unmapped.1_fastqc"
    ## 
    ## [[24]]
    ## [1] "Extracted overrepresented sequences for: 10749_S281_L004_unmapped.2_fastqc"
    ## 
    ## [[25]]
    ## [1] "Extracted overrepresented sequences for: 10751_S282_L004_unmapped.1_fastqc"
    ## 
    ## [[26]]
    ## [1] "Extracted overrepresented sequences for: 10751_S282_L004_unmapped.2_fastqc"
    ## 
    ## [[27]]
    ## [1] "Extracted overrepresented sequences for: 10754_S283_L004_unmapped.1_fastqc"
    ## 
    ## [[28]]
    ## [1] "Extracted overrepresented sequences for: 10754_S283_L004_unmapped.2_fastqc"
    ## 
    ## [[29]]
    ## [1] "Extracted overrepresented sequences for: 10755_S284_L004_unmapped.1_fastqc"
    ## 
    ## [[30]]
    ## [1] "Extracted overrepresented sequences for: 10755_S284_L004_unmapped.2_fastqc"
    ## 
    ## [[31]]
    ## [1] "Extracted overrepresented sequences for: 10763_S285_L004_unmapped.1_fastqc"
    ## 
    ## [[32]]
    ## [1] "Extracted overrepresented sequences for: 10763_S285_L004_unmapped.2_fastqc"
    ## 
    ## [[33]]
    ## [1] "Extracted overrepresented sequences for: 10764_S286_L004_unmapped.1_fastqc"
    ## 
    ## [[34]]
    ## [1] "Extracted overrepresented sequences for: 10764_S286_L004_unmapped.2_fastqc"
    ## 
    ## [[35]]
    ## [1] "Extracted overrepresented sequences for: 10765_S258_L002_unmapped.1_fastqc"
    ## 
    ## [[36]]
    ## [1] "Extracted overrepresented sequences for: 10765_S258_L002_unmapped.2_fastqc"
    ## 
    ## [[37]]
    ## [1] "Extracted overrepresented sequences for: 10767_S259_L002_unmapped.1_fastqc"
    ## 
    ## [[38]]
    ## [1] "Extracted overrepresented sequences for: 10767_S259_L002_unmapped.2_fastqc"
    ## 
    ## [[39]]
    ## [1] "Extracted overrepresented sequences for: 10769_S260_L002_unmapped.1_fastqc"
    ## 
    ## [[40]]
    ## [1] "Extracted overrepresented sequences for: 10769_S260_L002_unmapped.2_fastqc"
    ## 
    ## [[41]]
    ## [1] "Extracted overrepresented sequences for: 10770_S261_L002_unmapped.1_fastqc"
    ## 
    ## [[42]]
    ## [1] "Extracted overrepresented sequences for: 10770_S261_L002_unmapped.2_fastqc"
    ## 
    ## [[43]]
    ## [1] "Extracted overrepresented sequences for: Aso1_S262_L002_unmapped.1_fastqc"
    ## 
    ## [[44]]
    ## [1] "Extracted overrepresented sequences for: Aso1_S262_L002_unmapped.2_fastqc"
    ## 
    ## [[45]]
    ## [1] "Extracted overrepresented sequences for: Champ2_S263_L002_unmapped.1_fastqc"
    ## 
    ## [[46]]
    ## [1] "Extracted overrepresented sequences for: Champ2_S263_L002_unmapped.2_fastqc"
    ## 
    ## [[47]]
    ## [1] "Extracted overrepresented sequences for: Doc1_S264_L002_unmapped.1_fastqc"
    ## 
    ## [[48]]
    ## [1] "Extracted overrepresented sequences for: Doc1_S264_L002_unmapped.2_fastqc"
    ## 
    ## [[49]]
    ## [1] "Extracted overrepresented sequences for: Doc2_S265_L002_unmapped.1_fastqc"
    ## 
    ## [[50]]
    ## [1] "Extracted overrepresented sequences for: Doc2_S265_L002_unmapped.2_fastqc"
    ## 
    ## [[51]]
    ## [1] "Extracted overrepresented sequences for: Doc3_S266_L002_unmapped.1_fastqc"
    ## 
    ## [[52]]
    ## [1] "Extracted overrepresented sequences for: Doc3_S266_L002_unmapped.2_fastqc"
    ## 
    ## [[53]]
    ## [1] "Extracted overrepresented sequences for: Doc4_S287_L004_unmapped.1_fastqc"
    ## 
    ## [[54]]
    ## [1] "Extracted overrepresented sequences for: Doc4_S287_L004_unmapped.2_fastqc"
    ## 
    ## [[55]]
    ## [1] "Extracted overrepresented sequences for: Gron7_S267_L002_unmapped.1_fastqc"
    ## 
    ## [[56]]
    ## [1] "Extracted overrepresented sequences for: Gron7_S267_L002_unmapped.2_fastqc"
    ## 
    ## [[57]]
    ## [1] "Extracted overrepresented sequences for: Gron8_S268_L002_unmapped.1_fastqc"
    ## 
    ## [[58]]
    ## [1] "Extracted overrepresented sequences for: Gron8_S268_L002_unmapped.2_fastqc"
    ## 
    ## [[59]]
    ## [1] "Extracted overrepresented sequences for: Gron9_S269_L002_unmapped.1_fastqc"
    ## 
    ## [[60]]
    ## [1] "Extracted overrepresented sequences for: Gron9_S269_L002_unmapped.2_fastqc"
    ## 
    ## [[61]]
    ## [1] "Extracted overrepresented sequences for: Lamp1_S270_L002_unmapped.1_fastqc"
    ## 
    ## [[62]]
    ## [1] "Extracted overrepresented sequences for: Lamp1_S270_L002_unmapped.2_fastqc"
    ## 
    ## [[63]]
    ## [1] "Extracted overrepresented sequences for: Lamp2_S271_L002_unmapped.1_fastqc"
    ## 
    ## [[64]]
    ## [1] "Extracted overrepresented sequences for: Lamp2_S271_L002_unmapped.2_fastqc"
    ## 
    ## [[65]]
    ## [1] "Extracted overrepresented sequences for: Lamp4_S272_L002_unmapped.1_fastqc"
    ## 
    ## [[66]]
    ## [1] "Extracted overrepresented sequences for: Lamp4_S272_L002_unmapped.2_fastqc"
    ## 
    ## [[67]]
    ## [1] "Extracted overrepresented sequences for: Lamp5_S273_L002_unmapped.1_fastqc"
    ## 
    ## [[68]]
    ## [1] "Extracted overrepresented sequences for: Lamp5_S273_L002_unmapped.2_fastqc"
    ## 
    ## [[69]]
    ## [1] "Extracted overrepresented sequences for: Lamp6_S274_L002_unmapped.1_fastqc"
    ## 
    ## [[70]]
    ## [1] "Extracted overrepresented sequences for: Lamp6_S274_L002_unmapped.2_fastqc"
    ## 
    ## [[71]]
    ## [1] "Extracted overrepresented sequences for: Lamp8_S275_L002_unmapped.1_fastqc"
    ## 
    ## [[72]]
    ## [1] "Extracted overrepresented sequences for: Lamp8_S275_L002_unmapped.2_fastqc"
    ## 
    ## [[73]]
    ## [1] "Extracted overrepresented sequences for: Quill1_S276_L002_unmapped.1_fastqc"
    ## 
    ## [[74]]
    ## [1] "Extracted overrepresented sequences for: Quill1_S276_L002_unmapped.2_fastqc"
    ## 
    ## [[75]]
    ## [1] "Extracted overrepresented sequences for: Tisza2_S277_L002_unmapped.1_fastqc"
    ## 
    ## [[76]]
    ## [1] "Extracted overrepresented sequences for: Tisza2_S277_L002_unmapped.2_fastqc"
    ## 
    ## [[77]]
    ## [1] "Extracted overrepresented sequences for: Tisza3_S278_L002_unmapped.1_fastqc"
    ## 
    ## [[78]]
    ## [1] "Extracted overrepresented sequences for: Tisza3_S278_L002_unmapped.2_fastqc"

## Align sequences to 18s rRNA sequence

I took the top 75 overrepresented sequences of the unmapped reads for a
couple samples and blasted them against the nt database at NCBI. Some
aligned with a few bacterial species (which makes sense, as we performed
RNA sequencing on wild mushrooms) but the majority of them aligned with
the following gene:

Sequence ID: AY550243.1

Length: 5800

Amanita bisporigera intergenic spacer 2, partial sequence; 18S ribosomal
RNA gene, internal transcribed spacer 1, 5.8S ribosomal RNA gene,
internal transcribed spacer 2, and 28S ribosomal RNA gene, complete
sequence; and intergenic spacer 1, partial sequence

They typically aligned to the last ~400 bp of this sequence. However,
fastqc only includes 50 bp in their overrepresented sequence analysis.
Each of my reads are 150 bp and paired end, so to be thorough, I decided
to align all the unmapped reads from each sample to this ITS sequence. I
ran the following code on the supercomputer at MSI:

``` bash
module load bwa
cd /scratch.global/songs005/aphal_rnaseq/its18s

bwa index its18s.fna

bwa mem -t 10 its18s.fna ../hisat2_output_msdins_array_unmapp/10708_S250_L002/10708_S250_L002_unmapped.1.fastq ../hisat2_output_msdins_array_unmapp/10708_S250_L002/10708_S250_L002_unmapped.2.fastq
#Real time: 244.883 sec; CPU: 2140.206 sec

module load picard

java -jar /common/software/install/migrated/picard/2.25.6/picard.jar SortSam I=mAF.sam O=alignment.sorted.bam SORT_ORDER=coordinate
# Elapsed time: 10.92 minutes.

java -jar /common/software/install/migrated/picard/2.25.6/picard.jar BuildBamIndex INPUT=alignment.sorted.bam

#Elapsed time: 1.01 minutes.
rm mAF.sam
```

I then investigated the alignment in IGV, which confirmed this:

![IGV screenshot of unmapped read alignment to the A. bisporigera
sequence.](/Users/songs005/Downloads/image.png) ![Zooming in to the last
400 bp of the A. bisporigera
sequence.](/Users/songs005/Downloads/image2.png) We next took this ~400
bp region of the A. bisporigera gene and blasted it against multiple A.
phalloides genomes (from Drott et al. 2023 and NCBI), and found that
this gene locus is present in each genome but is truncated (only 70%
present) in the assembly we used in the RNAseq analysis. We have
therefore concluded that, on average ~47% of the mushroom RNA sequences
map to the reference genome due to 1) natural variation in bacteria and
other species present in the wild mushrooms we sampled and 2) the
truncation of this ITS housekeeping gene in the reference genome we
used.

### Session info

``` r
sessionInfo()
```

    ## R version 4.5.0 (2025-04-11)
    ## Platform: x86_64-apple-darwin20
    ## Running under: macOS Sequoia 15.5
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.5-x86_64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.5-x86_64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/Chicago
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] compiler_4.5.0    fastmap_1.2.0     cli_3.6.5         tools_4.5.0      
    ##  [5] htmltools_0.5.8.1 rstudioapi_0.17.1 yaml_2.3.10       rmarkdown_2.29   
    ##  [9] knitr_1.50        xfun_0.52         digest_0.6.37     rlang_1.1.6      
    ## [13] evaluate_1.0.4
