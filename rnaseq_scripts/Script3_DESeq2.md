# Extract mapping efficiency for hisat2

``` r
# import data and load packages
suppressMessages(library(DESeq2))
suppressMessages(library(rmarkdown))
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:Biobase':
    ## 
    ##     combine

    ## The following object is masked from 'package:matrixStats':
    ## 
    ##     count

    ## The following objects are masked from 'package:GenomicRanges':
    ## 
    ##     intersect, setdiff, union

    ## The following object is masked from 'package:GenomeInfoDb':
    ## 
    ##     intersect

    ## The following objects are masked from 'package:IRanges':
    ## 
    ##     collapse, desc, intersect, setdiff, slice, union

    ## The following objects are masked from 'package:S4Vectors':
    ## 
    ##     first, intersect, rename, setdiff, setequal, union

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     combine, intersect, setdiff, setequal, union

    ## The following object is masked from 'package:generics':
    ## 
    ##     explain

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(tidyr)
```

    ## 
    ## Attaching package: 'tidyr'

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     expand

``` r
library(readr)

dir <- getwd()
dir
```

    ## [1] "/Volumes/Liv/Lab_projects/Amanita_phalloides/bin/Scripts_for_publication"

``` r
dir2="/Volumes/Liv/Lab_projects/Amanita_phalloides"
# make output directories
subDir <- paste0(dir,"/Output_",Sys.Date())
dir.create(path=subDir, showWarnings = FALSE)
dir.create(paste0(subDir,"/Plots"), showWarnings = FALSE)
dir.create(paste0(subDir,"/NormCount-Graphs"), showWarnings = FALSE)
dir.create(paste0(subDir,"/NormCount-Graphs/Data"), showWarnings = FALSE)
```

``` r
# path to dir containing slurm-x.out files
# slurm_out_path <- "D:/Aphalloides_rnaseq_for_leaderless/hisat2_output_msdins_array/"
slurm_out_path <- paste0(dir2,"/MSI_outputs/hisat_output")
# get list of file names
out_files <- list.files(path=slurm_out_path, pattern=".out",full.names=TRUE,recursive=TRUE) 

# write a function to extract mapping efficiency

extract_hisat_mapping <- function(file_name) {
  # for debugging:
  #file_name = out_files[2]
  
  # open filename
  library(readr)
  temp <- read_tsv(file_name,show_col_types = FALSE)
  
  # Find the row index that contains "overall alignment rate"
  # this is typically third from the bottom
  index = nrow(temp)-2
  #temp2<-temp[grepl("overall alignment rate",temp[,1]),]
  mapping <- unlist(temp[index,1])
  
  
  
  # if this string contains %, then split by %, else print error
 if (grep("%", mapping)==TRUE) {
   
   # extract percent
   perc <- unlist(strsplit(mapping,split="%"))[1]
   
   # return this value to the console
   return(perc)
   
 } else print("Error: mapping rate not found")
  
}

# test it
extract_hisat_mapping(out_files[3])
```

    ## 2551749 reads; of these:1 
    ##                   "43.57"

``` r
# it works! now sapply to extract all the values

# mapping <- data.frame(Input = out_files, perc_mapping = sapply(out_files, extract_hisat_mapping))

# that didn't work because one or more files is not correctly formatted... use a loop
mapping <- data.frame(Array = list.files(path=slurm_out_path, pattern=".out",full.names=FALSE,recursive=TRUE) , perc_mapping = "")

# for (i in 1:length(out_files)) {
for (i in c(1:length(out_files))) {
  mapping[i,2] <- as.numeric(extract_hisat_mapping(out_files[i]))
} 

# find array index number for each sample
mapping$Array <- gsub("\\.out", "", mapping$Array)

# keep only second half of the Array name
mapping$Array <- sub(".*_", "", mapping$Array)

# Convert Array column to numeric and reorder
mapping <- mapping %>%
  mutate(Array = as.numeric(Array)) %>%
  arrange(Array)

# now match this to the metadata file
metadata <- readxl::read_xlsx(paste0(dir,"/metadata.xlsx"),sheet="metadata")
metadata$perc_mapping <- as.numeric(mapping$perc_mapping)
metadata[,c("Sample_ID","Population","Origin","perc_mapping")]
```

    ## # A tibble: 39 × 4
    ##    Sample_ID Population Origin perc_mapping
    ##    <chr>     <chr>      <chr>         <dbl>
    ##  1 Drake1    Drake      CA             37.0
    ##  2 Drake2    Drake      CA             56.4
    ##  3 Drake3    Drake      CA             44.7
    ##  4 Drake4    Drake      CA             31.6
    ##  5 Drake5    Drake      CA             67.2
    ##  6 Drake6    Drake      CA             44.6
    ##  7 Drake7    Drake      CA             42.2
    ##  8 Drake8    Drake      CA             43.6
    ##  9 Drake9    Drake      CA             46.2
    ## 10 Pet1      Pet        CA             43.5
    ## # ℹ 29 more rows

``` r
# find min, max, avg mapping rate
min(metadata$perc_mapping)
```

    ## [1] 4.31

``` r
max(metadata$perc_mapping)
```

    ## [1] 67.16

``` r
mean(metadata$perc_mapping)
```

    ## [1] 47.47615

``` r
sd(metadata$perc_mapping)
```

    ## [1] 10.72104

``` r
# print mapping
rmarkdown::paged_table(as.data.frame(metadata[,c("Sample_ID","Population","Origin","tot_reads_mill_trim","perc_mapping")]))
```

<script data-pagedtable-source type="application/json">
{"columns":[{"label":["Sample_ID"],"name":[1],"type":["chr"],"align":["left"]},{"label":["Population"],"name":[2],"type":["chr"],"align":["left"]},{"label":["Origin"],"name":[3],"type":["chr"],"align":["left"]},{"label":["tot_reads_mill_trim"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["perc_mapping"],"name":[5],"type":["dbl"],"align":["right"]}],"data":[{"1":"Drake1","2":"Drake","3":"CA","4":"74.4","5":"36.97"},{"1":"Drake2","2":"Drake","3":"CA","4":"56.8","5":"56.36"},{"1":"Drake3","2":"Drake","3":"CA","4":"54.4","5":"44.67"},{"1":"Drake4","2":"Drake","3":"CA","4":"48.6","5":"31.56"},{"1":"Drake5","2":"Drake","3":"CA","4":"76.6","5":"67.16"},{"1":"Drake6","2":"Drake","3":"CA","4":"60.0","5":"44.60"},{"1":"Drake7","2":"Drake","3":"CA","4":"75.2","5":"42.24"},{"1":"Drake8","2":"Drake","3":"CA","4":"58.2","5":"43.57"},{"1":"Drake9","2":"Drake","3":"CA","4":"38.2","5":"46.18"},{"1":"Pet1","2":"Pet","3":"CA","4":"31.2","5":"43.52"},{"1":"Pet1-redo","2":"Pet","3":"CA","4":"5.2","5":"43.57"},{"1":"Pet2","2":"Pet","3":"CA","4":"34.6","5":"38.62"},{"1":"Pet3","2":"Pet","3":"CA","4":"32.8","5":"31.37"},{"1":"Pet4","2":"Pet","3":"CA","4":"32.6","5":"50.14"},{"1":"Pet5","2":"Pet","3":"CA","4":"34.2","5":"4.31"},{"1":"Picnic1","2":"Picnic","3":"CA","4":"36.2","5":"52.81"},{"1":"Picnic2","2":"Picnic","3":"CA","4":"28.4","5":"41.24"},{"1":"Picnic3","2":"Picnic","3":"CA","4":"61.8","5":"40.67"},{"1":"Picnic4","2":"Picnic","3":"CA","4":"63.6","5":"43.72"},{"1":"Picnic5","2":"Picnic","3":"CA","4":"57.8","5":"45.35"},{"1":"Picnic6","2":"Picnic","3":"CA","4":"68.0","5":"52.84"},{"1":"Aso1","2":"Aso","3":"EU","4":"62.0","5":"41.16"},{"1":"Champ2","2":"Champ","3":"EU","4":"60.4","5":"53.60"},{"1":"Doc1","2":"Doc","3":"EU","4":"66.8","5":"44.12"},{"1":"Doc2","2":"Doc","3":"EU","4":"64.4","5":"50.78"},{"1":"Doc3","2":"Doc","3":"EU","4":"63.4","5":"50.87"},{"1":"Doc4","2":"Doc","3":"EU","4":"37.6","5":"55.68"},{"1":"Gron7","2":"Gron","3":"EU","4":"78.6","5":"55.69"},{"1":"Gron8","2":"Gron","3":"EU","4":"57.4","5":"57.22"},{"1":"Gron9","2":"Gron","3":"EU","4":"59.2","5":"57.83"},{"1":"Lamp1","2":"Lamp","3":"EU","4":"64.4","5":"59.69"},{"1":"Lamp2","2":"Lamp","3":"EU","4":"64.4","5":"53.22"},{"1":"Lamp4","2":"Lamp","3":"EU","4":"58.8","5":"39.31"},{"1":"Lamp5","2":"Lamp","3":"EU","4":"56.0","5":"56.41"},{"1":"Lamp6","2":"Lamp","3":"EU","4":"61.0","5":"56.49"},{"1":"Lamp8","2":"Lamp","3":"EU","4":"72.0","5":"55.78"},{"1":"Quill1","2":"Quill","3":"EU","4":"62.8","5":"50.53"},{"1":"Tisza2","2":"Tisza","3":"EU","4":"67.0","5":"55.58"},{"1":"Tisza3","2":"Tisza","3":"EU","4":"61.4","5":"56.14"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

``` r
# save as csv
write.csv(metadata,paste0(subDir,"/metadata_plus_mapping.csv"),row.names=FALSE)
```

# Import featurecounts data & filter genes that are not expressed

``` r
# set factor level; this is important for comparisons later
# group levels that are listed earlier in the list
# will be prioritized as the denominator when you do DESEq2 comparisons
# IE - the "wild type" or "control"
metadata$Origin <- factor(metadata$Origin, levels = c("EU", "CA"))


# paged_table(metadata)

# make sure all files are present
all(file.exists(metadata$Count_file))
```

    ## [1] TRUE

``` r
# if necessary - subset to only include samples where the file exists
metadata <- metadata[file.exists(metadata$Count_file),]
metadata[,c("Sample_ID","Population","Origin","perc_mapping")]
```

    ## # A tibble: 39 × 4
    ##    Sample_ID Population Origin perc_mapping
    ##    <chr>     <chr>      <fct>         <dbl>
    ##  1 Drake1    Drake      CA             37.0
    ##  2 Drake2    Drake      CA             56.4
    ##  3 Drake3    Drake      CA             44.7
    ##  4 Drake4    Drake      CA             31.6
    ##  5 Drake5    Drake      CA             67.2
    ##  6 Drake6    Drake      CA             44.6
    ##  7 Drake7    Drake      CA             42.2
    ##  8 Drake8    Drake      CA             43.6
    ##  9 Drake9    Drake      CA             46.2
    ## 10 Pet1      Pet        CA             43.5
    ## # ℹ 29 more rows

``` r
# find min, max, avg mapping rate
min(metadata$perc_mapping)
```

    ## [1] 4.31

``` r
max(metadata$perc_mapping)
```

    ## [1] 67.16

``` r
mean(metadata$perc_mapping)
```

    ## [1] 47.47615

``` r
all(file.exists(metadata$Count_file))
```

    ## [1] TRUE

``` r
# exclude pet5 (low mapping rate, 4%) and pet1-redo
metadata <- metadata[!metadata$Sample_ID %in% c("Pet1-redo","Pet5"),]
# import data - first get file names
files <- metadata$Count_file

# associate those with sample ID
names(files) <- paste(metadata$Sample_ID)

# import the featurecounts files
import_counts <- function(file_name) {
  # first line is for debugging
  # file_name = files[1]
  # extract sample prefix
  File_prefix <- unlist(strsplit(file_name,"/"))[8]
  
  temp <- read_tsv(file_name,comment="#",show_col_types = FALSE)
  # rename column to File_prefix
  colnames(temp)[7] <- File_prefix
  # output the column of counts
  output <- as.data.frame(temp[,7])
  rownames(output) <- temp$Geneid
  return(output)
}

# now run import_counts and cbind results together for all filenames in files
output_list <- lapply(files, import_counts)

# bind them together
cts <- as.matrix(do.call(cbind, output_list))

# rename rows to clean it up 
# when I originally ran Hisat2, my gff listed GPVFFA as GPVFFAY
# I also originally indexed the duplicated core sequences as core and core_1 instead of core_1 and core_2 so I will fix these here

old_msdin <- c("GPVFFAY","AWLATCP","AWLATCP_1","IWGIGCDP","IWGIGCDP_1","LIQRPFAP","LIQRPFAP_1")

# sort these alphabetical to match row order in cts
index <- row.names(cts) %in% old_msdin
row.names(cts)[index]
```

    ## [1] "AWLATCP_1"  "LIQRPFAP"   "AWLATCP"    "GPVFFAY"    "IWGIGCDP_1"
    ## [6] "IWGIGCDP"   "LIQRPFAP_1"

``` r
# now make sure new name order is the same
new_msdin <- c("AWLATCP_2","LIQRPFAP_1","AWLATCP_1","GPVFFA","IWGIGCDP_2","IWGIGCDP_1","LIQRPFAP_2")
row.names(cts)[index] <-new_msdin

# also fix messy gene IDs
row.names(cts) <- gsub("^[^-]+-|\\.m01$", "", row.names(cts))

head(cts,2)
```

    ##              10708_S250_L002 10709_S251_L002 10710_S252_L002 10711_S253_L002
    ## Ap.00g000010               0               0               0               0
    ## Ap.00g000020               0               0               0               0
    ##              10713_S254_L002 10716_S255_L002 10717_S256_L002 10720_S257_L002
    ## Ap.00g000010               0               0               0               0
    ## Ap.00g000020               0               0               0               0
    ##              10721_S279_L004 10745_S280_L004 10749_S281_L004 10751_S282_L004
    ## Ap.00g000010               0               0               0               0
    ## Ap.00g000020               0               0               0               0
    ##              10754_S283_L004 10763_S285_L004 10764_S286_L004 10765_S258_L002
    ## Ap.00g000010               0               0               0               0
    ## Ap.00g000020               0               0               0               0
    ##              10767_S259_L002 10769_S260_L002 10770_S261_L002 Aso1_S262_L002
    ## Ap.00g000010               0               0               0              0
    ## Ap.00g000020               0               0               0              0
    ##              Champ2_S263_L002 Doc1_S264_L002 Doc2_S265_L002 Doc3_S266_L002
    ## Ap.00g000010                0              0              0              0
    ## Ap.00g000020                0              0              0              0
    ##              Doc4_S287_L004 Gron7_S267_L002 Gron8_S268_L002 Gron9_S269_L002
    ## Ap.00g000010              0               0               0               0
    ## Ap.00g000020              0               0               0               0
    ##              Lamp1_S270_L002 Lamp2_S271_L002 Lamp4_S272_L002 Lamp5_S273_L002
    ## Ap.00g000010               0               0               0               0
    ## Ap.00g000020               0               0               0               0
    ##              Lamp6_S274_L002 Lamp8_S275_L002 Quill1_S276_L002 Tisza2_S277_L002
    ## Ap.00g000010               0               0                0                0
    ## Ap.00g000020               0               0                0                0
    ##              Tisza3_S278_L002
    ## Ap.00g000010                0
    ## Ap.00g000020                0

``` r
# check that the row order of metadata and column order of cts match

coldata <- as.data.frame(metadata[,c("Origin","Population")])
rownames(coldata) <- metadata$File_prefix

coldata$Origin <- factor(coldata$Origin)
coldata$Population <- factor(coldata$Population)

coldata
```

    ##                  Origin Population
    ## 10708_S250_L002      CA      Drake
    ## 10709_S251_L002      CA      Drake
    ## 10710_S252_L002      CA      Drake
    ## 10711_S253_L002      CA      Drake
    ## 10713_S254_L002      CA      Drake
    ## 10716_S255_L002      CA      Drake
    ## 10717_S256_L002      CA      Drake
    ## 10720_S257_L002      CA      Drake
    ## 10721_S279_L004      CA      Drake
    ## 10745_S280_L004      CA        Pet
    ## 10749_S281_L004      CA        Pet
    ## 10751_S282_L004      CA        Pet
    ## 10754_S283_L004      CA        Pet
    ## 10763_S285_L004      CA     Picnic
    ## 10764_S286_L004      CA     Picnic
    ## 10765_S258_L002      CA     Picnic
    ## 10767_S259_L002      CA     Picnic
    ## 10769_S260_L002      CA     Picnic
    ## 10770_S261_L002      CA     Picnic
    ## Aso1_S262_L002       EU        Aso
    ## Champ2_S263_L002     EU      Champ
    ## Doc1_S264_L002       EU        Doc
    ## Doc2_S265_L002       EU        Doc
    ## Doc3_S266_L002       EU        Doc
    ## Doc4_S287_L004       EU        Doc
    ## Gron7_S267_L002      EU       Gron
    ## Gron8_S268_L002      EU       Gron
    ## Gron9_S269_L002      EU       Gron
    ## Lamp1_S270_L002      EU       Lamp
    ## Lamp2_S271_L002      EU       Lamp
    ## Lamp4_S272_L002      EU       Lamp
    ## Lamp5_S273_L002      EU       Lamp
    ## Lamp6_S274_L002      EU       Lamp
    ## Lamp8_S275_L002      EU       Lamp
    ## Quill1_S276_L002     EU      Quill
    ## Tisza2_S277_L002     EU      Tisza
    ## Tisza3_S278_L002     EU      Tisza

``` r
all(rownames(coldata) %in% colnames(cts))
```

    ## [1] TRUE

``` r
all(rownames(coldata) == colnames(cts))
```

    ## [1] TRUE

``` r
# both of these are true!
# in case one of these was false, use the following logic:
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))
```

    ## [1] TRUE

``` r
# finally - generate dds object
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Origin)
```

    ## converting counts to integer mode

``` r
dds
```

    ## class: DESeqDataSet 
    ## dim: 8779 37 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(8779): Ap.00g000010 Ap.00g000020 ... Ap.00g068940 Ap.00g068950
    ## rowData names(0):
    ## colnames(37): 10708_S250_L002 10709_S251_L002 ... Tisza2_S277_L002
    ##   Tisza3_S278_L002
    ## colData names(2): Origin Population

``` r
# save cts object
write.csv(cts,paste0(subDir,"/all_gene_raw_counts.csv"))

# plot expression counts for all msdins quick
msdins <- read_tsv("/Volumes/Liv/Lab_projects/Amanita_phalloides/aphalloides_reff/reff/cds_for_liv.bed", comment = "#", col_names = FALSE,show_col_types = FALSE)

msdin_cores <- msdins$X4

msdin_cts <- cts[rownames(cts) %in% msdin_cores,]
msdin_cts <- msdin_cts[order(rownames(msdin_cts)), ]
write.csv(msdin_cts,paste0(subDir,"/msdin_raw_counts.csv"))

# get average for each row
msdin_means <- rowMeans(msdin_cts)
msdin_means
```

    ##    AWLATCP_1    AWLATCP_2      AWLVDCP   FFFPPFFIPP    FFPIVFSPP   FIFPPFFIPP 
    ## 1.048378e+02 0.000000e+00 2.210811e+01 5.567568e+00 1.135135e+01 1.864865e+00 
    ##       FMPLAP   FNILPFMLPP    FNLFRFPYP       GPVFFA      GVILIIP     HFASFIPP 
    ## 6.621622e+00 6.578378e+01 1.119027e+03 3.985341e+04 2.081081e+00 1.791892e+01 
    ##    IFLAFPIPP     IFWFIYFP     IIGILLPP   IRLPPLFLPP     ISDPTAYP   IWGIGCDP_1 
    ## 9.378378e+00 1.189189e+01 9.189189e+00 1.132432e+01 2.543243e+01 1.164865e+01 
    ##   IWGIGCDP_2     IWGIGCNP   LFFWFWFLWP     LGRPESLP   LILLAALGIP   LIQRPFAP_1 
    ## 4.254054e+01 6.591892e+01 2.829730e+01 1.470270e+01 4.532162e+02 2.764865e+01 
    ##   LIQRPFAP_2   LPILPIPPLP   LRLPPFMIPP     MLGFLVLP    MLPGMVAFS   MYNPPYFLPP 
    ## 1.251081e+02 1.081081e-01 2.181081e+01 4.535189e+04 1.378378e+00 4.216216e+00 
    ##      SFFFPIP    TIYYLYFIP     VQKPWSRP 
    ## 4.225981e+04 2.648649e+00 5.810811e+00

``` r
# get number of samples with a count of zeros for each msdin
num_zeros <- rowSums(msdin_cts == 0)
num_zeros
```

    ##  AWLATCP_1  AWLATCP_2    AWLVDCP FFFPPFFIPP  FFPIVFSPP FIFPPFFIPP     FMPLAP 
    ##          0         37          1          6          0         13          1 
    ## FNILPFMLPP  FNLFRFPYP     GPVFFA    GVILIIP   HFASFIPP  IFLAFPIPP   IFWFIYFP 
    ##          0          0          1         34          2          0          0 
    ##   IIGILLPP IRLPPLFLPP   ISDPTAYP IWGIGCDP_1 IWGIGCDP_2   IWGIGCNP LFFWFWFLWP 
    ##         21          6          1          1          1          0          0 
    ##   LGRPESLP LILLAALGIP LIQRPFAP_1 LIQRPFAP_2 LPILPIPPLP LRLPPFMIPP   MLGFLVLP 
    ##          2          0          0          0         36          9          0 
    ##  MLPGMVAFS MYNPPYFLPP    SFFFPIP  TIYYLYFIP   VQKPWSRP 
    ##         13          7          0          8          2

``` r
# print this:
summary_df <- data.frame(
  mean = msdin_means,
  num_zeros = num_zeros
)

# remove low-count genes that have a mean count below 5
# keep <- rowSums(counts(dds)) >= 10
keep <- rowMeans(counts(dds)) >= 5

sum(keep)
```

    ## [1] 8096

``` r
# total genes - counting msdins
length(dds)
```

    ## [1] 8779

``` r
# total genes - excluding msdins
length(dds) - length(msdin_cores)
```

    ## [1] 8746

``` r
# count number of msdins (out of 33) that are above 5
sum(msdin_means >= 5 )
```

    ## [1] 26

``` r
# percent
sum(msdin_means >= 5 ) / length(msdin_cores) * 100
```

    ## [1] 78.78788

``` r
# total genes detected - excluding msdins
sum(keep) - sum(msdin_means >= 5 )
```

    ## [1] 8070

``` r
# percent
(sum(keep) - sum(msdin_means >= 5 )) / (length(dds) - length(msdin_cores)) * 100
```

    ## [1] 92.27075

``` r
# find percentage of genes that have been kept after applying filter for low counts
sum(keep) / length(dds) * 100
```

    ## [1] 92.22007

``` r
# now subset your dds object to keep only "keep" rows
dds <- dds[keep,]
```

# Run DESeq2

``` r
dds <- DESeq(dds)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 93 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

``` r
res = results(dds, alpha=0.05, independentFiltering=FALSE, cooksCutoff=TRUE)
summary(res)
```

    ## 
    ## out of 8095 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 2490, 31%
    ## LFC < 0 (down)     : 2842, 35%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
resultsNames(dds)
```

    ## [1] "Intercept"       "Origin_CA_vs_EU"

``` r
# save results
resOrdered <- res[order(res$pvalue),]
# add column to indicate significance
resOrdered$significant_hits <- with(resOrdered, abs(log2FoldChange) > 1 & padj < 0.05)


write.csv(as.data.frame(resOrdered), 
          file=paste0(subDir,"/DESeq_output.csv"))
```

# QC results and generate first-pass plots

## Volcano plots

``` r
# summarize results from dds
# group1 will be the numerator (ie, treatment)
# group2 will be denominator (ie, control)
# gene_list is the list of genes you want to annotate on the volcano plot

make_volcano <- function(group1,group2,gene_list = msdins$X4) {
  # group1 = "CA"
  # group2 = "EU"
  # gene_list = msdins$X4
  print(paste0(group1," vs ",group2))
  
  res <- results(dds, contrast=c("Origin", group1, group2), alpha=0.05, independentFiltering=FALSE, cooksCutoff=TRUE)
  
  # now save deseq results as a dataframe
  deseqoutput <- as.data.frame(res)
  
  # plot it
  library(EnhancedVolcano)
  EnhancedVolcano(deseqoutput,
                  lab = rownames(res),
                  selectLab = gene_list,
                  title = paste0(group1," vs. ",group2),
                  subtitle = bquote(
                    italic("")),
                  caption = "",
                  x = 'log2FoldChange',
                  y = 'padj',
                  pCutoff = 0.05,
                  FCcutoff = 1,
                  ylab = bquote(~-Log[10] ~ italic("Adj. P value")),
                  legendPosition = 'bottom',
                  typeConnectors ="closed",
                  drawConnectors = TRUE,
                  widthConnectors = 1,
                  arrowheads = FALSE,
                  labSize = 4,
                  max.overlaps = 15,
                  lengthConnectors = unit(0.02, "npc")) +
    ggplot2::coord_cartesian(xlim=c(-6, 6), ylim=c(0,30)) +
    ggplot2::scale_x_continuous(breaks=seq(-5,5, 1))
  
  ggsave(paste0(subDir,"/Plots/",group1,"_vs_",group2,"_volcano.png"),
         units ="in",
         width = 10,
         height = 10)
  
  summary(res)
  
  # # shrink log fold changes
  # resLFC <- lfcShrink(dds, coef=paste0("Group_",group2,"_vs_",group1), type="apeglm")
  # resLFC
  # 
  # # use MA plot to visualize log2 fold changes
  # plotMA(res,ylim=c(-2,2))
  # 
  # # plot the shrunken LFC which should remove noise from low count genes
  # plotMA(resLFC, ylim=c(-2,2))
}
```

``` r
make_volcano(group2 = "EU",
             group1 = "CA")
```

    ## [1] "CA vs EU"

    ## Loading required package: ggplot2

    ## Loading required package: ggrepel

    ## Coordinate system already present. Adding new coordinate system, which will
    ## replace the existing one.
    ## Scale for x is already present.
    ## Adding another scale for x, which will replace the existing scale.

    ## Warning: ggrepel: 12 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

    ## 
    ## out of 8095 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 2490, 31%
    ## LFC < 0 (down)     : 2842, 35%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

## Plot counts for single genes

``` r
# make a plot for every gene and save it in a folder...
make_count_plot <- function(geneID) {
  # geneID = paste0(msdins$X4[10],".m01")
  gene <- plotCounts(dds, gene = geneID, intgroup = "Origin", returnData = TRUE)

  # Define color palette for Origin
  gene$Origin <- factor(gene$Origin, levels = c("CA",
                                                    "EU"))
  
  # also add pop information
  gene$Population <- metadata$Population   

  gene$Population <- factor(gene$Population, levels = c("Drake",
                                                            "Pet",
                                                            "Picnic",
                                                            "Aso",
                                                            "Champ",
                                                            "Doc",
                                                            "Gron",
                                                            "Lamp",
                                                            "Quill",
                                                            "Tisza"))
                          
  # genotype_colors <- c("#717568","#3F4739", "#FF01FB", "#0CCA4A","grey40","grey60")
  # colors for origin:
  # genotype_colors <- c("#FE4A49", "#2292A4")

genotype_colors <- c("#723D00","#fb8500","#ffb703",
                     "#023047","#0C556E","#177995","#219EBC","#5CB7CE","#96CFE0","#D1E8F2")

ggplot(data=gene,aes(x=Origin,y=count)) +
  geom_point(aes(colour = Population), size = 4, position = position_jitter(w = 0.3, h = 0)) +
  # add mean line
  stat_summary(fun= mean, fun.min=mean, fun.max=mean, geom="crossbar", width=0.5, linewidth=1,color="black") +
  # add error bars
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), geom="errorbar", width=0.2, linewidth=1,color="black") +
  theme_classic()  +
  scale_color_manual(values = genotype_colors) +
  theme(legend.position = "none",
        text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14)) +
  ylab("Normalized count") +
  ggtitle(paste0(geneID))

  # save graph
  ggsave(paste0(geneID,"_normalized_counts.png"),
         device="png",dpi="print",
         units = "in", width = 5, height = 3,
         path=paste0(subDir,"/NormCount-Graphs"))
  # save data table 
  write.csv(gene,paste0(subDir,"/NormCount-Graphs/Data/",geneID,".csv"))
}
```

``` r
# MSDINs, leadered and leaderless:
list <- msdins$X4

# first see which of these genes have any counts at all
# get counts matrix for all genes
counts <- counts(dds, normalized=TRUE)
# find list of RiPPs also in count matrix
list <- intersect(list, rownames(counts))

# loop through gene list
for(i in 1:length(list)) {
  genename = list[i]
  make_count_plot(geneID=genename) # takes ~1.5 hr to run on 9000 genes
}

# in the future we could also plot a handful of housekeeping genes

# next: generate facet_wrap graph for all MSDINS

# also: violin plot for vsd expression of all
```

## Variance stabilizing transformations

More info can be found at [DESeq2 Bioconductor
tutorial](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#exploring-and-exporting-results):

These figures plot the standard deviation of transformed data. A flat
curve of the square root of variance over the mean is ideal.

``` r
# Extracting transformed values
vsd <- vst(dds, blind=FALSE)
ntd <- normTransform(dds)

# extract the matrix of normalized values.
head(assay(vsd), 3)
```

    ##              10708_S250_L002 10709_S251_L002 10710_S252_L002 10711_S253_L002
    ## Ap.00g000050        3.931201        5.572686        4.894666        4.955907
    ## Ap.00g000070        4.385691        4.583703        4.452363        4.896566
    ## Ap.00g000080        3.573526        4.285382        3.887994        3.967977
    ##              10713_S254_L002 10716_S255_L002 10717_S256_L002 10720_S257_L002
    ## Ap.00g000050        5.840707        3.892063        4.118567        3.669032
    ## Ap.00g000070        4.684207        3.824621        4.801341        4.888208
    ## Ap.00g000080        5.095938        3.824621        3.702245        3.557715
    ##              10721_S279_L004 10745_S280_L004 10749_S281_L004 10751_S282_L004
    ## Ap.00g000050        3.139013        3.485095        4.743158        5.164954
    ## Ap.00g000070        4.225513        4.386186        4.675340        4.515019
    ## Ap.00g000080        3.139013        3.835854        3.948250        3.068117
    ##              10754_S283_L004 10763_S285_L004 10764_S286_L004 10765_S258_L002
    ## Ap.00g000050        6.278687        3.598048        4.785243        4.772664
    ## Ap.00g000070        4.518498        4.422308        4.307726        4.853514
    ## Ap.00g000080        4.164186        4.079780        4.154621        4.074898
    ##              10767_S259_L002 10769_S260_L002 10770_S261_L002 Aso1_S262_L002
    ## Ap.00g000050        4.477160        3.784712        4.042488       4.449782
    ## Ap.00g000070        4.657381        4.695893        4.781370       4.301738
    ## Ap.00g000080        4.017005        4.521756        3.947587       5.301503
    ##              Champ2_S263_L002 Doc1_S264_L002 Doc2_S265_L002 Doc3_S266_L002
    ## Ap.00g000050         3.968394       4.001065       5.165233       5.423174
    ## Ap.00g000070         4.741225       4.068286       4.336842       5.515404
    ## Ap.00g000080         4.553707       4.250250       3.900342       5.969678
    ##              Doc4_S287_L004 Gron7_S267_L002 Gron8_S268_L002 Gron9_S269_L002
    ## Ap.00g000050       4.338955        5.966887        6.790239        7.678984
    ## Ap.00g000070       5.083751        5.143560        4.769615        4.509287
    ## Ap.00g000080       5.662030        4.737057        4.360035        4.104067
    ##              Lamp1_S270_L002 Lamp2_S271_L002 Lamp4_S272_L002 Lamp5_S273_L002
    ## Ap.00g000050        4.537753        4.261712        3.851713        4.575055
    ## Ap.00g000070        4.867496        4.916421        4.129432        4.451162
    ## Ap.00g000080        4.656917        4.261712        3.851713        5.292489
    ##              Lamp6_S274_L002 Lamp8_S275_L002 Quill1_S276_L002 Tisza2_S277_L002
    ## Ap.00g000050        3.646494        8.055757         3.954890         4.509872
    ## Ap.00g000070        5.285198        5.044768         3.011301         4.436541
    ## Ap.00g000080        5.196493        6.014828         4.370828         3.875502
    ##              Tisza3_S278_L002
    ## Ap.00g000050         4.739591
    ## Ap.00g000070         4.359752
    ## Ap.00g000080         5.035299

``` r
# plots
library(vsn)
msd_dds <- meanSdPlot(assay(dds))
```

``` r
msd_vsd <- meanSdPlot(assay(vsd))
```

``` r
msd_dds$gg + ggtitle("Untransformed data") + ylim(-1,100)
```

    ## Warning: Removed 5807 rows containing non-finite outside the scale range
    ## (`stat_binhex()`).

    ## Warning: Removed 28 rows containing missing values or values outside the scale range
    ## (`geom_line()`).

``` r
msd_vsd$gg + ggtitle("Variance stabilized transformation")
```

<img src="Script3_DESeq2_files/figure-markdown_github/variance_plots-1.png" width="50%" /><img src="Script3_DESeq2_files/figure-markdown_github/variance_plots-2.png" width="50%" />

## Supervised clustering

Check the heatmap of a subset of the data to see if replicates look
similar to each other. This type of heatmap cluster is supervised, since
we are picking the 20 genes with the biggest changes.

It looks like 3 of the MSDINs are included in this list of top 20 genes
with biggest changes, including the putative leaderless!

``` r
# find 15 genes with highest expression levels from dds
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:15]

# make new data frame with only Population stage info
df <- as.data.frame(colData(dds)[,c("Population")])
colnames(df) <- c("Population")
#samples <- colnames(assay(ntd))
samples <- metadata$Sample_ID
rownames(df) <- samples

# heatmap of variance stabilized data
library(pheatmap)

# extract data
mat <- assay(vsd)[select,]
top15 <- row.names(mat)
colnames(mat) <- metadata$Sample_ID

# remove the random prefix and suffixes
row.names(mat) <- gsub("^[^-]+-|\\.m01$", "", row.names(mat))
row.names(mat)
```

    ##  [1] "Ap.00g078360" "Ap.00g054880" "Ap.00g019870" "Ap.00g054640" "Ap.00g057760"
    ##  [6] "MLGFLVLP"     "Ap.00g068810" "Ap.00g039900" "SFFFPIP"      "Ap.00g008070"
    ## [11] "GPVFFA"       "Ap.00g078370" "Ap.00g011790" "Ap.00g001030" "Ap.00g063180"

``` r
pheatmap(mat, cluster_rows=FALSE, 
         show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)
```

![](Script3_DESeq2_files/figure-markdown_github/supervised_clustering-1.png)

``` r
# also make heatmap of genes with custom colors
genotype_colors <- c("#723D00","#fb8500","#ffb703",
                     "#023047","#0C556E","#177995","#219EBC","#5CB7CE","#96CFE0","#D1E8F2")

annotation_colors <- list(Population = setNames(genotype_colors, unique(df$Population)))

# Define your custom color scale (from dark blue to aqua)
my_colors <- colorRampPalette(c("grey90","#03045e"))(100)


pheatmap(mat, cluster_rows=FALSE, 
         show_rownames=TRUE,
         cluster_cols=TRUE,
         annotation_col=df,
         annotation_colors = annotation_colors,
         legend=TRUE,
         color=my_colors,
         border_color = "transparent",
         fontsize_row = 10,  # Row label font size
         fontsize_col = 10)
```

![](Script3_DESeq2_files/figure-markdown_github/supervised_clustering-2.png)
##Do we know any functional information about the highly expressed
transcripts?

I used the following code to extract all the protein sequences that
match the highly expressed genes, and then I used the online NCBI blastp
to search human and Aspergillus fumigatus for closely related proteins
with known functions or names.

Next let’s figure out which genes are highly expressed…

``` r
library(rtracklayer)
library(Biostrings)
```

    ## Loading required package: XVector

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

``` r
library(GenomicFeatures)
```

    ## Loading required package: AnnotationDbi

    ## 
    ## Attaching package: 'AnnotationDbi'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

``` r
# Define file paths
genome_file <- "/Volumes/Liv/Lab_projects/Amanita_phalloides/aphalloides_reff/reff/10511_Aphal_PT_AllpathsLG_LINKS_jelly_pilon.fna"

# make sure chromosomes are assigned correctly
library(Rsamtools)

# Load Genome Sequence
genome_fa <- FaFile(genome_file)
class(genome)
```

    ## [1] "standardGeneric"
    ## attr(,"package")
    ## [1] "methods"

``` r
names(genome)
```

    ## NULL

``` r
# Load GFF3
gff <- paste0(dir,"/aphal_reff_with_msdins.gff3")
gff_sorted <- read_tsv(gff, comment = "#", col_names = FALSE,show_col_types = FALSE)

colnames(gff_sorted) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")

sum(is.na(gff_sorted$start))  # Count missing start positions
```

    ## [1] 0

``` r
sum(is.na(gff_sorted$end))
```

    ## [1] 0

``` r
# Create TxDb from GFF3
options(scipen = 999)
txdb <- makeTxDbFromGFF(gff, format="gff3")
```

    ## Warning in call_fun_in_txdbmaker("makeTxDbFromGFF", ...): makeTxDbFromGFF() has moved from GenomicFeatures to the txdbmaker package,
    ##   and is formally deprecated in GenomicFeatures >= 1.59.1. Please call
    ##   txdbmaker::makeTxDbFromGFF() to get rid of this warning.

    ## Import genomic features from the file as a GRanges object ...

    ## OK

    ## Prepare the 'metadata' data frame ... OK
    ## Make the TxDb object ...

    ## Warning in .makeTxDb_normarg_chrominfo(chrominfo): genome version information
    ## is not available for this TxDb object

    ## OK

``` r
# save txdb object
saveDb(txdb, file = paste0(dir,"/txdb.sqlite"))
```

    ## TxDb object:
    ## # Db type: TxDb
    ## # Supporting package: GenomicFeatures
    ## # Data source: /Volumes/Liv/Lab_projects/Amanita_phalloides/bin/Scripts_for_publication/aphal_reff_with_msdins.gff3
    ## # Organism: NA
    ## # Taxonomy ID: NA
    ## # miRBase build ID: NA
    ## # Genome: NA
    ## # Nb of transcripts: 8812
    ## # Db created by: txdbmaker package from Bioconductor
    ## # Creation time: 2025-08-11 11:31:37 -0500 (Mon, 11 Aug 2025)
    ## # txdbmaker version at creation time: 1.4.2
    ## # RSQLite version at creation time: 2.4.2
    ## # DBSCHEMAVERSION: 1.2

``` r
# Extract CDS sequences from genome
transcripts <- exonsBy(txdb, by="tx", use.names=TRUE)

cds_seqs <- extractTranscriptSeqs(genome_fa, transcripts)

# find cds that might not contain complete codons (length 3)
invalid_cds <- names(cds_seqs)[width(cds_seqs) %% 3 != 0]
invalid_cds
```

    ##   [1] "26557-Ap.00g000010.m01" "26557-Ap.00g002410.m01" "26557-Ap.00g069340.m01"
    ##   [4] "26557-Ap.00g069850.m01" "26557-Ap.00g074640.m01" "26557-Ap.00g074690.m01"
    ##   [7] "26557-Ap.00g028120.m01" "26557-Ap.00g077690.m01" "26557-Ap.00g078010.m01"
    ##  [10] "AWLATCP_2"              "AWLATCP_2.m01"          "FNLFRFPYP"             
    ##  [13] "FNLFRFPYP.m01"          "26557-Ap.00g078580.m01" "26557-Ap.00g078590.m01"
    ##  [16] "LIQRPFAP_1"             "LIQRPFAP_1.m01"         "26557-Ap.00g028540.m01"
    ##  [19] "26557-Ap.00g078640.m01" "26557-Ap.00g079210.m01" "26557-Ap.00g079670.m01"
    ##  [22] "26557-Ap.00g080450.m01" "26557-Ap.00g081210.m01" "SFFFPIP"               
    ##  [25] "SFFFPIP.m01"            "26557-Ap.00g031270.m01" "LFFWFWFLWP"            
    ##  [28] "LFFWFWFLWP.m01"         "26557-Ap.00g030190.m01" "26557-Ap.00g082400.m01"
    ##  [31] "26557-Ap.00g082490.m01" "26557-Ap.00g082500.m01" "26557-Ap.00g082510.m01"
    ##  [34] "26557-Ap.00g082650.m01" "FFFPPFFIPP"             "FFFPPFFIPP.m01"        
    ##  [37] "MYNPPYFLPP"             "MYNPPYFLPP.m01"         "VQKPWSRP"              
    ##  [40] "VQKPWSRP.m01"           "LGRPESLP"               "LGRPESLP.m01"          
    ##  [43] "IRLPPLFLPP"             "IRLPPLFLPP.m01"         "LRLPPFMIPP"            
    ##  [46] "LRLPPFMIPP.m01"         "FIFPPFFIPP"             "FIFPPFFIPP.m01"        
    ##  [49] "26557-Ap.00g083540.m01" "26557-Ap.00g083580.m01" "26557-Ap.00g083810.m01"
    ##  [52] "26557-Ap.00g083820.m01" "26557-Ap.00g032470.m01" "26557-Ap.00g033090.m01"
    ##  [55] "26557-Ap.00g083830.m01" "26557-Ap.00g083850.m01" "26557-Ap.00g033600.m01"
    ##  [58] "26557-Ap.00g084160.m01" "MLPGMVAFS"              "MLPGMVAFS.m01"         
    ##  [61] "26557-Ap.00g084380.m01" "26557-Ap.00g084500.m01" "26557-Ap.00g003000.m01"
    ##  [64] "26557-Ap.00g005340.m01" "AWLVDCP"                "AWLVDCP.m01"           
    ##  [67] "26557-Ap.00g033680.m01" "26557-Ap.00g084700.m01" "26557-Ap.00g084970.m01"
    ##  [70] "26557-Ap.00g034010.m01" "26557-Ap.00g085350.m01" "26557-Ap.00g085650.m01"
    ##  [73] "26557-Ap.00g085660.m01" "26557-Ap.00g085760.m01" "LILLAALGIP"            
    ##  [76] "LILLAALGIP.m01"         "FNILPFMLPP"             "FNILPFMLPP.m01"        
    ##  [79] "GVILIIP"                "GVILIIP.m01"            "LPILPIPPLP"            
    ##  [82] "LPILPIPPLP.m01"         "AWLATCP_1"              "AWLATCP_1.m01"         
    ##  [85] "IIGILLPP"               "IIGILLPP.m01"           "26557-Ap.00g085780.m01"
    ##  [88] "26557-Ap.00g085810.m01" "FFPIVFSPP"              "FFPIVFSPP.m01"         
    ##  [91] "26557-Ap.00g035390.m01" "GPVFFA"                 "GPVFFA.m01"            
    ##  [94] "26557-Ap.00g036150.m01" "26557-Ap.00g038290.m01" "26557-Ap.00g086220.m01"
    ##  [97] "MLGFLVLP"               "MLGFLVLP.m01"           "26557-Ap.00g038790.m01"
    ## [100] "26557-Ap.00g086250.m01" "26557-Ap.00g006790.m01" "26557-Ap.00g038890.m01"
    ## [103] "26557-Ap.00g038900.m01" "26557-Ap.00g086340.m01" "26557-Ap.00g086350.m01"
    ## [106] "26557-Ap.00g039780.m01" "26557-Ap.00g086430.m01" "26557-Ap.00g041060.m01"
    ## [109] "26557-Ap.00g041080.m01" "ISDPTAYP"               "ISDPTAYP.m01"          
    ## [112] "26557-Ap.00g041800.m01" "26557-Ap.00g086830.m01" "26557-Ap.00g086840.m01"
    ## [115] "26557-Ap.00g087100.m01" "26557-Ap.00g087150.m01" "26557-Ap.00g008720.m01"
    ## [118] "26557-Ap.00g013880.m01" "26557-Ap.00g049740.m01" "26557-Ap.00g087170.m01"
    ## [121] "IWGIGCNP"               "IWGIGCNP.m01"           "IWGIGCDP_1"            
    ## [124] "IWGIGCDP_1.m01"         "IWGIGCDP_2"             "IWGIGCDP_2.m01"        
    ## [127] "26557-Ap.00g050370.m01" "26557-Ap.00g049970.m01" "26557-Ap.00g050380.m01"
    ## [130] "26557-Ap.00g050390.m01" "26557-Ap.00g050400.m01" "26557-Ap.00g051340.m01"
    ## [133] "FMPLAP"                 "FMPLAP.m01"             "IFLAFPIPP"             
    ## [136] "IFLAFPIPP.m01"          "HFASFIPP"               "HFASFIPP.m01"          
    ## [139] "TIYYLYFIP"              "TIYYLYFIP.m01"          "26557-Ap.00g087420.m01"
    ## [142] "26557-Ap.00g013920.m01" "26557-Ap.00g014520.m01" "26557-Ap.00g053140.m01"
    ## [145] "26557-Ap.00g055300.m01" "26557-Ap.00g014670.m01" "26557-Ap.00g014680.m01"
    ## [148] "26557-Ap.00g018860.m01" "26557-Ap.00g055490.m01" "26557-Ap.00g057600.m01"
    ## [151] "LIQRPFAP_2"             "LIQRPFAP_2.m01"         "26557-Ap.00g058170.m01"
    ## [154] "26557-Ap.00g023560.m01" "26557-Ap.00g066240.m01" "26557-Ap.00g068680.m01"
    ## [157] "IFWFIYFP"               "IFWFIYFP.m01"

``` r
# check for any invalid bases in the CDS sequences
invalid_bases <- grepl("[^ATCG]", as.character(cds_seqs))
invalid_seqs <- cds_seqs[invalid_bases]
invalid_seqs
```

    ## DNAStringSet object of length 8:
    ##     width seq                                               names               
    ## [1]   633 ATGCAGAATCAACCTTCCGGTTC...TCTTAATTGATTTGATGCTATAA 26557-Ap.00g02881...
    ## [2]  1416 ATGTCCACGCTCGACGCGCTCTT...AGTTTGGGATCGTTTATGGCTAA 26557-Ap.00g08106...
    ## [3]   165 ATGGGTTCGTATTCTTCCAAATA...CTACCAAAGGCACTGATTTCTAG 26557-Ap.00g03087...
    ## [4]  1176 ATGACCGAAACCCCCCATGGCAG...AAGCAATGTTCTGGAAATATTAG 26557-Ap.00g08395...
    ## [5]   861 ATGGTTGATCCTGTTCCCAATGT...GAGATCTAATAANCAAAATATAG 26557-Ap.00g04095...
    ## [6]   957 ATGGATGTTCCCGTTGCACAGCA...GTCGNGATGCTGAGAGCCAATAG 26557-Ap.00g00932...
    ## [7]   915 ATGAACACCAAACGTGTCACAGA...AACGACCGTCCATGGACGGTTAA 26557-Ap.00g05106...
    ## [8]  2004 ATGACGCAGCGACCAGACGCTCT...CCGTAAAATATCTGTTGAGTTGA 26557-Ap.00g01398...

``` r
# remove seqs with invalid bases 
cds_seqs <- cds_seqs[!invalid_bases]

# Translate to protein sequences
protein_seqs <- translate(cds_seqs)
```

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[1]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[259]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[357]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[382]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[887]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[893]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[1456]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[1472]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[1501]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[1505]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[1506]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[1558]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[1559]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[1562]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[1563]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[1573]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[1574]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[1590]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[1628]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[1674]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[1742]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[1783]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[1940]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2019]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2020]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2098]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2124]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2125]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2142]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2262]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2271]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2272]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2273]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2285]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2398]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2399]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2451]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2452]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2453]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2454]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2455]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2456]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2460]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2461]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2462]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2463]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2464]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2465]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2509]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2525]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2536]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2537]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2538]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2600]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2601]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2603]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2685]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2687]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2699]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2700]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2711]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2718]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2728]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2848]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2935]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2936]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2967]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[2976]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3044]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3092]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3119]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3161]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3162]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3172]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3183]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3184]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3186]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3187]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3188]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3189]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3190]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3191]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3200]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3201]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3229]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3230]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3246]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3251]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3315]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3316]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3329]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3477]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3478]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3480]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3610]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3635]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3657]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3658]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3668]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3669]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3744]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3829]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3833]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3838]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3839]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[3932]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[4001]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[4058]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[4095]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[4269]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[4270]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[4325]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[4382]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[4383]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[4667]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[4813]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[5141]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[5407]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[5707]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[5710]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[5718]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[5719]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[5720]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[5721]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[5733]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[5734]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[5765]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[5766]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[5803]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[5804]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[5806]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[5848]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[5855]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[5856]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[5859]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[5860]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[5861]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[5862]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[5912]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[5913]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[5914]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[6033]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[6058]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[6156]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[6320]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[6381]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[6382]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[6951]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[7106]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[7317]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[7495]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[7496]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[7553]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[8366]]': last base was ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[8531]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[8775]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[8800]]': last 2 bases were ignored

    ## Warning in .Call2("DNAStringSet_translate", x, skip_code,
    ## dna_codes[codon_alphabet], : in 'x[[8801]]': last 2 bases were ignored

``` r
head(names(protein_seqs))
```

    ## [1] "26557-Ap.00g000010.m01" "26557-Ap.00g000020.m01" "26557-Ap.00g000030.m01"
    ## [4] "26557-Ap.00g000050.m01" "26557-Ap.00g000060.m01" "26557-Ap.00g000100.m01"

``` r
# fix names - remove "26557-" and ".m*"
names(protein_seqs) <- sub("^26557-", "", names(protein_seqs))
names(protein_seqs) <- sub("\\.m.*$", "", names(protein_seqs))

head(names(protein_seqs))
```

    ## [1] "Ap.00g000010" "Ap.00g000020" "Ap.00g000030" "Ap.00g000050" "Ap.00g000060"
    ## [6] "Ap.00g000100"

``` r
# Write to a FASTA file
writeXStringSet(protein_seqs, paste0(dir,"/aphal_proteins.fasta"))

cat("Protein sequences saved to aphal_proteins.fasta\n")
```

    ## Protein sequences saved to aphal_proteins.fasta

### Now extract highly expressed genes

``` r
top15
```

    ##  [1] "Ap.00g078360" "Ap.00g054880" "Ap.00g019870" "Ap.00g054640" "Ap.00g057760"
    ##  [6] "MLGFLVLP"     "Ap.00g068810" "Ap.00g039900" "SFFFPIP"      "Ap.00g008070"
    ## [11] "GPVFFA"       "Ap.00g078370" "Ap.00g011790" "Ap.00g001030" "Ap.00g063180"

``` r
# extract from fasta
# if reading in a new fasta, use the following:
# fasta <- paste0(dir,"/aphal_proteins.fasta")
# #BiocManager::install("Biostrings")
# 
# library(Biostrings)
# fasta_seqs <- readAAStringSet(fasta)

# Match names in FASTA with entries in top15
matched_seqs <- unique(protein_seqs[names(protein_seqs) %in% top15])
# remove duplicates - msdins are showing up twice ...
matched_seqs
```

    ## AAStringSet object of length 15:
    ##      width seq                                              names               
    ##  [1]   127 MRDRDTGRSRGFGFVTYSTNDEA...QGGYGGGGYSQGYQQGSGFGY* Ap.00g001030
    ##  [2]   468 MVASKTSILVAAVALGAASVYAA...SLGLDEDYFTKRFYDEDVSEF* Ap.00g078360
    ##  [3]   477 MVASKTSIIVAAVALGAATAYAA...SLGLDLDYFTKRFYDEDGLEY* Ap.00g078370
    ##  [4]    35 MSDINAARLPSFFFPIPCISDDIEMVLTRGERPLL              SFFFPIP
    ##  [5]    34 MSDVNTIRIPGPVFFAYVGDEVDNVLRSGERPFL               GPVFFA
    ##  ...   ... ...
    ## [11]   184 MVKFTSVFFVFATLALSVVASPV...AMRDKETVDAAFKQVIDYYSH* Ap.00g054880
    ## [12]   324 MPSPTVLRLLALASLALLAALGP...NRLLGSNGQPTSLGWYYVNQY* Ap.00g019870
    ## [13]   294 MQKSLILSVLLAFSLVAMCAPAV...TGFDEKITAVSTYWGYLRKST* Ap.00g057760
    ## [14]   233 MFFNRFFVLFAIVLYSLALITPH...DNLADEAVGRLFNIHWHWLFV* Ap.00g063180
    ## [15]   779 MYSLHYFLFLCTILSYTGVVFGL...AGVINQSGNANGYASTVTAWY* Ap.00g068810

``` r
# save
writeXStringSet(matched_seqs, filepath = paste0(dir,"/top15_sequences.fasta"))
```

## PCA / unsupervised clustering

A great interactive explanation of PCA plots is
[here](https://bioboot.github.io/bggn213_F20/class-material/pca/)

``` r
# find distances between samples
sampleDists <- dist(t(assay(vsd)))

# reformat distances
sampleDistMatrix <- as.matrix(sampleDists)

# make row and column names the same as sample IDs
rownames(sampleDistMatrix) <- rownames(df)
colnames(sampleDistMatrix) <- rownames(df)

# generate some nice colors
library(RColorBrewer)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

# plot it
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
```

![](Script3_DESeq2_files/figure-markdown_github/cluster_heatmap_and_pca-1.png)

``` r
# next lets generate a PCA plot on variance stabilized data
vsd_out <- assay(vsd)

# run PCA analysis
pca <- prcomp(t(vsd_out))

# plot scree
library(factoextra)
```

    ## Welcome! Want to learn more? See two factoextra-related books at https://goo.gl/ve3WBa

``` r
fviz_eig(pca)
```

![](Script3_DESeq2_files/figure-markdown_github/cluster_heatmap_and_pca-2.png)

``` r
# Create data frame with metadata and PC3 and PC4 values for input to ggplot
pcaData <- cbind(metadata[,c(2:5,8,13,15,18)], pca$x)

# Define color palette for genotypes
pcaData$Origin <- factor(pcaData$Origin, levels = c("CA",
                                                    "EU"))
genotype_colors <- c("#723D00","#fb8500","#ffb703",
                     "#023047","#0C556E","#177995","#219EBC","#5CB7CE","#96CFE0","#D1E8F2")

genotype_colors <- c("#fb8500", "#5CB7CE")

# find percent variance for PC1 and PC2
# summary of variance
pcasummary <- summary(pca)
pcasummary
```

    ## Importance of components:
    ##                            PC1     PC2      PC3      PC4      PC5      PC6
    ## Standard deviation     49.7044 25.8066 19.13221 14.78744 13.52954 11.80978
    ## Proportion of Variance  0.4717  0.1272  0.06989  0.04175  0.03495  0.02663
    ## Cumulative Proportion   0.4717  0.5989  0.66875  0.71051  0.74546  0.77209
    ##                             PC7     PC8     PC9    PC10    PC11   PC12    PC13
    ## Standard deviation     10.29650 9.94409 9.35800 9.14385 8.16263 8.0259 7.62119
    ## Proportion of Variance  0.02024 0.01888 0.01672 0.01596 0.01272 0.0123 0.01109
    ## Cumulative Proportion   0.79233 0.81121 0.82793 0.84389 0.85661 0.8689 0.88000
    ##                           PC14    PC15    PC16    PC17    PC18    PC19   PC20
    ## Standard deviation     7.31282 7.06540 6.88098 6.48477 6.33989 6.11415 5.9695
    ## Proportion of Variance 0.01021 0.00953 0.00904 0.00803 0.00767 0.00714 0.0068
    ## Cumulative Proportion  0.89021 0.89975 0.90879 0.91681 0.92449 0.93163 0.9384
    ##                           PC21    PC22    PC23    PC24    PC25    PC26    PC27
    ## Standard deviation     5.92059 5.63675 5.56342 5.07650 4.98911 4.78476 4.60256
    ## Proportion of Variance 0.00669 0.00607 0.00591 0.00492 0.00475 0.00437 0.00404
    ## Cumulative Proportion  0.94512 0.95119 0.95710 0.96202 0.96677 0.97114 0.97519
    ##                           PC28   PC29    PC30    PC31    PC32    PC33    PC34
    ## Standard deviation     4.53083 4.4004 4.21240 4.02068 3.99395 3.62851 3.34748
    ## Proportion of Variance 0.00392 0.0037 0.00339 0.00309 0.00305 0.00251 0.00214
    ## Cumulative Proportion  0.97911 0.9828 0.98619 0.98928 0.99233 0.99484 0.99698
    ##                           PC35   PC36                PC37
    ## Standard deviation     2.91502 2.7063 0.00000000000004633
    ## Proportion of Variance 0.00162 0.0014 0.00000000000000000
    ## Cumulative Proportion  0.99860 1.0000 1.00000000000000000

``` r
pc1var <- round(pcasummary$importance[2],digits = 4) * 100
pc2var <- round(pcasummary$importance[5],digits = 4) * 100
pc3var <- round(pcasummary$importance[8],digits = 4) * 100
pc4var <- round(pcasummary$importance[11],digits = 4) * 100


# plot PC1 and 2
ggplot(pcaData,aes(x=PC1, y=PC2, color = Origin, shape = Origin)) +
  geom_point(size=3) + 
  theme_classic() +
  scale_color_manual(values = genotype_colors) +
  xlab(paste0("PC1 (",pc1var,"% variance)")) +
  ylab(paste0("PC2 (",pc2var,"% variance)"))
```

![](Script3_DESeq2_files/figure-markdown_github/cluster_heatmap_and_pca-3.png)

``` r
ggsave("PCA_plot_pc12_origin.png",
       path = paste0(subDir,"/Plots/"),
       units="in",
       width = 6, height = 5)

# plot PC3 and 4
ggplot(pcaData,aes(x=PC3, y=PC4, color = Origin, size = Origin)) +
  geom_point(size=3) + 
  theme_classic() +
  scale_color_manual(values = genotype_colors) +
  xlab(paste0("PC3 (",pc3var,"% variance)")) +
  ylab(paste0("PC4 (",pc4var,"% variance)"))
```

![](Script3_DESeq2_files/figure-markdown_github/cluster_heatmap_and_pca-4.png)

``` r
ggsave("PCA_plot_pc34_origin.png",
       path = paste0(subDir,"/Plots/"),
       units="in",
       width = 6, height = 5)

# also plot where colors indicate population & shape indicates Origin

# deduplicate
t <- unique(metadata[,c(4,8)])
t
```

    ## # A tibble: 10 × 2
    ##    Population Origin
    ##    <chr>      <fct> 
    ##  1 Drake      CA    
    ##  2 Pet        CA    
    ##  3 Picnic     CA    
    ##  4 Aso        EU    
    ##  5 Champ      EU    
    ##  6 Doc        EU    
    ##  7 Gron       EU    
    ##  8 Lamp       EU    
    ##  9 Quill      EU    
    ## 10 Tisza      EU

``` r
pcaData$Population <- factor(pcaData$Population, levels = c("Drake",
                                                            "Pet",
                                                            "Picnic",
                                                            "Aso",
                                                            "Champ",
                                                            "Doc",
                                                            "Gron",
                                                            "Lamp",
                                                            "Quill",
                                                            "Tisza"))
# first 3 are CA, last 7 are EU
genotype_colors <- c("#723D00","#fb8500","#ffb703",
                     "#023047","#0C556E","#177995","#219EBC","#5CB7CE","#96CFE0","#D1E8F2")


ggplot(pcaData,aes(x=PC1, y=PC2, color = Population, shape = Origin)) +
  geom_point(size=4) + 
  theme_classic() +
  scale_color_manual(values = genotype_colors) +
stat_ellipse(aes(group = Origin), colour = "black")+                 
  xlab(paste0("PC1 (",pc1var,"% variance)")) +
  ylab(paste0("PC2 (",pc2var,"% variance)"))
```

![](Script3_DESeq2_files/figure-markdown_github/cluster_heatmap_and_pca-5.png)

``` r
ggsave("PCA_plot_pc12_population_origin.png",
       path = paste0(subDir,"/Plots/"),
       units="in",
       width = 6, height = 5)

ggplot(pcaData,aes(x=PC3, y=PC4, color = Population, shape = Origin)) +
  geom_point(size=3) + 
  theme_classic() +
  scale_color_manual(values = genotype_colors) +
  xlab(paste0("PC3 (",pc3var,"% variance)")) +
  ylab(paste0("PC4 (",pc4var,"% variance)"))
```

![](Script3_DESeq2_files/figure-markdown_github/cluster_heatmap_and_pca-6.png)

``` r
ggsave("PCA_plot_pc34_population_origin.png",
       path = paste0(subDir,"/Plots/"),
       units="in",
       width = 6, height = 5)

## now plot PCA where color indicates mapping rate
ggplot(pcaData,aes(x=PC1, y=PC2, color = perc_mapping, shape = Origin)) +
  geom_point(size=3) + 
  theme_classic() +
  scale_color_viridis_c() +
  #paletteer::scale_color_paletteer_c("viridis::plasma") +
  #scale_color_manual(values = genotype_colors) +
  xlab(paste0("PC1 (",pc1var,"% variance)")) +
  ylab(paste0("PC2 (",pc2var,"% variance)"))
```

![](Script3_DESeq2_files/figure-markdown_github/cluster_heatmap_and_pca-7.png)

``` r
ggsave("PCA_plot_pc12_percmapping_origin.png",
       path = paste0(subDir,"/Plots/"),
       units="in",
       width = 6, height = 5)

# plot mapping rates for both pops
ggplot(pcaData, aes(x=Origin,y=perc_mapping, color=Population,shape=Origin)) +
  geom_point(size = 4, position = position_jitter(w = 0.3, h = 0)) +
  # add mean line
  stat_summary(fun= mean, fun.min=mean, fun.max=mean, geom="crossbar", width=0.5, linewidth=1,color="black") +
  # add error bars
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), geom="errorbar", width=0.2, linewidth=1,color="black") + 
  scale_color_manual(values = genotype_colors) +
  theme_classic()
```

![](Script3_DESeq2_files/figure-markdown_github/cluster_heatmap_and_pca-8.png)

``` r
# do a statistical test
# Perform t-test comparing perc_mapping values between EU and CA
t_test_result <- t.test(perc_mapping ~ Origin, data = pcaData)

# View t-test results
print(t_test_result)
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  perc_mapping by Origin
    ## t = -3.3114, df = 31.932, p-value = 0.002312
    ## alternative hypothesis: true difference in means between group CA and group EU is not equal to 0
    ## 95 percent confidence interval:
    ##  -12.691346  -3.023741
    ## sample estimates:
    ## mean in group CA mean in group EU 
    ##         44.92579         52.78333

``` r
ggsave("percmapping_origin_population.png",
       path = paste0(subDir,"/Plots/"),
       units="in",
       width = 4, height = 5)

# export pca data for separate graphing
write.csv(pcaData[,1:10],paste0(subDir,"/pca_data.csv"),row.names=FALSE)


# also plot norm counts for leaderless with population colored...
gene <- plotCounts(dds, gene = "MLGFLVLP", intgroup = "Origin", returnData = TRUE)

  # Define color palette for Origin
  gene$Origin <- factor(gene$Origin, levels = c("CA",
                                                    "EU"))
  gene$Population <- metadata$Population   

  gene$Population <- factor(gene$Population, levels = c("Drake",
                                                            "Pet",
                                                            "Picnic",
                                                            "Aso",
                                                            "Champ",
                                                            "Doc",
                                                            "Gron",
                                                            "Lamp",
                                                            "Quill",
                                                            "Tisza"))
ggplot(data=gene,aes(x=Origin,y=count)) +
  geom_point(aes(colour = Population), size = 4, position = position_jitter(w = 0.3, h = 0)) +
  # add mean line
  stat_summary(fun= mean, fun.min=mean, fun.max=mean, geom="crossbar", width=0.5, linewidth=1,color="black") +
  # add error bars
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), geom="errorbar", width=0.2, linewidth=1,color="black") +
  theme_classic()  +
  scale_color_manual(values = genotype_colors) +
  theme(legend.position = "none",
        text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14)) +
  ylab("Normalized count") +
  ggtitle("MLGFLVLP")
```

![](Script3_DESeq2_files/figure-markdown_github/cluster_heatmap_and_pca-9.png)

``` r
  # save graph
  ggsave(paste0("MLGFLVLP","_pop_normalized_counts.png"),
         device="png",dpi="print",
         units = "in", width = 5, height = 3,
         path=paste0(subDir,"/NormCount-Graphs"))
```

## additional QC - Histogram of overall gene expression for EU vs CA

``` r
# get norm counts for all genes
norm_counts <- counts(dds, normalized=TRUE)

# aggregate by origin (CA or EU)
# Step 1: Read metadata
#metadata <- readxl::read_xlsx(paste0(dir,"/metadata.xlsx"),sheet="metadata")

# Step 2: Make sure metadata and norm_counts align
all(colnames(norm_counts) %in% metadata$File_prefix)  # Should return TRUE
```

    ## [1] TRUE

``` r
# Step 3: Get Origin for each sample
origin_vector <- metadata$Origin[match(colnames(norm_counts), metadata$File_prefix)]

# Step 4: Aggregate norm_counts by Origin
library(dplyr)

# Transpose, bind Origin, then group and summarize
norm_counts_t <- as.data.frame(t(norm_counts))
norm_counts_t$Origin <- origin_vector

aggregated <- norm_counts_t %>%
  group_by(Origin) %>%
  summarise(across(.cols = where(is.numeric), .fns = mean))

# Step 5: Transpose back to original gene-by-sample format
norm_counts_by_origin <- as.data.frame(t(aggregated[,-1]))
colnames(norm_counts_by_origin) <- aggregated$Origin

# plot histogram
library(ggplot2)
library(ggbreak)
```

    ## ggbreak v0.1.4 Learn more at https://yulab-smu.top/

    ## If you use ggbreak in published research, please cite the following
    ## paper:
    ## 
    ## S Xu, M Chen, T Feng, L Zhan, L Zhou, G Yu. Use ggbreak to effectively
    ## utilize plotting space to deal with large datasets and outliers.
    ## Frontiers in Genetics. 2021, 12:774846. doi: 10.3389/fgene.2021.774846

``` r
# Step 1: Reshape your data for plotting (long format)
norm_counts_long <- tidyr::pivot_longer(norm_counts_by_origin, 
                                         cols = c(EU, CA), 
                                         names_to = "Origin", 
                                         values_to = "Count")

# Step 2: Plot histogram
ggplot(norm_counts_long, aes(x = Count, fill = Origin)) +
  geom_histogram(binwidth = 1000, alpha = 0.7, position = "dodge") +
  scale_fill_manual(values = c("EU" = '#5CB7CE', "CA" = '#fb8500')) +
  labs(title = "Histogram of Counts for EU and CA", x = "Count", y = "Frequency") +
  theme_classic() +
  theme(legend.title = element_blank()) + # Optional: Remove legend title
  xlim(0, 20000) +
  scale_y_break(c(400, 1000), scales = 0.5)  # breaks the y-axis between 200 and 19500
```

![](Script3_DESeq2_files/figure-markdown_github/compare_all_gene_expr-1.png)

``` r
# calculate rank sum p value
wilcox.test(Count ~ Origin, data = norm_counts_long, subset = Origin %in% c("CA", "EU"))
```

    ## 
    ##  Wilcoxon rank sum test with continuity correction
    ## 
    ## data:  Count by Origin
    ## W = 32549607, p-value = 0.4534
    ## alternative hypothesis: true location shift is not equal to 0

There is no significant difference in the relative amounts of gene
expression between both EU and CA mushrooms.

# Plots for publication

## Extract DEseq2 results

The next bit of code will generate a new table from DESeq2 results,
keeping only log2FC and padj for each genotype and comparison.

I saved it as a function to mini-fy each DESeq2 output table

group1 or group2 = comparison from deseq you want to make mini

``` r
make_mini <- function(group1,group2,prefix="",keepmean=TRUE) {
  # group1 = "control"
  # group2 = "treatment"
  # prefix = "wildtype" # description you want used in your output minitable
  temp <- as.data.frame(results(dds,
                                contrast=c("Origin", group1, group2),
                                alpha=0.05,
                                , independentFiltering=FALSE, cooksCutoff=TRUE))
  temp$ID <- row.names(temp)
  
  output <- temp[,c("ID","baseMean","log2FoldChange","padj")]
  colnames(output) <- c("ID","baseMean",
                        paste0(prefix,"_L2FC"),
                        paste0(prefix,"_padj"))
  
  # if specified above, do not keep the baseMean column
  if(keepmean==FALSE) {
    output <- output[,c("ID",
                        paste0(prefix,"_L2FC"),
                        paste0(prefix,"_padj"))]
  }
  
  return(output)
}

deseq_results <- make_mini("CA","EU","EU_vs_CA")
deseq_results$ID <- gsub("^[^-]+-|\\.m01$", "", deseq_results$ID )

rmarkdown::paged_table(deseq_results[msdins$X4,], options = NULL)
```

<script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["ID"],"name":[1],"type":["chr"],"align":["left"]},{"label":["baseMean"],"name":[2],"type":["dbl"],"align":["right"]},{"label":["EU_vs_CA_L2FC"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["EU_vs_CA_padj"],"name":[4],"type":["dbl"],"align":["right"]}],"data":[{"1":"NA","2":"NA","3":"NA","4":"NA","_rn_":"__NA__"},{"1":"MLGFLVLP","2":"36535.675733","3":"3.31686642","4":"0.0000022006430009844596697144","_rn_":"MLGFLVLP"},{"1":"NA","2":"NA","3":"NA","4":"NA","_rn_":"NA.1"},{"1":"AWLATCP_1","2":"92.854224","3":"-0.40610347","4":"0.4046918782477719189749620909","_rn_":"AWLATCP_1"},{"1":"NA","2":"NA","3":"NA","4":"NA","_rn_":"NA.2"},{"1":"LIQRPFAP_1","2":"20.668640","3":"0.49033914","4":"0.0500381167695445397414211186","_rn_":"LIQRPFAP_1"},{"1":"LGRPESLP","2":"12.584763","3":"-1.99170674","4":"0.0000020660466072784794207989","_rn_":"LGRPESLP"},{"1":"LRLPPFMIPP","2":"19.783357","3":"-1.67571124","4":"0.0163291508128482094897027110","_rn_":"LRLPPFMIPP"},{"1":"IIGILLPP","2":"7.858281","3":"-3.84483686","4":"0.0023820415569868742637626990","_rn_":"IIGILLPP"},{"1":"AWLVDCP","2":"20.302731","3":"-1.30345978","4":"0.0012478509968585629278520210","_rn_":"AWLVDCP"},{"1":"IWGIGCNP","2":"40.206656","3":"-0.93389685","4":"0.1200156259843014883159639794","_rn_":"IWGIGCNP"},{"1":"IWGIGCDP_1","2":"10.534327","3":"-0.33787391","4":"0.4589252256084378478462326711","_rn_":"IWGIGCDP_1"},{"1":"IWGIGCDP_2","2":"38.222431","3":"1.18189743","4":"0.0177344067669819860366686726","_rn_":"IWGIGCDP_2"},{"1":"NA","2":"NA","3":"NA","4":"NA","_rn_":"NA.3"},{"1":"IFLAFPIPP","2":"8.528001","3":"0.17329788","4":"0.5174349796352358810125338096","_rn_":"IFLAFPIPP"},{"1":"IFWFIYFP","2":"8.942074","3":"-0.71643510","4":"0.0982624255744215802099716939","_rn_":"IFWFIYFP"},{"1":"LIQRPFAP_2","2":"107.061596","3":"0.63604031","4":"0.0440681353566374778285563707","_rn_":"LIQRPFAP_2"},{"1":"FNLFRFPYP","2":"969.661524","3":"-1.42097853","4":"0.0004269446460458155962434945","_rn_":"FNLFRFPYP"},{"1":"VQKPWSRP","2":"4.686841","3":"0.03511039","4":"0.9397560968337383968673748313","_rn_":"VQKPWSRP"},{"1":"IRLPPLFLPP","2":"9.756578","3":"-0.71644859","4":"0.1792088559790198376564518412","_rn_":"IRLPPLFLPP"},{"1":"FFFPPFFIPP","2":"5.580039","3":"-0.56461639","4":"0.3385093128667375772522518673","_rn_":"FFFPPFFIPP"},{"1":"FNILPFMLPP","2":"52.650357","3":"-0.82364245","4":"0.0985105245490786735063082347","_rn_":"FNILPFMLPP"},{"1":"FFPIVFSPP","2":"8.924417","3":"0.13507615","4":"0.7157211245901549112602424429","_rn_":"FFPIVFSPP"},{"1":"HFASFIPP","2":"12.030143","3":"-0.57397882","4":"0.0894781992911917178901148873","_rn_":"HFASFIPP"},{"1":"ISDPTAYP","2":"23.045584","3":"-1.47396718","4":"0.0022608450475629922153675277","_rn_":"ISDPTAYP"},{"1":"LFFWFWFLWP","2":"23.827455","3":"1.67747592","4":"0.0000920016734282405016127904","_rn_":"LFFWFWFLWP"},{"1":"SFFFPIP","2":"33476.783578","3":"4.64881849","4":"0.0000000622953708836750715519","_rn_":"SFFFPIP"},{"1":"FMPLAP","2":"5.640720","3":"-0.45150602","4":"0.2047581883090846543993279738","_rn_":"FMPLAP"},{"1":"LILLAALGIP","2":"349.713793","3":"2.58168622","4":"0.0000000000000000000001065662","_rn_":"LILLAALGIP"},{"1":"GPVFFA","2":"30589.347000","3":"4.20111767","4":"0.0000198672988550490839778154","_rn_":"GPVFFA"},{"1":"NA","2":"NA","3":"NA","4":"NA","_rn_":"NA.4"},{"1":"NA","2":"NA","3":"NA","4":"NA","_rn_":"NA.5"},{"1":"NA","2":"NA","3":"NA","4":"NA","_rn_":"NA.6"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

## Heatmap

``` r
# also make heatmap of genes with custom colors
genotype_colors <- c("#723D00","#fb8500","#ffb703",
                     "#023047","#0C556E","#177995","#219EBC","#5CB7CE","#96CFE0","#D1E8F2")

annotation_colors <- list(Population = setNames(genotype_colors, unique(df$Population)))

# Define your custom color scale (from dark blue to aqua)
my_colors <- colorRampPalette(c("grey90","#03045e"))(100)

# log transform count data
mat_log <- counts(dds, normalized=TRUE)[select,]
mat_log <- log10(mat_log)

top15 <- row.names(mat_log)
colnames(mat_log) <- metadata$Sample_ID

# remove the random prefix and suffixes
row.names(mat_log) <- gsub("^[^-]+-|\\.m01$", "", row.names(mat_log))
row.names(mat_log)
```

    ##  [1] "Ap.00g078360" "Ap.00g054880" "Ap.00g019870" "Ap.00g054640" "Ap.00g057760"
    ##  [6] "MLGFLVLP"     "Ap.00g068810" "Ap.00g039900" "SFFFPIP"      "Ap.00g008070"
    ## [11] "GPVFFA"       "Ap.00g078370" "Ap.00g011790" "Ap.00g001030" "Ap.00g063180"

``` r
# remove log10(0) values that are -Inf
mat_log[mat_log == -Inf] <- 0

library(grid)
```

    ## 
    ## Attaching package: 'grid'

    ## The following object is masked from 'package:Biostrings':
    ## 
    ##     pattern

``` r
heatmap_plot <- pheatmap(mat_log, cluster_rows=FALSE, 
         show_rownames=TRUE,
         cluster_cols=TRUE,
         treeheight_col = 30,
         annotation_col=df,
         annotation_colors = annotation_colors,
         legend=TRUE,
         color=my_colors,
         legend_width = unit(0.2, "in"),
         legend_position = "right",
         border_color = "transparent")
```

![](Script3_DESeq2_files/figure-markdown_github/pretty_heatmap-1.png)

``` r
gtable_obj <- heatmap_plot$gtable

# Modify the legend's vertical position
# Find the legend in the gtable and move it up
gtable_obj$grobs[[which(gtable_obj$layout$name == "legend")]]$grobs[[1]]$gp <- gpar(fontsize = 20)  # Set the font size to 10


# Modify the y position of the annotation_colors legend
# Adjust the y-coordinate to move the legend up
annotation_legend_index <- which(gtable_obj$layout$name == "annotation_legend")
gtable_obj$grobs[[annotation_legend_index]]$grobs[[1]]$gp <- gpar(fontsize = 20)  # Set the font size to 10



# Draw the modified heatmap with legend moved up
pdf(paste0(subDir,"/Plots/heatmap.pdf"), width = 8, height = 5, family="ArialMT")  # Specify desired width and height

grid.draw(gtable_obj)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
# Transpose the matrix
mat_rotated <- t(mat_log)

pdf(paste0(subDir,"/Plots/heatmap_transposed.pdf"), width = 5, height = 8, family="ArialMT")  # Specify desired width and height

heatmap_plot<-pheatmap(mat_rotated, cluster_rows=TRUE, 
         show_rownames=TRUE,
         cluster_cols=FALSE,
         treeheight_col = 30,
         annotation_col=df,
         annotation_colors = annotation_colors,
         legend=TRUE,
         color=my_colors,
         legend_width = unit(0.2, "in"),
         legend_position = "right",
         border_color = "transparent")
gtable_obj <- heatmap_plot$gtable
grid.draw(gtable_obj)
```

![](Script3_DESeq2_files/figure-markdown_github/pretty_heatmap-2.png)

``` r
dev.off()
```

    ## pdf 
    ##   3

## Supplemental figure: log10 counts for MSDINs

``` r
######### FACET_WRAP log10 COUNTS
# extract variance stabilized gene count data for all msdin genes

# list of all msdins
list
```

    ##  [1] "MLGFLVLP"   "AWLATCP_1"  "LIQRPFAP_1" "LGRPESLP"   "LRLPPFMIPP"
    ##  [6] "IIGILLPP"   "AWLVDCP"    "IWGIGCNP"   "IWGIGCDP_1" "IWGIGCDP_2"
    ## [11] "IFLAFPIPP"  "IFWFIYFP"   "LIQRPFAP_2" "FNLFRFPYP"  "VQKPWSRP"  
    ## [16] "IRLPPLFLPP" "FFFPPFFIPP" "FNILPFMLPP" "FFPIVFSPP"  "HFASFIPP"  
    ## [21] "ISDPTAYP"   "LFFWFWFLWP" "SFFFPIP"    "FMPLAP"     "LILLAALGIP"
    ## [26] "GPVFFA"

``` r
missing_genes <- setdiff(list, rownames(dds))
missing_genes
```

    ## character(0)

``` r
list <- intersect(list, rownames(dds))


# generate facet_wrap for these samples

# extract data
mat2 <- as.data.frame(counts(dds,normalized=TRUE)[list,])
# log transform
mat2 <- log10(mat2)
# remove log10(0) values that are -Inf
mat2[mat2 == -Inf] <- 0
colnames(mat2) <- metadata$Sample_ID

# Reshape the data from wide to long format
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ lubridate 1.9.4     ✔ tibble    3.3.0
    ## ✔ purrr     1.0.4     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ lubridate::%within%()   masks IRanges::%within%()
    ## ✖ Biostrings::collapse()  masks dplyr::collapse(), IRanges::collapse()
    ## ✖ dplyr::combine()        masks Biobase::combine(), BiocGenerics::combine()
    ## ✖ purrr::compact()        masks XVector::compact()
    ## ✖ dplyr::count()          masks matrixStats::count()
    ## ✖ dplyr::desc()           masks IRanges::desc()
    ## ✖ tidyr::expand()         masks S4Vectors::expand()
    ## ✖ dplyr::filter()         masks stats::filter()
    ## ✖ dplyr::first()          masks S4Vectors::first()
    ## ✖ dplyr::lag()            masks stats::lag()
    ## ✖ ggplot2::Position()     masks BiocGenerics::Position(), base::Position()
    ## ✖ purrr::reduce()         masks GenomicRanges::reduce(), IRanges::reduce()
    ## ✖ dplyr::rename()         masks S4Vectors::rename()
    ## ✖ lubridate::second()     masks S4Vectors::second()
    ## ✖ lubridate::second<-()   masks S4Vectors::second<-()
    ## ✖ AnnotationDbi::select() masks dplyr::select()
    ## ✖ XVector::slice()        masks dplyr::slice(), IRanges::slice()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
mat2_long <- mat2 %>%
  rownames_to_column("MSDIN") %>%        # Convert row names into a column
  pivot_longer(cols = -MSDIN,            # Pivot columns except for "Row"
               names_to = "Sample_ID",    # New column for column names
               values_to = "Value")    # New column for values

# add a new column to facet over:
# msdin type
msdin_type <- c("Known","Known SFFFPIP","Leaderless","Prolineless")
noncanonical <- c("MLGFLVLP","GPVFFA")

# Add the "Type" column
mat2_long$Type <- ifelse(mat2_long$MSDIN %in% noncanonical, "Noncanonical", "Known")
mat2_long$Type <- ifelse(mat2_long$MSDIN %in% "SFFFPIP", "Known SFFFPIP",mat2_long$Type)
mat2_long$Type <- ifelse(mat2_long$MSDIN %in% "MLGFLVLP", "Leaderless",mat2_long$Type)
mat2_long$Type <- ifelse(mat2_long$MSDIN %in% "GPVFFA", "Prolineless",mat2_long$Type)

# add additional metadata:
# Population
metadata2 <- metadata[,c("Sample_ID","Population","Origin")]
mat2_long <- mat2_long %>%
  left_join(metadata2, by = "Sample_ID")

# set order
mat2_long$Origin <- factor(mat2_long$Origin, levels = c("CA",
                                                    "EU"))

mat2_long$Population <- factor(mat2_long$Population, levels = c("Drake",
                                                            "Pet",
                                                            "Picnic",
                                                            "Aso",
                                                            "Champ",
                                                            "Doc",
                                                            "Gron",
                                                            "Lamp",
                                                            "Quill",
                                                            "Tisza"))

# extract p values for each MSDIN
deseq_msdin <- deseq_results[list,]

subset(deseq_msdin,EU_vs_CA_L2FC > 1 & EU_vs_CA_padj < 0.05)[,1]
```

    ## [1] "MLGFLVLP"   "IWGIGCDP_2" "LFFWFWFLWP" "SFFFPIP"    "LILLAALGIP"
    ## [6] "GPVFFA"

``` r
subset(deseq_msdin,EU_vs_CA_L2FC < -1 & EU_vs_CA_padj < 0.05)[,1]
```

    ## [1] "LGRPESLP"   "LRLPPFMIPP" "IIGILLPP"   "AWLVDCP"    "FNLFRFPYP" 
    ## [6] "ISDPTAYP"

``` r
# round p values
options(scipen = 999)

deseq_msdin$EU_vs_CA_padj <- signif(deseq_msdin$EU_vs_CA_padj, 3)
deseq_msdin$EU_vs_CA_padj <- formatC(deseq_msdin$EU_vs_CA_padj, format = "g", digits = 3)

colnames(deseq_msdin)[1] <-"MSDIN"
# merge this onto mat2_long
mat2_long <- merge(mat2_long,deseq_msdin,by="MSDIN",all.x=TRUE)
mat2_long$EU_vs_CA_padj <- as.numeric(mat2_long$EU_vs_CA_padj)

# now plot it and facet by MSDIN
pdf(paste0(subDir,"/Plots/log10_msdin.pdf"), width = 12, height = 12, family="ArialMT")  # Specify desired width and height

ggplot(data=mat2_long,aes(x=Origin,y=Value)) +
  geom_point(aes(colour = Population), size = 4, position = position_jitter(w = 0.3, h = 0)) +
  # add mean line
  stat_summary(fun= mean, fun.min=mean, fun.max=mean, geom="crossbar", width=0.5, linewidth=1,color="black") +
  # add error bars
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), geom="errorbar", width=0.2, linewidth=1,color="black") +
  theme_classic()  +
  scale_color_manual(values = genotype_colors) +
  theme(legend.position = "none",
        text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  ylab("Log10(normalized count)") +
  facet_wrap(Type~MSDIN)+
  # Add p-value labels on each facet
  geom_text(aes(x = 1.5, y = max(Value) * 0.9, label = sprintf("p = %.3f", EU_vs_CA_padj)), size = 5, color = "black", inherit.aes = FALSE)

dev.off()
```

    ## quartz_off_screen 
    ##                 2

## Another version; Normalized counts for each MSDIN

``` r
# also do this with normalized counts
# extract data
mat3 <- as.data.frame(counts[list,])
colnames(mat3) <- metadata$Sample_ID

# Reshape the data from wide to long format
mat3_long <- mat3 %>%
  rownames_to_column("MSDIN") %>%        # Convert row names into a column
  pivot_longer(cols = -MSDIN,            # Pivot columns except for "Row"
               names_to = "Sample_ID",    # New column for column names
               values_to = "Value")    # New column for values

# add a new column to facet over:
# msdin type
msdin_type <- c("canonical","new","leaderless","prolineless")
noncanonical <- c("MLGFLVLP","GPVFFA")

# Add the "Type" column
mat3_long$Type <- ifelse(mat3_long$MSDIN %in% noncanonical, "noncanonical", "canonical")
mat3_long$Type <- ifelse(mat3_long$MSDIN %in% "SFFFPIP", "new",mat3_long$Type)
mat3_long$Type <- ifelse(mat3_long$MSDIN %in% "MLGFLVLP", "leaderless",mat3_long$Type)
mat3_long$Type <- ifelse(mat3_long$MSDIN %in% "GPVFFA", "prolineless",mat3_long$Type)

# add additional metadata:
# Population
metadata2 <- metadata[,c("Sample_ID","Population","Origin")]
mat3_long <- mat3_long %>%
  left_join(metadata2, by = "Sample_ID")

# set order
mat3_long$Origin <- factor(mat3_long$Origin, levels = c("CA",
                                                    "EU"))

mat3_long$Population <- factor(mat3_long$Population, levels = c("Drake",
                                                            "Pet",
                                                            "Picnic",
                                                            "Aso",
                                                            "Champ",
                                                            "Doc",
                                                            "Gron",
                                                            "Lamp",
                                                            "Quill",
                                                            "Tisza"))



# now plot it and facet by MSDIN
pdf(paste0(subDir,"/Plots/normcount_msdin.pdf"), width = 12, height = 12, family="ArialMT")  # Specify desired width and height

ggplot(data=mat3_long,aes(x=Origin,y=Value)) +
  geom_point(aes(colour = Population), size = 4, position = position_jitter(w = 0.3, h = 0)) +
  # add mean line
  stat_summary(fun= mean, fun.min=mean, fun.max=mean, geom="crossbar", width=0.5, linewidth=1,color="black") +
  # add error bars
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), geom="errorbar", width=0.2, linewidth=1,color="black") +
  theme_classic()  +
  scale_color_manual(values = genotype_colors) +
  theme(legend.position = "none",
        text = element_text(size = 16),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  ylab("Normalized count") +
  facet_wrap(~MSDIN, scales = "free_y")

dev.off()
```

    ## quartz_off_screen 
    ##                 2

## Violin plots and statistical tests

``` r
# adjust type names for the graph
mat4_long<- mat2_long %>%
  mutate(Type = recode(Type, 
                       "Known" = "Known MSDINs",
                       "Known SFFFPIP" = "Known SFFFPIP",
                       "Prolineless" = "Prolineless GPVFFA", 
                       "Leaderless" = "Leaderless MLGFLVLP"))

genotype_colors <- c("#723D00","#fb8500","#ffb703",
                     "#023047","#0C556E","#177995","#219EBC","#5CB7CE","#96CFE0","#D1E8F2")


# set order
mat4_long$Type <- factor(mat4_long$Type, levels = c("Known MSDINs",
                                                    "Known SFFFPIP",
                                                    "Leaderless MLGFLVLP",
                                                    "Prolineless GPVFFA"))

pdf(paste0(subDir,"/Plots/violin.pdf"), width = 4.5, height = 2.5, family="ArialMT")

ggplot(data=mat4_long,aes(x=Origin,y=Value,fill=Origin)) +
geom_violin(trim = FALSE, adjust = 0.6, draw_quantiles = c(0.25, 0.75), linetype = "dashed") +
  geom_violin(trim = FALSE, adjust = 0.6, draw_quantiles = 0.5, fill="transparent") +
  theme_classic()  +
  scale_fill_manual(values = c("#fb8500", "#0C556E")) +
  ylab("Log10(normalized count)") +
    coord_cartesian(ylim = c(-0.5, 7)) +
  #geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 1) +  # Black line at y=0
  facet_wrap(~Type, ncol=4) +
  theme(panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "none",
        text = element_text(color = "black",size = 10))

dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
##### determine if normal distribution

# Group-wise Shapiro-Wilk test
normality_tests <- mat4_long %>%
  group_by(Origin, Type) %>%
  summarize(
    n = n(),
    mean = mean(Value, na.rm = TRUE),
    sem = sd(Value, na.rm = TRUE) / sqrt(n),
    shapiro_p = shapiro.test(Value)$p.value,
    is_normal = shapiro_p > 0.05
  )
```

    ## `summarise()` has grouped output by 'Origin'. You can override using the
    ## `.groups` argument.

``` r
print(normality_tests)
```

    ## # A tibble: 8 × 7
    ## # Groups:   Origin [2]
    ##   Origin Type                    n  mean    sem shapiro_p is_normal
    ##   <fct>  <fct>               <int> <dbl>  <dbl>     <dbl> <lgl>    
    ## 1 CA     Known MSDINs          437  1.14 0.0349  6.94e-10 FALSE    
    ## 2 CA     Known SFFFPIP          19  4.60 0.191   1.37e- 7 FALSE    
    ## 3 CA     Leaderless MLGFLVLP    19  4.44 0.201   3.96e- 4 FALSE    
    ## 4 CA     Prolineless GPVFFA     19  4.55 0.173   1.23e- 6 FALSE    
    ## 5 EU     Known MSDINs          414  1.18 0.0347  1.16e- 9 FALSE    
    ## 6 EU     Known SFFFPIP          18  1.95 0.264   2.63e- 3 FALSE    
    ## 7 EU     Leaderless MLGFLVLP    18  3.00 0.233   3.66e- 1 TRUE     
    ## 8 EU     Prolineless GPVFFA     18  1.54 0.306   1.56e- 4 FALSE

``` r
#### also plot known and new/noncanonical separately
sub_canon <- subset(mat4_long,Type=="Known MSDINs")
pdf(paste0(subDir,"/Plots/violin_canon.pdf"), width = 1.5, height = 2.5, family="ArialMT")

ggplot(data=sub_canon,aes(x=Origin,y=Value,fill=Origin)) +
geom_violin(trim = FALSE, adjust = 0.6, draw_quantiles = c(0.25, 0.75), linetype = "dashed") +
  geom_violin(trim = FALSE, adjust = 0.6, draw_quantiles = 0.5, fill="transparent") +
  theme_classic()  +
  scale_fill_manual(values = c("#fb8500", "#0C556E")) +
  ylab("Log10(normalized count)") +
    coord_cartesian(ylim = c(-0.5, 7)) +
  #geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 1) +  # Black line at y=0
  facet_wrap(~Type, ncol=4) +
  theme(panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "none",
        text = element_text(color = "black",size = 10))

dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
sub_noncanon <- subset(mat4_long,Type!="Known MSDINs")
pdf(paste0(subDir,"/Plots/violin_noncanon.pdf"), width = 3.5, height = 2.5, family="ArialMT")

ggplot(data=sub_noncanon,aes(x=Origin,y=Value,fill=Origin)) +
  geom_violin(trim = FALSE, adjust = 0.6, draw_quantiles = c(0.25, 0.75), linetype = "dashed") +
  geom_violin(trim = FALSE, adjust = 0.6, draw_quantiles = 0.5, fill="transparent") +
  theme_classic()  +
  scale_fill_manual(values = c("#fb8500", "#0C556E")) +
  ylab("Log10(normalized count)") +
  coord_cartesian(ylim = c(-0.5, 7)) +
  #geom_hline(yintercept = 0, color = "black", linetype = "solid", size = 1) +  # Black line at y=0
  facet_wrap(~Type, ncol=4) +
  theme(panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "none",
        text = element_text(color = "black",size = 10))

dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
library(dplyr)

# Run ANOVA for each Type
mat4_long$Var <- paste0(mat4_long$Origin,"_",mat4_long$Type)
# get log10:
# apply pseudocount for zeros
mat4_long$Value_log10 <- log10(mat4_long$Value +1)
df <- mat4_long[,c("Var","Value_log10")]

# make a plot of values
df$Var<-as.factor(df$Var)

df%>%ggplot(aes(x=Value_log10, fill=Var))+geom_density(alpha=0.5)
```

![](Script3_DESeq2_files/figure-markdown_github/pretty_violin-1.png)

``` r
# anova
model1<-lm(Value_log10~Var, data=df)
anova(model1)
```

    ## Analysis of Variance Table
    ## 
    ## Response: Value_log10
    ##            Df Sum Sq Mean Sq F value                Pr(>F)    
    ## Var         7 11.036 1.57662  74.285 < 0.00000000000000022 ***
    ## Residuals 954 20.248 0.02122                                  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Print results
# Tukey multiple comparisons
library(multcomp)
```

    ## Loading required package: mvtnorm
    ## Loading required package: survival
    ## Loading required package: TH.data
    ## Loading required package: MASS
    ## 
    ## Attaching package: 'MASS'
    ## 
    ## The following object is masked from 'package:AnnotationDbi':
    ## 
    ##     select
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select
    ## 
    ## 
    ## Attaching package: 'TH.data'
    ## 
    ## The following object is masked from 'package:MASS':
    ## 
    ##     geyser

``` r
anova_results<-summary(glht(model1, mcp(Var="Tukey")))
```

    ## Warning in RET$pfunction("adjusted", ...): Completion with error > abseps
    ## Warning in RET$pfunction("adjusted", ...): Completion with error > abseps

``` r
anova_results
```

    ## 
    ##   Simultaneous Tests for General Linear Hypotheses
    ## 
    ## Multiple Comparisons of Means: Tukey Contrasts
    ## 
    ## 
    ## Fit: lm(formula = Value_log10 ~ Var, data = df)
    ## 
    ## Linear Hypotheses:
    ##                                                       Estimate Std. Error
    ## CA_Known SFFFPIP - CA_Known MSDINs == 0               0.436471   0.034141
    ## CA_Leaderless MLGFLVLP - CA_Known MSDINs == 0         0.424007   0.034141
    ## CA_Prolineless GPVFFA - CA_Known MSDINs == 0          0.434134   0.034141
    ## EU_Known MSDINs - CA_Known MSDINs == 0                0.011693   0.009992
    ## EU_Known SFFFPIP - CA_Known MSDINs == 0               0.140664   0.035038
    ## EU_Leaderless MLGFLVLP - CA_Known MSDINs == 0         0.284605   0.035038
    ## EU_Prolineless GPVFFA - CA_Known MSDINs == 0          0.058656   0.035038
    ## CA_Leaderless MLGFLVLP - CA_Known SFFFPIP == 0       -0.012465   0.047266
    ## CA_Prolineless GPVFFA - CA_Known SFFFPIP == 0        -0.002337   0.047266
    ## EU_Known MSDINs - CA_Known SFFFPIP == 0              -0.424779   0.034181
    ## EU_Known SFFFPIP - CA_Known SFFFPIP == 0             -0.295808   0.047918
    ## EU_Leaderless MLGFLVLP - CA_Known SFFFPIP == 0       -0.151867   0.047918
    ## EU_Prolineless GPVFFA - CA_Known SFFFPIP == 0        -0.377815   0.047918
    ## CA_Prolineless GPVFFA - CA_Leaderless MLGFLVLP == 0   0.010127   0.047266
    ## EU_Known MSDINs - CA_Leaderless MLGFLVLP == 0        -0.412314   0.034181
    ## EU_Known SFFFPIP - CA_Leaderless MLGFLVLP == 0       -0.283343   0.047918
    ## EU_Leaderless MLGFLVLP - CA_Leaderless MLGFLVLP == 0 -0.139402   0.047918
    ## EU_Prolineless GPVFFA - CA_Leaderless MLGFLVLP == 0  -0.365351   0.047918
    ## EU_Known MSDINs - CA_Prolineless GPVFFA == 0         -0.422442   0.034181
    ## EU_Known SFFFPIP - CA_Prolineless GPVFFA == 0        -0.293471   0.047918
    ## EU_Leaderless MLGFLVLP - CA_Prolineless GPVFFA == 0  -0.149529   0.047918
    ## EU_Prolineless GPVFFA - CA_Prolineless GPVFFA == 0   -0.375478   0.047918
    ## EU_Known SFFFPIP - EU_Known MSDINs == 0               0.128971   0.035077
    ## EU_Leaderless MLGFLVLP - EU_Known MSDINs == 0         0.272912   0.035077
    ## EU_Prolineless GPVFFA - EU_Known MSDINs == 0          0.046964   0.035077
    ## EU_Leaderless MLGFLVLP - EU_Known SFFFPIP == 0        0.143941   0.048561
    ## EU_Prolineless GPVFFA - EU_Known SFFFPIP == 0        -0.082007   0.048561
    ## EU_Prolineless GPVFFA - EU_Leaderless MLGFLVLP == 0  -0.225949   0.048561
    ##                                                      t value Pr(>|t|)    
    ## CA_Known SFFFPIP - CA_Known MSDINs == 0               12.784  < 0.001 ***
    ## CA_Leaderless MLGFLVLP - CA_Known MSDINs == 0         12.419  < 0.001 ***
    ## CA_Prolineless GPVFFA - CA_Known MSDINs == 0          12.716  < 0.001 ***
    ## EU_Known MSDINs - CA_Known MSDINs == 0                 1.170  0.92471    
    ## EU_Known SFFFPIP - CA_Known MSDINs == 0                4.015  0.00136 ** 
    ## EU_Leaderless MLGFLVLP - CA_Known MSDINs == 0          8.123  < 0.001 ***
    ## EU_Prolineless GPVFFA - CA_Known MSDINs == 0           1.674  0.66020    
    ## CA_Leaderless MLGFLVLP - CA_Known SFFFPIP == 0        -0.264  0.99999    
    ## CA_Prolineless GPVFFA - CA_Known SFFFPIP == 0         -0.049  1.00000    
    ## EU_Known MSDINs - CA_Known SFFFPIP == 0              -12.427  < 0.001 ***
    ## EU_Known SFFFPIP - CA_Known SFFFPIP == 0              -6.173  < 0.001 ***
    ## EU_Leaderless MLGFLVLP - CA_Known SFFFPIP == 0        -3.169  0.02764 *  
    ## EU_Prolineless GPVFFA - CA_Known SFFFPIP == 0         -7.885  < 0.001 ***
    ## CA_Prolineless GPVFFA - CA_Leaderless MLGFLVLP == 0    0.214  1.00000    
    ## EU_Known MSDINs - CA_Leaderless MLGFLVLP == 0        -12.063  < 0.001 ***
    ## EU_Known SFFFPIP - CA_Leaderless MLGFLVLP == 0        -5.913  < 0.001 ***
    ## EU_Leaderless MLGFLVLP - CA_Leaderless MLGFLVLP == 0  -2.909  0.05977 .  
    ## EU_Prolineless GPVFFA - CA_Leaderless MLGFLVLP == 0   -7.624  < 0.001 ***
    ## EU_Known MSDINs - CA_Prolineless GPVFFA == 0         -12.359  < 0.001 ***
    ## EU_Known SFFFPIP - CA_Prolineless GPVFFA == 0         -6.124  < 0.001 ***
    ## EU_Leaderless MLGFLVLP - CA_Prolineless GPVFFA == 0   -3.121  0.03226 *  
    ## EU_Prolineless GPVFFA - CA_Prolineless GPVFFA == 0    -7.836  < 0.001 ***
    ## EU_Known SFFFPIP - EU_Known MSDINs == 0                3.677  0.00502 ** 
    ## EU_Leaderless MLGFLVLP - EU_Known MSDINs == 0          7.780  < 0.001 ***
    ## EU_Prolineless GPVFFA - EU_Known MSDINs == 0           1.339  0.85840    
    ## EU_Leaderless MLGFLVLP - EU_Known SFFFPIP == 0         2.964  0.05185 .  
    ## EU_Prolineless GPVFFA - EU_Known SFFFPIP == 0         -1.689  0.64998    
    ## EU_Prolineless GPVFFA - EU_Leaderless MLGFLVLP == 0   -4.653  < 0.001 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## (Adjusted p values reported -- single-step method)

``` r
# Extract coefficients and p-values into a data frame
tukey_df <- as.data.frame(anova_results$test[c("coefficients", "sigma", "tstat", "pvalues")])

# Add row names as a column for pairwise comparison names
tukey_df$comparison <- rownames(tukey_df)
rownames(tukey_df) <- NULL

# Reorder columns to have comparison first
tukey_df <- tukey_df[, c("comparison", "coefficients", "sigma", "tstat", "pvalues")]
write.csv(tukey_df, paste0(subDir,"/violin_anova_results.csv"),row.names=FALSE)


# List of comparisons
comparisons <- list(
  c("CA_Leaderless MLGFLVLP", "CA_Known MSDINs"),
  c("CA_Prolineless GPVFFA", "CA_Known MSDINs"),
  c("CA_Known SFFFPIP", "CA_Known MSDINs"),
  
  c("EU_Known MSDINs", "CA_Known MSDINs"),
  c("EU_Known SFFFPIP", "CA_Known SFFFPIP"),
  c("EU_Leaderless MLGFLVLP", "CA_Leaderless MLGFLVLP"),
  c("EU_Prolineless GPVFFA", "CA_Prolineless GPVFFA"),
  
  c("EU_Leaderless MLGFLVLP", "EU_Known MSDINs"),
  c("EU_Prolineless GPVFFA", "EU_Known MSDINs"),
  c("EU_Known SFFFPIP", "EU_Known MSDINs"),
  
  c("EU_Prolineless GPVFFA", "EU_Leaderless MLGFLVLP")
)

# Function to run t-test
run_t_test <- function(group1, group2) {
  # Subset the data for each group
  group1_data <- df[df$Var == group1, "Value_log10"]
  group2_data <- df[df$Var == group2, "Value_log10"]
  
  # Perform t-test
  test_result <- t.test(group1_data, group2_data)
  
  # Return the result (comparison name, p-value, and test statistics)
  return(data.frame(
    Comparison = paste(group1, "vs", group2),
    p_value = test_result$p.value,
    statistic = test_result$statistic,
    conf_low = test_result$conf.int[1],
    conf_high = test_result$conf.int[2]
  ))
}

# Apply the function to all comparisons and combine the results
t_test_results <- do.call(rbind, lapply(comparisons, function(x) run_t_test(x[1], x[2])))

# Write the results to a CSV file
write.csv(t_test_results, paste0(subDir,"/violin_t_test_results.csv"), row.names = FALSE)

# Display the results
print(t_test_results)
```

    ##                                           Comparison                   p_value
    ## t          CA_Leaderless MLGFLVLP vs CA_Known MSDINs 0.00000000000000089643056
    ## t1          CA_Prolineless GPVFFA vs CA_Known MSDINs 0.00000000000000002401309
    ## t2               CA_Known SFFFPIP vs CA_Known MSDINs 0.00000000000000353749137
    ## t3                EU_Known MSDINs vs CA_Known MSDINs 0.24978689322915040293971
    ## t4              EU_Known SFFFPIP vs CA_Known SFFFPIP 0.00000006145031541353583
    ## t5  EU_Leaderless MLGFLVLP vs CA_Leaderless MLGFLVLP 0.00021149978053678548961
    ## t6    EU_Prolineless GPVFFA vs CA_Prolineless GPVFFA 0.00000006686394210939543
    ## t7         EU_Leaderless MLGFLVLP vs EU_Known MSDINs 0.00000000429768457139605
    ## t8          EU_Prolineless GPVFFA vs EU_Known MSDINs 0.30901525839360577396775
    ## t9               EU_Known SFFFPIP vs EU_Known MSDINs 0.00179213718340842913940
    ## t10  EU_Prolineless GPVFFA vs EU_Leaderless MLGFLVLP 0.00015310780722544751965
    ##     statistic    conf_low   conf_high
    ## t   19.623348  0.37928647  0.46872717
    ## t1  21.924323  0.39325862  0.47500956
    ## t2  18.965621  0.38876041  0.48418246
    ## t3   1.151656 -0.00823501  0.03162007
    ## t4  -7.218920 -0.37962891 -0.21198684
    ## t5  -4.172728 -0.20742178 -0.07138226
    ## t6  -7.827818 -0.47476350 -0.27619236
    ## t7   9.950150  0.21561505  0.33020948
    ## t8   1.047048 -0.04730971  0.14123696
    ## t9   3.644155  0.05475250  0.20318956
    ## t10 -4.379130 -0.33167715 -0.12022012

``` r
# now do mann whitney u
run_mann_whitney <- function(group1, group2) {
  # Subset the data for each group
  group1_data <- df[df$Var == group1, "Value_log10"]
  group2_data <- df[df$Var == group2, "Value_log10"]
  
  # Perform Wilcoxon rank-sum test (Mann–Whitney U test)
  test_result <- wilcox.test(group1_data, group2_data, exact = FALSE)
  
  # Return the result (comparison name, p-value, and W statistic)
  return(data.frame(
    Comparison = paste(group1, "vs", group2),
    p_value = test_result$p.value,
    statistic = test_result$statistic
  ))
}

mann_whitney_results <- do.call(rbind, lapply(comparisons, function(comp) {
  run_mann_whitney(comp[1], comp[2])
}))
mann_whitney_results$p_value <- format(mann_whitney_results$p_value, scientific = TRUE)

print(mann_whitney_results)
```

    ##                                           Comparison      p_value statistic
    ## W          CA_Leaderless MLGFLVLP vs CA_Known MSDINs 6.535147e-13    8194.0
    ## W1          CA_Prolineless GPVFFA vs CA_Known MSDINs 5.887915e-13    8202.0
    ## W2               CA_Known SFFFPIP vs CA_Known MSDINs 1.305982e-12    8140.5
    ## W3                EU_Known MSDINs vs CA_Known MSDINs 1.048072e-01   96272.0
    ## W4              EU_Known SFFFPIP vs CA_Known SFFFPIP 1.701975e-06      13.0
    ## W5  EU_Leaderless MLGFLVLP vs CA_Leaderless MLGFLVLP 6.445666e-05      39.0
    ## W6    EU_Prolineless GPVFFA vs CA_Prolineless GPVFFA 5.773672e-07       6.0
    ## W7         EU_Leaderless MLGFLVLP vs EU_Known MSDINs 2.468038e-10    7008.0
    ## W8          EU_Prolineless GPVFFA vs EU_Known MSDINs 7.234004e-01    3910.0
    ## W9               EU_Known SFFFPIP vs EU_Known MSDINs 1.792991e-03    5345.5
    ## W10  EU_Prolineless GPVFFA vs EU_Leaderless MLGFLVLP 3.717844e-04      49.0

## Supplemental volcano plot

``` r
group1 = "CA"
group2 = "EU"
gene_list = msdins$X4 # list of MSDINs
print(paste0(group1," vs ",group2))
```

    ## [1] "CA vs EU"

``` r
res <- results(dds, contrast=c("Origin", group1, group2), alpha=0.05, independentFiltering=FALSE, cooksCutoff=TRUE)

# now save deseq results as a dataframe
deseqoutput <- as.data.frame(res)

# anything past -5 or 5 log2FoldChange, clip to -5 or 5 respectively
deseqoutput$log2FC_clip <- pmax(pmin(deseqoutput$log2FoldChange, 5), -5)

# also clip high pvals at 20
deseqoutput$padj_clip <- pmax(deseqoutput$padj, 0.00000000000000000001)

#c("#fb8500", "#5CB7CE")
# Create a custom color vector
keyvals <- ifelse(deseqoutput$log2FoldChange > 1 & deseqoutput$padj < 0.05, '#fb8500',
                  ifelse(deseqoutput$log2FoldChange < -1 & deseqoutput$padj < 0.05, '#5CB7CE', 'grey70'))

names(keyvals) <- rownames(deseqoutput)

deseqoutput$neglog10_padj <- -log10(deseqoutput$padj_clip)

label_data <- subset(deseqoutput, rownames(deseqoutput) %in% gene_list)

# plot it
library(ggplot2)

pdf(paste0(subDir,"/Plots/volcano.pdf"), width = 7, height = 5, family="ArialMT")  # Specify desired width and height

ggplot(deseqoutput, aes(x = log2FC_clip, y = neglog10_padj)) +
  geom_point(aes(color = keyvals), size = 2, alpha=0.5) +
  scale_color_identity() +  # uses the colors directly from keyvals
  geom_text_repel(data = label_data,
                  aes(label = rownames(label_data)),
                  size = 4,
                  max.overlaps = 10,
                  box.padding = 0.3,       # increase label box size
                  point.padding = 0.4,     # increase spacing between point and label
                  segment.color = 'grey50',
                  nudge_y = 0.5,           # move labels upward
                  force = 2) +                # stronger repulsion force
labs(
  #title = paste0(group1, " vs. ", group2),
  subtitle = "",  # italic("") is effectively blank
  caption = "",
  x = "Log2(Fold Change)",
  y = expression(-Log[10] ~ italic("Adj. P value"))
) +
  geom_point(data = label_data, color = "black", size = 2) +
  
  coord_cartesian(xlim = c(-5, 5), ylim = c(0, 20)) +
  scale_x_continuous(breaks = seq(-5, 5, 1)) +
  theme_classic() +
  theme(legend.position = "bottom")
```

    ## Warning: ggrepel: 15 unlabeled data points (too many overlaps). Consider
    ## increasing max.overlaps

``` r
dev.off()
```

    ## quartz_off_screen 
    ##                 2

``` r
# ggsave(paste0(subDir,"/Plots/",group1,"_vs_",group2,"_pretty_volcano.png"),
#        units ="in",
#        width = 8,
#        height = 6)

summary(res)
```

    ## 
    ## out of 8095 with nonzero total read count
    ## adjusted p-value < 0.1
    ## LFC > 0 (up)       : 2490, 31%
    ## LFC < 0 (down)     : 2842, 35%
    ## outliers [1]       : 0, 0%
    ## low counts [2]     : 0, 0%
    ## (mean count < 0)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

## Main principal component analysis figure (PCA)

``` r
##### PCA PLOT

pdf(paste0(subDir,"/Plots/pca.pdf"), width = 4.5, height = 2.5, family="ArialMT")  # Specify desired width and height

ggplot(pcaData,aes(x=PC1, y=PC2, color = Population, shape=Origin)) +
  geom_point(size=3) + 
  theme_classic() +
  scale_color_manual(values = genotype_colors) +
stat_ellipse(aes(group = Origin), colour = "black")+                 
  xlab(paste0("PC1 (",pc1var,"% variance)")) +
  ylab(paste0("PC2 (",pc2var,"% variance)")) +
  theme(panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        text = element_text(color = "black",size = 10),
        legend.key.height = unit(0.15, 'in'))+ 
  guides(color = guide_legend(override.aes=list(shape = c(rep(16,3),rep(17,7)))))  # Customize shape of color legend


dev.off()
```

    ## quartz_off_screen 
    ##                 2

## MSDIN expression vs mapping rate

``` r
# finally: also plot mapping rate for samples vs. counts for leaderless
dat <- read.csv(paste0(subDir,"/NormCount-Graphs/Data/MLGFLVLP.csv"))
colnames(dat)[1] <- "File_prefix"
# reread metadata
metadata <- read.csv(paste0(subDir,"/metadata_plus_mapping.csv"))

mini_meta <- metadata[,c("File_prefix","perc_mapping")]
mini_meta2 <- metadata[,c("Sample_ID","perc_mapping")]

# merge
dat2 <- merge(dat,mini_meta,by="File_prefix",all.x=TRUE)

dat2$Origin <- factor(dat2$Origin, levels = c("CA",
                                                    "EU"))

dat2$Population <- factor(dat2$Population, levels = c("Drake",
                                                            "Pet",
                                                            "Picnic",
                                                            "Aso",
                                                            "Champ",
                                                            "Doc",
                                                            "Gron",
                                                            "Lamp",
                                                            "Quill",
                                                            "Tisza"))


library(ggplot2)

# Fit the linear model
# for all points:
lm_model <- lm(count ~ perc_mapping, data = dat2)

# based on origin:
ca <- subset(dat2,Origin == "CA")
eu <- subset(dat2,Origin == "EU")

lm_model_ca <- lm(count ~ perc_mapping, data = ca)
lm_model_eu <- lm(count ~ perc_mapping, data = eu)

eqn_ca <- substitute(italic(y) == a + b %.% italic(x) * "," ~~ italic(R)^2 ~ "=" ~ r2, 
                  list(a = as.numeric(format(coef(lm_model_ca)[1], digits = 3)),
                       b = as.numeric(format(coef(lm_model_ca)[2], digits = 3)),
                       r2 = format(summary(lm_model_ca)$r.squared, digits = 3)))

eqn_eu <- substitute(italic(y) == a + b %.% italic(x) * "," ~~ italic(R)^2 ~ "=" ~ r2, 
                  list(a = as.numeric(format(coef(lm_model_eu)[1], digits = 3)),
                       b = as.numeric(format(coef(lm_model_eu)[2], digits = 3)),
                       r2 = format(summary(lm_model_eu)$r.squared, digits = 3)))

# Extract coefficients and R-squared
eqn <- substitute(italic(y) == a + b %.% italic(x) * "," ~~ italic(R)^2 ~ "=" ~ r2, 
                  list(a = as.numeric(format(coef(lm_model)[1], digits = 3)),
                       b = as.numeric(format(coef(lm_model)[2], digits = 3)),
                       r2 = format(summary(lm_model)$r.squared, digits = 3)))
genotype_colors <- c("#723D00","#fb8500","#ffb703",
                     "#023047","#0C556E","#177995","#219EBC","#5CB7CE","#96CFE0","#D1E8F2")

ggplot(dat2, aes(x=perc_mapping,y=count,color=Population, shape=Origin)) +
  geom_point(size=3) +
  scale_color_manual(values = genotype_colors) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add linear regression line
  annotate("text", x = min(dat2$perc_mapping), y = max(dat2$count), 
           label = as.character(as.expression(eqn)), parse = TRUE, hjust = -0.3) +  # Add equation
  theme_classic()

ggsave("mapping_vs_leaderless_counts.png",
       path = paste0(subDir,"/Plots/"),
       units="in",
       width = 6, height = 5)

ggplot(dat2, aes(x=perc_mapping,y=count,color=Population, shape=Origin)) +
  geom_point(size=3) +
  scale_color_manual(values = genotype_colors) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add linear regression line
  annotate("text", x = min(dat2$perc_mapping), y = max(dat2$count), 
           label = as.character(as.expression(eqn_eu)), parse = TRUE, hjust = -0.3) +  # Add equation
  theme_classic()
ggsave("mapping_vs_leaderless_counts_EUlm.png",
       path = paste0(subDir,"/Plots/"),
       units="in",
       width = 6, height = 5)


ggplot(dat2, aes(x=perc_mapping,y=count,color=Population, shape=Origin)) +
  geom_point(size=3) +
  scale_color_manual(values = genotype_colors) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add linear regression line
  annotate("text", x = min(dat2$perc_mapping), y = max(dat2$count), 
           label = as.character(as.expression(eqn_ca)), parse = TRUE, hjust = -0.3) +  # Add equation
  theme_classic()

ggsave("mapping_vs_leaderless_counts_CAlm.png",
       path = paste0(subDir,"/Plots/"),
       units="in",
       width = 6, height = 5)
```

Write function to access counts for each msdin and plot expr
vs. maturity

``` r
# first plot maturity for CA vs EU samples
mini_meta <- metadata[,c("File_prefix","Sample_ID","Origin","Population","Maturity")]

# t test for maturity
group1 <- "CA"
group2 <- "EU"
group1_data <- mini_meta[mini_meta$Origin == group1, "Maturity"]
group2_data <- mini_meta[mini_meta$Origin == group2, "Maturity"]
mean(group1_data)
```

    ## Warning in mean.default(group1_data): argument is not numeric or logical:
    ## returning NA

    ## [1] NA

``` r
plotrix::std.error(group1_data)
```

    ##  Maturity 
    ## 0.3845909

``` r
mean(group2_data)
```

    ## Warning in mean.default(group2_data): argument is not numeric or logical:
    ## returning NA

    ## [1] NA

``` r
plotrix::std.error(group2_data)
```

    ##  Maturity 
    ## 0.2914915

``` r
# Perform t-test
test_result <- t.test(group1_data, group2_data)
  
  # Return the result (comparison name, p-value, and test statistics)
  data.frame(
    Comparison = paste(group1, "vs", group2),
    p_value = test_result$p.value,
    statistic = test_result$statistic,
    conf_low = test_result$conf.int[1],
    conf_high = test_result$conf.int[2]
  )
```

    ##   Comparison    p_value statistic  conf_low conf_high
    ## t   CA vs EU 0.09964065 -1.694132 -1.799273 0.1641848

``` r
  # extract text for plot
  p_value = paste0("p = ",round(test_result$p.value,4))

mini_meta$Origin <- factor(mini_meta$Origin, levels = c("CA",
                                                "EU"))
  
mini_meta$Population <- factor(mini_meta$Population, levels = c("Drake",
                                                        "Pet",
                                                        "Picnic",
                                                        "Aso",
                                                        "Champ",
                                                        "Doc",
                                                        "Gron",
                                                        "Lamp",
                                                        "Quill",
                                                        "Tisza"))

genotype_colors <- c("#723D00","#fb8500","#ffb703",
                       "#023047","#0C556E","#177995","#219EBC","#5CB7CE","#96CFE0","#D1E8F2")


pdf(paste0(subDir,"/Plots/maturity.pdf"), width = 4, height = 4, family="ArialMT")

ggplot(mini_meta, aes(x=Origin,y=Maturity,color=Population)) +
  geom_point(aes(colour = Population), size = 4, position = position_jitter(w = 0.4, h = 0)) +
  scale_color_manual(values = genotype_colors) +
  # add mean line
  stat_summary(fun= mean, fun.min=mean, fun.max=mean, geom="crossbar", width=0.5, linewidth=1,color="black") +
  # add error bars
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), geom="errorbar", width=0.2, linewidth=1,color="black") +
  annotate("text", x = "CA", y = 8, 
           label = p_value, hjust = -0.3) +  # add p value
  coord_cartesian(ylim = c(0, 8)) +
  theme_classic()
dev.off()
```

    ## quartz_off_screen 
    ##                 2

## MSDIN expression vs maturity of mushroom

The maturity of each mushroom was scored prior to sampling:

<figure>
<img src="/Users/songs005/Documents/Screenshot.png"
alt="A. phalloides maturity" />
<figcaption aria-hidden="true">A. phalloides maturity</figcaption>
</figure>

``` r
plot_maturity <- function (msdin) {
  #msdin <- "MLGFLVLP"
  dat <- read.csv(paste0(subDir,"/NormCount-Graphs/Data/",msdin,".csv"))
  colnames(dat)[1] <- "File_prefix"
  mini_meta <- metadata[,c("File_prefix","Maturity")]
  mini_meta2 <- metadata[,c("Sample_ID","Maturity")]
  
  # merge
  dat2 <- merge(dat,mini_meta,by="File_prefix",all.x=TRUE)
  
  dat2$Origin <- factor(dat2$Origin, levels = c("CA",
                                                "EU"))
  
  dat2$Population <- factor(dat2$Population, levels = c("Drake",
                                                        "Pet",
                                                        "Picnic",
                                                        "Aso",
                                                        "Champ",
                                                        "Doc",
                                                        "Gron",
                                                        "Lamp",
                                                        "Quill",
                                                        "Tisza"))
  
  
  library(ggplot2)

  
  # Fit the linear model
  # for all points:
  lm_model <- lm(count ~ Maturity, data = dat2)
  
  # based on origin:
  ca <- subset(dat2,Origin == "CA")
  eu <- subset(dat2,Origin == "EU")
  
  # calculate p value with spearman ranked correlation
  corr_ca <- cor.test(x=ca$Maturity, y=ca$count, method = "kendall")
  corr_ca
  corr_ca2 <- corr_ca[["p.value"]]
  
  corr_eu <- cor.test(x=eu$Maturity, y=eu$count, method = "kendall")
  corr_eu2 <- corr_eu[["p.value"]]
  
  # do linear regression
  lm_model_ca <- lm(count ~ Maturity, data = ca)
  lm_model_eu <- lm(count ~ Maturity, data = eu)
  
  eqn_ca <- substitute(italic(y) == a + b %.% italic(x) * "," ~~ italic(R)^2 ~ "=" ~ r2, 
                       list(a = as.numeric(format(coef(lm_model_ca)[1], digits = 3)),
                            b = as.numeric(format(coef(lm_model_ca)[2], digits = 3)),
                            r2 = format(summary(lm_model_ca)$r.squared, digits = 3)))
  
  eqn_eu <- substitute(italic(y) == a + b %.% italic(x) * "," ~~ italic(R)^2 ~ "=" ~ r2, 
                       list(a = as.numeric(format(coef(lm_model_eu)[1], digits = 3)),
                            b = as.numeric(format(coef(lm_model_eu)[2], digits = 3)),
                            r2 = format(summary(lm_model_eu)$r.squared, digits = 3)))
  
  # Extract coefficients and R-squared
  eqn <- substitute(italic(y) == a + b %.% italic(x) * "," ~~ italic(R)^2 ~ "=" ~ r2, 
                    list(a = as.numeric(format(coef(lm_model)[1], digits = 3)),
                         b = as.numeric(format(coef(lm_model)[2], digits = 3)),
                         r2 = format(summary(lm_model)$r.squared, digits = 3)))
  genotype_colors <- c("#723D00","#fb8500","#ffb703",
                       "#023047","#0C556E","#177995","#219EBC","#5CB7CE","#96CFE0","#D1E8F2")
  labeleqn <- paste(eqn_ca, eqn_eu, sep = "\n")
  
  ggplot(dat2, aes(x=Maturity,y=count,color=Population, shape=Origin)) +
    geom_point(size=3) +
    scale_color_manual(values = genotype_colors) +
    #geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add linear regression line
    geom_smooth(
    data = subset(dat2, Origin == "CA"),
    method = "lm", se = FALSE, color = "#fb8500"
  ) +

  # One regression line for EU
  geom_smooth(
    data = subset(dat2, Origin == "EU"),
    method = "lm", se = FALSE, color = "#0C556E"
  ) +
      annotate("text", 
           x = min(dat2$Maturity), 
           y = max(dat2$count), 
           label = as.character(as.expression(eqn_ca)), 
           parse = TRUE, 
           hjust = 0,
         color = "#fb8500") +
  annotate("text", 
           x = min(dat2$Maturity), 
           y = max(dat2$count) * 0.9,  # slightly lower
           label = as.character(as.expression(eqn_eu)), 
           parse = TRUE, 
           hjust = 0,
         color = "#0C556E") +
    theme_classic()
  
  ggsave(paste0("maturity_vs_",msdin,"_counts.png"),
         path = paste0(subDir,"/Plots/"),
         units="in",
         width = 5, height = 4)
  
  # also get values for each plot and return that
  ca <- c(a = as.numeric(format(coef(lm_model_ca)[1], digits = 3)),
                            b = as.numeric(format(coef(lm_model_ca)[2], digits = 3)),
                            r2 = format(summary(lm_model_ca)$r.squared, digits = 3))
  eu <- c(a = as.numeric(format(coef(lm_model_eu)[1], digits = 3)),
                            b = as.numeric(format(coef(lm_model_eu)[2], digits = 3)),
                            r2 = format(summary(lm_model_eu)$r.squared, digits = 3))
  output <- as.data.frame(rbind(ca,eu))
  output$MSDIN <- msdin
  output$Origin <- row.names(output)
  output$pval <- c(corr_ca2,corr_eu2)
  return(output)
}

# now run on each msdin and save the output
# make sure to only run on list of msdins 
# MSDINs, leadered and leaderless:
list <- msdins$X4

# first see which of these genes have any counts at all
# get counts matrix for all genes
counts <- counts(dds, normalized=TRUE)
# find list of RiPPs also in count matrix
list <- intersect(list, rownames(counts))

results_df <- do.call(rbind, lapply(list, function(genename) {
  plot_maturity(msdin = genename)
}))
```

    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties

    ## Warning in cor.test.default(x = eu$Maturity, y = eu$count, method = "kendall"):
    ## Cannot compute exact p-value with ties

    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties
    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties

    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties
    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties

    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties
    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties

    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties
    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties

    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties
    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties

    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties
    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties

    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties
    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties

    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties
    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties

    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties
    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties

    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties
    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties

    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties
    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties

    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties
    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties

    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties
    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties

    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties
    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties

    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties
    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties

    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties
    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties

    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties
    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties

    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties
    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties

    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties
    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties

    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties
    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties

    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties
    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties

    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties
    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties

    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties
    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties

    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties
    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties

    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties
    ## Warning in cor.test.default(x = ca$Maturity, y = ca$count, method = "kendall"):
    ## Cannot compute exact p-value with ties

    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

``` r
rmarkdown::paged_table(results_df, options = NULL)
```

<script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["a"],"name":[1],"type":["chr"],"align":["left"]},{"label":["b"],"name":[2],"type":["chr"],"align":["left"]},{"label":["r2"],"name":[3],"type":["chr"],"align":["left"]},{"label":["MSDIN"],"name":[4],"type":["chr"],"align":["left"]},{"label":["Origin"],"name":[5],"type":["chr"],"align":["left"]},{"label":["pval"],"name":[6],"type":["dbl"],"align":["right"]}],"data":[{"1":"-53934","2":"26331","3":"0.617","4":"MLGFLVLP","5":"ca","6":"0.0004314902","_rn_":"ca"},{"1":"1905","2":"865","3":"0.00685","4":"MLGFLVLP","5":"eu","6":"0.5576087627","_rn_":"eu"},{"1":"184","2":"-22.9","3":"0.171","4":"AWLATCP_1","5":"ca","6":"0.3370504430","_rn_":"ca1"},{"1":"261","2":"-28.9","3":"0.124","4":"AWLATCP_1","5":"eu","6":"0.2911900393","_rn_":"eu1"},{"1":"30.7","2":"-1.41","3":"0.0524","4":"LIQRPFAP_1","5":"ca","6":"0.3370504430","_rn_":"ca2"},{"1":"32.3","2":"-2.71","3":"0.149","4":"LIQRPFAP_1","5":"eu","6":"0.0461794057","_rn_":"eu2"},{"1":"1.7","2":"0.852","3":"0.0663","4":"LGRPESLP","5":"ca","6":"0.1547241742","_rn_":"ca3"},{"1":"31","2":"-1.86","3":"0.0377","4":"LGRPESLP","5":"eu","6":"0.4116714486","_rn_":"eu3"},{"1":"-1.17","2":"2.46","3":"0.109","4":"LRLPPFMIPP","5":"ca","6":"0.0705364206","_rn_":"ca4"},{"1":"4.41","2":"5.02","3":"0.0493","4":"LRLPPFMIPP","5":"eu","6":"0.4565528101","_rn_":"eu4"},{"1":"1.81","2":"-0.0561","3":"0.00322","4":"IIGILLPP","5":"ca","6":"0.6204221886","_rn_":"ca5"},{"1":"32.3","2":"-3.15","3":"0.0185","4":"IIGILLPP","5":"eu","6":"0.7376579728","_rn_":"eu5"},{"1":"14.7","2":"-0.513","3":"0.0163","4":"AWLVDCP","5":"ca","6":"0.6957122528","_rn_":"ca6"},{"1":"25.3","2":"0.811","3":"0.00119","4":"AWLVDCP","5":"eu","6":"0.9688161358","_rn_":"eu6"},{"1":"74.1","2":"-10.1","3":"0.517","4":"IWGIGCNP","5":"ca","6":"0.0039763929","_rn_":"ca7"},{"1":"716","2":"-114","3":"0.325","4":"IWGIGCNP","5":"eu","6":"0.1273508641","_rn_":"eu7"},{"1":"20.5","2":"-2.35","3":"0.269","4":"IWGIGCDP_1","5":"ca","6":"0.1262882472","_rn_":"ca8"},{"1":"49.6","2":"-7.01","3":"0.346","4":"IWGIGCDP_1","5":"eu","6":"0.0088125945","_rn_":"eu8"},{"1":"108","2":"-12.2","3":"0.228","4":"IWGIGCDP_2","5":"ca","6":"0.1262882472","_rn_":"ca9"},{"1":"93.4","2":"-13.1","3":"0.136","4":"IWGIGCDP_2","5":"eu","6":"0.5063170780","_rn_":"eu9"},{"1":"7.92","2":"0.329","3":"0.0265","4":"IFLAFPIPP","5":"ca","6":"0.5455442266","_rn_":"ca10"},{"1":"10","2":"-0.261","3":"0.00975","4":"IFLAFPIPP","5":"eu","6":"0.7843521213","_rn_":"eu10"},{"1":"8.38","2":"-0.264","3":"0.00395","4":"IFWFIYFP","5":"ca","6":"0.8034448288","_rn_":"ca11"},{"1":"27.2","2":"-2.89","3":"0.234","4":"IFWFIYFP","5":"eu","6":"0.0785450746","_rn_":"eu11"},{"1":"12.7","2":"26","3":"0.337","4":"LIQRPFAP_2","5":"ca","6":"0.0039763929","_rn_":"ca12"},{"1":"-31.7","2":"21.6","3":"0.177","4":"LIQRPFAP_2","5":"eu","6":"0.0210829976","_rn_":"eu12"},{"1":"934","2":"-88.5","3":"0.101","4":"FNLFRFPYP","5":"ca","6":"0.0814665572","_rn_":"ca13"},{"1":"1772","2":"-64.2","3":"0.00826","4":"FNLFRFPYP","5":"eu","6":"0.6113055670","_rn_":"eu13"},{"1":"-0.273","2":"1.17","3":"0.305","4":"VQKPWSRP","5":"ca","6":"0.0055213489","_rn_":"ca14"},{"1":"2.72","2":"0.499","3":"0.0238","4":"VQKPWSRP","5":"eu","6":"0.8450284144","_rn_":"eu14"},{"1":"4.31","2":"0.77","3":"0.0377","4":"IRLPPLFLPP","5":"ca","6":"0.5686096359","_rn_":"ca15"},{"1":"24.9","2":"-2.26","3":"0.0931","4":"IRLPPLFLPP","5":"eu","6":"0.2557909986","_rn_":"eu15"},{"1":"6.81","2":"-0.399","3":"0.0555","4":"FFFPPFFIPP","5":"ca","6":"0.1262882472","_rn_":"ca16"},{"1":"4.26","2":"0.554","3":"0.00574","4":"FFFPPFFIPP","5":"eu","6":"0.8124777480","_rn_":"eu16"},{"1":"23.2","2":"3.5","3":"0.03","4":"FNILPFMLPP","5":"ca","6":"0.5937992115","_rn_":"ca17"},{"1":"138","2":"-13.1","3":"0.0259","4":"FNILPFMLPP","5":"eu","6":"0.5576087627","_rn_":"eu17"},{"1":"8.27","2":"0.3","3":"0.00993","4":"FFPIVFSPP","5":"ca","6":"0.5455442266","_rn_":"ca18"},{"1":"1.22","2":"1.5","3":"0.0746","4":"FFPIVFSPP","5":"eu","6":"0.9066386495","_rn_":"eu18"},{"1":"65.2","2":"-9.48","3":"0.0871","4":"HFASFIPP","5":"ca","6":"0.8588961211","_rn_":"ca19"},{"1":"35.1","2":"-3.77","3":"0.202","4":"HFASFIPP","5":"eu","6":"0.0461794057","_rn_":"eu19"},{"1":"3.51","2":"2.05","3":"0.0519","4":"ISDPTAYP","5":"ca","6":"0.4993186301","_rn_":"ca20"},{"1":"1.06","2":"6.35","3":"0.0424","4":"ISDPTAYP","5":"eu","6":"0.2255555914","_rn_":"eu20"},{"1":"26.5","2":"2.18","3":"0.0148","4":"LFFWFWFLWP","5":"ca","6":"0.4552603429","_rn_":"ca21"},{"1":"7.06","2":"0.855","3":"0.0103","4":"LFFWFWFLWP","5":"eu","6":"0.7843521213","_rn_":"eu21"},{"1":"44117","2":"4142","3":"0.0453","4":"SFFFPIP","5":"ca","6":"0.3740580788","_rn_":"ca22"},{"1":"-1863","2":"819","3":"0.0284","4":"SFFFPIP","5":"eu","6":"0.9066386495","_rn_":"eu22"},{"1":"3.16","2":"0.481","3":"0.0783","4":"FMPLAP","5":"ca","6":"0.3740580788","_rn_":"ca23"},{"1":"7.13","2":"-0.0248","3":"0.0000706","4":"FMPLAP","5":"eu","6":"0.9066386495","_rn_":"eu23"},{"1":"530","2":"13.1","3":"0.00997","4":"LILLAALGIP","5":"ca","6":"0.4552603429","_rn_":"ca24"},{"1":"183","2":"-15.8","3":"0.0871","4":"LILLAALGIP","5":"eu","6":"0.1712308140","_rn_":"eu24"},{"1":"10470","2":"10227","3":"0.223","4":"GPVFFA","5":"ca","6":"0.0359229494","_rn_":"ca25"},{"1":"-2138","2":"978","3":"0.0265","4":"GPVFFA","5":"eu","6":"0.0382714053","_rn_":"eu25"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>

``` r
write.csv(results_df,paste0(subDir,"/maturity_msdins.csv"))

###########
# also make nice leaderless plot
msdin <- "MLGFLVLP"
  dat <- read.csv(paste0(subDir,"/NormCount-Graphs/Data/",msdin,".csv"))
  colnames(dat)[1] <- "File_prefix"
  mini_meta <- metadata[,c("File_prefix","Maturity")]
  mini_meta2 <- metadata[,c("Sample_ID","Maturity")]
  
  # merge
  dat2 <- merge(dat,mini_meta,by="File_prefix",all.x=TRUE)
  
  dat2$Origin <- factor(dat2$Origin, levels = c("CA",
                                                "EU"))
  
  dat2$Population <- factor(dat2$Population, levels = c("Drake",
                                                        "Pet",
                                                        "Picnic",
                                                        "Aso",
                                                        "Champ",
                                                        "Doc",
                                                        "Gron",
                                                        "Lamp",
                                                        "Quill",
                                                        "Tisza"))
  
  
  library(ggplot2)
  
  # Fit the linear model
  # for all points:
  lm_model <- lm(count ~ Maturity, data = dat2)
  
  # based on origin:
  ca <- subset(dat2,Origin == "CA")
  eu <- subset(dat2,Origin == "EU")
  
  lm_model_ca <- lm(count ~ Maturity, data = ca)
  lm_model_eu <- lm(count ~ Maturity, data = eu)
  
  eqn_ca <- substitute(italic(y) == a + b %.% italic(x) * "," ~~ italic(R)^2 ~ "=" ~ r2, 
                       list(a = as.numeric(format(coef(lm_model_ca)[1], digits = 3)),
                            b = as.numeric(format(coef(lm_model_ca)[2], digits = 3)),
                            r2 = format(summary(lm_model_ca)$r.squared, digits = 3)))
  
  eqn_eu <- substitute(italic(y) == a + b %.% italic(x) * "," ~~ italic(R)^2 ~ "=" ~ r2, 
                       list(a = as.numeric(format(coef(lm_model_eu)[1], digits = 3)),
                            b = as.numeric(format(coef(lm_model_eu)[2], digits = 3)),
                            r2 = format(summary(lm_model_eu)$r.squared, digits = 3)))
  
  # Extract coefficients and R-squared
  eqn <- substitute(italic(y) == a + b %.% italic(x) * "," ~~ italic(R)^2 ~ "=" ~ r2, 
                    list(a = as.numeric(format(coef(lm_model)[1], digits = 3)),
                         b = as.numeric(format(coef(lm_model)[2], digits = 3)),
                         r2 = format(summary(lm_model)$r.squared, digits = 3)))
  genotype_colors <- c("#723D00","#fb8500","#ffb703",
                       "#023047","#0C556E","#177995","#219EBC","#5CB7CE","#96CFE0","#D1E8F2")
  labeleqn <- paste(eqn_ca, eqn_eu, sep = "\n")
  
  
pdf(paste0(subDir,"/Plots/maturity_leaderless.pdf"), width = 5, height = 4, family="ArialMT")
  
  ggplot(dat2, aes(x=Maturity,y=count,color=Population, shape=Origin)) +
    geom_point(size=3) +
    scale_color_manual(values = genotype_colors) +
    #geom_smooth(method = "lm", se = FALSE, color = "black") +  # Add linear regression line
    geom_smooth(
    data = subset(dat2, Origin == "CA"),
    method = "lm", se = FALSE, color = "#fb8500"
  ) +

  # One regression line for EU
  geom_smooth(
    data = subset(dat2, Origin == "EU"),
    method = "lm", se = FALSE, color = "#0C556E"
  ) +
      annotate("text", 
           x = min(dat2$Maturity), 
           y = max(dat2$count), 
           label = as.character(as.expression(eqn_ca)), 
           parse = TRUE, 
           hjust = 0,
         color = "#fb8500") +
  annotate("text", 
           x = min(dat2$Maturity), 
           y = max(dat2$count) * 0.9,  # slightly lower
           label = as.character(as.expression(eqn_eu)), 
           parse = TRUE, 
           hjust = 0,
         color = "#0C556E") +
    theme_classic()
```

    ## `geom_smooth()` using formula = 'y ~ x'
    ## `geom_smooth()` using formula = 'y ~ x'

``` r
  dev.off()
```

    ## quartz_off_screen 
    ##                 2

## PCA of MSDINs only

``` r
# subset dds to only include msdins

# generate dds object
library(DESeq2)
setClassUnion("ExpData", c("matrix", "SummarizedExperiment"))
dds2 <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Origin)
```

    ## converting counts to integer mode

``` r
# keep only the msdins
list <- msdins$X4

# first see which of these genes have any counts at all
# Check which entries in 'list' are NOT in the rownames
missing <- setdiff(list, rownames(dds2))
missing
```

    ## character(0)

``` r
valid_list <- intersect(list, rownames(dds2))
# this should be all msdins
valid_list
```

    ##  [1] "MYNPPYFLPP" "MLGFLVLP"   "MLPGMVAFS"  "AWLATCP_1"  "AWLATCP_2" 
    ##  [6] "LIQRPFAP_1" "LGRPESLP"   "LRLPPFMIPP" "IIGILLPP"   "AWLVDCP"   
    ## [11] "IWGIGCNP"   "IWGIGCDP_1" "IWGIGCDP_2" "TIYYLYFIP"  "IFLAFPIPP" 
    ## [16] "IFWFIYFP"   "LIQRPFAP_2" "FNLFRFPYP"  "VQKPWSRP"   "IRLPPLFLPP"
    ## [21] "FFFPPFFIPP" "FNILPFMLPP" "FFPIVFSPP"  "HFASFIPP"   "ISDPTAYP"  
    ## [26] "LFFWFWFLWP" "SFFFPIP"    "FMPLAP"     "LILLAALGIP" "GPVFFA"    
    ## [31] "GVILIIP"    "FIFPPFFIPP" "LPILPIPPLP"

``` r
dds2 <- dds2[valid_list, ]

dds2 <- DESeq(dds2)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## -- note: fitType='parametric', but the dispersion trend was not well captured by the
    ##    function: y = a/x + b, and a local regression fit was automatically substituted.
    ##    specify fitType='local' or 'mean' to avoid this message next time.

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 2 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

``` r
vsd2 <- varianceStabilizingTransformation(dds2, blind=FALSE)

# find distances between samples
sampleDists <- dist(t(assay(vsd2)))

# reformat distances
sampleDistMatrix <- as.matrix(sampleDists)

# make new data frame with only Population stage info
#metadata <- readxl::read_xlsx(paste0(dir,"/metadata.xlsx"),sheet="metadata")

#metadata <- metadata[!metadata$Sample_ID %in% c("Pet1-redo","Pet5"),]

samples <- metadata$Sample_ID

# make row and column names the same as sample IDs
rownames(sampleDistMatrix) <- samples
colnames(sampleDistMatrix) <- samples

# generate some nice colors
library(RColorBrewer)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

# plot it
library(pheatmap)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
```

![](Script3_DESeq2_files/figure-markdown_github/pca_msdins_only-1.png)

``` r
# next lets generate a PCA plot on variance stabilized data
vsd_out <- assay(vsd2)

# run PCA analysis
pca <- prcomp(t(vsd_out))

# plot scree
library(factoextra)
fviz_eig(pca)
```

![](Script3_DESeq2_files/figure-markdown_github/pca_msdins_only-2.png)

``` r
# Create data frame with metadata and PC3 and PC4 values for input to ggplot
pcaData <- cbind(metadata[,c(2:5,8,13,15)], pca$x)

# Define color palette for genotypes
pcaData$Origin <- factor(pcaData$Origin, levels = c("CA",
                                                    "EU"))
# find percent variance for PC1 and PC2
# summary of variance
pcasummary <- summary(pca)
pcasummary
```

    ## Importance of components:
    ##                            PC1     PC2     PC3    PC4     PC5     PC6     PC7
    ## Standard deviation     14.8640 5.99708 5.67979 4.3618 3.98125 3.71512 3.27998
    ## Proportion of Variance  0.5435 0.08848 0.07936 0.0468 0.03899 0.03395 0.02647
    ## Cumulative Proportion   0.5435 0.63200 0.71137 0.7582 0.79716 0.83112 0.85758
    ##                            PC8     PC9    PC10    PC11    PC12    PC13    PC14
    ## Standard deviation     2.98778 2.82122 2.62042 2.24451 2.18749 1.92729 1.83987
    ## Proportion of Variance 0.02196 0.01958 0.01689 0.01239 0.01177 0.00914 0.00833
    ## Cumulative Proportion  0.87954 0.89912 0.91602 0.92841 0.94018 0.94932 0.95765
    ##                           PC15   PC16    PC17    PC18    PC19    PC20    PC21
    ## Standard deviation     1.78200 1.5997 1.42650 1.31553 1.19035 1.11217 1.08699
    ## Proportion of Variance 0.00781 0.0063 0.00501 0.00426 0.00349 0.00304 0.00291
    ## Cumulative Proportion  0.96546 0.9718 0.97676 0.98102 0.98450 0.98755 0.99045
    ##                           PC22   PC23    PC24    PC25    PC26   PC27    PC28
    ## Standard deviation     1.03047 0.8542 0.79616 0.65940 0.61565 0.4919 0.41226
    ## Proportion of Variance 0.00261 0.0018 0.00156 0.00107 0.00093 0.0006 0.00042
    ## Cumulative Proportion  0.99307 0.9949 0.99642 0.99749 0.99842 0.9990 0.99944
    ##                           PC29    PC30   PC31    PC32                 PC33
    ## Standard deviation     0.33936 0.26002 0.2005 0.07742 0.000000000000001112
    ## Proportion of Variance 0.00028 0.00017 0.0001 0.00001 0.000000000000000000
    ## Cumulative Proportion  0.99972 0.99989 1.0000 1.00000 1.000000000000000000

``` r
pc1var <- round(pcasummary$importance[2],digits = 4) * 100
pc2var <- round(pcasummary$importance[5],digits = 4) * 100
pc3var <- round(pcasummary$importance[8],digits = 4) * 100
pc4var <- round(pcasummary$importance[11],digits = 4) * 100


pcaData$Population <- factor(pcaData$Population, levels = c("Drake",
                                                            "Pet",
                                                            "Picnic",
                                                            "Aso",
                                                            "Champ",
                                                            "Doc",
                                                            "Gron",
                                                            "Lamp",
                                                            "Quill",
                                                            "Tisza"))



pdf(paste0(subDir,"/Plots/pca_msdins_only.pdf"), width = 4.5, height = 2.5, family="ArialMT")  # Specify desired width and height

ggplot(pcaData,aes(x=PC1, y=PC2, color = Population, shape=Origin)) +
  geom_point(size=3) + 
  theme_classic() +
  scale_color_manual(values = genotype_colors) +
stat_ellipse(aes(group = Origin), colour = "black")+                 
  xlab(paste0("PC1 (",pc1var,"% variance)")) +
  ylab(paste0("PC2 (",pc2var,"% variance)")) +
  theme(panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        text = element_text(color = "black",size = 10),
        legend.key.height = unit(0.15, 'in'))+ 
  guides(color = guide_legend(override.aes=list(shape = c(rep(16,3),rep(17,7)))))  # Customize shape of color legend


dev.off()
```

    ## quartz_off_screen 
    ##                 2

## PCA of all non-MSDIN transcripts

``` r
# subset dds to exclude msdins
# generate dds object
library(DESeq2)
dds2 <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ Origin)
```

    ## converting counts to integer mode

``` r
# exclude the msdins
list <- msdins$X4

valid_list <- !rownames(dds2) %in% list
dds2 <- dds2[valid_list, ]

dds2 <- DESeq(dds2)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 112 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

``` r
vsd2 <- varianceStabilizingTransformation(dds2, blind=FALSE)

# find distances between samples
sampleDists <- dist(t(assay(vsd2)))

# reformat distances
sampleDistMatrix <- as.matrix(sampleDists)

# make new data frame with only Population stage info
df <- as.data.frame(colData(dds2)[,c("Population")])
colnames(df) <- c("Population")
#samples <- colnames(assay(ntd))
samples <- metadata$Sample_ID

# make row and column names the same as sample IDs
rownames(sampleDistMatrix) <- samples
colnames(sampleDistMatrix) <- samples

# generate some nice colors
library(RColorBrewer)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

# plot it
library(pheatmap)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
```

![](Script3_DESeq2_files/figure-markdown_github/pca_no_msdins-1.png)

``` r
# next lets generate a PCA plot on variance stabilized data
vsd_out <- assay(vsd2)

# run PCA analysis
pca <- prcomp(t(vsd_out))

# plot scree
library(factoextra)
fviz_eig(pca)
```

![](Script3_DESeq2_files/figure-markdown_github/pca_no_msdins-2.png)

``` r
# Create data frame with metadata and PC3 and PC4 values for input to ggplot
pcaData <- cbind(metadata[,c(2:5,8,13,15)], pca$x)

# Define color palette for genotypes
pcaData$Origin <- factor(pcaData$Origin, levels = c("CA",
                                                    "EU"))
# find percent variance for PC1 and PC2
# summary of variance
pcasummary <- summary(pca)
pcasummary
```

    ## Importance of components:
    ##                            PC1     PC2      PC3      PC4      PC5      PC6
    ## Standard deviation     47.3575 24.9699 18.45869 14.18347 12.92323 11.21752
    ## Proportion of Variance  0.4691  0.1304  0.07127  0.04208  0.03493  0.02632
    ## Cumulative Proportion   0.4691  0.5995  0.67081  0.71289  0.74783  0.77415
    ##                            PC7     PC8     PC9    PC10   PC11    PC12    PC13
    ## Standard deviation     9.88720 9.54770 8.95558 8.70501 7.8231 7.63303 7.24637
    ## Proportion of Variance 0.02045 0.01907 0.01678 0.01585 0.0128 0.01219 0.01098
    ## Cumulative Proportion  0.79459 0.81366 0.83044 0.84629 0.8591 0.87128 0.88226
    ##                           PC14    PC15    PC16    PC17   PC18   PC19    PC20
    ## Standard deviation     6.90728 6.67042 6.49784 6.08704 6.0258 5.7869 5.68652
    ## Proportion of Variance 0.00998 0.00931 0.00883 0.00775 0.0076 0.0070 0.00676
    ## Cumulative Proportion  0.89224 0.90155 0.91038 0.91813 0.9257 0.9327 0.93950
    ##                           PC21    PC22    PC23    PC24    PC25    PC26    PC27
    ## Standard deviation     5.59398 5.34926 5.28243 4.81172 4.72775 4.53029 4.36940
    ## Proportion of Variance 0.00655 0.00599 0.00584 0.00484 0.00468 0.00429 0.00399
    ## Cumulative Proportion  0.94604 0.95203 0.95786 0.96271 0.96738 0.97167 0.97567
    ##                           PC28    PC29    PC30    PC31    PC32    PC33    PC34
    ## Standard deviation     4.27551 4.15915 3.99409 3.80822 3.76587 3.41928 3.17808
    ## Proportion of Variance 0.00382 0.00362 0.00334 0.00303 0.00297 0.00245 0.00211
    ## Cumulative Proportion  0.97949 0.98311 0.98645 0.98948 0.99245 0.99489 0.99701
    ##                           PC35    PC36                PC37
    ## Standard deviation     2.77734 2.56994 0.00000000000004139
    ## Proportion of Variance 0.00161 0.00138 0.00000000000000000
    ## Cumulative Proportion  0.99862 1.00000 1.00000000000000000

``` r
pc1var <- round(pcasummary$importance[2],digits = 4) * 100
pc2var <- round(pcasummary$importance[5],digits = 4) * 100
pc3var <- round(pcasummary$importance[8],digits = 4) * 100
pc4var <- round(pcasummary$importance[11],digits = 4) * 100


pcaData$Population <- factor(pcaData$Population, levels = c("Drake",
                                                            "Pet",
                                                            "Picnic",
                                                            "Aso",
                                                            "Champ",
                                                            "Doc",
                                                            "Gron",
                                                            "Lamp",
                                                            "Quill",
                                                            "Tisza"))



pdf(paste0(subDir,"/Plots/pca_no_msdins.pdf"), width = 4.5, height = 2.5, family="ArialMT")  # Specify desired width and height

ggplot(pcaData,aes(x=PC1, y=PC2, color = Population, shape=Origin)) +
  geom_point(size=3) + 
  theme_classic() +
  scale_color_manual(values = genotype_colors) +
stat_ellipse(aes(group = Origin), colour = "black")+                 
  xlab(paste0("PC1 (",pc1var,"% variance)")) +
  ylab(paste0("PC2 (",pc2var,"% variance)")) +
  theme(panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        text = element_text(color = "black",size = 10),
        legend.key.height = unit(0.15, 'in'))+ 
  guides(color = guide_legend(override.aes=list(shape = c(rep(16,3),rep(17,7)))))  # Customize shape of color legend


dev.off()
```

    ## quartz_off_screen 
    ##                 2

# Session info

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
    ## [1] grid      stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] multcomp_1.4-28             TH.data_1.1-3              
    ##  [3] MASS_7.3-65                 survival_3.8-3             
    ##  [5] mvtnorm_1.3-3               lubridate_1.9.4            
    ##  [7] forcats_1.0.0               stringr_1.5.1              
    ##  [9] purrr_1.0.4                 tibble_3.3.0               
    ## [11] tidyverse_2.0.0             ggbreak_0.1.4              
    ## [13] factoextra_1.0.7            RColorBrewer_1.1-3         
    ## [15] Rsamtools_2.24.0            GenomicFeatures_1.60.0     
    ## [17] AnnotationDbi_1.70.0        Biostrings_2.76.0          
    ## [19] XVector_0.48.0              rtracklayer_1.68.0         
    ## [21] pheatmap_1.0.13             vsn_3.76.0                 
    ## [23] EnhancedVolcano_1.26.0      ggrepel_0.9.6              
    ## [25] ggplot2_3.5.2               readr_2.1.5                
    ## [27] tidyr_1.3.1                 dplyr_1.1.4                
    ## [29] rmarkdown_2.29              DESeq2_1.48.1              
    ## [31] SummarizedExperiment_1.38.1 Biobase_2.68.0             
    ## [33] MatrixGenerics_1.20.0       matrixStats_1.5.0          
    ## [35] GenomicRanges_1.60.0        GenomeInfoDb_1.44.0        
    ## [37] IRanges_2.42.0              S4Vectors_0.46.0           
    ## [39] BiocGenerics_0.54.0         generics_0.1.4             
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] rstudioapi_0.17.1        jsonlite_2.0.0           magrittr_2.0.3          
    ##   [4] farver_2.1.2             fs_1.6.6                 BiocIO_1.18.0           
    ##   [7] ragg_1.4.0               vctrs_0.6.5              memoise_2.0.1           
    ##  [10] RCurl_1.98-1.17          rstatix_0.7.2            htmltools_0.5.8.1       
    ##  [13] S4Arrays_1.8.1           progress_1.2.3           plotrix_3.8-4           
    ##  [16] curl_6.4.0               broom_1.0.8              cellranger_1.1.0        
    ##  [19] gridGraphics_0.5-1       Formula_1.2-5            SparseArray_1.8.0       
    ##  [22] sandwich_3.1-1           httr2_1.1.2              zoo_1.8-14              
    ##  [25] cachem_1.1.0             GenomicAlignments_1.44.0 lifecycle_1.0.4         
    ##  [28] pkgconfig_2.0.3          Matrix_1.7-3             R6_2.6.1                
    ##  [31] fastmap_1.2.0            GenomeInfoDbData_1.2.14  aplot_0.2.7             
    ##  [34] digest_0.6.37            colorspace_2.1-1         patchwork_1.3.1         
    ##  [37] textshaping_1.0.1        RSQLite_2.4.2            ggpubr_0.6.0            
    ##  [40] filelock_1.0.3           labeling_0.4.3           timechange_0.3.0        
    ##  [43] mgcv_1.9-1               httr_1.4.7               abind_1.4-8             
    ##  [46] compiler_4.5.0           bit64_4.6.0-1            withr_3.0.2             
    ##  [49] backports_1.5.0          BiocParallel_1.42.1      carData_3.0-5           
    ##  [52] DBI_1.2.3                hexbin_1.28.5            ggsignif_0.6.4          
    ##  [55] biomaRt_2.64.0           rappdirs_0.3.3           DelayedArray_0.34.1     
    ##  [58] rjson_0.2.23             tools_4.5.0              glue_1.8.0              
    ##  [61] restfulr_0.0.15          nlme_3.1-168             gtable_0.3.6            
    ##  [64] tzdb_0.5.0               preprocessCore_1.70.0    hms_1.1.3               
    ##  [67] car_3.1-3                xml2_1.3.8               utf8_1.2.6              
    ##  [70] pillar_1.10.2            yulab.utils_0.2.0        vroom_1.6.5             
    ##  [73] limma_3.64.1             splines_4.5.0            BiocFileCache_2.16.1    
    ##  [76] lattice_0.22-6           bit_4.6.0                tidyselect_1.2.1        
    ##  [79] locfit_1.5-9.12          knitr_1.50               xfun_0.52               
    ##  [82] statmod_1.5.0            stringi_1.8.7            UCSC.utils_1.4.0        
    ##  [85] ggfun_0.1.9              yaml_2.3.10              evaluate_1.0.4          
    ##  [88] codetools_0.2-20         BiocManager_1.30.26      ggplotify_0.1.2         
    ##  [91] cli_3.6.5                affyio_1.78.0            systemfonts_1.2.3       
    ##  [94] dichromat_2.0-0.1        Rcpp_1.0.14              readxl_1.4.5            
    ##  [97] dbplyr_2.5.0             png_0.1-8                XML_3.99-0.18           
    ## [100] parallel_4.5.0           blob_1.2.4               prettyunits_1.2.0       
    ## [103] bitops_1.0-9             txdbmaker_1.4.2          viridisLite_0.4.2       
    ## [106] scales_1.4.0             affy_1.86.0              crayon_1.5.3            
    ## [109] rlang_1.1.6              KEGGREST_1.48.1
