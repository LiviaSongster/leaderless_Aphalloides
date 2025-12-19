# Leaderless MSDINs in A. phalloides.
## Github repository contents:
- [Software S1](https://github.com/LiviaSongster/leaderless_Aphalloides/tree/main/Software%20S1): Scripts for MSDIN finding from Park et al. 2025 (leaderless MSDINs) and Drott et al. 2023 (leadered MSDINs); Author: Milton T. Drott 
- [rnaseq_scripts](https://github.com/LiviaSongster/leaderless_Aphalloides/tree/main/rnaseq_scripts): Code related to RNA Seq data analysis and figure generation; Author: Livia D. S. Oster
- [all_msdin_cds.bed](https://github.com/LiviaSongster/leaderless_Aphalloides/blob/main/all_msdin_cds.bed)): Genome annotation file of canonical and noncanonical MSDINs in A. phalloides reference genome 10511
- [aphal_msdin_annotations.gff3](https://github.com/LiviaSongster/leaderless_Aphalloides/blob/main/aphal_msdin_annotations.gff3): Genome annotation file of canonical and noncanonical MSDINs in A. phalloides reference genome 10511
- [fasta_files.zip](https://github.com/LiviaSongster/leaderless_Aphalloides/blob/main/fasta_files.zip): fasta files used for MSDIN phylogenies in Park et al. 2025

## Detailed descriptions for [Software S1](https://github.com/LiviaSongster/leaderless_Aphalloides/tree/main/Software%20S1)
Author: Milton T. Drott

This directory contains MSDIN finding scripts from Park et al. 2025 (leaderless MSDINs) and Drott et al. 2023 (leadered MSDINs).
- Please refer to the nested [README_new.txt](https://github.com/LiviaSongster/leaderless_Aphalloides/blob/main/Software%20S1/README_new.txt) for detailed information on how to use these scripts and example usage.
- [leadered.zip](https://github.com/LiviaSongster/leaderless_Aphalloides/blob/main/Software%20S1/leadered.zip): These scripts were developed by Milton T. Drott for the manuscript entitled <i>[Pangenomics of the death cap mushroom Amanita phalloides, and of Agaricales, reveals dynamic evolution of toxin genes in an invasive range](https://doi.org/10.1038/s41396-023-01432-x)</i>. Please cite accordingly if using this pipeline.
- [leaderless.zip](https://github.com/LiviaSongster/leaderless_Aphalloides/blob/main/Software%20S1/leaderless.zip): Scripts have been updated to accommodate leaderless MSDINs by Milton T. Drott for the manuscript entitled <i>[Leaderless RiPPs expand the repertoire of fungal secondary metabolites](https://www.pnas.org/)</i>. Please cite accordingly if using the leaderless portion of the pipeline 

## Detailed descriptions for [rnaseq_scripts](https://github.com/LiviaSongster/leaderless_Aphalloides/tree/main/rnaseq_scripts)
Author: Livia D. S. Oster

Hyperlinks will connect to .md markdown files that are rendered in Github. For proper rendering & easier reviewing on your own computer, I recommend you download the [html](https://github.com/LiviaSongster/leaderless_Aphalloides/tree/main/rnaseq_scripts/html) files locally and open them in your .html viewer/browser.
- Script 1 [Add MSDINs to .gff3 gene annotation file](https://github.com/LiviaSongster/leaderless_Aphalloides/blob/main/rnaseq_scripts/Script1_gff.md)
- Script 2 [Read mapping with Hisat2 and Featurecounts](https://github.com/LiviaSongster/leaderless_Aphalloides/blob/main/rnaseq_scripts/Script2_RNAseq_genome_mapping.md)
- Script 3 [Differential expression analysis with DESeq2](https://github.com/LiviaSongster/leaderless_Aphalloides/blob/main/rnaseq_scripts/Script3_DESeq2.md)
- Script 4 [BAM alignment of reads to a leaderless MSDIN gene](https://github.com/LiviaSongster/leaderless_Aphalloides/blob/main/rnaseq_scripts/Script4_leaderless_bam_alignment.md)
