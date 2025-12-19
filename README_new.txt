README.txt

Please read to the end before attempting to use scripts

Note that the "ReadME" variable in the MSDIN_blast_new.sh script has to be changed from "NO" to "YES" for this script to work.

The following scripts were developed by Milton T. Drott for the manuscript entitled ?angenomics of the death cap mushroom Amanita phalloides, and of Agaricales, reveals dynamic evolution of toxin genes in an invasive range?please cite accordingly if using this pipeline.
Scripts have been updated to accommodate leaderless MSDINs by Milton T. Drott for the manuscript entitled "Leaderless RiPPs expand the chemical repertoire of the invasive death cap mushroom" please cite accordingly if using the leaderless portion of the pipeline 

The example input contains all of the scripts you will need to find MSDIN genes. Additionally, a genome (A. Phalloides 10511) is included to facilitate testing of scripts before using it on your own data. Scripts are included as .sh files. Each script starts with hashed-out instructions of what the user needs to have done before using the script (e.g., installing certain programs like BLAST and setting user-defined variables). Default values were used in manuscripts mentioned above, but may need to be changed for future efforts. If content with using past parameters, see ?uick start?section below for a quick?start?

In the example_input scripts have already been run and their output can be found in the corresponding folders (as defined at the header of the scripts). You can delete these folders if you? like a fresh start using the following command:

rm -r putative_msdin_prots/ preprocessed_msdins/ processed_msdin/ mast_output/ deduped_and_gene/ final_msdin/

############ 
Some necessary files are included from the associated manuscript but could be changed to reflect future discoveries:
############ 
MSDIN_genes_all_lib.txt    --    a tsp library of all known MSDIN sequences, 
final_motif.txt            --    and an associated MEME motif files 
name_list                  --    a list of all accessions that you will be processing -- naming must match genome name exactly


############
Quick start:
############ 
-Install bedtools, EMBOSS, meme suite, samtools, coreutils, and BLAST and add them to PATH. 
-Confirm that the above are in PATH.
-change Path_to_mast variable in MSDIN_motif.sh to be blank between quotes
-change ReadME variable in MSDIN_blast.sh to read "YES"
-run scripts:

Please note that for simplicity, the example inputs have been separated into leadered and leaderless subfolders. Each subfolder already contains the relevant scripts. For consistency, I recommend running these sets of scripts in separate folders are exemplified. For quick start, navigate to the respective folders and run:
scripts may throw "rm, no such file or directory" errors if you remove the existing folders, they should still work fine -- no need for alarm.

For leaderless:

bash MSDIN_blast_new.sh
bash MSDIN_motif_runmodes.sh
bash MSDIN_motif_process_runmodes_leaderless.sh
bash make_msdin_bed_and_parse_motif_leaderless.sh
bash MSDIN_score_and_locus_leaderless.sh


For leadered:

bash MSDIN_blast_new.sh
bash MSDIN_motif.sh                           #if running both leaders and leaderless analysis, reference MSDIN_blast_new.sh output for both analyses and avoid having to run this twice
bash MSDIN_motif_process_runmodes.sh
bash make_msdin_bed_and_parse_motif_large.sh
and
bash MSDIN_score_and_locus.sh
         -- or --
bash MSDIN_score_and_locus_corepro.sh.        # not used in Park et al. to keep scoring consistent with Drott et al. 2023


-See header in MSDIN_score_and_locus.sh file for information about final outputs
-If using your own input data please see information in MDSIN_blast_new.sh for input formatting. You will need to make a new name list in the folder you?e running scripts.
-relevant outputs are put into a folder called final_msdin. Files are defined in the header of MSDIN_score_and_locus.sh 

############ 
Troubleshooting:
############ 
 if the scripts appear to not be working please verify that all of the required programs have been installed and are exported to PATH. Also confirm that variables at the top of scripts are defined correctly. If problems persist compare your input/output files?format to those that are present in the example_input folder when it is first uncompressed (note that the example input also contains example outputs for all steps ?as indicated above).

Previous validation of original MSDIN finding pipeline:
Scripts were tested on OSX, Catalina using a range of user-defined parameters with the following versions:
EMBOSS v6.6.0 
Meme 5.3.3 
Bedtools 2.30.0 
Blast 2.12.0 

Scripts were tested on Linux, CentOS using only the default user-defined parameters with the following versions:
EMBOSS v6.5.7 
Meme 5.5.1 
Bedtools 2.30.0 
Blast 2.13.0 

Validation of new non canonical MSDIN finding pipeline, high-throughput msdin bed parser, and scoring with core pro -- only on OSX Sonoma:
EMBOSS v6.6.0 
Meme 5.5.6
Bedtools 2.31.1
Blast 2.16.0 



In trying to find EMBOSS again during testing I had some trouble. The binaries seem to still be available as of this writing from: ftp://emboss.open-bio.org/pub/EMBOSS/old/6.5.0/EMBOSS-6.5.7.tar.gz

