# SARS_CoV_2_haplotypes_Tajima_D_2020_time_series
@Author: Arnaud N'Guessan

## Overview
This repository contains the scripts and data used for analyzing the spread of SARS-CoV-2 haplotypes worldwide during the first and second waves (2020). 

## Dependencies
R (version 3.5.2+) packages: "ggplot2", "seqinr", "grid", "RColorBrewer", "randomcoloR", "gplots", "lmPerm", "ggpubr", "gridExtra", "RColorBrewer", "tidyr", "dendextend", "Cairo", "UpSetR", "parallel", "foreach", "doParallel", "infotheo", "igraph", "FD", "vegan", "indicspecies", "plyr", "lme4", "lmerTest", "MuMIn", "AICcmodavg", "EnvStats", "session"

## Order of script execution
The scripts "first_wave_Tajima_D_per_Continent_top3_haplotypes.r" and "second_wave_Tajima_D_per_Continent_top3_haplotypes.r", which calculate Tajima's D with a control for differences in sample size across haplotypes and time, should be executed before "Taj_D_time_series_SARS-CoV-2_top3_continent_Haplotypes.r", which allows to visualize the time series. 

1."first_wave_Tajima_D_per_Continent_top3_haplotypes.r"

You should run this script using the command "Rscript first_wave_Tajima_D_per_Continent_top3_haplotypes.r WORKSPACE NB_CORES THE_CONTINENT".


a) Arguments:


-->WORKSPACE: The absolute path of the repertory containing the data 


-->NB_CORES: Number of cpus for Tajima's D analysis. It corresponds to the number of months that will be analyzed in parallel. Thus, it should be <=7 for Wave 1.


-->THE_CONTINENT: The continent analyzed


b) Inputs (files you need to copy into the workspace): 


-->"Hussingroup_Inter_db_consensus_sequences.fasta": This file contains the consensus sequences in the database (This file is too large to be uploaded in Github. Contact the main author for access)


-->"HussinGroup_Interdb_metadata.csv": Metadata of the consensus sequences


c) Output: 


-->"Table_time_series_Taj_D_with_resamplings_first_wave_{THE_CONTINENT}_top3_haplotypes.csv"


-->"Taj_D_first_wave_{THE_CONTINENT}_top3_haplotypes_RSession.Rda"

2."second_wave_Tajima_D_per_Continent_top3_haplotypes.r"

You should run this script using the command "Rscript second_wave_Tajima_D_per_Continent_top3_haplotypes.r WORKSPACE NB_CORES THE_CONTINENT".

a) Arguments:

-->WORKSPACE: The absolute path of the repertory containing the data 

--> NB_CORES: Number of cpus for Tajima's D analysis. It corresponds to the number of months that will be analyzed in parallel. Thus, it should be <=5 for Wave 2.

--> THE_CONTINENT: The continent analyzed


b) Inputs (files you need to copy into the workspace): 

--> "Hussingroup_Inter_db_consensus_sequences.fasta": This file contains the consensus sequences in the database (This file is too large to be uploaded in Github. Contact the main author for access)

--> "HussinGroup_Interdb_metadata.csv": Metadata of the consensus sequences

c) Output: 

--> "Table_time_series_Taj_D_with_resamplings_second_wave_{THE_CONTINENT}_top3_haplotypes.csv" 

--> "Taj_D_second_wave_{THE_CONTINENT}_top3_haplotypes_RSession.Rda"

3."Taj_D_time_series_SARS-CoV-2_top3_continent_Haplotypes.r"

You should run this script using the command "Rscript Taj_D_time_series_SARS-CoV-2_top3_continent_Haplotypes.r WORKSPACE".

a) Arguments:

-->WORKSPACE: The absolute path of the repertory containing the data (same workspace as the previous scripts)

b) Inputs (files you need to copy into the workspace): 

-->The 12 Output files from "first_wave_Tajima_D_per_Continent_top3_haplotypes.r" and "second_wave_Tajima_D_per_Continent_top3_haplotypes.r"

-->"HussinGroup_Interdb_metadata.csv": Metadata of the consensus sequences

c) Output: 

-->The time series plots (Tajima's D and number of sequences per haplotype per month per continent)

-->the Rsession file
