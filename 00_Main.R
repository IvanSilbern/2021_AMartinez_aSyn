

# Sets of scripts in support to the paper
# "Low-expressing synucleinopathy mouse models based on oligomer formation and C-terminal truncation of α-synuclein"
# https://www.frontiersin.org/articles/10.3389/fnins.2021.643391/full
# by A. Martinez & I. Silbern et al
# 
# run the scripts one after another
# to obtain quantifications of alpha-synuclein peptides
# input quantification tables from Skyline are located in /Skyline_reults folder
# intermediate tables are saved in /temp folder
# figures will be saved in /plots folder
# result table is saved in /tables folder

source("01_Analysis_Hippoc.R")
source("02_Analysis_Str_OB_SN.R")
source("03_Dilutions_aSyn15N.R")
source("04_prepare_figures_Hippoc.R")
source("05_prepare_figures_Str_OB_SN.R")
source("06_prepare_report.R")