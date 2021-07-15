
#' Title: "A database of pediatric drug effects to evaluate ontogenic mechanisms from child growth and development" study
#' 
#' Author details: Nicholas Giangreco
#' 
#' This script generates the table of Pediatric FAERS characteristics

# Purpose -----------------------------------------------------------------

#' Make table1 of pediatric FAERS data used for database generation
#'


# Setup -------------------------------------------------------------------


#install.packages("pacman")
pacman::p_load(tidyverse,data.table,tableone)

data_dir <- "../../data/"
img_dir <- "../../docs/imgs/"
seed = 0
set.seed(seed)

basename="database_generation_table1"
out_dir <- paste0(data_dir,basename,"/")


# Load data and covariates ------------------------------------------------

source("database_generation_load_data.R")


# Table1 ------------------------------------------------------------------

#https://cran.r-project.org/web/packages/tableone/vignettes/introduction.html

raw_data[,c(drug_col,rx_col),with=F] %>% unique() %>% nrow()

vars <- c(id_col,sex_col,stage_col,"reporter_qualification","polypharmacy")

sub <- 
raw_data[,
         vars,
         with=F
] %>% 
    unique() %>% 
    merge(atc_covariates_full_bin,by=id_col)

t1 <- 
    CreateTableOne(
        data = sub,
        vars = c(colnames(atc_covariates_full_bin)[2:ncol(atc_covariates_full_bin)],
                 setdiff(vars,id_col)),
        factorVars = colnames(atc_covariates_full_bin)[2:ncol(atc_covariates_full_bin)]
        )

t1Mat <- print(t1, quote = FALSE, noSpaces = TRUE, printToggle = FALSE,contDigits = 4)
## Save to a CSV file
write.csv(t1Mat, file = paste0(data_dir,basename,".csv"))

