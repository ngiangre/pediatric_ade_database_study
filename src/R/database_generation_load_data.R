#' Title: "A database of pediatric drug effects to evaluate ontogenic mechanisms from child growth and development" study
#' 
#' Script author details: Nicholas Giangreco
#' 
#' This script loads the data and covariates for drug-events
#' 
# PURPOSE -----------------------------------------------------------------

#' Load data and covariates for database generation
#' 

# Setup -------------------------------------------------------------------

pacman::p_load(tidyverse,data.table)

# load data ---------------------------------------------------------------


data_dir <- "../../data/"

drug_col = "atc_concept_id"
drug_col_name = "atc_concept_name"
rx_col = "meddra_concept_id"
rx_col_name = "meddra_concept_name"
stage_col = "nichd"
age_col = "age"
sex_col="sex"
dz_col = "drug_indication"
id_col="safetyreportid"

raw_data <- 
    data.table::fread(paste0(data_dir,"database_generation_pediatric_FAERS.csv.gz"))

raw_data[,`:=`(ade = paste0(get(drug_col),"_",get(rx_col)))]

# adding covariates -------------------------------------------------------

category_levels = list()

category_levels[[stage_col]] =
    c("term_neonatal","infancy",
      "toddler","early_childhood",
      "middle_childhood","early_adolescence",
      "late_adolescence")

raw_data[,stage_col] <- 
  factor(raw_data[,get(stage_col)],levels=category_levels[[stage_col]])

category_levels[[sex_col]] =
    c("Male","Female")

category_levels[["reporter_qualification"]] =
    c("Consumer or non-health professional",
      "Lawyer","Other health professional",
      "Pharmacist","Physician")

dates <- raw_data$receive_date %>% unique()

category_levels[["receive_date"]] <-
    dates[order(dates)]


# atc covariates level 4-------------------------------------------------------


atc_covariates <- 
  merge(
        raw_data[,.(safetyreportid,atc_concept_id)] %>% 
            unique(),
        data.table::fread(paste0(data_dir,"compound_drugbank05/drug_atc_codes_rxnorm_joined.csv")) %>% 
            .[,.(atc_concept_id,code_4)] %>% 
            unique()
    ) %>% 
    data.table::dcast(safetyreportid ~ code_4,value.var="code_4",
                      fun.aggregate = length)
colnames(atc_covariates)[2:ncol(atc_covariates)] <- 
    paste0("X",colnames(atc_covariates)[2:ncol(atc_covariates)])

atc_covariates_bin <- 
  merge(
    raw_data[,.(safetyreportid,atc_concept_id)] %>% 
      unique(),
    data.table::fread(paste0(data_dir,"compound_drugbank05/drug_atc_codes_rxnorm_joined.csv")) %>% 
      .[,.(atc_concept_id,code_4)] %>% 
      unique()
  ) %>% 
  data.table::dcast(safetyreportid ~ code_4,value.var="code_4",
                    fun.aggregate = function(x){as.integer(length(x)>0)})
colnames(atc_covariates_bin)[2:ncol(atc_covariates_bin)] <- 
  paste0("X",colnames(atc_covariates_bin)[2:ncol(atc_covariates_bin)])

atc_covariates_full <- 
  merge(
    raw_data[,.(safetyreportid,atc_concept_id)] %>% 
      unique(),
    data.table::fread(paste0(data_dir,"compound_drugbank05/drug_atc_codes_rxnorm_joined.csv")) %>% 
      .[,.(atc_concept_id,level_4)] %>% 
      unique()
  ) %>% 
  data.table::dcast(safetyreportid ~ level_4,value.var="level_4",
                    fun.aggregate = length)

atc_covariates_full_bin <- 
  merge(
    raw_data[,.(safetyreportid,atc_concept_id)] %>% 
      unique(),
    data.table::fread(paste0(data_dir,"compound_drugbank05/drug_atc_codes_rxnorm_joined.csv")) %>% 
      .[,.(atc_concept_id,level_4)] %>% 
      unique()
  ) %>% 
  data.table::dcast(safetyreportid ~ level_4,value.var="level_4",
                    fun.aggregate = function(x){as.integer(length(x)>0)})


# atc covariates level 3-------------------------------------------------------


atc_covariates_l3 <- 
  merge(
    raw_data[,.(safetyreportid,atc_concept_id)] %>% 
      unique(),
    data.table::fread(paste0(data_dir,"compound_drugbank05/drug_atc_codes_rxnorm_joined.csv")) %>% 
      .[,.(atc_concept_id,code_3)] %>% 
      unique()
  ) %>% 
  data.table::dcast(safetyreportid ~ code_3,value.var="code_3",
                    fun.aggregate = length)
colnames(atc_covariates_l3)[2:ncol(atc_covariates_l3)] <- 
  paste0("X",colnames(atc_covariates_l3)[2:ncol(atc_covariates_l3)])

atc_covariates_l3_bin <- 
  merge(
    raw_data[,.(safetyreportid,atc_concept_id)] %>% 
      unique(),
    data.table::fread(paste0(data_dir,"compound_drugbank05/drug_atc_codes_rxnorm_joined.csv")) %>% 
      .[,.(atc_concept_id,code_3)] %>% 
      unique()
  ) %>% 
  data.table::dcast(safetyreportid ~ code_3,value.var="code_3",
                    fun.aggregate = function(x){as.integer(length(x)>0)})
colnames(atc_covariates_l3_bin)[2:ncol(atc_covariates_l3_bin)] <- 
  paste0("X",colnames(atc_covariates_l3_bin)[2:ncol(atc_covariates_l3_bin)])

atc_covariates_l3_full <- 
  merge(
    raw_data[,.(safetyreportid,atc_concept_id)] %>% 
      unique(),
    data.table::fread(paste0(data_dir,"compound_drugbank05/drug_atc_codes_rxnorm_joined.csv")) %>% 
      .[,.(atc_concept_id,level_3)] %>% 
      unique()
  ) %>% 
  data.table::dcast(safetyreportid ~ level_3,value.var="level_3",
                    fun.aggregate = length)

atc_covariates_l3_full_bin <- 
  merge(
    raw_data[,.(safetyreportid,atc_concept_id)] %>% 
      unique(),
    data.table::fread(paste0(data_dir,"compound_drugbank05/drug_atc_codes_rxnorm_joined.csv")) %>% 
      .[,.(atc_concept_id,level_3)] %>% 
      unique()
  ) %>% 
  data.table::dcast(safetyreportid ~ level_3,value.var="level_3",
                    fun.aggregate = function(x){as.integer(length(x)>0)})


# polypharmacy ------------------------------------------------------------

raw_data <- 
  merge(
    raw_data,
    raw_data[,
             c(id_col,drug_col),
             with=F
             ] %>% 
      unique() %>% 
      .[,.(polypharmacy = .N),id_col],
    by=id_col
    )
