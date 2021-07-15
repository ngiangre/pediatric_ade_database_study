
#' Title: "A database of pediatric drug effects to evaluate ontogenic mechanisms from child growth and development" study
#' 
#' Author details: Nicholas Giangreco
#' 
#' This script generates the random sample of drug-events for evaluating GAM
#' fit and generalization

# Purpose -----------------------------------------------------------------

#' To generate random drug-events for sampling GAM statistical fits and generalizations
#' 

# Setup -------------------------------------------------------------------

pacman::p_load(data.table,tidyverse)

data_dir <- "../../data/"
img_dir <- "../../docs/imgs/"
seed = 0
set.seed(seed)

basename <- paste0("database_generation_random_ade_candidates")

theme_set(theme_bw(base_size=16) + theme(text = element_text(face="bold")))


# Load data and covariates ------------------------------------------------

source("database_generation_load_data.R")

stage_knots <- 
    seq_along(category_levels[[stage_col]])

# choose and output candidate ADEs ---------------------------------------------------

candidate_ades <- 
    raw_data[,.(atc_concept_id,meddra_concept_id,atc_concept_name,meddra_concept_name)] %>% unique()

candidate_ades %>% 
    fwrite(paste0(data_dir,"database_generation_candidate_ades.csv"))

# random candidate ADEs ----------------------------------------------------------

n <- 2000
gt <- 50


ade_counts <- 
    raw_data %>% 
    .[,c(id_col,drug_col,rx_col,drug_col_name,rx_col_name),with=F] %>% 
    unique() %>% 
    .[,.N,c(drug_col,rx_col,drug_col_name,rx_col_name)]

candidate_ades_fit <- 
    ade_counts %>% 
    .[,c(drug_col,rx_col,drug_col_name,rx_col_name,"N"),with=F] %>% 
    sample_n(n)

candidate_ades_gen_setdiff <- 
    fsetdiff(
        ade_counts,
        candidate_ades[N>gt]
        ) %>% 
        .[N>gt] %>% 
        sample_n(n-nrow(candidate_ades[N>gt]))

candidate_ades_gen <- 
    bind_rows(
        candidate_ades_gen_setdiff,
        candidate_ades[N>gt]
    )



# Output random candidates ------------------------------------------------

candidate_ades_fit %>% 
    fwrite(paste0(data_dir,basename,"_for_fits.csv"))

candidate_ades_gen %>% 
    fwrite(paste0(data_dir,basename,"_for_generalizations.csv"))


# plot distributions ------------------------------------------------------

g <- candidate_ades_fit %>% 
    ggplot(aes(N)) + 
    geom_histogram(binwidth = 0.1) +
    scale_x_continuous(trans="log10") +
    xlab("Number of Reports") +
    ylab("Number of ADEs")
ggsave(paste0(img_dir,basename,"_fits_distribution.png"),g,width=5,height=4)

g <- candidate_ades_gen %>% 
    ggplot(aes(N)) + 
    geom_histogram(binwidth = 0.1) +
    scale_x_continuous(trans="log10") +
    xlab("Number of Reports") +
    ylab("Number of ADEs")
ggsave(paste0(img_dir,basename,"_generalizations_distribution.png"),g,width=5,height=4)
