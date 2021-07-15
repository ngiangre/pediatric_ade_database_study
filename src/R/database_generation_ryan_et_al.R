#' Title: "A database of pediatric drug effects to evaluate ontogenic mechanisms from child growth and development" study
#' 
#' Author details: Nicholas Giangreco
#' 
#' This script generates the enrichment data of our
#' identified significant drug-events in an ADE reference set

# Purpose -----------------------------------------------------------------

#' To find in pediatric FAERS the AEs and drugs (maspped to ATC) within the Ryan et al. reference set
#' 


# Setup -------------------------------------------------------------------

pacman::p_load(tidyverse,data.table,doParallel)

data_dir <- "../../data/"
img_dir <- "../../docs/imgs/"
basename <- "database_generation_ryan_et_al_"
seed = 0
set.seed(seed)

registerDoParallel(cores=5)

theme_set(theme_bw(base_size=16) + theme(text = element_text(face="bold")))

# Load data and covariates ------------------------------------------------

source("database_generation_load_data.R")


# Load GAM data ------------------------------------------------

source("database_generation_load_GAM_data.R")

# load raw ryan et al data ----------------------------------------------------

meddra_map <- 
    fread(paste0(data_dir,"ryan_to_meddra_concept_map.csv"))

meddra_to_rxnorm_map <- 
    fread(paste0(data_dir,"ryan_to_rxnorm_concept_map.csv"))


# map rxnorm to atc and join-------------------------------------------------------

atc_map <- 
    fread(paste0(data_dir,"compound_drugbank05/drug_atc_codes_rxnorm_joined.csv"))

merged <- 
merge(
    meddra_to_rxnorm_map[,
                         .(condition_name,
                           ryan_concept_id,
                           rxnorm_concept_id = as.character(rxnorm_concept_id),
                           control = ifelse(positive_control==1,"positive","negative"))
                         ],
    atc_map[,
            .(rxnorm_concept_id = as.character(rxnorm_concept_id),
               atc_concept_id)
            ] %>% 
        unique(),
    by="rxnorm_concept_id"
    ) %>% 
    merge(
        meddra_map,
        by="ryan_concept_id",
        allow.cartesian = T
    ) %>% 
    .[,.(condition_name,
         meddra_concept_id,
         atc_concept_id,
         control)
      ] %>% 
    unique()


# Output ------------------------------------------------------------------

merged %>% 
    fwrite(paste0(data_dir,"database_generation_ryan_et_al.csv"))


# load processed ryan et al data ------------------------------------------

merged <-  
  fread(paste0(data_dir,"database_generation_ryan_et_al.csv"))

# summarizing -------------------------------------------------------------

nrow(merged)
merged[,.N,.(control)]
merged[,.N,.(condition_name)]
merged[,.N,.(condition_name,control)]

raw_data_merged <- 
raw_data[,.(safetyreportid,atc_concept_id,meddra_concept_id)] %>% 
    .[,.N,c(rx_col,drug_col)] %>% 
    merge(
        merged,
        by=c("atc_concept_id","meddra_concept_id")
        )

raw_data_merged[,.N,.(control)]
raw_data_merged[,.N,.(condition_name,control)]


g <- raw_data_merged %>% 
    ggplot(aes(N)) +
    geom_histogram(fill="gray",color="black",binwidth = 1) +
    scale_x_continuous(breaks=c(1,5,10,25,50,75,100)) +
    scale_y_continuous(breaks=c(1,5,10,15,20,25,30)) +
    facet_wrap(~condition_name,ncol=2) +
    xlab("Number of reports") +
    ylab("Number of adverse drug events")
ggsave(paste0(img_dir,"ryan_et_al_childhood_ade_report_distribution.png"),width=8,height=4)

# average GAM scores per stage per positive/negative controls ------------

raw_data_merged$ade <- 
  paste0(raw_data_merged$atc_concept_id,"_",raw_data_merged$meddra_concept_id)

ades <- 
  dts[database=="covariate_adjusted" & gam_score>0,ade]

tmp <- 
merge(
  dts[ade %in% ades &
    database=="covariate_adjusted",
    .(ade,nichd,gam_score,gam_score_90mse)
    ],
  raw_data_merged,
  by="ade",
  allow.cartesian = T
) %>% 
  dcast(ade + nichd + condition_name ~ control,
        value.var = "gam_score") 

test_dts <- NULL
for(st in category_levels$nichd){
  a <- tmp[nichd==st,na.omit(positive)] %>% as.numeric()
  b <- tmp[nichd==st,na.omit(negative)] %>% as.numeric()
  test <- ks.test(a,b,alternative = "greater")
  dt <- data.table(nichd = st, pvalue = test$p.value,
                   positive_mean = mean(a),negative_mean = mean(b)
                   )
  test_dts <- 
    rbind(test_dts,dt)
}

tmp <- 
  merge(
    dts[ade %in% ades &
          database=="covariate_adjusted",
        .(ade,nichd,gam_score,gam_score_90mse)
    ],
    raw_data_merged,
    by="ade",
    allow.cartesian = T
  )

tmp[gam_score_90mse>0,
    .N,
    .(control,condition_name)
    ] %>% 
  dcast(condition_name ~ control,value.var = "N")

prop.test(249,38+249)$p.value
prop.test(9,4+9)$p.value
prop.test(27,38+27)$p.value

prop.test(27+9+249,38+4+38+27+9+249)$p.value

# Sig vs total positive vs negative ADEs ----------------------------------

raw_data_merged$ade <- 
  paste0(raw_data_merged$atc_concept_id,"_",raw_data_merged$meddra_concept_id)

tmp <- 
  dts[
    database=="covariate_adjusted",
    .(ade,nichd,gam_score_90mse)
  ] %>% 
  unique() %>% 
  merge(
    raw_data_merged,
    by="ade"
  )

tmp[,length(unique(ade))]

merge(
  tmp[,.(total = length(unique(ade))),.(control)], 
  tmp[gam_score_90mse>0,
      .(sig = length(unique(ade))),
      .(control)]
)

a = 74
b = 238
c = 20
d = 74
or = (a / b) / (c / d)
se = sqrt( (1/a) + (1/b) + (1/c) + (1/d))
or 
exp( log(or) - 1.645*se )
exp( log(or) + 1.645*se )
