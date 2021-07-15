
# PURPOSE -----------------------------------------------------------------

#' To load database generation data
#' 

# Setup -------------------------------------------------------------------

pacman::p_load(data.table,tidyverse,doParallel)

data_dir <- "../../data/"
img_dir <- "../../docs/imgs/"
seed = 0
set.seed(seed)

stages <- 
  c("term_neonatal","infancy",
    "toddler","early_childhood",
    "middle_childhood","early_adolescence",
    "late_adolescence")

# Load data (dts) ---------------------------------------------------------------

candidate_ades <- fread(paste0(data_dir,"database_generation_candidate_ades.csv"))

database_base="database_generation_base/"
database_covariate_adjusted="database_generation_with_covariates/"

dbs=list("base" = database_base,"covariate_adjusted" = database_covariate_adjusted)

dts <- NULL
for(d in seq_along(dbs)){
  tmp_dir <- paste0(data_dir,dbs[[d]])
  files <- list.files(tmp_dir)
  dt <- foreach(x=files,
                .combine = "rbind",
                .errorhandling = "remove") %dopar% {
                  fread(paste0(tmp_dir,x))
                }
  dt$database <- names(dbs)[d]
  dts <- rbind(dts,dt)
}

# db_base <-
#     fread(paste0(data_dir,"database_generation_base.csv"),nThread=cores)
# db_base$database <- "base"
# db_ca <-
#     fread(paste0(data_dir,"database_generation_with_covariates.csv"),nThread=cores)
# db_ca$database <- "covariate_adjusted"
# 
# dts <-
#     bind_rows(
#         db_base,
#         db_ca
#     )

common_ades <- 
    dts[,.(database,ade)] %>% 
    unique() %>% 
    .[,.N,ade] %>% 
    .[N==2,unique(ade)]

sample_ades <- sample(common_ades,min(10000,length(common_ades)),replace=F)

dts$nichd <- factor(dts$nichd,levels=stages)

narrow_ades <- 
  dts %>% 
  .[ade %in% common_ades,
    .(ade,database,nichd,gam_score_90mse,gam_score_90pse)] %>% 
  unique() %>% 
  .[gam_score_90pse<5 &
      gam_score_90mse> -5,
    .N,.(database,ade)] %>% 
  .[N==7,.(database,ade)] %>% 
  fsetdiff(dts[gam_score<0,
               .N,.(database,ade)] %>% 
             .[N==7,.(database,ade)] %>% 
             unique(),
           all=T) %>% 
  fsetdiff(dts[gam_score<1 &
                 gam_score> -1,
               .N,.(database,ade)] %>% 
             .[N==7,.(database,ade)] %>% 
             unique(),
           all=T) %>% 
  fintersect(dts[gam_score_90mse>1,
                 .(database,ade)] %>% 
               unique(),
             all=T) %>% 
  fintersect(dts[gam_score>1,
                 .(database,ade)] %>% 
               unique(),
             all=T) %>% 
  .[,unique(ade)]

sig_ades <- dts[database=="covariate_adjusted" & gam_score_90mse>0,unique(ade)]

null_dts <- fread(paste0(data_dir,"database_generation_with_covariates_null_simulation.csv"))

null_dts[,.(null_99 = quantile(gam_score,c(0.99))),nichd] %>% 
  .[,.(nichd = factor(nichd,levels=stages),null_99)] %>% 
  .[order(nichd)]

sig_null_ades <- 
merge(
  dts[database=="covariate_adjusted"],
  null_dts[,.(null_99 = quantile(gam_score,c(0.99))),nichd],
  by="nichd"
) %>% 
  .[gam_score_90mse>null_99,unique(ade)]
