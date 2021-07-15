#' Title: "A database of pediatric drug effects to evaluate ontogenic mechanisms from child growth and development" study
#' 
#' Author details: Nicholas Giangreco
#' 
#' This script generates the data for the chosen covariate GAM specification on 
#' a sample of random drug-events

# Purpose -----------------------------------------------------------------

#' 
#' To trial different GAM specifications on the pediatric reference set to prototype the output for the database. 

# Setup -------------------------------------------------------------------


#install.packages("pacman")
pacman::p_load(tidyverse,data.table,mgcv,doParallel,progress,ROCR)

data_dir <- "../../data/"
img_dir <- "../../docs/imgs/"
seed = 0
set.seed(seed)

registerDoParallel(cores=5)

basename="database_generation_with_covariates_null_simulation"


# Load data and covariates ------------------------------------------------

source("database_generation_load_data.R")

stage_knots <- 
  seq_along(category_levels[[stage_col]])


# randomize drug and event reporting --------------------------------------

drugs <- raw_data[,unique(get(drug_col))]
events <- raw_data[,unique(get(rx_col))]

base_drug_dt <- 
  raw_data[,c(id_col,drug_col),with=F] %>% unique()
#base_drug_dt[,.N,id_col] %>% .[,mean(N)]
base_drug_dt[,(drug_col) := sample(get(drug_col),nrow(base_drug_dt),replace=F)]
#base_drug_dt[,.N,id_col] %>% .[,mean(N)]

base_rx_dt <- 
  raw_data[,c(id_col,rx_col),with=F] %>% unique()
#base_rx_dt[,.N,id_col] %>% .[,mean(N)]
base_rx_dt[,(rx_col) := sample(get(rx_col),nrow(base_rx_dt),replace=F)]
#base_rx_dt[,.N,id_col] %>% .[,mean(N)]

base_drug_rx_dt <- 
  merge(base_drug_dt,base_rx_dt,allow.cartesian = T) %>% 
  merge(raw_data[,
                 c(drug_col,rx_col,drug_col_name,rx_col_name),with=F] %>% 
          unique(),
        by=c(drug_col,rx_col))

raw_data <- 
merge(
  base_drug_rx_dt,
  raw_data[,c(stage_col,sex_col,"reporter_qualification",
              "receive_date",names(atc_covariates),
              "polypharmacy"),with=F] %>% 
    unique(),
  by=id_col,
  allow.cartesian = T
)

# make candidate ADEs -----------------------------------------------------

candidate_ades <- 
  raw_data[,c(drug_col,rx_col,drug_col_name,rx_col_name),with=F] %>% 
  unique() %>% 
  sample_n(1e4)

# functions ---------------------------------------------------------------

source("database_generation_functions.R")

# GAM modeling ------------------------------------------------------------


x_cols=c(names(category_levels),"polypharmacy")
ade_data <- 
  make_ade_data(
    x_cols=x_cols,
    category_levels = category_levels
  ) %>% 
  .[order(get(stage_col))]

dts <- NULL
minrange = 1
maxrange = 500
cat(maxrange," candidates\n")
t0 = Sys.time()
cat("Start:",format(Sys.time(), "%a %b %d %X %Y"),"\n")

range <- minrange:maxrange
step=0.02
start_cuts <-  seq(step,1-step,step)
end_cuts <-  seq(step,1,step)
starts <- sapply(start_cuts,function(x)range[floor(x*length(range))])
starts <- c(min(range),starts)
ends <- sapply(end_cuts,function(x)range[floor(x*length(range))])
for(i in seq_along(end_cuts)){
  
  start <- ifelse(starts[i]==min(range),starts[i],starts[i]+1)
  end <- ends[i]
  cat("Indices: ",start," to ",end,"\n")
  
  system.time({dt = foreach(x=start:end,.combine="rbind",.errorhandling='remove') %dopar% {

    p <- candidate_ades[x]

    drug_id <- p[[drug_col]]
    rx_id <- p[[rx_col]]
    drug_name <- p[[drug_col_name]]
    rx_name <- p[[rx_col_name]]

    dat <- set_ade_pair(ade_data,drug_id,rx_id,drug_name,rx_name)

    time_ <- system.time(bam_mod <- compute_bam_with_reporter_date_sex_effects_atc(dat))
    bam_cs <- process_gam_model(
      bam_mod,
      dat,
      drug_id,rx_id,drug_name,rx_name,
      knots=stage_knots,
      coef_str="s(nichd):D.",time_[3]
    )

    merge(
      bam_cs,
      compute_prr(drug_id,
                  rx_id,
                  drug_name,
                  rx_name,dat)[,
                               .(atc_concept_id,meddra_concept_id,nichd,
                                 PRR,PRR_90mse,PRR_95mse,PRR_99mse,PRR_999mse,
                                 PRR_90pse,PRR_95pse,PRR_99pse,PRR_999pse)
                               ]
      )
  }})
  
  dts <- bind_rows(dts,dt)
  
  cat((i/length(end_cuts))*100,"% complete\n",format(Sys.time(), "%a %b %d %X %Y"),"\n")
}
t1 = Sys.time()
cat("End:",format(Sys.time(), "%a %b %d %X %Y"),"\n")
cat("\n",round(as.numeric(difftime(t1,t0,units="mins")),2)," minutes\n")

dts %>% 
  fwrite(paste0(data_dir,basename,".csv"))