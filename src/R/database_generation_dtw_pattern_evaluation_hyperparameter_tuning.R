#' Title: "A database of pediatric drug effects to evaluate ontogenic mechanisms from child growth and development" study
#' 
#' Author details: Nicholas Giangreco
#' 
#' This script generates the data of our
#' data-driven clustering hyperparameter tuning

# Purpose -----------------------------------------------------------------

#' To quantify pattern homogeneity using synthetic drug-event patterns
#' on munnin (multicore machine)
#' 

# Setup -------------------------------------------------------------------

pacman::p_load(data.table,tidyverse,doParallel,dtwclust)

data_dir <- "../../data/"
img_dir <- "../../docs/imgs/"
basename <- "database_generation_dtw_bootstrap_pattern_evaluation_"
seed = 0
set.seed(seed)

registerDoParallel(cores=50)

stages <- 
    c("term_neonatal","infancy",
      "toddler","early_childhood",
      "middle_childhood","early_adolescence",
      "late_adolescence")

# Load functions ----------------------------------------------------------

source("database_generation_functions.R")

# Load synthetic drug-events with predetermined patterns ------------------

syn_inc <- 
    fread(paste0(data_dir,"sim_positive_nichd_increase_data.csv")) %>% 
    .[,.(ade,nichd,gam_score,spikein)] %>% 
    unique()

syn_dec <- 
    fread(paste0(data_dir,"sim_positive_nichd_decrease_data.csv")) %>% 
    .[,.(ade,nichd,gam_score,spikein)] %>% 
    unique()

syn_plat <- 
    fread(paste0(data_dir,"sim_positive_nichd_plateau_data.csv")) %>% 
    .[,.(ade,nichd,gam_score,spikein)] %>% 
    unique()

syn_data <- 
    bind_rows(
        syn_inc,
        syn_dec,
        syn_plat
    )

syn_data$nichd <- factor(syn_data$nichd,levels=stages)
norm_syn_data <- 
    normalize_data(syn_data,groups=c("spikein"))

syn_unif <- 
    fread(paste0(data_dir,"sim_positive_nichd_uniform_data.csv")) %>% 
    .[,.(ade,nichd,gam_score,spikein)] %>% 
    unique()

syn_lst=list("increase"=syn_inc,"decrease"=syn_dec,"plateau"=syn_plat)
ddiff_dts <- NULL
for(i in seq_along(syn_lst)){
    dat=syn_lst[[i]]
    spikein = dat[,unique(spikein)]
    dat_merge = merge(dat,syn_unif,by=c("ade","nichd"))
    ddiff <- 
        dat_merge[,.(ade,nichd,ddiff = gam_score.x - gam_score.y)] 
    ddiff_dts <- 
        bind_rows(
            ddiff_dts,
            merge(
                dat,
                merge(
                    syn_unif[,.(ade,nichd,gam_score_unif = gam_score)],
                    ddiff,by=c("ade","nichd")
                ),
                by=c("ade","nichd")
            )
        )
    
}

ddiff_dts$nichd <- factor(ddiff_dts$nichd,levels=stages)

norm_ddiff <- 
  bind_rows(
    normalize_data(ddiff_dts[spikein=="increase"],score="ddiff",groups=c("spikein")),
    normalize_data(ddiff_dts[spikein=="decrease"],score="ddiff",groups=c("spikein")),
    normalize_data(ddiff_dts[spikein=="plateau"],score="ddiff",groups=c("spikein"))
  ) %>% 
  na.omit() # many are 0

norm_ddiff$nichd <- factor(norm_ddiff$nichd,levels=stages)

norm_ddiff[,rank:=frank(norm),.(ade,spikein)]
norm_ddiff[,str_pattern := paste0(rank,collapse = ""),.(ade,spikein)]

#syn_data[ade %in% norm_ddiff[is.na(norm),ade],.(S = sum(DE)),.(ade,spikein)] %>% .[,mean(S),spikein]

plat_ades <- 
  norm_ddiff[
    spikein=="plateau",
    .(ade,nichd,rank)
  ] %>% 
  unique() %>% 
  dcast(ade ~ nichd) %>% 
  .[(term_neonatal<infancy) & (early_adolescence>late_adolescence),ade]

clean_norm_ddiff <- 
  bind_rows(
    norm_ddiff[str_pattern=="1234567" & spikein=="increase"],
    norm_ddiff[str_pattern=="7654321" & spikein=="decrease"],
    norm_ddiff[
      (ade %in% plat_ades & spikein=="plateau")]
  )


# load GAM data -----------------------------------------------------------

source("database_generation_load_GAM_data.R")

# give synthetic data unique IDs for evaluation ---------------------------

tmp <- 
  clean_norm_ddiff %>% 
  .[,.(ade,spikein)] %>% 
  unique()

norm_syn_ades <- 
  tmp %>% 
  .[,.(ade,spikein,id = paste0("X",1:nrow(tmp)))] %>% 
  merge(
    clean_norm_ddiff[,.(ade,nichd,norm,spikein)],
    by=c("ade","spikein"),
    allow.cartesian = T
  ) %>% 
  .[,.(ade = id,nichd,norm,spikein)] %>% 
  .[order(ade)]

# Normalize ade gam scores -----------------------------

norm_ades <- 
  normalize_data(dts[database=="covariate_adjusted" &
                       ade %in% common_ades],score="gam_score",groups=c())

# hyperparameter tuning DTW -----------------------------------------------

cat("\n## hyperparameter tuning DTW ## \n")
t0 = Sys.time()
cat("Start:",format(Sys.time(), "%a %b %d %X %Y"),"\n")

distances=c("dtw_basic","lbi","sbd","gak")
centroids=c("mean","shape","dba","pam","fcm","fcmdd")
types=c("partitional","fuzzy")
ks = c(2,3,4,5,6)
param_grid <- NULL
for(di in distances){
    for(ce in centroids){
        for(ty in types){
            for(k in ks){
                dt <- data.table(distance = di, centroid = ce,type = ty,k=k)
                param_grid <- 
                    bind_rows(
                        param_grid,
                        dt
                    )
            }
        }
    }
}
param_grid$param_grid_ind <- 1:nrow(param_grid)

dts <- NULL

Nsamp = 1e4
N=10

for(i in 1:nrow(param_grid)){
  
  cat("Parameters:\n")
  print(param_grid[i])
  t1 = Sys.time()
  cat(format(Sys.time(), "%a %b %d %X %Y"),"\n")
  
  boots <- 
    foreach(j=1:N,
            .combine = "rbind",
            .errorhandling = "remove") %dopar% {
              
              set.seed(j)
              RcppParallel::setThreadOptions(1L)
              
              X <- 
                bind_rows(
                  norm_ades %>% 
                    .[ade %in% sample(sig_null_ades,
                                      min(Nsamp,length(sig_null_ades)),
                                          replace=F)],
                  norm_syn_ades
                ) %>% 
                dcast(ade ~ nichd,value.var="norm") %>% 
                merge(norm_syn_ades[,.(ade,spikein)] %>% unique(),by="ade",all.x=T)
              
              time_ <- system.time(
                dtw <- 
                  tsclust(X[,stages,with=F],
                          type=param_grid[i,type],
                          centroid=param_grid[i,centroid],
                          distance=param_grid[i,distance],
                          preproc = NULL,
                          k = param_grid[i,k],
                          seed=0)
              )
              X[,"cluster":=dtw@cluster]
              X[,"param_grid_ind":=i]
              
              param_dt <- 
                evaluate_patterns_in_parameter_sets(X)
              param_dt[,"time":=time_[3]]
              param_dt[,"bootstrap":=j]
              
              param_dt
              
            }
  
  dts <- 
    bind_rows(dts,
              boots
    )
}

t1 = Sys.time()
cat("End:",format(Sys.time(), "%a %b %d %X %Y"),"\n")
cat("\n",round(as.numeric(difftime(t1,t0,units="mins")),2)," minutes\n")

param_grid %>% 
  fwrite(paste0(data_dir,basename,"hyperparameter_grid.csv"))
dts %>% 
  fwrite(paste0(data_dir,basename,"tuning_results.csv"))
