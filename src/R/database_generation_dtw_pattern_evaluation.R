#' Title: "A database of pediatric drug effects to evaluate ontogenic mechanisms from child growth and development" study
#' 
#' Author details: Nicholas Giangreco
#' 
#' This script generates the data and plots from our
#' data-driven clustering of drug-events

# Purpose -----------------------------------------------------------------

#' To strategize and quantify pattern homogeneity using synthetic drug-event patterns
#' 

# Setup -------------------------------------------------------------------

pacman::p_load(data.table,tidyverse,doParallel,dtwclust)

data_dir <- "../../data/"
img_dir <- "../../docs/imgs/"
basename <- "database_generation_dtw_pattern_evaluation_"
seed = 0
set.seed(seed)

registerDoParallel(cores=5)

theme_set(theme_bw(base_size=16) + theme(text = element_text(face="bold")))

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
    .[,.(ade,nichd,gam_score,DE,Tij,spikein)] %>% 
    unique()

syn_dec <- 
    fread(paste0(data_dir,"sim_positive_nichd_decrease_data.csv")) %>% 
    .[,.(ade,nichd,gam_score,DE,Tij,spikein)] %>% 
    unique()
syn_dec %>% dcast(ade ~ factor(nichd,levels=stages),value.var="gam_score")

syn_plat <- 
    fread(paste0(data_dir,"sim_positive_nichd_plateau_data.csv")) %>% 
    .[,.(ade,nichd,gam_score,DE,Tij,spikein)] %>% 
    unique()
syn_plat %>% dcast(ade ~ factor(nichd,levels=stages),value.var="gam_score")

syn_data <- 
    bind_rows(
        syn_inc,
        syn_dec,
        syn_plat
    )

syn_data$nichd <- factor(syn_data$nichd,levels=stages)
norm_syn_data <- 
    normalize_data(syn_data,groups=c("spikein"))

g <- norm_syn_data %>% 
    ggplot(aes(
        factor(nichd,levels=stages),norm,group=ade)
    ) +
    geom_path(alpha=0.1) +
    facet_grid(spikein~.) +
    xlab("") +
    ylab("Normalized score") +
    theme(
        axis.text.x = element_text(angle=45,vjust=1,hjust=1)
    )
ggsave(paste0(img_dir,basename,"synthetic_ade_norm_scores.png"),g,height=7,width=5)


syn_unif <- 
    fread(paste0(data_dir,"sim_positive_nichd_uniform_data.csv")) %>% 
    .[,.(ade,nichd,gam_score,spikein)] %>% 
    unique()

syn_lst=list("increase"=syn_inc,"decrease"=syn_dec,"plateau"=syn_plat)
ddiff_dts <- NULL
for(i in seq_along(syn_lst)){
    dat=syn_lst[[i]]
    spikein = names(syn_lst)[i]
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

g <- clean_norm_ddiff %>% 
    ggplot(aes(
    factor(nichd,levels=stages),norm,group=ade)) +
    geom_path(alpha=0.1) +
    facet_grid(spikein~.) +
    xlab("") +
    ylab("Normalized difference from uniform score") +
    theme(
      strip.background = element_blank(),
        axis.text.x = element_text(angle=45,vjust=1,hjust=1)
    )
ggsave(paste0(img_dir,basename,"synthetic_ade_diff_from_uniform_norm_scores.png"),g,height=7,width=5)


# Load GAM data -----------------------------------------------------------

source("database_generation_load_GAM_data.R")


# Normalizing GAM data ----------------------------------------------------

norm_ades <- 
  normalize_data(dts[database=="covariate_adjusted" &
                       ade %in% common_ades],score="gam_score",groups=c())


# GAM scores across stages ------------------------------------------------

g <- norm_ades[ade %in% sample_ades] %>%
  ggplot(aes(nichd,norm,group=ade)) +
  geom_line(alpha=0.1) +
  xlab("") +
  ylab("Normalized GAM score") +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"gam_score_minmax_score_patterns_across_stages.png"),g,
       width=6,height=5)

# Proportion of distinct patterns - increase, decrease, other -------------

norm_ades$nichd <- factor(norm_ades$nichd,levels=stages)

norm_ades[,rank:=frank(norm),.(ade)]

rank <- 
  norm_ades %>% 
  .[order(nichd)] %>% 
  dcast(ade~nichd,value.var="rank")
rank[, str_pattern:=do.call(paste0,.SD), .SDcols=-1]

tmp = data.table(
  dynamic = c("increase","decrease","other"),
  N = c(
    rank[,.(ade,str_pattern)][,.N,str_pattern][str_pattern=="1234567",N],
    rank[,.(ade,str_pattern)][,.N,str_pattern][str_pattern=="7654321",N],
    nrow(rank) - rank[,.(ade,str_pattern)][,.N,str_pattern][1:2,sum(N)]
  )
)
tmp$prop <- tmp$N/nrow(rank)

g <- tmp %>% 
  ggplot(aes(forcats::fct_reorder(dynamic,prop,.desc = T),prop)) +
  geom_bar(stat="identity",color="black",fill="gray") +
  xlab("") +
  ylab("Percent of drug-events") +
  scale_y_continuous(labels=scales::percent) +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"percent_pattern_types.png"),g,width=5,height=5)

# Give synthetic data unique IDs for evaluation ---------------------------

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

norm_syn_ades[,length(unique(ade)),spikein]

# Make X for DTW hyperparameter tuning -----------------------------

Nsamp = 1e4

X <- 
  bind_rows(
    norm_ades %>% 
      .[ade %in% sample(common_ades,Nsamp,replace=F)],
    norm_syn_ades
  ) %>% 
  dcast(ade ~ nichd,value.var="norm") %>% 
  merge(norm_syn_ades[,.(ade,spikein)] %>% unique(),by="ade",all.x=T)

X[,.N,spikein]

# SAMPLE hyperparameter tuning DTW -----------------------------------------------

cat("\n## hyperparameter tuning DTW ## \n")
t0 = Sys.time()
cat("Start:",format(Sys.time(), "%a %b %d %X %Y"),"\n")

distances=c("dtw_basic","sbd","gak")
centroids=c("mean","pam","dba","shape")
types=c("partitional","hierarchical")
ks = c(2,3,4,5,6,7,8,9,10)
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
minrange = 1
maxrange = nrow(param_grid)
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
  
eval_dt <- 
    foreach(j=start:end,
            .combine = "rbind",
            .errorhandling = "remove") %dopar% {
                time_ <- system.time(
                    dtw <- 
                        tsclust(X[,norm_ades[,levels(nichd)],with=F],
                                type=param_grid[j,type],
                                centroid=param_grid[j,centroid],
                                distance=param_grid[j,distance],
                                preproc = NULL,
                                k = param_grid[j,k],
                                seed=seed)
                )
                X[,"cluster":=dtw@cluster]
                X[,"time":=time_[3]]
                X[,"param_grid_ind":=j]
                
                X
            }

dts <- 
  bind_rows(dts,
            eval_dt
            )
}

param_grid %>% 
    fwrite(paste0(data_dir,basename,"hyperparameter_grid.csv"))
dts %>% 
    fwrite(paste0(data_dir,basename,"tuning_cluster_results.csv.gz"))

t1 = Sys.time()
cat("End:",format(Sys.time(), "%a %b %d %X %Y"),"\n")
cat("\n",round(as.numeric(difftime(t1,t0,units="mins")),2)," minutes\n")

# Load hyperparameter tuning results -----------------------------------------------

param_grid <- 
    fread(paste0(data_dir,basename,"hyperparameter_grid.csv"))
eval_dt <- 
    fread(paste0(data_dir,basename,"tuning_cluster_results.csv.gz"))
    
eval_param_dt <- 
eval_dt %>% 
    .[spikein!=""] %>% 
    merge(
        param_grid,
        by="param_grid_ind"
    ) 


# Hyperparameter tuning - #ADEs per spikein --------------------------------------------

nades_in_dynamic_dt <- 
  eval_param_dt[,
                .(param_grid_ind,spikein,ade)
  ] %>% 
  unique() %>% 
  .[,.(total_dynamic_ades = .N),.(param_grid_ind,spikein)]

# Hyperparameter tuning - #ADEs per cluster --------------------------------------------

nades_in_cluster_dt <- 
  eval_param_dt[,
                .(param_grid_ind,cluster,ade)
  ] %>% 
  unique() %>% 
  .[,.(total_cluster_ades = .N),.(param_grid_ind,cluster)]

# Hyperparameter tuning - #ADEs per spikein,cluster --------------------------------------------


nades_in_cluster_dynamic_dt <- 
eval_param_dt[,
              .(param_grid_ind,cluster,spikein,ade)
              ] %>% 
    unique() %>% 
    .[,.N,.(param_grid_ind,cluster,spikein)] %>%
  merge(nades_in_dynamic_dt,
        by=c("param_grid_ind","spikein")
  ) %>% 
  merge(nades_in_cluster_dt,
        by=c("param_grid_ind","cluster")
  ) %>% 
  .[,.(param_grid_ind,spikein,cluster,N,
       frac_cluster = N/total_cluster_ades,
       frac_dynamic = N/total_dynamic_ades)]

nades_in_cluster_dynamic_dt[
    param_grid_ind==74
    ] %>% 
    .[spikein!=""] %>% 
    dcast(spikein ~ cluster,value.var="N")


# Example params: pg ------------------------------------------------------

pg=74

param_grid[
  param_grid_ind==74
]

# Example of cluster purity ----------------------------------------

g <- nades_in_cluster_dynamic_dt[
  param_grid_ind==pg
] %>% 
  ggplot(aes(factor(cluster),spikein,color=frac_cluster)) +
  geom_text(aes(label=round(frac_cluster,2))) +
  geom_point(size=15,pch=21) +
  scale_color_gradient(high="red",low="blue",name="Fraction of \ncanonical drug-events") +
  xlab("Cluster") +
  ylab("") +
  ggtitle(paste0("Cluster purity"))
ggsave(paste0(img_dir,basename,"example_of_cluster_purity.png"),g,width=7,height=3)

# Hyperparameter tuning - cluster purity --------------------------------------------

ch_dt <- 
    nades_in_cluster_dynamic_dt %>% 
    .[spikein!=""] %>% 
    .[,
      .(
          cluster_homog = Rfast::nth(frac_cluster,1,descending = T) / 
              Rfast::nth(frac_cluster,2,descending = T)
      ),
      .(param_grid_ind,cluster)
    ] %>% 
    .[,.(cluster_homog_sum = mean(cluster_homog)),param_grid_ind]

g <- ch_dt %>% 
    merge(param_grid) %>% 
  .[order(cluster_homog_sum,decreasing = T)] %>% 
    ggplot(aes(forcats::fct_reorder(factor(param_grid_ind),cluster_homog_sum,.desc = T),
               cluster_homog_sum,fill=factor(k))) +
    geom_bar(stat="identity") +
    scale_fill_brewer(palette = "Blues") +
    guides(fill=guide_legend(title="K")) +
    xlab("DTW model") +
    ylab("Cluster purity") +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    )
ggsave(paste0(img_dir,basename,"cluster_purity.png"),g,width=5,height=4)

# Example of dynamics localization ----------------------------------------

g <- nades_in_cluster_dynamic_dt[
  param_grid_ind==pg
] %>% 
  ggplot(aes(factor(cluster),spikein,color=frac_dynamic)) +
  geom_text(aes(label=round(frac_dynamic,2))) +
  geom_point(size=15,pch=21) +
  scale_color_gradient(high="red",low="blue",name="Fraction of \ncanonical drug-events") +
  xlab("Cluster") +
  ylab("") +
  ggtitle(paste0("Dynamics localization"))
ggsave(paste0(img_dir,basename,"example_of_dynamics_localization.png"),g,width=7,height=3)

# hyperparameter tuning - dynamics localization --------------------------------------------

sl_dt <- 
    nades_in_cluster_dynamic_dt %>% 
    .[,
      .(
          spikein_homog = Rfast::nth(frac_dynamic,1,descending = T) / 
              Rfast::nth(frac_dynamic,2,descending = T)
      ),
      .(param_grid_ind,spikein)
    ] %>% 
    .[,.(spikein_homog_sum = mean(spikein_homog)),param_grid_ind]

g <- sl_dt %>% 
    merge(param_grid) %>% 
  .[order(spikein_homog_sum,decreasing = T)] %>% 
    ggplot(aes(forcats::fct_reorder(factor(param_grid_ind),spikein_homog_sum,.desc = T),
               spikein_homog_sum,fill=factor(k))) +
    geom_bar(stat="identity") +
    scale_fill_brewer(palette = "Blues") +
    guides(fill=guide_legend(title="K")) +
    xlab("DTW model") +
    ylab("Dynamics localization") +
    theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
    )

ggsave(paste0(img_dir,basename,"dynamics_localization.png"),g,width=5,height=4)

# hyperparameter tuning - example of cluster purity vs dynamics localization --------------------------------------------

tmp <- 
  nades_in_cluster_dynamic_dt %>% 
  .[param_grid_ind==pg] %>% 
  merge(param_grid,by="param_grid_ind")

tmp %>% dcast(spikein ~ cluster,value.var = "frac_cluster")
tmp %>% dcast(spikein ~ cluster,value.var = "frac_dynamic")
tmp[,.(m = frac_cluster * frac_dynamic),.(spikein,cluster)] %>% dcast(spikein ~ cluster,value.var = "m")
g <- tmp[,.(m = frac_cluster * frac_dynamic),.(spikein,cluster)] %>% 
  ggplot(aes(cluster,spikein,fill=m)) +
  geom_tile() +
  xlab("") +
  ylab("") +
  theme(
    legend.position = "none"
  )
ggsave(paste0(img_dir,basename,"example_of_dynamic_localization_vs_cluster_purity_heatmap.png"),g,width=2,height=1.5)
g <- 
  tmp %>% 
  ggplot(aes(frac_cluster,frac_dynamic,color=factor(cluster),shape=spikein)) +
  geom_point(size=4) +
  guides(color=guide_legend(title="Cluster"),
         shape=guide_legend(title="Dynamic")) +
  xlab("Dynamics localization") +
  ylab("Cluster purity") 
ggsave(paste0(img_dir,basename,"example_of_dynamic_localization_vs_cluster_purity.png"),g,width=5,height=3)

# hyperparameter tuning - example aggregation of cluster purity vs dynamics localization --------------------------------------------

g <- 
  merge(
    nades_in_cluster_dynamic_dt %>% 
      .[param_grid_ind==pg] %>% 
      .[,.(max_frac_cluster = max(frac_cluster)),.(param_grid_ind,spikein)] %>% 
      .[,.(mean_max_frac_cluster = mean(max_frac_cluster)),.(param_grid_ind)],
    nades_in_cluster_dynamic_dt %>% 
      .[param_grid_ind==pg] %>% 
      .[,.(max_frac_dynamic = max(frac_dynamic)),.(param_grid_ind,cluster)] %>% 
      .[,.(mean_max_frac_dynamic = mean(max_frac_dynamic)),.(param_grid_ind)],
    by="param_grid_ind"
  ) %>% 
  merge(param_grid,by="param_grid_ind") %>% 
  .[,.(mean_max_frac_cluster,mean_max_frac_dynamic,
       param_grid_ind,label="optimized score")] %>% 
  ggplot(aes(mean_max_frac_cluster,mean_max_frac_dynamic)) +
  geom_point(size=4) +
  ggrepel::geom_text_repel(aes(label=label),force=80,nudge_x = 0.1,nudge_y=-0.1) +
  xlim(0,1) +
  ylim(0,1) +
  xlab("Dynamics localization") +
  ylab("Cluster purity") 
ggsave(paste0(img_dir,basename,"example_of_aggregate_dynamic_localization_vs_cluster_purity.png"),g,width=4,height=3)

# hyperparameter tuning - aggregation of cluster purity vs dynamics localization --------------------------------------------

clusters_for_purity <- 
  nades_in_cluster_dynamic_dt %>% 
  .[order(frac_cluster,decreasing = T)] %>% 
  .[,.SD[1],.(param_grid_ind,spikein)] %>% 
  .[,.(param_grid_ind,spikein,cluster)] %>% 
  unique()

tmp = merge(
  nades_in_cluster_dynamic_dt %>% 
    .[,.(max_frac_cluster = max(frac_cluster)),.(param_grid_ind,spikein)] %>% 
    .[,.(mean_max_frac_cluster = mean(max_frac_cluster)),.(param_grid_ind)],
  nades_in_cluster_dynamic_dt %>% 
    .[,.(max_frac_dynamic = max(frac_dynamic)),.(param_grid_ind,cluster)] %>% 
    merge(clusters_for_purity,by=c("param_grid_ind","cluster")) %>% 
    .[,.(mean_max_frac_dynamic = mean(max_frac_dynamic)),.(param_grid_ind)],
  by="param_grid_ind"
) %>% 
  merge(param_grid,by="param_grid_ind") %>% 
  merge(
    eval_dt[,.(param_grid_ind,time)] %>% 
      unique(),
    by="param_grid_ind"
  )

g <- 
  tmp %>% 
  .[k==3] %>% 
  ggplot(aes(
    mean_max_frac_cluster,
    mean_max_frac_dynamic,
    color=time,
    shape=centroid)) +
  geom_point(size=6) +
  scale_color_distiller(
    "Time in seconds",
    palette = "Spectral",
    breaks=c(1,100,200,300)
  ) +
  scale_shape_manual(values=15:18) +
  facet_grid(factor(k) ~ distance,scales="free") +
  xlab("Dynamics localization score") +
  ylab("Cluster purity score") 
ggsave(paste0(img_dir,basename,"aggregate_dynamic_localization_vs_cluster_purity_cluster3.png"),g,width=10,height=5)

g <- 
  tmp %>% 
  ggplot(aes(
    mean_max_frac_cluster,
    mean_max_frac_dynamic,
    color=time,
    shape=centroid)) +
  geom_point(size=3) +
  scale_color_distiller(
    "Time in seconds",
    palette = "Spectral",
    breaks=c(1,100,200,300)
  ) +
  scale_shape_manual(values=15:18) +
  facet_grid(distance ~ factor(k)) +
  ylim(0,1) +
  geom_abline(intercept=0.9,slope=0,col="red",linetype="dashed") +
  geom_vline(xintercept = 0.9,col="red",linetype="dashed") +
  xlab("Dynamics localization score") +
  ylab("Cluster purity score") 
ggsave(paste0(img_dir,basename,"aggregate_dynamic_localization_vs_cluster_purity.png"),g,width=15,height=6)

tmp2 <- 
  tmp[,.(mean_max_frac_cluster,mean_max_frac_dynamic,opt = mean_max_frac_cluster*mean_max_frac_dynamic,time,param_grid_ind)] %>% 
  merge(param_grid)
setorderv(tmp2,
          cols = c("opt","time"),order=c(-1,-1)
)
best_params <- 
  tmp2 

best_params %>% 
  fwrite(paste0(data_dir,basename,"hyperparameter_set_scores.csv"))
  

# pattern evaluation function and max metrics---------------------------------------------

#' Input: eval_dt
#' Output: cluster metric, dynamics metric, and cluster*dynamics metric
evaluate_patterns_in_parameter_sets <- function(dt,param_grid){

  # remove drug-events that are not canonical from the parameter results
  dt <- eval_dt %>% .[spikein!=""]
  
  # count number of canonical drug-events per dynamic
  # should be the same for each parameter set
  dynamic_dt <- 
    dt[,
            .(param_grid_ind,spikein,ade)
    ] %>% 
    unique() %>% 
    .[,.(total_dynamic_ades = .N),.(param_grid_ind,spikein)]
  
  # count number of canonical drug-events in each cluster
  cluster_dt <- 
    dt[,
            .(param_grid_ind,cluster,ade)
            ] %>% 
    unique() %>% 
    .[,.(total_cluster_ades = .N),.(param_grid_ind,cluster)]
  
  # count number of canonical drug-events per dynamic in clusters
  # calculate fraction across clusters and across dynamics
  cluster_dynamic_dt <- 
    dt[,
                  .(param_grid_ind,cluster,spikein,ade)
    ] %>% 
    unique() %>% 
    .[,.N,.(param_grid_ind,cluster,spikein)] %>%
    merge(dynamic_dt,
          by=c("param_grid_ind","spikein")
    ) %>% 
    merge(cluster_dt,
          by=c("param_grid_ind","cluster")
    ) %>% 
    .[,.(param_grid_ind,spikein,cluster,N,
         frac_cluster = N/total_cluster_ades,
         frac_dynamic = N/total_dynamic_ades)]
  
  clusters_for_purity <- 
    cluster_dynamic_dt %>% 
    .[order(frac_cluster,decreasing = T)] %>% 
    .[,.SD[1],.(param_grid_ind,spikein)] %>% 
    .[,.(param_grid_ind,spikein,cluster)] %>% 
    unique()
  
  tmp = merge(
    cluster_dynamic_dt %>% 
      .[,.(max_frac_cluster = max(frac_cluster)),.(param_grid_ind,spikein)] %>% 
      .[,.(mean_max_frac_cluster = mean(max_frac_cluster)),.(param_grid_ind)],
    cluster_dynamic_dt %>% 
      .[,.(max_frac_dynamic = max(frac_dynamic)),.(param_grid_ind,cluster)] %>% 
      merge(clusters_for_purity,by=c("param_grid_ind","cluster")) %>% 
      .[,.(mean_max_frac_dynamic = mean(max_frac_dynamic)),.(param_grid_ind)],
    by="param_grid_ind"
  ) %>% 
    merge(param_grid,by="param_grid_ind") %>% 
    merge(
      dt[,.(param_grid_ind,time)] %>% 
        unique(),
      by="param_grid_ind"
    )
  
  tmp %>% 
    .[,.(param_grid_ind,
         mean_max_frac_dynamic,
         mean_max_frac_cluster,
         opt = mean_max_frac_dynamic * mean_max_frac_cluster)] %>% 
    .[order(opt,decreasing = T)]
}

evaluate_patterns_in_parameter_sets(eval_dt,param_grid)


# str_pattern evaluaton ---------------------------------------------------

eval_dt_norm_ades <- 
eval_dt %>% 
  melt(id.vars=c("ade"),measure.vars=stages,
       variable.name="nichd",value.name="norm") %>% 
  unique()

eval_dt_norm_ades[,rank:=frank(norm),.(ade)]

eval_dt_rank <- 
  eval_dt_norm_ades %>% 
  .[order(nichd)] %>% 
  dcast(ade~nichd,value.var="rank")

eval_dt_rank[, str_pattern:=do.call(paste0,.SD), .SDcols=-1]

merge(
  eval_dt,
  eval_dt_rank[,.(ade,str_pattern)],
  by="ade"
) %>% 
  .[spikein!="",.(ade,spikein,cluster,str_pattern)] %>% 
  unique() %>% 
  .[,.(spikein,str_pattern)] %>% 
  unique()

# hyperparameter tuning - all patterns in clusters for best set -----------

for(i in 1:5){
  
  print(best_params[i])
  
  g <- norm_ades %>% 
    merge(
      eval_dt[
        param_grid_ind==best_params[i,param_grid_ind]
        ,.(ade,cluster)
      ],
      by="ade"
    ) %>% 
    ggplot(aes(factor(nichd,levels=stages),norm,color=factor(cluster),group=ade)) +
    geom_line(alpha=0.1) +
    guides(color=guide_legend(title="Cluster",override.aes = list(alpha = 1))) +
    xlab("") +
    ylab("Normalized score") +
    facet_grid(cluster~.) +
    theme(
      axis.text.x = element_text(angle=45,vjust=1,hjust=1)
    )
  ggsave(paste0(img_dir,basename,"drug_event_risk_patterns_in_clusters_",i,"best_hyperparameter_set.png"),g,width=6,height=6)
  
}

# Bootstrap hyperparameter tuned DTW from munnin -----------------------------------------------

boot_eval_dt <- 
  fread(paste0(data_dir,"database_generation_dtw_bootstrap_pattern_evaluation_tuning_results.csv"))


boot_param_grid_dt <- 
  fread(paste0(data_dir,"database_generation_dtw_bootstrap_pattern_evaluation_hyperparameter_grid.csv"))

inds <- boot_param_grid_dt[
  type=="partitional" & 
    distance %in% c("dtw_basic","gak","sbd") &
    centroid %in% c("mean","dba","pam","shape") &
    k %in% seq(3,6),
  param_grid_ind
  ]

boot_rank <- 
boot_eval_dt[param_grid_ind %in% inds,
             .(param_grid_ind,rank = frank(-opt)),
             .(bootstrap)
] %>% 
  .[,.(mrank = mean(rank),sdrank = sd(rank)),param_grid_ind] %>% 
  .[order(mrank)]

merge(
  boot_eval_dt %>% 
    .[,.(avg_time = mean(time),
         avg_cp = mean(mean_max_frac_cluster),
         avg_dl = mean(mean_max_frac_dynamic),
         avg_opt = mean(opt)),
      param_grid_ind],
  boot_param_grid_dt,
  by="param_grid_ind"
) %>% 
  merge(boot_rank) %>%
  .[order(avg_opt,decreasing = T)] %>% 
  fwrite(paste0(data_dir,basename,"bootstrap_hyperparameter_set_scores.csv"))

svals <- 15:20

cvals <- RColorBrewer::brewer.pal(4,"Dark2")

g <- merge(
  boot_rank,
  boot_param_grid_dt,
  by="param_grid_ind"
) %>% 
  ggplot(aes(sdrank,mrank,color=factor(k),
             shape=centroid)) +
  geom_point(size=5) +
  guides(color=guide_legend(title="Cluster")) +
  scale_color_manual(values=cvals) +
  scale_shape_manual(values=svals) +
  facet_wrap(~distance) +
  xlab("Rank standard error") +
  ylab("Average rank") +
  theme(
    strip.background = element_blank()
  )
ggsave(paste0(img_dir,basename,"bootstrap_aggregate_sd_vs_avg_rank_color_cluster.png"),g,width=8,height=4)

g <- merge(
  boot_eval_dt %>% 
    .[,.(avg_time = mean(time),
         avg_cp = mean(mean_max_frac_cluster),
         avg_dl = mean(mean_max_frac_dynamic),
         avg_opt = mean(opt)),
      param_grid_ind],
  boot_param_grid_dt,
  by="param_grid_ind"
) %>% 
  merge(boot_rank) %>% 
  .[order(avg_opt,decreasing = T)] %>% 
  ggplot(aes(avg_dl,avg_cp,color=log10(avg_time),
             shape=centroid)) +
  geom_point(size=5) +
  scale_color_distiller(
    "Time in seconds",
    palette = "PRGn",
    breaks=c(log10(3),1,2,3),
    labels=round(10^c(log10(3),1,2,3),1)
  ) +
  scale_shape_manual(values=svals) +
  facet_grid(distance ~ factor(k)) +
  geom_abline(intercept=0.8,slope=0,col="red",linetype="dashed") +
  geom_vline(xintercept = 0.8,col="red",linetype="dashed") +
  xlab("Dynamics localization score") +
  ylab("Cluster purity score") +
  theme(
    strip.background = element_blank()
  )
ggsave(paste0(img_dir,basename,"bootstrap_aggregate_dynamic_localization_vs_cluster_purity.png"),g,width=12,height=6)


g <- merge(
  boot_eval_dt %>% 
    .[,.(avg_time = mean(time),
         avg_cp = mean(mean_max_frac_cluster),
         avg_dl = mean(mean_max_frac_dynamic),
         avg_opt = mean(opt)),
      param_grid_ind],
  boot_param_grid_dt,
  by="param_grid_ind"
) %>% 
  merge(boot_rank) %>% 
  .[order(avg_opt,decreasing = T)] %>% 
  ggplot(aes(avg_opt,mrank,color=log10(avg_time),
             shape=centroid)) +
  geom_point(size=5) +
  scale_color_distiller(
    "Time in seconds",
    palette = "PRGn",
    breaks=c(log10(3),1,2,3),
    labels=round(10^c(log10(3),1,2,3),1)
  ) +
  scale_shape_manual(values=svals) +
  facet_grid(distance ~ factor(k)) +
  xlab("Average clustering optimization score") +
  ylab("Average rank") 
ggsave(paste0(img_dir,basename,"bootstrap_aggregate_clustering_opt_score_by_rank.png"),g,width=12,height=6)


g <- merge(
  boot_eval_dt %>% 
    .[,.(avg_time = mean(time),
         avg_cp = mean(mean_max_frac_cluster),
         avg_dl = mean(mean_max_frac_dynamic),
         avg_opt = mean(opt)),
      param_grid_ind],
  boot_param_grid_dt,
  by="param_grid_ind"
) %>% 
  merge(boot_rank) %>% 
  .[order(avg_opt,decreasing = T)] %>% 
  ggplot(aes(avg_time,avg_opt,color=factor(k),
             shape=centroid)) +
  geom_point(size=5) +
  scale_x_continuous(trans="log10",labels=scales::comma) +
  guides(color=guide_legend(title="Cluster")) +
  scale_color_manual(values=cvals) +
  scale_shape_manual(values=svals) +
  facet_wrap(~distance) +
  xlab("Average time in seconds") +
  ylab("Clustering optimization score")
ggsave(paste0(img_dir,basename,"bootstrap_aggregate_clustering_opt_score_by_time_color_cluster.png"),g,
       width=8,height=4)

g <- merge(
  boot_eval_dt %>% 
    .[,.(avg_time = mean(time),
         avg_cp = mean(mean_max_frac_cluster),
         avg_dl = mean(mean_max_frac_dynamic),
         avg_opt = mean(opt)),
      param_grid_ind],
  boot_param_grid_dt,
  by="param_grid_ind"
) %>% 
  merge(boot_rank) %>% 
  .[order(avg_opt,decreasing = T)] %>% 
  ggplot(aes(avg_dl,avg_cp,color=factor(k),
             shape=centroid)) +
  geom_point(size=5) +
  guides(color=guide_legend(title="Cluster")) +
  scale_color_manual(values=cvals) +
  scale_shape_manual(values=svals) +
  facet_wrap(~distance) +
  xlab("Dynamics localization score") +
  ylab("Cluster purity score") +
  theme(
    strip.background = element_blank()
  )
ggsave(paste0(img_dir,basename,"bootstrap_aggregate_clustering_dynamic_localization_vs_cluster_purity_color_cluster.png"),g,
       width=8,height=4)

g <- merge(
  boot_eval_dt %>% 
    .[,.(avg_time = mean(time),
         avg_cp = mean(mean_max_frac_cluster),
         avg_dl = mean(mean_max_frac_dynamic),
         avg_opt = mean(opt)),
      param_grid_ind],
  boot_param_grid_dt,
  by="param_grid_ind"
) %>% 
  merge(boot_rank) %>% 
  .[order(avg_opt,decreasing = T)] %>% 
  .[,.(
    dl_lwr = quantile(avg_dl,c(0.025)),
    dl_mean = mean(avg_dl),
    dl_upr = quantile(avg_dl,c(0.975)),
    cp_lwr = quantile(avg_cp,c(0.025)),
    cp_mean = mean(avg_cp),
    cp_upr = quantile(avg_cp,c(0.975))
    ),.(k)] %>% 
  ggplot(aes(x=k)) +
  geom_point(aes(y=dl_mean,color="Dynamics localization")) +
  geom_errorbar(aes(ymin=dl_lwr,ymax=dl_upr,color="Dynamics localization"),width=0.2) +
  geom_point(aes(y=cp_mean,color="Cluster Purity"),
             position = position_nudge(x=0.1)) +
  geom_errorbar(aes(ymin=cp_lwr,ymax=cp_upr,color="Cluster Purity"),width=0.2,
                position = position_nudge(x=0.1)) +
  scale_color_manual(values = c("Dynamics localization" = "red","Cluster Purity" = "black")) +
  guides(color=guide_legend(title="Metric")) +
  ylab("Metric") +
  xlab("K") +
  theme(
    legend.position = "bottom"
  )
ggsave(paste0(img_dir,basename,"bootstrap_aggregate_clustering_metric_by_K.png"),g,
       width=6,height=4)

g <- merge(
  boot_eval_dt %>% 
    .[,.(avg_time = mean(time),
         avg_cp = mean(mean_max_frac_cluster),
         avg_dl = mean(mean_max_frac_dynamic),
         avg_opt = mean(opt)),
      param_grid_ind],
  boot_param_grid_dt,
  by="param_grid_ind"
) %>% 
  merge(boot_rank) %>% 
  .[order(avg_opt,decreasing = T)] %>% 
  .[,.(
    dl_lwr = quantile(avg_dl,c(0.025)),
    dl_mean = mean(avg_dl),
    dl_upr = quantile(avg_dl,c(0.975)),
    cp_lwr = quantile(avg_cp,c(0.025)),
    cp_mean = mean(avg_cp),
    cp_upr = quantile(avg_cp,c(0.975))
  ),.(distance)] %>% 
  ggplot(aes(x=distance)) +
  geom_point(aes(y=dl_mean,color="Dynamics localization")) +
  geom_errorbar(aes(ymin=dl_lwr,ymax=dl_upr,color="Dynamics localization"),width=0.2) +
  geom_point(aes(y=cp_mean,color="Cluster Purity"),
             position = position_nudge(x=0.1)) +
  geom_errorbar(aes(ymin=cp_lwr,ymax=cp_upr,color="Cluster Purity"),width=0.2,
                position = position_nudge(x=0.1)) +
  scale_color_manual(values = c("Dynamics localization" = "red","Cluster Purity" = "black")) +
  guides(color=guide_legend(title="Metric")) +
  ylab("Metric") +
  xlab("Distance") +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1),
    legend.position = "bottom"
  )
ggsave(paste0(img_dir,basename,"bootstrap_aggregate_clustering_metric_by_distance.png"),g,
       width=6,height=4)
