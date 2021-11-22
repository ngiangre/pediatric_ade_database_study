#' Title: "A database of pediatric drug effects to evaluate ontogenic mechanisms from child growth and development" study
#' 
#' Script author details: Nicholas Giangreco
#' 
#' This script generates the plots of statistic distributions
#' between GAM fit and/or generalization for the different
#' GAM specifications per a random drug-event sample

# Purpose -----------------------------------------------------------------

#' To compare the statistical fit with the ML generalization of the GAMs
#' 

# Setup -------------------------------------------------------------------

pacman::p_load(data.table,tidyverse,doParallel)

data_dir <- "../../data/"
img_dir <- "../../docs/imgs/"

seed = 0
set.seed(seed)

basename <- "database_generation_generalization_fit_viz_"
gout_dir <- paste0(data_dir,"database_generation_generalization/")
fout_dir <- paste0(data_dir,"database_generation_fit/")

method="fREML"

theme_set(theme_bw(base_size=16) + theme(text = element_text(face="bold")))

# Load data and covariates ------------------------------------------------

source("database_generation_load_data.R")

stage_knots <- 
    seq_along(category_levels[[stage_col]])

# Functions ---------------------------------------------------------------

source("database_generation_functions.R")


# Model colors ------------------------------------------------------------

model_colors <- rainbow(length(names(fcn_lst)))
names(model_colors) <- names(fcn_lst)
model_colors_no_base <- 
    model_colors[2:length(model_colors)]

# load generalization and fit data ---------------------------------------------------------------


files <- list.files(gout_dir)
gdt <- lapply(files,
             function(x){
                 fread(paste0(gout_dir,x))
             }) %>% 
    bind_rows()

files <- list.files(fout_dir)
fdt <- lapply(files,
              function(x){
                  fread(paste0(fout_dir,x))
              }) %>% 
    bind_rows()




# ade count ---------------------------------------------------------------

ade_count <- 
  raw_data[,.(safetyreportid,ade)] %>% 
  unique() %>% 
  .[,.N,ade]

# Getting base aic ---------------------------------------------------

fdt_base <- 
  fread(paste0(data_dir,"database_generation_base.csv")) %>% 
  .[ade %in% fdt[,unique(ade)]] %>% 
  .[,.(ade,model_aic = AIC,model_time = time)] %>% 
  unique()

fdt_base$model <- "Base"

# model stats datatable ---------------------------------------------------

stats_dt <- data.table()

# statistical fit - time vs aic ----------------------

stats_dt <- 
  fdt[,
      .(ade,model,model_aic,model_time)
  ] %>% 
  unique() %>% 
  bind_rows(
    fdt_base
  ) %>% 
  .[,.(
    aic = mean(model_aic) %>% round(3),
    time = mean(model_time) %>% round(3)
  ),model]

g <- fdt[,
         .(ade,model,model_aic,model_time)
] %>% 
  unique() %>% 
  bind_rows(
    fdt_base
  ) %>% 
  .[,.(
    aic = mean(model_aic),
    time = mean(model_time)
  ),model] %>% 
  ggplot(aes(time,aic,fill=model)) +
  geom_point(color="black",size=4,pch=21) +
  scale_fill_manual(values=model_colors) +
  guides(fill=guide_legend(title="GAM",title.position = "top",ncol=1)) +
  xlab("Time in seconds") +
  ylab("AIC") +
  theme(
    legend.position = "none"
  )
ggsave(paste0(img_dir,basename,"aic_by_time.png"),
       g,width=3,height=3)
g <- g +
  theme(
    legend.position = "bottom"
  )
ggsave(paste0(img_dir,basename,"aic_by_time_with_legend.png"),
       g,width=5,height=6)


# statistical fit - low number ADE reports vs AIC  -------------------------------

tmp <- 
  merge(
    fdt[,.(ade,model,model_aic)] %>% 
      bind_rows(
        fdt_base[,.(ade,model,model_aic)]
      ),
    ade_count,
    by="ade"
  ) %>% 
  unique()

tmp$N_cut <- cut(tmp$N,breaks=c(1,5,10,50))
tmp[is.na(N_cut),"N_cut"] <- paste0("(50,",max(tmp$N),"]")

g <- tmp[,
         .(
           lwr = quantile(model_aic,c(0.025)),
           m = mean(model_aic),
           upr = quantile(model_aic,c(0.975))
         ),
         .(model,N_cut)
] %>% 
  ggplot(aes(N_cut,m,color=forcats::fct_reorder(model,m,.fun=mean,.desc=F))) + 
  geom_point(position = position_dodge(width=0.6),size=0.5) +
  geom_errorbar(aes(ymin=lwr,ymax=upr),
                position = position_dodge(width=0.6)) +
  scale_color_manual(values=model_colors) +
  guides(color=guide_legend(title="Number of model coefficients")) +
  scale_y_continuous(trans="log10",labels=scales::comma) +
  xlab("Number of\ndrug-event reports") +
  ylab("AIC") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"model_aic_by_nreports.png"),g,width=3,height=3)

# statistical fit - time vs case probs ----------------------

stats_dt <- 
  merge(
    stats_dt,
    gdt[,
        .(ade,model,model_mean_true_probability,model_mean_false_probability,model_time)
    ] %>% 
      unique() %>% 
      .[,.(
        drugevent_probability_fold_change = (mean(model_mean_true_probability) /
          gdt[model=="Base",mean(model_mean_true_probability)]) %>% round(3)
      ),model],
    by="model"
  )

g <- gdt[,
    .(ade,model,model_mean_true_probability,model_mean_false_probability,model_time)
    ] %>% 
  unique() %>% 
    .[,.(
        control_prob = mean(model_mean_false_probability) /
          gdt[model=="Base",mean(model_mean_false_probability)],
        case_prob = mean(model_mean_true_probability) /
          gdt[model=="Base",mean(model_mean_true_probability)],
        time = mean(model_time)
    ),model] %>% 
  .[order(case_prob)] %>% 
    ggplot(aes(time,case_prob,fill=model)) +
    geom_point(color="black",size=4,pch=21) +
  geom_point(aes(time,control_prob),fill="gray",color="black",size=1,pch=22) +
  scale_fill_manual(values=model_colors) +
    xlab("Time in seconds") +
    ylab("Drug-event probability\nfold change") +
    theme(
        legend.position = "none"
    )
ggsave(paste0(img_dir,basename,"drug_event_probability_by_time.png"),
       g,width=3,height=3)
  
# generalization - test vs train auroc -----------------------------------

stats_dt <- 
  merge(
    stats_dt,
    gdt %>% 
      .[,.(ade,model,train_auc,test_auc)] %>% 
      unique() %>% 
      .[,.(mean_training_auroc = mean(train_auc) %>% round(3),
           mean_testing_auroc = mean(test_auc) %>% round(3)),
        .(model)],
    by="model"
  )

g <- 
  gdt %>% 
  .[,.(ade,model,train_auc,test_auc)] %>% 
  unique() %>% 
  .[,.(mean_train = mean(train_auc),
       mean_test = mean(test_auc)),
    .(model)] %>% 
  ggplot(aes(mean_train,mean_test,fill=model)) +
  geom_point(color="black",size=4,pch=21) +
  scale_fill_manual(values=model_colors) +
  geom_abline(intercept=0,slope=1,col="red",linetype=2) +
  xlab("Train AUROC") +
  ylab("Test AUROC") +
  theme(
    legend.position = "none"
  )
ggsave(paste0(img_dir,basename,"train_by_test_auroc_model_performance.png"),g,width=3,height=3)

# generalization - train vs test auroc dumbells ------------------------------------------

stats_dt <- 
  merge(
    stats_dt,
    gdt %>% 
      .[,.(ade,model,train_auc,test_auc)] %>% 
      unique() %>% 
      .[,.(mean_train = mean(train_auc),
           mean_test = mean(test_auc)),
        .(model)] %>% 
      .[,.(model,difference_from_training_auroc = (mean_test - mean_train) %>% round(3))],
    by="model"
    )

g <- gdt %>% 
  .[,.(ade,model,train_auc,test_auc)] %>% 
  unique() %>% 
  .[,.(mean_train = mean(train_auc),
       mean_test = mean(test_auc)),
    .(model)] %>% 
  .[,.(model,overfit = mean_test - mean_train)] %>% 
  .[order(overfit)] %>% 
  ggplot(aes(overfit,forcats::fct_inorder(model))) +
  geom_segment( aes(y=forcats::fct_inorder(model), 
                    yend=forcats::fct_inorder(model), 
                    x=0, xend=overfit)) +
  geom_point(size=5, aes(fill=model),color="black",pch=21) +
  scale_fill_manual(values=model_colors) +
  xlab("Difference from\ntraining AUROC") +
  ylab("") +
  theme(
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
ggsave(paste0(img_dir,basename,"overfit_model_performance.png"),g,width=3,height=3)

# generalization - test auroc by time -------------------------------------

g <- gdt %>% 
    .[,.(ade,test_auc,train_auc,model_time,model)] %>% 
    unique() %>% 
    .[,.(
        mean_auc = mean(test_auc),
        mean_time = mean(model_time)
    ),model] %>% 
    ggplot(aes(mean_time,mean_auc,fill=model)) +
    geom_point(pch=21,color="black",size=4) +
  scale_fill_manual(values=model_colors) +
    guides(fill=guide_legend(title.position = "top",ncol=2)) +
    scale_x_continuous(trans="log10",labels=scales::comma) +
    xlab("Time in seconds") +
    ylab("Test AUROC") +
    theme(
        legend.position = "none"
    )
ggsave(paste0(img_dir,basename,"average_test_performance_by_average_time.png"),g,width=3,height=3)



# output - model statistics -----------------------------------------------

stats_dt %>% 
  fwrite(paste0(data_dir,basename,"dGAM_statistics.csv"))
