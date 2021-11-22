#' Title: "A database of pediatric drug effects to evaluate ontogenic mechanisms from child growth and development" study
#' 
#' Script author details: Nicholas Giangreco
#' 
#' This script generates the data and plots to evaluate the bias and confounding mitigation by our dGAMs

# PURPOSE -----------------------------------------------------------------

#' To evaluate the correlation between drug-event GAM scores and 
#' covariate prevalence for pediatric-specific events. 
#' 

# Setup -------------------------------------------------------------------

pacman::p_load(tidyverse,data.table,doParallel)

data_dir <- "../../data/"
img_dir <- "../../docs/imgs/"
basename <- "database_generation_pediatric_ae_"
seed = 0
set.seed(seed)

cores=4
registerDoParallel(cores=cores)

stages <- 
  c("term_neonatal","infancy",
    "toddler","early_childhood",
    "middle_childhood","early_adolescence",
    "late_adolescence")

theme_set(theme_bw(base_size=16) + theme(text = element_text(face="bold")))

# load meddra pediatric aes -----------------------------------------------

ped_aes <- 
  fread(paste0(data_dir,"paediatric_term_list_19-0_concept_joined.csv"))


# PROCESSING -- load data -------------------------------------------------

source("database_generation_load_data.R")

# PROCESSING -- Load GAM data -----------------------------------------------------------

source("database_generation_load_GAM_data.R")

# PROCESSING -- calculate AE report counts per stage ----------------------

stage_ae_counts <- 
  raw_data[,
           .(safetyreportid,meddra_concept_id)
  ] %>% 
  unique() %>% 
  .[,.N,.(meddra_concept_id)]

ae_counts <- 
  raw_data[,
         .(safetyreportid,nichd,meddra_concept_id)
         ] %>% 
  unique() %>% 
  .[,.N,.(meddra_concept_id,nichd)]

ped_ae_counts <- 
ae_counts[
  meddra_concept_id %in% ped_aes$concept_id
  ] %>% 
  dcast(
    meddra_concept_id ~ factor(nichd,levels=stages),
    value.var="N",
    fill = 0
    )

stage_only_logical <- 
  rowSums(ped_ae_counts[,stages,with=F])==apply(ped_ae_counts[,stages,with=F],1,max)

stage_only_aes <- 
  ped_ae_counts[stage_only_logical]$meddra_concept_id


# PROCESSING -- calculate drug reports probability per stage --------------

drug_counts <- 
  raw_data[,
           .(safetyreportid,nichd,atc_concept_id)
  ] %>% 
  unique() %>% 
  .[,.(total_drug_stage_N = .N),.(atc_concept_id,nichd)] %>% 
  merge(
    raw_data[,
             .(safetyreportid,atc_concept_id)
    ] %>% 
      unique() %>% 
      .[,.(total_drug_N = .N),.(atc_concept_id)]
  )


# PROCESSING -- perform stage-chi squared tests ---------------------------

chi_dt_new <- 
  lapply(
    1:nrow(ped_ae_counts),
    function(i){
      tab <- 
        ped_ae_counts[i,stages,with=F] %>% unlist %>% unname
      mat <- 
        matrix(
          c(tab,
            rep(
              stage_ae_counts[meddra_concept_id==ped_ae_counts[i,meddra_concept_id],N],
                length(tab)
              )
            ),
          nrow=length(tab))
      chi <- chisq.test(mat)
      data.table(
        meddra_concept_id = ped_ae_counts[i,meddra_concept_id],
        chi_statistic = chi$statistic %>% unname,
        chi_pvalue = chi$p.value %>% unname,
        nichd_max_reporting = stages[which.max(tab)],
        max_reporting_or = 
          (tab[which.max(tab)]/stage_ae_counts[meddra_concept_id==ped_ae_counts[i,meddra_concept_id],N])/((sum(tab)-tab[which.max(tab)])/stage_ae_counts[meddra_concept_id==ped_ae_counts[i,meddra_concept_id],N]),
        nichd_max_reporting_N = tab[which.max(tab)],
        nichd_max_reporting_percentage = tab[which.max(tab)]/sum(tab)
      )
    }) %>% bind_rows()

chi_dt_new$one_stage_only <- chi_dt_new$meddra_concept_id %in% stage_only_aes
chi_dt_new$fdr <- p.adjust(chi_dt_new$chi_pvalue,method = "fdr")


# PROCESSING -- join chi and risk data ------------------------------------

plot_data <- 
dts[,
    .(database,meddra_concept_id,atc_concept_id,nichd,D,E,DE,gam_score)
    ] %>% 
  unique() %>% 
  merge(
    chi_dt_new,
    by="meddra_concept_id",
    allow.cartesian = T
  ) %>%
  merge(
    drug_counts,
    by=c("atc_concept_id","nichd")
  ) %>% 
  merge(
    dts[database=='base' &
      meddra_concept_id %in% chi_dt$meddra_concept_id,
      .SD[which.max(gam_score)],
      .(meddra_concept_id,atc_concept_id)
    ] %>% 
      .[,
        .(nichd_max_base_risk = nichd,meddra_concept_id,atc_concept_id)
      ],
    by=c("meddra_concept_id","atc_concept_id")
  ) %>% 
  merge(
    dts[database=='base' &
      meddra_concept_id %in% chi_dt$meddra_concept_id,
      .SD[which.max(E)],
      .(meddra_concept_id,atc_concept_id)
    ] %>% 
      .[,
        .(nichd_max_base_reporting = nichd,meddra_concept_id,atc_concept_id)
      ],
    by=c("meddra_concept_id","atc_concept_id")
  )

plot_data$percent_drug <- 
  plot_data$total_drug_stage_N/plot_data$total_drug_N

# Visualization --------

fdr_=0.05
tmp <- plot_data %>%  
  .[fdr<fdr_ &
      nichd_max_base_reporting==nichd &
      nichd_max_base_risk==nichd
  ] %>% unique()
tmp <- tmp[order(percent_drug)]
tmp$bin <- 0
seq_=floor(seq(1,nrow(tmp),length=100))
for(i in seq_along(seq_)){
  if(i==100){
    next
  }else{
    inds=seq(seq_[i],seq_[i+1],1)
    tmp[inds,bin:= i]
  }
}


ca_lm <- lm(gam_score ~ percent_drug,data=tmp[database=="covariate_adjusted"])
b_lm <- lm(gam_score ~ percent_drug,data=tmp[database=="base"])

g <- tmp[,
    .(
      m = mean(gam_score),
      v = sd(gam_score),
      avg = mean(percent_drug)
      ),
    .(database,bin)
    ] %>% 
  ggplot(aes(avg,m,color=database)) +
  geom_point() +
  geom_errorbar(aes(ymin=m-v,ymax=m+v),width=0) +
  geom_smooth(data=tmp,aes(percent_drug,gam_score),method="lm",color="darkgrey") +
  scale_color_brewer(palette="Set2") +
  scale_x_continuous(
    breaks=c(
      seq(0.1,1,length=10)
      ),
    minor_breaks = tmp[,.(m = mean(percent_drug)),.(bin)]$m,
    labels=scales::label_percent(accuracy=1)
    ) +
  xlab("Drug reporting percent in stage") +
  ylab('Log odds risk across\npediatric adverse events') +
  facet_wrap(~database) +
  theme(
    legend.position = "none",
    strip.background = element_blank()
  ) +
  ggtitle(
    paste0(
      "FDR < ",
      fdr_,
      "\nN_events = ",
      tmp[,length(unique(meddra_concept_id))],
      "\nBase int: ",round(coef(b_lm)[1],2),
      " [",round(confint(b_lm)[1,1],2),", ",round(confint(b_lm)[1,2],2),"]",
      " CA int: ",round(coef(ca_lm)[1],2),
      " [",round(confint(ca_lm)[1,1],2),", ",round(confint(ca_lm)[1,2],2),"]",
      "\nBase slope: ",round(coef(b_lm)[2],2),
      " [",round(confint(b_lm)[2,1],2),", ",round(confint(b_lm)[2,2],2),"]",
      " CA slope: ",round(coef(ca_lm)[2],2),
      " [",round(confint(ca_lm)[2,1],2),", ",round(confint(ca_lm)[2,2],2),"]"
      )
    )

ggsave(paste0(img_dir,basename,"base_vs_covariated_adjusted_correction.png"),g,width=12,height=5)
