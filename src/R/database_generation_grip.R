#' Title: "A database of pediatric drug effects to evaluate ontogenic mechanisms from child growth and development" study
#' 
#' Script author details: Nicholas Giangreco
#' 
#' This script generates the enrichment data and plots of our
#' identified significant drug-events in an pediatric ADE reference set

# Purpose -----------------------------------------------------------------

#' To visualize and compare risks within the GRiP ADE reference set
#' 

# Setup -------------------------------------------------------------------

pacman::p_load(tidyverse,data.table,doParallel)

data_dir <- "../../data/"
img_dir <- "../../docs/imgs/"
basename <- "database_generation_grip_"
seed = 0
set.seed(seed)

registerDoParallel(cores=5)

theme_set(theme_bw(base_size=16) + theme(text = element_text(face="bold")))

# Load data and covariates ------------------------------------------------

source("database_generation_load_data.R")

# Load GAM data ------------------------------------------------

source("database_generation_load_GAM_data.R")


# load GRiP reference set -------------------------------------------------

pediatric_reference <- 
    fread("https://raw.githubusercontent.com/ngiangre/GRiP_pediatric_ADE-reference_set/master/data/GRiP_Pediatric_ADE_Gold_Standard_List_minimal_joined.csv")

pediatric_reference$ade <- 
    pediatric_reference[,paste0(ATC_concept_id,"_",MedDRA_concept_id)]

pediatric_reference$ade_name <- 
    pediatric_reference[,paste0(ATC_concept_name," and ",MedDRA_concept_name)]

pediatric_reference_with_negatives <- 
    fread("https://raw.githubusercontent.com/ngiangre/GRiP_pediatric_ADE-reference_set/master/data/GRiP_Pediatric_ADE_with_generated_negatives.csv")

pediatric_reference_with_negatives$ade <- 
    pediatric_reference_with_negatives[,paste0(ATC_concept_id,"_",MedDRA_concept_id)]

pediatric_reference_with_negatives$ade_name <- 
    pediatric_reference_with_negatives[,paste0(ATC_concept_name," and ",MedDRA_concept_name)]

pediatric_reference_with_negatives[,table(Control)]

ades_child_both <- 
    pediatric_reference[Population %in% c("B","C"),unique(ade)]


# dGAM performance metrics to detect GRiP ADEs ---------------------------------

tmp <- 
    dts[
        database=="covariate_adjusted",
        .(ade,nichd,gam_score_90mse,gam_score)
    ] %>% 
    unique() %>% 
    merge(
        pediatric_reference_with_negatives[,
                                           .(ade,Control)
        ] %>% unique(),
        by="ade"
    ) 

nomsigany_avg <- 
    tmp[ade %in% tmp[gam_score_90mse>0,unique(ade)],
        .(Control,
          gam_score = mean(gam_score),
          gam_score_90mse = mean(gam_score_90mse)),
        ade] %>% unique()



metric_boot <- 
    function(dat,metric="auc",boot=100){
        vec <- sapply(1:boot,function(i){
            set.seed(i)
            ind <- sample(1:nrow(dat),nrow(dat),replace=T)
            truth <- 
                dat[ind,Control] %>% factor(levels=c("N","P"))
            estimate <- 
                dat[ind,gam_score_90mse]
            ROCR::performance(
                ROCR::prediction(estimate,truth),
                metric
            )@y.values[[1]]
        })
        return(
            c(
                "lwr" = quantile(vec,c(0.025)) %>% unname,
                "mean" = mean(vec),
                "upr" = quantile(vec,c(0.975)) %>% unname
            )
        )
    }

metric_boot(nomsigany_avg,"auc")
metric_boot(nomsigany_avg,"aucpr")

metric_boot(nomsigany_avg,"rec")
metric_boot(nomsigany_avg,"aucpr")

# PRR performance metrics to detect GRiP ADEs ---------------------------------

tmp <- 
    dts[
        database=="covariate_adjusted",
        .(ade,nichd,PRR_90mse,PRR)
    ] %>% 
    unique() %>% 
    merge(
        pediatric_reference_with_negatives[,
                                           .(ade,Control)
        ] %>% unique(),
        by="ade"
    )

nomsigany_avg <- 
    tmp[ade %in% tmp[PRR_90mse>0,unique(ade)],
        .(Control,
          PRR = mean(PRR,na.rm=T),
          PRR_90mse = mean(PRR_90mse,na.rm=T)),
        ade] %>% unique()



metric_boot <- 
    function(dat,metric="auc",boot=100){
        vec <- sapply(1:boot,function(i){
            set.seed(i)
            ind <- sample(1:nrow(dat),nrow(dat),replace=T)
            truth <- 
                dat[ind,Control] %>% factor(levels=c("N","P"))
            estimate <- 
                dat[ind,PRR_90mse]
            ROCR::performance(
                ROCR::prediction(estimate,truth),
                metric
            )@y.values[[1]]
        })
        return(
            c(
                "lwr" = quantile(vec,c(0.025)) %>% unname,
                "mean" = mean(vec),
                "upr" = quantile(vec,c(0.975)) %>% unname
            )
        )
    }

metric_boot(nomsigany_avg,"auc")
metric_boot(nomsigany_avg,"aucpr")

# dGAM statistically significant drug event risks vs total positive vs negative ADEs ----------------------------------


tmp <- 
    dts[
        database=="covariate_adjusted",
        .(ade,nichd,gam_score_90mse)
    ] %>% 
    unique() %>% 
    merge(
        pediatric_reference_with_negatives[,
                                           .(ade,Control)
        ] %>% unique(),
        by="ade"
    )

tmp[,length(unique(ade))]
tmp[Control=="P",length(unique(ade))]
tmp[Control=="N",length(unique(ade))]

merge(
    tmp[,.(total = length(unique(ade))),.(Control)], 
    tmp[gam_score_90mse>0,
    .(sig = length(unique(ade))),
    .(Control)]
)

a = 23
b = 75
c = 18
d = 112
or = (a / b) / (c / d)
se = sqrt( (1/a) + (1/b) + (1/c) + (1/d))

or
exp( log(or) - 1.96*se )
exp( log(or) + 1.96*se )

tp = a
fp = c
fn = (b-a)
tn = (d-c)
prec = tp/(tp+fp)
prec
rec = tp/b
rec
f1 = tp/(tp+0.5*(fp+fn))
f1

# PRR statistically significant drug event risks vs total positive vs negative ADEs ----------------------------------


tmp <- 
    dts[
        database=="covariate_adjusted",
        .(ade,nichd,PRR_90mse)
    ] %>% 
    unique() %>% 
    merge(
        pediatric_reference_with_negatives[,
                                           .(ade,Control)
        ] %>% unique(),
        by="ade"
    )

tmp[,length(unique(ade))]
tmp[Control=="P",length(unique(ade))]
tmp[Control=="N",length(unique(ade))]

merge(
    tmp[,.(total = length(unique(ade))),.(Control)], 
    tmp[PRR_90mse>0,
        .(sig = length(unique(ade))),
        .(Control)]
)

a = 71
b = 75
c = 110
d = 112
tp = a
fp = c
fn = (b-a)
tn = (d-c)
or = (a/b)/(c/d)
se = sqrt( (1/a) + (1/b) + (1/c) + (1/d))

or
exp( log(or) - 1.645*se )
exp( log(or) + 1.645*se )

tp = a
fp = c
fn = (b-a)
tn = (d-c)
prec = tp/(tp+fp)
prec
rec = tp/b
rec
f1 = tp/(tp+0.5*(fp+fn))
f1

# significant drug-event risks and reporting -----------------------------------

length(ades_child_both)

pediatric_reference[
    ade %in% sig_ades & 
        ade %in% ades_child_both &
        Control=="P",
    .(ade,ade_name)
] %>% 
    unique() %>% 
    merge(dts[database=="covariate_adjusted"][,.(ade,nichd,Drug = D,`Drug-event`=DE,Event = E,gam_score,gam_score_90mse,gam_score_90pse)],by=c("ade")) %>% 
    pivot_longer(cols=c("Drug","Event","Drug-event")) %>%
    mutate(
        value = log10(value + 1),
        nichd = factor(nichd,levels=category_levels$nichd)
    ) %>% select(ade) %>% distinct() %>% nrow()

22/75

g <- pediatric_reference[
    ade %in% sig_ades & 
        ade %in% ades_child_both &
        Control=="P",
    .(ade,ade_name)
] %>% 
    unique() %>% 
    merge(dts[database=="covariate_adjusted"][,.(ade,nichd,Drug = D,`Drug-event`=DE,Event = E,gam_score,gam_score_90mse,gam_score_90pse)],by=c("ade")) %>% 
    pivot_longer(cols=c("Drug","Event","Drug-event")) %>%
    mutate(
        value = log10(value + 1),
        nichd = factor(nichd,levels=category_levels$nichd)
    ) %>%
    ggplot(aes(factor(nichd,levels=category_levels$nichd))) +
    geom_bar(aes(nichd,value,fill=name),color="black",
             stat="identity",position="dodge") +
    geom_point(aes(y=gam_score),size=1,color="black",fill="gray") +
    geom_errorbar(aes(ymin=gam_score_90mse,ymax=gam_score_90pse),
                  color="black",width=0.2,size=0.5) +
    scale_y_continuous(sec.axis=sec_axis(~.,name="Risk of ADE\n(dGAM score)")) +
    guides(fill=guide_legend(title="Report type")) +
    xlab("") +
    ylab("log10(Number of reports)") +
    geom_abline(intercept = 0,slope=0,color="red",linetype="dashed") +
    facet_wrap(~ade_name,ncol=5,labeller = label_wrap_gen(width=20),scales="free_y") +
    theme(
        legend.position = "left",
        strip.background = element_blank(),
        axis.text.x = element_text(angle=45,vjust=1,hjust=1)
    )
ggsave(paste0(img_dir,basename,"significant_drug_events.png"),g,width=18,height=12)





