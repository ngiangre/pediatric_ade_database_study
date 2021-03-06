
#' Title: "A database of pediatric drug effects to evaluate ontogenic mechanisms from child growth and development" study
#' 
#' Script author details: Nicholas Giangreco
#' 
#' This script generates the plots of covariate distributions for
#' drug-events in Pediatric FAERS

# Purpose -----------------------------------------------------------------


#' Generate covariate distribution plots
#' 


# load libraries and set variables ----------------------------------------------------------

pacman::p_load(data.table,tidyverse)

data_dir <- "../../data/"
img_dir <- "../../docs/imgs/"
basename <- "database_generation_openFDA_summary_"

theme_set(theme_bw(base_size=16) + theme(text = element_text(face="bold")))

`%notin%` <- Negate(`%in%`)

# Load data and covariates ------------------------------------------------

source("database_generation_load_data.R")


# polypharmacy 95% CI --------------------------------------------------

raw_data[,.(safetyreportid,polypharmacy)] %>% 
  unique() %>% 
  .[,mean(polypharmacy)]
  
raw_data[,.(safetyreportid,polypharmacy)] %>% 
  unique() %>% 
  .[,quantile(polypharmacy,c(.025,0.975))]

# Reporting of ADEs --------------------------------------------------------------------

df <- raw_data[,
         .(safetyreportid,atc_concept_id,meddra_concept_id,nichd)
         ] %>% 
    unique() %>% 
    .[,
      .N,
      .(atc_concept_id,meddra_concept_id,nichd)
      ]

num_reports_vs_ades <- 
    df %>% 
    distinct(atc_concept_id,
             meddra_concept_id,N) %>% 
    filter(N!=0) %>% 
    group_by(N) %>% 
    count() %>% 
    rename(
        Num_Reports = N,
        Num_ADEs = n
    ) %>% 
    arrange(desc(Num_ADEs))

bins <- cut(num_reports_vs_ades$Num_Reports,
            breaks=c(0,1,2,3,4,seq(5,100,5),200),
            right = F)


num_reports_vs_ades$nreport = bins

num_reports_vs_ades_binned <- 
    num_reports_vs_ades %>% 
    group_by(nreport) %>% 
    summarise(
        nADEs = sum(Num_ADEs)
    )

num_reports_vs_ades_binned_nreport <- 
    as.character(num_reports_vs_ades_binned$nreport)

num_reports_vs_ades_binned_nreport[
    num_reports_vs_ades_binned_nreport=="[0,1)"
] <- "0"
num_reports_vs_ades_binned_nreport[
    num_reports_vs_ades_binned_nreport=="[1,2)"
] <- "1"
num_reports_vs_ades_binned_nreport[
    num_reports_vs_ades_binned_nreport=="[2,3)"
] <- "2"
num_reports_vs_ades_binned_nreport[
    num_reports_vs_ades_binned_nreport=="[3,4)"
] <- "3"
num_reports_vs_ades_binned_nreport[
    num_reports_vs_ades_binned_nreport=="[4,5)"
] <- "4"
num_reports_vs_ades_binned_nreport[
    is.na(num_reports_vs_ades_binned_nreport)
] <- ">200"

num_reports_vs_ades_binned$nreport <- 
    factor(num_reports_vs_ades_binned_nreport,
           levels = num_reports_vs_ades_binned_nreport)

g <- num_reports_vs_ades_binned %>% 
    ggplot(aes(x=nreport,y=nADEs/sum(nADEs))) +
    geom_bar(stat="identity") +
    geom_text(aes(
        label=paste0(round(nADEs/sum(nADEs),4)*100,"%"),
        vjust=-.5,fontface="bold"
    )) +
    scale_y_continuous(labels=scales::percent) +
    ylab("Proportion of drug-events") +
    xlab("Number of Reports") +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle=45,vjust=.7),
        text = element_text(face="bold",size=16)
    )
ggsave(paste0(img_dir,basename,"number_pediatric_reports_versus_proportion_of_ADEs_PT.png"),g,
       width=20,height=6)

# Reporting of ADEs about covariates ---------------------------------------

raw_data$ade <- paste0(raw_data$atc_concept_id,"_",raw_data$meddra_concept_id)

reports_in_stages <- raw_data[,
                            .(safetyreportid,nichd)
] %>% 
  unique() %>% 
  .[,
    .(total = .N),
    .(nichd)
  ] 

g <- raw_data[
  ,.(safetyreportid,nichd)
  ] %>% 
  unique() %>% 
  .[,
    .(frac_reports = .N / raw_data[,length(unique(safetyreportid))]),
    .(nichd)
  ]  %>% 
  ggplot(aes(factor(nichd,levels=category_levels$nichd),frac_reports)) +
  geom_bar(stat="identity",fill="gray",color="black") +
  xlab("") +
  scale_y_continuous(labels=scales::percent) +
  ylab("Percent of reports") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"reports_in_stages.png"),g,width=8,height=6)

g <- raw_data[
  ,.(safetyreportid,sex,nichd)
] %>% 
  unique() %>% 
  .[,
    .(N = .N),
    .(nichd,sex)
  ]  %>% 
  merge(reports_in_stages) %>% 
  .[,.(prop = N/total,nichd,sex)] %>% 
  ggplot(aes(factor(nichd,levels=category_levels$nichd),prop,fill=sex)) +
  geom_bar(stat="identity",position="dodge",color="black") +
  xlab("") +
  scale_y_continuous(labels=scales::percent) +
  guides(fill=guide_legend(title="Sex")) +
  ylab("Percent of reports") +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"report_sex_in_stages.png"),g,width=9,height=6)

reporter_in_stages <- 
  raw_data[
    ,.(safetyreportid,reporter_qualification)
  ] %>% 
  unique() %>% 
  .[,
    .(total = .N),
    .(reporter_qualification)
  ]

g <- raw_data[
  ,.(safetyreportid,reporter_qualification,nichd)
] %>% 
  unique() %>% 
  .[,
    .(N = .N),
    .(nichd,reporter_qualification)
  ]  %>% 
  merge(reporter_in_stages) %>% 
  .[,.(prop = N/total,nichd,reporter_qualification)] %>% 
  ggplot(aes(factor(nichd,levels=category_levels$nichd),prop,fill=reporter_qualification)) +
  geom_bar(stat="identity",position="dodge",color="black") +
  xlab("") +
  facet_grid(~reporter_qualification,scales="free",labeller = label_wrap_gen(width=25)) +
  guides(fill=guide_legend(title="Reporter qualification",title.position = "top",nrow=1)) +
  scale_y_continuous(labels=scales::percent) +
  ylab("Percent of reports") +
  theme(
    strip.background = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"reporter_in_stages.png"),g,width=16,height=6)


g <- raw_data[,
         .(receive_date,nichd,safetyreportid)
         ] %>% 
  unique() %>% 
  .[,.N,.(receive_date,nichd)] %>% 
  merge(reports_in_stages) %>% 
  .[order(receive_date)] %>% 
  .[,.(prop = N/total,nichd,receive_date)] %>% 
  ggplot(aes(receive_date,prop)) +
  geom_line(color="darkgray") +
  scale_y_continuous(labels=scales::percent) +
  facet_wrap(~factor(nichd,category_levels$nichd),scales="free",ncol=1) +
  xlab("Time") +
  ylab("Percent of reports") +
  theme(
    strip.background = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )
ggsave(paste0(img_dir,basename,"reports_over_time_in_stages.png"),g,width=8,height=16)

total <- 
  atc_covariates_full %>% 
  melt(id.vars=c(id_col)) %>% 
  .[,.(total = sum(value)),.(variable)]

g <- merge(
  raw_data[,.(safetyreportid,nichd)] %>% unique(),
  atc_covariates_full,by=id_col
  ) %>% 
  melt(id.vars=c(id_col,stage_col)) %>% 
  .[,.(S = sum(value)),.(nichd,variable)] %>% 
  merge(total) %>% 
  .[,.(prop = S/total,nichd,variable)] %>% 
  ggplot(aes(prop,variable)) +
  geom_bar(stat="identity",position="dodge",color="black",fill="gray") +
  facet_wrap(~factor(nichd,category_levels$nichd),scales="free_x",nrow=2) +
  scale_x_continuous(labels=scales::percent) +
  ylab("") +
  xlab("Percent of reports") +
  theme(
    strip.background = element_blank(),
    legend.position = "none"
  )
ggsave(paste0(img_dir,basename,"reports_in_atc_classes_in_stages.png"),g,width=15,height=10)
