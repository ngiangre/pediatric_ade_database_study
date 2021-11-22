#' Title: "A database of pediatric drug effects to evaluate ontogenic mechanisms from child growth and development" study
#' 
#' Script author details: Nicholas Giangreco
#' 
#' This script generates the plots of the distributions and characteristics for dGAMs and their comparison to the PRR

# Purpose -----------------------------------------------------------------

#' Tabulate drug-events and risk distributions for dGAMs
#' 

# Setup -------------------------------------------------------------------

pacman::p_load(data.table,tidyverse,doParallel)

data_dir <- "../../data/"
img_dir <- "../../docs/imgs/"
basename <- "database_generation_with_covariates_"
seed = 0
set.seed(seed)

registerDoParallel(cores=5)

theme_set(theme_bw(base_size=16) + theme(text = element_text(face="bold")))

# Load data and covariates ------------------------------------------------

source("database_generation_load_data.R")

# Load GAM data ------------------------------------------------

source("database_generation_load_GAM_data.R")


# Functions ---------------------------------------------------------------

source("database_generation_functions.R")

# Number of drugs, aes, and ades  ------------------------------------

dts[database=="covariate_adjusted",length(unique(ade))]
dts[database=="covariate_adjusted",length(unique(meddra_concept_id))]
dts[database=="covariate_adjusted",length(unique(atc_concept_id))]

dts[database=="covariate_adjusted" & gam_score_90mse>0,length(unique(ade))]

merge(
  dts[database=="covariate_adjusted"],
  null_dts[,.(null_99 = quantile(gam_score,c(0.99))),nichd],
  by="nichd"
) %>% 
  .[gam_score_90mse>null_99,length(unique(ade))]

# Comparing PRR and GAM scores correlating to scores at stages  ------------------------

sub <- dts[
  database=="covariate_adjusted" &
    ade %in% sample(sig_null_ades,2e3,replace=F), 
  .(ade,nichd,gam_score,PRR)
]

sub[is.na(PRR),"PRR"] <- 0

cor_dts <- NULL
for(st1 in stages){
  for(st2 in stages){
    cor_ = cor.test(sub[nichd==st1,gam_score],sub[nichd==st2,gam_score],method="spearman")
    cor_dts <- 
      bind_rows(
        cor_dts,
        data.table(nichd_1 = st1, nichd_2 = st2,cor_estimate = cor_$estimate,cor_pvalue = cor_$p.value,score="GAM")
      )
  }
}
for(st1 in stages){
  for(st2 in stages){
    cor_ = cor.test(sub[nichd==st1,PRR],sub[nichd==st2,PRR],method="spearman")
    cor_dts <- 
      bind_rows(
        cor_dts,
        data.table(nichd_1 = st1, nichd_2 = st2,cor_estimate = cor_$estimate,cor_pvalue = cor_$p.value,score="PRR")
      )
  }
}

cor_dts[score=="GAM",score:= list("dGAM")]
g <- cor_dts %>% 
  ggplot(aes(factor(nichd_1,levels=stages),
             factor(nichd_2,levels=stages),
             fill=cor_estimate)) +
  geom_tile() +
  scale_fill_viridis_c(guide = guide_legend()) +
  facet_wrap(~score) +
  guides(fill = guide_colourbar(barwidth = 12,title.position = "top",title="Spearman correlation")) +
  xlab("") +
  ylab("") +
  theme(
    strip.background = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"PRR_and_GAM_score_stage_correlation.png"),g,width=7,height=5)


# Comparing PRR and GAM dynamics across stages -----------------------------------------

sub <- 
  dts[database=="covariate_adjusted"]
sub[is.na(PRR)]$PRR <- 0
norm_ades_gam <- 
  normalize_data(sub,
                 score="gam_score",strata="nichd",groups = c())
norm_ades_gam$type <- "dGAM"
norm_ades_prr <- 
  normalize_data(sub,
                 score="PRR",strata="nichd",groups = c())
norm_ades_prr$type <- "PRR"

g <- 
  bind_rows(
    norm_ades_prr,
    norm_ades_gam
  ) %>%
  .[,.(ade,nichd,norm,type)] %>% 
  unique() %>% 
  .[ade %in% sample(sig_null_ades,2e3,replace=F)] %>%  
  ggplot(aes(nichd,norm)) +
  geom_line(aes(group=ade),alpha=0.1) +
  facet_grid(type~.) +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  ) +
  xlab("") +
  ylab("Normalized score") +
  theme(
    strip.background = element_blank()
  )
ggsave(paste0(img_dir,basename,"gam_versus_prr_norm_score_across_stages.png"),g,width=6,height=6)

# ADE risk distributions per stage ---------------------------------------

sub <- 
  dts[database=="covariate_adjusted"]

seq_= seq(-10,10,1)
Nsamp=10000
Nsim = 100
tmp <- foreach(i=1:Nsim,.combine='rbind') %dopar%{
  set.seed(i)
  sinds <- sample(1:nrow(sub),Nsamp,replace=F)
  tmp <- dts[sinds] 
  tmp %>% 
    .[,
      .(m = mean(gam_score)),
      .(nichd,ade)
    ] %>% 
    .[,.(nichd,ade,m,
         m_cat = cut(m,breaks=seq_,right=F,ordered_result=T),
         bootstrap =  i)] %>% 
    na.omit() %>% 
    .[,.(N = .N),.(nichd,m_cat)] %>% 
    merge(
      tmp[,.(total = length(unique(ade))),.(nichd)]
    ) %>% 
    .[,.(N = N/total,nichd,m_cat)]
}

#gam_score
g <- tmp[,
         .(lwr = quantile(N,c(0.025)),
           m = mean(N),
           upr = quantile(N,c(0.975))),
         .(nichd,m_cat)
] %>% 
  ggplot(aes(m_cat,m)) +
  geom_bar(color="black",fill="grey",stat="identity",
           width = 1, position = position_nudge(x = 0.5)) +
  geom_errorbar(aes(ymin=lwr,ymax=upr),
                position = position_nudge(x = 0.5),width=0.2) +
  scale_fill_brewer(palette = "Dark2") +
  facet_wrap(~factor(nichd,levels=category_levels$nichd)) +
  scale_y_continuous(labels=scales::percent) +
  geom_vline(aes(xintercept=11),color="red", linetype="dashed", size=1) +
  xlab("dGAM score categories") +
  ylab("Percent of drug-events") +
  theme(
    strip.background = element_blank(),
    legend.position = 'none',
    text = element_text(size=22),
    axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=12)
  )
ggsave(paste0(img_dir,basename,"gam_score_distributions_-10to10_boot.png"),g,width=15,height=8)


# Percent nominally significant risks across stages ---------------------------------

g <- dts[
  database=="covariate_adjusted" & 
    gam_score_90mse>0,
  .(N = length(unique(ade))),
  nichd] %>% 
  merge(
    dts[
      database=="covariate_adjusted",
      .(total = length(unique(ade))),
      nichd
    ],
    by="nichd"
  ) %>% 
  .[,.(nichd,N,frac = N/total)] %>% 
  ggplot(aes(nichd,frac)) +
  geom_bar(stat="identity",color="black",fill="grey") +
  xlab("") +
  ylab("Percent of nominally significant\ndrug-events") +
  scale_y_continuous(labels=scales::percent) +
  theme(
    text = element_text(size=20),
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,'significant_risks_across_stages.png'),g,width=6,height=7)

# Significance from report shuffling - for defining putative ADEs --------------------------------------

tmp <- merge(
  dts[database=="covariate_adjusted"],
  null_dts[,.(null_99 = quantile(gam_score,c(0.99))),nichd],
  by="nichd"
) %>% 
  .[,.(ade,sig_null=as.integer(gam_score_90mse>null_99),sig_risk = as.integer(gam_score_90mse>0)),nichd] %>% 
  .[,.N,.(nichd,sig_de_risk = paste0(sig_null,sig_risk))] 

g <- tmp %>% 
  .[order(factor(nichd,levels=stages))] %>% 
  .[sig_de_risk!="00"] %>% 
  dcast(nichd ~ sig_de_risk,value.var="N") %>% 
  .[,.(nichd,frac = `11`/`01`)] %>% 
  ggplot(aes(factor(nichd,levels=stages),frac)) +
  geom_bar(stat="identity",position="dodge",color="black",fill="gray") +
  scale_y_continuous(labels=scales::percent) +
  xlab("") +
  ylab("Percent of significant drug-event risks\nnominally and by the null model") +
  theme(
    text = element_text(size=20),
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,'null_simulation_percent_significant_risks_across_stages.png'),g,width=6,height=7)

g <- merge(
  dts[database=="covariate_adjusted"],
  null_dts[,.(null_99 = quantile(gam_score,c(0.99))),nichd],
  by="nichd"
) %>% 
  .[gam_score_90mse>null_99,.(N = length(unique(ade))),nichd] %>% 
  .[,.(nichd = factor(nichd,levels=category_levels$nichd),N)] %>% 
  .[order(nichd)] %>% 
  merge(
    dts[
      database=="covariate_adjusted",
      .(total = length(unique(ade))),
      nichd
    ],
    by="nichd"
  ) %>% 
  .[,.(nichd,N,frac = N/total)] %>% 
  ggplot(aes(nichd,frac)) +
  geom_bar(stat="identity",color="black",fill="grey") +
  xlab("") +
  ylab("Percent of significant drug-events\nby the null model") +
  scale_y_continuous(labels=scales::percent) +
  theme(
    text = element_text(size=20),
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,'null_simulation_significant_risks_across_stages.png'),g,width=6,height=7)


# putative ADEs across stages by socs --------------------------------------

meddra_relationships <- 
  fread(paste0(data_dir,"standard_reactions_meddra_relationships.csv"))

g <- merge(
  dts[database=="covariate_adjusted"],
  null_dts[,.(null_99 = quantile(gam_score,c(0.99))),nichd],
  by="nichd"
) %>% 
  .[gam_score_90mse>null_99] %>% 
  merge(
    meddra_relationships[,.(meddra_concept_id = meddra_concept_id_1,
                            soc = meddra_concept_name_4)],
    by="meddra_concept_id",
    allow.cartesian = T
  ) %>% 
  .[gam_score_90mse>null_99,.(N = length(unique(ade))),.(soc,nichd)] %>% 
  .[,.(nichd = factor(nichd,levels=category_levels$nichd),soc,N)] %>% 
  .[order(nichd)] %>% 
  merge(
    merge(
      dts[database=="covariate_adjusted"],
      null_dts[,.(null_99 = quantile(gam_score,c(0.99))),nichd],
      by="nichd"
    ) %>% 
      .[gam_score_90mse>null_99] %>% 
      merge(
        meddra_relationships[,.(meddra_concept_id = meddra_concept_id_1,
                                soc = meddra_concept_name_4)],
        by="meddra_concept_id",
        allow.cartesian = T
      ) %>% 
      .[gam_score_90mse>null_99,.(total = length(unique(ade))),.(soc)],
    by="soc"
  ) %>% 
  .[,.(nichd,soc,N,frac = N/total)] %>% 
  ggplot(aes(factor(nichd,levels=stages),frac)) +
  geom_bar(stat="identity",color="black",fill="grey") +
  xlab("") +
  ylab("Percent of putative ADEs") +
  scale_y_continuous(labels=scales::percent) +
  facet_wrap(~soc,nrow=3,labeller = label_wrap_gen(width=20)) +
  theme(
    strip.background = element_blank(),
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,'null_simulation_significant_risks_across_stages_by_soc.png'),
       g,width=20,height=10)

# putative ADEs across stages by atc1s --------------------------------------

drugbank_atc <- 
  fread(paste0(data_dir,"compound_drugbank05/drug_atc_codes_rxnorm_joined.csv"))

g <- 
  merge(
    dts[database=="covariate_adjusted"],
    null_dts[,.(null_99 = quantile(gam_score,c(0.99))),nichd],
    by="nichd"
  ) %>% 
  .[gam_score_90mse>null_99] %>% 
  merge(
    drugbank_atc[,.(atc_concept_id,atc1 = level_4)],
    by="atc_concept_id",
    allow.cartesian = T
  ) %>% 
  .[gam_score_90mse>null_99,.(N = length(unique(ade))),.(atc1,nichd)] %>% 
  .[,.(nichd = factor(nichd,levels=category_levels$nichd),atc1,N)] %>% 
  .[order(nichd)] %>% 
  merge(
    merge(
      dts[database=="covariate_adjusted"],
      null_dts[,.(null_99 = quantile(gam_score,c(0.99))),nichd],
      by="nichd"
    ) %>% 
      .[gam_score_90mse>null_99] %>% 
      merge(
        drugbank_atc[,.(atc_concept_id,atc1 = level_4)],
        by="atc_concept_id",
        allow.cartesian = T
      ) %>% 
      .[gam_score_90mse>null_99,.(total = length(unique(ade))),.(atc1)],
    by="atc1"
  ) %>% 
  .[,.(nichd,atc1,N,frac = N/total)] %>% 
  ggplot(aes(nichd,frac)) +
  geom_bar(stat="identity",color="black",fill="grey") +
  xlab("") +
  ylab("Percent of putative ADEs") +
  scale_y_continuous(labels=scales::percent) +
  facet_wrap(~atc1,nrow=2,labeller = label_wrap_gen(width=25)) +
  theme(
    strip.background = element_blank(),
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,'null_simulation_significant_risks_across_stages_by_atc1.png'),
       g,width=20,height=6)


# putative ADEs across stages by CYP enzymes ----------------------------------------------

drugbank_atc_cyp_substrates <- 
  fread(paste0(data_dir,"drugbank_atc_cyp_substrates.csv"))

drugbank_atc <- 
  fread(paste0(data_dir,"compound_drugbank05/drug_atc_codes_rxnorm_joined.csv"))

meddra_relationships <- 
  fread(paste0(data_dir,"standard_reactions_meddra_relationships.csv"))

sub <- merge(
  dts[database=="covariate_adjusted"],
  null_dts[,.(null_99 = quantile(gam_score,c(0.99))),nichd],
  by="nichd"
) %>% 
  .[gam_score_90mse>null_99] %>% 
  merge(
    meddra_relationships[,
                         .(meddra_concept_id = meddra_concept_id_1,
                           soc = meddra_concept_name_4)
    ],
    by="meddra_concept_id",
    allow.cartesian=T
  ) %>% 
  merge(
    drugbank_atc[,
                 .(atc_concept_id,
                   atc1 = level_4,
                   atc2 = level_3,
                   atc3 = level_2,
                   atc4 = level_1)
    ] %>% unique(),
    by="atc_concept_id"
  ) %>% 
  merge(
    drugbank_atc_cyp_substrates[,
                                .(atc_concept_id,enzyme_name)
    ] %>% unique(),
    by="atc_concept_id",
    allow.cartesian = T
  )

total <- 
  sub[
    ,.(ade,enzyme_name)
  ] %>% 
  unique() %>% 
  .[,
    .(
      total = length(unique(ade))
    ),
    .(enzyme_name)
  ] 

tmp <- 
  sub[gam_score_90mse>null_99,
      .(
        N = length(unique(ade))
      ),
      .(nichd,enzyme_name)
  ] %>% 
  merge(
    total,
    by="enzyme_name"
  ) %>% 
  .[,.(nichd,enzyme_name,frac_risk = N / total)]

setorderv(tmp,c("nichd","frac_risk"),c(1,-1))
g <- tmp %>% 
  ggplot(aes(factor(nichd,levels=stages),frac_risk)) +
  geom_bar(stat="identity",color="black",fill="grey") +
  xlab("") +
  ylab("Percent of significant drug-events") +
  scale_y_continuous(labels=scales::percent) +
  facet_wrap(~enzyme_name,ncol=5,labeller = label_wrap_gen(width=25)) +
  theme(
    strip.background = element_blank(),
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"cyp_drug_risks_across_stages.png"),g,width=15,height=12)


# Outputting clinically significant ades using sider ---------------------------------

drugbank_atc <- 
  fread(paste0(data_dir,"compound_drugbank05/drug_atc_codes_rxnorm_joined.csv"))

sider_drug_atc <- 
  fread(paste0(data_dir,"sider/drug_atc.tsv"),sep="\t")
colnames(sider_drug_atc) <- c("stitch_id","atc_concept_code")

sider_drug_atc_id <- 
  merge(
    sider_drug_atc,
    drugbank_atc[,.(atc_concept_code = atc_code,atc_concept_id,atc_concept_name)] %>% 
      unique(),
    by="atc_concept_code"
  )

meddra_all_label_se <- 
  fread(paste0(data_dir,"sider/meddra_all_label_se.tsv"),sep="\t")
colnames(meddra_all_label_se) <- c("source_label","stitch_id","umls_source_value","medgen_id","meddra_concept_class_id","medgen_id","meddra_concept_name")

sider <- 
  merge(
    meddra_all_label_se[
      meddra_concept_class_id=="PT",
      .(stitch_id,medgen_id,meddra_concept_name)] %>% 
      unique(),
    fread(paste0(data_dir,"standard_reactions_meddra_relationships.csv")) %>% 
      .[,
        .(
          meddra_concept_name = meddra_concept_name_1,
          meddra_concept_id = meddra_concept_id_1,
          soc = meddra_concept_name_4
        )
      ] %>% unique(),
    by="meddra_concept_name",
    allow.cartesian = T
  ) %>% 
  merge(
    sider_drug_atc_id,
    by="stitch_id",
    allow.cartesian = T
  )

sider$ade <- paste0(sider$atc_concept_id,"_",sider$meddra_concept_id)
sider$ade_name <- paste0(sider$atc_concept_name," and ",sider$meddra_concept_name)

merge(
  sider,
  dts[,.(ade)] %>% unique(),
  by=c("ade")
  ) %>% 
  fwrite(paste0(data_dir,"database_generation_sider_standardized_named_data.csv"))

