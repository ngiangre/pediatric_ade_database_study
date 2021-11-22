#' Title: "A database of pediatric drug effects to evaluate ontogenic mechanisms from child growth and development" study
#' 
#' Author details: Nicholas Giangreco
#' 
#' This script generates the plots for the stage enrichment data

# Purpose -----------------------------------------------------------------

#' To evaluate stage enrichment of drug-events from categories
#' 


# Setup -------------------------------------------------------------------

pacman::p_load(data.table,tidyverse)

data_dir <- "../../data/"
img_dir <- "../../docs/imgs/"
basename <- "database_generation_stage_enrichment_"
seed = 0
set.seed(seed)

stages <- 
    c("term_neonatal","infancy",
      "toddler","early_childhood",
      "middle_childhood","early_adolescence",
      "late_adolescence")

theme_set(theme_bw(base_size=16) + theme(text = element_text(face="bold")))


# Load functions ----------------------------------------------------------

source("database_generation_functions.R")


# load stage enrichment data ---------------------------------------------

category_test_dts <- 
  fread(paste0(data_dir,"database_generation_stage_enrichment_results.csv"))
category_test_dts$atc_concept_class_id <- 
  sapply(str_split(category_test_dts$category,"_"),function(x){str_to_upper(x[2])})
category_test_dts$meddra_concept_class_id <- 
  sapply(str_split(category_test_dts$category,"_"),function(x){str_to_upper(x[1])})
category_test_dts$atc_concept_name <- 
  sapply(str_split(category_test_dts$subcategory," --- "),function(x){x[2]})
category_test_dts$meddra_concept_name <- 
  sapply(str_split(category_test_dts$subcategory," --- "),function(x){x[1]})

category_test_dts[,table(category)]

stage_colors <- RColorBrewer::brewer.pal(length(stages),name = "Dark2")
names(stage_colors) <- stages

# Tally classes -------------------------------------------------------

category_test_dts[
  a>0 &
    category %in% c("atc1","atc2","atc3","atc4","atc5") & 
    subcategory!="",length(unique(subcategory))
]
category_test_dts[
  a>0 &
    category %in% c("pt","hlt","hlgt","soc") & 
    subcategory!="",length(unique(subcategory))
]

category_test_dts[
  a>0 &
  category %in% c("soc_atc1","soc_atc2","soc_atc3","soc_atc4","soc_atc5",
                  "hlgt_atc1","hlgt_atc2","hlgt_atc3","hlgt_atc4","hlgt_atc5",
                  "hlt_atc1","hlt_atc2","hlt_atc3","hlt_atc4","hlt_atc5",
                  "pt_atc1","pt_atc2","pt_atc3","pt_atc4"),
  length(unique(subcategory))]


# drug class enrichments ---------------------------------------------------------

#Number
category_test_dts[
  a>0 & odds_ratio>1 & fdr<0.05 &
    category %in% c("atc1","atc2","atc3","atc4","atc5") & 
    subcategory!="",.(subcategory,meddra_concept_class_id)
] %>% 
  unique() %>% 
  nrow()

#Number by stage
category_test_dts[
  a>0 & odds_ratio>1 & fdr<0.05 &
    category %in% c("atc1","atc2","atc3","atc4","atc5") & 
    subcategory!="",
  .(subcategory,meddra_concept_class_id,nichd)
  ] %>% 
  unique() %>% 
  .[,.(N = length(unique(subcategory))),
  .(nichd,meddra_concept_class_id)
] %>% 
  unique() %>%
  .[,.(N = sum(N)),nichd] %>% 
  .[order(factor(nichd,levels=stages))]
#Number in adolescence
category_test_dts[
  a>0 & odds_ratio>1 & fdr<0.05 &
    category %in% c("atc1","atc2","atc3","atc4","atc5") & 
    (nichd=="early_adolescence" | nichd=="late_adolescence") &
    subcategory!="",.(N = length(unique(subcategory))),
  .(nichd,meddra_concept_class_id)
][,.(N = sum(N))]

#Table
category_test_dts[
  a>0 & odds_ratio>1 & is.finite(odds_ratio) & fdr<0.05 &
    category %in% c("atc1","atc2","atc3","atc4","atc5") & 
    subcategory!=""] %>% 
  .[order(odds_ratio,decreasing = T)] %>% 
  .[,.SD[1],.(nichd,meddra_concept_class_id)] %>% 
  .[,
    .(drug_class = meddra_concept_name,
      nichd,a,
      lwr = round(lwr,2),
      odds_ratio = round(odds_ratio,2),
      upr = round(upr,2),
      pvalue = format(pvalue,scientific=T),
      fdr = format(fdr,scientific=T),
      atc_concept_class_id = meddra_concept_class_id)
    ] %>% 
  .[order(factor(nichd,stages))] %>% 
  fwrite(paste0(data_dir,basename,"top_drug_classes_in_stages.csv"))

category_test_dts[
  a>0 & odds_ratio>1 & fdr<0.05 &
    category %in% c("atc1","atc2","atc3","atc4","atc5") & 
    subcategory!="",length(unique(subcategory)),meddra_concept_class_id
]

g <- category_test_dts[
  category=="atc1" & subcategory!=""
] %>% 
  .[,.(meddra_concept_name,nichd = factor(nichd,levels=stages),
       lwr,odds_ratio,upr,fdr)
  ] %>% 
  ggplot(aes(nichd,odds_ratio,group=meddra_concept_name)) +
  geom_bar(stat="identity",color="black",fill="gray") +
  geom_errorbar(aes(ymin=lwr,ymax=upr),width=0.2) +
  geom_hline(yintercept=1,color="red",linetype="dashed") +
  facet_wrap(~meddra_concept_name,labeller = label_wrap_gen(width=20),scales="free_y",nrow=2) +
  xlab("") +
  ylab("Odds ratio") +
  theme(
    strip.background = element_blank(),
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"atc1_risk_across_stages.png"),g,width=20,height=6)

# adverse event enrichments ---------------------------------------------------------

#Number
category_test_dts[
  a>0 & odds_ratio>1 & fdr<0.05 &
    category %in% c("pt","hlt","hlgt","soc") & 
    subcategory!="",.(subcategory,meddra_concept_class_id)
] %>% 
  unique() %>% 
  nrow()

#Number by stage
category_test_dts[
  a>0 & odds_ratio>1 & fdr<0.05 &
    category %in% c("pt","hlt","hlgt","soc") & 
    subcategory!="",
  .(subcategory,meddra_concept_class_id,nichd)
] %>% 
  unique() %>% 
  .[,.(N = length(unique(subcategory))),
    .(nichd,meddra_concept_class_id)
  ] %>% 
  unique() %>%
  .[,.(N = sum(N)),nichd] %>% 
  .[order(factor(nichd,levels=stages))]
(category_test_dts[
  a>0 & odds_ratio>1 & fdr<0.05 &
    category %in% c("pt","hlt","hlgt","soc") & 
    subcategory!="",
  .(subcategory,meddra_concept_class_id,nichd)
] %>% 
    unique() %>% 
    .[,.(N = length(unique(subcategory))),
      .(nichd,meddra_concept_class_id)
    ] %>% 
    unique() %>%
    .[,.(N = sum(N)),nichd] %>% 
    .[order(factor(nichd,levels=stages))] %>% 
    .[,N])/302

#Table
category_test_dts[
  a>0 & odds_ratio>1 & is.finite(odds_ratio) & fdr<0.05 &
    category %in% c("pt","hlt","hlgt","soc") & 
    subcategory!=""] %>% 
  .[order(odds_ratio,decreasing = T)] %>% 
  .[,.SD[1],.(nichd,meddra_concept_class_id)] %>% 
  .[,
    .(adverse_event_class = meddra_concept_name,
      nichd,a,
      lwr = round(lwr,2),
      odds_ratio = round(odds_ratio,2),
      upr = round(upr,2),
      pvalue = format(pvalue,scientific=T),
      fdr = format(fdr,scientific=T),
      meddra_concept_class_id)
  ] %>% 
  .[order(factor(meddra_concept_class_id,c("PT","HLT","HLGT","SOC")))] %>% 
  .[order(factor(nichd,stages))] %>% 
  fwrite(paste0(data_dir,basename,"top_adverse_event_classes_in_stages.csv"))

category_test_dts[
  a>0 & odds_ratio>1 & fdr<0.05 &
    category %in% c("pt","hlt","hlgt","soc") & 
    subcategory!="",length(unique(subcategory)),meddra_concept_class_id
]

g <- category_test_dts[
  category=="soc" & meddra_concept_name!=""
  ] %>% 
  .[,.(meddra_concept_name,nichd = factor(nichd,levels=stages),
    lwr,odds_ratio,upr,fdr)
    ] %>% 
  ggplot(aes(nichd,odds_ratio,group=meddra_concept_name)) +
  geom_bar(stat="identity",color="black",fill="gray") +
  geom_errorbar(aes(ymin=lwr,ymax=upr),width=0.2) +
  geom_hline(yintercept=1,color="red",linetype="dashed") +
  facet_wrap(~meddra_concept_name,labeller = label_wrap_gen(width=20),scales="free_y",nrow=3) +
  xlab("") +
  ylab("Odds ratio") +
  theme(
    strip.background = element_blank(),
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"soc_risk_across_stages.png"),g,width=20,height=10)

# ADE class combination enrichments ---------------------------------------------------

category_test_dts[
  a>0 & odds_ratio>1 & fdr<0.05 &
    category %in% c("soc_atc1","soc_atc2","soc_atc3","soc_atc4","soc_atc5",
                    "hlgt_atc1","hlgt_atc2","hlgt_atc3","hlgt_atc4","hlgt_atc5",
                    "hlt_atc1","hlt_atc2","hlt_atc3","hlt_atc4","hlt_atc5",
                    "pt_atc1","pt_atc2","pt_atc3","pt_atc4"),
  length(unique(subcategory))]

g <- category_test_dts[
  a>0 & odds_ratio>1 & fdr<0.05 &
    category %in% c("soc_atc1","soc_atc2","soc_atc3","soc_atc4","soc_atc5",
                    "hlgt_atc1","hlgt_atc2","hlgt_atc3","hlgt_atc4","hlgt_atc5",
                    "hlt_atc1","hlt_atc2","hlt_atc3","hlt_atc4","hlt_atc5",
                    "pt_atc1","pt_atc2","pt_atc3","pt_atc4"),
  .(N = length(unique(subcategory))),
  .(atc_concept_class_id,
    meddra_concept_class_id)
  ] %>% 
  bind_rows(
    data.table(
      atc_concept_class_id = c("ATC1","ATC2","ATC3","ATC4","ATC5"),
      meddra_concept_class_id = rep("PT",5),
      N = c(rep(0,4),NA)
    )
  ) %>% 
  .[,.(N,
       atc_concept_class_id,
       meddra_concept_class_id = 
         factor(meddra_concept_class_id,level=c("PT","HLT","HLGT","SOC")))] %>% 
  ggplot(aes(atc_concept_class_id,meddra_concept_class_id)) +
  geom_tile(aes(fill=N)) +
  colorspace::scale_fill_continuous_divergingx() +
  geom_text(aes(label=N),fontface="bold",size=16) +
  xlab("") +
  ylab("") +
  theme(
    legend.position = "none",
    line = element_blank(),
    axis.line = element_line(size=0)
  )
ggsave(paste0(img_dir,basename,"ADE_class_enrichments_heatmap.png"),g,width=10,height=5)

# atc4 enrichments ---------------------------------------------

atc4_cats <- 
  category_test_dts %>% 
  .[category %in% c("soc_atc4")] %>% 
  .[fdr<0.05 & odds_ratio>1,
    .(atc_concept_name,meddra_concept_name,
      nichd,a,
      lwr = round(lwr,2),
      odds_ratio = round(odds_ratio,2),
      upr = round(upr,2),fdr = round(fdr,4),
      atc_concept_class_id,meddra_concept_class_id)]
atc4_cats$nichd <- factor(atc4_cats$nichd,levels=stages)
setorderv(atc4_cats,c("nichd","fdr"))
atc4_cats[,length(unique(atc_concept_name))]
atc4_cats[fdr<0.05,length(unique(atc_concept_name)),nichd][order(factor(nichd,levels=stages))]

atc4_cats %>% 
  fwrite(paste0(data_dir,basename,"significant_atc4_soc_enrichment.csv"))

tmp <- 
  atc4_cats[,
            .(atc_concept_name,meddra_concept_name)
  ] %>% 
  unique() %>% 
  merge(
    category_test_dts %>% 
      .[grepl("_atc4",category)],
    by=c("atc_concept_name","meddra_concept_name")
  ) %>% 
  .[,.(atc_concept_name,
       nichd = factor(nichd,levels=stages),
       subcategory,
       lwr,odds_ratio,upr,
       class = meddra_concept_name)
  ]

mycolors <- 
  colorRampPalette(
    RColorBrewer::brewer.pal(name="Dark2", n = 8)
  )(nrow(unique(tmp[,.(atc_concept_name,class)])))

g <- 
  tmp %>% 
  ggplot(
    aes(factor(nichd,levels=stages),lwr,
        group=atc_concept_name,
        color=atc_concept_name)
  ) +
  geom_line(size=2) +
  facet_wrap(~class,labeller = label_wrap_gen(25),scales="free_y") +
  ggrepel::geom_label_repel(
    data=tmp[,.SD[which.max(lwr)],.(atc_concept_name,class)],
    aes(label=str_wrap(atc_concept_name, width=30)),
    na.rm = F,
    force=80,fontface="bold"
  ) +
  scale_color_manual(values=mycolors) +
  geom_hline(yintercept=1,color="red",linetype="dashed") +
  xlab("") +
  ylab("Odds ratio 95% lower bound") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle=45,vjust=1,hjust=1),
    strip.background = element_blank()
  ) 
ggsave(paste0(img_dir,basename,"odds_lwr_enriched_atc4_across_stages_colored.png"),g,width=20,height=18)

# atc5 enrichments ---------------------------------------------

atc5_cats <- 
  category_test_dts %>% 
  .[grepl("_atc5",category) & fdr<0.05 & odds_ratio>1,
    .(atc_concept_name,meddra_concept_name,
      nichd,a,
      lwr = round(lwr,2),
      odds_ratio = round(odds_ratio,2),
      upr = round(upr,2),fdr = round(fdr,4),
      atc_concept_class_id,meddra_concept_class_id)]
atc5_cats$nichd <- factor(atc5_cats$nichd,levels=stages)
setorderv(atc5_cats,c("nichd","fdr"))
atc5_cats[,length(unique(atc_concept_name))]
atc5_cats %>% 
  fwrite(paste0(data_dir,basename,"significant_atc5_soc_enrichment.csv"))

tmp <- 
  atc5_cats[,
            .(atc_concept_name,meddra_concept_name)
  ] %>% 
  unique() %>% 
  merge(
    category_test_dts %>% 
      .[grepl("_atc5",category)],
    by=c("atc_concept_name","meddra_concept_name")
  ) %>% 
  .[,.(atc_concept_name,
       nichd = factor(nichd,levels=stages),
       subcategory,
       lwr,odds_ratio,upr,
       class = meddra_concept_name)
  ]

mycolors <- 
  colorRampPalette(
    RColorBrewer::brewer.pal(name="Dark2", n = 8)
    )(nrow(unique(tmp[,.(atc_concept_name,class)])))

g <- 
  tmp %>% 
  ggplot(
    aes(factor(nichd,levels=stages),lwr,
        group=atc_concept_name,
        color=atc_concept_name)
  ) +
  geom_line(size=2) +
  facet_wrap(~class,labeller = label_wrap_gen(20),scales="free_y",nrow=3) +
  ggrepel::geom_label_repel(
    data=tmp[,.SD[which.max(lwr)],.(atc_concept_name,class)],
    aes(label=str_wrap(atc_concept_name, width=30)),
    na.rm = F,
    force=20,fontface="bold"
  ) +
  scale_color_manual(values=mycolors) +
  geom_hline(yintercept=1,color="red",linetype="dashed") +
  xlab("") +
  ylab("Odds ratio 95% lower bound") +
  theme(
    legend.position = "none",
    axis.text = element_text(color="black"),
    axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=12),
    strip.background = element_blank(),
    strip.text = element_text(size=18)
  ) 
ggsave(paste0(img_dir,basename,"odds_lwr_enriched_drugs_across_stages_colored.png"),g,width=20,height=10)

# CYP enrichments ---------------------------------------------------------

g <- category_test_dts[
  category=="cyp_substrate" & meddra_concept_name!="" &
    meddra_concept_name %in% category_test_dts[category=="cyp_substrate" & meddra_concept_name!="" & fdr<0.05,unique(meddra_concept_name)]
] %>% 
  .[,.(meddra_concept_name,nichd = factor(nichd,levels=stages),
       lwr,odds_ratio,upr,fdr)
  ] %>% 
  ggplot(aes(nichd,odds_ratio,group=meddra_concept_name)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=lwr,ymax=upr),width=0.2) +
  geom_hline(yintercept=1,color="red",linetype="dashed") +
  facet_wrap(~meddra_concept_name,labeller = label_wrap_gen(width=25),scales="free_y") +
  xlab("") +
  ylab("Odds ratio") +
  theme(
    strip.background = element_blank(),
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"cyp_risk_across_stages.png"),g,width=15,height=12)

category_test_dts[
  category=="cyp_substrate" & meddra_concept_name!="" &
    meddra_concept_name %in% category_test_dts[category=="cyp_substrate" & meddra_concept_name!="" & fdr<0.05,unique(meddra_concept_name)] & odds_ratio>1 & fdr<0.05
] %>% 
  .[,.(meddra_concept_name,nichd = factor(nichd,levels=stages),
       lwr = round(lwr,2),
       odds = round(odds_ratio,2),
       upr = round(upr,2),
       FDR = round(fdr,2)
  )
  ] %>% 
  .[order(factor(nichd,stages))]
