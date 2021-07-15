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

# SOC enrichments ---------------------------------------------------------

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

category_test_dts[
  category=="soc" & meddra_concept_name!=""
] %>% 
  .[,.(meddra_concept_name,nichd = factor(nichd,levels=stages),
       lwr,odds_ratio,upr,fdr)
  ] %>% .[order(lwr,decreasing=T)] %>% .[fdr<0.05 & odds_ratio>1,.SD,nichd]

category_test_dts[
  category=="soc" & meddra_concept_name!=""
] %>% 
  .[,.(meddra_concept_name,nichd = factor(nichd,levels=stages),
       lwr,odds_ratio,upr,fdr)
  ] %>% .[order(lwr,decreasing=T)] %>% .[fdr<0.05 & odds_ratio>1,.SD,nichd] %>% dcast(meddra_concept_name ~ factor(nichd,levels=stages),value.var="lwr")

# ATC1 enrichments ---------------------------------------------------------

g <- category_test_dts[
  category=="atc1" & meddra_concept_name!=""
] %>% 
  .[,.(meddra_concept_name,nichd = factor(nichd,levels=stages),
       lwr,odds_ratio,upr,fdr)
  ] %>% 
  ggplot(aes(nichd,odds_ratio,group=meddra_concept_name)) +
  geom_bar(stat="identity",color="black",fill="gray") +
  geom_errorbar(aes(ymin=lwr,ymax=upr),width=0.2) +
  geom_hline(yintercept=1,color="red",linetype="dashed") +
  facet_wrap(~meddra_concept_name,labeller = label_wrap_gen(width=20),scales="free_y",ncol=5) +
  xlab("") +
  ylab("Odds ratio") +
  theme(
    strip.background = element_blank(),
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"atc1_risk_across_stages.png"),g,width=15,height=10)

category_test_dts[
  category=="atc1" & meddra_concept_name!=""
] %>% 
  .[,.(meddra_concept_name,nichd = factor(nichd,levels=stages),
       lwr,odds_ratio,upr,fdr)
  ] %>% .[order(lwr,decreasing=T)] %>% .[fdr<0.05 & odds_ratio>1,.SD,nichd]

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

tmp <- 
  category_test_dts[
    category %in% c("pt_cyp","hlt_cyp","hlgt_cyp","soc_cyp")
  ]

tmp[,length(unique(subcategory))]
tmp[fdr<0.05 & odds_ratio>1,length(unique(subcategory))]
tmp[fdr<0.05 & odds_ratio>1,.(subcategory,nichd)]

g <- tmp %>% 
  ggplot(aes(odds_ratio,-log10(pvalue))) + 
  geom_point(pch=21,color="black",fill="grey") +
  geom_point(
    data=category_test_dts[
      category %in% c("sider_cyp")
    ],
    pch=21,color="black",fill="red"
  ) +
  geom_point(
    data=category_test_dts[
      category %in% c("siderlabel_cyp")
    ],
    pch=21,color="black",fill="blue"
  ) +
  geom_point(
    data=category_test_dts[
      category %in% c("cyp_substrate")
    ],
    pch=21,color="black",fill="green"
  ) +
  scale_fill_manual(values=stage_colors) +
  scale_x_continuous(trans="log10",labels=scales::comma) +
  scale_y_continuous(trans="sqrt") +
  xlab("Odds ratio") +
  theme(
    legend.position = "bottom"
  )
ggsave(paste0(img_dir,basename,"cyp_class_volcano_plot.png"),
       g,width=5,height=5)

g <- tmp[odds_ratio>1 & fdr<0.05,.(N = length(unique(subcategory))),nichd] %>% 
  ggplot(aes(factor(nichd,levels=stages),N)) +
  geom_bar(stat="identity",color="black",fill="gray") +
  xlab("") +
  ylab("Number of significant\nstage enrichments") +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  ) 
ggsave(paste0(img_dir,basename,"number_of_significant_cyp_classes_at_stages.png"),g,width=6,height=5)

g <- 
  category_test_dts[
    category=="cyp_substrate" & 
      odds_ratio>1 & fdr<0.05,
    .(N = length(unique(subcategory))),nichd] %>% 
  ggplot(aes(factor(nichd,levels=stages),N)) +
  geom_bar(stat="identity",color="black",fill="gray") +
  xlab("") +
  ylab("Number of significant\nstage enrichments") +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  ) 
ggsave(paste0(img_dir,basename,"number_of_significant_cyp_risks_at_stages.png"),g,width=6,height=5)

g <- category_test_dts[
  subcategory %in% category_test_dts[category=="cyp_substrate" &
                                       odds_ratio>1 & fdr<0.05,unique(subcategory)]
] %>% 
  ggplot(aes(factor(nichd,levels=stages),lwr,group=subcategory)) +
  geom_line(size=3) +
  xlab("") +
  ylab("Odds ratio 95% lower bound") +
  facet_wrap(~meddra_concept_name,scales="free_y",labeller = label_wrap_gen(width=30)) +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1),
    strip.background = element_blank()
  ) 
ggsave(paste0(img_dir,basename,"cyp_enrichments_across_stages.png"),g,width=8,height=6)

# Transporter enrichments -------------------------------------------------

g <- category_test_dts[
  category=="transporter" & meddra_concept_name!="" &
    meddra_concept_name %in% category_test_dts[category=="transporter" & meddra_concept_name!="" & fdr<0.05,unique(meddra_concept_name)]
] %>% 
  .[,.(meddra_concept_name,nichd = factor(nichd,levels=stages),
       lwr,odds_ratio,upr,fdr)
  ] %>% 
  ggplot(aes(nichd,odds_ratio,group=meddra_concept_name)) +
  geom_bar(stat="identity",color="black",fill="gray") +
  geom_errorbar(aes(ymin=lwr,ymax=upr),width=0.2) +
  geom_hline(yintercept=1,color="red",linetype="dashed") +
  facet_wrap(~meddra_concept_name,labeller = label_wrap_gen(width=25),scales="free_y") +
  xlab("") +
  ylab("Odds ratio") +
  theme(
    strip.background = element_blank(),
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"transporter_risk_across_stages.png"),g,width=15,height=12)

tmp <- 
  category_test_dts[
    category %in% c("pt_transporter","hlt_transporter","hlgt_transporter","soc_transporter")
  ]

tmp[,length(unique(subcategory))]
tmp[fdr<0.05 & odds_ratio>1,length(unique(subcategory))]
tmp[fdr<0.05 & odds_ratio>1,.(subcategory,nichd,percent_events_in_stage = a/(a+c),lwr,odds_ratio,upr,fdr)]

g <- tmp %>% 
  ggplot(aes(odds_ratio,-log10(pvalue))) + 
  geom_point(pch=21,color="black",fill="grey") +
  scale_fill_manual(values=stage_colors) +
  scale_x_continuous(trans="log10",labels=scales::comma) +
  scale_y_continuous(trans="sqrt") +
  xlab("Odds ratio") +
  theme(
    legend.position = "bottom"
  )
ggsave(paste0(img_dir,basename,"transporter_class_volcano_plot.png"),
       g,width=5,height=5)

tmp[odds_ratio>1 & fdr<0.05]


tmp <- 
  category_test_dts[
    category %in% c('transporter')
  ]

tmp[,length(unique(subcategory))]
tmp[odds_ratio>1 & fdr<0.05,length(unique(subcategory))]
tmp[odds_ratio>1 & fdr<0.05]

g <- 
  category_test_dts[
    category=="transporter" & 
      odds_ratio>1 & fdr<0.05,
    .(N = length(unique(subcategory))),nichd] %>% 
  ggplot(aes(factor(nichd,levels=stages),N)) +
  geom_bar(stat="identity",color="black",fill="gray") +
  xlab("") +
  ylab("Number of significant\nstage enrichments") +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  ) 
ggsave(paste0(img_dir,basename,"number_of_significant_transporter_risks_at_stages.png"),g,width=6,height=5)

g <- category_test_dts[
  subcategory %in% tmp[odds_ratio>1 & fdr<0.05,unique(subcategory)]
] %>% 
  ggplot(aes(factor(nichd,levels=stages),lwr,group=subcategory)) +
  geom_line(size=3) +
  xlab("") +
  ylab("Odds ratio 95% lower bound") +
  facet_wrap(~meddra_concept_name,scales="free_y",labeller = label_wrap_gen(width=30)) +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1),
    strip.background = element_blank()
  ) 
ggsave(paste0(img_dir,basename,"transporter_enrichments_across_stages.png"),g,width=18,height=10)

# ADE enrichments ---------------------------------------------------

tmp <- 
  category_test_dts %>% 
  .[category %in% c("soc","hlgt","hlt","pt","atc1","atc2","atc3","atc4")]

tmp[grepl("ATC[12345]",meddra_concept_class_id),length(unique(meddra_concept_name))]
tmp[meddra_concept_class_id %in% c("PT","HLT","HLGT","SOC"),length(unique(meddra_concept_name))]
tmp[,length(unique(subcategory))]

tmp %>% 
  fwrite(paste0(data_dir,basename,"ade_data.csv"))
tmp[,length(unique(subcategory))]
tmp[subcategory!="",.(category,subcategory)] %>% unique() %>% .[,.N,category]

g <- tmp %>% 
  ggplot(aes(odds_ratio,-log10(pvalue))) + 
  geom_point(pch=21,color="black",fill="grey") +
  scale_fill_manual(values=stage_colors) +
  scale_x_continuous(trans="log10",labels=scales::comma) +
  scale_y_continuous(trans="sqrt") +
  xlab("Odds ratio") +
  theme(
    legend.position = "bottom"
  )
ggsave(paste0(img_dir,basename,"ade_volcano_plot.png"),
       g,width=5,height=5)

g <- tmp %>% 
  ggplot(aes(a,odds_ratio,fill=-log10(pvalue))) + 
  geom_point(pch=21,color="black") +
  scale_fill_gradient2(low="blue",mid = "white", high="red", 
                       midpoint = (tmp[,min(-log10(pvalue))] + tmp[,max(-log10(pvalue))])/2) +
  xlab("Number of significant drug-events") +
  ylab("Odds ratio") +
  scale_x_continuous(trans="log10",labels=scales::comma) +
  scale_y_continuous(trans="log10",labels=scales::comma) +
  facet_wrap(~factor(nichd,levels=stages),nrow=2)
ggsave(paste0(img_dir,basename,"number_and_odds_correlation_ade.png"),
       g,width=12,height=5)


tmp[odds_ratio>1 & fdr<0.05,length(unique(subcategory))]
tmp[odds_ratio>1 & fdr<0.05,unique(subcategory)]

tmp[odds_ratio>1 & fdr<0.05,length(unique(subcategory))] /
  tmp[,length(unique(subcategory))]

tmp[a>5 & fdr<0.05 & lwr>1 &
      grepl("ATC[12345]",meddra_concept_class_id)
    ] %>% 
  .[order(fdr)] %>% 
  .[,.SD[1],nichd] %>% 
  .[,.(subcategory,
       percent_events_in_stage = round(a/(a+c),2),
       nichd = factor(nichd,levels=stages),
       lwr = round(lwr,2),
       odds_ratio = round(odds_ratio,2),
       upr = round(upr,2),
       fdr)] %>% 
  .[order(nichd)]

tmp[a>5 & fdr<0.05 & lwr>1 &
      meddra_concept_class_id %in% c("PT","HLT","HLGT","SOC")
    ] %>% 
  .[order(fdr)] %>% 
  .[,.SD[1],nichd] %>% 
  .[,.(subcategory,
       percent_events_in_stage = round(a/(a+c),2),
       nichd = factor(nichd,levels=stages),
       lwr = round(lwr,2),
       odds_ratio = round(odds_ratio,2),
       upr = round(upr,2),
       fdr)] %>% 
  .[order(nichd)]


g <- 
  tmp[
    fdr<0.05 & 
      odds_ratio>1 &
      grepl("ATC[12345]",meddra_concept_class_id)
    ,.(N = length(unique(subcategory))),.(nichd,category)] %>% 
  ggplot(aes(factor(nichd,levels=stages),N)) +
  geom_bar(stat="identity",color="black",fill="gray") +
  xlab("") +
  ylab("Number of significant\nstage enrichments") +
  facet_grid(category~.,scales="free_y",labeller = label_wrap_gen(30)) +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  ) 
ggsave(paste0(img_dir,basename,"number_of_significant_drug_classes_at_stages.png"),g,width=5,height=10)

g <- 
  tmp[
    fdr<0.05 & 
      odds_ratio>1 &
      meddra_concept_class_id %in% c("PT","HLT","HLGT","SOC")
    ,.(N = length(unique(subcategory))),.(nichd,category)] %>% 
  ggplot(aes(factor(nichd,levels=stages),N)) +
  geom_bar(stat="identity",color="black",fill="gray") +
  xlab("") +
  ylab("Number of significant\nstage enrichments") +
  facet_grid(factor(category,levels=c("soc","hlgt","hlt","pt"))~.,scales="free_y",labeller = label_wrap_gen(30)) +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  ) 
ggsave(paste0(img_dir,basename,"number_of_significant_event_classes_at_stages.png"),g,width=5,height=10)

g <- 
  atc5_cats[,
            .(atc_concept_name,meddra_concept_name)
  ] %>% 
  unique() %>% 
  merge(
    tmp,
    by=c("atc_concept_name","meddra_concept_name")
  ) %>% 
  .[,.(nichd = factor(nichd,levels=stages),subcategory,
       lwr,soc = meddra_concept_name)] %>% 
  ggplot(aes(factor(nichd,levels=stages),lwr,group=subcategory,color=soc)) +
  geom_line(size=3) +
  facet_wrap(~subcategory,labeller = label_wrap_gen(30)) +
  xlab("") +
  ylab("Odds ratio 95% lower bound") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle=45,vjust=1,hjust=1),
    strip.background = element_blank()
  ) 
ggsave(paste0(img_dir,basename,"odds_lwr_enriched_drugs_SOC_across_stages.png"),g,width=20,height=15)


# ADE class enrichments ---------------------------------------------------

tmp <- 
  category_test_dts %>% 
  .[category %in% c("soc_atc1","soc_atc2","soc_atc3","soc_atc4","soc_atc5",
                    "hlgt_atc1","hlgt_atc2","hlgt_atc3","hlgt_atc4","hlgt_atc5",
                    "hlt_atc1","hlt_atc2","hlt_atc3","hlt_atc4","hlt_atc5",
                    "pt_atc1","pt_atc2","pt_atc3","pt_atc4")]

tmp[atc_concept_class_id!="ATC5",length(unique(atc_concept_name))]
tmp[meddra_concept_class_id!="PT",length(unique(atc_concept_name))]
tmp[,length(unique(subcategory))]

tmp %>% 
  fwrite(paste0(data_dir,basename,"ade_class_data.csv"))
tmp[subcategory!="",length(unique(subcategory))]
tmp[subcategory!="",.(category,subcategory)] %>% unique() %>% .[,.N,category]

g <- tmp %>% 
  ggplot(aes(odds_ratio,-log10(pvalue))) + 
  geom_point(pch=21,color="black",fill="grey") +
  scale_fill_manual(values=stage_colors) +
  scale_x_continuous(trans="log10",labels=scales::comma) +
  scale_y_continuous(trans="sqrt") +
  xlab("Odds ratio") +
  theme(
    legend.position = "bottom"
  )
ggsave(paste0(img_dir,basename,"ade_class_volcano_plot.png"),
       g,width=5,height=5)

g <- tmp %>% 
  ggplot(aes(a,odds_ratio,fill=-log10(pvalue))) + 
  geom_point(pch=21,color="black") +
  scale_fill_gradient2(low="blue",mid = "white", high="red", 
                       midpoint = (tmp[,min(-log10(pvalue))] + tmp[,max(-log10(pvalue))])/2) +
  xlab("Number of significant drug-events") +
  ylab("Odds ratio") +
  scale_x_continuous(trans="log10",labels=scales::comma) +
  scale_y_continuous(trans="log10",labels=scales::comma) +
  facet_wrap(~factor(nichd,levels=stages),nrow=2)
ggsave(paste0(img_dir,basename,"number_and_odds_correlation_ade_class.png"),
       g,width=12,height=5)


tmp[odds_ratio>1 & fdr<0.05,length(unique(subcategory))]
tmp[odds_ratio>1 & fdr<0.05,unique(subcategory)]

tmp[odds_ratio>1 & fdr<0.05,length(unique(subcategory))] /
  tmp[,length(unique(subcategory))]

tmp[subcategory=="Endocrine disorders --- DIRECT ACTING ANTIVIRALS" & nichd=="infancy"]

tmp[odds_ratio>1 & fdr<0.05 & lwr>1] %>% 
  .[order(fdr)] %>% 
  .[,.SD[1],nichd] %>% 
  .[,.(subcategory,
       percent_events_in_stage = round(a/(a+c),2),
       nichd = factor(nichd,levels=stages),
       lwr = round(lwr,2),
       odds_ratio = round(odds_ratio,2),
       upr = round(upr,2),
       fdr)] %>% 
  .[order(nichd)]

atc5_cats <- 
  tmp[grepl("atc5",category) & fdr<0.05 & odds_ratio>1,
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

atc_cats <- 
  tmp[fdr<0.05 & odds_ratio>1,
      .(atc_concept_name,meddra_concept_name,
        nichd,a,
        lwr = round(lwr,2),
        odds_ratio = round(odds_ratio,2),
        upr = round(upr,2),fdr = fdr,
        atc_concept_class_id,meddra_concept_class_id)]
atc_cats$nichd <- factor(atc_cats$nichd,levels=stages)
setorderv(atc_cats,c("nichd","fdr","odds_ratio"))
atc_cats[,length(unique(atc_concept_name))]
atc_cats %>% 
  fwrite(paste0(data_dir,basename,"significant_atc_enrichment.csv"))

g <- tmp[fdr<0.05 & odds_ratio>1,.(N = length(unique(subcategory))),nichd] %>% 
  ggplot(aes(factor(nichd,levels=stages),N)) +
  geom_bar(stat="identity",color="black",fill="gray") +
  xlab("") +
  ylab("Number of significant\nstage enrichments") +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  ) 
ggsave(paste0(img_dir,basename,"number_of_significant_ades_at_stages.png"),g,width=6,height=5)

g <- tmp[fdr<0.05 & odds_ratio>1 & atc_concept_class_id=="ATC5",
         .(N = length(unique(subcategory))),nichd] %>% 
  ggplot(aes(factor(nichd,levels=stages),N)) +
  geom_bar(stat="identity",color="black",fill="gray") +
  geom_label(aes(label=N)) +
  xlab("") +
  ylab("Number of significant\nstage enrichments") +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  ) 

ggsave(paste0(img_dir,basename,"number_of_significant_atc5_ades_at_stages.png"),g,width=6,height=5)

g <- tmp[fdr<0.05 & odds_ratio>1,.(N = length(unique(subcategory))),
    .(nichd,atc_concept_class_id,meddra_concept_class_id)] %>% 
  ggplot(aes(factor(nichd,levels=stages),N)) +
  geom_bar(stat="identity",color="black",fill="gray") +
  xlab("") +
  ylab("Number of significant\nstage enrichments") +
  scale_y_continuous(trans="sqrt") +
  facet_grid(meddra_concept_class_id~atc_concept_class_id) +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  ) 
ggsave(paste0(img_dir,basename,"number_of_significant_ades_at_stages_by_classes.png"),g,width=10,height=6)

g <- 
  atc5_cats[,
            .(atc_concept_name,meddra_concept_name)
            ] %>% 
  unique() %>% 
  merge(
    tmp,
    by=c("atc_concept_name","meddra_concept_name")
  ) %>% 
  .[,.(nichd = factor(nichd,levels=stages),subcategory,
       lwr,odds_ratio,upr,
       soc = meddra_concept_name)] %>% 
  ggplot(aes(factor(nichd,levels=stages),lwr,group=subcategory)) +
  geom_line(size=3) +
  facet_wrap(~subcategory,labeller = label_wrap_gen(30)) +
  geom_hline(yintercept=1,color="red",linetype="dashed") +
  xlab("") +
  ylab("Odds ratio 95% lower bound") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle=45,vjust=1,hjust=1),
    strip.background = element_blank()
  ) 
ggsave(paste0(img_dir,basename,"odds_lwr_enriched_drugs_across_stages.png"),g,width=20,height=15)


# drug enrichments --------------------------------------------------------

drugbank_atc <- 
  fread(paste0(data_dir,"compound_drugbank05/drug_atc_codes_rxnorm_joined.csv"))


risky_drugs <- category_test_dts[
  grepl("^atc",category) & 
    odds_ratio>1 & 
    fdr<0.05,length(unique(meddra_concept_name))]
total <- category_test_dts[
  grepl("^atc",category),length(unique(meddra_concept_name))]

risky_drugs
total
risky_drugs / total


g <- category_test_dts[
  grepl("^atc",category) & 
    odds_ratio>1 & 
    fdr<0.05,.N,nichd] %>% 
  ggplot(aes(factor(nichd,levels=stages),N)) +
  geom_bar(stat="identity",color="black",fill="gray") +
  xlab("") +
  ylab("Number of risky drugs and drug classes") +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"risky_drugs_and_classes_across_stages.png"),g,width=8,height=7)

atc5_risky_drugs <- 
  category_test_dts[
    category=="atc5",
    .(atc_concept_name = meddra_concept_name,nichd,odds_ratio,fdr)
  ] %>% 
  unique() %>% 
  merge(
    drugbank_atc,
    by="atc_concept_name",
    allow.cartesian = T
  ) %>% 
  .[,
    .(
      nichd,
      odds_ratio,fdr,
      atc5 = atc_concept_name,
      atc4 = level_1,atc3 = level_2,
      atc2 = level_3,atc1 = level_4
    )
  ] %>% 
  .[odds_ratio>1 & fdr<0.05] %>% 
  unique()

atc5_risky_drugs[,length(unique(atc5))]
atc5_risky_drugs[,length(unique(atc4))]
atc5_risky_drugs[,length(unique(atc3))]
atc5_risky_drugs[,length(unique(atc2))]
atc5_risky_drugs[,length(unique(atc1))]

category_test_dts[
  category=="atc4" & 
    odds_ratio>1 & 
    fdr<0.05,
  .(level_1 = meddra_concept_name)
] %>% unique()

category_test_dts[
  category=="atc3" & 
    odds_ratio>1 & 
    fdr<0.05,
  .(level_2 = meddra_concept_name)
] %>% unique() 

category_test_dts[
  category=="atc2" & 
    odds_ratio>1 & 
    fdr<0.05,
  .(level_3 = meddra_concept_name)
] %>% unique()

atc1_risky_classes <- 
  setdiff(
    category_test_dts[
      category=="atc1" & 
        lwr>1 & 
        fdr<0.05,meddra_concept_name
    ] %>% unique(),category_test_dts[
      category=="atc1" & 
        lwr<1 & 
        nichd=="late_adolescence",meddra_concept_name
    ])

g <- category_test_dts[
  subcategory %in% atc1_risky_classes
] %>% 
  ggplot(aes(factor(nichd,levels=stages),lwr,group=subcategory)) +
  geom_line(size=3) +
  xlab("") +
  ylab("Odds ratio 95% lower bound") +
  facet_wrap(~meddra_concept_name,scales="free_y",labeller = label_wrap_gen(width=30)) +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1),
    strip.background = element_blank()
  ) 
ggsave(paste0(img_dir,basename,"odds_lwr_enriched_atc1_risk_across_stages.png"),g,width=8,height=5)

atc2_risky_classes <- 
  setdiff(
    category_test_dts[
      category=="atc2" & 
        lwr>1 & 
        fdr<0.05,meddra_concept_name
    ] %>% unique(),category_test_dts[
      category=="atc2" & 
        lwr<1 & 
        nichd=="late_adolescence",meddra_concept_name
    ])

g <- category_test_dts[
  subcategory %in% atc2_risky_classes
] %>% 
  ggplot(aes(factor(nichd,levels=stages),lwr,group=subcategory)) +
  geom_line(size=3) +
  xlab("") +
  ylab("Odds ratio 95% lower bound") +
  facet_wrap(~meddra_concept_name,scales="free_y",labeller = label_wrap_gen(width=30)) +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1),
    strip.background = element_blank()
  ) 
ggsave(paste0(img_dir,basename,"odds_lwr_enriched_atc2_risk_across_stages.png"),g,width=15,height=10)

atc3_risky_classes <- 
  setdiff(
    category_test_dts[
      category=="atc3" & 
        lwr>1 & 
        fdr<0.05,meddra_concept_name
    ] %>% unique(),category_test_dts[
      category=="atc3" & 
        lwr<1 & 
        nichd=="late_adolescence",meddra_concept_name
    ])

g <- category_test_dts[
  category=="atc3" & subcategory %in% atc3_risky_classes
] %>% 
  ggplot(aes(factor(nichd,levels=stages),lwr,group=subcategory)) +
  geom_line(size=3) +
  xlab("") +
  ylab("Odds ratio 95% lower bound") +
  facet_wrap(~meddra_concept_name,scales="free_y",labeller = label_wrap_gen(width=20)) +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1),
    strip.background = element_blank()
  ) 
ggsave(paste0(img_dir,basename,"odds_lwr_enriched_atc3_risk_across_stages.png"),g,width=15,height=10)

atc4_risky_classes <- 
  setdiff(
    category_test_dts[
      category=="atc4" & 
        lwr>1 & 
        fdr<0.05,meddra_concept_name
    ] %>% unique(),category_test_dts[
      category=="atc4" & 
        lwr<1 & 
        nichd=="late_adolescence",meddra_concept_name
    ])

g <- category_test_dts[
  category=="atc4" & subcategory %in% atc4_risky_classes
] %>% 
  ggplot(aes(factor(nichd,levels=stages),lwr,group=subcategory)) +
  geom_line(size=3) +
  xlab("") +
  ylab("Odds ratio 95% lower bound") +
  facet_wrap(~meddra_concept_name,scales="free_y",labeller = label_wrap_gen(width=30)) +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1),
    strip.background = element_blank()
  ) 
ggsave(paste0(img_dir,basename,"odds_lwr_enriched_atc4_risk_across_stages.png"),g,width=25,height=15)

atc5_risky_classes <- 
  setdiff(
    category_test_dts[
      category=="atc5" & 
        lwr>1 & 
        fdr<0.05,meddra_concept_name
    ] %>% unique(),category_test_dts[
      category=="atc5" & 
        lwr<1 & 
        nichd=="late_adolescence",meddra_concept_name
    ])

g <- category_test_dts[
  category=="atc5" & subcategory %in% atc5_risky_classes
] %>% 
  ggplot(aes(factor(nichd,levels=stages),lwr,group=subcategory)) +
  geom_line(size=3) +
  xlab("") +
  ylab("Odds ratio 95% lower bound") +
  facet_wrap(~meddra_concept_name,scales="free_y",labeller = label_wrap_gen(width=30)) +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1),
    strip.background = element_blank()
  ) 
ggsave(paste0(img_dir,basename,"odds_lwr_enriched_atc5_risk_across_stages.png"),g,width=25,height=20)

# Valproic acid investigations --------------------------------------------

drugbank_atc_cyp_substrates <- 
  fread(paste0(data_dir,"drugbank_atc_cyp_substrates.csv"))

drugbank_atc <- 
  fread(paste0(data_dir,"compound_drugbank05/drug_atc_codes_rxnorm_joined.csv"))

drugbank_atc[atc_concept_name=="valproic acid",.(atc_concept_name,level_1,level_2,level_3,level_4)] %>% unique()
meddra_relationships[meddra_concept_name_4=="Congenital, familial and genetic disorders",
                     sample(unique(meddra_concept_name_1),10)]
events <- meddra_relationships[meddra_concept_name_4=="Congenital, familial and genetic disorders",
                     unique(meddra_concept_name_1)]
category_test_dts[
  meddra_concept_class_id=="PT" & 
    meddra_concept_name %in% events & 
    atc_concept_name=="NERVOUS SYSTEM" & fdr<0.05]

id_ <- drugbank_atc[atc_concept_name=="valproic acid",unique(atc_concept_id)]
drugbank_atc_cyp_substrates[atc_concept_id==id_]
es_ <- drugbank_atc_cyp_substrates[atc_concept_id==id_,enzyme_name]
category_test_dts[
  atc_concept_name %in% es_ & 
    meddra_concept_name=="Congenital, familial and genetic disorders" &
    nichd=="term_neonatal",
  .(subcategory,a,lwr,odds_ratio,fdr)
  ]
category_test_dts[
  atc_concept_name %in% es_ & 
    grepl("siderlabel",category) &
    nichd=="term_neonatal",
  .(subcategory,a,lwr,odds_ratio,fdr)
]


# enrichment of CYP metabolized drugs and liver related events?  NO --------

category_test_dts[
  meddra_concept_name=="Hepatobiliary disorders" &
    atc_concept_class_id=="CYP"] %>% 
  .[order(fdr)] %>% 
  .[fdr<0.05 & odds_ratio>1]

# SIDER adverse events and CYP enzymes ------------------------------------

g <- category_test_dts[
  category %in% c("sider_cyp"),
  .(atc_concept_name,meddra_concept_name,nichd,score=odds_ratio)
  ] %>% 
  .[
    is.finite(score)
    ] %>% 
  ggplot(aes(factor(nichd,levels=stages),score)) +
  geom_boxplot() +
  geom_hline(yintercept=1,color="red",linetype="dashed") +
  xlab("") +
  facet_wrap(~atc_concept_name,labeller = label_wrap_gen(width=30)) +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"sider_by_cyp_odds_across_stages.png"),g,width=20,height=15)

