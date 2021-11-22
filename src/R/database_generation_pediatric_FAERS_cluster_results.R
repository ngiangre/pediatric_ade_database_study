#' Title: "A database of pediatric drug effects to evaluate ontogenic mechanisms from child growth and development" study
#' 
#' Script author details: Nicholas Giangreco
#' 
#' This script generates the plots from our
#' data-driven clustering of drug-events in Pediatric FAERS 

# Purpose -----------------------------------------------------------------

#' Assign cluster categories to pediatric FAERS' drug-events
#' 


# Setup -------------------------------------------------------------------

pacman::p_load(data.table,tidyverse,doParallel)

data_dir <- "../../data/"
img_dir <- "../../docs/imgs/"
basename <- "database_generation_pediatric_FAERS_cluster_results_"
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

# load GAM data -----------------------------------------------------------

source("database_generation_load_GAM_data.R")

# load clustering data -----------------------------------------------------------

bp=5
X <-  
    fread(paste0(data_dir,
                 "database_generation_dtw_clustering_tuned_cluster_results_",
                 bp,"best_hyperparameter_set.csv.gz"))
cluster_names <- c("")
X_melt <- 
X %>% 
    melt(id.vars=c("ade","cluster","time","type","centroid","distance","preproc","k"),
         variable.name="nichd",value.name = "score")

norm_ades <- 
  normalize_data(dts[database=="covariate_adjusted" &
                       ade %in% common_ades],score="gam_score",groups=c())

clusters <- 
  1:length(X_melt[,unique(cluster)])
cluster_colors <- RColorBrewer::brewer.pal(length(clusters),name = "Dark2")
names(cluster_colors) <- clusters  

cluster_names = c("1" = "Plateau", "2" = "Increase", 
                  "3" = "Inverse Plateau", "4" = "Decrease")

# Load reference sets -----------------------------------------------------

pediatric_reference <- 
  fread("https://raw.githubusercontent.com/ngiangre/GRiP_pediatric_ADE-reference_set/master/data/GRiP_Pediatric_ADE_Gold_Standard_List_minimal_joined.csv")

pediatric_reference$ade <- 
  pediatric_reference[,paste0(ATC_concept_id,"_",MedDRA_concept_id)]

pediatric_reference$ade_name <- 
  pediatric_reference[,paste0(ATC_concept_name," and ",MedDRA_concept_name)]

pediatric_aes <- 
  fread(paste0(data_dir,"paediatric_term_list_19-0_concept_joined.csv"))

meddra_relationships <- 
  fread(paste0(data_dir,"standard_reactions_meddra_relationships.csv"))

pediatric_aes <- 
  merge(
    pediatric_aes,
    meddra_relationships,
    by.x="concept_id",
    by.y="meddra_concept_id_1"
  )

ryan_et_al <- 
  fread(paste0(data_dir,"database_generation_ryan_et_al.csv"))

ryan_et_al$ade <- 
  paste0(ryan_et_al$atc_concept_id,"_",ryan_et_al$meddra_concept_id)

drugbank_atc_cyp_substrates <- 
  fread(paste0(data_dir,"drugbank_atc_cyp_substrates.csv"))

drugbank_atc <- 
  fread(paste0(data_dir,"compound_drugbank05/drug_atc_codes_rxnorm_joined.csv"))

# Plot drug event clusters ------------------------------------------------

g <- norm_ades %>% 
  merge(
    X[,.(ade,cluster)] %>% 
      sample_n(1e4),
    by="ade"
  ) %>% 
  ggplot(aes(factor(nichd,levels=stages),norm,color=factor(cluster),group=ade)) +
  geom_line(alpha=0.1) +
  guides(color=guide_legend(title="Cluster\ncategory",override.aes = list(alpha = 1))) +
  xlab("") +
  ylab("Normalized dGAM score") +
  scale_color_manual(values=cluster_colors,labels=cluster_names) +
  facet_wrap(~cluster,labeller = labeller(cluster = cluster_names)) +
  theme(
    strip.background = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"sample_risk_patterns_in_clusters_",bp,"best_hyperparameter_set.png"),g,width=8,height=6)

# Proportion of significant drug-events in clusters ----------------------------------------

tmp <- X
tmp$sig_null <- tmp[,ade %in% sig_null_ades]

cluster_names
tmp[,.(N = .N),cluster]
tmp[,.(N = .N/nrow(tmp)),cluster]
tmp[,.(frac = sum(sig_null)),cluster]
tmp[,.(frac = sum(sig_null)/.N),cluster]

g <- tmp[,.(frac = sum(sig_null)/.N),cluster] %>% 
  ggplot(aes(factor(cluster),frac,fill=factor(cluster))) +
  geom_bar(stat="identity",color="black") +
  xlab("") +
  ylab("Proportoin of\nputative ADEs") +
  scale_fill_manual(values=cluster_colors,labels=cluster_names) +
  guides(fill=guide_legend(title="Cluster category",title.position = "top")) +
  scale_x_discrete(labels=cluster_names) +
  scale_y_continuous(labels=scales::percent) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"prop_significant_drug_events_in_clusters_",bp,"best_hyperparameter_set.png"),g,width=6,height=6)


# Talley of plateau putative ADEs at stages -------------------------------

tmp <- X[cluster==names(cluster_names[unname(cluster_names)=="Plateau"])]
tmp$sig_null <- tmp[,ade %in% sig_null_ades]

g <- tmp[
  sig_null==T
  ] %>% 
  pivot_longer(cols=stages) %>% 
  data.table() %>% 
  .[,.SD[which.max(value)],ade] %>% 
  .[,.N,name] %>% 
  .[order(factor(name,stages))] %>% 
  ggplot(aes(factor(name,stages),N)) +
  geom_bar(color="black",fill=cluster_colors[unname(cluster_names)=="Plateau"],
           stat="identity") +
  geom_text(aes(y=N,label=N),nudge_y = 100,fontface="bold") +
  xlab("") +
  ylab("Number of putative ADEs\nwith plateau dynamic") +
  scale_y_continuous(labels=scales::comma) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"prop_significant_plateau_drug_events_in_stages.png"),g,width=6,height=5)

# Cluster assignments for ADE class enrichments --------------------------

enriched <- 
  fread(paste0(data_dir,"database_generation_stage_enrichment_significant_atc5_enrichment.csv"))

enriched[18]

drug_class_map <- 
bind_rows(
  drugbank_atc[,
               .(atc_concept_id,
                 atc_concept_name = atc_concept_name
               )
  ] %>% unique(),
  drugbank_atc[,
               .(atc_concept_id,
                 atc_concept_name = level_1
               )
  ] %>% unique(),
  drugbank_atc[,
               .(atc_concept_id,
                 atc_concept_name = level_2
               )
  ] %>% unique(),
  drugbank_atc[,
               .(atc_concept_id,
                 atc_concept_name = level_3
               )
  ] %>% unique(),
  drugbank_atc[,
               .(atc_concept_id,
                 atc_concept_name = level_4
               )
  ] %>% unique()
)

event_class_map <- 
bind_rows(
  meddra_relationships[,.(
    meddra_concept_name = meddra_concept_name_1,
    meddra_concept_id = meddra_concept_id_1
  )] %>% unique(),
  meddra_relationships[,.(
    meddra_concept_name = meddra_concept_name_2,
    meddra_concept_id = meddra_concept_id_1
  )] %>% unique(),
  meddra_relationships[,.(
    meddra_concept_name = meddra_concept_name_3,
    meddra_concept_id = meddra_concept_id_1
  )] %>% unique(),
  meddra_relationships[,.(
    meddra_concept_name = meddra_concept_name_4,
    meddra_concept_id = meddra_concept_id_1
  )] %>% unique()
)

tmp <- 
merge(
  enriched,
  drug_class_map,
  by="atc_concept_name",
  allow.cartesian = T
) %>% 
  merge(
    event_class_map,
    by="meddra_concept_name",
    allow.cartesian = T
  ) %>% 
  .[,.(ade = paste0(atc_concept_id,"_",meddra_concept_id),nichd,lwr,odds_ratio,upr,fdr,
       atc_concept_name,meddra_concept_name)] %>% 
  unique() %>% 
  merge(
    X[ade %in% sig_null_ades,.(ade,cluster)] %>% 
      unique(),
    by="ade"
    )

tmp$cluster_name <- 
  sapply(tmp$cluster,function(x){cluster_names[x]})
  
frac_clusters <- merge(
  tmp[,.(ade,atc_concept_name,meddra_concept_name,nichd,cluster_name)] %>% 
    unique() %>% 
    .[,.N,.(atc_concept_name,meddra_concept_name,nichd,cluster_name)],
  tmp[,.(ade,atc_concept_name,meddra_concept_name,nichd)] %>% 
    unique() %>% 
    .[,.(total = .N),.(atc_concept_name,meddra_concept_name,nichd)]
) %>% 
  .[,.(atc_concept_name,meddra_concept_name,nichd,cluster_name,frac = N / total)] %>% 
  .[order(frac,decreasing = T)] %>%
  .[,.SD[1],.(atc_concept_name,meddra_concept_name,nichd)]

g <- merge(
  tmp[,.(ade_name = paste0(atc_concept_name," --- ",meddra_concept_name),
         ade,cluster)] %>% 
    unique(),
  norm_ades,
  by=c("ade"),
  allow.cartesian=T
) %>% 
  .[,.(nichd = factor(nichd,levels=stages),ade,cluster,ade_name,norm)] %>% 
  .[order(nichd)] %>% 
  .[grepl("montelukast",ade_name)] %>% 
  ggplot(aes(factor(nichd,levels=stages),norm,color=factor(cluster))) +
  geom_line(aes(group=ade),size=1) +
  scale_color_manual(values=cluster_colors,labels=cluster_names) +
  facet_wrap(~ade_name,
             labeller = label_wrap_gen(width=20),ncol=8) +
  xlab("") +
  ylab("Normalized dGAM score") +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"cluster_assignments_montelukast_psychiatric_disorders.png"),g,width=6,height=5)

merge(
  tmp[,.(ade_name = paste0(atc_concept_name," --- ",meddra_concept_name),
         ade,cluster_name)] %>% 
    unique(),
  norm_ades,
  by=c("ade"),
  allow.cartesian=T
) %>% 
  .[,.(nichd = factor(nichd,levels=stages),ade,cluster_name,ade_name,norm)] %>% 
  .[order(nichd)] %>% 
  .[,.N,cluster_name]

g <- merge(
  tmp[,.(ade_name = paste0(atc_concept_name," --- ",meddra_concept_name),
         ade,cluster)] %>% 
    unique(),
  norm_ades,
  by=c("ade"),
  allow.cartesian=T
) %>% 
  .[,.(nichd = factor(nichd,levels=stages),ade,cluster,ade_name,norm)] %>% 
  .[order(nichd)] %>% 
  ggplot(aes(factor(nichd,levels=stages),norm,color=factor(cluster))) +
  geom_line(aes(group=ade),size=1) +
  scale_color_manual(values=cluster_colors,labels=cluster_names) +
  facet_wrap(~ade_name,
             labeller = label_wrap_gen(width=20),ncol=8) +
  guides(color=guide_legend(title="Cluster",title.position = "top")) +
  xlab("") +
  ylab("Normalized dGAM score") +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"cluster_assignments_enriched_drug_risks.png"),g,width=20,height=15)
