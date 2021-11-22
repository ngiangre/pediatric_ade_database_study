#' Title: "A database of pediatric drug effects to evaluate ontogenic mechanisms from child growth and development" study
#' 
#' Script author details: Nicholas Giangreco
#' 
#' This script generates the data and database for this study

# Purpose -----------------------------------------------------------------

#' To make tables that will make up the final database for this study
#' 


# Setup -------------------------------------------------------------------

#install.packages("pacman")
pacman::p_load(tidyverse,data.table,mgcv,doParallel,progress,ROCR)

data_dir <- "../../data/"
img_dir <- "../../docs/imgs/"
seed = 0
set.seed(seed)

registerDoParallel(cores=5)

basename="database_generation_er_tables"
out_dir <- paste0(data_dir,basename,"/")

ifelse(
    !dir.exists(file.path(out_dir)),
    dir.create(file.path(out_dir)),
    FALSE)

# Load functions ----------------------------------------------------------

source("database_generation_functions.R")

# Load data -----------------------------------------------------------

source("database_generation_load_data.R")

# Load GAM data -----------------------------------------------------------

source("database_generation_load_GAM_data.R")

# Load clustering data -----------------------------------------------------------

bp=5
X <-  
    fread(paste0(data_dir,
                 "database_generation_dtw_clustering_tuned_cluster_results_",
                 bp,"best_hyperparameter_set.csv.gz"))

clusters <- 
    1:length(X[,unique(cluster)])
cluster_colors <- RColorBrewer::brewer.pal(length(clusters),name = "Dark2")
names(cluster_colors) <- clusters  

cluster_names = c("1" = "Plateau", "2" = "Increase", 
                  "3" = "Inverse Plateau", "4" = "Decrease")


# Load drugbank data ------------------------------------------------------

drugbank_atc <- 
  fread(paste0(data_dir,"compound_drugbank05/drug_atc_codes_rxnorm_joined.csv"))

# Load sider side effects -------------------------------------------------

sider <- 
  fread(paste0(data_dir,"database_generation_sider_standardized_named_data.csv"))

# Load meddra relationships -----------------------------------------------

meddra_relationships <- 
  fread(paste0(data_dir,"standard_reactions_meddra_relationships.csv"))

# er_tables list ----------------------------------------------------------

er_tables <- list()

# Primary keys ------------------------------------------------------------

primary_keys <- common_ades

# DRUG --------------------------------------------------------------------

tmp <- 
  raw_data[,
           .(atc_concept_id,safetyreportid)
  ] %>% 
  unique() %>% 
  .[,.(ndrugreports = .N),atc_concept_id] %>% 
  merge(
    drugbank_atc[,
                 .(atc_concept_id,atc_concept_name,atc_concept_code = atc_code,
                   atc4_concept_name = level_1,atc4_concept_code = code_1,
                   atc3_concept_name = level_2,atc3_concept_code = code_2,
                   atc2_concept_name = level_3,atc2_concept_code = code_3,
                   atc1_concept_name = level_4,atc1_concept_code = code_4
                 )] %>% 
      unique(),
    by=drug_col,
    all.x=T
  )

#er_tables[["drug"]] <- tmp


# DRUG --- recent ATC vocabulary ------------------------------------------

fread("../../vocabulary_download_v5_RxNorm_ATC_20210825/CONCEPT.csv") %>% 
  .[concept_class_id=='ATC 5th',length(unique(concept_id))]

atc5_concept_table <- 
  fread("../../vocabulary_download_v5_RxNorm_ATC_20210825/CONCEPT.csv") %>% 
  .[concept_class_id=='ATC 5th' & concept_id %in% tmp[,atc_concept_id],
    .(atc_concept_id = concept_id,
      atc_concept_name = concept_name,
      atc_concept_code = concept_code)
  ] %>% 
  merge(
    raw_data[,
             .(atc_concept_id,safetyreportid)
    ] %>% 
      unique() %>% 
      .[,.(ndrugreports = .N),atc_concept_id],
    by="atc_concept_id",
    all.y=T
  )

atc5_concept_table[,.N,atc_concept_name][N>1] #3 names mapping to 2 concept ids

atc_anc_desc <- 
  fread("../../vocabulary_download_v5_RxNorm_ATC_20210825/CONCEPT_ANCESTOR.csv") %>% 
  .[descendant_concept_id %in% atc5_concept_table[,atc_concept_id],
    .(ancestor_concept_id,
      atc_concept_id = descendant_concept_id)
  ]

atc_anc_desc_long <- 
  atc_anc_desc %>% 
  merge(
    atc5_concept_table,
    by="atc_concept_id",
    all.y=T
  ) %>% 
  merge(
    fread("../../vocabulary_download_v5_RxNorm_ATC_20210825/CONCEPT.csv") %>% 
      .[concept_class_id %in% c("ATC 4th","ATC 3rd","ATC 2nd","ATC 1st"),
        .(concept_id,concept_name,concept_code,concept_class_id)],
    all.x=T,
    by.x="ancestor_concept_id",
    by.y="concept_id"
  )

tmp <- 
  atc5_concept_table %>% 
  merge(
    atc_anc_desc_long[
      concept_class_id=="ATC 4th",
      .(atc_concept_id,atc_concept_name,atc_concept_code,
        atc4_concept_name = concept_name,
        atc4_concept_code = concept_code,
        ndrugreports)
    ],
    all.x = T,
    by=c("atc_concept_id","atc_concept_name","atc_concept_code","ndrugreports")
  ) %>% 
  merge(
    atc_anc_desc_long[
      concept_class_id=="ATC 3rd",
      .(atc_concept_id,atc_concept_name,atc_concept_code,
        atc3_concept_name = concept_name,
        atc3_concept_code = concept_code,
        ndrugreports)
    ],
    all.x = T,
    by=c("atc_concept_id","atc_concept_name","atc_concept_code","ndrugreports")
  ) %>% 
  merge(
    atc_anc_desc_long[
      concept_class_id=="ATC 2nd",
      .(atc_concept_id,atc_concept_name,atc_concept_code,
        atc2_concept_name = concept_name,
        atc2_concept_code = concept_code,
        ndrugreports)
    ],
    all.x = T,
    by=c("atc_concept_id","atc_concept_name","atc_concept_code","ndrugreports")
  ) %>% 
  merge(
    atc_anc_desc_long[
      concept_class_id=="ATC 1st",
      .(atc_concept_id,atc_concept_name,atc_concept_code,
        atc1_concept_name = concept_name,
        atc1_concept_code = concept_code,
        ndrugreports)
    ],
    all.x = T,
    by=c("atc_concept_id","atc_concept_name","atc_concept_code","ndrugreports")
  )

er_tables[["drug"]] <- tmp


# DRUG - corresponding to rxnorm ids - not for db -------------------------

tmp <- 
  er_tables[["drug"]] %>% 
  .[,.(atc_concept_id,atc_concept_name)] %>% 
  merge(
    fread("../../vocabulary_download_v5_RxNorm_ATC_20210825/CONCEPT_RELATIONSHIP.csv") %>% 
      .[,.(atc_concept_id = concept_id_1,
           concept_id = concept_id_2,
           relationship_id)],
    by="atc_concept_id"
  ) %>% 
  merge(
    fread("../../vocabulary_download_v5_RxNorm_ATC_20210825/CONCEPT.csv") %>% 
      .[,.(concept_id,concept_name,concept_class_id,vocabulary_id)],
    by="concept_id"
  )

#tmp[grepl("hydrocortisone",atc_concept_name)] %>% View()

# EVENT -------------------------------------------------------------------

soc_category <-  fread(paste0(data_dir,"database_generation_gene_expression_evaluation_soc_categories.csv"))

tmp <- 
  raw_data[,
           .(meddra_concept_id,safetyreportid)
  ] %>% 
  unique() %>% 
  .[,.(neventreports = .N),meddra_concept_id] %>% 
  merge(
    meddra_relationships,
    by.x=rx_col,
    by.y=paste0(rx_col,"_1"),
    all.x=T
  ) %>% 
  merge(
    soc_category[,.(meddra_concept_name_4 = soc,soc_category = category)],
    by="meddra_concept_name_4",
    all.x = T
  )

ped_aes <- 
  fread(paste0(data_dir,"paediatric_term_list_19-0_concept_joined.csv"))

tmp$pediatric_adverse_event <- 
  as.integer(tmp$meddra_concept_id %in% ped_aes$concept_id)

er_tables[["event"]] <- tmp


# ADE_RAW -----------------------------------------------------------------

tmp <- 
  raw_data[,
           c("safetyreportid","ade",
             "atc_concept_id","meddra_concept_id",
             "nichd","sex","reporter_qualification",
             "receive_date",
             names(atc_covariates)[2:ncol(atc_covariates)],
             "polypharmacy"),with=F] %>% 
  unique()
tmp$receive_date <- 
  as.character(tmp$receive_date)

er_tables[["ade_raw"]] <- tmp

# ATC_RAW_MAP -----------------------------------------------------------------

tmp <- er_tables[["drug"]] %>% 
  .[,.(atc1_concept_name,raw_code = paste0("X",atc1_concept_code))] %>% 
  .[order(raw_code,decreasing = F)] %>% 
  unique()

er_tables[["atc_raw_map"]] <- tmp

# ADE_NULL ----------------------------------------------------------------

tmp <- 
null_dts[,.(null_99 = quantile(gam_score,c(0.99))),nichd] %>% 
  .[,.(nichd = factor(nichd,levels=stages),null_99)] %>% 
  .[order(nichd)]

er_tables[["ade_null"]] <- tmp


# ADE_NULL_DISTRIBUTION ---------------------------------------------------

tmp <- null_dts[,
         .(
           nichd,
           gam_score,
           ade=paste0(atc_concept_id,"_",meddra_concept_id)
           )
         ]

er_tables[["ade_null_distribution"]] <- tmp

# ADE ---------------------------------------------------------------------

tmp <- X[
  ade %in% primary_keys
] %>% 
  merge(
    dts[,.(ade,atc_concept_id,meddra_concept_id)] %>% 
      unique(),
    by="ade"
  ) %>% .[,
          .(
            ade,
            atc_concept_id,
            meddra_concept_id,
            cluster_id = cluster,
            gt_null_statistic = ade %in% sig_ades,
            gt_null_99 = ade %in% sig_null_ades
          )
  ] %>% 
  merge(
    dts[
      database=="covariate_adjusted",
      .SD[which.max(gam_score)],
      ade
    ] %>% 
      .[,.(ade,max_score_nichd = nichd)],
    by="ade"
  )

tmp$cluster_name <- 
  sapply(tmp$cluster,function(x){cluster_names[[x]]})

tmp <- 
tmp %>% 
  merge(
    raw_data[,.(safetyreportid,ade)] %>% unique() %>% .[,.(ade_nreports = .N),ade],
    by="ade"
  )

er_tables[["ade"]] <- tmp

# ADE_NICHD ---------------------------------------------------------------

sub <- dts[database=="covariate_adjusted" &
               ade %in% primary_keys]

norm_ades <- 
    normalize_data(sub,score="gam_score",groups=c())

tmp <- 
  norm_ades[,
            .(
              ade,
              nichd,
              gam_score,
              norm
            )] %>% 
  merge(
    dts[database=="covariate_adjusted",
        .(ade,atc_concept_id,meddra_concept_id,
          nichd,gam_score_se,
          gam_score_90mse,gam_score_90pse,D,E,DE)] %>% 
      unique(),
    by=c("ade","nichd")
  )

er_tables[["ade_nichd"]] <-
  merge(
    tmp,
    merge(
      merge(
        tmp[,.(meddra_concept_id,atc_concept_id)] %>% unique(),
        er_tables[["drug"]][,.(atc_concept_id,atc_concept_name)],
        by="atc_concept_id",
        allow.cartesian=T
      ),
      er_tables[['event']][,.(meddra_concept_id,meddra_concept_name_1)],
      by="meddra_concept_id",
      allow.cartesian=T
    ) %>% 
      .[,.(atc_concept_id,
           meddra_concept_id,
           ade_name = paste0(atc_concept_name," and ",meddra_concept_name_1))
      ],
    by=c("atc_concept_id","meddra_concept_id"),
    allow.cartesian=T
  ) %>% 
  unique()

# ADE_NICHD_ENRICHMENT ----------------------------------------------------

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

er_tables[["ade_nichd_enrichment"]] <- 
  category_test_dts[
    a>0 &
    category %in% c("atc1","atc2","atc3","atc4","atc5",
                    "pt","hlt","hlgt","soc",
                    "soc_atc1","soc_atc2","soc_atc3","soc_atc4","soc_atc5",
                    "hlgt_atc1","hlgt_atc2","hlgt_atc3","hlgt_atc4","hlgt_atc5",
                    "hlt_atc1","hlt_atc2","hlt_atc3","hlt_atc4","hlt_atc5",
                    "pt_atc1","pt_atc2","pt_atc3","pt_atc4"),
    .(category,
      atc_concept_name,meddra_concept_name,nichd,
      atc_concept_class_id,meddra_concept_class_id,
      a,b,c,d,lwr,odds_ratio,upr,pvalue,fdr = fdr_all)
    ]

# DRUG_GENE ---------------------------------------------------------------

tmp <- 
  bind_rows(
    fread(paste0(data_dir,"drug_carriers_actions_rvest.csv")) %>% 
      .[action!="",
        .(id = carrier_id,
          type="carrier",
          action,
          uniprot_id = UNIPROT,
          entrez_id = ENTREZID,
          gene_symbol = SYMBOL,
          drugbank_id = parent_key)
      ] %>% 
      unique(),
    fread(paste0(data_dir,"drug_transporters_actions_rvest.csv")) %>% 
      .[action!="",
        .(id = transporter_id,
          type="transporter",
          action,
          uniprot_id = UNIPROT,
          entrez_id = ENTREZID,
          gene_symbol = SYMBOL,
          drugbank_id = parent_key)
      ] %>% 
      unique(),
    fread(paste0(data_dir,"drug_enzymes_actions_rvest.csv")) %>% 
      .[action!="",
        .(id = enzyme_id,
          type="enzyme",
          action,
          uniprot_id = UNIPROT,
          entrez_id = ENTREZID,
          gene_symbol = SYMBOL,
          drugbank_id = parent_key)
      ] %>% 
      unique(),
    fread(paste0(data_dir,"drug_targets_actions_rvest.csv")) %>% 
      .[action!="",
        .(id = target_id,
          type="target",
          action,
          uniprot_id = UNIPROT,
          entrez_id = ENTREZID,
          gene_symbol = SYMBOL,
          drugbank_id = parent_key)
      ] %>% 
      unique()
  ) %>% 
  merge(
    drugbank_atc[,.(drugbank_id = `drugbank-id`,atc_concept_id)],
    allow.cartesian = T
  ) %>% 
  unique()

er_tables[["drug_gene"]] <- tmp

# SIDER -------------------------------------------------------------------

tmp <- 
  sider

er_tables[["sider"]] <- tmp

# RYAN --------------------------------------------------------------------

merged <-  
  fread(paste0(data_dir,"database_generation_ryan_et_al.csv"))

er_tables[["ryan"]] <- merged

# GRiP --------------------------------------------------------------------

pediatric_reference_with_negatives <- 
  fread("https://raw.githubusercontent.com/ngiangre/GRiP_pediatric_ADE-reference_set/master/data/GRiP_Pediatric_ADE_with_generated_negatives.csv")

pediatric_reference_with_negatives$ade <- 
  pediatric_reference_with_negatives[,paste0(ATC_concept_id,"_",MedDRA_concept_id)]

pediatric_reference_with_negatives$ade_name <- 
  pediatric_reference_with_negatives[,paste0(ATC_concept_name," and ",MedDRA_concept_name)]

er_tables[["grip"]] <- pediatric_reference_with_negatives

# GENE --------------------------------------------------------------------

tmp <- 
  fread(
    paste0(
      data_dir,
      "database_generation_gene_expression_evaluation_stage_association_to_residuals_samples_significance.csv"
      )
    ) %>% 
  .[,.(gene_symbol = gene,probe)] %>% 
  unique()

er_tables[["gene"]] <- 
  tmp

# GENE_EXPRESSION ---------------------------------------------------------

tmp <- 
  fread(
    paste0(data_dir,
           "database_generation_gene_expression_evaluation_stage_association_to_residuals_samples_significance.csv")) %>% 
  .[,
    .(sample,nichd,probe,gene_symbol = gene,actual,prediction,residual,fdr,f_statistic,f_pvalue)]

er_tables[["gene_expression"]] <- 
  tmp

# CYP_GENE_EXPRESSION_SUBSTRATE_RISK_INFORMATION ---------------------------------------------------

all_gene_dt <- 
  fread(paste0(
    data_dir,
    "database_generation_gene_expression_evaluation_cyp_expression_soc_risk_information_summaries.csv")) %>% 
  rename(
    gene_symbol = gene
  )

er_tables[["cyp_gene_expression_substrate_risk_information"]] <- 
  all_gene_dt

# DICTIONARY --------------------------------------------------------------

tmp <- NULL
for(i in 1:length(er_tables)){
  tablename <- names(er_tables)[i]
  cols_ <- er_tables[[i]] %>% colnames()
  tmp <- 
    bind_rows(
      tmp,
      lapply(cols_,function(x){data.table(table = tablename,field = x,description="",type="character")}) %>% bind_rows()
    )
}

##ade_raw##
tmp[table=="ade_raw" & field=="ade","description"] = 
  "Primary key. This is the unique identifier of an adverse drug event (drug-event). It is a combination of the atc_concept_id and the meddra_concept_id."
tmp[table=="ade_raw" & field=="safetyreportid","description"] = 
  "The unique identifier for the report."
tmp[table=="ade_raw" & field=="atc_concept_id","description"] = 
  "The ATC 5th level OMOP concept identifier."
tmp[table=="ade_raw" & field=="meddra_concept_id","description"] = 
  "The MedDRA preferred term OMOP concept identifier. In the event table, this would be equivalent in 'meddra_concept_id_1'."
tmp[table=="ade_raw" & field=="safetyreportid","description"] = 
  "The unique identifier for the report."
tmp[table=="ade_raw" & field=="nichd","description"] = 
  "This is the NICHD-defined child development stage. Defined in https://doi.org/10.1542/peds.2012-0055I."
tmp[table=="ade_raw" & field=="sex","description"] = 
  "The reported sex."
tmp[table=="ade_raw" & field=="polypharmacy","description"] = 
  "The number of drugs reported."
tmp[table=="ade_raw" & field=="reporter_qualification","description"] = 
  "The type of reporter."
tmp[table=="ade_raw" & field=="receive_date","description"] = 
  "The date the report was first submitted."
tmp[table=="ade_raw" & field=="XA","description"] =
  "GAM covariate name for the ATC 1st level concept name 'ALIMENTARY TRACT AND METABOLISM'"
tmp[table=="ade_raw" & field=="XB","description"] =
  "GAM covariate name for the ATC 1st level concept name 'BLOOD AND BLOOD FORMING ORGANS'"
tmp[table=="ade_raw" & field=="XC","description"] =
  "GAM covariate name for the ATC 1st level concept name 'CARDIOVASCULAR SYSTEM'"
tmp[table=="ade_raw" & field=="XJ","description"] =
  "GAM covariate name for the ATC 1st level concept name 'ANTIINFECTIVES FOR SYSTEMIC USE'"
tmp[table=="ade_raw" & field=="XL","description"] =
  "GAM covariate name for the ATC 1st level concept name 'ANTINEOPLASTIC AND IMMUNOMODULATING AGENTS'"
tmp[table=="ade_raw" & field=="XD","description"] =
  "GAM covariate name for the ATC 1st level concept name 'DERMATOLOGICALS'"
tmp[table=="ade_raw" & field=="XG","description"] =
  "GAM covariate name for the ATC 1st level concept name 'GENITO URINARY SYSTEM AND SEX HORMONES'"
tmp[table=="ade_raw" & field=="XH","description"] =
  "GAM covariate name for the ATC 1st level concept name 'SYSTEMIC HORMONAL PREPARATIONS, EXCL. SEX HORMONES AND INSULINS'"
tmp[table=="ade_raw" & field=="XM","description"] =
  "GAM covariate name for the ATC 1st level concept name 'MUSCULO-SKELETAL SYSTEM'"
tmp[table=="ade_raw" & field=="XR","description"] =
  "GAM covariate name for the ATC 1st level concept name 'RESPIRATORY SYSTEM'"
tmp[table=="ade_raw" & field=="XS","description"] =
  "GAM covariate name for the ATC 1st level concept name 'SENSORY ORGANS'"
tmp[table=="ade_raw" & field=="XN","description"] =
  "GAM covariate name for the ATC 1st level concept name 'NERVOUS SYSTEM'"
tmp[table=="ade_raw" & field=="XP","description"] =
  "GAM covariate name for the ATC 1st level concept name 'ANTIPARASITIC PRODUCTS, INSECTICIDES AND REPELLENTS'"
tmp[table=="ade_raw" & field=="XV","description"] =
  "GAM covariate name for the ATC 1st level concept name 'VARIOUS'"
tmp[table=="ade_raw" & field=="age","description"] =
  "The age (normalized to year units) value given for the report's subject"

##atc_raw_map##
tmp[table=="atc_raw_map" & field=="atc1_concept_name","description"] =
  "The ATC 1st level OMOP concept name."
tmp[table=="atc_raw_map" & field=="raw_code","description"] =
  "GAM covariate name for the ATC 1st level concept names - seen in the ade_raw table"

##ade##

tmp[table=="ade" & field=="ade","description"] = 
  "Primary key. This is the unique identifier of an adverse drug event (drug-event). It is a combination of the atc_concept_id and the meddra_concept_id."
tmp[table=="ade" & field=="atc_concept_id","description"] = 
  "The ATC 5th level OMOP concept identifier."
tmp[table=="ade" & field=="meddra_concept_id","description"] = 
  "The MedDRA preferred term OMOP concept identifier. In the event table, this would be equivalent in 'meddra_concept_id_1'."
tmp[table=="ade" & field=="cluster_id","description"] = 
  "The identifier for the cluster group assigned to a drug-event by our data-driven clustering approach. See the manuscript's methods for details."
tmp[table=="ade" & field=="cluster_name","description"] = 
  "The dynamics name given to the identfier of a cluster group. This is descriptive of the risk trend across stages, from birth through adolescence."
tmp[table=="ade" & field=="gt_null_statistic","description"] = 
  "The boolean value indicating whether at least one stage's score was greater than nominal significance (the 90 percent confidence interval was above 0)."
tmp[table=="ade" & field=="gt_null_99","description"] = 
  "The boolean value indicating whether at least one stage's score was greater than significance by the null model, as referenced in the paper (the score was greater than the 99th percentile of the null distribution of randomly co-reported drugs and events)."
tmp[table=="ade" & field=="max_score_nichd","description"] = 
  "The child development stage that had the highest risk score for the drug-event."
tmp[table=="ade" & field=="ade_nreports","description"] = 
  "The number of reports of the drug and event co-occurring"

##ade_nichd##

tmp[table=="ade_nichd" & field=="ade","description"] = 
  "Primary key. This is the unique identifier of an adverse drug event (drug-event). It is a combination of the atc_concept_id and the meddra_concept_id."
tmp[table=="ade_nichd" & field=="atc_concept_id","description"] = 
  "The ATC 5th level OMOP concept identifier."
tmp[table=="ade_nichd" & field=="meddra_concept_id","description"] = 
  "The MedDRA preferred term OMOP concept identifier."
tmp[table=="ade_nichd" & field=="ade_name","description"] = 
  "The named identifier of an adverse drug event (drug-event). It is a combination of the atc_concept_name and the meddra_concept_name."
tmp[table=="ade_nichd" & field=="nichd","description"] = 
  "This is the NICHD-defined child development stage. Defined in https://doi.org/10.1542/peds.2012-0055I."
tmp[table=="ade_nichd" & field=="gam_score","description"] = 
  "The risk coefficient from a drug-event GAM, given to each nichd stage. It is the log odds risk of event occurrence given the data as specified in the manuscript."
tmp[table=="ade_nichd" & field=="norm","description"] = 
  "The normalized risk coefficient, between 0 and 1, across stages for a drug-event. This preserves the risk trend but constrains the range of the risk scores between 0 and 1."
tmp[table=="ade_nichd" & field=="gam_score_se","description"] = 
  "The standard deviation of the risk coefficient."
tmp[table=="ade_nichd" & field=="gam_score_90mse","description"] = 
  "The 90 percent lower bounded risk score using the formula gam_score - (1.645*gam_score_se)."
tmp[table=="ade_nichd" & field=="gam_score_90pse","description"] = 
  "The 90 percent upper bounded risk score using the formula gam_score + (1.645*gam_score_se)."
tmp[table=="ade_nichd" & field=="D","description"] = 
  "The number of reports of the drug at the child development stage"
tmp[table=="ade_nichd" & field=="E","description"] = 
  "The number of reports of the event at the child development stage"
tmp[table=="ade_nichd" & field=="DE","description"] = 
  "The number of reports of the drug & event at the child development stage"

##ade_nichd_enrichment##

tmp[table=="ade_nichd_enrichment" & field=="category","description"] = 
  "The category on enrichment. Either a MedDRA adverse event class, ATC drug class, or a combination of ATC and MedDRA classes. These categories are included in the manuscript results associated to this database."
tmp[table=="ade_nichd_enrichment" & field=="atc_concept_name","description"] = 
  "The ATC concept identifier."
tmp[table=="ade_nichd_enrichment" & field=="meddra_concept_name","description"] = 
  "The MedDRA concept identifier."
tmp[table=="ade_nichd_enrichment" & field=="nichd","description"] = 
  "This is the NICHD-defined child development stage. Defined in https://doi.org/10.1542/peds.2012-0055I."
tmp[table=="ade_nichd_enrichment" & field=="atc_concept_class_id","description"] = 
  "The ATC concept class identifier."
tmp[table=="ade_nichd_enrichment" & field=="meddra_concept_class_id","description"] = 
  "The MedDRA concept class identifier."
tmp[table=="ade_nichd_enrichment" & field=="a","description"] = 
  "The number of significant, by the null model, drug-events in both the stage and ATC/MedDRA concept category."
tmp[table=="ade_nichd_enrichment" & field=="b","description"] = 
  "The number of significant, by the null model, drug-events in the stage and not in the ATC/MedDRA concept category."
tmp[table=="ade_nichd_enrichment" & field=="c","description"] = 
  "The number of significant, by the null model, drug-events not in the stage but in the ATC/MedDRA concept category."
tmp[table=="ade_nichd_enrichment" & field=="d","description"] = 
  "The number of significant, by the null model, drug-events not in the stage and not in the ATC/MedDRA concept category."
tmp[table=="ade_nichd_enrichment" & field=="lwr","description"] = 
  "The 95% lower bound of the odds ratio."
tmp[table=="ade_nichd_enrichment" & field=="odds_ratio","description"] = 
  "The odds ratio for the category and stage enrichment."
tmp[table=="ade_nichd_enrichment" & field=="upr","description"] = 
  "The 95% lower bound of the odds ratio."
tmp[table=="ade_nichd_enrichment" & field=="pvalue","description"] = 
  "The p-value from the fisher exact test."
tmp[table=="ade_nichd_enrichment" & field=="fdr","description"] = 
  "The FDR corrected pvalue."

##ade_null##

tmp[table=="ade_null" & field=="ade","description"] = 
  "Primary key. This is the unique identifier of an adverse drug event (drug-event). It is a combination of the atc_concept_id and the meddra_concept_id."
tmp[table=="ade_null" & field=="nichd","description"] = 
  "This is the NICHD-defined child development stage. Defined in https://doi.org/10.1542/peds.2012-0055I."
tmp[table=="ade_null" & field=="null_99","description"] = 
  "The gam risk coefficient value at the 99th percentile of the risk coefficient distribution from random drug and event GAMs. If a risk coefficient at at least one stage for a drug-event GAM is higher than this value, then the drug-event GAM is significant by the null model."

##ade_null_distribution##
tmp[table=="ade_null_distribution" & field=="nichd","description"] = 
  "This is the NICHD-defined child development stage. Defined in https://doi.org/10.1542/peds.2012-0055I."
tmp[table=="ade_null_distribution" & field=="gam_score","description"] = 
  "The gam risk score coefficient from calculating the association between a random drug and event."
tmp[table=="ade_null_distribution" & field=="ade","description"] = 
  "The random drug and event concept identifier pair. These areATC 5th level and MedDRA PT OMOP concept identifiers that were randomly paired together."

##drug##

tmp[table=="drug" & field=="atc_concept_name","description"] = 
  "The ATC 5th level OMOP concept name. In the ade_nichd_enrichment table, this ATC concept is from any level in the hierarchy."
tmp[table=="drug" & field=="atc_concept_id","description"] = 
  "The ATC 5th level OMOP concept identifier."
tmp[table=="drug" & field=="atc_concept_code","description"] = 
  "The ATC 5th level OMOP concept code."
tmp[table=="drug" & field=="atc4_concept_code","description"] = 
  "The ATC 4th level OMOP concept code."
tmp[table=="drug" & field=="atc3_concept_code","description"] = 
  "The ATC 3rd level OMOP concept code."
tmp[table=="drug" & field=="atc2_concept_code","description"] = 
  "The ATC 2nd level OMOP concept code."
tmp[table=="drug" & field=="atc1_concept_code","description"] = 
  "The ATC 1st level OMOP concept code."
tmp[table=="drug" & field=="atc4_concept_name","description"] = 
  "The ATC 4th level OMOP concept name."
tmp[table=="drug" & field=="atc3_concept_name","description"] = 
  "The ATC 3rd level OMOP concept name."
tmp[table=="drug" & field=="atc2_concept_name","description"] = 
  "The ATC 2nd level OMOP concept name."
tmp[table=="drug" & field=="atc1_concept_name","description"] = 
  "The ATC 1st level OMOP concept name."
tmp[table=="drug" & field=="ndrugreports","description"] = 
  "The number of reports of the drug in Pediatric FAERS."

##drug_gene##

tmp[table=="drug_gene" & field=="drugbank_id","description"] = 
  "The identifier for the drug on drugbank's website."
tmp[table=="drug_gene" & field=="id","description"] = 
  "The identifier of the gene/protein on drugbank's website."
tmp[table=="drug_gene" & field=="type","description"] = 
  "The action of the drug onto the gene/protein. This was scraped off of drugbank webpages - see R script associated with the manuscript for scraping details."
tmp[table=="drug_gene" & field=="action","description"] = 
  "The action value."
tmp[table=="drug_gene" & field=="uniprot_id","description"] = 
  "The protein identifier on uniprot's website."
tmp[table=="drug_gene" & field=="entrez_id","description"] = 
  "The gene identifier from joining the uniprot identifier to the entrez identifer from the microarray platform database package within Bioconductor."
tmp[table=="drug_gene" & field=="gene_symbol","description"] = 
  "The gene symbol identifier from joining the uniprot identifier to the entrez identifer from the microarray platform database package within Bioconductor."
tmp[table=="drug_gene" & field=="atc_concept_id","description"] = 
  "The ATC 5th level OMOP concept identifier."

##event##

tmp[table=="event" & field=="meddra_concept_name_4","description"] = 
  "The MedDRA system organ class concept name."
tmp[table=="event" & field=="meddra_concept_id","description"] = 
  "The MedDRA preferred term OMOP concept identifier."
tmp[table=="event" & field=="meddra_concept_id_2","description"] = 
  "The MedDRA higher level concept identifier."
tmp[table=="event" & field=="meddra_concept_id_3","description"] = 
  "The MedDRA higher level greater term concept identifier."
tmp[table=="event" & field=="meddra_concept_id_4","description"] = 
  "The MedDRA system organ class concept identifier."
tmp[table=="event" & field=="neventreports","description"] = 
  "The number of adverse event reports in Pediatric FAERS."
tmp[table=="event" & field=="meddra_concept_class_id_1","description"] = 
  "The MedDRA preferred term concept class identifier."
tmp[table=="event" & field=="meddra_concept_class_id_2","description"] = 
  "The MedDRA higher level concept class identifier."
tmp[table=="event" & field=="meddra_concept_class_id_3","description"] = 
  "The MedDRA higher level greater term concept class identifier."
tmp[table=="event" & field=="meddra_concept_class_id_4","description"] = 
  "The MedDRA system organ class concept class identifier."
tmp[table=="event" & field=="meddra_concept_code_1","description"] = 
  "The MedDRA preferred term concept code identifier."
tmp[table=="event" & field=="meddra_concept_code_2","description"] = 
  "The MedDRA higher level concept code identifier."
tmp[table=="event" & field=="meddra_concept_code_3","description"] = 
  "The MedDRA higher level greater term concept code identifier."
tmp[table=="event" & field=="meddra_concept_code_4","description"] = 
  "The MedDRA system organ class concept code identifier."
tmp[table=="event" & field=="meddra_concept_name_2","description"] = 
  "The MedDRA higher level concept name."
tmp[table=="event" & field=="meddra_concept_name_3","description"] = 
  "The MedDRA higher level greater term concept name."
tmp[table=="event" & field=="meddra_concept_name_1","description"] = 
  "The MedDRA preferred term concept name. Same as 'meddra_concept_name'"
tmp[table=="event" & field=="relationship_id_12","description"] = 
  "The relationship identifier between columns *1 and *2; should be 'Is a' denoting 1-to-1 mapping."
tmp[table=="event" & field=="relationship_id_23","description"] = 
  "The relationship identifier between columns *2 and *3; should be 'Is a' denoting 1-to-1 mapping."
tmp[table=="event" & field=="relationship_id_34","description"] = 
  "The relationship identifier between columns *3 and *4; should be 'Is a' denoting 1-to-1 mapping."
tmp[table=="event" & field=="soc_category","description"] = 
  "The customized category to represent meddra_concept_name_4 events more broadly as used in the manuscript. Developed in consultation with https://admin.new.meddra.org/sites/default/files/guidance/file/intguide_21_0_english.pdf."
tmp[table=="event" & field=="pediatric_adverse_event","description"] = 
  "Whether this event concept (meddra concept id) was defined by MedDRA 19th edition vocabulary as a pediatric-specific adverse event. One (1) indicates yes and zero (0) indicates no. The list of events were curated from this site: https://www.meddra.org/paediatric-and-gender-adverse-event-term-lists."

##sider##

tmp[table=="sider" & field=="ade","description"] = 
  "Primary key. This is the unique identifier of an adverse drug event (drug-event). It is a combination of the atc_concept_id and the meddra_concept_id."
tmp[table=="sider" & field=="stitch_id","description"] = 
  "The STITCH (http://stitch.embl.de/) concept identifier for the ATC 5th level drug."
tmp[table=="sider" & field=="meddra_concept_name","description"] = 
  "The MedDRA preferred term concept name."
tmp[table=="sider" & field=="medgen_id","description"] = 
  "The identifier of the adverse event concept on NCBI's human medical genetiics database"
tmp[table=="sider" & field=="meddra_concept_id","description"] = 
  "The MedDRA preferred term OMOP concept identifier."
tmp[table=="sider" & field=="soc","description"] = 
  "The MedDRA system organ class concept name. Abbreviated to SOC"
tmp[table=="sider" & field=="atc_concept_id","description"] = 
  "The ATC 5th level OMOP concept identifier."
tmp[table=="sider" & field=="atc_concept_name","description"] = 
  "The ATC 5th level concept name."
tmp[table=="sider" & field=="atc_concept_code","description"] = 
  "The ATC 5th level concept code."
tmp[table=="sider" & field=="ade_name","description"] = 
  "The named identifier of an adverse drug event (drug-event). It is a combination of the atc_concept_name and the meddra_concept_name."

##gene##

tmp[table=="gene" & field=="gene_symbol","description"] = 
  "The gene symbol identifier from joining the uniprot identifier to the entrez identifer from the microarray platform database package within Bioconductor."
tmp[table=="gene" & field=="probe","description"] = 
  "The probe identifier from the affymetrix microarrays used in the manuscript."

##gene_expression##

tmp[table=="gene_expression" & field=="sample","description"] = 
  "The GEO sample identifier used in the GSE datasets."
tmp[table=="gene_expression" & field=="nichd","description"] = 
  "This is the NICHD-defined child development stage. Defined in https://doi.org/10.1542/peds.2012-0055I."
tmp[table=="gene_expression" & field=="probe","description"] = 
  "The probe identifier on the affymetrix gene chip."
tmp[table=="gene_expression" & field=="gene_symbol","description"] = 
  "The gene symbol identifier from joining the uniprot identifier to the entrez identifer from the microarray platform database package within Bioconductor."
tmp[table=="gene_expression" & field=="actual","description"] = 
  "The sample value from the stage-association GLM. See the manuscript for details."
tmp[table=="gene_expression" & field=="prediction","description"] = 
  "The sample predicted value from the stage-association GLM. See the manuscript for details."
tmp[table=="gene_expression" & field=="residual","description"] = 
  "The sample residual (actual - predicted) value from the stage-association GLM. See the manuscript for details."
tmp[table=="gene_expression" & field=="fdr","description"] = 
  "The F test FDR corrected pvalue."
tmp[table=="gene_expression" & field=="f_statistic","description"] = 
  "The F test, as summarized from the glm, statistic."
tmp[table=="gene_expression" & field=="f_pvalue","description"] = 
  "The F test, as summarized from the glm, pvalue."

##cyp_gene_expression_substrate_risk_information##

tmp[table=="cyp_gene_expression_substrate_risk_information" & field=="gene_symbol","description"] = 
  "The gene symbol identifier from joining the uniprot identifier to the entrez identifer from the microarray platform database package within Bioconductor."
tmp[table=="cyp_gene_expression_substrate_risk_information" & field=="type","description"] = 
  "The pharmacological action this gene has on drugs."
tmp[table=="cyp_gene_expression_substrate_risk_information" & field=="soc","description"] = 
  "The MedDRA system organ class concept name. These concepts are a subset of what is represented in the event table."
tmp[table=="cyp_gene_expression_substrate_risk_information" & field=="auroc","description"] = 
  "The normalized W statistic from the Mann Whitney test evaluating the rank difference between the mutual information distributions of drug substrates versus drug non-substrates. See the manuscript methods for details."
tmp[table=="cyp_gene_expression_substrate_risk_information" & field=="wt_pvalue","description"] = 
  "The Mann Whitney p-value from evaluating the rank difference between the mutual information distributions of drug substrates versus drug non-substrates. See the manuscript methods for details."
tmp[table=="cyp_gene_expression_substrate_risk_information" & field=="ttest_statistic","description"] = 
  "The Student's t-test statistic from evaluating the average difference between the mutual information distributions of drug substrates versus drug non-substrates. See the manuscript methods for details."
tmp[table=="cyp_gene_expression_substrate_risk_information" & field=="ttest_pvalue","description"] = 
  "The Student's t-test p-value from evaluating the average difference between the mutual information distributions of drug substrates versus drug non-substrates. See the manuscript methods for details."

tmp <- 
tmp %>% 
  bind_rows(
    data.table(table=NA,field=NA,type=NA,description="The ade_nichd_covariate_summary table fields are the covariates and their categories in constructing the drug-event GAMs. These can be used for describing the report covariates for a drug-event at a particular stage. This covariate summary does not include the controls that made up each drug-event GAM, only the cases.")
  )

fieldtypes <- c("int","tinyint","bigint","float","double","date","character","varchar","text")

tmp[grepl('concept_id',field),"type"] = "int"
tmp[grepl('polypharmacy',field),"type"] = "int"
tmp[grepl('ndrugreports',field),"type"] = "int"
tmp[grepl('neventreports',field),"type"] = "int"
tmp[grepl('polypharmacy',field),"type"] = "int"
tmp[grepl('date',field),"type"] = "date"
tmp[grepl('score',field,perl = T),"type"] = "float"
tmp[grepl('norm',field,perl = T),"type"] = "float"
tmp[grepl('null_99',field),"type"] = "float"
tmp[grepl('pvalue',field),"type"] = "float"
tmp[grepl('statistic',field),"type"] = "float"
tmp[grepl('auroc',field),"type"] = "float"
tmp[grepl('residual',field),"type"] = "float"
tmp[grepl('actual',field),"type"] = "float"
tmp[grepl('prediction',field),"type"] = "float"
tmp[grepl('fdr',field),"type"] = "float"
tmp[grepl('odds',field),"type"] = "float"
tmp[grepl('X[A-Z]',field,perl = T),"type"] = "float"
tmp[grepl('^D$',field),"type"] = "int"
tmp[grepl('^E$',field),"type"] = "int"
tmp[grepl('^DE$',field,perl = T),"type"] = "int"
tmp[grepl('^a$',field,perl = T),"type"] = "int"
tmp[grepl('^b$',field,perl = T),"type"] = "int"
tmp[grepl('^c$',field,perl = T),"type"] = "int"
tmp[grepl('^d$',field,perl = T),"type"] = "int"
tmp[grepl('lwr',field),"type"] = "float"
tmp[grepl('odds ratio',field),"type"] = "float"
tmp[grepl('upr',field),"type"] = "float"
tmp[grepl('f_pvalue',field),"type"] = "float"
tmp[grepl('f_statistic',field),"type"] = "float"
tmp[grepl('pediatric_adverse_event',field),"type"] = "int"


er_tables[["dictionary"]] = tmp 


# Output to directory -----------------------------------------------------


for(i in seq_along(er_tables)){
    tmp <- er_tables[[i]]
    tmp_name <- names(er_tables)[i]
    tmp %>% 
        #arrow::write_parquet(paste0(out_dir,tmp_name,".parquet"))
      fwrite(paste0(out_dir,tmp_name,".csv.gz"))
}

# Create sqlite database --------------------------------------------------

db <- 
  RSQLite::dbConnect(
    RSQLite::SQLite(), 
    dbname = "../../data/database_generation_er_tables.sqlite"
    )
for(i in seq_along(er_tables)){
  tmp <- er_tables[[i]]
  tmp_name <- names(er_tables)[i]
  fieldtypes <- er_tables[["dictionary"]][table==tmp_name,type]
  names(fieldtypes) <- er_tables[["dictionary"]][table==tmp_name,field]
  RSQLite::dbWriteTable(db, tmp_name, tmp,
                        overwrite=T,
                        row.names=FALSE,
                        field.types=fieldtypes)
}
RSQLite::dbDisconnect(db)

db <- 
  RSQLite::dbConnect(
    RSQLite::SQLite(), 
    dbname = "../../data/database_generation_er_tables_sample.sqlite"
  )
for(i in seq_along(er_tables)){
  tmp <- er_tables[[i]]
  tmp_name <- names(er_tables)[i]
  if(!(tmp_name %in% c("sider")))next
  fieldtypes <- er_tables[["dictionary"]][table==tmp_name,type]
  names(fieldtypes) <- er_tables[["dictionary"]][table==tmp_name,field]
  RSQLite::dbWriteTable(db, tmp_name, tmp,
                        overwrite=T,
                        row.names=FALSE,
                        field.types=fieldtypes)
}
RSQLite::dbDisconnect(db)

# Create sql database -----------------------------------------------------

system("sqlite3 ../../data/database_generation_er_tables.sqlite .dump > ../../data/database_generation_er_tables.sql")
system("sqlite3 ../../data/database_generation_er_tables_sample.sqlite .dump > ../../data/database_generation_er_tables_sample.sql")
