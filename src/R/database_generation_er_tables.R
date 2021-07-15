#' Title: "A database of pediatric drug effects to evaluate ontogenic mechanisms from child growth and development" study
#' 
#' Author details: Nicholas Giangreco
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

# Load GAM data -----------------------------------------------------------

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

# ADE_NULL ----------------------------------------------------------------

tmp <- 
null_dts[,.(null_99 = quantile(gam_score,c(0.99))),nichd] %>% 
  .[,.(nichd = factor(nichd,levels=stages),null_99)] %>% 
  .[order(nichd)]

er_tables[["ade_null"]] <- tmp

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
        .(ade,ade_name,atc_concept_id,meddra_concept_id,
          nichd,gam_score_se,
          gam_score_90mse,gam_score_90pse,D,E,DE)] %>% 
      unique(),
    by=c("ade","nichd")
  )

er_tables[["ade_nichd"]] <- tmp

# ADE_NICHD_ENRICHMENT ----------------------------------------------------

tmp <- 
  fread(paste0(data_dir,"database_generation_stage_enrichment_ade_class_data.csv"))

er_tables[["ade_nichd_enrichment"]] <- 
  tmp[,.(atc_concept_name,meddra_concept_name,nichd,
         atc_concept_class_id,meddra_concept_class_id,
         a,b,c,d,lwr,odds_ratio,upr,pvalue,fdr = fdr_all)]

# ADE_NICHD_COVARIATE_SUMMARY ---------------------------------------------

sub <- 
  raw_data

tmp <- 
  sub[,
      c("ade","nichd",colnames(atc_covariates)),
      with=F
  ] %>% 
  unique()

colnames(tmp) <- c("ade","nichd",colnames(atc_covariates_full_bin))


ade_nichd_cov_summ <- 
sub[,
    .(ade,nichd,safetyreportid)
] %>% 
  unique() %>% 
  .[,
    .(
      nreports = length(unique(safetyreportid))
    ),
    .(ade,nichd)
  ] %>% 
  .[order(ade)] %>% 
  merge(
    sub[,
        .(ade,nichd,safetyreportid,polypharmacy)
    ] %>% 
      unique() %>%  
      .[,
        .(
          avg_polypharmacy = mean(polypharmacy)
        ),
        .(ade,nichd)
      ] %>% 
      .[order(ade)]
    ,
    by=c("ade","nichd")
  ) %>% 
  merge(
    sub[,
        .(ade,nichd,safetyreportid,sex)
    ] %>% 
      unique() %>% 
      .[,
        .(
          avg_male = mean(sex=="Male")
        ),
        .(ade,nichd)
      ] %>% 
      .[order(ade)],
    by=c("ade","nichd")
  ) %>% 
  merge(
    sub[,
        .(ade,nichd,safetyreportid,reporter_qualification)
    ] %>% 
      unique() %>% 
      .[,
        .(nreports = length(unique(safetyreportid))),
        .(ade,nichd,reporter_qualification)
      ] %>% 
      dcast(ade + nichd~ reporter_qualification,fill=0,value.var="nreports") %>% 
      .[order(ade)],
    by=c("ade","nichd")
  ) %>% 
  merge(
    tmp %>% 
      melt(
        id.vars=c("ade","nichd",id_col),
        variable.name="atc1_drug_class"
      ) %>% 
      .[,
        .(nreports = sum(value==1)),
        .(ade,nichd,atc1_drug_class)
      ] %>% 
      dcast(ade + nichd~ atc1_drug_class,fill=0,value.var="nreports") %>% 
      .[order(ade)],
    by=c("ade","nichd")
  ) %>% 
  merge(
    raw_data[,
             .(
               ade,nichd,atc_concept_id,meddra_concept_id,
               atc_concept_name,meddra_concept_name
             )
    ] %>% 
      unique(),
    by=c("ade","nichd")
  )

er_tables[["ade_nichd_covariate_summary"]] <- 
  ade_nichd_cov_summ

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

er_tables[["drug"]] <- tmp

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
  

er_tables[["event"]] <- tmp

# SIDER -------------------------------------------------------------------

tmp <- 
  sider

er_tables[["sider"]] <- tmp

# RYAN --------------------------------------------------------------------


# GRIP --------------------------------------------------------------------


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

# CYP_GENE_EXPRESSION_SUBSTRATE_RISK_INFORMATION ---------------------------------------------------

all_gene_dt <- 
  fread(paste0(
    data_dir,
    "database_generation_gene_expression_evaluation_cyp_expression_soc_risk_information_summaries.csv"))

er_tables[["cyp_gene_expression_substrate_risk_information"]] <- 
  all_gene_dt

# DICTIONARY --------------------------------------------------------------

tmp <- 
data.table(
    ade = "Primary key. This is the unique identifier of an adverse drug event (drug-event). It is a combination of the atc_concept_id and the meddra_concept_id.",
    nichd = "This is the NICHD-defined child development stage. Defined in https://doi.org/10.1542/peds.2012-0055I.",
    
    ##ade##
    atc_concept_id = "The ATC 5th level OMOP concept identifier.",
    meddra_concept_id = "The MedDRA preferred term OMOP concept identifier. In the event table, this would be equivalent in 'meddra_concept_id_1'.",
    cluster_id = "The identifier for the cluster group assigned to a drug-event by our data-driven clustering approach. See the manuscript's methods for details.",
    cluster_name = "The dynamics name given to the identfier of a cluster group. This is descriptive of the risk trend across stages, from birth through adolescence.",
    gt_null_statistic = "The boolean value indicating whether at least one stage's score was greater than nominal significance (the 90 percent confidence interval was above 0).",
    gt_null_99 = "The boolean value indicating whether at least one stage's score was greater than significance by the null model, as referenced in the paper (the score was greater than the 99th percentile of the null distribution of randomly co-reported drugs and events).",
    max_score_nichd = "The child development stage that had the highest risk score for the drug-event.",
    
    ##ade_nichd##
    gam_score = "The risk coefficient from a drug-event GAM, given to each nichd stage. It is the log odds risk of event occurrence given the data as specified in the manuscript.",
    norm = "The normalized risk coefficient, between 0 and 1, across stages for a drug-event. This preserves the risk trend but constrains the range of the risk scores between 0 and 1.",
    gam_score_se = "The standard dfeviation of the risk coefficient.",
    gam_score_90mse = "The 90 percent lower bounded risk score using the formula gam_score - (1.645*gam_score_se).",
    gam_score_90pse = "The 90 percent upper bounded risk score using the formula gam_score + (1.645*gam_score_se).",
    D = "The number of reports of the drug at the child development stage",
    E = "The number of reports of the event at the child development stage",
    DE = "The number of reports of the drug & event at the child development stage",
    
    ##ade_null##
    null_99 = "The gam risk coefficient value at the 99th percentile of the risk coefficient distribution from random drug and event GAMs. If a risk coefficient at at least one stage for a drug-event GAM is higher than this value, then the drug-event GAM is significant by the null model.",
    
    ##drug##
    atc_concept_name = "The ATC 5th level OMOP concept name. In the ade_nichd_enrichment table, this ATC concept is from any level in the hierarchy.",
    atc_concept_code = "The ATC 5th level OMOP concept code.",
    atc4_concept_name = "The ATC 4th level OMOP concept name.",
    atc4_concept_code = "The ATC 4th level OMOP concept code.",
    atc3_concept_name = "The ATC 3rd level OMOP concept name.",
    atc3_concept_code = "The ATC 3rd level OMOP concept code.",
    atc2_concept_name = "The ATC 2nd level OMOP concept name.",
    atc2_concept_code = "The ATC 2nd level OMOP concept code.",
    atc1_concept_name = "The ATC 1st level OMOP concept name.",
    atc1_concept_code = "The ATC 1st level OMOP concept code.",
    
    ##drug_gene##
    drugbank_id = "The identifier for the drug on drugbank's website.",
    id = "The identifier of the gene/protein on drugbank's website.",
    type = "The action of the drug onto the gene/protein. This was scraped off of drugbank webpages - see R script associated with the manuscript for scraping details.",
    action = "The action value.",
    uniprot_id = "The protein identifier on uniprot's website.",
    entrez_id = "The gene identifier from joining the uniprot identifier to the entrez identifer from the microarray platform database package within Bioconductor.",
    gene_symbol = "The gene symbol identifier from joining the uniprot identifier to the entrez identifer from the microarray platform database package within Bioconductor.",
    
    ##gene##
    probe = "The probe identifier from the affymetrix microarrays used in the manuscript.",
    
    ##event##
    meddra_concept_id_1 = "The MedDRA preferred term OMOP concept identifier.",
    meddra_concept_id_2 = "The MedDRA higher level OMOP concept identifier.",
    meddra_concept_id_3 = "The MedDRA higher level greater term OMOP concept identifier.",
    meddra_concept_id_4 = "The MedDRA system organ class OMOP concept identifier.",
    meddra_concept_code_1 = "The MedDRA preferred term concept code identifier.",
    meddra_concept_code_2 = "The MedDRA higher level concept code identifier.",
    meddra_concept_code_3 = "The MedDRA higher level greater term concept code identifier.",
    meddra_concept_code_4 = "The MedDRA system organ class concept code identifier.",
    meddra_concept_name = "The MedDRA preferred term concept name. In the ade_nichd_enrichment table, this MedDRA concept is from any level in the hierarchy.",
    meddra_concept_name_2 = "The MedDRA higher level concept name.",
    meddra_concept_name_3 = "The MedDRA higher level greater term concept name.",
    meddra_concept_name_4 = "The MedDRA system organ class concept name.",
    meddra_concept_class_id_1 = "The MedDRA preferred term concept class identifier.",
    meddra_concept_class_id_2 = "The MedDRA higher level concept class identifier.",
    meddra_concept_class_id_3 = "The MedDRA higher level greater term concept class identifier.",
    meddra_concept_class_id_4 = "The MedDRA system organ class concept class identifier.",
    relationship_id_12 = "The relationship identifier between columns *1 and *2; should be 'Is a' denoting 1-to-1 mapping.",
    relationship_id_23 = "The relationship identifier between columns *2 and *3; should be 'Is a' denoting 1-to-1 mapping.",
    relationship_id_34 = "The relationship identifier between columns *3 and *4; should be 'Is a' denoting 1-to-1 mapping.",
    soc_category = "The customized category to represent meddra_concept_name_4 events more broadly as used in the manuscript. Developed in consultation with https://admin.new.meddra.org/sites/default/files/guidance/file/intguide_21_0_english.pdf.",
    
    ##sider##
    meddra_umls_id = 'The identifier to the UMLS (https://www.nlm.nih.gov/research/umls/index.html) concept identifier for the MedDRA term.',
    stitch_id_flat = "The STITCH (http://stitch.embl.de/) concept identifier for the ATC 5th level drug.",
    method_of_detection = "The method of detection on the drug label for a MedDRA preferred term.",
    meddra_concept_name = "The MedDRA preferred term concept name.",
    ade_name = "The concatenation of the ATC 5th level and MedDRA preferred term concept names.",
    
    ##cyp_gene_expression_substrate_risk_information##
    gene = "The cytochrome P450 gene symbol representing a drug's enzyme. A drug is the substrate for this enzyme so action=='substrate'. These genes are considered in our analysis and are a subset of the available gene symbols in the drug_gene and gene tables.",
    soc = "The MedDRA system organ class concept name. These concepts are a subset of what is represented in the event table.",
    auroc = "The normalized W statistic from the Mann Whitney test evaluating the rank difference between the mutual information distributions of drug substrates versus drug non-substrates. See the manuscript methods for details.",
    wt_pvalue = "The Mann Whitney p-value from evaluating the rank difference between the mutual information distributions of drug substrates versus drug non-substrates. See the manuscript methods for details.",
    ttest_statistic = "The Student's t-test statistic from evaluating the average difference between the mutual information distributions of drug substrates versus drug non-substrates. See the manuscript methods for details.", 
    ttestt_pvalue = "The Student's t-test p-value from evaluating the average difference between the mutual information distributions of drug substrates versus drug non-substrates. See the manuscript methods for details.",
    
    ##ade_nichd_enrichment##
    atc_concept_class_id = "The MedDRA concept class identifier.",
    meddra_concept_class_id = "The MedDRA concept class identifier.",
    a = "The number of significant, by the null model, drug-events in both the stage and ATC/MedDRA concept category.",
    b = "The number of significant, by the null model, drug-events in the stage and not in the ATC/MedDRA concept category.",
    c = "The number of significant, by the null model, drug-events not in the stage but in the ATC/MedDRA concept category.",
    d = "The number of significant, by the null model, drug-events not in the stage and not in the ATC/MedDRA concept category.",
    lwr = "The 95% lower bound of the odds ratio.",
    odds_ratio = "The odds ratio for the category and stage enrichment.",
    upr= "The 95% lower bound of the odds ratio.",
    pvalue = "The p-value from the fisher exact test.",
    fdr = "The FDR corrected pvalue.",
    
    tmp2 = "description"
) %>% melt(id.vars="tmp2",variable.name="field",value.name = "description")

tmp[,tmp2:=NULL] 

tmp$database_notes <- 
  c(
    "The ade_nichd_covariate_summary table fields are the covariates and their categories in constructing the drug-event GAMs. These can be used for describing the report covariates for a drug-event at a particular stage. This covariate summary does not include the controls that made up each drug-event GAM, only the cases.",
    rep("NA",nrow(tmp)-1)
    )

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
  RSQLite::dbWriteTable(db, tmp_name, tmp,overwrite=T)
}
RSQLite::dbDisconnect(db)

# Create sql database -----------------------------------------------------

system("sqlite3 ../../data/database_generation_er_tables.sqlite .dump > ../../data/database_generation_er_tables.sql")
