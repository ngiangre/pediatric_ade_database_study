#' Title: "A database of pediatric drug effects to evaluate ontogenic mechanisms from child growth and development" study
#' 
#' Author details: Nicholas Giangreco
#' 
#' This script generates the data for evaluating enrichment of drug-events
#' in stages when grouped by categories

# Purpose -----------------------------------------------------------------

#' To evaluate stage enrichment of drug-events from categories
#' 


# Setup -------------------------------------------------------------------

pacman::p_load(data.table,tidyverse,doParallel)

data_dir <- "../../data/"
img_dir <- "../../docs/imgs/"
basename <- "database_generation_stage_enrichment_"
seed = 0
set.seed(seed)

registerDoParallel(cores=50)

stages <- 
    c("term_neonatal","infancy",
      "toddler","early_childhood",
      "middle_childhood","early_adolescence",
      "late_adolescence")

# Load functions ----------------------------------------------------------

source("database_generation_functions.R")

# load GAM data -----------------------------------------------------------

source("database_generation_load_GAM_data.R")

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

drugbank_atc <- 
  fread(paste0(data_dir,"compound_drugbank05/drug_atc_codes_rxnorm_joined.csv"))

drugbank_atc_cyp_substrates <- 
    fread(paste0(data_dir,"drugbank_atc_cyp_substrates.csv"))

drugbank_atc_transporter_substrates <- 
  fread(paste0(data_dir,"drugbank_atc_transporters_substrates.csv"))
  

sider <- 
  fread(paste0(data_dir,"database_generation_sider_standardized_data.csv"))

sider_pt <- 
merge(
  sider[,.(meddra_concept_id)] %>% unique(),
  meddra_relationships[,
                       .(meddra_concept_id = meddra_concept_id_1,
                         sider = meddra_concept_name_1,
                         sider_label = "sider_adverse_event")] %>% 
    unique(),
  by="meddra_concept_id"
)
# Stage enrichment - high throughput approach ---------------------------

cat("\n# stage enrichment #\n")

#' Stage and Category enrichments
#' Categories have subcategories
#' Subcategories are assigned to drug-events
#' Intersect/Setdiff to get 2x2 table
#' Test
#' Output in foreach loop

DT <- 
    dts[
        database=="covariate_adjusted" &
          ade %in% sig_null_ades,
        .(ade,atc_concept_id,meddra_concept_id,nichd,gam_score_90mse)
    ] %>% 
  merge(
    null_dts[,.(null_99 = quantile(gam_score,c(0.99))),nichd],
    by="nichd"
  ) %>% 
  .[gam_score_90mse>null_99] %>% 
    merge(
        meddra_relationships[,.(meddra_concept_id = meddra_concept_id_1,
                                pt = meddra_concept_name_1,
                                hlt = meddra_concept_name_2,
                                hlgt = meddra_concept_name_3,
                                soc = meddra_concept_name_4)],
        by="meddra_concept_id",
        allow.cartesian = T,
        all.x=T
    ) %>% 
    merge(
        drugbank_atc_cyp_substrates[,
                                    .(atc_concept_id,cyp_substrate = enzyme_name)
        ] %>% unique(),
        by="atc_concept_id",
        allow.cartesian = T,
        all.x=T
    ) %>% 
  merge(
    drugbank_atc_transporter_substrates[,
                                .(atc_concept_id,transporter = transporter_name)
    ] %>% unique(),
    by="atc_concept_id",
    allow.cartesian = T,
    all.x=T
  ) %>% 
    merge(
        drugbank_atc[,
                     .(atc_concept_id,
                       atc5 = atc_concept_name,
                       atc4 = level_1,
                       atc3 = level_2,atc2 = level_3,
                       atc1 = level_4)
        ] %>% unique(),
        by="atc_concept_id",
        allow.cartesian = T,
        all.x=T
    ) %>% 
    merge(
        pediatric_aes[,
                      .(meddra_concept_id = concept_id,
                        pediatric_ae = concept_name)
        ],
        by="meddra_concept_id",
        allow.cartesian = T,
        all.x=T
    ) %>% 
  merge(
    sider_pt,
    by="meddra_concept_id",
    allow.cartesian = T,
    all.x=T
  )
DT[,sider_atc1 := paste0(sider," --- ",atc1)]
DT[,sider_atc2 := paste0(sider," --- ",atc2)]
DT[,sider_atc3 := paste0(sider," --- ",atc3)]
DT[,sider_atc4 := paste0(sider," --- ",atc4)]
DT[,pt_atc1 := paste0(pt," --- ",atc1)]
DT[,pt_atc2 := paste0(pt," --- ",atc2)]
DT[,pt_atc3 := paste0(pt," --- ",atc3)]
DT[,pt_atc4 := paste0(pt," --- ",atc4)]
DT[,hlt_atc1 := paste0(hlt," --- ",atc1)]
DT[,hlt_atc2 := paste0(hlt," --- ",atc2)]
DT[,hlt_atc3 := paste0(hlt," --- ",atc3)]
DT[,hlt_atc4 := paste0(hlt," --- ",atc4)]
DT[,hlt_atc5 := paste0(hlt," --- ",atc5)]
DT[,hlgt_atc1 := paste0(hlgt," --- ",atc1)]
DT[,hlgt_atc2 := paste0(hlgt," --- ",atc2)]
DT[,hlgt_atc3 := paste0(hlgt," --- ",atc3)]
DT[,hlgt_atc4 := paste0(hlgt," --- ",atc4)]
DT[,hlgt_atc5 := paste0(hlgt," --- ",atc5)]
DT[,soc_atc1 := paste0(soc," --- ",atc1)]
DT[,soc_atc2 := paste0(soc," --- ",atc2)]
DT[,soc_atc3 := paste0(soc," --- ",atc3)]
DT[,soc_atc4 := paste0(soc," --- ",atc4)]
DT[,soc_atc5 := paste0(soc," --- ",atc5)]
DT[,siderlabel_atc1 := paste0(sider_label," --- ",atc1)]
DT[,siderlabel_atc2 := paste0(sider_label," --- ",atc2)]
DT[,siderlabel_atc3 := paste0(sider_label," --- ",atc3)]
DT[,siderlabel_atc4 := paste0(sider_label," --- ",atc4)]
DT[,siderlabel_atc5 := paste0(sider_label," --- ",atc5)]


DT[,siderlabel_cyp := paste0(sider_label," --- ",cyp_substrate)]
DT[,sider_cyp := paste0(sider," --- ",cyp_substrate)]
DT[,pediatricae_cyp := paste0(pediatric_ae," --- ",cyp_substrate)]
DT[,pt_cyp := paste0(pt," --- ",cyp_substrate)]
DT[,hlt_cyp := paste0(hlt," --- ",cyp_substrate)]
DT[,hlgt_cyp := paste0(hlgt," --- ",cyp_substrate)]
DT[,soc_cyp := paste0(soc," --- ",cyp_substrate)]

DT[,siderlabel_transporter := paste0(sider_label," --- ",transporter)]
DT[,sider_transporter := paste0(sider," --- ",transporter)]
DT[,pediatricae_transporter := paste0(pediatric_ae," --- ",transporter)]
DT[,pt_transporter := paste0(pt," --- ",transporter)]
DT[,hlt_transporter := paste0(hlt," --- ",transporter)]
DT[,hlgt_transporter := paste0(hlgt," --- ",transporter)]
DT[,soc_transporter := paste0(soc," --- ",transporter)]

DT[,pediatricae_atc1 := paste0(pediatric_ae," --- ",atc1)]
DT[,pediatricae_atc2 := paste0(pediatric_ae," --- ",atc2)]
DT[,pediatricae_atc3 := paste0(pediatric_ae," --- ",atc3)]
DT[,pediatricae_atc4 := paste0(pediatric_ae," --- ",atc4)]

totaladesdt <- DT[,unique(ade)]

category <- 
  c("soc","hlgt","hlt","pt",
  "atc5","atc4","atc3","atc2","atc1",
  "soc_atc1","soc_atc2","soc_atc3","soc_atc4","soc_atc5",
    "hlgt_atc1","hlgt_atc2","hlgt_atc3","hlgt_atc4","hlgt_atc5",
    "hlt_atc1","hlt_atc2","hlt_atc3","hlt_atc4","hlt_atc5",
    "pt_atc1","pt_atc2","pt_atc3","pt_atc4",
    "pediatricae_atc1","pediatricae_atc2",
    "pediatricae_atc3","pediatricae_atc4",
    "cyp_substrate","sider_cyp","pediatricae_cyp",
    "pt_cyp","hlt_cyp","hlgt_cyp","soc_cyp","siderlabel_cyp",
    "transporter","sider_transporter","pediatricae_transporter",
    "pt_transporter","hlt_transporter","hlgt_transporter",
    "soc_transporter","siderlabel_transporter",
    "sider_atc1","sider_atc2","sider_atc3","sider_atc4")

DT[,
   lapply(.SD,function(x){length(na.omit(unique(x)))}),
   .SDcols=category]

rm(dts)
rm(null_dts)

category_test_dts <- NULL

t0 = Sys.time()
cat("Start:",format(Sys.time(), "%a %b %d %X %Y"),"\n")
for(category_ in category){
    
    grid <- 
        DT[,c("nichd",category_),with=F] %>% 
        unique() %>% 
        .[!grepl("^NA",get(category_)) & !grepl("NA$",get(category_))]
    setorderv(grid,c("nichd",category_),c(1,1))
    
    totalades <- 
      DT[!grepl("^NA",get(category_)) & !grepl("NA$",get(category_)),unique(ade)]
    
    cat("Processing: ",category_," at ",format(Sys.time(), "%a %b %d %X %Y"),"\n")
    
    test_dts <- foreach(i=1:nrow(grid),
                        .combine = "rbind",
                        .errorhandling = "remove") %dopar% {
                            
                            row <- grid[i]
                            stagevar <- row[,nichd]
                            catvar <- row[,..category_] %>% unlist %>% unname
                            catades <- DT[get(category_)==catvar,unique(ade)]
                            stagesigades <- DT[nichd==stagevar,unique(ade)]
                            #         Category
                            #                Subcategory     Not 
                            # dev
                            #         stage       a           b
                            #         not         c           d 
                            a <- intersect(stagesigades,catades) %>% length()
                            b <- setdiff(stagesigades,catades) %>% length()
                            c <- setdiff(catades,stagesigades) %>% length()
                            d <- setdiff(totalades,c(catades,stagesigades)) %>% 
                                length()
                            mat <- matrix(c(a,b,c,d),nrow = 2,byrow = T)
                            test <- fisher.test(mat,conf.level=0.95)
                            dt <- data.table(
                                nichd = stagevar,
                                category = category_,
                                subcategory = catvar,
                                a=a,
                                b=b,
                                c=c,
                                d=d,
                                odds_ratio = test$estimate,
                                lwr = test$conf.int[1],
                                upr = test$conf.int[2],
                                pvalue = test$p.value
                            )
                            dt[,prr:=( a / (a + b) )/( c / (c + d) )]
                            dt[,prr_se:=sqrt( (1 / a) + (1 / c) - (1 / (a + b)) - (1/(c + d)) )]
                            dt[,prr_lwr:=dt[,prr]/exp(1.96*dt[,prr])]
                            dt[,prr_upr:=dt[,prr]*exp(1.96*dt[,prr])]
                            dt
                        }
    test_dts$fdr <- p.adjust(test_dts$pvalue,"fdr")
    category_test_dts <- 
        bind_rows(
            category_test_dts,
            test_dts
        )
    
}
category_test_dts$fdr_all <- p.adjust(category_test_dts$pvalue,"fdr")
t1 = Sys.time()
cat("End:",format(Sys.time(), "%a %b %d %X %Y"),"\n")
cat("\n",round(as.numeric(difftime(t1,t0,units="mins")),2)," minutes\n")

category_test_dts %>% 
    fwrite(paste0(data_dir,basename,"results.csv"))