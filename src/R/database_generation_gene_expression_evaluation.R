#' Title: "A database of pediatric drug effects to evaluate ontogenic mechanisms from child growth and development" study
#' 
#' Author details: Nicholas Giangreco
#' 
#' This script generates the data and plots to evaluate the information between dynamics of drug substrate risks and gene expression across childhood 

# Purpose -----------------------------------------------------------------

#' To evaluate significant temporal expression for gene products and evaluate correlation with ADE risks by drugs that are substrates for those gene products
#' 


# Setup -------------------------------------------------------------------

pacman::p_load(tidyverse,data.table,doParallel)

data_dir <- "../../data/"
img_dir <- "../../docs/imgs/"
basename <- "database_generation_gene_expression_evaluation_"
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


# PROCESSING -- Load GAM data -----------------------------------------------------------

source("database_generation_load_GAM_data.R")

rm(dt)
rm(candidate_ades)
rm(common_ades)
rm(sample_ades)
rm(sig_ades)
rm(narrow_ades)

# PROCESSING -- Load functions ----------------------------------------------------------

source("database_generation_functions.R")

# PROCESSING -- Load sider side effects -------------------------------------------------

sider <- fread(paste0(data_dir,"database_generation_sider_standardized_named_data.csv"))

#sider <-  fread(paste0(data_dir,"database_generation_sider_standardized_named_data.csv"))

dts_sub <- 
  dts[
    ade %in% sig_null_ades & 
      database=="covariate_adjusted",
    .(atc_concept_id,meddra_concept_id,ade,
      nichd,gam_score,gam_score_se,
      gam_score_90mse,gam_score_90pse,
      ade_name,D,E,DE,database)
  ] %>% .[order(nichd)]

rm(dts)

labeled_sig_ades <- 
  merge(
    dts_sub[,c("atc_concept_id","meddra_concept_id","ade"),with=F] %>% unique(),
    sider[,c("atc_concept_id","meddra_concept_id"),with=F] %>% unique(),
    by=c("atc_concept_id","meddra_concept_id")
  ) %>% unique() %>% .[,ade]

# PROCESSING -- Load drugbank data ------------------------------------------------------

drugbank_atc <- 
  fread(paste0(data_dir,"compound_drugbank05/drug_atc_codes_rxnorm_joined.csv"))

drugbank_all <- 
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
      ],
    fread(paste0(data_dir,"drug_transporters_actions_rvest.csv")) %>% 
      .[action!="",
        .(id = transporter_id,
          type="transporter",
          action,
          uniprot_id = UNIPROT,
          entrez_id = ENTREZID,
          gene_symbol = SYMBOL,
          drugbank_id = parent_key)
      ],
    fread(paste0(data_dir,"drug_enzymes_actions_rvest.csv")) %>% 
      .[action!="",
        .(id = enzyme_id,
          type="enzyme",
          action,
          uniprot_id = UNIPROT,
          entrez_id = ENTREZID,
          gene_symbol = SYMBOL,
          drugbank_id = parent_key)
      ],
    fread(paste0(data_dir,"drug_targets_actions_rvest.csv")) %>% 
      .[action!="",
        .(id = target_id,
          type="target",
          action,
          uniprot_id = UNIPROT,
          entrez_id = ENTREZID,
          gene_symbol = SYMBOL,
          drugbank_id = parent_key)
      ]
  ) %>% 
  merge(
    drugbank_atc[,.(drugbank_id = `drugbank-id`,atc_concept_id)],
    allow.cartesian = T
  ) %>% 
  unique()

drugbank_substrates <- drugbank_all[action=="substrate"]

drugbank_substrates[
  atc_concept_id %in% 
    (sapply(labeled_sig_ades,function(x){str_split(x,"_")[[1]][1]}) %>% unique()),
  length(unique(atc_concept_id))
  ]

gene_type_substrates <- 
  drugbank_substrates[
    atc_concept_id %in% 
      (sapply(labeled_sig_ades,function(x){str_split(x,"_")[[1]][1]}) %>% unique()),
    .(N = length(unique(atc_concept_id))),
    .(gene_symbol,type)
  ]

drugbank_substrates[,
                    .(
                      Ndrugbank_substrate_drugs = 
                        length(unique(atc_concept_id))
                    ),
                    .(gene_symbol,type)
] %>% 
  merge(
    drugbank_substrates[
      atc_concept_id %in% 
        (sapply(labeled_sig_ades,function(x){str_split(x,"_")[[1]][1]}) %>% unique()),
      .(
        Nsignificant_substrate_drugs = 
          length(unique(atc_concept_id))
      ),
      .(gene_symbol,type)
    ],
    all.x=T
  ) %>% 
  fwrite(paste0(data_dir,basename,"drugbank_gene_Ndrugs_Nsigdrugs.csv"))

# PROCESSING -- Load gene expression data -----------------------------------------------

gdata <- 
    fread(
      paste0(
        data_dir,
        "database_generation_stevens_et_al_raw_cel_processed_data.csv"
        )
      ) %>% 
  .[gse %in% c("GSE9006","GSE26440","GSE11504","TABM666")]

gdata$nichd <- factor(gdata$nichd,levels=stages)

gdata_stages <- intersect(stages,gdata[,unique(nichd)])

gdata %>% 
  .[,.(N = length(unique(sample))),.(gse,assay)] %>% 
  fwrite(paste0(data_dir,basename,"Nsamples_gse_assays.csv"))

gdata %>% 
  .[,.(N = length(unique(sample))),.(gse,nichd)] %>% 
  .[,.(N,gse,nichd = factor(nichd,levels=stages))] %>% 
  na.omit() %>% 
  dcast(gse ~ nichd,value.var="N",fill = 0) %>% 
  fwrite(paste0(data_dir,basename,"Nsamples_across_stages.csv"))


# PROCESSING -- Load genome-wide PCA --------------------------------------------------

pca_data <- 
  fread(paste0(data_dir,"database_generation_gene_expression_validation_pca_data.csv"))

# PROCESSING -- Generate cyp expression using glm residuals ------------------------------

probes <- gdata[grep("^CYP",SYMBOL),unique(PROBEID)]
resid_dt <- lapply(
  probes,function(x){
    dat <-
      gdata[
        PROBEID==x
      ] %>%
      .[,
        .(sample,years,
          value)
      ] %>% 
      unique() %>% 
      merge(
        pca_data[,
                 .(sample,
                   PC1,PC2,PC3,PC4,PC5,PC6,
                   gse,sex,nichd,nichd_int = as.integer(factor(nichd,levels=stages)))
        ],
        by="sample"
      )
    if(nrow(dat)==0){
      NA
    }else{
      test <- glm(value ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + gse,data=dat)
      dat$residual <- test$fitted.values
      ftest <- anova(glm(residual ~ nichd_int,data=dat),test = "F")
      
      data.table(
        probe = x,
        sample = dat$sample,
        nichd = dat$nichd,
        actual = dat$value,
        prediction = test$fitted.values,
        residual = test$residuals,
        f_statistic = 
          ftest$`F`[
            rownames(ftest)=="nichd_int"
          ],
        f_pvalue = 
          ftest$`Pr(>F)`[
            rownames(ftest)=="nichd_int"
          ]
      )
    }
  }
) %>% bind_rows()

resid_dt_gene <- 
  resid_dt %>% 
  merge(
    gdata[,.(probe = PROBEID,sample,gene = SYMBOL)] %>% unique(),
    by=c("sample","probe"),
    allow.cartesian=T
  )

resid_dt_gene$fdr <- p.adjust(resid_dt_gene$f_pvalue,method="fdr")

resid_dt_gene %>% 
  fwrite(paste0(data_dir,basename,"cyp_stage_association_to_residuals_samples_significance.csv"))

rm(resid_dt)
rm(resid_dt_gene)

# PROCESSING -- Compute relationship between pairwise within-cyp probe expression ----------------------

resid_gene_dt <- 
  fread(paste0(data_dir,basename,"cyp_stage_association_to_residuals_samples_significance.csv"))

genes_ <- resid_gene_dt[fdr<0.1,unique(gene)]

nprobes_dt <- 
  resid_gene_dt[gene %in% genes_,
                .(Nprobes = length(unique(probe))),
                gene] %>% 
  .[order(gene)]

cor_mat_melts <- NULL
for(gene_ in genes_){
  sub <- 
    resid_gene_dt[
      fdr<0.1 & gene==gene_
    ] %>% 
    na.omit() %>% 
    .[,.(probe,sample,residual,nichd)] %>% 
    unique() %>% 
    .[,
      .(m = mean(residual)),
      .(nichd,probe)
    ] %>% 
    dcast(probe~ factor(nichd,levels=stages),value.var="m") %>% 
    transpose()
  
  probes <- sub[1] %>% unlist %>% unname
  res <- Hmisc::rcorr(as.matrix(sub[-c(1)]),type="spearman")
  cor_mat <- res$r
  colnames(cor_mat) <- probes
  rownames(cor_mat) <- probes
  g <- ggcorrplot::ggcorrplot(cor_mat) + ggtitle(paste0(gene_," probe correlations"))
  #ggsave(paste0(img_dir,basename,gene_,"_probe_correlations.png"),g,width=7,height=7)
  cor_mat <- cor_mat %>% data.table()
  cor_mat$PROBEID <- probes
  cor_mat_melt <- 
    cor_mat %>% 
    melt(id.vars="PROBEID",variable.name="PROBEID2") %>%  
    merge(
      data.table(
        id = 1:length(probes),
        PROBEID = probes
      ),
      by="PROBEID"
    ) %>% 
    merge(
      data.table(
        id2 = 1:length(probes),
        PROBEID2 = probes
      ),
      by="PROBEID2"
    ) %>% 
    .[PROBEID!=PROBEID2]
  cor_mat_melt$id12 <- cor_mat_melt[,id*id2]
  cor_mat_melt$gene <- gene_
  cor_mat_melt <- 
    cor_mat_melt %>% 
    .[!duplicated(cor_mat_melt[,paste0(gene,"_",id12)])]
  cor_mat_melts <- 
    bind_rows(
      cor_mat_melts,
      cor_mat_melt
    )
}

cor_mat_melts %>% 
  fwrite(paste0(data_dir,basename,"cyp_gene_probe_correlations.csv"))



# PROCESSING -- Tabulate cyp dataset----------------------------------------------------------------

resid_value_gene_dt <- 
  fread(paste0(data_dir,basename,"cyp_stage_association_to_residuals_samples_significance.csv"))

cor_mat_melts <- 
  fread(paste0(data_dir,basename,"cyp_gene_probe_correlations.csv"))

discordant_probes <- 
  union(
    cor_mat_melts[,.(m = mean(value)),PROBEID][m<0,PROBEID],
    cor_mat_melts[,.(m = mean(value)),PROBEID2][m<0,PROBEID2]
  )

nconcordantprobes_dt <-
  bind_rows(
    cor_mat_melts[
      (!PROBEID %in% discordant_probes &
         !PROBEID2 %in% discordant_probes),
      .(PROBEID,gene)
    ] %>% unique(),
    cor_mat_melts[
      (!PROBEID %in% discordant_probes &
         !PROBEID2 %in% discordant_probes),
      .(PROBEID = PROBEID2,gene)
    ] %>% unique()
  ) %>% unique() %>% 
  .[,.(N_concordant_probes = .N),gene]

concordant_probes <- 
  union(
    cor_mat_melts %>% 
      .[(!PROBEID %in% discordant_probes &
           !PROBEID2 %in% discordant_probes),unique(PROBEID)],
    cor_mat_melts %>% 
      .[(!PROBEID %in% discordant_probes &
           !PROBEID2 %in% discordant_probes),unique(PROBEID2)]
  )

concordant_and_single_probes <- 
  union(
    concordant_probes,
    resid_value_gene_dt[gene %in% 
            (resid_value_gene_dt[fdr<0.1,.(gene,probe)] %>% unique() %>% .[,.N,gene] %>% .[N==1,gene]),
          unique(probe)]
  )

genes_concordant_probes <- 
  cor_mat_melts %>% 
  .[(!PROBEID %in% discordant_probes &
       !PROBEID2 %in% discordant_probes),
    unique(gene)]

genes_concordant_and_single_probes <- 
  union(
    genes_concordant_probes,
    (resid_value_gene_dt[fdr<0.1,.(gene,probe)] %>% unique() %>% .[,.N,gene] %>% .[N==1,gene])
  )

resid_value_gene_dt[,length(unique(gene))]
resid_value_gene_dt[fdr<0.1,length(unique(gene))]

length(genes_concordant_and_single_probes)
length(concordant_and_single_probes) 

length(sig_null_ades)
merge(
  dts_sub[
    ade %in% sig_null_ades
  ] %>% 
    .[,.(atc_concept_id,meddra_concept_id)] %>% unique(),
  sider[,.(atc_concept_id,meddra_concept_id)] %>% unique(),
  by=c("atc_concept_id","meddra_concept_id")
) %>% nrow()
sapply(labeled_sig_ades,function(x){str_split(x,"_")[[1]][1]}) %>% unique() %>% length()
data.table(
  drug = sapply(labeled_sig_ades,function(x){str_split(x,"_")[[1]][1]}),
  ade = labeled_sig_ades
) %>% .[,.N,drug] %>% .[,summary(N)]

intersect(
  drugbank_substrates[,unique(atc_concept_id)],
  sapply(labeled_sig_ades,function(x){str_split(x,"_")[[1]][1]}) %>% unique()
) %>% length()
intersect(
  drugbank_substrates[,unique(atc_concept_id)],
  sapply(labeled_sig_ades,function(x){str_split(x,"_")[[1]][1]}) %>% unique()
) %>% length() / sapply(labeled_sig_ades,function(x){str_split(x,"_")[[1]][1]}) %>% unique() %>% length()
drugbank_substrates[,length(unique(gene_symbol))]

gene_type_soc_ade <-
  dts_sub[
    ade %in% labeled_sig_ades,
    .(atc_concept_id,meddra_concept_id,ade,
      nichd,gam_score,gam_score_se,ade_name)
  ] %>% 
  .[order(nichd)] %>% 
  merge(
    sider[,.(meddra_concept_id,soc)] %>% unique(),
    by="meddra_concept_id",
    allow.cartesian = T,
    all.x=T
  ) %>% 
  merge(
    drugbank_substrates[
      gene_symbol %in% genes_concordant_and_single_probes,
      .(atc_concept_id,type,gene_symbol)
    ] %>% 
      unique(),
    by="atc_concept_id",
    allow.cartesian = T,
    all.x=T
  )

gene_type_soc_ade[!is.na(gene_symbol),length(unique(gene_symbol))]
gene_type_soc_ade[!is.na(gene_symbol),length(unique(gene_symbol)),type]
gene_type_soc_ade[,length(unique(soc))]
gene_type_soc_ade[,length(unique(ade))]
gene_type_soc_ade[,length(unique(atc_concept_id))]

gene_type_soc_ade %>% 
  fwrite(paste0(data_dir,basename,"cyp_gene_type_soc_ade.csv"))

# PROCESSING -- Delete unused variables -------------------------------------------------

rm(gdata)
rm(sig_null_ades)
rm(cor_mat_melts)

# ANALYSIS -- CYP substrate vs. nonsubstrate MI test ---------------------------

resid_value_gene_dt <- 
  fread(paste0(data_dir,basename,"cyp_stage_association_to_residuals_samples_significance.csv"))
resid_value_gene_dt[,length(unique(gene))]

gene_type_soc_ade <- 
  fread(paste0(data_dir,basename,"cyp_gene_type_soc_ade.csv"))
setdiff(gene_type_soc_ade[,unique(ade)],gene_type_soc_ade[gene_symbol=="",unique(ade)]) %>% length()
setdiff(gene_type_soc_ade[,unique(atc_concept_id)],gene_type_soc_ade[gene_symbol=="",unique(atc_concept_id)]) %>% length()

registerDoParallel(cores=50)

iter <- 
  gene_type_soc_ade[!is.na(gene_symbol) & gene_symbol!="",
                    .(
                      Ndrugs = length(unique(atc_concept_id)),
                      Nades = length(unique(ade))
                    ),
                    .(gene_symbol,type)
  ]
iter[,length(unique(gene_symbol))]
resid_value_gene_dt[gene %in% iter[,unique(gene_symbol)] & fdr<0.1,unique(probe)]

Nsamp = 100
zs = rnorm(Nsamp)
all_gene_dt <- NULL
all_full_gene_dt <- NULL
all_gsubs <- NULL
all_overall_risks_comparisons <- NULL
system.time({
  for(i in 1:nrow(iter)){
    
    gene_=iter[i,gene_symbol]
    type_=iter[i,type]
    cat("\nIteration",i,":",gene_,type_,"influencing drug risks\n")
    
    probes <- resid_value_gene_dt[gene==gene_ & fdr<0.1,unique(probe)]
    if(length(probes)==0)next
    
    gsub_obs <- resid_value_gene_dt[probe %in% probes & gene==gene_]
    all_gsubs <- 
      bind_rows(
        all_gsubs,
        gsub_obs
      )
    
    substrates <- 
      gene_type_soc_ade[
        gene_symbol==gene_ &
          type==type_,
        unique(ade)
      ]
    nonsubstrates <- 
      setdiff(
        gene_type_soc_ade[
          gene_symbol!=gene_,
          unique(ade)
        ],
        substrates
      )
    
    
    if(length(substrates)<2)next
    if(length(nonsubstrates)<2)next
    
    substrate_risks <- 
      dts_sub[
        ade %in% substrates
      ]
    substrate_risks$substrate_type <- "substrate drug-events"
    nonsubstrate_risks <- 
      dts_sub[
        ade %in% nonsubstrates
      ]
    nonsubstrate_risks$substrate_type <- "nonsubstrate drug-events"
    overall_risks_comparison <- bind_rows(substrate_risks,nonsubstrate_risks)
    overall_risks_comparison$gene <- gene_
    overall_risks_comparison$type <- type_
    all_overall_risks_comparisons <- 
      bind_rows(
        all_overall_risks_comparisons,
        overall_risks_comparison
      )
    
    substrate_drugs <- 
      sapply(substrates,function(x){str_split(x,"_")[[1]][1]}) %>% unique()
    nonsubstrate_drugs <- 
      sapply(nonsubstrates,function(x){str_split(x,"_")[[1]][1]}) %>% unique()  
    
    expand.grid(ade = substrates,probe=probes) %>% data.table() %>% nrow() %>% cat("\n\t\t",.)
    expand.grid(ade = nonsubstrates,probe=probes) %>% data.table() %>% nrow() %>% cat("\n\t\t",.)
    
    gene_dt <- 
      foreach(boot_=1:length(zs),.combine = "rbind",.errorhandling = "remove") %dopar% {
        z_ = zs[boot_]
        tmp <- data.table()
        
        grid <-
          bind_rows(
            expand.grid(ade = substrates,probe=probes,substrate="substrate") %>% data.table(),
            expand.grid(ade = nonsubstrates,probe=probes,substrate="nonsubstrate") %>% data.table()
          )
        
        tmp <- 
          lapply(1:nrow(grid),function(i){
            ade_ <- grid[i,ade]
            probe_ <- grid[i,probe]
            stype_ <- grid[i,substrate]
            mi_ <- 
              merge(
                dts_sub[ade==ade_,.(ade,nichd,score = (z_*gam_score_se)+gam_score)] %>% unique(),
                gsub_obs[probe==probe_,.(m = mean(residual)),nichd],
                by="nichd"
              ) %>% .[,maigesPack::MI(m,score)]
            cor_ <- 
              merge(
                dts_sub[ade==ade_,.(ade,nichd,score = (z_*gam_score_se)+gam_score)] %>% unique(),
                gsub_obs[probe==probe_,.(m = mean(residual)),nichd],
                by="nichd"
              ) %>% .[,cor(m,score,method="spearman")]
            data.table(ade = ade_,probe = probe_,mi=mi_,cor=cor_,substrate=stype_)
          }) %>% bind_rows()
        
        tmp$gene = gene_
        tmp$type = type_
        tmp$z = z_
        tmp$boot=boot_
        
        tmp
        
      }    
    
    smi = gene_dt[substrate=="substrate",mean(mi),.(ade,probe)][,`V1`]
    nsmi = gene_dt[substrate=="nonsubstrate",mean(mi),.(ade,probe)][,`V1`]
    denom = (length(smi)*length(nsmi))
    wtest <- 
      wilcox.test(smi,nsmi,alternative="g")
    
    smi = gene_dt[substrate=="substrate"][,mi]
    nsmi = gene_dt[substrate=="nonsubstrate"][,mi]
    ttest <- 
      t.test(smi,nsmi,alternative="g")
    
    all_gene_dt <- 
      bind_rows(
        all_gene_dt,
        gene_dt[,
                .(
                  gene = gene_,
                  type = type_,
                  auroc = 
                    (wtest$statistic/
                       denom),
                  wt_pvalue = 
                    wtest$p.value,
                  ttest_statistic = 
                    ttest$statistic,
                  ttest_pvalue = 
                    ttest$p.value
                )
        ]
      )
    
    
    all_full_gene_dt <- 
      bind_rows(
        all_full_gene_dt,
        gene_dt
      )
    
  }
})

all_gene_dt %>% 
  fwrite(paste0(data_dir,basename,"cyp_expression_risk_information_summaries.csv"))

all_full_gene_dt %>% 
  fwrite(paste0(data_dir,basename,"cyp_expression_risk_information.csv"))

all_overall_risks_comparisons %>% 
  fwrite(paste0(data_dir,basename,"cyp_expression_risk_substrate_nonsubstrate_risks.csv"))

all_gsubs %>% 
  fwrite(paste0(data_dir,basename,"cyp_expression_risk_probe_sample_residuals.csv"))

# ANALYSIS -- CYP substrate vs. nonsubstrate MI test for systemic disorders ---------------------------

resid_value_gene_dt <- 
  fread(paste0(data_dir,basename,"cyp_stage_association_to_residuals_samples_significance.csv"))

gene_type_soc_ade <- 
  fread(paste0(data_dir,basename,"cyp_gene_type_soc_ade.csv"))

registerDoParallel(cores=4)

iter <- 
  gene_type_soc_ade[!is.na(gene_symbol),
                    .(
                      Ndrugs = length(unique(atc_concept_id)),
                      Nades = length(unique(ade))
                    ),
                    .(gene_symbol,soc,type)
  ] 

Nsamp = 100
zs = rnorm(Nsamp)
all_gene_dt <- NULL
all_full_gene_dt <- NULL
all_gsubs <- NULL
all_overall_risks_comparisons <- NULL
system.time({
  for(i in 1:nrow(iter)){
    
    #i=4
    gene_=iter[i,gene_symbol]
    type_=iter[i,type]
    soc_=iter[i,soc]
    cat("\nIteration",i,":",gene_,type_,"influencing drug risks for",soc_,"\n")
    probes <- resid_value_gene_dt[gene==gene_ & fdr<0.1,unique(probe)]
    if(length(probes)==0)next
    
    gsub_obs <- resid_value_gene_dt[probe %in% probes & gene==gene_]
    all_gsubs <- 
      bind_rows(
        all_gsubs,
        gsub_obs
      )
    
    substrates <- 
      gene_type_soc_ade[
        gene_symbol==gene_ &
          type==type_ &
          soc==soc_,
        unique(ade)
      ]
    nonsubstrates <- 
      setdiff(
        gene_type_soc_ade[
          gene_symbol!=gene_ &
            soc==soc_,
          unique(ade)
        ],
        substrates
      )
    
    if(length(substrates)<2)next
    if(length(nonsubstrates)<2)next
    
    substrate_risks <- 
      dts_sub[
        ade %in% substrates
      ]
    substrate_risks$substrate_type <- "substrate drug-events"
    nonsubstrate_risks <- 
      dts_sub[
        ade %in% nonsubstrates
      ]
    nonsubstrate_risks$substrate_type <- "nonsubstrate drug-events"
    overall_risks_comparison <- bind_rows(substrate_risks,nonsubstrate_risks)
    overall_risks_comparison$gene <- gene_
    overall_risks_comparison$type <- type_
    overall_risks_comparison$soc <- soc_
    all_overall_risks_comparisons <- 
      bind_rows(
        all_overall_risks_comparisons,
        overall_risks_comparison
      )
    
    substrate_drugs <- 
      sapply(substrates,function(x){str_split(x,"_")[[1]][1]}) %>% unique()
    nonsubstrate_drugs <- 
      sapply(nonsubstrates,function(x){str_split(x,"_")[[1]][1]}) %>% unique()  
    
    expand.grid(ade = substrates,probe=probes) %>% data.table() %>% nrow() %>% cat("\n\t\t",.)
    expand.grid(ade = nonsubstrates,probe=probes) %>% data.table() %>% nrow() %>% cat("\n\t\t",.)
    
    if(length(substrates)==0 | length(nonsubstrates)==0)next
    
    gene_dt <- 
      foreach(boot_=1:length(zs),.combine = "rbind",.errorhandling = "remove") %dopar% {
        z_ = zs[boot_]
        tmp <- data.table()
        
        grid <-
          bind_rows(
            expand.grid(ade = substrates,probe=probes,substrate="substrate") %>% data.table(),
            expand.grid(ade = nonsubstrates,probe=probes,substrate="nonsubstrate") %>% data.table()
          )
        
        tmp <- 
          lapply(1:nrow(grid),function(i){
            ade_ <- grid[i,ade]
            probe_ <- grid[i,probe]
            stype_ <- grid[i,substrate]
            mi_ <- 
              merge(
                dts_sub[ade==ade_,.(ade,nichd,score = (z_*gam_score_se)+gam_score)] %>% unique(),
                gsub_obs[probe==probe_,.(m = mean(residual)),nichd],
                by="nichd"
              ) %>% .[,maigesPack::MI(m,score)]
            cor_ <- 
              merge(
                dts_sub[ade==ade_,.(ade,nichd,score = (z_*gam_score_se)+gam_score)] %>% unique(),
                gsub_obs[probe==probe_,.(m = mean(residual)),nichd],
                by="nichd"
              ) %>% .[,cor(m,score,method="spearman")]
            data.table(ade = ade_,probe = probe_,mi=mi_,cor=cor_,substrate=stype_)
          }) %>% bind_rows()
        
        tmp$gene = gene_
        tmp$type = type_
        tmp$soc = soc_
        tmp$z = z_
        tmp$boot=boot_
        
        tmp
        
      }    
    
    smi = gene_dt[substrate=="substrate",median(mi),.(ade,probe)][,`V1`]
    nsmi = gene_dt[substrate=="nonsubstrate",median(mi),.(ade,probe)][,`V1`]
    denom = (length(smi)*length(nsmi))
    wtest <- 
      wilcox.test(smi,nsmi,paired = F,alternative="g")
    
    smi = gene_dt[substrate=="substrate"][,mi]
    nsmi = gene_dt[substrate=="nonsubstrate"][,mi]
    ttest <- 
      t.test(smi,nsmi,paired=F,alternative="g")
    
    all_gene_dt <- 
      bind_rows(
        all_gene_dt,
        gene_dt[,
                .(
                  gene = gene_,
                  type = type_,
                  soc=soc_,
                  auroc = 
                    (wtest$statistic/
                       denom),
                  wt_pvalue = 
                    wtest$p.value,
                  ttest_statistic = 
                    ttest$statistic,
                  ttest_pvalue = 
                    ttest$p.value
                )
        ]
      )
    
    
    all_full_gene_dt <- 
      bind_rows(
        all_full_gene_dt,
        gene_dt
      )
    
  }
})

all_gene_dt %>% 
  fwrite(paste0(data_dir,basename,"cyp_expression_soc_risk_information_summaries.csv"))

all_full_gene_dt %>% 
  fwrite(paste0(data_dir,basename,"cyp_expression_soc_risk_information.csv"))

all_overall_risks_comparisons %>% 
  fwrite(paste0(data_dir,basename,"cyp_expression_soc_risk_substrate_nonsubstrate_risks.csv"))

all_gsubs %>% 
  fwrite(paste0(data_dir,basename,"cyp_expression_soc_risk_probe_sample_residuals.csv"))

# ANALYSIS -- CYP substrate vs. nonsubstrate cluster enrichment ---------------------------

gene_type_soc_ade <- 
  fread(paste0(data_dir,"database_generation_cyp_gene_type_soc_ade.csv"))

iter <- 
  gene_type_soc_ade[!is.na(gene_symbol) & gene_symbol!="",
                    .(
                      Ndrugs = length(unique(atc_concept_id)),
                      Nades = length(unique(ade))
                    ),
                    .(gene_symbol,type)
  ] 

ade_er_table <- 
  fread(paste0(data_dir,"database_generation_er_tables/ade.csv.gz"))
cluster_ades <- ade_er_table[,.(cluster_name)] %>% unique()
cluster_enrichments <- NULL
for(i in 1:nrow(iter)){
  
  #i=4
  gene_=iter[i,gene_symbol]
  type_=iter[i,type]
  
  substrates <- 
    gene_type_soc_ade[
      gene_symbol==gene_ &
        type==type_,
      unique(ade)
    ]
  nonsubstrates <- 
    setdiff(
      gene_type_soc_ade[
        gene_symbol!=gene_,
        unique(ade)
      ],
      substrates
    )
  
  
  if(length(substrates)<2)next
  if(length(nonsubstrates)<2)next
  
  tmp <- 
    ade_er_table[
      ade %in% substrates,.(cluster_name,ade,substrate="substrate")
    ] %>% .[,.N,cluster_name] %>% transpose()
  colnames(tmp) <- tmp[1,] %>% unlist() %>% unname()
  
  cluster_dt <- 
    cluster_ades %>%
    merge(
      ade_er_table[
        ade %in% substrates,.(cluster_name,ade,substrate="substrate")
      ]
    )  %>%
    bind_rows(
      ade_er_table[
        ade %in% nonsubstrates,.(cluster_name,ade,substrate="nonsubstrate")
      ]
    ) %>% 
    dcast(cluster_name ~ substrate,fun.aggregate = length)
  
  ctest <- cluster_dt[,-c("cluster_name")] %>% as.matrix() %>% chisq.test()
  
  cluster_enrichments <- 
    bind_rows(
      cluster_enrichments,
      data.table(gene=gene_,type=type_,chi = ctest$statistic,pvalue = ctest$p.value) %>% 
        cbind(tmp[-c(1)])
    )
  
  
}

cluster_enrichments %>% .[order(pvalue)]
cluster_enrichments[gene=="CYP17A1"]
tmp <- 
  cluster_enrichments %>% 
  melt(
    id.vars=c("gene","type","chi","pvalue"),
    variable.name="cluster_name")
tmp[is.na(value),"value"] <- 0
tmp %>% 
  fwrite(paste0(data_dir,basename,"cyp_expression_risk_cluster_enrichments.csv"))

# ANALYSIS -- CYP substrate vs. nonsubstrate cluster enrichment for systemic disorders ---------------------------

gene_type_soc_ade <- 
  fread(paste0(data_dir,"database_generation_cyp_gene_type_soc_ade.csv"))

iter <- 
  gene_type_soc_ade[!is.na(gene_symbol) & gene_symbol!="",
                    .(
                      Ndrugs = length(unique(atc_concept_id)),
                      Nades = length(unique(ade))
                    ),
                    .(gene_symbol,soc,type)
  ] 

ade_er_table <- 
  fread(paste0(data_dir,"database_generation_er_tables/ade.csv.gz"))
cluster_ades <- ade_er_table[,.(cluster_name)] %>% unique()
cluster_enrichments <- NULL
for(i in 1:nrow(iter)){
    
    #i=4
    gene_=iter[i,gene_symbol]
    type_=iter[i,type]
    soc_=iter[i,soc]

    substrates <- 
      gene_type_soc_ade[
        gene_symbol==gene_ &
          type==type_ &
          soc==soc_,
        unique(ade)
      ]
    nonsubstrates <- 
      setdiff(
        gene_type_soc_ade[
          gene_symbol!=gene_ &
            soc==soc_,
          unique(ade)
        ],
        substrates
      )

    
    if(length(substrates)<2)next
    if(length(nonsubstrates)<2)next
    
    tmp <- 
      ade_er_table[
        ade %in% substrates,.(cluster_name,ade,substrate="substrate")
      ] %>% .[,.N,cluster_name] %>% transpose()
    colnames(tmp) <- tmp[1,] %>% unlist() %>% unname()
    
    cluster_dt <- 
      cluster_ades %>%
      merge(
        ade_er_table[
          ade %in% substrates,.(cluster_name,ade,substrate="substrate")
        ]
      )  %>%
      bind_rows(
        ade_er_table[
          ade %in% nonsubstrates,.(cluster_name,ade,substrate="nonsubstrate")
        ]
      ) %>% 
      dcast(cluster_name ~ substrate,fun.aggregate = length)

    ctest <- cluster_dt[,-c("cluster_name")] %>% as.matrix() %>% chisq.test()
    
    cluster_enrichments <- 
      bind_rows(
        cluster_enrichments,
        data.table(gene=gene_,type=type_,soc=soc_,chi = ctest$statistic,pvalue = ctest$p.value) %>% 
          cbind(tmp[-c(1)])
      )
    
    
    }

cluster_enrichments %>% .[order(pvalue)]
cluster_enrichments[gene=="CYP17A1"]
tmp <- 
  cluster_enrichments %>% 
  melt(
    id.vars=c("gene","soc","type","chi","pvalue"),
    variable.name="cluster_name")
tmp[is.na(value),"value"] <- 0
tmp %>% 
  fwrite(paste0(data_dir,basename,"cyp_expression_soc_risk_cluster_enrichments.csv"))


# VISUALIZATION -- CYP substrate vs. nonsubstrate MI test viz -----------------------------

gene_type_soc_ade <- 
  fread(paste0(data_dir,basename,"cyp_gene_type_soc_ade.csv"))

cluster_enrichments <- 
  fread(paste0(data_dir,basename,"cyp_expression_risk_cluster_enrichments.csv"))

all_gene_dt <- 
  fread(paste0(data_dir,basename,"cyp_expression_risk_information_summaries.csv"))

all_full_gene_dt <- 
  fread(paste0(data_dir,basename,"cyp_expression_risk_information.csv"))

all_overall_risks_comparisons <- 
  fread(paste0(data_dir,basename,"cyp_expression_risk_substrate_nonsubstrate_risks.csv"))

all_gsubs <- 
  fread(paste0(data_dir,basename,"cyp_expression_risk_probe_sample_residuals.csv"))

all_gene_dt$fdr <- p.adjust(all_gene_dt$ttest_pvalue,method="fdr")

all_gene_dt %>% 
  merge(
    all_full_gene_dt[,.(avg_mi = mean(mi)),.(gene,type,substrate)],
    by=c("gene","type")
  ) %>% 
  .[substrate=="substrate"] %>% 
  merge(
    cluster_enrichments[,.(gene,type,cluster_chisq_enrichment = chi)] %>% unique(),
    by=c("gene","type")
  ) %>% 
  .[order(ttest_pvalue)] %>% 
  .[,
    .(
      gene,
      cluster_chisq_enrichment = round(cluster_chisq_enrichment,2),
      ttest_statistic = round(ttest_statistic,2),
      fdr = scales::scientific(fdr,digits=3)
    )
  ] %>% 
  fwrite(paste0(data_dir,basename,"cyp_expression_risk_information_table.csv"))

all_gene_dt[,.N,type]
all_gene_dt[ttest_pvalue<0.05]

tmp <- 
  all_overall_risks_comparisons %>% 
  merge(
    null_dts[,.(null_99 = quantile(gam_score,c(0.99))),nichd],
    by="nichd"
  ) %>% 
  .[
    gene %in% all_gene_dt[ttest_pvalue<0.05,gene] &
      gam_score_90mse>null_99 & DE>0 &
      substrate_type=="substrate drug-events",
    .(gene,ade_name)
  ] %>% unique()

tmp[gene %in% all_gene_dt[ttest_pvalue<0.05,gene],length(unique(ade_name)),gene]
tmp[gene %in% all_gene_dt[ttest_pvalue<0.05,gene],.(gene,ade_name)] %>% unique()

# VISUALIZATION -- CYP substrate vs. nonsubstrate MI soc test viz -----------------------------

gene_type_soc_ade <- 
  fread(paste0(data_dir,basename,"cyp_gene_type_soc_ade.csv"))

soc_category <-  fread(paste0(data_dir,basename,"soc_categories.csv"))

cluster_enrichments <- 
  fread(paste0(data_dir,basename,"cyp_expression_soc_risk_cluster_enrichments.csv"))

all_gene_dt <- 
  fread(paste0(data_dir,basename,"cyp_expression_soc_risk_information_summaries.csv"))

all_full_gene_dt <- 
  fread(paste0(data_dir,basename,"cyp_expression_soc_risk_information.csv"))

all_overall_risks_comparisons <- 
  fread(paste0(data_dir,basename,"cyp_expression_soc_risk_substrate_nonsubstrate_risks.csv"))

all_gsubs <- 
  fread(paste0(data_dir,basename,"cyp_expression_soc_risk_probe_sample_residuals.csv"))

all_gene_dt$fdr <- p.adjust(all_gene_dt$ttest_pvalue,method="fdr")

mi_auroc_gene_dt <- 
  all_gene_dt %>% 
  merge(
    all_full_gene_dt[,.(avg_mi = mean(mi)),.(gene,type,soc,substrate)],
    by=c("gene","soc","type")
  ) 

all_gene_dt %>% 
  merge(
    all_full_gene_dt[,.(avg_mi = mean(mi)),.(gene,type,soc,substrate)],
    by=c("gene","soc","type")
  ) %>% 
  .[substrate=="substrate" & fdr<0.05,
    .(gene,soc)] %>% 
  merge(
    all_overall_risks_comparisons[
      substrate_type=="substrate",
      .(gene,soc,ade_name)],
    by=c("gene","soc")
  ) %>% 
  .[,length(unique(gene))]
all_gene_dt %>% 
  merge(
    all_full_gene_dt[,.(avg_mi = mean(mi)),.(gene,type,soc,substrate)],
    by=c("gene","soc","type")
  ) %>% 
  .[substrate=="substrate" & fdr<0.05,
    .(gene,soc)] %>% 
  merge(
    all_overall_risks_comparisons[
      substrate_type=="substrate",
      .(gene,soc,ade_name)],
    by=c("gene","soc")
  ) %>% 
  .[,length(unique(ade_name))]
all_gene_dt %>% 
  merge(
    all_full_gene_dt[,.(avg_mi = mean(mi)),.(gene,type,soc,substrate)],
    by=c("gene","soc","type")
  ) %>% 
  .[substrate=="substrate" & fdr<0.05 &
      gene %in% c("CYP2C18","CYP27B1"),
    .(gene,soc)] %>% 
  merge(
    all_overall_risks_comparisons[
      substrate_type=="substrate",
      .(gene,soc,ade_name)],
    by=c("gene","soc")
  ) %>% 
  .[,length(unique(ade_name))]

g <- all_gene_dt %>% 
  merge(
    all_full_gene_dt[,.(avg_mi = mean(mi)),.(gene,type,soc,substrate)],
    by=c("gene","soc","type")
  ) %>% 
  .[substrate=="substrate"] %>% 
  ggplot(aes(avg_mi,-log10(fdr),fill=soc)) +
  geom_point(size=3,pch=21) +
  geom_hline(yintercept = -log10(0.05),color="red",linetype="dashed") +
  colorspace::scale_fill_discrete_qualitative(palette="Dark 3") +
  xlab("Average mutual information") +
  ylab("-log10(FDR)") +
  guides(fill=guide_legend(ncol=1,title="System Organ Class")) +
  facet_wrap(~gene) +
  theme(
    strip.background = element_blank(),
    legend.position="right"
  )
ggsave(paste0(img_dir,basename,"cyp_expression_soc_risk_mi_versus_significance_by_gene.png"),g,width=18,height=7)

all_gene_dt %>% 
  merge(
    all_full_gene_dt[,.(avg_mi = mean(mi)),.(gene,type,soc,substrate)],
    by=c("gene","soc","type")
  ) %>% 
  .[fdr<0.05,length(unique(gene))]


all_overall_risks_comparisons[substrate_type=="substrate drug-events","substrate_type"] <- "substrate"
all_overall_risks_comparisons[substrate_type=="nonsubstrate drug-events","substrate_type"] <- "nonsubstrate"

g <- all_overall_risks_comparisons[
  nichd %in% gdata_stages & substrate_type=="substrate"
] %>% 
  merge(
    all_gene_dt %>% 
      .[fdr<0.01],
    by=c("gene","soc","type")
  ) %>% 
  .[,.(lwr = mean(gam_score) - sd(gam_score),
       avg = mean(gam_score),
       upr = mean(gam_score) + sd(gam_score)),
    .(nichd,substrate_type,soc,gene)
  ] %>% 
  merge(
    soc_category
  ) %>% 
  ggplot(aes(factor(nichd,levels=stages),
           avg,
           color=category)) +
  geom_path(size=1,group=1) +
  xlab("") +
  ylab("dGAM risk") +
  colorspace::scale_color_discrete_qualitative() +
  guides(color=guide_legend(title="Broad SOC category",title.position = "top",ncol=2)) +
  facet_wrap(~gene,scales="free_y",ncol=4) +
  theme(
    strip.background = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"cyp_expression_soc_risk_all_substrate_to_nonsubstrate.png"),g,width=12,height=10)

g <- all_gsubs[
  fdr<0.1,
  .(sample,probe,gene,nichd,residual)
] %>% 
  merge(
    all_gene_dt %>% 
      .[order(ttest_pvalue)] %>% 
      .[fdr<0.01,.SD[1],gene],
    by="gene"
  ) %>% 
  .[,
    .(avg = mean(residual)),
    .(nichd,probe,gene)
  ] %>% 
  .[nichd!=""] %>% 
  ggplot(aes(factor(nichd,levels=stages[2:6]),avg,group=probe)) +
  geom_point() +
  geom_line() +
  xlab("") +
  ylab("Probe residual expression") +
  facet_wrap(~gene,scales="free_y") +
  theme(
    strip.background = element_blank(),
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"cyp_expression_soc_risk_residual_probe_expression.png"),g,width=12,height=9)

cyp2c_ades <- all_full_gene_dt %>% 
  merge(
    all_gene_dt[grepl("^CYP2C",gene) & !(gene=="CYP2C18") & fdr<0.05]
  ) %>% .[substrate=="substrate",unique(ade)]
cyp2c18_ades <- all_full_gene_dt %>% 
  merge(
    all_gene_dt[(gene=="CYP2C18") & fdr<0.05]
  ) %>% .[substrate=="substrate",unique(ade)]

all_overall_risks_comparisons[
  ade %in% setdiff(cyp2c18_ades,cyp2c_ades) &
    substrate_type=="substrate",.(ade_name,soc)
  ] %>% unique()

g <- 
  plot_drug_events(
    dts_sub[
      ade %in% setdiff(cyp2c18_ades,cyp2c_ades)],
    color_name="ade"
    ) + 
  theme(legend.position = "none")
ggsave(paste0(img_dir,basename,"cyp2c18_specific_ades.png"),g,width=7,height=5)


