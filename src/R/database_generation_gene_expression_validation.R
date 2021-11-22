#' Title: "A database of pediatric drug effects to evaluate ontogenic mechanisms from child growth and development" study
#' 
#' Script author details: Nicholas Giangreco
#' 
#' This script generates the data and plots to validate the
#' gene expression across childhood dataset 

# Purpose -----------------------------------------------------------------

#' To evaluate significant temporal expression for gene products and evaluate correlation with ADE risks by drugs that are substrates for those gene products
#' 


# Setup -------------------------------------------------------------------

pacman::p_load(tidyverse,data.table,doParallel)

data_dir <- "../../data/"
img_dir <- "../../docs/imgs/"
basename <- "database_generation_gene_expression_validation_"
seed = 0
set.seed(seed)

registerDoParallel(cores=4)

stages <- 
  c("term_neonatal","infancy",
    "toddler","early_childhood",
    "middle_childhood","early_adolescence",
    "late_adolescence")

theme_set(theme_bw(base_size=16) + theme(text = element_text(face="bold")))


# load gene expression data -----------------------------------------------

gdata <- 
  fread(paste0(data_dir,"database_generation_stevens_et_al_raw_cel_processed_data.csv"))

gdata$value <- as.numeric(gdata$value)
gdata$nichd <- factor(gdata$nichd,levels=stages)

gdata_stages <- intersect(stages,gdata[,unique(nichd)])
gene_all <- gdata[,unique(SYMBOL)]
Nsamp=100
gene_random <- sample(gdata[,unique(SYMBOL)],Nsamp,replace=F)

# make stevens' age groups ------------------------------------------------

gdata <- 
  gdata %>% 
  merge(
    data.table(
      years_rounded = seq(0,17,1),
      stevens_group = c(2,rep(seq(2,16,2),each=2),17)
    ),
    by="years_rounded"
  )

# all expression across stages --------------------------------------------

g <- gdata %>% 
  .[,.(sample,nichd,PROBEID,value)] %>% 
  unique() %>% 
  na.omit() %>% 
  ggplot(aes(factor(nichd,levels=stages),log2(value))) +
  geom_boxplot(outlier.shape = NA) +
  xlab("") +
  ylab("Log2 gene expression") +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"normalized_expression_all_genes_across_stages.png"),g,width=8,height=6)

# single gene expression across stages ------------------------------------

supp2 <- fread(paste0(data_dir,"stevens_et_al_supp2_probes.csv"))

# Data validation -- PCA --------------------------------------------------

tmp <- gdata[,
             .(sample,PROBEID,value)
] %>% 
  .[is.finite(value)] %>% 
  unique() %>% 
  na.omit() %>% 
  dcast(PROBEID ~ sample,value.var = "value",fill=0)

mat <- 
  tmp[,-c("PROBEID")]  %>% 
  transpose()

dim(mat)

prcomp_ <- prcomp(mat,center = T,scale=T)

g <- data.table(
  PC1 = prcomp_$x[,1],
  PC2 = prcomp_$x[,2],
  sample = colnames(tmp)[2:ncol(tmp)]
) %>% 
  merge(
    gdata[,.(sample,gse)] %>% unique(),
    by="sample"
  ) %>% 
  ggplot(aes(PC1,PC2,fill=gse)) +
  geom_jitter(pch=21,color="black") +
  guides(fill=guide_legend(title="GEO Series"))
ggsave(paste0(img_dir,basename,"genomewide_pca_gse.png"),g,width=6,height=4)

pca_data <- 
  data.table(
    PC1 = prcomp_$x[,1],
    PC2 = prcomp_$x[,2],
    PC3 = prcomp_$x[,3],
    PC4 = prcomp_$x[,4],
    PC5 = prcomp_$x[,5],
    PC6 = prcomp_$x[,6],
    PC7 = prcomp_$x[,7],
    PC8 = prcomp_$x[,8],
    sample = colnames(tmp)[2:ncol(tmp)]
  ) %>% 
  merge(
    gdata[,.(sample,gse,sex,nichd)] %>% unique(),
    by="sample"
  )

pca_data %>% 
  fwrite(paste0(data_dir,basename,"pca_data.csv"))

# Data validation - uploading significant findings from expression --------

validation_data <- list()

validation_data[["Stevens"]] <- 
  fread(paste0(data_dir,"stevens_et_al_supp2_probes.csv")) %>% 
  .[,
    .(
      probe = `Probeset ID`, cluster = Cluster,fdr = `qvalue(Age Group)`
    )
  ] %>% 
  .[cluster!="Adult"] %>% 
  merge(
    gdata[,.(probe=PROBEID,gene=SYMBOL)] %>% unique(),
    all.x=T
  )


# Data validation -- enrichment and plots --------------------------------------------------

probes <- gdata[,unique(PROBEID)]
coef_pvalues <- lapply(
  probes,function(x){
    dat <-
      gdata[
        gse %in% c("GSE9006","GSE26440","GSE11504","TABM666") &
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
                   gse,sex,nichd = as.integer(factor(nichd,levels=stages)))
        ],
        by="sample"
      )
    if(nrow(dat)==0){
      NA
    }else{
      test <- glm(value ~ nichd + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + gse - 1,data=dat)
      data.table(
        probe = x,
        nichd_coef = summary(test)$coefficients[,1][rownames(summary(test)$coefficients)=="nichd"] %>% as.numeric(),
        pvalue = summary(test)$coefficients[,4][rownames(summary(test)$coefficients)=="nichd"] %>% as.numeric()
      )
    }
  }
) %>% bind_rows()

coef_pvalues$fdr <- 
  p.adjust(coef_pvalues$pvalue,method="fdr")

probe_gene_qvalue <- 
  coef_pvalues %>% 
  merge(
    gdata[,.(gene = SYMBOL,probe = PROBEID)] %>% unique(),
    by="probe"
  ) %>% 
  .[,.(probe,gene,nichd_coef,pvalue,fdr)]

probe_gene_qvalue %>% 
  fwrite(paste0(data_dir,basename,"our_age_associated_data.csv"))
validation_data$Stevens %>% 
  fwrite(paste0(data_dir,basename,"stevens_age_associated_data.csv"))

qs=probe_gene_qvalue[seq(100,probe_gene_qvalue[,length(probe)],100),fdr]
enrichment_data <- NULL
all_genes <- probe_gene_qvalue[,unique(probe)]
for(q_ in qs){
  
  our_sig_data <- probe_gene_qvalue[fdr<q_,unique(probe)]
  our_nonsig_data <- setdiff(all_genes,our_sig_data)
  their_sig_data <- validation_data[["Stevens"]][,unique(probe)]
  their_notsig_data <- setdiff(all_genes,their_sig_data)
  
  a <-
    intersect(
      our_sig_data,
      their_sig_data
    ) %>% length()
  b <-
    length(our_sig_data) - a
  c <-
    intersect(
      our_nonsig_data,
      their_sig_data
    ) %>% length()
  d <-
    length(our_nonsig_data) - c
  
  mat <- matrix(c(a,b,c,d),byrow = F,nrow=2,
                dimnames = list(c('our_qvalues_ltq_','not'),c('their_sig_qvalues','not')))
  
  test <- fisher.test(mat)
  
  enrichment_data <- 
    bind_rows(
      enrichment_data,
      data.table(
        a = a,
        b = b,
        c = c,
        d = d,
        a_plus_b = a+b,
        c_plus_d = c+d,
        our_data_set = union(our_sig_data,our_nonsig_data) %>% length(),
        validation_set = length(their_sig_data),
        odds_lwr = exp( log((a/b)/(c/d)) - ( 1.96*sqrt((1/a) + (1/b) + (1/c) + (1/d)) ) ),
        odds = (a/b)/(c/d),
        odds_upr = exp( log((a/b)/(c/d)) + ( 1.96*sqrt((1/a) + (1/b) + (1/c) + (1/d)) ) ),
        pvalue = test$p.value,
        alpha=q_,
        study="Stevens"
      )
    )
  
}

enrichment_data %>% 
  fwrite(paste0(data_dir,basename,"validation_enrichment.csv"))

g <- enrichment_data %>% 
  ggplot(aes(alpha,odds)) +
  geom_line() +
  geom_point() +
  xlab("Alpha significance level") +
  ylab("Enrichment odds ratio") +
  geom_hline(yintercept=1,color="red",linetype="dashed")
ggsave(paste0(img_dir,basename,"stevens_enrichment.png"),g,width=6,height=4)
