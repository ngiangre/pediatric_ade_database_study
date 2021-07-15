#' Title: "A database of pediatric drug effects to evaluate ontogenic mechanisms from child growth and development" study
#' 
#' Author details: Nicholas Giangreco
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

g <- gdata %>% 
  .[,.(gse,sample,nichd,PROBEID,value)] %>% 
  unique() %>% 
  na.omit() %>% 
  ggplot(aes(factor(nichd,levels=stages),log2(value))) +
  geom_boxplot(outlier.shape = NA) +
  xlab("") +
  ylab("Log2 gene expression") +
  facet_wrap(~gse,nrow=1) +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"normalized_expression_all_genes_across_stages_by_gse.png"),g,width=20,height=6)

g <- gdata %>% 
  .[,.(gse,assay,sample,nichd,PROBEID,value)] %>% 
  unique() %>% 
  na.omit() %>% 
  ggplot(aes(factor(nichd,levels=stages),log2(value))) +
  geom_boxplot(outlier.shape = NA) +
  xlab("") +
  ylab("Log2 gene expression") +
  facet_wrap(~assay) +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"normalized_expression_all_genes_across_stages_by_assay.png"),g,width=12,height=6)

# single gene expression across stages ------------------------------------

supp2 <- fread(paste0(data_dir,"stevens_et_al_supp2_probes.csv"))

genes <- c("KAT2B","DTX1","VEGFA","SKP1","CAMK2D","LEF1")
for(gene in genes){
  g <- gdata %>% 
    .[SYMBOL==gene] %>% 
    .[gse %in% c("GSE9006","GSE26440","GSE11504","TABM666"),
      .(gse,sample,nichd,stevens_group,PROBEID,value)] %>% 
    unique() %>% 
    na.omit() %>% 
    ggplot(aes(factor(nichd,levels=stages),log2(value))) +
    geom_boxplot() +
    geom_jitter(alpha=0.5) +
    xlab("") +
    ylab("Log2 expression values") +
    ggtitle(gene) +
    theme(
      axis.text.x = element_text(angle=45,vjust=1,hjust=1)
    )
  ggsave(paste0(img_dir,basename,gene,"_expression_across_stages.png"),g,width=7,height=6)
}

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
  guides(fill=guide_legend(title="GEO Series")) +
ggsave(paste0(img_dir,basename,"genomewide_pca_gse.png"),g,width=6,height=4)

g <- data.table(
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
    gdata[,.(sample,gse)] %>% unique(),
    by="sample"
  ) %>% 
  ggplot(aes(PC7,PC8,fill=gse)) +
  geom_jitter(pch=21,color="black")
ggsave(paste0(img_dir,basename,"genomewide_pca78_gse.png"),g,width=6,height=5)

summ_ <- summary(prcomp_)
data.table(
  pc = colnames(summ_$importance),
  var_prop = summ_$importance[2,]
  ) %>% 
  head(10) %>% 
  ggplot(aes(forcats::fct_inorder(pc),var_prop)) +
  geom_bar(stat="identity") +
  ylab("Proportion of variance") +
  xlab("")

g <- data.table(
  PC1 = prcomp_$x[,1],
  PC2 = prcomp_$x[,2],
  sample = colnames(tmp)[2:ncol(tmp)]
) %>% 
  merge(
    gdata[,.(sample,gse,sex)] %>% unique(),
    by="sample"
  ) %>% 
  ggplot(aes(PC1,PC2,fill=sex)) +
  geom_jitter(pch=21,color="black")
ggsave(paste0(img_dir,basename,"genomewide_pca_sex.png"),g,width=6,height=5)

g <- data.table(
  PC1 = prcomp_$x[,1],
  PC2 = prcomp_$x[,2],
  sample = colnames(tmp)[2:ncol(tmp)]
) %>% 
  merge(
    gdata[,.(sample,gse,nichd)] %>% unique(),
    by="sample"
  ) %>% 
  ggplot(aes(PC1,PC2,fill=factor(nichd,levels=stages))) +
  geom_jitter(pch=21,color="black") +
  guides(fill=guide_legend(title="NICHD",ncol=2,title.position = "top")) +
  theme(
    legend.position = "bottom"
  )
ggsave(paste0(img_dir,basename,"genomewide_pca_nichd.png"),g,width=6,height=5)

g <- data.table(
  PC1 = prcomp_$x[,1],
  PC2 = prcomp_$x[,2],
  sample = colnames(tmp)[2:ncol(tmp)]
) %>% 
  merge(
    gdata[,.(sample,gse,years)] %>% unique(),
    by="sample"
  ) %>% 
  ggplot(aes(PC1,PC2,fill=years)) +
  geom_jitter(pch=21,color="black") +
  colorspace::scale_fill_continuous_sequential() +
  theme(
    legend.position = "bottom"
  )
ggsave(paste0(img_dir,basename,"genomewide_pca_years.png"),g,width=6,height=5)

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

# Data validation -- Genome-wide expression ttest across stages ------------------

tests <- NULL
for(st1 in 1:length(gdata_stages)){
  for(st2 in 1:length(gdata_stages)){
    if((st2-st1)==1){
      a = gdata[nichd==gdata_stages[st1],.(sample,PROBEID,value)] %>% unique() %>% .[,value]
      b = gdata[nichd==gdata_stages[st2],.(sample,PROBEID,value)] %>% unique() %>% .[,value]
      test = t.test(sample(log2(a),10000,replace=T),sample(log2(b),10000,replace=T))
      tests <- 
        bind_rows(
          tests,
          data.table(nichd1=gdata_stages[st1],nichd2=gdata_stages[st2],pvalue=test$p.value)
        )
    }
  }
}

tests$fdr <- p.adjust(tests$pvalue,method="fdr")
tests[pvalue<0.05]
tests[fdr<0.05]


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


# Data validation -- enrichment --------------------------------------------------

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

# Data validation - enrichment plots --------------------------------------

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

g <- enrichment_data[alpha>0.01] %>% 
  ggplot(aes(alpha,odds)) +
  geom_line() +
  geom_point() +
  xlab("Alpha significance level") +
  ylab("Odds ratio") +
  geom_hline(yintercept=1,color="red",linetype="dashed")
ggsave(paste0(img_dir,basename,"stevens_enrichment_reduced_range.png"),g,width=6,height=5)


# Data validation -- example gene dynamics using residuals ------------------------

gene_="CYP2C19"
probes <- gdata[SYMBOL==gene_,unique(PROBEID)]
resid_dt <- lapply(
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
                   gse,sex,nichd)
        ],
        by="sample"
      )
    if(nrow(dat)==0){
      NA
    }else{
      test <- glm(value ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + gse,data=dat)
      data.table(
        probe = x,
        sample = dat$sample,
        nichd = dat$nichd,
        actual = dat$value,
        prediction = test$fitted.values,
        residual = test$residuals
      )
    }
  }
) %>% bind_rows()

resid_dt %>% 
  ggplot(aes(factor(nichd,levels=stages),residual)) +
  geom_boxplot() +
  geom_jitter(alpha=0.1) +
  xlab("") +
  ylab("Residual values") +
  ggtitle(gene_) +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )

# Data validation -- Generating correlation of stage association between values and residuals ------------------------------

probes <- gdata[,unique(PROBEID)]
resid_value_dt <- lapply(
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
      tmp <- glm(value ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + gse,data=dat)
      dat$residual <- tmp$residuals
      test <- glm(value ~ nichd + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + gse,data=dat)
      rtest <- glm(residual ~ nichd,data=dat)
      data.table(
        probe = x,
        nichd_coef = 
          summary(test)$coefficients[,1][rownames(summary(test)$coefficients)=="nichd"] %>% as.numeric(),
        pvalue = 
          summary(test)$coefficients[,4][rownames(summary(test)$coefficients)=="nichd"] %>% as.numeric(),
        nichd_rcoef = 
          summary(rtest)$coefficients[,1][rownames(summary(rtest)$coefficients)=="nichd"] %>% as.numeric(),
        rpvalue = 
          summary(rtest)$coefficients[,4][rownames(summary(rtest)$coefficients)=="nichd"] %>% as.numeric() 
      )
    }
  }
) %>% bind_rows()

resid_value_dt$rfdr <- 
  p.adjust(resid_value_dt$rpvalue,method="fdr")

resid_value_dt %>% 
  fwrite(paste0(data_dir,basename,"stage_association_by_residuals_valus.csv"))

# Data validation -- Viewing correlation of stage association between values and residuals ------------------------------

resid_value_dt <-  
  fread(paste0(data_dir,basename,"stage_association_by_residuals_valus.csv"))

sig_probes <- 
  resid_value_dt[rfdr<0.1,probe]

resid_value_dt %>% 
  ggplot(aes(pvalue,rpvalue)) +
  geom_point(alpha=0.5) +
  geom_abline(slope=1,intercept=0,color="red",linetype="dashed") +
  xlab("Stage-association of probe values") +
  ylab("Stage-association of probe residuals")

# Generating CYP expression using glm residuals ------------------------------

resid_value_dt <-  
  fread(paste0(data_dir,basename,"stage_association_by_residuals_valus.csv"))

cyp_genes <- gene_all[grepl("^CYP",gene_all)]
probes <- gdata[SYMBOL %in% cyp_genes,unique(PROBEID)]
resid_dt <- lapply(
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
                   gse,sex,nichd,nichd_int = as.integer(factor(nichd,levels=stages)))
        ],
        by="sample"
      )
    if(nrow(dat)==0){
      NA
    }else{
      test <- glm(value ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + gse,data=dat)
      dat$residual <- test$fitted.values
      
      data.table(
        probe = x,
        sample = dat$sample,
        nichd = dat$nichd,
        actual = dat$value,
        prediction = test$fitted.values,
        residual = test$residuals,
        f_statistic = 
          anova(glm(residual ~ nichd_int,data=dat),test = "F")$F[
            rownames(anova(glm(residual ~ nichd_int,data=dat),test = "F"))=="nichd_int"
            ],
        f_pvalue = 
          anova(glm(residual ~ nichd_int,data=dat),test = "F")$`Pr(>F)`[
            rownames(anova(glm(residual ~ nichd_int,data=dat),test = "F"))=="nichd_int"
            ]
      )
    }
  }
) %>% bind_rows()

resid_dt_gene <- 
  resid_dt %>% 
  merge(
    gdata[SYMBOL %in% cyp_genes,.(probe = PROBEID,sample,gene = SYMBOL)] %>% unique(),
    by=c("sample","probe"),
    allow.cartesian=T
  )

resid_dt_gene$fdr <- p.adjust(resid_dt_gene$f_pvalue,method="fdr")

resid_dt_gene$value = resid_dt_gene$residual

resid_dt_gene %>% 
  fwrite(paste0(data_dir,basename,"cyp_stage_association_to_residuals_samples_significance.csv"))

# Viewing CYP expression using glm residuals ------------------------------

resid_value_cyp_gene_dt <- 
  fread(paste0(data_dir,basename,"cyp_stage_association_to_residuals_samples_significance.csv"))

resid_value_cyp_gene_dt[pvalue<0.05,unique(gene)]

gene_="CYP2C8"
resid_value_cyp_gene_dt[gene==gene_] %>% 
  .[,.(sample,nichd,residual,gene)] %>% 
  unique() %>% 
  .[,.(
    lwr = quantile(residual,c(0.025)),
    avg = mean(residual),
    upr = quantile(residual,c(0.975))
  ),nichd] %>% 
  ggplot(aes(factor(nichd,levels=stages),avg)) +
  geom_point() +
  geom_errorbar(aes(ymin=lwr,ymax=upr),width=0.1) +
  xlab("") +
  ylab("Residual values") +
  ggtitle(paste0(gene_,"; Pvalue < ",resid_value_cyp_gene_dt[gene==gene_,round(unique(pvalue),18)])) +
  theme(
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )

# Relationship between CYP probe expression within genes ----------------------

resid_gene_dt <- 
  fread(paste0(data_dir,basename,"cyp_stage_association_to_residuals_samples_significance.csv"))

cyp_genes <- resid_gene_dt[f_pvalue<0.05,unique(gene)]

nprobes_dt <- 
  resid_dt_gene[gene %in% cyp_genes,
        .(Nprobes = length(unique(probe))),
        gene] %>% 
  .[order(gene)]

cor_mat_melts <- NULL
for(gene_ in cyp_genes){
  sub <- 
    resid_dt_gene[
      gene==gene_
    ] %>% 
    na.omit() %>% 
    .[,.(probe,sample,value,nichd)] %>% 
    unique() %>% 
    .[,
      .(m = mean(value)),
      .(nichd,probe)
    ] %>% 
    dcast(probe~ factor(nichd,levels=stages),value.var="m") %>% 
    transpose()
  
  probes <- sub[1] %>% unlist %>% unname
  res <- Hmisc::rcorr(as.matrix(sub[-c(1)]),type="pearson")
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

g <- cor_mat_melts %>% 
  ggplot(aes(value,gene)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  ylab("") +
  xlab("Pearson correlation between\nprobes in a CYP gene")
ggsave(paste0(img_dir,basename,"cyp_gene_correlation_boxplots.png"),g,width=6,height=6)

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

g <- cor_mat_melts %>% 
  .[(!PROBEID %in% discordant_probes &
       !PROBEID2 %in% discordant_probes)] %>% 
  .[order(gene)] %>% 
  ggplot(aes(value,forcats::fct_inorder(gene))) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  ylab("") +
  xlab("Pearson correlation between\nprobes in a CYP gene")
ggsave(paste0(img_dir,basename,"cyp_gene_concordant_correlation_boxplots.png"),g,width=6,height=6)

merge(
  nprobes_dt,
  nconcordantprobes_dt,
  all.x=T
) %>% 
  fwrite(paste0(data_dir,basename,"cyp_gene_Nprobes_Nconcordantprobes.csv"))

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
    gdata[SYMBOL %in% nprobes_dt[Nprobes==1,gene],unique(PROBEID)]
  )

genes_concordant_probes <- 
  cor_mat_melts %>% 
  .[(!PROBEID %in% discordant_probes &
       !PROBEID2 %in% discordant_probes),unique(gene)]

genes_concordant_and_single_probes <- 
  union(
    genes_concordant_probes,
    nprobes_dt[Nprobes==1,gene]
  )

# load GAM data -----------------------------------------------------------

source("database_generation_load_GAM_data.R")

sub <- 
  dts[
    ade %in% sig_null_ades & 
      database=="covariate_adjusted",
    .(atc_concept_id,meddra_concept_id,ade,
      nichd,gam_score,gam_score_se)
  ] %>% .[order(nichd)]

# load functions ----------------------------------------------------------

source("database_generation_functions.R")

# load sider side effects -------------------------------------------------

sider <- fread(paste0(data_dir,"database_generation_sider_standardized_named_data.csv"))

#sider <-  fread(paste0(data_dir,"database_generation_sider_standardized_named_data.csv"))

labeled_sig_ades <- 
  merge(
    sub[,c("atc_concept_id","meddra_concept_id","ade"),with=F] %>% unique(),
    sider[,c("atc_concept_id","meddra_concept_id"),with=F] %>% unique(),
    by=c("atc_concept_id","meddra_concept_id")
  ) %>% unique() %>% .[,ade]

# load drugbank data ------------------------------------------------------

drugbank_atc <- 
  fread(paste0(data_dir,"compound_drugbank05/drug_atc_codes_rxnorm_joined.csv"))

drugbank_cyp_substrates <- 
  fread(paste0(data_dir,"drugbank_atc_cyp_substrates.csv"))

genes_w_gt_5_sig_drug_substrates <- 
  drugbank_cyp_substrates[
    atc_concept_id %in% (sapply(labeled_sig_ades,function(x){str_split(x,"_")[[1]][1]}) %>% unique()),
    .(N = length(unique(atc_concept_id))),
    SYMBOL][N>5,SYMBOL]

drugbank_cyp_substrates[,
                        .(
                          Ndrugbank_cyp_substrate_drugs = 
                            length(unique(atc_concept_id))
                        ),
                        SYMBOL] %>% 
  merge(
    drugbank_cyp_substrates[
      atc_concept_id %in% (sapply(labeled_sig_ades,function(x){str_split(x,"_")[[1]][1]}) %>% unique()),
                            .(
                              Nsignificant_cyp_substrate_drugs = 
                                length(unique(atc_concept_id))
                            ),
                            SYMBOL],
    all.x=T
  ) %>% 
  fwrite(paste0(data_dir,basename,"drugbank_cyp_gene_Ndrugs_Nsigdrugs.csv"))

# tabulate ----------------------------------------------------------------

resid_value_cyp_gene_dt <- 
  fread(paste0(data_dir,basename,"cyp_stage_association_to_residuals_samples_significance.csv"))


resid_value_cyp_gene_dt[,length(unique(gene))]

length(genes_concordant_and_single_probes)
length(concordant_and_single_probes) 


length(sig_null_ades)
merge(
  sub[
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
  drugbank_cyp_substrates[,unique(atc_concept_id)],
  sapply(labeled_sig_ades,function(x){str_split(x,"_")[[1]][1]}) %>% unique()
) %>% length()
intersect(
  drugbank_cyp_substrates[,unique(atc_concept_id)],
  sapply(labeled_sig_ades,function(x){str_split(x,"_")[[1]][1]}) %>% unique()
) %>% length() / sapply(labeled_sig_ades,function(x){str_split(x,"_")[[1]][1]}) %>% unique() %>% length()
drugbank_cyp_substrates[,length(unique(SYMBOL))]

genes_for_risk_expression_comparison <- 
  intersect(genes_w_gt_5_sig_drug_substrates,genes_concordant_and_single_probes)
length(genes_for_risk_expression_comparison)

eval_dt <- 
  resid_value_cyp_gene_dt[
    fdr<0.05 &
    gene %in% genes_for_risk_expression_comparison &
      probe %in% concordant_and_single_probes
  ] %>% 
  .[,.(sample,nichd,value,gene,probe,f_statistic,f_pvalue,fdr)] %>% 
  unique() %>% 
  .[,.(m = mean(value)),.(probe,gene,nichd,f_pvalue,f_statistic,fdr)] %>% 
  merge(
    merge(
      sub[
        ade %in% labeled_sig_ades,
        .(ade,atc_concept_id,meddra_concept_id,ade_name,
          nichd,gam_score,gam_score_se)
        ],
      drugbank_cyp_substrates[,.(atc_concept_id,gene = SYMBOL)] %>% unique(),
      by=("atc_concept_id"),
      allow.cartesian = T
      ),
    by=c("nichd","gene"),
    allow.cartesian = T
  ) %>% 
  .[,.(ade,ade_name,nichd,gam_score,gam_score_se,meddra_concept_id,atc_concept_id,gene,probe,m,f_statistic,f_pvalue,fdr)] %>% 
  merge(
    sider[,.(meddra_concept_id,soc)] %>% unique(),
    by="meddra_concept_id",
    allow.cartesian = T
  )

eval_dt[,length(unique(gene))]
eval_dt[,length(unique(probe))]
eval_dt[,unique(gene)]

# CYP substrate vs. nonsubstrate MI test ---------------------------

resid_value_cyp_gene_dt <- 
  fread(paste0(data_dir,basename,"cyp_stage_association_to_residuals_samples_significance.csv"))

sub <- 
  dts[
    ade %in% sig_null_ades & 
      database=="covariate_adjusted",
    .(atc_concept_id,meddra_concept_id,ade,
      nichd,gam_score,gam_score_se,ade_name)
  ] %>% .[order(nichd)]

Nsamp = 100
zs = rnorm(Nsamp)
gene_dts <- NULL
for(gene_ in genes_for_risk_expression_comparison){
  
  #gene_="CYP2C18"
  probes <- eval_dt[gene==gene_,unique(probe)]
  gsub_obs <- resid_value_cyp_gene_dt[probe %in% probes]
  gsub_obs_mean <- 
    gsub_obs[,.(m = mean(residual)),nichd] %>% 
    na.omit()
  
  substrates <- 
    eval_dt[gene==gene_,unique(ade)]
  
  substrate_drugs <- 
    sapply(substrates,function(x){str_split(x,"_")[[1]][1]}) %>% unique()
  nonsubstrates <- 
    setdiff(labeled_sig_ades,substrates)

  nonsubstrate_drugs <- 
    sapply(nonsubstrates,function(x){str_split(x,"_")[[1]][1]}) %>% unique()  
  
  gene_dt <- 
    foreach(boot_=1:length(zs),.combine = "rbind",.errorhandling = "remove") %dopar% {
      z_ = zs[boot_]
      
      grid <- 
        expand.grid(ade = substrates,probe=probes) %>% data.table()
      dist_ <- c()
      for(i in 1:nrow(grid)){ 
        ade_ <- grid[i,ade]
        probe_ <- grid[i,probe]

        dist_ <- 
          c(dist_,
            merge(
              sub[ade==ade_,.(ade,nichd,score = (z_*gam_score_se)+gam_score)] %>% unique(),
              gsub_obs[probe==probe_,.(m = mean(residual)),nichd],
              by="nichd"
            ) %>% .[,maigesPack::MI(m,score)])
      }
      substrate_mis <- dist_
      
      grid <- 
        expand.grid(ade = nonsubstrates,probe=probes) %>% data.table()
      dist_ <- c()
      for(i in 1:nrow(grid)){ 
        ade_ <- grid[i,ade]
        probe_ <- grid[i,probe]

        dist_ <- 
          c(dist_,
            merge(
              sub[ade==ade_,.(ade,nichd,score = (z_*gam_score_se)+gam_score)] %>% unique(),
              gsub_obs[probe==probe_,.(m = mean(residual)),nichd],
              by="nichd"
            ) %>% .[,maigesPack::MI(m,score)])
      }
      nonsubstrate_mis <- dist_
      
      test <- wilcox.test(substrate_mis,nonsubstrate_mis,alternative = 'g')
      
      auroc = test$statistic/(length(substrate_mis)*length(nonsubstrate_mis))
      
      data.table(gene = gene_,Nprobes = length(probes),z = z_,
                 Nsubstrates = length(substrates),
                 Nsubstratedrugs = length(substrate_drugs),
                 Nnonsubstrates = length(nonsubstrates),
                 Nnonsubstratedrugs = length(nonsubstrate_drugs),
                 Nsubstratedist = length(substrate_mis),
                 Nnonsubstratedist = length(nonsubstrate_mis),
                 substrate_mi = mean(substrate_mis),
                 nonsubstrate_mi = mean(nonsubstrate_mis),
                 pvalue = test$p.value,
                 auroc = auroc,
                 U = test$statistic,
                 bootstrap = boot_
      )
      
    }
  
  g <- gene_dt %>% 
    ggplot() + 
    geom_histogram(
      aes(nonsubstrate_mi,fill="nonsubstrate"),
      color="black",alpha=0.5,bins = 30
      ) + 
    geom_histogram(
      aes(substrate_mi,fill="substrate"),
      color="black",alpha=0.5,bins = 30
    ) +
    guides(fill=guide_legend(title="substrate")) +
    scale_fill_manual(values=c("blue","red")) +
    xlab("Bootstrap mutual information")
  ggsave(paste0(img_dir,basename,gene_,"_bootstrap_mis.png"),g,width=7,height=4)
  
  g <- gene_dt %>% 
    ggplot(aes(pvalue)) +
    geom_histogram(color="black",fill="gray",bins = 30)
  ggsave(paste0(img_dir,basename,gene_,"_bootstrap_pvalues.png"),g,width=5,height=4)
  
  g <- gene_dt %>% 
    ggplot(aes(auroc)) +
    geom_histogram(color="black",fill="gray",bins = 30)
  ggsave(paste0(img_dir,basename,gene_,"_bootstrap_aurocs.png"),g,width=5,height=4)
  
  g <- gene_dt %>% 
    ggplot(aes(substrate_mi,nonsubstrate_mi,fill=auroc)) + 
    geom_jitter(color="black",pch=21,size=3) +
    geom_abline(slope=1,color="red",linetype="dashed") +
    colorspace::scale_fill_continuous_diverging(
      palette = "Blue-Red 3",name="AUROC",mid=(gene_dt[,max(auroc)] + gene_dt[,min(auroc)])/2
    ) +
    xlab("Substrate mutual information") +
    ylab("Nonsubstrate mutual information")
  ggsave(paste0(img_dir,basename,gene_,"_bootstrap_substrate_vs_nonsubstrate_information.png"),g,width=6,height=4)
  
  g <- gene_dt %>% 
    ggplot(aes(substrate_mi,nonsubstrate_mi,fill=z)) + 
    geom_jitter(color="black",pch=21,size=3) +
    geom_abline(slope=1,color="red",linetype="dashed") +
    colorspace::scale_fill_continuous_diverging(
      palette = "Blue-Red 3",name="Z score",mid=(gene_dt[,max(auroc)] + gene_dt[,min(auroc)])/2
    ) +
    xlab("Substrate mutual information") +
    ylab("nonsubstrate mutual information")
  ggsave(paste0(img_dir,basename,gene_,"_bootstrap_substrate_vs_nonsubstrate_information_zscore.png"),g,width=6,height=4)
  
  g <- gene_dt %>% 
    ggplot(aes(substrate_mi,nonsubstrate_mi,fill=z)) + 
    geom_jitter(color="black",pch=21,size=3) +
    geom_abline(slope=1,color="red",linetype="dashed",size=2) +
    facet_wrap(~factor(gene,levels=ord_genes),scales="free",ncol=1) +
    colorspace::scale_fill_continuous_diverging(
      palette = "Tofino",name="Z score",mid=(gene_dts[,max(z)] + gene_dts[,min(z)])/2,
      guide=guide_colorbar(barwidth = 15)
    ) +
    xlab("Substrate mutual information") +
    ylab("Nonsubstrate mutual information") +
    theme(
      strip.background = element_blank(),
      legend.position = "bottom"
    )
  ggsave(paste0(img_dir,basename,gene_,"_bootstrap_substrate_vs_nonsubstrate_information_zscore.png"),g,width=5,height=10)
  
  g <- sub[
    ade %in% substrates
  ] %>% 
    ggplot(aes(factor(nichd,levels=stages),gam_score)) + 
    geom_path(aes(group=ade)) + 
    xlab("") + 
    ylab("GAM risk score") +
    theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))
  ggsave(paste0(img_dir,basename,gene_,"_substrate_risks.png"),g,width=7,height=6)
  
  g <- sub[
    ade %in% substrates,
    .(
      lwr = mean(gam_score) - sd(gam_score),
      avg=mean(gam_score),
      upr = mean(gam_score) + sd(gam_score)
      ),
    nichd
  ] %>% 
    ggplot(aes(factor(nichd,levels=stages),avg)) + 
    geom_point() +
    geom_errorbar(aes(ymin=lwr,ymax=upr),width=0.1) +
    xlab("") + 
    ylab("GAM risk score") +
    theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))
  ggsave(paste0(img_dir,basename,gene_,"_substrate_risk_summary.png"),g,width=7,height=6)
  
  g <- sub[
    ade %in% nonsubstrates
  ] %>% 
    ggplot(aes(factor(nichd,levels=stages),gam_score)) + 
    geom_path(aes(group=ade),alpha=0.5) +
    xlab("") + 
    ylab("GAM risk score") +
    theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))
  ggsave(paste0(img_dir,basename,gene_,"_nonsubstrate_risks.png"),g,width=7,height=6)
  
  g <- sub[
    ade %in% nonsubstrates,
    .(
      lwr = mean(gam_score) - sd(gam_score),
      avg=mean(gam_score),
      upr = mean(gam_score) + sd(gam_score)
    ),
    nichd
  ] %>% 
    ggplot(aes(factor(nichd,levels=stages),avg)) + 
    geom_point() +
    geom_errorbar(aes(ymin=lwr,ymax=upr),width=0.1) +
    xlab("") + 
    ylab("GAM risk score") +
    theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))
  ggsave(paste0(img_dir,basename,gene_,"_nonsubstrate_risk_summary.png"),g,width=7,height=6)
  
  g <- 
    gsub_obs %>% 
    .[,
      .(
        lwr = mean(residual) - sd(residual),
        avg = mean(residual),
        upr = mean(residual) + sd(residual)
      ),
      .(probe,nichd)
    ] %>% 
    merge(
      gsub_obs_mean,
      by="nichd"
    ) %>% 
    ggplot(aes(factor(nichd,levels=stages),avg)) + 
    geom_point(aes(group=probe),position = position_dodge(width=0.5)) +
    geom_errorbar(aes(ymin=lwr,ymax=upr,group=probe),width=0.1,position = position_dodge(width=0.5)) +
    geom_line(aes(factor(nichd,levels=stages),m,group=probe),size=3) + 
    xlab("") + 
    ylab("Residual expression") +
    theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))
  ggsave(paste0(img_dir,basename,gene_,"_expression.png"),g,width=7,height=6)
  
  g <- 
    gsub_obs_mean %>% 
    ggplot(aes(factor(nichd,levels=stages),m)) + 
    geom_line(aes(group=1),size=3) + 
    xlab("") + 
    ylab("Mean residual expression") +
    theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))
  ggsave(paste0(img_dir,basename,gene_,"_mean_expression.png"),g,width=7,height=6)
  
  gene_dts <- 
    bind_rows(
      gene_dts,
      gene_dt
    )
  
}

gene_dts %>% 
  fwrite(paste0(data_dir,basename,"cyp_substrate_mean_MI_test_results.csv"))

gene_dts <- 
  fread(paste0(data_dir,basename,"cyp_substrate_mean_MI_test_results.csv"))

eval_dt[,.(gene,probe,
           f_statistic = round(f_statistic,3),
           f_pvalue = scales::scientific(f_pvalue,3),
           fdr = scales::scientific(fdr,3))] %>% 
  unique() %>% 
  merge(
    gene_dts[,
              .(
                Nsubstrates %>% unique(),
                Nnonsubstrates %>% unique(),
                auroc = 
                  (wilcox.test(substrate_mi,nonsubstrate_mi,alternative="g")$statistic/
                     (length(substrate_mi)*length(nonsubstrate_mi))) %>% round(3),
                wt_pvalue = 
                  wilcox.test(substrate_mi,nonsubstrate_mi,alternative="g")$p.value,
                ttest_statistic = 
                  t.test(substrate_mi,nonsubstrate_mi,alternative="g")$statistic %>% round(3),
                ttest_pvalue = 
                  t.test(substrate_mi,nonsubstrate_mi,alternative="g")$p.value),
              .(gene)],
    by="gene",
    allow.cartesian = T
    ) %>% 
  .[order(ttest_pvalue)] %>% 
  fwrite(paste0(data_dir,basename,"substrate_mi_comparison_ttest.csv"))

ord_genes <- 
  gene_dts[,
           .(
             ttest_pvalue = 
               t.test(substrate_mi,nonsubstrate_mi,alternative="g")$p.value
             ),
           gene] %>% 
  .[order(ttest_pvalue,decreasing = T)] %>% 
  .[ttest_pvalue<0.05,gene]

g <- 
  eval_dt %>% 
  .[gene %in% ord_genes] %>% 
  .[,
    .(
      lwr = quantile(gam_score,c(0.025)),
      score = mean(gam_score),
      upr = quantile(gam_score,c(0.975))
      ),
    .(gene,nichd)
    ] %>% 
  ggplot(aes(factor(nichd,levels=stages),score)) +
  geom_point() +
  geom_errorbar(aes(ymin=lwr,ymax=upr),width=0.1)+
  facet_wrap(~factor(gene,levels=ord_genes),scales = "free_y",ncol=1) +
  xlab("") +
  ylab("GAM risk score") +
  theme(
    strip.background = element_blank(),
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"cyp_drug_risks_across_stages.png"),g,width=5,height=10)

g <- 
  eval_dt %>% 
  .[,.(gene,probe,nichd,m)] %>% 
  unique() %>% 
  .[gene %in% ord_genes] %>% 
  ggplot(aes(factor(nichd,level=stages),m)) +
  geom_point() +
  geom_line(aes(group=probe)) +
  facet_wrap(~factor(gene,levels=ord_genes),scales = "free_y",ncol=1) +
  xlab("") +
  ylab("Probe expression") +
  theme(
    strip.background = element_blank(),
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"cyp_expression_from_random_across_stages.png"),g,width=5,height=10)

gene_dts[,
         .(
           avg_substrate_MI = mean(substrate_mi),
           avg_nonsubstrate_MI = mean(nonsubstrate_mi),
           greater_substrate_information = sum(substrate_mi>nonsubstrate_mi)/100
           ),
         gene
         ] %>% unique() %>% .[order(greater_substrate_information)]

g <- gene_dts %>% 
  .[gene %in% ord_genes] %>% 
  ggplot(aes(substrate_mi %>% round(2),nonsubstrate_mi %>% round(2),fill=auroc)) + 
  geom_jitter(color="black",pch=21,size=3) +
  geom_abline(slope=1,color="red",linetype="dashed",size=2) +
  facet_wrap(~factor(gene,levels=ord_genes),scales="free",ncol=1) +
  colorspace::scale_fill_continuous_diverging(
    palette = "Blue-Red 3",name="AUROC",mid=(gene_dts[,max(auroc)] + gene_dts[,min(auroc)])/2,
    guide=guide_colorbar(barwidth = 15)
    ) +
  xlab("Substrate mutual information") +
  ylab("Nonsubstrate mutual information") +
  theme(
    strip.background = element_blank(),
    legend.position = "bottom"
  )
ggsave(paste0(img_dir,basename,"bootstrap_substrate_vs_nonsubstrate_information.png"),g,width=5,height=10)

g <- gene_dts %>% 
  .[gene %in% ord_genes] %>% 
  ggplot(aes(substrate_mi %>% round(2),z)) + 
  geom_jitter(pch=21,size=3,color="black",fill="gray") +
  geom_smooth(method="lm",size=2) +
  facet_wrap(~factor(gene,levels=ord_genes),scales="free") +
  ylab("Z score") +
  xlab("Substrate mutual information")
ggsave(paste0(img_dir,basename,"bootstrap_substrate_information_vs_zscore.png"),g,width=11,height=6)

g <- gene_dts %>% 
  .[gene %in% ord_genes] %>% 
  ggplot(aes(substrate_mi %>% round(2),nonsubstrate_mi %>% round(2),fill=z)) + 
  geom_jitter(color="black",pch=21,size=3) +
  geom_abline(slope=1,color="red",linetype="dashed",size=2) +
  facet_wrap(~factor(gene,levels=ord_genes),scales="free",ncol=1) +
  colorspace::scale_fill_continuous_diverging(
    palette = "Tofino",name="Z score",mid=(gene_dts[,max(z)] + gene_dts[,min(z)])/2,
    guide=guide_colorbar(barwidth = 15)
  ) +
  xlab("Substrate mutual information") +
  ylab("Nonsubstrate mutual information") +
  theme(
    strip.background = element_blank(),
    legend.position = "bottom"
  )
ggsave(paste0(img_dir,basename,"bootstrap_substrate_vs_nonsubstrate_information_zscore.png"),g,width=5,height=10)

g <- gene_dts %>% 
  .[gene %in% ord_genes] %>% 
  ggplot(aes(substrate_mi,factor(gene,levels=ord_genes),fill=auroc)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(color="black",pch=21) +
  colorspace::scale_fill_continuous_sequential(
    name="AUROC",palette="Plasma",
    guide=guide_colorbar(barwidth = 10)
  ) +
  xlab("CYP-substrate\nmutual information") +
  ylab("") +
  theme(
    strip.background = element_blank(),
    legend.position = "bottom"
  )
ggsave(paste0(img_dir,basename,"bootstrap_substrate_information_by_auroc.png"),g,width=5,height=6)

g <- gene_dts %>% 
  .[gene %in% ord_genes] %>% 
  ggplot() + 
  geom_histogram(
    aes(substrate_mi,fill="substrate"),
    color="black",alpha=0.5
  ) + 
  geom_histogram(
    aes(nonsubstrate_mi,fill="nonsubstrate"),
    color="black",alpha=0.5
  ) +
  scale_fill_manual(values=c("blue","red"), 
                    name="",
                    breaks=c("nonsubstrate","substrate"),
                    labels=c("Nonsubstrate","Substrate")) +
  facet_wrap(~factor(gene,levels=ord_genes),scales="free",ncol=1) +
  guides(fill=guide_legend(title="substrate")) +
  xlab("Bootstrap mutual information") +
  ylab("Frequency") +
  theme(
    strip.background = element_blank(),
    legend.position = "bottom"
  )
ggsave(paste0(img_dir,basename,"bootstrap_mi.png"),g,width=5,height=10)

summ_zscore <- 
  gene_dts %>% 
  .[gene %in%
      ord_genes
    ] %>%
  .[
    substrate_mi>nonsubstrate_mi,
    .(
      lwr_zscore = quantile(z,c(0.025)),
      avg_zscore = mean(z),
      upr_zscore = quantile(z,c(0.975))),
    .(gene)
    ]

g <- eval_dt[,
        .(ade,nichd,gam_score,gam_score_se,gene)
        ] %>% 
  unique() %>% 
  merge(
    summ_zscore,
    by='gene'
  ) %>% 
  .[,
    .(
      lwr = mean(lwr_zscore*gam_score_se+gam_score),
      risk = mean(avg_zscore*gam_score_se+gam_score),
      upr = mean(upr_zscore*gam_score_se+gam_score)
    ),
      .(nichd,gene)
    ] %>% 
  ggplot(aes(factor(nichd,levels=stages),risk)) +
  geom_point() +
  geom_errorbar(aes(ymin=lwr,ymax=upr),width=0.1) +
  facet_wrap(~factor(gene,levels=ord_genes),ncol=1) +
  xlab("") +
  ylab("GAM risk score") +
  theme(
    strip.background = element_blank(),
    axis.text.x = element_text(angle=45,vjust=1,hjust=1)
  )
ggsave(paste0(img_dir,basename,"subtrate_gt_nonsubstrate_mi_zscore_risk.png"),g,width=5,height=10)

# CYP substrate vs. nonsubstrate MI test for systemic disorders ---------------------------

resid_value_cyp_gene_dt <- 
  fread(paste0(data_dir,basename,"cyp_stage_association_to_residuals_samples_significance.csv"))

sub <- 
  dts[
    ade %in% sig_null_ades & 
      database=="covariate_adjusted",
    .(atc_concept_id,meddra_concept_id,ade,
      nichd,gam_score,gam_score_se,ade_name)
  ] %>% .[order(nichd)]

Nsamp = 100
zs = rnorm(Nsamp)
gene_dts <- NULL
all_gsubs <- NULL
all_overall_risks_comparisons <- NULL
for(gene_ in genes_for_risk_expression_comparison){
  
  #gene_="CYP2D6"
  cat(gene_,"\n")
  probes <- eval_dt[gene==gene_,unique(probe)]
  gsub_obs <- resid_value_cyp_gene_dt[probe %in% probes]
  all_gsubs <- 
    bind_rows(
      all_gsubs,
      gsub_obs
    )
  gsub_obs_mean <- 
    gsub_obs[,.(m = mean(residual)),nichd] %>% 
    na.omit()
  socs_ = eval_dt[gene==gene_,.(N=length(unique(ade))),soc][,soc]
  
  if(length(socs_)==0){next}
  for(soc_ in socs_){
    cat("\t",soc_,"\n")
    
    events=sider[ade %in% labeled_sig_ades & soc==soc_,unique(meddra_concept_id)]
    
    substrates <- 
      sub[
        ade %in% labeled_sig_ades & 
          atc_concept_id %in% eval_dt[gene==gene_,unique(atc_concept_id)] &
          meddra_concept_id %in% events,
        unique(ade)
      ]
    nonsubstrates <- 
      sub[
        ade %in% labeled_sig_ades & 
          atc_concept_id %in% eval_dt[gene!=gene_,unique(atc_concept_id)] &
          meddra_concept_id %in% events,
        unique(ade)
      ] %>% setdiff(substrates)
    
    if(length(substrates)==0 | length(nonsubstrates)==0){next}
    
    substrate_risks <- 
      sub[
        ade %in% substrates
      ]
    substrate_risks$type <- "substrate drug-events"
    nonsubstrate_risks <- 
      sub[
        ade %in% nonsubstrates
      ]
    nonsubstrate_risks$type <- "nonsubstrate drug-events"
    overall_risks_comparison <- bind_rows(substrate_risks,nonsubstrate_risks)
    overall_risks_comparison$gene <- gene_
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
    
    gene_dt <- 
      foreach(boot_=1:length(zs),.combine = "rbind",.errorhandling = "remove") %dopar% {
        z_ = zs[boot_]
        
        grid <- 
          expand.grid(ade = substrates,probe=probes) %>% data.table()
        dist_ <- c()
        for(i in 1:nrow(grid)){ 
          ade_ <- grid[i,ade]
          probe_ <- grid[i,probe]
          
          dist_ <- 
            c(dist_,
              merge(
                sub[ade==ade_,.(ade,nichd,score = (z_*gam_score_se)+gam_score)] %>% unique(),
                gsub_obs[probe==probe_,.(m = mean(residual)),nichd],
                by="nichd"
              ) %>% .[,maigesPack::MI(m,score)])
        }
        substrate_mis <- dist_
        
        grid <- 
          expand.grid(ade = nonsubstrates,probe=probes) %>% data.table()
        dist_ <- c()
        for(i in 1:nrow(grid)){ 
          ade_ <- grid[i,ade]
          probe_ <- grid[i,probe]
          
          dist_ <- 
            c(dist_,
              merge(
                sub[ade==ade_,.(ade,nichd,score = (z_*gam_score_se)+gam_score)] %>% unique(),
                gsub_obs[probe==probe_,.(m = mean(residual)),nichd],
                by="nichd"
              ) %>% .[,maigesPack::MI(m,score)])
        }
        nonsubstrate_mis <- dist_
        
        test <- wilcox.test(substrate_mis,nonsubstrate_mis,alternative = 'g')
        
        auroc = test$statistic/(length(substrate_mis)*length(nonsubstrate_mis))
        
        data.table(gene = gene_,Nprobes = length(probes),
                   soc = soc_,z = z_,
                   Nsubstrates = length(substrates),
                   Nsubstratedrugs = length(substrate_drugs),
                   Nnonsubstrates = length(nonsubstrates),
                   Nnonsubstratedrugs = length(nonsubstrate_drugs),
                   Nsubstratedist = length(substrate_mis),
                   Nnonsubstratedist = length(nonsubstrate_mis),
                   substrate_mi = mean(substrate_mis),
                   nonsubstrate_mi = mean(nonsubstrate_mis),
                   z_pvalue = test$p.value,
                   z_auroc = auroc,
                   z_U = test$statistic,
                   bootstrap = boot_
        )
        
      }
    
    gene_dts <- 
      bind_rows(
        gene_dts,
        gene_dt %>% 
          merge(
            gene_dt[,
                    .(
                      auroc = 
                        (wilcox.test(substrate_mi,nonsubstrate_mi,alternative="g")$statistic/
                           (length(substrate_mi)*length(nonsubstrate_mi))),
                      wt_pvalue = 
                        wilcox.test(substrate_mi,nonsubstrate_mi,alternative="g")$p.value,
                      ttest_statistic = 
                        t.test(substrate_mi,nonsubstrate_mi,alternative="g")$statistic,
                      ttest_pvalue = 
                        t.test(substrate_mi,nonsubstrate_mi,alternative="g")$p.value
                      ),
                    .(gene,soc)
            ],
            by=c("gene","soc")
          )
      )
    
    gene_dt %>% 
      ggplot() + 
      geom_histogram(
        aes(nonsubstrate_mi,fill="nonsubstrate"),
        color="black",alpha=0.5,bins = 30
      ) + 
      geom_histogram(
        aes(substrate_mi,fill="substrate"),
        color="black",alpha=0.5,bins = 30
      ) +
      guides(fill=guide_legend(title="substrate")) +
      scale_fill_manual(values=c("blue","red")) +
      xlab("Bootstrap mutual information") +
      ggtitle(paste0(gene_," substrates and risks for\n",soc_)) 
    
    gene_dt %>% 
      ggplot(aes(substrate_mi,nonsubstrate_mi,fill=z_auroc)) + 
      geom_jitter(color="black",pch=21,size=3) +
      geom_abline(slope=1,color="red",linetype="dashed") +
      colorspace::scale_fill_continuous_diverging(
        palette = "Blue-Red 3",name="AUROC",mid=(gene_dt[,max(auroc)] + gene_dt[,min(auroc)])/2
      ) +
      xlab("Substrate mutual information") +
      ylab("Nonsubstrate mutual information") +
      ggtitle(paste0(gene_," substrates and risks for\n",soc_)) 
    
    gene_dt %>% 
      ggplot(aes(substrate_mi,nonsubstrate_mi,fill=z)) + 
      geom_jitter(color="black",pch=21,size=3) +
      geom_abline(slope=1,color="red",linetype="dashed") +
      colorspace::scale_fill_continuous_diverging(
        palette = "Blue-Red 3",name="Z score",mid=(gene_dt[,max(auroc)] + gene_dt[,min(auroc)])/2
      ) +
      xlab("Substrate mutual information") +
      ylab("nonsubstrate mutual information") +
      ggtitle(paste0(gene_," substrates and risks for\n",soc_)) 
    
    gene_dt %>% 
      ggplot(aes(substrate_mi,nonsubstrate_mi,fill=z)) + 
      geom_jitter(color="black",pch=21,size=3) +
      geom_abline(slope=1,color="red",linetype="dashed",size=2) +
      colorspace::scale_fill_continuous_diverging(
        palette = "Tofino",name="Z score",mid=(gene_dts[,max(z)] + gene_dts[,min(z)])/2,
        guide=guide_colorbar(barwidth = 15)
      ) +
      xlab("Substrate mutual information") +
      ylab("Nonsubstrate mutual information") +
      theme(
        strip.background = element_blank(),
        legend.position = "bottom"
      ) +
      ggtitle(paste0(gene_," substrates and risks for\n",soc_)) 
    
  }
  
}

gene_dts %>% 
  fwrite(paste0(data_dir,basename,"cyp_substrate_soc_risk_MI_test_results.csv"))

gene_dts <- 
  fread(paste0(data_dir,basename,"cyp_substrate_soc_risk_MI_test_results.csv"))

gene_dts_summary <- 
  gene_dts[,
           .(
             gene,
             soc,
             Nsubstrates,Nnonsubstrates,
             auroc,ttest_statistic,
             ttest_pvalue
           )] %>% unique() %>% .[order(ttest_pvalue)] 

gene_dts_summary[,length(unique(soc))]

gene_dts_summary$ttest_fdr <- 
  p.adjust(gene_dts_summary$ttest_pvalue,method="fdr")

gene_dts_summary[auroc>0.5 & ttest_fdr<0.05,length(unique(gene))]
gene_dts_summary[auroc>0.5 & ttest_fdr<0.05] %>% nrow() / gene_dts_summary %>% nrow()

gene_dts_summary %>% 
  fwrite(paste0(data_dir,basename,"substrate_mi_soc_comparison_ttest.csv"))


g <- gene_dts %>% 
  .[,.(substrate_mi = mean(substrate_mi)),.(gene,soc,auroc)] %>% 
  ggplot(aes(substrate_mi,auroc,fill=gene)) +
  geom_point(color="black",pch=21,size=4) +
  xlab("CYP-substrate\nmutual information") +
  ylab("AUROC") +
  guides(fill=guide_legend(ncol=4,title="")) +
  colorspace::scale_fill_discrete_qualitative(palette="Dynamic") +
  theme(
    strip.background = element_blank(),
    legend.position = "bottom"
  )
ggsave(paste0(img_dir,basename,"substrate_mi_soc_auroc.png"),g,width=7,height=5)

# Deep dive on drug,gene,event similar direction hypothesis ---------------------------

resid_value_cyp_gene_dt <- 
  fread(paste0(data_dir,basename,"cyp_stage_association_to_residuals_samples_significance.csv"))

sub <- 
  dts[
    ade %in% sig_null_ades & 
      database=="covariate_adjusted",
    .(atc_concept_id,meddra_concept_id,ade,
      nichd,gam_score,gam_score_se,ade_name)
  ] %>% .[order(nichd)]

Nsamp = 100
zs = rnorm(Nsamp)

gene_="CYP2D6"
probes <- eval_dt[gene==gene_,unique(probe)]
gsub_obs <- resid_value_cyp_gene_dt[probe %in% probes]
gsub_obs_mean <- 
  gsub_obs[,.(m = mean(residual)),nichd] %>% 
  na.omit()
  
gsub_obs %>% 
  .[,
    .(
      lwr = mean(residual) - sd(residual),
      avg = mean(residual),
      upr = mean(residual) + sd(residual)
    ),
    .(probe,nichd)
  ] %>% 
  ggplot(aes(factor(nichd,levels=stages),avg)) + 
  geom_point(aes(group=probe),position = position_dodge(width=0.5)) +
  geom_line(aes(group=probe)) +
  facet_wrap(~probe,ncol=1) +
  xlab("") + 
  ylab("Residual expression") +
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))

gsub_obs_mean %>% 
  ggplot(aes(factor(nichd,levels=stages),m)) + 
  geom_line(aes(group=1),size=3) + 
  xlab("") + 
  ylab("Mean residual expression") +
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))


ade_ = eval_dt[gene==gene_ & ade_name=="phenytoin and Drug Reaction With Eosinophilia And Systemic Symptoms",unique(ade)]
soc_=sider[ade==ade_,unique(soc)]
drug_="phenytoin"
drug_id_="21604399"
event_="Drug Reaction With Eosinophilia And Systemic Symptoms"
event_id_="43053854"
soc_=sider[meddra_concept_id==event_id_,unique(soc)][1]
events=sider[ade %in% labeled_sig_ades & soc==soc_,unique(meddra_concept_id)]

substrates <- 
  sub[
    ade %in% labeled_sig_ades & 
      atc_concept_id %in% eval_dt[gene==gene_,unique(atc_concept_id)] &
      meddra_concept_id %in% events,
    unique(ade)
    ]
nonsubstrates <- 
  sub[
    ade %in% labeled_sig_ades & 
      atc_concept_id %in% eval_dt[gene!=gene_,unique(atc_concept_id)],
    unique(ade)
  ] %>% setdiff(substrates)

grid <- 
  expand.grid(ade = substrates,probe=probes) %>% data.table()
dist_ <- c()
for(i in 1:nrow(grid)){ 
  ade_ <- grid[i,ade]
  probe_ <- grid[i,probe]
  
  dist_ <- 
    c(dist_,
      merge(
        sub[ade==ade_,.(ade,nichd,score = (z_*gam_score_se)+gam_score)] %>% unique(),
        gsub_obs[probe==probe_,.(m = mean(residual)),nichd],
        by="nichd"
      ) %>% .[,maigesPack::MI(m,score)])
}
substrate_mis <- dist_

grid <- 
  expand.grid(ade = nonsubstrates,probe=probes) %>% data.table()
dist_ <- c()
for(i in 1:nrow(grid)){ 
  ade_ <- grid[i,ade]
  probe_ <- grid[i,probe]
  
  dist_ <- 
    c(dist_,
      merge(
        sub[ade==ade_,.(ade,nichd,score = (z_*gam_score_se)+gam_score)] %>% unique(),
        gsub_obs[probe==probe_,.(m = mean(residual)),nichd],
        by="nichd"
      ) %>% .[,maigesPack::MI(m,score)])
}
nonsubstrate_mis <- dist_

data.table(
  mi = c(nonsubstrate_mis,substrate_mis),
  s = c(rep("nonsubstrate",length(nonsubstrate_mis)),rep("substrate",length(substrate_mis)))
  ) %>% 
  ggplot(aes(mi,fill=s)) +
  geom_histogram()
  

gene_dt <- 
    foreach(boot_=1:length(zs),.combine = "rbind",.errorhandling = "remove") %dopar% {
      z_ = zs[boot_]
      
      grid <- 
        expand.grid(ade = substrates,probe=probes) %>% data.table()
      dist_ <- c()
      for(i in 1:nrow(grid)){ 
        ade_ <- grid[i,ade]
        probe_ <- grid[i,probe]
        
        dist_ <- 
          c(dist_,
            merge(
              sub[ade==ade_,.(ade,nichd,score = (z_*gam_score_se)+gam_score)] %>% unique(),
              gsub_obs[probe==probe_,.(m = mean(residual)),nichd],
              by="nichd"
            ) %>% .[,maigesPack::MI(m,score)])
      }
      substrate_mis <- dist_
      
      grid <- 
        expand.grid(ade = nonsubstrates,probe=probes) %>% data.table()
      dist_ <- c()
      for(i in 1:nrow(grid)){ 
        ade_ <- grid[i,ade]
        probe_ <- grid[i,probe]
        
        dist_ <- 
          c(dist_,
            merge(
              sub[ade==ade_,.(ade,nichd,score = (z_*gam_score_se)+gam_score)] %>% unique(),
              gsub_obs[probe==probe_,.(m = mean(residual)),nichd],
              by="nichd"
            ) %>% .[,maigesPack::MI(m,score)])
      }
      nonsubstrate_mis <- dist_
      
      test <- wilcox.test(substrate_mis,nonsubstrate_mis,alternative = 'g')
      
      auroc = test$statistic/(length(substrate_mis)*length(nonsubstrate_mis))
      
      data.table(gene = gene_,Nprobes = length(probes),z = z_,
                 Nsubstrates = length(substrates),
                 Nsubstratedrugs = length(substrate_drugs),
                 Nnonsubstrates = length(nonsubstrates),
                 Nnonsubstratedrugs = length(nonsubstrate_drugs),
                 Nsubstratedist = length(substrate_mis),
                 Nnonsubstratedist = length(nonsubstrate_mis),
                 substrate_mi = mean(substrate_mis),
                 nonsubstrate_mi = mean(nonsubstrate_mis),
                 pvalue = test$p.value,
                 auroc = auroc,
                 U = test$statistic,
                 bootstrap = boot_
      )
      
    }
  
gene_dt %>% 
  ggplot() + 
  geom_histogram(
    aes(nonsubstrate_mi,fill="nonsubstrate"),
    color="black",alpha=0.5,bins = 30
  ) + 
  geom_histogram(
    aes(substrate_mi,fill="substrate"),
    color="black",alpha=0.5,bins = 30
  ) +
  guides(fill=guide_legend(title="substrate")) +
  scale_fill_manual(values=c("blue","red")) +
  xlab("Bootstrap mutual information")

gene_dt %>% 
    ggplot(aes(substrate_mi,nonsubstrate_mi,fill=auroc)) + 
    geom_jitter(color="black",pch=21,size=3) +
    geom_abline(slope=1,color="red",linetype="dashed") +
    colorspace::scale_fill_continuous_diverging(
      palette = "Blue-Red 3",name="AUROC",mid=(gene_dt[,max(auroc)] + gene_dt[,min(auroc)])/2
    ) +
    xlab("Substrate mutual information") +
    ylab("Nonsubstrate mutual information")

  
gene_dt %>% 
    ggplot(aes(substrate_mi,nonsubstrate_mi,fill=z)) + 
    geom_jitter(color="black",pch=21,size=3) +
    geom_abline(slope=1,color="red",linetype="dashed") +
    colorspace::scale_fill_continuous_diverging(
      palette = "Blue-Red 3",name="Z score",mid=(gene_dt[,max(auroc)] + gene_dt[,min(auroc)])/2
    ) +
    xlab("Substrate mutual information") +
    ylab("nonsubstrate mutual information")
  
gene_dt %>% 
    ggplot(aes(substrate_mi,nonsubstrate_mi,fill=z)) + 
    geom_jitter(color="black",pch=21,size=3) +
    geom_abline(slope=1,color="red",linetype="dashed",size=2) +
    facet_wrap(~gene,scales="free",ncol=1) +
    colorspace::scale_fill_continuous_diverging(
      palette = "Tofino",name="Z score",mid=(gene_dts[,max(z)] + gene_dts[,min(z)])/2,
      guide=guide_colorbar(barwidth = 15)
    ) +
    xlab("Substrate mutual information") +
    ylab("Nonsubstrate mutual information") +
    theme(
      strip.background = element_blank(),
      legend.position = "bottom"
    )

sub[
    ade %in% substrates,
    .(
      lwr = mean(gam_score) - sd(gam_score),
      avg=mean(gam_score),
      upr = mean(gam_score) + sd(gam_score)
    ),
    nichd
  ] %>% 
    ggplot(aes(factor(nichd,levels=stages),avg)) + 
    geom_point() +
    geom_errorbar(aes(ymin=lwr,ymax=upr),width=0.1) +
    xlab("") + 
    ylab("GAM risk score") +
  ggtitle("Substrate risk summary") +
    theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))

sub[
    ade %in% nonsubstrates,
    .(
      lwr = mean(gam_score) - sd(gam_score),
      avg=mean(gam_score),
      upr = mean(gam_score) + sd(gam_score)
    ),
    nichd
  ] %>% 
    ggplot(aes(factor(nichd,levels=stages),avg)) + 
    geom_point() +
    geom_errorbar(aes(ymin=lwr,ymax=upr),width=0.1) +
    xlab("") + 
    ylab("GAM risk score") +
  ggtitle("Nonsubstrate risk summary") +
    theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))


  