#' Title: "A database of pediatric drug effects to evaluate ontogenic mechanisms from child growth and development" study
#' 
#' Author details: Nicholas Giangreco
#' 
#' This script generates the data and plots evaluating adjustment 
#' in drug-event detection bias

# Purpose -----------------------------------------------------------------

#' To visualize and compare significant risks before and after covariate 
#' adjustment for ADEs with pediatric AEs
#' 

# Setup -------------------------------------------------------------------

pacman::p_load(tidyverse,data.table,doParallel)

data_dir <- "../../data/"
img_dir <- "../../docs/imgs/"
basename <- "database_generation_bias_"
seed = 0
set.seed(seed)

registerDoParallel(cores=5)

theme_set(theme_bw(base_size=16) + theme(text = element_text(face="bold")))

# Load data and covariates ------------------------------------------------

source("database_generation_load_data.R")

# Load GAM data ------------------------------------------------

source("database_generation_load_GAM_data.R")


# load pediatric aes -------------------------------------------------

pediatric_aes <- 
    fread(paste0(data_dir,"paediatric_term_list_19-0_concept_joined.csv"))


# join pediatric aes to gam data ------------------------------------------


sub <- 
dts[,
    .(ade,nichd,gam_score_90mse,meddra_concept_id,database)
    ] %>% 
    merge(
        null_dts[,.(null_99 = quantile(gam_score,c(0.99))),nichd],
        by="nichd"
    )

sub$pediatric_ae <- 
    sub$meddra_concept_id %in% pediatric_aes$concept_id

# enrichment for ade being nominally significant and pediatric ae ---------

dt <- data.table()

a <- 
    sub[database=="base" & pediatric_ae==T & gam_score_90mse>0,unique(ade)]
b <- 
    setdiff(sub[database=="base" & pediatric_ae==T & gam_score_90mse<=0,unique(ade)],a)
c <- 
    sub[database=="base" & pediatric_ae==F & gam_score_90mse>0,unique(ade)]
d <- 
    setdiff(sub[database=="base" & pediatric_ae==F & gam_score_90mse<=0,unique(ade)],c)

mat <- 
    matrix(
        c(
            length(a),
            length(b),
            length(c),
            length(d)
        ),
        byrow = T,nrow=2)
test <- fisher.test(mat)

dt$database <- "base"
dt$lwr <- test$conf.int[1]
dt$upr <- test$conf.int[2]
dt$odds <- test$estimate
dt$pvalue <- test$p.value

base_dt <- dt

dt <- data.table()

a <- 
    sub[database=="covariate_adjusted" & pediatric_ae==T & gam_score_90mse>0,unique(ade)]
b <- 
    setdiff(sub[database=="covariate_adjusted" & pediatric_ae==T & gam_score_90mse<=0,unique(ade)],a)
c <- 
    sub[database=="covariate_adjusted" & pediatric_ae==F & gam_score_90mse>0,unique(ade)]
d <- 
    setdiff(sub[database=="covariate_adjusted" & pediatric_ae==F & gam_score_90mse<=0,unique(ade)],c)

mat <- 
    matrix(
        c(
            length(a),
            length(b),
            length(c),
            length(d)
        ),
        byrow = T,nrow=2)
test <- fisher.test(mat)

dt$database <- "covariate_adjusted"
dt$lwr <- test$conf.int[1]
dt$upr <- test$conf.int[2]
dt$odds <- test$estimate
dt$pvalue <- test$p.value

ca_dt <- dt

g <- 
    rbind(
        base_dt,
        ca_dt
    ) %>% 
    ggplot(aes(database,odds)) +
    geom_point() +
    geom_errorbar(aes(ymin=lwr,ymax=upr),width=0.1) +
    xlab("") +
    ylab("Odds") +
    theme(
        axis.text.x = element_text(angle=45,vjust=1,hjust=1)
    )
ggsave(paste0(img_dir,basename,"nominally_significant_ade_enrichment.png"),g,width=4,height=4)

# enrichment for ade being nominally significant and pediatric ae ---------

dt <- data.table()

a <- 
    sub[database=="base" & pediatric_ae==T & gam_score_90mse>null_99,unique(ade)]
b <- 
    setdiff(sub[database=="base" & pediatric_ae==T & gam_score_90mse<=null_99,unique(ade)],a)
c <- 
    sub[database=="base" & pediatric_ae==F & gam_score_90mse>null_99,unique(ade)]
d <- 
    setdiff(sub[database=="base" & pediatric_ae==F & gam_score_90mse<=null_99,unique(ade)],c)

mat <- 
    matrix(
        c(
            length(a),
            length(b),
            length(c),
            length(d)
        ),
        byrow = T,nrow=2)
test <- fisher.test(mat)

dt$database <- "base"
dt$lwr <- test$conf.int[1]
dt$upr <- test$conf.int[2]
dt$odds <- test$estimate
dt$pvalue <- test$p.value

base_dt <- dt

dt <- data.table()

a <- 
    sub[database=="covariate_adjusted" & pediatric_ae==T & gam_score_90mse>null_99,unique(ade)]
b <- 
    setdiff(sub[database=="covariate_adjusted" & pediatric_ae==T & gam_score_90mse<=null_99,unique(ade)],a)
c <- 
    sub[database=="covariate_adjusted" & pediatric_ae==F & gam_score_90mse>null_99,unique(ade)]
d <- 
    setdiff(sub[database=="covariate_adjusted" & pediatric_ae==F & gam_score_90mse<=null_99,unique(ade)],c)

mat <- 
    matrix(
        c(
            length(a),
            length(b),
            length(c),
            length(d)
            ),
        byrow = T,nrow=2)
test <- fisher.test(mat)

dt$database <- "covariate_adjusted"
dt$lwr <- test$conf.int[1]
dt$upr <- test$conf.int[2]
dt$odds <- test$estimate
dt$pvalue <- test$p.value

ca_dt <- dt

g <- 
rbind(
    base_dt,
    ca_dt
) %>% 
    ggplot(aes(database,odds)) +
    geom_point() +
    geom_errorbar(aes(ymin=lwr,ymax=upr),width=0.1) +
    xlab("") +
    ylab("Odds") +
    theme(
        axis.text.x = element_text(angle=45,vjust=1,hjust=1)
    )
ggsave(paste0(img_dir,basename,"significant_ade_enrichment.png"),g,width=4,height=4)

# bias reduction using null model vs nominal significance -----------------

a <- 
    sub[database=="covariate_adjusted" & pediatric_ae==T & gam_score_90mse>null_99,unique(ade)]
b <- 
    setdiff(sub[database=="covariate_adjusted" & pediatric_ae==T & gam_score_90mse>0,unique(ade)],a)
c <- 
    sub[database=="covariate_adjusted" & pediatric_ae==F & gam_score_90mse>null_99,unique(ade)]
d <- 
    setdiff(sub[database=="covariate_adjusted" & pediatric_ae==F & gam_score_90mse>0,unique(ade)],c)

mat <- 
    matrix(
        c(
            length(a),
            length(b),
            length(c),
            length(d)
        ),
        byrow = T,nrow=2)
fisher.test(mat)