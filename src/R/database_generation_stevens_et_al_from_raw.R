#' Title: "A database of pediatric drug effects to evaluate ontogenic mechanisms from child growth and development" study
#' 
#' Script author details: Nicholas Giangreco
#' 
#' This script processes and generates the data for our
#' gene expression across childhood dataset

# Purpose -----------------------------------------------------------------

#' To download the raw (CEL) gene expression from GEO during childhood from Stevens et al. and then process
#' 

# Setup -------------------------------------------------------------------
# http://rstudio-pubs-static.s3.amazonaws.com/14438_729bb8430be242fba106c8ae3968458f.html
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("affy",force=T)
library(affy)
#BiocManager::install("org.Hs.eg.db",force=T)
library(org.Hs.eg.db)
#BiocManager::install("hgu133a.db",force=T)
library(hgu133a.db)
#BiocManager::install("hgu133b.db",force=T)
library(hgu133b.db)
#BiocManager::install("hgu133plus2.db",force=T)
library(hgu133plus2.db)
library(AnnotationDbi)
#BiocManager::install("illuminaHumanv3.db",force=T)
library(illuminaHumanv3.db)

pacman::p_load(tidyverse,data.table,doParallel,RCurl,GEOquery)

devtools::install_github("ngiangre/ROMOPOmics",force=T,dep=F)
library(ROMOPOmics)

data_dir <- "../../data/"
img_dir <- "../../docs/imgs/"
basename <- "database_generation_stevens_et_al_"
seed = 0
set.seed(seed)

registerDoParallel(cores=50)

stages <- 
    c("term_neonatal","infancy",
      "toddler","early_childhood",
      "middle_childhood","early_adolescence",
      "late_adolescence")

theme_set(theme_bw(base_size=16) + theme(text = element_text(face="bold")))

# list GSEs of interest and urls (put CEL files in GSE directories in an overall directory) --------------------------------------------------------

# Use the ftp site to download the CEL files
gse_ids <- 
    c("GSE9006", "GSE26440", "GSE11504", "TABM666", 
      "GSE6011", "GSE37721", "GSE20307", "GSE20436")

get_geo_cel_url  <- function(id_input = "GSE15896"){
    url   <- case_when(
        grepl("^GSE",id_input) ~ 
            paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/",
                   substr(id_input,1,nchar(id_input)-3),"nnn/GSE",
                   substr(id_input,4,stop = nchar(id_input)),
                   paste0("/suppl/",id_input,"_RAW.tar")),
        grepl("^GDS",id_input) ~
            paste0("https://ftp.ncbi.nlm.nih.gov/geo/datasets/",
                   substr(id_input,1,nchar(id_input)-3),"nnn/GDS",
                   substr(id_input,4,stop = nchar(id_input)),
                   paste0("/suppl/",id_input,"_RAW.tar"))
    )
    return(ifelse(RCurl::url.exists(url),url,NA))
}

urls <- 
    sapply(gse_ids,
           get_geo_cel_url
    )


# identify assays/metadata of GSE CEL data ------------------------------------------------

affy::bgcorrect.methods()
affy::normalize.AffyBatch.methods()
affy::pmcorrect.methods()
affy::express.summary.stat.methods()
# check out paper for combination evaluation
# https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gnj010

exprSet_files <- NULL
for(id_ in gse_ids){
  cat("Processing ",id_," at ",format(Sys.time())," ... \n")
  if(grepl("TAB",id_)){
    files = 
      sapply(paste0("../../data/stevens_et_al/E-TABM-666.raw.",c(1,2,3)),
             function(x){list.files(x,pattern=".CEL$",
                                    full.names = T)}) %>% unlist %>% unname()
  }else{
    files = list.files(
      paste0("../../data/stevens_et_al/",id_,"_RAW"),
      pattern="CEL.gz$",full.names=T)
  }
  exprSetdt <- 
    foreach(
      file=files,
      .combine = "rbind",
      .errorhandling = "remove"
    ) %dopar% {
      affy.data <- affy::ReadAffy(filenames = file,verbose = F)
      dt <- data.table()
      dt$sample <- str_split(colnames(Biobase::exprs(affy.data)),".CEL")[[1]][1]
      dt$file <- file
      dt$assay <- Biobase::annotation(affy.data)
      dt
    }
  exprSetdt$gse <- id_
  exprSet_files <- 
    bind_rows(
      exprSet_files,
      exprSetdt %>% data.table() %>% na.omit()
    )
}

exprSet_files[,.(Nsamples = length(unique(sample))),.(gse,assay)]

# process GSE CEL data ------------------------------------------------

assays <- exprSet_files[,unique(assay)] %>% na.omit()
exprSets <- NULL
for(assay_ in assays){
  files = exprSet_files[assay==assay_,unique(file)]
  affy.data <- affy::ReadAffy(filenames = files,verbose = F)
  eset.mas5 <- 
    expresso(affy.data,
             bgcorrect.method = "rma",
             normalize.method="quantiles",
             pmcorrect.method="pmonly",
             summary.method="avgdiff", 
             verbose = F)
  mat <- exprs(eset.mas5) 
  dt <- data.table(mat)
  dt$PROBEID <- rownames(mat)
  dt$assay <- assay_
  colnames(dt) <- sapply(str_split(colnames(dt),".CEL"),function(x){x[1]})
  dt_melt = dt %>% melt(id.vars=c("PROBEID","assay"),variable.name="sample")
  exprSets <- 
    bind_rows(
      exprSets,
      dt_melt
    )
  
}

merge(
  exprSets,
  exprSet_files[,.(file,sample,gse)],
  by=c("sample")
) %>% 
  fwrite(paste0(data_dir,basename,"raw_cel_all_preprocessed_data.csv"))

# post-process ------------------------------------------

exprSets <- 
    fread(paste0(data_dir,basename,"raw_cel_all_preprocessed_data.csv"))


exprSets[,length(unique(sample))]
exprSets[,length(unique(PROBEID))]
exprSets[,length(unique(sample)),assay]
exprSets[,length(unique(sample)),.(gse,assay)]

# download geo data via romopomics -------------------------------------------------------

stevens_gse_lst <- fetch_geo_series(gse_ids,data_dir = tempdir())

lst <- stevens_gse_lst[grep("GSE",names(stevens_gse_lst))]

# process phenotype data --------------------------------------------------

pdata <- NULL

proc <- 
  data.table(pData(lst$GSE9006$GSE$GSE$`GSE9006-GPL96_series_matrix.txt.gz`))
proc$age <- 
  sapply(str_split(proc$`Age:ch1`," [a-z]{3}"),function(x){as.numeric(x[1])})
proc[,'years' := fifelse(age>18,age/12,age)]
pdata <- bind_rows(
  pdata,
  proc[,.(sample = geo_accession,years,sex = proc$`Gender:ch1`,gse="GSE9006")]
)
proc <- 
  data.table(pData(lst$GSE9006$GSE$GSE$`GSE9006-GPL97_series_matrix.txt.gz`))
proc$age <- 
  sapply(str_split(proc$`Age:ch1`," [a-z]{3}"),function(x){as.numeric(x[1])})
proc[,'years' := fifelse(age>18,age/12,age)]
pdata <- bind_rows(
  pdata,
  proc[,.(sample = geo_accession,years,sex = proc$`Gender:ch1`,gse="GSE9006")]
)

proc <- 
  data.table(pData(lst$GSE26440$GSE$GSE$GSE26440_series_matrix.txt.gz))
proc$years <- as.numeric(proc$`age (years):ch1`)
pdata <- bind_rows(
  pdata,
  proc[,.(sample = geo_accession,years,sex = NA,gse="GSE26440")]
)

proc <- 
  data.table(pData(lst$GSE11504$GSE$GSE$GSE11504_series_matrix.txt.gz))
proc$years <- 
  sapply(str_split(proc$title," months"),
         function(x){
           as.numeric(str_split(x,"human age ")[[1]][2])/12
         })
pdata <- bind_rows(
  pdata,
  proc[,.(sample = geo_accession,years,sex = NA,gse="GSE11504")]
)

proc <- 
  data.table(pData(lst$GSE6011$GSE$GSE$GSE6011_series_matrix.txt.gz))
proc$years <- 
  sapply(str_split(proc$`Age:ch1`," months"),
         function(x){
           as.numeric(x[1])/12
         })
pdata <- bind_rows(
  pdata,
  proc[,.(sample = geo_accession,years,sex = ifelse(proc$`sex:ch1`=="Male","M","F"),gse = "GSE6011")]
)

proc <- 
  data.table(pData(lst$GSE37721$GSE$GSE$GSE37721_series_matrix.txt.gz))
proc$age <- 
  sapply(str_split(proc$title,"age "),
         function(x){
           str_split(x[2]," rep[0-9]")[[1]][1]
         })
proc$years <- 
  sapply(proc$age,
         function(x){
           gsub(" years", ".",x) %>% 
             trimws("left")
         }) %>% 
  unname() %>% 
  sapply(function(x){
    str_split(x," months")[[1]][1] %>% 
      gsub(" ","",.) %>% 
      gsub("and","",.) %>% 
      as.numeric()
  }) %>% 
  unname()
pdata <- bind_rows(
  pdata,
  proc[,.(sample = geo_accession,years,sex = ifelse(proc$`Sex:ch1`=="male","M","F"),gse="GSE37721")]
)


proc <- 
  data.table(pData(lst$GSE20307$GSE$GSE$GSE20307_series_matrix.txt.gz))
proc$years <- 
  sapply(proc$`age_at_onset:ch1`,
         function(x){
           if(grepl("GTE",x)){
             runif(n=1,min=6,max=18)
           }else{
             runif(n=1,min=0,max=6)
           }
         }) %>% unname
pdata <- bind_rows(
  pdata,
  proc[,.(sample = geo_accession,years,sex = ifelse(proc$`sex:ch1`=="Male","M","F"),gse="GSE20307")]
)


proc <- 
  data.table(pData(lst$GSE20436$GSE$GSE$GSE20436_series_matrix.txt.gz))
proc$years <- as.numeric(proc$`age:ch1`)
pdata <- bind_rows(
  pdata,
  proc[,.(sample = geo_accession,years,sex = ifelse(proc$`gender:ch1`=="Male","M","F"),gse="GSE20436")]
)

proc <- 
  fread(paste0(data_dir,"stevens_et_al/E-TABM-666.sdrf.txt")) %>% 
  .[,.(sample = `Scan Name`,years = `Characteristics [age]`,sex = ifelse(`Characteristics [sex]`=="male","M","F"),gse="TABM666")]
pdata <- bind_rows(
  pdata,
  proc
)

pdata[,length(unique(sample))]
pdata[,.N,gse]

expr_pdata <- 
  merge(
    exprSets,
    pdata,
    by=c('gse','sample')
  )

expr_pdata[,length(unique(sample))]
expr_pdata[,length(unique(sample)),gse]
expr_pdata[,length(unique(sample)),sex]
expr_pdata[,length(unique(sample)),.(assay,gse)]
expr_pdata[,length(unique(sample)),.(assay,gse,sex)]

rm(exprSets)
rm(stevens_gse_lst)
rm(lst)

# join feature data -------------------------------------------------------


hgu133a_dt <- 
  AnnotationDbi::select(hgu133a.db,
                        keys=expr_pdata[assay=="hgu133a",unique(PROBEID)],
                        columns = c("PROBEID","UNIPROT","ENTREZID","SYMBOL","GENENAME"),
                        keytype = "PROBEID") %>% data.table()
hgu133b_dt <- 
  AnnotationDbi::select(hgu133b.db,
                        keys=expr_pdata[assay=="hgu133b",unique(PROBEID)],
                        columns = c("PROBEID","UNIPROT","ENTREZID","SYMBOL","GENENAME"),
                        keytype = "PROBEID") %>% data.table()
plus2 <- 
  AnnotationDbi::select(hgu133plus2.db,
                        keys=expr_pdata[assay=="hgu133plus2",unique(PROBEID)],
                        columns = c("PROBEID","UNIPROT","ENTREZID","SYMBOL","GENENAME"),
                        keytype = "PROBEID") %>% data.table()
illum <- 
  AnnotationDbi::select(illuminaHumanv3.db,
                        keys=expr_pdata[assay=="illuminaHumanv3",unique(PROBEID)],
                        columns = c("PROBEID","UNIPROT","ENTREZID","SYMBOL","GENENAME"),
                        keytype = "PROBEID") %>% data.table()


expr_pfdata <- 
bind_rows(
  merge(
    expr_pdata[assay=="hgu133a"],
    hgu133a_dt,
    by="PROBEID",
    allow.cartesian=T,
    all.x=T
    ) %>% unique() %>% .[!is.na(UNIPROT)],
  merge(
    expr_pdata[assay=="hgu133b"],
    hgu133b_dt,
  by="PROBEID",
  allow.cartesian=T,
  all.x=T
  ) %>% unique() %>% .[!is.na(UNIPROT)],
  merge(
    expr_pdata[assay=="hgu133plus2" & !gse %in% c("GSE11504","GSE26440")],
    plus2,
    by="PROBEID",
    allow.cartesian=T,
    all.x=T
  ) %>% unique() %>% .[!is.na(UNIPROT)],
  merge(
    expr_pdata[assay=="illuminaHumanv3"],
    illum,
    by="PROBEID",
    allow.cartesian=T,
    all.x=T
  ) %>% unique() %>% .[!is.na(UNIPROT)],
  merge(
    expr_pdata[gse=="GSE26440"],
    plus2,
    by="PROBEID",
    allow.cartesian=T,
    all.x=T
  ) %>% unique() %>% .[!is.na(UNIPROT)],
  merge(
    expr_pdata[gse=="GSE11504"],
    plus2,
    by="PROBEID",
    allow.cartesian=T,
    all.x=T
  ) %>% unique() %>% .[!is.na(UNIPROT)],
  merge(
    expr_pdata[assay=="illuminaHumanv3"],
    illum,
    by="PROBEID",
    allow.cartesian=T,
    all.x=T
  ) %>% unique() %>% .[!is.na(UNIPROT)]
)

expr_pdata[,length(unique(sample)),.(assay,gse)]
expr_pfdata[,length(unique(sample)),.(assay,gse)]
expr_pfdata[,length(unique(sample))]

rm(expr_pdata)

# add stages --------------------------------------------------------------

dts <- expr_pfdata

dts[dts$years>0 & dts$years<=(1/12),"nichd"] <- "term_neonatal"
dts[dts$years>(1/12) & dts$years<=1,"nichd"] <- "infancy"
dts[dts$years>1 & dts$years<=2,"nichd"] <- "toddler"
dts[dts$years>2 & dts$years<=5,"nichd"] <- "early_childhood"
dts[dts$years>5 & dts$years<=11,"nichd"] <- "middle_childhood"
dts[dts$years>11 & dts$years<=18,"nichd"] <- "early_adolescence"
dts[dts$years>18 & dts$years<=21,"nichd"] <- "late_adolescence"
dts$years_rounded <- round(dts$years)

expr_pfdata <- dts

rm(dts)

# output ------------------------------------------------------------------


expr_pfdata %>% 
  .[,.(N = length(unique(sample))),.(gse,nichd)] %>% 
  .[,.(N,gse,nichd = factor(nichd,levels=stages))] %>% 
  na.omit() %>% 
  dcast(gse ~ nichd,value.var="N",fill = 0)

expr_pfdata %>% 
  .[,.(N = length(unique(sample)))]
expr_pfdata %>% 
  .[,.(N = length(unique(sample))),gse]


expr_pfdata %>% 
    fwrite(paste0(data_dir,basename,"raw_cel_processed_data.csv"))

