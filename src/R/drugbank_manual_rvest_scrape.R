# PURPOSE -----------------------------------------------------------------

#' manual sacrape drugbank pages to get correst drug-target-action mapping
#'

# load libraries and set variables ----------------------------------------


#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("org.Hs.eg.db",force=T)
library(AnnotationDbi)
library(org.Hs.eg.db)

pacman::p_load(rvest,purrr,data.table,doParallel,tidyverse,Biobase)

data_dir <- "../../data/"
drugbank_dir <- paste0(data_dir,"compound_drugbank05/")

registerDoParallel(cores=5)

# load drugbank data ------------------------------------------------------

drug_carriers <- fread(paste0(drugbank_dir,"drug_carriers.csv"))
drug_transporters <- fread(paste0(drugbank_dir,"drug_transporters.csv"))
drug_enzymes <- fread(paste0(drugbank_dir,"drug_enzymes.csv"))
drug_targets <- fread(paste0(drugbank_dir,"drug_targets.csv"))
drug_targets_actions <- fread(paste0(drugbank_dir,"drug_targets_actions.csv"))

# extract example data ----------------------------------------------------

zipra <- read_html("https://go.drugbank.com/drugs/DB00246")

zipra %>%
  html_element("#targets") %>%
  html_element("div.bond-list") %>%
  html_element("#BE0000145") %>%
  html_element("div.badge.badge-pill.badge-action") %>%
  html_text()

zipra %>%
  html_element("#targets") %>%
  html_element("div.bond-list") %>%
  html_element("#BE0000145") %>% 
  html_nodes("a") %>% 
  html_attr("href") %>% 
  .[3] %>% 
  str_split("uniprot/") %>% 
  .[[1]] %>% 
  .[2]

targets <-
  zipra %>%
  html_element("#targets") %>%
  html_element("div.bond-list") %>%
  html_children() %>%
  map_chr(function(x){
    x %>%
      html_attr("id")
  }) %>% 
  na.omit()

dts <- NULL
for(target in targets){
  action <- zipra %>%
    html_element("#targets") %>%
    html_element(paste0("#",target)) %>%
    html_element("div.badge.badge-pill.badge-action") %>%
    html_text()

  tmp <- zipra %>%
    html_element("#targets") %>%
    html_element("div.bond-list") %>%
    html_element(paste0("#",target)) %>% 
    html_nodes("a") %>% 
    html_attr("href")
  uniprot <- str_split(tmp[grepl("www.uniprot.org",tmp)],"uniprot/")[[1]][2]
  
  
  dts <- 
    bind_rows(dts,
              data.table(
                target = target, 
                action = action
                ) %>% 
                cbind(
                  AnnotationDbi::select(org.Hs.eg.db,
                         keys=uniprot,columns = c("ENTREZID","SYMBOL","GENENAME"),
                         keytype = "UNIPROT")
                )
              )
}
dts

# protein to action to drug scrapping/mapping --------------------------------------------------------------------

## TARGETS ##
drugs <- drug_targets[,unique(parent_key)]

action_protein_drug <-
  foreach(drug=drugs,.combine="rbind",.errorhandling = "remove") %dopar% {

    page <- read_html(paste0("https://go.drugbank.com/drugs/",drug))

    targets <-
      page %>%
      html_element("#targets") %>%
      html_element("div.bond-list") %>%
      html_children() %>%
      map_chr(function(x){
        x %>%
          html_attr("id")
      }) %>% 
      na.omit()

    dts <- NULL
    for(target in targets){
      action <- page %>%
        html_element("#targets") %>%
        html_element(paste0("#",target)) %>%
        html_element("div.badge.badge-pill.badge-action") %>%
        html_text()

      tmp <- page %>%
        html_element("#targets") %>%
        html_element("div.bond-list") %>%
        html_element(paste0("#",target)) %>% 
        html_nodes("a") %>% 
        html_attr("href")
      uniprot <- str_split(tmp[grepl("www.uniprot.org",tmp)],"uniprot/")[[1]][2]
      
      dts <-
        bind_rows(dts,
                  data.table(
                    target_id = target, 
                    action = stringr::str_to_lower(action)
                  ) %>% 
                    cbind(
                      AnnotationDbi::select(org.Hs.eg.db,
                             keys=uniprot,columns = c("ENTREZID","SYMBOL","GENENAME"),
                             keytype = "UNIPROT")
                    )
                  )
    }

    dts$parent_key <- drug

    dts

}

action_protein_drug %>%
  fwrite(paste0(drugbank_dir,"drug_targets_actions_rvest.csv"))
action_protein_drug %>%
  fwrite(paste0(data_dir,"drug_targets_actions_rvest.csv"))
## ##

## ENZYMES ##
drugs <- drug_enzymes[,unique(parent_key)]

action_protein_drug <-
  foreach(drug=drugs,.combine="rbind",.errorhandling = "remove") %dopar% {
    
    page <- read_html(paste0("https://go.drugbank.com/drugs/",drug))
    
    ps <-
      page %>%
      html_element("#enzymes") %>%
      html_element("div.bond-list") %>%
      html_children() %>%
      map_chr(function(x){
        x %>%
          html_attr("id")
      }) %>% 
      na.omit()
    
    dts <- NULL
    for(p in ps){
      action <- page %>%
        html_element("#enzymes") %>%
        html_element(paste0("#",p)) %>%
        html_element("div.badge.badge-pill.badge-action") %>%
        html_text()
      
      tmp <- page %>%
        html_element("#enzymes") %>%
        html_element("div.bond-list") %>%
        html_element(paste0("#",p)) %>% 
        html_nodes("a") %>% 
        html_attr("href")
      uniprot <- str_split(tmp[grepl("www.uniprot.org",tmp)],"uniprot/")[[1]][2]
      
      dts <-
        bind_rows(dts,
                  data.table(
                    enzyme_id = p, 
                    action = stringr::str_to_lower(action)
                  ) %>% 
                    cbind(
                      AnnotationDbi::select(org.Hs.eg.db,
                             keys=uniprot,columns = c("ENTREZID","SYMBOL","GENENAME"),
                             keytype = "UNIPROT")
                    )
        )
    }
    
    dts$parent_key <- drug
    
    dts
    
  }

action_protein_drug %>%
  fwrite(paste0(drugbank_dir,"drug_enzymes_actions_rvest.csv"))
action_protein_drug %>%
  fwrite(paste0(data_dir,"drug_enzymes_actions_rvest.csv"))

## ##

## TRANSPORTERS ##
drugs <- drug_transporters[,unique(parent_key)]

action_protein_drug <-
  foreach(drug=drugs,.combine="rbind",.errorhandling = "remove") %dopar% {
    
    page <- read_html(paste0("https://go.drugbank.com/drugs/",drug))
    
    ps <-
      page %>%
      html_element("#transporters") %>%
      html_element("div.bond-list") %>%
      html_children() %>%
      map_chr(function(x){
        x %>%
          html_attr("id")
      }) %>% 
      na.omit()
    
    dts <- NULL
    for(p in ps){
      action <- page %>%
        html_element("#transporters") %>%
        html_element(paste0("#",p)) %>%
        html_element("div.badge.badge-pill.badge-action") %>%
        html_text()
      
      tmp <- page %>%
        html_element("#transporters") %>%
        html_element("div.bond-list") %>%
        html_element(paste0("#",p)) %>% 
        html_nodes("a") %>% 
        html_attr("href")
      uniprot <- str_split(tmp[grepl("www.uniprot.org",tmp)],"uniprot/")[[1]][2]
      
      dts <-
        bind_rows(dts,
                  data.table(
                    transporter_id = p, 
                    action = stringr::str_to_lower(action)
                  ) %>% 
                    cbind(
                      AnnotationDbi::select(org.Hs.eg.db,
                             keys=uniprot,columns = c("ENTREZID","SYMBOL","GENENAME"),
                             keytype = "UNIPROT")
                    )
        )
      }
    
    dts$parent_key <- drug
    
    dts
    
  }

action_protein_drug %>%
  fwrite(paste0(drugbank_dir,"drug_transporters_actions_rvest.csv"))
action_protein_drug %>%
  fwrite(paste0(data_dir,"drug_transporters_actions_rvest.csv"))

## ##

## CARRIERS ##
drugs <- drug_carriers[,unique(parent_key)]

action_protein_drug <-
  foreach(drug=drugs,.combine="rbind",.errorhandling = "remove") %dopar% {
    
    page <- read_html(paste0("https://go.drugbank.com/drugs/",drug))
    
    ps <-
      page %>%
      html_element("#carriers") %>%
      html_element("div.bond-list") %>%
      html_children() %>%
      map_chr(function(x){
        x %>%
          html_attr("id")
      }) %>% 
      na.omit()
    
    dts <- NULL
    for(p in ps){
      action <- page %>%
        html_element("#carriers") %>%
        html_element(paste0("#",p)) %>%
        html_element("div.badge.badge-pill.badge-action") %>%
        html_text()
      
      tmp <- page %>%
        html_element("#carriers") %>%
        html_element("div.bond-list") %>%
        html_element(paste0("#",p)) %>% 
        html_nodes("a") %>% 
        html_attr("href")
      uniprot <- str_split(tmp[grepl("www.uniprot.org",tmp)],"uniprot/")[[1]][2]
      
      dts <-
        bind_rows(dts,
                  data.table(
                    carrier_id = p, 
                    action = stringr::str_to_lower(action)
                  ) %>% 
                    cbind(
                      AnnotationDbi::select(org.Hs.eg.db,
                             keys=uniprot,columns = c("ENTREZID","SYMBOL","GENENAME"),
                             keytype = "UNIPROT")
                    )
        )
    }
    
    dts$parent_key <- drug
    
    dts
    
  }

action_protein_drug %>%
  fwrite(paste0(drugbank_dir,"drug_carriers_actions_rvest.csv"))
action_protein_drug %>%
  fwrite(paste0(data_dir,"drug_carriers_actions_rvest.csv"))

## ##

# load drugbank actions rvest data and save cyp substrates ----------------------------------------

drugbank_dir <- paste0(data_dir,"compound_drugbank05/")
atc_db <- 
  fread(paste0(drugbank_dir,"drug_atc_codes_rxnorm_joined.csv")) %>% 
  .[,.(atc_concept_id,parent_key = `drugbank-id`)] %>% 
  unique()
atc_db_enzymes <- 
  merge(
    atc_db,
    merge(
      fread(paste0(drugbank_dir,"drug_enzymes_actions_rvest.csv")),
      fread(paste0(drugbank_dir,"drug_enzymes.csv")) %>% 
        .[,.(enzyme_id = id,enzyme_name = name,parent_key)] %>% 
        unique(),by=c("parent_key","enzyme_id")
    ),
    by="parent_key"
  )

atc_db_cyp_substrates <- 
  atc_db_enzymes[grepl("CYP",SYMBOL) &
                   action=="substrate"]

atc_db_cyp_substrates %>% 
  fwrite(paste0(data_dir,"drugbank_atc_cyp_substrates.csv"))

# load drugbank actions rvest data and save cyp substrates ----------------------------------------

drugbank_dir <- paste0(data_dir,"compound_drugbank05/")
atc_db <- 
  fread(paste0(drugbank_dir,"drug_atc_codes_rxnorm_joined.csv")) %>% 
  .[,.(atc_concept_id,parent_key = `drugbank-id`)] %>% 
  unique()
atc_db_transporters <- 
  merge(
    atc_db,
    merge(
      fread(paste0(drugbank_dir,"drug_transporters_actions_rvest.csv")),
      fread(paste0(drugbank_dir,"drug_transporters.csv")) %>% 
        .[,.(transporter_id = id,transporter_name = name,parent_key)] %>% 
        unique(),by=c("parent_key","transporter_id")
    ),
    by="parent_key"
  )

atc_db_transporter_substrates <- 
  atc_db_transporters[action=="substrate"]

atc_db_transporter_substrates %>% 
  fwrite(paste0(data_dir,"drugbank_atc_transporters_substrates.csv"))
