
# PURPOSE -----------------------------------------------------------------

#' To parse and clean sider data to use later
#' 

# Setup -------------------------------------------------------------------

pacman::p_load(data.table,tidyverse,doParallel)

data_dir <- "../../data/"
img_dir <- "../../docs/imgs/"
basename <- "database_generation_sider_"
seed = 0
set.seed(seed)

registerDoParallel(cores=5)

theme_set(theme_bw(base_size=16) + theme(text = element_text(face="bold")))


# download and store sider data -------------------------------------------------------

meddra_relationships <- 
    fread(paste0(data_dir,"standard_reactions_meddra_relationships.csv"))

db <- "effect_sider4"

con <- DBI::dbConnect(RMySQL::MySQL(),
                      user = readr::read_tsv("../../.my.cnf")$u,
                      password = readr::read_tsv("../../.my.cnf")$pw,
                      host = "127.0.0.1",
                      port=3307,
                      dbname=db)

meddra <- 
    tbl(con,"medDRA") %>% 
    collect() %>% 
    data.table()

atc_levels <- 
    tbl(con,"atc_levels") %>% 
    collect() %>% 
    data.table()


drug_atc <- 
    tbl(con,"drug_atc") %>% 
    collect() %>% 
    data.table()

meddra_all_label_indications <- 
    tbl(con,"medDRA_all_label_indications") %>% 
    collect() %>% 
    data.table()

DBI::dbDisconnect(con)

drugbank_atc <- 
    fread(paste0(data_dir,"compound_drugbank05/drug_atc_codes_rxnorm_joined.csv"))

meddra_all_label_indications[,
                             .(stitch_id_flat,method_of_detection,
                               meddra_umls_id = MedDRA_umls_id)
                             ]
meddra[MedDRA_concept_type=="PT",
       .(meddra_umls_id = umls_id,meddra_concept_id = MedDRA_id)
       ]

sider_merge <- 
merge(
    merge(
        meddra_all_label_indications[,
                                     .(stitch_id_flat,method_of_detection,
                                       meddra_umls_id = MedDRA_umls_id)
                                     ],
        merge(
            meddra[MedDRA_concept_type=="PT",
                   .(meddra_umls_id = umls_id,meddra_concept_code = MedDRA_id)
                   ],
            meddra_relationships[,.(meddra_concept_code = meddra_concept_code_1,
                                    meddra_concept_id = meddra_concept_id_1)] %>% 
                unique(),
            by="meddra_concept_code"
        ),
        by="meddra_umls_id"
    ),
    drug_atc[,.(stitch_id_flat,atc_concept_code = atc)],
    by="stitch_id_flat",
    allow.cartesian = T
)

merge(
    sider_merge,
    drugbank_atc[,.(atc_concept_code = atc_code,atc_concept_id)] %>% unique(),
    by="atc_concept_code"
) %>% 
    .[,.(meddra_umls_id,stitch_id_flat,meddra_concept_id,method_of_detection,atc_concept_id)] %>% 
    unique() %>% 
    fwrite(paste0(data_dir,basename,"standardized_data.csv"))
