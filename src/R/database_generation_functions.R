#' DATABASE GENERATION STUDY FUNCTIONS
#' 
#' Functions sourced by database_generation* scripts
#' 

#' Construct the drug-event dataset for GAMs 
#' 
#' Calls 'raw_data' in global environment and uses the function-defined drug, event, and drug-event column names to construct the dataset. Optional additional covariates and their ordered factors are also taken into consideration
#'
make_ade_data <- 
    function(id_col="safetyreportid",e_col="E",d_col="D",de_col="DE",x_cols=c("nichd"),category_levels=list()){
    # Extract IDs
    ids <- raw_data[,unique(get(id_col))]
    dt <- data.table()
    dt[,id_col] <- ids
    dt[,e_col] <- 0
    dt[,d_col] <- 0
    dt[,de_col] <- 0
    
    # Map stages to reports
    group_cols <- c(id_col,x_cols)
    reports_case_control_dat <-
        data.table::merge.data.table(dt,
                                     raw_data[,group_cols,with=F] %>% unique(),
                                     by=id_col
        ) %>% unique
    #Order reports by age
    reports_case_control_dat <- reports_case_control_dat[match(raw_data[,unique(safetyreportid)],get(id_col))]
    gam_data <- reports_case_control_dat 
    # Convert stage categories to integers
    group_cols <- c(d_col,e_col,de_col,x_cols)
    for(x_col in group_cols){
        x <- gam_data[,x_col,with=F] %>%
            unlist %>% unname
        if(x_col %in% names(category_levels)){
            gam_data[,paste0(x_col,"_name")] <- gam_data[,c(x_col),with=F]
            gam_data[,x_col] <- match(x,category_levels[[x_col]])
        }
    }
    
    # Output
    return(gam_data)
    
}

#' Set the drug and event indicator for the dataset
#'
set_ade_pair <- 
    function(ade_data,drug_ids,rx_ids,drug_name,rx_name,id_col="safetyreportid",e_col="E",d_col="D",de_col="DE"){
    #drug_ids=21604757;rx_ids=36919145;drug_name="methylphenidate";rx_name="Acute psychosis"
    # Extract event, drug, and not drug reports
    drug_reports <- raw_data[get(drug_col) %in% drug_ids,unique(get(id_col))]
    rx_reports <- raw_data[get(rx_col) %in% rx_ids,unique(get(id_col))]
    
    ade_data[,(e_col) := as.integer(get(id_col) %in% rx_reports)]
    ade_data[,(d_col) := as.integer(get(id_col) %in% drug_reports)]
    ade_data[,(de_col) := get(d_col) * get(e_col)]
    
    ade_data[,drug_name := drug_name]
    ade_data[,rx_name := rx_name]
    
    # Output
    return(ade_data)
    
}

#'Make a balanced train-test dataset split
#'
train_test_balanced_split <- 
    function(dat,test_split=0.2,response="E",seed=0){
    
    X <- dat %>% select(-!!as.name(response))
    y <- dat %>% select(!!as.name(response))
    
    case_inds <- which(y==1)
    ctrl_inds <- which(y==0)
    
    if((length(case_inds)+length(ctrl_inds))!=nrow(y)){errorCondition("unequal case-ctrl split")}
    
    set.seed(seed)
    test_case_inds <- 
        sample(case_inds,ceiling(length(case_inds)*test_split),replace=F)
    
    set.seed(seed)
    test_ctrl_inds <- 
        sample(ctrl_inds,floor(length(ctrl_inds)*test_split),replace=F)
    
    #(length(test_case_inds)+length(test_ctrl_inds))==ceiling((test_split*nrow(y)))
    
    train_case_inds <- 
        setdiff(case_inds,test_case_inds)
    train_ctrl_inds <- 
        setdiff(ctrl_inds,test_ctrl_inds)
    
    #(length(train_case_inds)+length(train_ctrl_inds))==floor(((1-test_split)*nrow(y)))
    
    X_test <- 
        X[c(test_case_inds,test_ctrl_inds)]
    X_train <- 
        X[c(train_case_inds,train_ctrl_inds)]
    
    y_test <- 
        y[c(test_case_inds,test_ctrl_inds)]
    y_train <- 
        y[c(train_case_inds,train_ctrl_inds)]
    
    dat_test <- 
        dat[c(test_case_inds,test_ctrl_inds)]
    dat_train <- 
        dat[c(train_case_inds,train_ctrl_inds)]
    
    if((nrow(dat_test)+nrow(dat_train))!=nrow(dat)){
        errorCondition("unequal train-test split")
    }else{
        return(list("train"= dat_train,
                    "test" = dat_test,
                    "X_train"= X_train,
                    "y_train" = y_train, 
                    "X_test" = X_test,
                    "y_test" = y_test))
    }
    
}

#'Compute AUROC from model
#'
get_auc <- 
    function(mod){
    predictions <- predict(mod, type="response")
    pred <- ROCR::prediction(as.numeric(predictions), as.numeric(mod$y))
    auc <- ROCR::performance( pred, "auc")@y.values[[1]]
    return(auc)
}
    
#'Compute AUROC from model on new data
#'
get_newdata_auc <- 
    function(mod,newdata,response="E",id_col="safetyreportid"){
        if(
            all(
                colnames(atc_covariates)[2:ncol(atc_covariates)] %in% 
                rownames(attr(mod$terms,"factors"))[
                    !rownames(attr(mod$terms,"factors")) %in% colnames(newdata)
                    ]
                )
        ){
            newdata <- merge(newdata,atc_covariates,by=id_col)
        }
        if(
            all(
                colnames(atc_covariates_l3)[2:ncol(atc_covariates_l3)] %in% 
                rownames(attr(mod$terms,"factors"))[
                    !rownames(attr(mod$terms,"factors")) %in% colnames(newdata)
                ]
            )
        ){
            newdata <- merge(newdata,atc_covariates_l3,by=id_col)
        }
    newdata <- newdata[,rownames(attr(mod$terms,"factors")),with=F]
    predictions <- predict(mod, newdata = newdata, newdata.guaranteed = T,type="response")
    pred <- ROCR::prediction(as.numeric(predictions), as.numeric(newdata[[response]]))
    auc <- ROCR::performance( pred, "auc")@y.values[[1]]
    return(auc)
}

#' Compute proportional reporting ratio of ADE across childhood
#' 
compute_prr <- 
    function(drug_id,rx_id,drug_name,rx_name,dat_new,id_col="safetyreportid",e_col="E",d_col="D",de_col="DE"){
    t0=Sys.time()
    # Set GAM data columns as factors
    dat_new[,d_col] <- factor(dat_new[,get(d_col)],levels=c(0,1))
    dat_new[,e_col] <- factor(dat_new[,get(e_col)],levels=c(0,1))
    # Make 2x2 ADE tables at stages
    tab <- dat_new[,
                   .(tab = list(table(D,E))
                   ),by=c(paste0(stage_col,"_name"))
    ]
    # Extract 2x2 parameters
    #https://www.ema.europa.eu/en/documents/regulatory-procedural-guideline/draft-guideline-use-statistical-signal-detection-methods-eudravigilance-data-analysis-system_en.pdf
    # stage order
    # a => drug, event
    # b => drug, no event
    # c => no drug, event
    # d => no drug, no event
    stage <- tab[,paste0(stage_col,"_name"),with=F] %>% unlist %>% unname
    a <- sapply(tab$tab,function(x){x["1","1"]})
    b <- sapply(tab$tab,function(x){x["1","0"]})
    # check:
    # dat_new[which(dat_new$D==1 & dat_new$E==0)][,table(get(stage_col))]
    # names(b) <- stage; b
    c <- sapply(tab$tab,function(x){x["0","1"]})
    d <- sapply(tab$tab,function(x){x["0","0"]})
    
    dt <- data.table()
    dt[,"a":=a]
    dt[,"b":=b]
    dt[,"c":=c]
    dt[,"d":=d]
    # (drug, events / all drug) / (no drug, events / all other drugs)
    dt[,"PRR":=(a/(a+b))/(c/(c+d))]
    dt[,"PRR_se":=sqrt( (1 / a) + (1 / c) - (1 / (a + b)) - (1/(c + d)) )]
    dt[,(stage_col):=stage]
    dt[,(drug_col_name):=drug_name]
    dt[,(rx_col_name):=rx_name]
    dt[,"ade_name"] = paste0(drug_name," and ",rx_name)
    dt[,(drug_col):=drug_id]
    dt[,(rx_col):=rx_id]
    dt[,"ade"] = paste0(drug_id,"_",rx_id)
    dt[,"PRR_90mse":=dt[,"PRR"]/exp(1.645*dt[,"PRR_se"])]
    dt[,"PRR_95mse":=dt[,"PRR"]/exp(1.96*dt[,"PRR_se"])]
    dt[,"PRR_99mse":=dt[,"PRR"]/exp(2.567*dt[,"PRR_se"])]
    dt[,"PRR_999mse":=dt[,"PRR"]/exp(3.291*dt[,"PRR_se"])]
    dt[,"PRR_90pse":=dt[,"PRR"]*exp(1.645*dt[,"PRR_se"])]
    dt[,"PRR_95pse":=dt[,"PRR"]*exp(1.96*dt[,"PRR_se"])]
    dt[,"PRR_99pse":=dt[,"PRR"]*exp(2.567*dt[,"PRR_se"])]
    dt[,"PRR_999pse":=dt[,"PRR"]*exp(3.291*dt[,"PRR_se"])]
    dt[,"time"] = difftime(Sys.time(),t0,"secs") %>% as.numeric()
    
    dt
    
}

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' Base: interaction effect of drug exposure during child development stages
#' 
compute_bam_base <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T){
    
    #bs="cs" gives lower AIC and df - a better fit model
    bam_mod <-
        mgcv::bam(E ~ s(nichd,
                  by=D,
                  bs=bs_,
                  k=length(stage_knots)),
            data=dat_new,
            knots = list(x=stage_knots),
            family=binomial(link="logit"),
            discrete=discrete,method = method
        )
    
    bam_mod
    
}

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' Stage: main effect of drug exposure across child development stages
#' 
compute_bam_with_only_stage_effects <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T){
    
    base_form <- "E ~ "
    
    nichd_form <- 
        deparse(~s(nichd,
                   bs=bs_,
                   k=length(stage_knots))) %>% substring(2)
    
    form <- paste0(base_form,nichd_form) %>% 
        as.formula()
    
    #bs="cs" gives lower AIC and df - a better fit model
    bam_mod <-
        mgcv::bam(form,
            data=dat_new,
            family=binomial(link="logit"),
            discrete = discrete,method = method
        )
    
    bam_mod
    
}

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' Sex: interaction effect of Female/Male sex across child development stages
#' 
compute_bam_with_only_sex_effects <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T){
    
    base_form <- "E ~ "
    sex_form <- 
        deparse(~s(nichd,
                   by=sex,
                   bs=bs_,
                   k=length(stage_knots)))%>% substring(2)
    
    form <- paste0(base_form,sex_form) %>% 
        as.formula()
    
    #bs="cs" gives lower AIC and df - a better fit model
    bam_mod <-
        mgcv::bam(form,
            data=dat_new,
            family=binomial(link="logit"),
            discrete = discrete,method = method
        )
    
    bam_mod
    
}

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' ATC: Effect of the number of drugs within anatomical/pharmacological (ATC 1st) classes e.g. XA (Alimentary tract and metabolism drugs - 'X' is just for making the name valid in a model formula)
#' 
compute_bam_with_only_atc_effects <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T,id_col="safetyreportid"){
    
    dat_new <- merge(dat_new,atc_covariates,by=id_col)
    
    base_form <- "E ~ " 
    
    cov_form <- paste0(colnames(atc_covariates)[2:ncol(atc_covariates)],collapse=" + ")
    
    form <- paste0(base_form," + ",cov_form) %>% 
        as.formula()
    
    #bs="cs" gives lower AIC and df - a better fit model
    bam_mod <-
        mgcv::bam(form,
            data=dat_new,
            family=binomial(link="logit"),
            discrete = discrete,method = method
        )
    
    bam_mod
    
}

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' ATCbin: Effect of drugs within anatomical/pharmacological (ATC 1st) classes e.g. XA (Alimentary tract and metabolism drugs - 'X' is just for making the name valid in a model formula)
#' 
compute_bam_with_only_atcbin_effects <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T,id_col="safetyreportid"){
        
        dat_new <- merge(dat_new,atc_covariates_bin,by=id_col)
        
        base_form <- "E ~ " 
        
        cov_form <- paste0(colnames(atc_covariates_bin)[2:ncol(atc_covariates_bin)],collapse=" + ")
        
        form <- paste0(base_form," + ",cov_form) %>% 
            as.formula()
        
        #bs="cs" gives lower AIC and df - a better fit model
        bam_mod <-
            mgcv::bam(form,
                      data=dat_new,
                      family=binomial(link="logit"),
                      discrete = discrete,method = method
            )
        
        bam_mod
        
    }

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' ATC3: Effect of the number of drugs within therapeutic (ATC 2nd) classes e.g. XA0 (Drugs used in diabetes - 'X' is just for making the name valid in a model formula)
#' 
compute_bam_with_only_atc3_effects <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T,id_col="safetyreportid"){
        
        dat_new <- merge(dat_new,atc_covariates_l3,by=id_col)
        
        base_form <- "E ~ " 
        
        cov_form <- paste0(colnames(atc_covariates_l3)[2:ncol(atc_covariates_l3)],collapse=" + ")
        
        form <- paste0(base_form," + ",cov_form) %>% 
            as.formula()
        
        #bs="cs" gives lower AIC and df - a better fit model
        bam_mod <-
            mgcv::bam(form,
                      data=dat_new,
                      family=binomial(link="logit"),
                      discrete = discrete,method = method
            )
        
        bam_mod
        
    }

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' ATC3bin: Effect of drugs within therapeutic (ATC 2nd) classes e.g. XA0 (Drugs used in diabetes - 'X' is just for making the name valid in a model formula)
#' 
compute_bam_with_only_atc3bin_effects <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T,id_col="safetyreportid"){
        
        dat_new <- merge(dat_new,atc_covariates_l3_bin,by=id_col)
        
        base_form <- "E ~ " 
        
        cov_form <- paste0(colnames(atc_covariates_l3_bin)[2:ncol(atc_covariates_l3_bin)],collapse=" + ")
        
        form <- paste0(base_form," + ",cov_form) %>% 
            as.formula()
        
        #bs="cs" gives lower AIC and df - a better fit model
        bam_mod <-
            mgcv::bam(form,
                      data=dat_new,
                      family=binomial(link="logit"),
                      discrete = discrete,method = method
            )
        
        bam_mod
        
    }

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' Date: Effect of the date of reporting across the time period of report collection
#' 
compute_bam_with_only_date_effects <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T){
    
    base_form <- "E ~ "
    
    rd_form <- 
        deparse(
            ~s(receive_date,
               bs=bs_)) %>% substring(2)
    
    form <- paste0(base_form,rd_form) %>% 
        as.formula()
    
    #bs="cs" gives lower AIC and df - a better fit model
    bam_mod <-
        mgcv::bam(form,
            data=dat_new,
            family=binomial(link="logit"),
            discrete = discrete,method = method
        )
    
    bam_mod
    
}

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' Reporter: Effect of the type of reporter e.g. Physician
#' 
compute_bam_with_only_reporter_effects <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T){
    
    base_form <- "E ~ "
    
    rq_form <- 
        deparse(
            ~te(reporter_qualification,
                bs=bs_,
                k=length(category_levels$reporter_qualification))) %>% substring(2)
    
    form <- paste0(base_form,rq_form) %>% 
        as.formula()
    
    #bs="cs" gives lower AIC and df - a better fit model
    bam_mod <-
        mgcv::bam(form,
            data=dat_new,
            family=binomial(link="logit"),
            discrete = discrete,method = method
        )
    
    bam_mod
    
}

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' Base+Age: Base interaction and the main effect of child development stage
#' 
compute_bam_with_age_effects <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T){
    
    base_form <- 
        deparse(E ~ s(nichd,
                      by=D,
                      bs=bs_,
                      k=length(stage_knots)))
    nichd_form <- 
        deparse(~s(nichd,
                   bs=bs_,
                   k=length(stage_knots))) %>% substring(2)
    
    form <- paste0(base_form," + ",nichd_form) %>% 
        as.formula()
    
    #bs="cs" gives lower AIC and df - a better fit model
    bam_mod <-
        mgcv::bam(form,
            data=dat_new,
            family=binomial(link="logit"),
            discrete = discrete,method = method
        )
    
    bam_mod
    
}

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' Base+Sex: Base interaction and the sex interaction across child development stages
#' 
compute_bam_with_sex_effects <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T){
    
    base_form <- 
        deparse(E ~ s(nichd,
                      by=D,
                      bs=bs_,
                      k=length(stage_knots)))
    sex_form <- 
        deparse(~s(nichd,
                   by=sex,
                   bs=bs_,
                   k=length(stage_knots)))%>% substring(2)
    
    form <- paste0(base_form," + ",sex_form) %>% 
        as.formula()
    
    #bs="cs" gives lower AIC and df - a better fit model
    bam_mod <-
        mgcv::bam(form,
            data=dat_new,
            family=binomial(link="logit"),
            discrete = discrete,method = method
        )
    
    bam_mod
    
}

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' Base+NdrugsL: Base interaction and the smooth effect of the number of drugs
#' 
compute_bam_with_ndrugs_Leffects <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T){
        
        base_form <- 
            deparse(E ~ s(nichd,
                          by=D,
                          bs=bs_,
                          k=length(stage_knots)))
        
        form <- paste0(base_form," + polypharmacy") %>% 
            as.formula()
        
        #bs="cs" gives lower AIC and df - a better fit model
        bam_mod <-
            mgcv::bam(form,
                      data=dat_new,
                      family=binomial(link="logit"),
                      discrete = discrete,method = method
            )
        
        bam_mod
        
    }

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' Base+NdrugsS: Base interaction and the smooth effect of the number of drugs
#' 
compute_bam_with_ndrugs_Seffects <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T){
        
        base_form <- 
            deparse(E ~ s(nichd,
                          by=D,
                          bs=bs_,
                          k=length(stage_knots)))
        
        ndrugs_form <- 
            deparse(~s(polypharmacy,
                       bs=bs_)) %>% substring(2)    
        
        form <- paste0(base_form," + ",ndrugs_form) %>% 
            as.formula()
        
        #bs="cs" gives lower AIC and df - a better fit model
        bam_mod <-
            mgcv::bam(form,
                      data=dat_new,
                      family=binomial(link="logit"),
                      discrete = discrete,method = method
            )
        
        bam_mod
        
    }


#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' Base+ATC: Base interaction and the effect of the number of drugs within anatomical/pharmacological classes
#' 
compute_bam_with_atc_effects <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T,id_col="safetyreportid"){
    
    dat_new <- merge(dat_new,atc_covariates,by=id_col)
    
    base_form <- 
        deparse(E ~ s(nichd,
                      by=D,
                      bs=bs_,
                      k=length(stage_knots)))
    
    cov_form <- paste0(colnames(atc_covariates)[2:ncol(atc_covariates)],collapse=" + ")
    
    form <- paste0(base_form," + ",cov_form) %>% 
        as.formula()
    
    #bs="cs" gives lower AIC and df - a better fit model
    bam_mod <-
        mgcv::bam(form,
            data=dat_new,
            family=binomial(link="logit"),
            discrete = discrete,method = method
        )
    
    bam_mod
    
}

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' Base+ATCbin: Base interaction and the effect of drugs within anatomical/pharmacological classes
#' 
compute_bam_with_atcbin_effects <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T,id_col="safetyreportid"){
        
        dat_new <- merge(dat_new,atc_covariates_bin,by=id_col)
        
        base_form <- 
            deparse(E ~ s(nichd,
                          by=D,
                          bs=bs_,
                          k=length(stage_knots)))
        
        cov_form <- paste0(colnames(atc_covariates_bin)[2:ncol(atc_covariates_bin)],collapse=" + ")
        
        form <- paste0(base_form," + ",cov_form) %>% 
            as.formula()
        
        #bs="cs" gives lower AIC and df - a better fit model
        bam_mod <-
            mgcv::bam(form,
                      data=dat_new,
                      family=binomial(link="logit"),
                      discrete = discrete,method = method
            )
        
        bam_mod
        
    }

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' Base+ATC3: Base interaction and the effect of the number of drugs within therapeutic classes
#' 
compute_bam_with_atc3_effects <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T,id_col="safetyreportid"){
        
        dat_new <- merge(dat_new,atc_covariates_l3,by=id_col)
        
        base_form <- 
            deparse(E ~ s(nichd,
                          by=D,
                          bs=bs_,
                          k=length(stage_knots)))
        
        cov_form <- paste0(colnames(atc_covariates_l3)[2:ncol(atc_covariates_l3)],collapse=" + ")
        
        form <- paste0(base_form," + ",cov_form) %>% 
            as.formula()
        
        #bs="cs" gives lower AIC and df - a better fit model
        bam_mod <-
            mgcv::bam(form,
                      data=dat_new,
                      family=binomial(link="logit"),
                      discrete = discrete,method = method
            )
        
        bam_mod
        
    }

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' Base+ATC3bin: Base interaction and the effect of drugs within therapeutic classes
#' 
compute_bam_with_atc3bin_effects <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T,id_col="safetyreportid"){
        
        dat_new <- merge(dat_new,atc_covariates_l3_bin,by=id_col)
        
        base_form <- 
            deparse(E ~ s(nichd,
                          by=D,
                          bs=bs_,
                          k=length(stage_knots)))
        
        cov_form <- paste0(colnames(atc_covariates_l3_bin)[2:ncol(atc_covariates_l3_bin)],collapse=" + ")
        
        form <- paste0(base_form," + ",cov_form) %>% 
            as.formula()
        
        #bs="cs" gives lower AIC and df - a better fit model
        bam_mod <-
            mgcv::bam(form,
                      data=dat_new,
                      family=binomial(link="logit"),
                      discrete = discrete,method = method
            )
        
        bam_mod
        
    }

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' Base+Date: Base interaction and the effect of the date of reporting
#' 
compute_bam_with_date_effects <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T){
    
    base_form <- 
        deparse(E ~ s(nichd,
                      by=D,
                      bs=bs_,
                      k=length(stage_knots)))
    rd_form <- 
        deparse(
            ~s(receive_date,
               bs=bs_)) %>% substring(2)
    
    form <- paste0(base_form," + ",rd_form) %>% 
        as.formula()
    
    #bs="cs" gives lower AIC and df - a better fit model
    bam_mod <-
        mgcv::bam(form,
            data=dat_new,
            family=binomial(link="logit"),
            discrete=discrete,method = method
        )
    
    bam_mod
    
}

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' Base+Reporter: Base interaction and the effect of the reporter type
#' 
compute_bam_with_reporter_effects <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T){
    
    base_form <- 
        deparse(E ~ s(nichd,
                      by=D,
                      bs=bs_,
                      k=length(stage_knots)))
    rq_form <- 
        deparse(
            ~te(reporter_qualification,
                bs=bs_,
                k=length(category_levels$reporter_qualification))) %>% substring(2)
    
    form <- paste0(base_form," + ",rq_form) %>% 
        as.formula()
    
    #bs="cs" gives lower AIC and df - a better fit model
    bam_mod <-
        mgcv::bam(form,
            data=dat_new,
            family=binomial(link="logit"),
            discrete = discrete,method = method
        )
    
    bam_mod
    
}

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' Base+Reporter+Date: Base interaction and the effect of the reporter type and date of report
#' 
compute_bam_with_reporter_date_effects <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T){
        
        base_form <- 
            deparse(E ~ s(nichd,
                          by=D,
                          bs=bs_,
                          k=length(stage_knots)))
        rq_form <- 
            deparse(
                ~te(reporter_qualification,
                    bs=bs_,
                    k=length(category_levels$reporter_qualification))) %>% substring(2)
        rd_form <- 
            deparse(
                ~s(receive_date,
                   bs=bs_)) %>% substring(2)
        
        form <- paste0(base_form," + ",rq_form," + ",rd_form) %>% 
            as.formula()
        
        #bs="cs" gives lower AIC and df - a better fit model
        bam_mod <-
            mgcv::bam(form,
                      data=dat_new,
                      family=binomial(link="logit"),
                      discrete = discrete,method = method
            )
        
        bam_mod
        
    }

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' Base+Stage+Sex: Base interaction, child stage main effect, and sex interaction across childhood
#' 
compute_bam_with_age_sex_effects <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T){
    
    base_form <- 
        deparse(E ~ s(nichd,
                      by=D,
                      bs=bs_,
                      k=length(stage_knots)))
    nichd_form <- 
        deparse(~s(nichd,
                   bs=bs_,
                   k=length(stage_knots))) %>% substring(2)
    sex_form <- 
        deparse(~s(nichd,
                   by=sex,
                   bs=bs_,
                   k=length(stage_knots)))%>% substring(2)
    
    form <- paste0(base_form," + ",nichd_form," + ",sex_form) %>% 
        as.formula()
    
    #bs="cs" gives lower AIC and df - a better fit model
    bam_mod <-
        mgcv::bam(form,
            data=dat_new,
            family=binomial(link="logit"),
            discrete = discrete,method = method
        )
    
    bam_mod
    
}

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' Base+Sex+NdrugsS: Base interaction, sex interaction across childhood, and the number of drugs
#' 
compute_bam_with_sex_ndrugsS_effects <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T){
        
        base_form <- 
            deparse(E ~ s(nichd,
                          by=D,
                          bs=bs_,
                          k=length(stage_knots)))
        
        ndrugs_form <- 
            deparse(~s(polypharmacy,
                       bs=bs_)) %>% substring(2)  
        sex_form <- 
            deparse(~s(nichd,
                       by=sex,
                       bs=bs_,
                       k=length(stage_knots)))%>% substring(2)
        
        form <- paste0(base_form," +",ndrugs_form," + ",sex_form) %>% 
            as.formula()
        
        #bs="cs" gives lower AIC and df - a better fit model
        bam_mod <-
            mgcv::bam(form,
                      data=dat_new,
                      family=binomial(link="logit"),
                      discrete = discrete,method = method
            )
        
        bam_mod
        
    }

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' Base+Sex+NdrugsL: Base interaction, sex interaction across childhood, and the number of drugs
#' 
compute_bam_with_sex_ndrugsL_effects <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T){
        
        base_form <- 
            deparse(E ~ s(nichd,
                          by=D,
                          bs=bs_,
                          k=length(stage_knots)))
        sex_form <- 
            deparse(~s(nichd,
                       by=sex,
                       bs=bs_,
                       k=length(stage_knots)))%>% substring(2)
        
        form <- paste0(base_form," + polypharmacy + ",sex_form) %>% 
            as.formula()
        
        #bs="cs" gives lower AIC and df - a better fit model
        bam_mod <-
            mgcv::bam(form,
                      data=dat_new,
                      family=binomial(link="logit"),
                      discrete = discrete,method = method
            )
        
        bam_mod
        
    }

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' Base+Sex+ATC: Base interaction, sex interaction across childhood, and anatomical/pharmacological number of drugs effect
#' 
compute_bam_with_sex_effects_atc <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T,id_col="safetyreportid"){
        
        dat_new <- merge(dat_new,atc_covariates,by=id_col)
        
        base_form <- 
            deparse(E ~ s(nichd,
                          by=D,
                          bs=bs_,
                          k=length(stage_knots)))
        sex_form <- 
            deparse(~s(nichd,
                       by=sex,
                       bs=bs_,
                       k=length(stage_knots)))%>% substring(2)
        
        cov_form <- paste0(colnames(atc_covariates)[2:ncol(atc_covariates)],collapse=" + ")
        
        form <- paste0(base_form," + ",sex_form,
                       " + ",cov_form) %>% 
            as.formula()
        
        #bs="cs" gives lower AIC and df - a better fit model
        bam_mod <-
            mgcv::bam(form,
                      data=dat_new,
                      family=binomial(link="logit"),
                      discrete = discrete,method = method
            )
        
        bam_mod
        
    }

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' Base+Sex+ATC3: Base interaction, sex interaction across childhood, and therapeutic number of drugs effect
#' 
compute_bam_with_sex_effects_atc3 <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T,id_col="safetyreportid"){
        
        dat_new <- merge(dat_new,atc_covariates_l3,by=id_col)
        
        base_form <- 
            deparse(E ~ s(nichd,
                          by=D,
                          bs=bs_,
                          k=length(stage_knots)))
        sex_form <- 
            deparse(~s(nichd,
                       by=sex,
                       bs=bs_,
                       k=length(stage_knots)))%>% substring(2)
        
        cov_form <- paste0(colnames(atc_covariates_l3)[2:ncol(atc_covariates_l3)],collapse=" + ")
        
        form <- paste0(base_form," + ",sex_form,
                       " + ",cov_form) %>% 
            as.formula()
        
        #bs="cs" gives lower AIC and df - a better fit model
        bam_mod <-
            mgcv::bam(form,
                      data=dat_new,
                      family=binomial(link="logit"),
                      discrete = discrete,method = method
            )
        
        bam_mod
        
    }

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' Base+Sex+NdrugsL+ATCbin: Base interaction, sex interaction across childhood, the mean of drugs, and anatomical/pharmacological drugs effect
#' 
compute_bam_with_sex_ndrugsL_effects_atcbin <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T,id_col="safetyreportid"){
        
        dat_new <- merge(dat_new,atc_covariates_bin,by=id_col)
        
        base_form <- 
            deparse(E ~ s(nichd,
                          by=D,
                          bs=bs_,
                          k=length(stage_knots)))
        
        sex_form <- 
            deparse(~s(nichd,
                       by=sex,
                       bs=bs_,
                       k=length(stage_knots)))%>% substring(2)
        
        cov_form <- paste0(colnames(atc_covariates_bin)[2:ncol(atc_covariates_bin)],collapse=" + ")
        
        form <- paste0(base_form," + polypharmacy + ",
                       sex_form," + ",cov_form) %>% 
            as.formula()
        
        #bs="cs" gives lower AIC and df - a better fit model
        bam_mod <-
            mgcv::bam(form,
                      data=dat_new,
                      family=binomial(link="logit"),
                      discrete = discrete,method = method
            )
        
        bam_mod
        
    }

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' Base+Sex+NdrugsL+ATC3bin: Base interaction, sex interaction across childhood, the mean of drugs, and anatomical/pharmacological drugs effect
#' 
compute_bam_with_sex_ndrugsL_effects_atc3bin <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T,id_col="safetyreportid"){
        
        dat_new <- merge(dat_new,atc_covariates_l3_bin,by=id_col)
        
        base_form <- 
            deparse(E ~ s(nichd,
                          by=D,
                          bs=bs_,
                          k=length(stage_knots)))
        
        sex_form <- 
            deparse(~s(nichd,
                       by=sex,
                       bs=bs_,
                       k=length(stage_knots)))%>% substring(2)
        
        cov_form <- paste0(colnames(atc_covariates_l3_bin)[2:ncol(atc_covariates_l3_bin)],collapse=" + ")
        
        form <- paste0(base_form," + polypharmacy + ",
                       sex_form," + ",cov_form) %>% 
            as.formula()
        
        #bs="cs" gives lower AIC and df - a better fit model
        bam_mod <-
            mgcv::bam(form,
                      data=dat_new,
                      family=binomial(link="logit"),
                      discrete = discrete,method = method
            )
        
        bam_mod
        
    }

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' Base+Sex+NdrugsS: Base interaction, sex interaction across childhood, the number of drugs, and anatomical/pharmacological drugs effect
#' 
compute_bam_with_sex_ndrugsS_effects_atcbin <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T,id_col="safetyreportid"){
        
        dat_new <- merge(dat_new,atc_covariates_bin,by=id_col)
        
        base_form <- 
            deparse(E ~ s(nichd,
                          by=D,
                          bs=bs_,
                          k=length(stage_knots)))
        
        ndrugs_form <- 
            deparse(~s(polypharmacy,
                       bs=bs_)) %>% substring(2)  
        
        sex_form <- 
            deparse(~s(nichd,
                       by=sex,
                       bs=bs_,
                       k=length(stage_knots)))%>% substring(2)
        
        cov_form <- paste0(colnames(atc_covariates_bin)[2:ncol(atc_covariates_bin)],collapse=" + ")
        
        form <- paste0(base_form," + ",ndrugs_form," + ",
                       sex_form," + ",cov_form) %>% 
            as.formula()
        
        #bs="cs" gives lower AIC and df - a better fit model
        bam_mod <-
            mgcv::bam(form,
                      data=dat_new,
                      family=binomial(link="logit"),
                      discrete = discrete,method = method
            )
        
        bam_mod
        
    }

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' Base+Sex+Date+Reporter+NdrugsS+ATC: Base interaction, sex interaction across childhood, the number of drugs, and anatomical/pharmacological drugs effect, reporter type and date of report
#' 
compute_bam_with_sex_ndrugsS_reporter_date_effects_atcbin <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T,id_col="safetyreportid"){
        
        dat_new <- merge(dat_new,atc_covariates_bin,by=id_col)
        
        base_form <- 
            deparse(E ~ s(nichd,
                          by=D,
                          bs=bs_,
                          k=length(stage_knots)))
        
        ndrugs_form <- 
            deparse(~s(polypharmacy,
                       bs=bs_)) %>% substring(2)  
        
        sex_form <- 
            deparse(~s(nichd,
                       by=sex,
                       bs=bs_,
                       k=length(stage_knots)))%>% substring(2)
        rq_form <- 
            deparse(
                ~te(reporter_qualification,
                    bs=bs_,
                    k=length(category_levels$reporter_qualification))) %>% substring(2)
        rd_form <- 
            deparse(
                ~s(receive_date,
                   bs=bs_)) %>% substring(2)
        
        cov_form <- paste0(colnames(atc_covariates_bin)[2:ncol(atc_covariates_bin)],collapse=" + ")
        
        form <- paste0(base_form," + ",ndrugs_form," + ",
                       sex_form," + ",rq_form," + ",rd_form,
                       " + ",cov_form) %>% 
            as.formula()
        
        #bs="cs" gives lower AIC and df - a better fit model
        bam_mod <-
            mgcv::bam(form,
                      data=dat_new,
                      family=binomial(link="logit"),
                      discrete = discrete,method = method
            )
        
        bam_mod
        
    }

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' Base+Sex+NdrugsS+ATC3: Base interaction, sex interaction across childhood, the number of drugs, and anatomical/pharmacological drugs effect
#' 
compute_bam_with_sex_ndrugsS_effects_atc3bin <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T,id_col="safetyreportid"){
        
        dat_new <- merge(dat_new,atc_covariates_l3_bin,by=id_col)
        
        base_form <- 
            deparse(E ~ s(nichd,
                          by=D,
                          bs=bs_,
                          k=length(stage_knots)))
        
        ndrugs_form <- 
            deparse(~s(polypharmacy,
                       bs=bs_)) %>% substring(2)  
        
        sex_form <- 
            deparse(~s(nichd,
                       by=sex,
                       bs=bs_,
                       k=length(stage_knots)))%>% substring(2)
        
        cov_form <- paste0(colnames(atc_covariates_l3_bin)[2:ncol(atc_covariates_l3_bin)],collapse=" + ")
        
        form <- paste0(base_form," + ",ndrugs_form," + ",
                       sex_form," + ",cov_form) %>% 
            as.formula()
        
        #bs="cs" gives lower AIC and df - a better fit model
        bam_mod <-
            mgcv::bam(form,
                      data=dat_new,
                      family=binomial(link="logit"),
                      discrete = discrete,method = method
            )
        
        bam_mod
        
    }

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' Base+Stage+Sex+ATC: Base interaction, child stage main effect, sex interaction across childhood, and anatomical/pharmacological number of drugs effect
#' 
compute_bam_with_age_sex_effects_atc <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T,id_col="safetyreportid"){
        
        dat_new <- merge(dat_new,atc_covariates,by=id_col)
    
    base_form <- 
        deparse(E ~ s(nichd,
                      by=D,
                      bs=bs_,
                      k=length(stage_knots)))
    sex_form <- 
        deparse(~s(nichd,
                   by=sex,
                   bs=bs_,
                   k=length(stage_knots)))%>% substring(2)
    nichd_form <- 
        deparse(~s(nichd,
                   bs=bs_,
                   k=length(stage_knots))) %>% substring(2)
    
    cov_form <- paste0(colnames(atc_covariates)[2:ncol(atc_covariates)],collapse=" + ")
    
    form <- paste0(base_form," + ",nichd_form," + ",sex_form,
                   " + ",cov_form) %>% 
        as.formula()
    
    #bs="cs" gives lower AIC and df - a better fit model
    bam_mod <-
        mgcv::bam(form,
            data=dat_new,
            family=binomial(link="logit"),
            discrete = discrete,method = method
        )
    
    bam_mod
    
}
#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' Base+Stage+Sex+ATC3: Base interaction, child stage main effect, sex interaction across childhood, and therapeutic number of drugs effect
#' 
compute_bam_with_age_sex_effects_atc3 <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T,id_col="safetyreportid"){
        
        dat_new <- merge(dat_new,atc_covariates_l3,by=id_col)
        
        base_form <- 
            deparse(E ~ s(nichd,
                          by=D,
                          bs=bs_,
                          k=length(stage_knots)))
        sex_form <- 
            deparse(~s(nichd,
                       by=sex,
                       bs=bs_,
                       k=length(stage_knots)))%>% substring(2)
        nichd_form <- 
            deparse(~s(nichd,
                       bs=bs_,
                       k=length(stage_knots))) %>% substring(2)
        
        cov_form <- paste0(colnames(atc_covariates_l3)[2:ncol(atc_covariates_l3)],collapse=" + ")
        
        form <- paste0(base_form," + ",nichd_form," + ",sex_form,
                       " + ",cov_form) %>% 
            as.formula()
        
        #bs="cs" gives lower AIC and df - a better fit model
        bam_mod <-
            mgcv::bam(form,
                      data=dat_new,
                      family=binomial(link="logit"),
                      discrete = discrete,method = method
            )
        
        bam_mod
        
    }

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' Base+Stage+ATC+Date+Reporter: Base interaction, child stage main effect, anatomical/pharmacological number of drugs effect, reporter type effect, and reporting date effect
#' 
compute_bam_with_reporter_date_age_effects_atc <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T,id_col="safetyreportid"){
        
        dat_new <- merge(dat_new,atc_covariates,by=id_col)
        
    base_form <- 
        deparse(E ~ s(nichd,
                      by=D,
                      bs=bs_,
                      k=length(stage_knots)))
    nichd_form <- 
        deparse(~s(nichd,
                   bs=bs_,
                   k=length(stage_knots))) %>% substring(2)
    rq_form <- 
        deparse(
            ~te(reporter_qualification,
                bs=bs_,
                k=length(category_levels$reporter_qualification))) %>% substring(2)
    rd_form <- 
        deparse(
            ~s(receive_date,
               bs=bs_)) %>% substring(2)
    
    cov_form <- paste0(colnames(atc_covariates)[2:ncol(atc_covariates)],collapse=" + ")
    
    form <- paste0(base_form," + ",nichd_form," + ",rq_form," + ",rd_form," + ",cov_form) %>% 
        as.formula()
    
    #bs="cs" gives lower AIC and df - a better fit model
    bam_mod <-
        mgcv::bam(form,
            data=dat_new,
            family=binomial(link="logit"),
            discrete=discrete,method = method
        )
    
    bam_mod
    
}

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' Base+Sex+ATC+Date+Reporter: Base interaction, sex interaction across childhood, anatomical/pharmacological number of drugs effect, reporter type effect, and reporting date effect
#' 
compute_bam_with_reporter_date_sex_effects_atc <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T,id_col="safetyreportid"){
        
        dat_new <- merge(dat_new,atc_covariates,by=id_col)
        
        base_form <- 
            deparse(E ~ s(nichd,
                          by=D,
                          bs=bs_,
                          k=length(stage_knots)))
        sex_form <- 
            deparse(~s(nichd,
                       by=sex,
                       bs=bs_,
                       k=length(stage_knots)))%>% substring(2)
        rq_form <- 
            deparse(
                ~te(reporter_qualification,
                    bs=bs_,
                    k=length(category_levels$reporter_qualification))) %>% substring(2)
        rd_form <- 
            deparse(
                ~s(receive_date,
                   bs=bs_)) %>% substring(2)
        
        cov_form <- paste0(colnames(atc_covariates)[2:ncol(atc_covariates)],collapse=" + ")
        
        form <- paste0(base_form," + "," + ",sex_form,
                       " + ",rq_form," + ",rd_form," + ",cov_form) %>% 
            as.formula()
        
        #bs="cs" gives lower AIC and df - a better fit model
        bam_mod <-
            mgcv::bam(form,
                      data=dat_new,
                      family=binomial(link="logit"),
                      discrete = discrete,method = method
            )
        
        bam_mod
        
    }

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' Base+Sex+ATC3+Date+Reporter: Base interaction, sex interaction across childhood, therapeutic number of drugs effect, reporter type effect, and reporting date effect
#' 
compute_bam_with_reporter_date_sex_effects_atc3 <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T,id_col="safetyreportid"){
        
        dat_new <- merge(dat_new,atc_covariates_l3,by=id_col)
        
        base_form <- 
            deparse(E ~ s(nichd,
                          by=D,
                          bs=bs_,
                          k=length(stage_knots)))
        sex_form <- 
            deparse(~s(nichd,
                       by=sex,
                       bs=bs_,
                       k=length(stage_knots)))%>% substring(2)
        rq_form <- 
            deparse(
                ~te(reporter_qualification,
                    bs=bs_,
                    k=length(category_levels$reporter_qualification))) %>% substring(2)
        rd_form <- 
            deparse(
                ~s(receive_date,
                   bs=bs_)) %>% substring(2)
        
        cov_form <- paste0(colnames(atc_covariates_l3)[2:ncol(atc_covariates_l3)],collapse=" + ")
        
        form <- paste0(base_form," + "," + ",sex_form,
                       " + ",rq_form," + ",rd_form," + ",cov_form) %>% 
            as.formula()
        
        #bs="cs" gives lower AIC and df - a better fit model
        bam_mod <-
            mgcv::bam(form,
                      data=dat_new,
                      family=binomial(link="logit"),
                      discrete = discrete,method = method
            )
        
        bam_mod
        
    }

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' Base+Stage+Sex+ATC+Date+Reporter: Base interaction, child stage main effect, sex interaction across childhood, anatomical/pharmacological number of drugs effect, reporter type effect, and reporting date effect
#' 
compute_bam_with_reporter_date_age_sex_effects_atc <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T,id_col="safetyreportid"){
        
        dat_new <- merge(dat_new,atc_covariates,by=id_col)
        
    base_form <- 
        deparse(E ~ s(nichd,
                      by=D,
                      bs=bs_,
                      k=length(stage_knots)))
    sex_form <- 
        deparse(~s(nichd,
                   by=sex,
                   bs=bs_,
                   k=length(stage_knots)))%>% substring(2)
    nichd_form <- 
        deparse(~s(nichd,
                   bs=bs_,
                   k=length(stage_knots))) %>% substring(2)
    rq_form <- 
        deparse(
            ~te(reporter_qualification,
                bs=bs_,
                k=length(category_levels$reporter_qualification))) %>% substring(2)
    rd_form <- 
        deparse(
            ~s(receive_date,
               bs=bs_)) %>% substring(2)
    
    cov_form <- paste0(colnames(atc_covariates)[2:ncol(atc_covariates)],collapse=" + ")
    
    form <- paste0(base_form," + ",nichd_form," + ",sex_form,
                   " + ",rq_form," + ",rd_form," + ",cov_form) %>% 
        as.formula()
    
    #bs="cs" gives lower AIC and df - a better fit model
    bam_mod <-
        mgcv::bam(form,
            data=dat_new,
            family=binomial(link="logit"),
            discrete = discrete,method = method
        )
    
    bam_mod
    
}

#' Compute GAM
#' 
#' Predicting adverse event occurrence
#' 
#' Base+Stage+Sex+ATC3+Date+Reporter: Base interaction, child stage main effect, sex interaction across childhood, therapeutic number of drugs effect, reporter type effect, and reporting date effect
#' 
compute_bam_with_reporter_date_age_sex_effects_atc3 <- 
    function(dat_new,bs_="cs",method="fREML",discrete=T,id_col="safetyreportid"){
        
        dat_new <- merge(dat_new,atc_covariates_l3,by=id_col)
        
        base_form <- 
            deparse(E ~ s(nichd,
                          by=D,
                          bs=bs_,
                          k=length(stage_knots)))
        sex_form <- 
            deparse(~s(nichd,
                       by=sex,
                       bs=bs_,
                       k=length(stage_knots)))%>% substring(2)
        nichd_form <- 
            deparse(~s(nichd,
                       bs=bs_,
                       k=length(stage_knots))) %>% substring(2)
        rq_form <- 
            deparse(
                ~te(reporter_qualification,
                    bs=bs_,
                    k=length(category_levels$reporter_qualification))) %>% substring(2)
        rd_form <- 
            deparse(
                ~s(receive_date,
                   bs=bs_)) %>% substring(2)
        
        cov_form <- paste0(colnames(atc_covariates_l3)[2:ncol(atc_covariates_l3)],collapse=" + ")
        
        form <- paste0(base_form," + ",nichd_form," + ",sex_form,
                       " + ",rq_form," + ",rd_form," + ",cov_form) %>% 
            as.formula()
        
        #bs="cs" gives lower AIC and df - a better fit model
        bam_mod <-
            mgcv::bam(form,
                      data=dat_new,
                      family=binomial(link="logit"),
                      discrete = discrete,method = method
            )
        
        bam_mod
        
    }

#' Process GAM model data
#' 
#' Construct datatable with Base interaction only GAM parameter statistics and drug-event information
#' 
process_gam_model <- 
    function(bam_mod,dat_new,drug_id,rx_id,drug_name,rx_name,knots,coef_str="s(nichd):D.",time_=7,pad="",id_col="safetyreportid",e_col="E",d_col="D",de_col="DE"){
    
    # Extract GAM coefficients
    # mapping spline knots to category names
    # add overall GAM score - weighted average across stages
    #https://www.andrewheiss.com/blog/2016/04/25/convert-logistic-regression-standard-errors-to-odds-ratios-with-r/
    #https://stats.stackexchange.com/questions/364568/extracting-significance-of-gam-smoothing-terms
    
    cnames=names(coef(bam_mod))
    inds = grepl(coef_str,cnames,fixed=T)
    cknots = str_replace(cnames[inds],fixed(coef_str),"") %>% as.integer()
    scores_gam_dt <- data.table(stage = knots)
    scores_gam_dt[match(cknots,scores_gam_dt$stage),"gam_score"] <- coef(bam_mod,complete=T)[inds]
    scores_gam_dt[match(cknots,scores_gam_dt$stage),"gam_score_se"] <- summary(bam_mod,freq=F)$se[inds]
    scores_gam_dt$N <- dat_new[,.N,by=stage_col]$N
    setnames(scores_gam_dt,c("stage"),c(stage_col))
    scores_gam_dt[[stage_col]] <- as.factor(category_levels[[stage_col]])
    
    # Mapping GAM data to categories
    # for summarizing column counts
    dat_new[,(stage_col)] <-
        dat_new[,lapply(.SD,function(x){category_levels[[stage_col]][x]}),.SDcols=stage_col]
    dat_agg <-
        dat_new[,
                lapply(.SD,function(x){sum(x==1)}),
                by=stage_col,
                .SDcols=c(e_col,d_col,de_col)
        ] %>%
        pivot_longer(-!!as.name(stage_col)) %>%
        data.table()
    dat_agg[[stage_col]] <- factor(dat_agg[[stage_col]] %>% unlist %>% unname,
                                   levels=category_levels[[stage_col]])
    
    # Computing confidence (credible) interval scores
    scores_gam_dt$gam_score_90mse <-
        scores_gam_dt[,.(gam_score_90mse = gam_score - (1.645*(gam_score_se)))]$gam_score_90mse
    scores_gam_dt$gam_score_95mse <-
        scores_gam_dt[,.(gam_score_95mse = gam_score - (1.96*(gam_score_se)))]$gam_score_95mse
    scores_gam_dt$gam_score_99mse <-
        scores_gam_dt[,.(gam_score_99mse = gam_score - (2.576*(gam_score_se)))]$gam_score_99mse
    scores_gam_dt$gam_score_999mse <-
        scores_gam_dt[,.(gam_score_999mse = gam_score - (3.291*(gam_score_se)))]$gam_score_999mse
    scores_gam_dt$gam_score_90pse <-
        scores_gam_dt[,.(gam_score_90pse = gam_score + (1.645*(gam_score_se)))]$gam_score_90pse
    scores_gam_dt$gam_score_95pse <-
        scores_gam_dt[,.(gam_score_95pse = gam_score + (1.96*(gam_score_se)))]$gam_score_95pse
    scores_gam_dt$gam_score_99pse <-
        scores_gam_dt[,.(gam_score_99pse = gam_score + (2.576*(gam_score_se)))]$gam_score_99pse
    scores_gam_dt$gam_score_999pse <-
        scores_gam_dt[,.(gam_score_999pse = gam_score + (3.291*(gam_score_se)))]$gam_score_999pse
    
    stage_dat_agg <- dat_agg[,.(value=sum(value,na.rm=T)),by=c(stage_col,"name")] %>% dcast(nichd ~ name)
    colnames(stage_dat_agg)[1] <- stage_col
    all_dat_agg <- dat_agg[,.(value=sum(value,na.rm=T)),by=c("name")]
    all_dat_agg[[stage_col]] <- "all"
    all_dat_agg <- all_dat_agg %>% dcast(nichd ~ name)
    colnames(all_dat_agg)[1] <- stage_col
    dat_agg <- bind_rows(stage_dat_agg,all_dat_agg)
    # Joining mapped GAM data
    # to GAM results
    joined <- inner_join(scores_gam_dt,dat_agg,by=stage_col) %>% data.table()
    summ <- summary(bam_mod)
    joined[,drug_col] = drug_id
    joined[,rx_col] = rx_id
    joined[,"ade"] = paste0(drug_id,"_",rx_id)
    joined[,drug_col_name] = drug_name
    joined[,rx_col_name] = rx_name
    joined[,"ade_name"] = paste0(drug_name," and ",rx_name)
    joined[,"formula"] = paste0(as.character(deparse(bam_mod$formula)),collapse="")
    joined[,"AIC"] = bam_mod$aic
    joined[,"gcv_ubre"] = bam_mod$gcv.ubre %>% unname
    
    y_true = bam_mod$y
    y_pred = bam_mod$fitted.values
    joined[,"mean_true_probability"] = mean(y_pred[y_true==1])
    joined[,"mean_false_probability"] = mean(y_pred[y_true==0])
    
    if(length(unique(y_true))==2){
        pred = ROCR::prediction(y_pred,y_true)
        pred = ROCR::performance(prediction.obj = pred,"auc")
        auc = pred@y.values[[1]]
        joined[,"auc"] = auc        
    }else{
        joined[,"auc"] = NA
    }
    
    joined[,"time"] = time_
    
    stringi::stri_sub(colnames(joined)[grepl("gam",colnames(joined))],4,3) <- pad
    joined
}

#' Extract GAM statistics
#' 
#' Add gam statistics to a data.table
#' 
extract_gam_statistics <- 
    function(mod,name,dt,time_=""){
        
        ss <- summary(mod)
        y_true = mod$y
        y_pred = mod$fitted.values
        
        dt[,"model"] = name
        dt[,"model_time"] = time_
        dt[,"formula"] = mod$formula %>% deparse() %>% paste0(collapse="")
        dt[,"model_gcv_ubre"] = mod$gcv.ubre
        dt[,"model_nparameters"] = ss$np-1
        dt[,"model_deviance"] = anova(mod)$dev.expl
        dt[,"model_rsq"] = ss$r.sq
        dt[,"model_aic"] = mod$aic
        dt[,"model_mean_probability"] = mean(y_pred,na.rm=T)
        dt[,"model_mean_true_probability"] = mean(y_pred[y_true==1],na.rm=T)
        dt[,"model_mean_false_probability"] = mean(y_pred[y_true==0],na.rm=T)
        
        dt
    }
#' Extract all GAM coefficients 
#'
#' Construct datatable with all GAM parameter statistics and drug-event information
#' 
extract_gam_coefficients <- 
    function(bam_mod,drug_id,rx_id,drug_name,rx_name){
    
    scores_gam_dt <- data.table(predictor = names(coefficients(bam_mod)),
                                gam_score = coefficients(bam_mod),
                                gam_score_se = diag(vcov(bam_mod)))
    
    scores_gam_dt$gam_score_90mse <-
        scores_gam_dt[,.(gam_score_90mse = gam_score - (1.645*(gam_score_se)))]$gam_score_90mse
    scores_gam_dt$gam_score_95mse <-
        scores_gam_dt[,.(gam_score_95mse = gam_score - (1.96*(gam_score_se)))]$gam_score_95mse
    scores_gam_dt$gam_score_99mse <-
        scores_gam_dt[,.(gam_score_99mse = gam_score - (2.576*(gam_score_se)))]$gam_score_99mse
    scores_gam_dt$gam_score_999mse <-
        scores_gam_dt[,.(gam_score_999mse = gam_score - (3.291*(gam_score_se)))]$gam_score_999mse
    scores_gam_dt$gam_score_90pse <-
        scores_gam_dt[,.(gam_score_90pse = gam_score + (1.645*(gam_score_se)))]$gam_score_90pse
    scores_gam_dt$gam_score_95pse <-
        scores_gam_dt[,.(gam_score_95pse = gam_score + (1.96*(gam_score_se)))]$gam_score_95pse
    scores_gam_dt$gam_score_99pse <-
        scores_gam_dt[,.(gam_score_99pse = gam_score + (2.576*(gam_score_se)))]$gam_score_99pse
    scores_gam_dt$gam_score_999pse <-
        scores_gam_dt[,.(gam_score_999pse = gam_score + (3.291*(gam_score_se)))]$gam_score_999pse
        
        
    scores_gam_dt[,drug_col] = drug_id
    scores_gam_dt[,rx_col] = rx_id
    scores_gam_dt[,drug_col_name] = drug_name
    scores_gam_dt[,rx_col_name] = rx_name
    
    scores_gam_dt
    
}

#' Min-max normalize GAM parameter statistics
#' 
#' Min-max normalize score values across strata with optional groups
#' 
normalize_data <- 
    function(dat,obs="ade",strata="nichd",score="gam_score",groups="database"){
        
        normalize <- function(x) {
            return ((x - min(x)) / (max(x) - min(x)))
        }
        
        sub <- 
            dat[,
                c(obs,strata,score,groups),
                with=F
            ] %>% 
            unique() %>% 
            .[order(get(strata))]
        
        sub[,norm := normalize(get(score)),c(obs,groups)]
        
        sub
    }

#' Make data from tSNE
#' 
make_tsne_data <- function(X,ids,subset=100,ndims=2,perplexity=3,check_duplicates=F,dimname="tsne",pca=T,partial_pca=F,max_iter=500,num_threads=0,seed=0,theta=1){
    
    set.seed(seed)
    
    inds = sample(1:nrow(X),min(nrow(X),subset),replace=F)
    
    tsne <- Rtsne(as.matrix(X[inds]),initial_dims=7,
                  check_duplicates = check_duplicates,
                  perplexity = perplexity,
                  theta=theta,
                  pca=pca,
                  partial_pca=partial_pca,
                  max_itrer=max_iter,
                  dims = ndims,
                  num_threads=num_threads)
    
    tsne_dt <- data.table(tsne$Y) %>% cbind(ades[inds])
    
    setnames(tsne_dt,paste0("V",seq(1,ndims)),paste0(dimname,seq(1,ndims)))
    
    tsne_dt$kl_divergence <- tsne$itercosts[length(tsne$itercosts)]
    
    tsne_dt
    
}

#' Make data from UMAP
#'
make_umap_data <- 
    function(X,ades,ndims=2,seed=0,
             subset=1000,metric="cosine",spread=spread,pca=T,a=NA,b=NA){
        
        set.seed(seed)
        
        inds = sample(1:nrow(X),min(nrow(X),subset),replace=F)
        
        params <- umap.defaults
        params$metric <- metric
        params$n_components <- ndims
        params$random_state <- seed
        params$transform_state <- seed
        params$spread <- spread
        params$a <- a
        params$b <- b
        
        if(pca){
            X <- prcomp(X[inds],scale=F,retx=T)$x
        }else{
            X <- X[inds]
        }
        
        datum.umap <- 
            umap(X,config=params)
        
        dt <- data.table(datum.umap$layout) %>% 
            cbind(ades[inds]) 
        
        setnames(dt,paste0("V",seq_along(1:ndims)),paste0("umap",seq_along(1:ndims)))
        
    }

#' Plot drug-event risk scores across childhood
#' 
plot_drug_events <- 
    function(dts,color_name="database",wrap_gen_width=20){
        dts %>% 
        pivot_longer(cols=c("D","E","DE")) %>% 
            filter(value!=0) %>% 
            mutate(
                value = log10(value+1)
            ) %>% 
            ggplot() +
            geom_bar(aes(nichd,value,fill=name),
                     stat="identity",position="dodge") +
            geom_errorbar(aes(x=nichd,y=gam_score,
                              ymin=gam_score_90mse,
                              ymax=gam_score_90pse,
                              color=!!as.name(color_name)),
                          position = position_dodge(width=0.5),
                          width=0.2,size=1) +
            geom_point(aes(nichd,gam_score,
                           color=!!as.name(color_name)),
                       position = position_dodge(width=0.5),size=2) +
            scale_y_continuous(sec.axis=sec_axis(~.,name="Risk of ADE\n(GAM score)")) + 
            scale_color_brewer(palette = "Dark2") +
            facet_wrap(~ade_name,labeller = label_wrap_gen(width = wrap_gen_width),scales="free_y") +
            theme(
                axis.text.x = element_text(angle=45,vjust=1,hjust=1),
                strip.background = element_blank(),
                legend.position = "bottom"
            ) +
            xlab("") +
            ylab("log10(Number of Reports)")
    }

#' Evaluate patterns from clustering with hyperparameter sets
#' 
#' Input: eval_dt
#' Output: cluster metric, dynamics metric, and cluster*dynamics metric
evaluate_patterns_in_parameter_sets <- function(dt){
    
    # remove drug-events that are not canonical from the parameter results
    dt <- dt %>% .[spikein!=""]
    
    # count number of canonical drug-events per dynamic
    # should be the same for each parameter set
    dynamic_dt <- 
        dt[,
           .(param_grid_ind,spikein,ade)
        ] %>% 
        unique() %>% 
        .[,.(total_dynamic_ades = .N),.(param_grid_ind,spikein)]
    
    # count number of canonical drug-events in each cluster
    cluster_dt <- 
        dt[,
           .(param_grid_ind,cluster,ade)
        ] %>% 
        unique() %>% 
        .[,.(total_cluster_ades = .N),.(param_grid_ind,cluster)]
    
    # count number of canonical drug-events per dynamic in clusters
    # calculate fraction across clusters and across dynamics
    cluster_dynamic_dt <- 
        dt[,
           .(param_grid_ind,cluster,spikein,ade)
        ] %>% 
        unique() %>% 
        .[,.N,.(param_grid_ind,cluster,spikein)] %>%
        merge(dynamic_dt,
              by=c("param_grid_ind","spikein")
        ) %>% 
        merge(cluster_dt,
              by=c("param_grid_ind","cluster")
        ) %>% 
        .[,.(param_grid_ind,spikein,cluster,N,
             frac_cluster = N/total_cluster_ades,
             frac_dynamic = N/total_dynamic_ades)]
    
    clusters_for_purity <- 
        cluster_dynamic_dt %>% 
        .[order(frac_cluster,decreasing = T)] %>% 
        .[,.SD[1],.(param_grid_ind,spikein)] %>% 
        .[,.(param_grid_ind,spikein,cluster)] %>% 
        unique()
    
    tmp = merge(
        cluster_dynamic_dt %>% 
            .[,.(max_frac_cluster = max(frac_cluster)),.(param_grid_ind,spikein)] %>% 
            .[,.(mean_max_frac_cluster = mean(max_frac_cluster)),.(param_grid_ind)],
        cluster_dynamic_dt %>% 
            .[,.(max_frac_dynamic = max(frac_dynamic)),.(param_grid_ind,cluster)] %>% 
            merge(clusters_for_purity,by=c("param_grid_ind","cluster")) %>% 
            .[,.(mean_max_frac_dynamic = mean(max_frac_dynamic)),.(param_grid_ind)],
        by="param_grid_ind"
    )
    
    tmp %>% 
        .[,.(param_grid_ind,
             mean_max_frac_dynamic,
             mean_max_frac_cluster,
             opt = mean_max_frac_dynamic * mean_max_frac_cluster)] %>% 
        .[order(opt,decreasing = T)]
}

fcn_lst = 
    list(
        "Base"= compute_bam_base,
        "Base+Sex" = compute_bam_with_sex_effects,
        "Base+Date+Reporter" = compute_bam_with_reporter_date_effects,
        "Base+Sex+ATC" = compute_bam_with_sex_effects_atc,
        "Base+Sex+ATC3" = compute_bam_with_sex_effects_atc3,
        "Base+Sex+NdrugsS" = compute_bam_with_sex_ndrugsS_effects,
        "Base+Sex+NdrugsS+ATCbin" = compute_bam_with_sex_ndrugsS_effects_atcbin,
        "Base+Sex+Date+Reporter+NdrugsS+ATCbin" = 
            compute_bam_with_sex_ndrugsS_reporter_date_effects_atcbin,
        "Base+Sex+Date+Reporter+ATC" = compute_bam_with_reporter_date_sex_effects_atc
    )
