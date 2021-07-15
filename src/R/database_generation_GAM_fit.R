
#' Title: "A database of pediatric drug effects to evaluate ontogenic mechanisms from child growth and development" study
#' 
#' Author details: Nicholas Giangreco
#' 
#' This script generates the GAM fit statistics for a random drug-event sample

# Purpose -----------------------------------------------------------------

#' To evaluate the fit and output of GAMs from the base specification to adding covariates
#' 


# Setup -------------------------------------------------------------------

pacman::p_load(data.table,tidyverse,doParallel,mgcv)

data_dir <- "../../data/"
img_dir <- "../../docs/imgs/"
seed = 0
set.seed(seed)

method="fREML"
discrete=T

registerDoParallel(cores=20)

basename <- paste0("database_generation_fit")
out_dir <- paste0(data_dir,basename,"/")

ifelse(
    !dir.exists(file.path(out_dir)),
    dir.create(file.path(out_dir)),
    FALSE)

theme_set(theme_bw(base_size=16) + theme(text = element_text(face="bold")))

# Load data and covariates ------------------------------------------------

source("database_generation_load_data.R")

stage_knots <- 
    seq_along(category_levels[[stage_col]])

# functions ---------------------------------------------------------------

source("database_generation_functions.R")

# candidate ADEs ----------------------------------------------------------

candidate_ades <- 
    fread(paste0(data_dir,"database_generation_random_ade_candidates_for_fits.csv"))


# Main --------------------------------------------------------------------

x_cols=c(names(category_levels),"polypharmacy")
ade_data <- 
    make_ade_data(
        x_cols=x_cols,
        category_levels = category_levels
    )

dts <- NULL
minrange = 1
maxrange = nrow(candidate_ades)
cat(maxrange," candidates\n")
t0 = Sys.time()
cat("Start:",format(Sys.time(), "%a %b %d %X %Y"),"\n")

range <- minrange:maxrange
step=0.02
start_cuts <-  seq(step,1-step,step)
end_cuts <-  seq(step,1,step)
starts <- sapply(start_cuts,function(x)range[floor(x*length(range))])
starts <- c(min(range),starts)
ends <- sapply(end_cuts,function(x)range[floor(x*length(range))])
for(i in seq(end_cuts)){
    
    start <- ifelse(starts[i]==min(range),starts[i],starts[i]+1)
    end <- ends[i]
    cat("Indices: ",start," to ",end,"\n")
        
    dt <-  foreach(x=start:end,.combine="rbind",.errorhandling='remove') %dopar% {
            
            row=candidate_ades[x]
            drug_id = row[,..drug_col][[1]]
            rx_id = row[,..rx_col][[1]]
            drug_name = row[,..drug_col_name][[1]]
            rx_name = row[,..rx_col_name][[1]]
            
            dat <- 
                set_ade_pair(ade_data,drug_id,rx_id,drug_name,rx_name)
            
            mod_dts <- NULL
            mod_lst <- list()
            for(j in seq_along(fcn_lst)){
                time_ <- system.time(
                    mod <- fcn_lst[[j]](dat,method=method,discrete=discrete)
                )
                mod_lst[[names(fcn_lst)[j]]] <- mod
                ss <- summary(mod)
                y_true = mod$y
                y_pred = mod$fitted.values
                
                dt <- broom::tidy(mod) %>% data.table()
                
                dt <- 
                    extract_gam_statistics(mod,names(fcn_lst)[j],dt,time_[3])
                
                mod_dts <- rbind(mod_dts,dt)
            }
            
            av_dts <- NULL
            for(ii in seq_along(mod_lst)){
                for(ij in seq_along(mod_lst)){
                    if(ii!=ij & ij>ii){
                        av <- anova(mod_lst[[ii]],mod_lst[[ij]],test="Chisq")
                        comp <- paste0("anova(",names(mod_lst)[ii],",",names(mod_lst)[ij],")")
                        dt <- data.table(
                            model_ground = names(mod_lst)[ii],
                            model_comparator = names(mod_lst)[ij],
                            anova_comparison = comp,
                            anova_deviance = av$Deviance[2],
                            anova_pvalue = av$`Pr(>Chi)`[2]
                            )
                        av_dts <- rbind(av_dts,dt)
                    }
                }
            }
            
            merged_dt <- 
                merge(
                    mod_dts,
                    av_dts,
                    by.x="model",
                    by.y="model_comparator",
                    allow.cartesian = T,
                    all.x=T
                    )
            
            merged_dt[,"ade"] <- paste0(drug_id,"_",rx_id)
            merged_dt[,drug_col] <- drug_id
            merged_dt[,rx_col] <- rx_id
            merged_dt[,drug_col_name] <- drug_name
            merged_dt[,rx_col_name] <- rx_name
            
            merged_dt
            
        }
        
    dt %>% 
        fwrite(paste0(out_dir,basename,"_",method,"_",start,"to",end,".csv"))
    
    dts <- bind_rows(dts,dt)
    
    cat((i/length(end_cuts))*100,"% complete\n",format(Sys.time(), "%a %b %d %X %Y"),"\n")
        
}
t1 = Sys.time()
cat("End:",format(Sys.time(), "%a %b %d %X %Y"),"\n")
cat("\n",round(as.numeric(difftime(t1,t0,units="mins")),2)," minutes\n")

dts %>% 
    fwrite(paste0(data_dir,basename,"_",method,".csv"))
