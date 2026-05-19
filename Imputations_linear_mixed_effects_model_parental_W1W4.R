
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#   Authors: M.H., Y.R. #
#   Imputations and mixed-effect models 
#   for parental-reported parenting behaviors over time (4 waves)
#   Code for adolescents-reports looks similar with 3 waves
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


# To install:
install.packages("devtools")
devtools::install_github(repo = "ryannick28/CustomFunctionsYrotha", ref = 'main') 

### Libraries:
library(lmerTest)
library(car)
library(mice)
library(nlme)
library(broom.mixed)
library(CustomFunctionsYrotha) # We use the function wideToLong: converts wide data format to long and prepares data for repeated-measures analyses / mixed models


### Define paths:
path_data <- ".csv" # define path to the data location 
path_save <- ".rds" # define path to save imputation outputs

# Prepare wide data
# Read in data (with missings)

dat_m0 <- read.csv(path_data, stringsAsFactors = TRUE)

### Take a selection:
# P means parent-reported; number reflects assessment wave
dat_mw <- dat_m0[,c("P1_Involvment_mean","P2_Involvment_mean","P3_Involvment_mean","P4_Involvment_mean",
                    "P1_Positive_Parenting_mean", "P2_Positive_Parenting_mean", 
                    "P3_Positive_Parenting_mean", "P4_Positive_Parenting_mean",
                    "P1_Poor_Monitoring_mean","P2_Poor_Monitoring_mean",
                    "P3_Poor_Monitoring_mean", "P4_Poor_Monitoring_mean",
                    "P1_Inconsistent_Discipline_mean", "P2_Inconsistent_Discipline_mean",
                    "P3_Inconsistent_Discipline_mean", "P4_Inconsistent_Discipline_mean",
                    "P1_Corporal_Punishment_mean", "P2_Corporal_Punishment_mean", 
                    "P3_Corporal_Punishment_mean", "P4_Corporal_Punishment_mean", 
                    "ADR_GenderChild")]

# Add ID variable
dat_mw$ID <- factor(paste0('S_', 1:nrow(dat_mw)))
# Make sure gender variable is a factor
dat_mw$ADR_GenderChild <- as.factor(dat_mw$ADR_GenderChild)

# Run MICE Imputation

# Set seed for reproducibility
set.seed(123) 

# Nr of imputations
m <- 30

# Mice imputation:
m.out <- mice(dat_mw, m = m, pred=quickpred(dat_mw)) # quickpred: predictor matrix with default values of r=.1 and cases of0.25
saveRDS(m.out, path_save)

# --- DIAGNOSTICS ---
plot(m.out) 
densityplot(m.out)

# FUNCTION

compute_mixed_models_T1_T4 <- function(m_out_obj, PARENTING_BEH, PARENTING_save_name, log ="",pot="", out_dir_main = "results/linear_mixed_models", out_dir_seed = "results/linear_mixed_models/seed_123", out_dir_formatted = "results/linear_mixed_models/formatted") {

  m <- m_out_obj$m  # Extract m from the object 
  
  # Collect completed datasets in list:
  comp.out <- vector('list', length = m) #create empty list, that will be filled with outputs
  for(i in 1:m){
    comp.out[[i]] <- complete(m_out_obj, action = i)
  }

  
  ### Turn each completed data set to long format & change reference of factor:
  comp.out_l <- lapply(comp.out, FUN = function(x){
    d_l <-  wideToLong(x, ind = 'P*_', nRep=4, repColnm = 'duration')
    ### Change reference:
    d_l$duration <- factor(d_l$duration, levels=c("P1_", "P2_", "P3_", "P4_"))
    return(d_l)
  })
  # Fit mixed model (optim needed because of convergence issues)
  imp_fit <<- lapply(comp.out_l, function(x) {
    lme(fixed = as.formula(paste(PARENTING_BEH, "~ ADR_GenderChild + duration")), # ADR_GenderChild = Sex
        data = x, random = ~1 | ID, control = lmeControl(opt = 'optim'))
  })
  ### Pool results from 30 imputations:
  summary(pool(imp_fit)) ##
  t1 <- summary(pool(imp_fit))
  t1 <- t1[1:3,] # keep relevant effects

  ## repeat above analyses with changing factor levels (P2_ =reference)
  ### Collect completed datasets in list:
  comp.out <- vector('list', length = m) #create empty list, that will be filled with outputs
  for(i in 1:m){
    comp.out[[i]] <- complete(m_out_obj, action = i)
  }
  comp.out_l <- lapply(comp.out, FUN = function(x){
    d_l <-  wideToLong(x, ind = 'P*_', nRep=4, repColnm = 'duration')
    ### Change reference:
    d_l$duration <- factor(d_l$duration, levels=c("P2_", "P3_", "P4_", "P1_"))
    return(d_l)
  })

  ### Fit mixed model (optim needed because of convergence issues)
  imp_fit <- lapply(comp.out_l, function(x){
    lme(fixed = as.formula(paste(PARENTING_BEH, "~ ADR_GenderChild + duration")),
        data = x, random = ~1 | ID, control = lmeControl(opt = 'optim'))
  })
  ## pool and add relevant effect to output table
  t2 <- summary(pool(imp_fit))
  t1 <- rbind(t1,t2[3,])
  
  ## repeat above analyses with changing factor levels (P3_ =reference)
  comp.out <- vector('list', length = m) #create empty list, that will be filled with outputs
  for(i in 1:m){
    comp.out[[i]] <- complete(m_out_obj, action = i)
  }
  comp.out_l <- lapply(comp.out, FUN = function(x){
    d_l <-  wideToLong(x, ind = 'P*_', nRep=4, repColnm = 'duration')
    ### Change reference:
    d_l$duration <- factor(d_l$duration, levels=c("P3_", "P4_", "P1_", "P2_"))
    return(d_l)
  })

  ### Fit mixed model (optim needed because of convergence issues)
  imp_fit <- lapply(comp.out_l, function(x){
    lme(fixed = as.formula(paste(PARENTING_BEH, "~ ADR_GenderChild + duration")),
        data = x, random = ~1 | ID, control = lmeControl(opt = 'optim'))
  })
  ## pool and add relevant effect to output table
  t2 <- summary(pool(imp_fit))
  t1 <- rbind(t1,t2[3,])
  t1$term <- c("T1", "Sex", "T1-T2", "T2-T3", "T3-T4")
  estimate_per_year <- t1$estimate / c(1,1,0.9,1,2.1) # calculate estimate per year (dividing by time interval)
  library(tibble)
  t1 <- add_column(t1, estimate_per_year, .after = "estimate")
  # save table
  write.csv(x=t1, file = paste0(out_dir_main, "/", PARENTING_save_name, "_untransformed.csv"), row.names = FALSE)
  
  ### check QQ Plots
  qqPlot(resid(imp_fit[[1]]))
  
  ########### MODEL WITH TRANSFORMED DV ----
  ### Turn each completed data set to long format:
  comp.out <- vector('list', length = m) #create empty list, that will be filled with outputs
  for(i in 1:m){
    comp.out[[i]] <- complete(m_out_obj, action = i)
  }
  comp.out_l <- lapply(comp.out, FUN = function(x){
    d_l <-  wideToLong(x, ind = 'P*_', nRep=4, repColnm = 'duration')
    ### Change reference:
    d_l$duration <- factor(d_l$duration, levels=c("P1_", "P2_", "P3_", "P4_"))
    return(d_l)
  })
  
  ### Fit mixed model (optim needed because of convergence issues) fit model with transformation of DV
  imp_fit.2 <- lapply(comp.out_l, function(x){
    lme(fixed = as.formula(paste(log,"(", PARENTING_BEH, pot,") ~ ADR_GenderChild + duration")), # If log = default of "" --> nothing pasted, if pot = default "" nothing added, leave one at default always when calling the function
        data = x, random = ~1 | ID, control = lmeControl(opt = 'optim'))
  })
  ### get ICC
  
  # Extract the first imputed model
  model <- imp_fit.2[[1]]
  
  # Extract the random effects variance
  # Extract the random effects variance
  random_var <- as.numeric(unlist(VarCorr(model))[1, 1])
  
  # Extract the residual variance
  residual_var <- as.numeric(unlist(VarCorr(model))[2, 1])
  
  # Calculate the ICC
  icc <- round(random_var / (random_var + residual_var), 2)
  
  # assign number of observations
  observations <- summary(model)[["dims"]][["N"]]
  
  N_of_subjects <- summary(model)[["dims"]][["ngrps"]][["ID"]]
  
  
  ### check QQ Plots if ^2 fits better
  qqPlot(resid(imp_fit.2[[1]]))
  
  ### Pool results:
  summary(pool(imp_fit.2)) ##
  k1 <- summary(pool(imp_fit.2))
  
  T1_4 <- k1[5,]
  
  k1 <- k1[1:3,] # keep relevant effects
  
  ## repeat above analyses with changing factor levels (T2_ =reference)
  ### Collect completed datasets in list:
  comp.out <- vector('list', length = m) #create empty list, that will be filled with outputs
  for(i in 1:m){
    comp.out[[i]] <- complete(m_out_obj, action = i)
  }
  comp.out_l <- lapply(comp.out, FUN = function(x){
    d_l <-  wideToLong(x, ind = 'P*_', nRep=4, repColnm = 'duration')
    ### Change reference:
    d_l$duration <- factor(d_l$duration, levels=c("P2_", "P3_", "P4_", "P1_"))
    return(d_l)
  })

  ### Fit mixed model (optim needed because of convergence issues)
  imp_fit.2 <- lapply(comp.out_l, function(x){
    lme(fixed = as.formula(paste(log,"(", PARENTING_BEH, pot,") ~ ADR_GenderChild + duration")),
        data = x, random = ~1 | ID, control = lmeControl(opt = 'optim'))
  })
  
  ## pool and add relevant effect to output table
  k2 <- summary(pool(imp_fit.2))
  k1 <- rbind(k1,k2[3,])
  ## repeat above analyses with changing factor levels (P3_ =reference)
  comp.out <- vector('list', length = m) #create empty list, that will be filled with outputs
  for(i in 1:m){
    comp.out[[i]] <- complete(m_out_obj, action = i)
  }
  comp.out_l <- lapply(comp.out, FUN = function(x){
    d_l <-  wideToLong(x, ind = 'P*_', nRep=4, repColnm = 'duration')
    ### Change reference:
    d_l$duration <- factor(d_l$duration, levels=c("P3_", "P4_", "P1_","P2_"))
    return(d_l)
  })

  ### Fit mixed model (optim needed because of convergence issues)
  imp_fit.2 <- lapply(comp.out_l, function(x){
    lme(fixed = as.formula(paste(log,"(", PARENTING_BEH, pot,") ~ ADR_GenderChild + duration")),
        data = x, random = ~1 | ID, control = lmeControl(opt = 'optim'))
  })
  ## pool and add relevant effect to output table
  k2 <- summary(pool(imp_fit.2))
  k1 <- rbind(k1,k2[3,])
  
  k1$term <- c("T1", "Sex", "T1-T2", "T2-T3", "T3-T4")
  
  estimate_per_year <- k1$estimate / c(1,1,0.9,1,2.1)
  
  library(tibble)
  k1 <- add_column(k1, estimate_per_year, .after = "estimate")
  k1
  
  
  T1_4 <- rbind(k1[1:2,-3],T1_4)
  T1_4 <- rbind(T1_4,k1[-c(1:2),-3])
  
  T1_4$term <- c("T1", "Sex","T1-T4", "T1-T2", "T2-T3", "T3-T4")
  estimate_per_year_T1_4 <- T1_4$estimate / c(1,1,4,0.9,1,2.1) # calculate estimate per year (dividing by time interval)
  
  T1_4 <- add_column(T1_4, estimate_per_year_T1_4, .after = "estimate")
  T1_4 <- T1_4%>%
    dplyr::rename(estimate_per_year = estimate_per_year_T1_4)
  
  T1_4 <- add_column(T1_4, icc)
  T1_4 <- add_column(T1_4, paste(N_of_subjects, "/", observations))
  T1_4 <- dplyr::rename(T1_4, "N of subjects / observations" = 9)
  T1_4
  # save table
  write.csv(x=k1, file = paste0(out_dir_seed, "/", PARENTING_save_name, "_pooled_30.csv"), row.names = FALSE)
  write.csv(x=T1_4, file = paste0(out_dir_seed, "/", PARENTING_save_name, "_pooled_30_T1_T4.csv"), row.names = FALSE)
  
  ######### format table ----
  # read in table
  table <- read.csv(paste0(out_dir_seed, "/",  PARENTING_save_name, "_pooled_30.csv"))
  table
  ## round
  table[,c(2:6)] <- round(table[,c(2:6)],2)
  table[,7] <- round(table[,7],3)
  table
  for (i in 1:5) {
    if (as.numeric(table$p.value[i]) < 0.001){
      table$p.value[i] <- "<.001"
    }
  }
  colnames(table) <- c("term", "estimate", "estimate per year", "SE", "statistic", "df", "p")
  table
  write.csv(x=table, file = paste0(out_dir_formatted,  "/",PARENTING_save_name, "_formatted.csv"), row.names = FALSE)
  
  # same with the T1-T4 effect included
  table <- read.csv(paste0(out_dir_seed, "/",  PARENTING_save_name, "_pooled_30_T1_T4.csv"))
  table
  ## round
  table[,c(2:6)] <- round(table[,c(2:6)],2)
  table[,7] <- round(table[,7],3)
  table
  for (i in 1:6) {
    if (as.numeric(table$p.value[i]) < 0.001){
      table$p.value[i] <- "<.001"
    }
  }
  colnames(table) <- c("Variable", "B", "B per year", "SE B", "t", "df", "p", "ICC", "N of subjects / observations")
  table
  write.csv(x=table, file = paste0(out_dir_formatted, "/",  PARENTING_save_name, "_formatted_incl_T1_T4.csv"), row.names = FALSE)
}

# FUNCTION CALLS

# info about arguments "log" and pot: If log not specified, the default of "" is taken--> nothing; if pot = default "" nothing added, leave one at default always when calling the function
# ev. change path with out_dir_main, out_dir_seed, out_dir_formatted

#Involvement
compute_mixed_models_T1_T4(m_out_obj = m.out, PARENTING_BEH = "Involvment_mean", "involvement_^2", pot="^2")

#Positive Parenting
compute_mixed_models_T1_T4(m_out_obj = m.out, PARENTING_BEH = "Positive_Parenting_mean","positive_parenting_^3", pot="^3")

#Poor Monitoring
compute_mixed_models_T1_T4(m_out_obj = m.out, PARENTING_BEH = "Poor_Monitoring_mean","poor_monitoring_log", log="log")

#Inconsistent Discipline
compute_mixed_models_T1_T4(m_out_obj = m.out, PARENTING_BEH = "Inconsistent_Discipline_mean","inconsistent_discipline")

#Corporal Punishment
compute_mixed_models_T1_T4(m_out_obj = m.out, PARENTING_BEH = "Corporal_Punishment_mean","corporal_punishment_log", log="log")
