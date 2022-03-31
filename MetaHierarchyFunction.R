MetaHierarchy <- function(dataf, study_id, minstudies, bootstrap_n, neg){
  if (!require(diagmeta)) install.packages('diagmeta')
  library(diagmeta)
  if (!require(mada)) install.packages('mada')
  library(mada) 
  if (!require(dmetatools)) install.packages('dmetatools')
  library(dmetatools)  
  if (!require(dplyr)) install.packages('dplyr')
  library(dplyr)  
  if (!require(lme4)) install.packages('lme4')
  library(lme4) 
  if (!require(tryCatchLog)) install.packages('tryCatchLog')
  library(tryCatchLog) 
  df_use <- dataf
  df_use <- df_use %>% mutate(threshold_use = case_when(
    neg==1 ~ threshold*-1,
    TRUE ~ threshold))
  k <- 2
  ind <- 8
  list_use <- vector(mode="list", k)
  for(i in seq(k)){
    list_use[[i]] <- NaN*seq(ind)
    for (j in seq(ind)){
      model_type = case_when(j==1 ~ "DIDS",
                             j==2 ~ "CIDS",
                             j==3 ~ "DICS",
                             j==4 ~ "CICS",
                             j==5 ~ "DS",
                             j==6 ~ "CS",
                             j==7 ~ "DI",
                             j==8 ~ "CI")
      model_result <- try(diagmeta(tp, fp, tn, fn, threshold_use, data=df_use, studlab = study_id, model=model_type, incr=0.5), silent=TRUE)
      
      list_use[[1]][[j]] <- tryCatch(model_result[[10]], error = function(e) paste("NA"))
      
      list_use[[2]][[j]] <- tryCatch(getME(model_result$result.lmer, name=c("devcomp"))[[1]][[7]], error = function(e) paste("NA"))
    }
  }
  model_resultsdf <- as.data.frame(list_use,col.names = c("modelfit", "criterion")) 
  
  if (sum(model_resultsdf$modelfit == "NA") < 8){
    model_resultsdf <- model_resultsdf %>% filter(modelfit != "NA") %>% slice_min(criterion)
    model_result <- diagmeta(tp, fp, tn, fn, threshold_use, data=df_use, studlab = study_id, model=model_resultsdf$modelfit, incr=0.5)
    model_result
  } else {
    limiteddf <- df_use %>% group_by(threshold_use) %>% filter(n() >= minstudies)
    for(i in seq(k)){
      list_use[[i]] <- NaN*seq(ind)
      for (j in seq(ind)){
        model_type = case_when(j==1 ~ "DIDS",
                               j==2 ~ "CIDS",
                               j==3 ~ "DICS",
                               j==4 ~ "CICS",
                               j==5 ~ "DS",
                               j==6 ~ "CS",
                               j==7 ~ "DI",
                               j==8 ~ "CI")
        model_result <- try(diagmeta(tp, fp, tn, fn, threshold_use, data=limiteddf, studlab = study_id, model=model_type, incr=0.5), silent=TRUE)
        
        list_use[[1]][[j]] <- tryCatch(model_result[[10]], error = function(e) paste("NA"))
        
        list_use[[2]][[j]] <- tryCatch(getME(model_result$result.lmer, name=c("devcomp"))[[1]][[7]], error = function(e) paste("NA"))
      }
    }
    model_resultsdf <- as.data.frame(list_use,col.names = c("modelfit", "criterion")) 
    if (sum(model_resultsdf$modelfit == "NA") < 8){
      model_resultsdf <- model_resultsdf %>% filter(modelfit != "NA") %>%  slice_min(as.numeric(criterion))
      model_result <- diagmeta(tp, fp, tn, fn, threshold_use, data=limiteddf, studlab = study_id, model=model_resultsdf$modelfit, incr=0.5)
      model_result
      .GlobalEnv$DFused <- limiteddf 
    } else {
      limiteddf <- df_use %>% group_by(threshold) %>% filter(n() >= minstudies)
      split_data <- split(limiteddf, f = limiteddf$threshold)
      list2env(split_data,envir = .GlobalEnv) 
      my_list <- list() 
      for(i in 1:length(split_data))
      {   my_list[[i]] <- reitsma(data=split_data[[i]], TP="tp", FN="fn", FP="fp", TN="tn")
      print("NEW MODEL: threshold listed below")
      print(split_data[[i]]$threshold[[1]]) 
      print(summary(my_list[[i]]))
      print(AUC_boot(split_data[[i]]$tp,split_data[[i]]$fp,split_data[[i]]$fn,split_data[[i]]$tn,B=bootstrap_n))
      }}
  }
}