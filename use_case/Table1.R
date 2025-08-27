library(Rcpp)
library(dplyr)
library(purrr)

# generate white noise
library(mvtnorm)

# to use start_with
library(tidyselect)

# to use dummy
library(fastDummies)

# to generate multivariate binary random variates
library(bindata)

# to calculate AUC
library(ROCR)

source("Fed_simulation_functions.R")

##################################################
###### Data pre-processing #######################
##################################################
## read in data
sample_data = read.csv("sample_data_to_run.csv")[,-1]

## site
sites = sample_data$new_site

## covariates matrix
design_matrix = sample_data[,-c(1,2,ncol(sample_data))]
new_ind <- which(apply(design_matrix,2,mean) > 0.1)
X = design_matrix[,new_ind] %>% as.matrix()

## outcome vector
treat = sample_data[,"outcome"] %>% as.matrix()

###################################################################################################
###### State-of-the-art to be compared with: when data can be pooled together #######################
###################################################################################################
####BHq###
fit_BHq_0<-BHq(Treatall=treat,
               Xall=X)
fit_BHq<-BH(fit_BHq_0$pvalue,0.1)



############################################################################
###### proposed method: when data can not be shared #######################
############################################################################
fit_proposed = fit_local_wide_lasso = fit_meta_BH = fit_SCAD = NULL
tryCatch({
  
  #### 1.0 --------- Our method:
  fit_proposed_raw<-local.glm.lasso(Treatall=treat,
                                    Xall=X,
                                    site=sites)
  fit_proposed<-MFDR(fit_proposed_raw$Theta,0.1)
  
  
  
  ##### other methods
  fit_proposed_0<-local.glm(Treatall=treat,
                            Xall=X,
                            site=sites)
  
  ##### 1.1 -------- local site use node wide : Fed_FDR_N
  fit_local_wide_lasso<-MFDR(fit_proposed_0$delasso,0.1)
  
  #### 1.2 ---------- metha+BH method
  fit_meta_BH<-BH(fit_proposed_0$pvalue,0.1)
  
  
  #### 1.3 ---------- Fed-FDR-SCAD method
  fit_SCAD<-local.glm.scad(Treatall=treat,
                           Xall=X,
                           site=sites)
  
  fit_SCAD_final<-MFDR(fit_SCAD$Theta,0.1)
  
  
},error=function(e){
  cat("ERROR :",conditionMessage(e), "\n")
})



############################################################################
###### Final comparison #######################
############################################################################
index_null = which(map_lgl(list(fit_BHq$FDR,
                                fit_local_wide_lasso$FDR,
                                fit_meta_BH$FDR,
                                fit_SCAD_final$FDR,
                                fit_proposed$FDR),is.null) == TRUE)
if (length(index_null) == 0){
  results_list = list(fit_BHq$FDR,
                      fit_local_wide_lasso$FDR,
                      fit_meta_BH$FDR,
                      fit_SCAD_final$FDR,
                      fit_proposed$FDR)
  results = matrix(unlist(results_list), 
                   ncol = length(fit_BHq$FDR), byrow = TRUE)
  rownames(results) = c("Gold_Standard_BHq",
                        "to_compared-local_wide_lasso",
                        "to_compared-meta_BH",
                        "to_compared-SCAD",
                        "proposed")
}else{
  results_list = list(fit_BHq$FDR,
                      fit_local_wide_lasso$FDR,
                      fit_meta_BH$FDR,
                      fit_SCAD_final$FDR,
                      fit_proposed$FDR)
  results = matrix(unlist(results_list[-index_null]), 
                   ncol = length(fit_BHq$FDR), byrow = TRUE)
  rownames(results) = c("Gold_Standard_BHq",
                        "to_compared-local_wide_lasso",
                        "to_compared-meta_BH",
                        "to_compared-SCAD",
                        "proposed")[-index_null]
}

colnames(results) = colnames(X)


output = subset(as.data.frame(results), 
                select = as.numeric(which(results[1,] == 1 & (apply(results,2,sum) != 1))))
if(ncol(output)  >= 2){
  output_values = apply(output, 1, sum)
}else{
  output_values = 0
}



if(length(output_values) == 5){
  output_new = subset(as.data.frame(results), 
                      select = as.numeric(which(apply(results,2,sum) >= 1)))
  print(paste("Output value is", output_values))
  print("------------------")
  print(output_new)
  print("------------------")
}



##### Save results
print("start_saving_results ---- ")
final_results = list(gold_standard = fit_BHq,
                     proposed = fit_proposed,
                     lasso = fit_local_wide_lasso,
                     SCAD = fit_SCAD_final,
                     meta = fit_meta_BH,
                     output = output,
                     output_values = output_values)

saveRDS(final_results, file = "output_table1.rds")

