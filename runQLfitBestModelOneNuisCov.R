source("readdata.R")
source("design_list.R")
source("fit_model.R")
source("simulation_criteria.R")
vnf <- c("Line", "Block", "Concb", 
          "RINa", "neut", "lymp", 
          "mono", "eosi", "baso")
vntest <- vnf

criteria <- 1
model_th <- 1
subdir <- "BestModelOneNuisCov"
repeat{
out_model <- fit_model(vnf, vntest, criteria, counts, model_th, subdir)
#out_model <- s[[i]]
assign(paste("ms_criteria", model_th, sep = "_" ),out_model)
ms_val <- get(paste("ms_criteria", model_th, sep = "_" ))
cov_del <- ms_val[1,criteria] # cov_del <- 14; i <- 1# model_th <- 2 # criteria <- 1

cov_set <- list_model(vnf, vntest)$test.mat # dim(cov_set)
res <- data.frame(criteria = colnames(ms_val)[criteria], 
                  model = model_th, 
                  cov_del = rownames(cov_set)[ cov_del])
if (cov_del ==1) break
vnf <- vntest <- vnf[-cov_del]
  model_th <- model_th +1
}

# write.csv(out_pairedend_cbc, "out_pairedend_cbc.csv", row.names = F)
proc.time() -pm1