fit_model_sim1 <- function(vnf, vntest, criteria, counts, model_th, subdir){ # model_th <- 1
  y <- counts
  log.offset <- log(apply(y, 2, quantile, 0.75))
  list_out <- list_model(vnf, vntest)
  design.list <- list_out$design.list
  test.mat <- list_out$test.mat
  fit <- QL.fit(y, design.list, test.mat, # dim(counts)
                 log.offset = log.offset, print.progress=TRUE,
                 Model = "NegBin", method = "optim")
  result<- QL.results(fit, Plot = FALSE)
  res_sel <- sel_criteria(result)
  k <- nrow(test.mat)
  name_model <- NULL 
  for (i in 1:k) name_model <- paste(name_model, row.names(test.mat)[i], sep =".")
  model_dir <- paste(resultdir,"/", subdir, "/",
                     colnames(res_sel)[criteria], 
                     "/Model",model_th,name_model, 
                     sep ="")
  dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)
  save(result, file = paste(model_dir,"/Model",model_th, "_result.RData", sep =""))
  save(fit, file = paste(model_dir,"/Model",model_th, "_fit.RData", sep =""))
  for(i in 1:(nrow(test.mat))){
    postscript(paste(model_dir,"/Model", 
                     model_th, row.names(test.mat)[i],".eps", sep =""))
    hist(result$P.values[[3]][,i], 
         main=row.names(test.mat)[i],
         xlab = "p-values", col = 'green',nclass=100)
    box()
    dev.off()
    
    pdf(paste(model_dir,"/Model", 
              model_th, row.names(test.mat)[i],".pdf", sep =""))
    hist(result$P.values[[3]][,i], 
         main=row.names(test.mat)[i],
         xlab = "p-values", col = 'green',nclass=100)
    box()
    dev.off()
  }
  print(paste("Model", model_th, sep = " "))
  
  return(res_sel)
}