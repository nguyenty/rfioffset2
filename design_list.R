list_model <- function(vnf, vntest){ # vnf is a vector of  all variable names
  n <- length(vntest)               # vntest is the vector of all variable names tested for significance   
  design.list <- vector("list", n+1)
  test.mat <- matrix(c(rep(1, n), 2:(n+1)), byrow = F, ncol = 2)
  row.names(test.mat) <- vntest
  if ((length(vntest) == 0) | (length(vnf) == 0)) stop ("vnf and vntest has to be nonempty")
  if (any(!vntest %in%vnf)) stop ("tested variables vntest have to be a subset of full set of variables vnf")
  fmf  <- as.formula(paste(" ~ ", paste(vnf, collapse= "+"))) 
  design.list[[1]] <- model.matrix(fmf)
  if (length(vnf) == 1){
    design.list[[2]] <- rep(1, length(get(vnf[1])) )
  } else {
    for(i in 1:n){
      fmr <- as.formula(paste(" ~ ", paste(setdiff(vnf, vntest[i]), collapse= "+"))) 
      design.list[[i+1]] <- model.matrix(fmr)
    }
  }
  return(list(design.list = design.list, test.mat = test.mat))
}
