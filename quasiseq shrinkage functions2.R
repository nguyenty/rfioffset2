myQL.fit <- function(counts, design.list, betavec=NULL, test.mat = NULL, log.offset = NULL, Model = "NegBin", method="glm", print.progress = TRUE, 
                     NBdisp = "trend", ...) {
  if(length(grep(method,c("optim","glm")))!=1) stop(paste("Supplied method argument =",method,"does not uniquely match either 'optim' or 'glm'."))
  if(is.data.frame(counts)) counts<-as.matrix(counts)
  ### Note: First element of design.list should pertain to overall full model.  This is the design used to obtain
  ### dispersion estimates for quasi-likelihood models.
  if(any(rowSums(counts)==0)) stop(cat(sum(rowSums(counts)==0)," genes have 0 counts across all samples. \n Please remove genes with zero total counts before analyzing.\n",sep=""))
  ### Check for errors
  if (any(round(counts) != counts)) 
    stop("Count data contains non-integers.")
  if (any(counts < 0)) 
    stop("Count data contains negative counts.")
  
  if (!Model %in% c("NegBin", "Poisson")) 
    stop(cat("Unidentified Model: Model must be either 'NegBin' or 'Poisson'.\n"))
  if (Model == "NegBin") {
    
    if (!NBdisp[1] %in% c("trend", "common") & length(NBdisp) != nrow(counts)) 
      stop(cat("Unidentified NegBin Dispersion: NBdisp must be set as 'trend' or 'common' to estimate negative binomial dispersion from data using GLM edgeR (McCarthy et al., 2012),\n or it must be a vector providing negative binomial dispersion parameter value to use for each gene.\n"))
    
    if (length(NBdisp) == nrow(counts) & !is.numeric(NBdisp)) 
      stop(cat("NBdisp contains non-numeric values.\n\tAll negative binomial dispersion parameters must be non-negative real numbers.\n"))
    
    if (length(NBdisp) == nrow(counts) & any(NBdisp < 0)) 
      stop(cat("NBdisp contains negative values.\nAll negative binomial dispersion parameters must be non-negative real numbers.\n"))
  }
  
  ### Fit model and evaluate deviance under each design provided in design.list
  deviance.list <- vector("list", length(design.list))
  p <- NULL
  n <- ncol(counts)  # p is used to store the d.f. for each model (it will be a vector) and n is the total number of samples (also the number of observations for each gene)
  
  for (jj in 1:length(design.list)) {
    design <- design.list[[jj]]
    
    if (is.vector(design)) {
      p <- c(p, length(unique(design)))  ## Record the d.f. for current model
      
      ### Check for errors if the current model design is specified as a vector
      if (p[jj] > p[1]) 
        stop(cat("Full model design must be first element in 'design.list'.\n'p' for element", jj, "is larger than 'p' for first element,\nindicating first element does not provide full model design.\n"))
      if (length(design) != n) 
        stop(cat("Element", jj, "in 'design.list' has length", length(design), ".\nDesign vectors must have length", 
                 n, "(to match number of columns in data).\n"))
    }
    
    if (is.matrix(design)) {
      p <- c(p, ncol(design))
      
      ### Check for errors if the current model design is specified as a matrix 
      ### if(prod(design[,1]==1)!=1)
      ### stop(paste('The first column of matrix in element',jj,'of 'design.list' is not a column of 1s for the
      ### intercept. Please include intercept.'))
      if (nrow(design) != n) 
        stop(cat("Element", jj, "in 'design.list' has", nrow(design), "rows.\nDesign matrices must have", 
                 n, "rows (to match number of columns in data).\n"))
      if (p[jj] > p[1]) 
        stop(cat("Full model design must be first element in 'design.list'.\n'p' for element", jj, "is larger than 'p' for first element,\nindicating first element does not provide full model design.\n"))
    }
    
    ### Analyze using quasi-negative binomial model, if chosen
    if (Model == "NegBin") {
      if (is.vector(design)) {
        if (length(unique(design)) > 1) 
          design <- model.matrix(~as.factor(design))
        if (length(unique(design)) == 1) 
          design <- matrix(1, ncol(counts), 1)
      }
      if (jj == 1) {
        if (is.null(log.offset)) 
          d <- DGEList(counts = counts, group = design[, 2], lib.size = rep(1, ncol(counts)))
        if (!is.null(log.offset)) 
          d <- DGEList(counts = counts, group = design[, 2], lib.size = exp(log.offset))
        #                d <- calcNormFactors(d)
        ### If requested, use gene-specific trended dispersion estimates from GLM edgeR (McCarthy et al., 2012).
        ### if(NBdisp=='trend')nb.disp<-estimateGLMTrendedDisp(d, design,...)$trended.dispersion
        if (NBdisp == "trend") 
          nb.disp <- myestimateGLMTrendedDispDGEList(d, design, betavec=betavec)$trended.dispersion
        
        ### If requested, use common dispersion estimate from GLM edgeR (McCarthy et al., 2012).
        if (NBdisp == "common") 
          nb.disp <- rep(estimateGLMCommonDisp(d, design, ...)$common.dispersion, nrow(counts))
        
        ### If provided, use prespecified dispersion estimates.
        if (length(NBdisp) == nrow(counts)) {
          if (is.numeric(NBdisp) & !any(NBdisp < 0)) 
            nb.disp <- NBdisp
        }
      }
      
      ### Analyze genes with positive dispersion parameters using quasi-negative binomial model
      if (any(nb.disp > 0)) 
        res <- myNBDev(counts[nb.disp > 0, ], design, betavec, log.offset, nb.disp[nb.disp > 0],method, print.progress)
      
      ### If present, analyze genes for which nb.disp==0 using quasi-Poisson model
      if (any(nb.disp == 0)) {
        res2 <- myPoisDev(counts[nb.disp == 0, ], design, betavec, log.offset, print.progress)
        means <- dev <- rep(NA, nrow(counts))
        parms <- matrix(NA, nrow(counts), ncol(design))
        means[nb.disp == 0] <- res2$means
        dev[nb.disp == 0] <- res2$dev
        parms[nb.disp == 0, ] <- res2$parms
        if (any(nb.disp > 0)) {
          means[nb.disp > 0] <- res$means
          dev[nb.disp > 0] <- res$dev
          parms[nb.disp > 0, ] <- res$parms
        }
        res <- list(dev = dev, means = means, parms = parms)
      }
    }
    
    ### Analyze using quasi-Poisson model, if chosen
    if (Model == "Poisson") res <- PoisDev(counts, design, log.offset, print.progress)
    ### Record means and parameter estimate from full model
    if (jj == 1) {
      means <- res$means
      parms <- res$parms
    }
    deviance.list[[jj]] <- res$dev
  }
  
  LRT <- num.df <- NULL
  
  if (length(design.list) > 1) {
    ### Compute likelihood ratio test statistics. If not otherwise specified, compare each model to the first model
    ### in design.list, which should be the full model
    
    if (is.null(test.mat)) {
      cat("Note: 'test.mat' not provided. Comparing each model \nfrom 'design.list' to first model in 'design.list', which must be the full model\n")
      test.mat <- cbind(1, 2:length(design.list))
      rownames(test.mat) <- paste("Design", 1, " vs Design", 2:length(design.list), sep = "")
    }
    
    for (i in 1:nrow(test.mat)) {
      i1 <- test.mat[i, 1]
      i2 <- test.mat[i, 2]
      num.df <- c(num.df, abs(p[i2] - p[i1]))
      LRT <- cbind(LRT, -(deviance.list[[i2]] - deviance.list[[i1]])/(p[i2] - p[i1]))
    }
    colnames(LRT) <- rownames(test.mat)
  }
  den.df <- (n - p[1])
  
  ### Compute deviance dispersion estimate
  phi.hat.dev <- deviance.list[[1]]/den.df
  
  ### Compute Pearson dispersion estimate
  if (Model == "NegBin") 
    phi.hat.pearson <- (means - counts)^2/(means + means^2 * nb.disp)
  if (Model == "Poisson") 
    phi.hat.pearson <- (means - counts)^2/means
  phi.hat.pearson[means == 0] <- 0
  phi.hat.pearson <- rowSums(phi.hat.pearson)/den.df
  
  if (Model == "Poisson") 
    return(list(LRT = LRT, phi.hat.dev = phi.hat.dev, phi.hat.pearson = phi.hat.pearson, mn.cnt = rowMeans(counts), 
                den.df = den.df, num.df = num.df, Model = Model, fitted.values = means, coefficients = parms))
  if (Model == "NegBin") 
    return(list(LRT = LRT, phi.hat.dev = phi.hat.dev, phi.hat.pearson = phi.hat.pearson, mn.cnt = rowMeans(counts), 
                den.df = den.df, num.df = num.df, Model = Model, NB.disp = nb.disp, fitted.values = means, coefficients = parms))
}

#Gets NegBin Dispersion; adjusted for 
mydispPearson <- function (y, design = NULL,betavec=NULL,offset = NULL, min.row.sum = 5, subset = 10000, 
    AveLogCPM = NULL, tol = 1e-06, trace = FALSE, initial.dispersion = 0.1) 
{
    y <- as.matrix(y)
    if (is.null(design)) {
        design <- matrix(1, ncol(y), 1)
        rownames(design) <- colnames(y)
        colnames(design) <- "Intercept"
    }
    else {
        design <- as.matrix(design)
    }
    if (is.null(offset)) 
        offset <- rep(0,ncol(y))
    small.row.sum <- which(rowSums(y) < min.row.sum)
    if (length(small.row.sum)) {
        y <- y[-small.row.sum, , drop = FALSE]
        offset <- offset[-small.row.sum, , drop = FALSE]
        betavec <- betavec[-small.row.sum,]
    }
    if (nrow(y) < 1) 
        stop("no data rows with required number of counts")
    if (!is.null(subset) && subset <= nrow(y)/2) {
        if (is.null(AveLogCPM)) 
            AveLogCPM <- aveLogCPM(y, offset = offset)
        else {
            if (length(small.row.sum)) 
                AveLogCPM <- AveLogCPM[-small.row.sum]
        }
        i <- systematicSubset(subset, AveLogCPM)
        y <- y[i, , drop = FALSE]
        offset <- offset[i, , drop = FALSE]
        betavec <- betavec[i,]
    }
    mu <- matrix(NA,nrow(y),ncol(y))
    for(i in 1:nrow(y)){
    	ADJoffset <- offset+betavec[i,]
   		fit <- glmFit(y = y[i,,drop=FALSE], design = design, dispersion = initial.dispersion, 
        	offset = ADJoffset, prior.count = 0)
		mu[i,] <- fit$fitted.values
    	}
    nlibs <- ncol(y)
    df.residual <- nlibs - ncol(design)
    one <- df.residual/nlibs
    phi <- 0
    iter <- 0
    pos <- mu > 0
    y <- y[pos]
    mu <- mu[pos]
    repeat {
        iter <- iter + 1
        s2 <- (y - mu)^2
        Q <- mean(s2/mu/(1 + phi * mu))
        dQ <- mean(s2/(1 + phi * mu)^2)
        dif <- (Q - one)/dQ
        if (dif < 0) 
            break
        phi <- phi + dif
        if (trace) 
            cat(iter, phi, Q, dQ, dif, "\n")
        if (dif < tol) 
            break
        if (iter > 100) {
            warning("iteration limit reached")
            break
        }
    }
    phi
}

myestimateGLMCommonDisp <- function (y, design = NULL, betavec=NULL, offset = NULL, method = "Pearson", 
    subset = 10000, AveLogCPM = NULL, verbose = FALSE, ...) 
{
    y <- as.matrix(y)
    if (is.null(design)) {
        design <- matrix(1, ncol(y), 1)
        rownames(design) <- colnames(y)
        colnames(design) <- "Intercept"
    }
    else {
        design <- as.matrix(design)
    }
    if (ncol(design) >= ncol(y)) {
        warning("No residual df: setting dispersion to NA")
        return(NA)
    }
    method <- match.arg(method, c("CoxReid", "Pearson", "Pearson2", 
        "deviance"))
    if (is.null(offset)) 
        offset <- log(colSums(y))
    if (is.null(AveLogCPM)) 
        AveLogCPM <- aveLogCPM(y)
    disp <- switch(method, CoxReid = dispCoxReid(y, design = design, 
        offset = offset, subset = subset, AveLogCPM = AveLogCPM, 
        ...), Pearson = mydispPearson(y, design = design, offset = offset, betavec=betavec, 
        subset = subset, AveLogCPM = AveLogCPM, ...), deviance = dispDeviance(y, 
        design = design, offset = offset, subset = subset, AveLogCPM = AveLogCPM, 
        ...))
    if (verbose) 
        cat("Disp =", round(disp, 5), ", BCV =", round(sqrt(disp), 
            4), "\n")
    disp
}

mydispBinTrend <- function (y, design = NULL, betavec=NULL, offset = NULL, df = 5, span = 0.3, 
    min.n = 400, method.bin = "Pearson", method.trend = "spline", 
    AveLogCPM = NULL, ...) 
{
    y <- as.matrix(y)
    nlibs <- ncol(y)
    ntags <- nrow(y)
    pos <- rowSums(y) > 0
    if (!any(pos)) 
        return(AveLogCPM = AveLogCPM, dispersion = rep(0, ntags))
    npostags <- sum(pos)
    if (is.null(design)) {
        design <- matrix(1, nlibs, 1)
    }
    else {
        design <- as.matrix(design)
    }
    if (is.null(offset)) 
    	offset <- log(colSums(y))
    method.bin <- match.arg(method.bin, c("CoxReid", "Pearson", 
        "deviance"))
    method.trend <- match.arg(method.trend, c("spline", "loess"))
    if (is.null(AveLogCPM)) 
        AveLogCPM <- aveLogCPM(y)
    group <- as.numeric(pos)
    if (npostags < 100) 
        nbins <- 1
    else {
        nbins <- floor(npostags^0.4)
        nbins <- min(nbins, 1000)
        min.n <- min(min.n, floor(npostags/nbins))
    }
    if (min.n < 50) {
        nbins <- floor(npostags/50)
        min.n <- 50
    }
    if (nbins > 1) 
        group[pos] <- cutWithMinN(AveLogCPM[pos], intervals = nbins, 
            min.n = min.n)$group
    bin.d <- bin.A <- rep(0, nbins)
    for (i in 1:nbins) {
        bin <- group == i
        bin.d[i] <- myestimateGLMCommonDisp(y[bin, ], design, betavec=betavec, method = method.bin, 
            offset, min.row.sum = 0, ...)
        bin.A[i] <- mean(AveLogCPM[bin])
    }
    if (nbins == 1) {
        dispersion <- rep.int(bin.d, ntags)
        return(list(AveLogCPM = AveLogCPM, dispersion = dispersion, 
            bin.AveLogCPM = bin.A, bin.dispersion = bin.d))
    }
    if (nbins < 7) {
        f <- approxfun(bin.A, sqrt(bin.d), rule = 2)
        dispersion <- f(AveLogCPM)^2
        return(list(AveLogCPM = AveLogCPM, dispersion = dispersion, 
            bin.AveLogCPM = bin.A, bin.dispersion = bin.d))
    }
    if (method.trend == "spline") {
        require("splines")
        p1 <- (1:(df - 1))/df
        knots1 <- quantile(bin.A, p = p1)
        r <- range(bin.A)
        knots2 <- r[1] + p1 * (r[2] - r[1])
        knots <- 0.3 * knots1 + 0.7 * knots2
        basisbins <- ns(bin.A, df = df, knots = knots, intercept = TRUE)
        beta <- coefficients(lm.fit(basisbins, sqrt(bin.d)))
        basisall <- predict(basisbins, newx = AveLogCPM)
        dispersion <- drop(basisall %*% beta)^2
    }
    if (method.trend == "loess") {
        fit <- loessFit(sqrt(bin.d), bin.A, span = span, iterations = 1)
        f <- approxfun(bin.A, fit$fitted, rule = 2)
        dispersion <- f(AveLogCPM)^2
    }
    list(AveLogCPM = AveLogCPM, dispersion = dispersion, bin.AveLogCPM = bin.A, 
        bin.dispersion = bin.d)
}

myestimateGLMTrendedDisp <- function (y, design = NULL, betavec=NULL, offset = NULL, AveLogCPM = NULL, 
    method = "auto", ...) 
{
    y <- as.matrix(y)
    ntags <- nrow(y)
    if (is.null(design)) {
        design <- matrix(1, ncol(y), 1)
        rownames(design) <- colnames(y)
        colnames(design) <- "Intercept"
    }
    else {
        design <- as.matrix(design)
    }
    if (ncol(design) >= ncol(y)) {
        warning("No residual df: cannot estimate dispersion")
        return(NA, ntags)
    }
    if (is.null(offset)) {
        lib.size <- colSums(y)
        offset <- log(lib.size)
    }
    if (is.null(AveLogCPM)) 
        AveLogCPM <- aveLogCPM(y, lib.size = exp(expandAsMatrix(offset,dim(y))))
    method <- match.arg(method, c("auto", "bin.spline", "bin.loess", 
        "power", "spline"))
    if (method == "auto") {
        if (ntags < 200) {
            method <- "power"
        }
        else {
            method <- "bin.spline"
        }
    }
    trend <- switch(method, bin.spline = mydispBinTrend(y, design, betavec=betavec,
        offset = offset, method.trend = "spline", AveLogCPM = AveLogCPM, 
        ...), bin.loess = dispBinTrend(y, design, offset = offset, 
        method.trend = "loess", AveLogCPM = AveLogCPM, ...), 
        power = dispCoxReidPowerTrend(y, design, offset = offset, 
            AveLogCPM = AveLogCPM, ...), spline = dispCoxReidSplineTrend(y, 
            design, offset = offset, AveLogCPM = AveLogCPM, ...))
    trend$design <- design
    trend$AveLogCPM <- AveLogCPM
    trend$trend.method <- method
    trend
}

myestimateGLMTrendedDispDGEList <- function (y, design = NULL, betavec=NULL, offset = NULL, AveLogCPM = NULL, 
    method = "auto", ...) 
{
    if (!is.null(AveLogCPM)) 
        y$AveLogCPM <- AveLogCPM
    if (is.null(y$AveLogCPM)) 
        y$AveLogCPM <- aveLogCPM(y)
    d <- myestimateGLMTrendedDisp(y = y$counts, design = design, betavec=betavec,
        offset = getOffset(y), AveLogCPM = y$AveLogCPM, method = method, 
        ...)
    y$trended.dispersion <- d$dispersion
    y$trend.method <- d$trend.method
    y$bin.dispersion <- d$bin.dispersion
    y$bin.AveLogCPM <- d$bin.AveLogCPM
    y$design <- d$design
    y
}

myNBDev <- function(counts, design, betavec, log.offset, nb.disp, method, print.progress = TRUE) {

    n <- ncol(counts)
    
    if (is.null(log.offset)) 
        log.offset <- rep(0, ncol(counts))
    est.offset <- exp(log.offset)

  SAT.LIKE <- function(counts, disp) {
        means <- counts
        like <- disp * log(disp/(disp + means))
        like[counts != 0] <- like[counts != 0] + counts[counts != 
            0] * log(means[counts != 0]/(disp + means[counts != 
            0]))
        -sum(like)
    }
    LIKE <- function(parms, design, counts, disp, est.offset) {
        means <- as.vector(exp(design %*% parms) * est.offset)
        like <- disp * log(disp/(disp + means))
        like[counts != 0] <- like[counts != 0] + counts[counts != 
            0] * log(means[counts != 0]/(disp + means[counts != 
            0]))
        -sum(like)
    }
    GRAD <- function(parms, design, counts, disp, est.offset) {
        means <- as.vector(exp(design %*% parms) * est.offset)
        colSums(-(counts - means * (disp + counts)/(disp + means)) * 
            design)
    }

    
    deviance.vector <- rep(NA, nrow(counts))
    means <- matrix(NA, nrow(counts), ncol(counts))
    parms <- matrix(NA, nrow(counts), ncol(design))

    ### For each gene and given design matrix, fit glm using provided negative binomial dispersion estimate
    for (gn in 1:nrow(counts)) {
        ### If wanted, provide running progress update (eventually once every 5000 genes)
        ADJlog.offset <- betavec[gn,]+log.offset
        if (gn %in% c(2, 10, 100, 500, 1000, 2500, 5000 * (1:200)) & print.progress) 
            print(paste("Analyzing Gene #", gn))
        
if(length(grep(method,"optim"))==1){
        if (ncol(design) > 1) 
            init.parms <- lm(log(counts[gn, ] + 1) ~ design[, 
                -1], offset = ADJlog.offset)$coefficients
        if (ncol(design) == 1) 
            init.parms <- lm(log(counts[gn, ] + 1) ~ 1, offset = ADJlog.offset)$coefficients
        opt <- optim(init.parms, fn = LIKE, gr = GRAD, method = "BFGS", 
            design = design, control = list(reltol = 1e-25, maxit = 1000), 
            counts = counts[gn, ], disp = 1/nb.disp[gn], est.offset = exp(ADJlog.offset))
        means[gn, ] <- as.vector(exp(design %*% opt$par) * exp(ADJlog.offset))
        parms[gn, ] <- opt$par
        deviance.vector[gn] <- 2 * (opt$value - SAT.LIKE(counts[gn, 
            ], 1/nb.disp[gn]))
    }

if(length(grep(method,"glm"))==1){
        #### For 2000 Fly genes, glm takes roughly 9 seconds (optim took roughly 21 seconds)
        glm.fitted <- glm(formula = counts[gn, ] ~ . - 1 + offset(ADJlog.offset), family = negbin.br("log", nb.disp[gn]), 
            data = as.data.frame(design), control = glm.control(epsilon = 1e-08, maxit = 100, trace = FALSE))
        parms[gn, ] <- coefficients(glm.fitted)
        
        ### Save optimized means (used in Pearson's dispersion estimator)
        means[gn, ] <- as.vector(exp(design %*% coefficients(glm.fitted)) * exp(ADJlog.offset))
        
        ### Save deviance (used to compute LRT for comparing models and also deviance dispersion estimator)
        deviance.vector[gn] <- glm.fitted$deviance
        
        ### When a group has only zero counts, the MLE doesn't exist.  For these genes, we use the method Kosmidis & Firth(2009)
        ### to moderate the parameter estimates.
        ### For 2000 Fly genes, fbr takes roughly 41 seconds
        if (any(counts[gn, ] == 0) && any(abs(coefficients(glm.fitted)) > 3)) {
            fbrflm.fit <- fbrglm(formula = counts[gn, ] ~ . - 1 + offset(ADJlog.offset), data = as.data.frame(design), start=coefficients(glm.fitted),
                family = negbin.br("log", nb.disp[gn]), control = glm.control(epsilon = 1e-08, maxit = 100, trace = FALSE))
            parms[gn, ] <- coefficients(fbrflm.fit)
        }
    }
      }
    return(list(dev = deviance.vector, means = means, parms = parms))
} 

myPoisDev <- function (counts, design, betavec, log.offset, print.progress = TRUE) 
{
    n <- ncol(counts)
    if (is.vector(design)) {
            offset <- expandAsMatrix((log.offset),dim(counts))
            offset <- exp(offset+betavec)
        counts <- as.matrix(counts)
        means <- counts
        parms <- NULL
        for (i in 1:length(unique(design))) {
            if (sum(design == unique(design)[i]) > 1) 
                means[, design == unique(design)[i]] <- rowSums(counts[, 
                  design == unique(design)[i]])/rowSums(offset[,design == 
                  unique(design)[i]])
            if (sum(design == unique(design)[i]) == 1) 
                means[, design == unique(design)[i]] <- counts[, 
                  design == unique(design)[i]]/offset[,design == 
                  unique(design)[i]]
            parms <- cbind(parms, means[, design == unique(design)[i]][, 
                1])
        }
        parms <- cbind(log(parms[, 1]), log(parms[, -1]/parms[, 
            1]))
        colnames(parms) <- c("(Intercept)", unique(design)[-1])
        means <- means * offset
        deviance <- means - counts
        deviance[counts != 0] <- deviance[counts != 0] + (counts * 
            log(counts/means))[counts != 0]
        deviance <- 2 * rowSums(deviance)
    }
    if (is.matrix(design)) {
        deviance <- rep(NA, nrow(counts))
        means <- matrix(NA, nrow(counts), ncol(counts))
        parms <- matrix(NA, nrow(counts), ncol(design))
        for (i in 1:nrow(counts)) {
            if (i %in% c(2, 10, 100, 500, 1000, 2500, 4000, 5000 * 
                (1:200)) & print.progress) 
                print(paste("Analyzing Gene #", i))
            res <- glm(counts[i, ] ~ design[, -1], family = "poisson", 
                offset = betavec[i,]+log.offset)
            means[i, ] <- res$fitted.values
            parms[i, ] <- res$coefficients
            deviance[i] <- res$deviance
        }
    }
    return(list(dev = deviance, means = means, parms = parms))
}

negbin.br=function( link = "log", overdisp = stop("'overdisp' must be specified"))
{
	mu.eta=variance=linkfun=function(...){return('Pass CRAN checker.')}
	ans=mgcv::negbin(1/overdisp,link)
	assign('variance', ans$variance, envir=environment(ans$variance))
	assign('mu.eta', ans$mu.eta, envir=environment(ans$variance))
	assign('linkfun', ans$linkfun, envir=environment(ans$variance))
	# if(ans$link=='log'){
		# tmp=function(eta)pmin(pmax(exp(eta), .Machine$double.eps), .Machine$double.xmax)
		# environment(tmp)=environment(ans$linkinv)
		# ans$linkinv = tmp
	# }
	ans$adjresponse=function(hat, mu, weight,...)
	{
		.5*hat*(1+variance(mu)/mu.eta(linkfun(mu)))
	}
	environment(ans$adjresponse)=environment(ans$getTheta)
	class(ans)=c('fbrfamily','family')
	ans
}

fbrglm = function (formula, family = negbin.br('log',0), data, weights, subset, 
    na.action, start = NULL, etastart, mustart, offset, control = list(...), 
    model = TRUE, method = "glm.fit", x = FALSE, y = TRUE, contrasts = NULL, 
    ...) 
{
    call <- match.call()
    if (is.character(family)) 
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) 
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
	if (!inherits(family,"fbrfamily"))	stop("'family' needs to inherit 'brfamily'.")

    if (missing(data)) 
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action", 
        "etastart", "mustart", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    if (identical(method, "model.frame")) 
        return(mf)
    if (!is.character(method) && !is.function(method)) 
        stop("invalid 'method' argument")
    if (identical(method, 'glm.fit')) {
		method='glm.fit0'
		tmp.mit=control$maxit
       control <- do.call("glm.control", control)
		if(is.null(tmp.mit)) control$maxit = 1e3L
	}
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm)) 
            names(Y) <- nm
    }
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(Y), 0L)
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights)) 
        stop("'weights' must be a numeric vector")
    if (!is.null(weights) && any(weights < 0)) 
        stop("negative weights not allowed")
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) != NROW(Y)) 
            stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
                length(offset), NROW(Y)), domain = NA)
    }
	keep.x=x; keep.y=y
    mustart <- model.extract(mf, "mustart")
    etastart <- model.extract(mf, "etastart")
	control0=control
	control0$maxit=1L
	yy=Y; sts=start; ests=etastart; musts=mustart; iter=1L

	repeat{		
		fit <- suppressWarnings(
			eval(call(if (is.function(method)) "method" else method, 
			x = X, y = yy, weights = weights, start = sts, etastart = if(iter==1L) ests else NULL, 
			mustart = musts, offset = offset, family = family, 
			control = control0, intercept = attr(mt, "intercept") > 0L)))
		adjresponse = family$adjresponse(rowSums(qr.Q(fit$qr)^2), fit$fitted.values, fit$prior.weights)
		# adjresponse = family$adjresponse(
			# rowSums(qr.Q(qr(family$mu.eta(fit$linear)/sqrt(family$var(fit$fitted))*X))^2),
			# fit$fitted.values
		# )
		convb=abs(max(-Inf, abs(sts-fit$coef),na.rm=TRUE))/max(c(.1, abs(fit$coef)),na.rm=TRUE)
		if(control$trace>0)cat("iteration:", iter, "\t:convb =", convb, '\n')
		if( convb < control$epsilon ) {
			fit$converge = TRUE
			fit$iter = iter
			break
		}
		if(iter >= control$maxit) {
			warning('did not converge.')
			fit$converge = FALSE
			fit$iter = iter
			break
		}
		yy=Y+adjresponse
		sts = fit$coef;					#sts = NULL
		ests = fit$linear.predictors;  #ests = NULL
		musts = fit$fitted.values
		iter = iter + 1L
	}
	fit$vcov=solve(as.matrix(Matrix::nearPD(crossprod(qr.R(fit$qr)))$mat))  ## inverse Fisher's information, not yet observed info
	# fit$vcov=tcrossprod(backsolve(qr.R(fit$qr)[,order(fit$qr$pivot),drop=FALSE], diag(1,NCOL(X)))) 

    if (FALSE && length(offset) && attr(mt, "intercept") > 0L) {
        fit2 <- eval(call(if (is.function(method)) "method" else method, 
            x = X[, "(Intercept)", drop = FALSE], y = Y, weights = weights, 
            offset = offset, family = family, control = control, 
            intercept = TRUE))
        if (!fit2$converged) 
            warning("fitting to calculate the null deviance did not converge -- increase 'maxit'?")
        fit$null.deviance <- fit2$deviance
    }
	fit$null.deviance = NA_real_
    if (model) 
        fit$model <- mf
    fit$na.action <- attr(mf, "na.action")
    if (x) 
        fit$x <- X
    if (!y) {
        fit$y <- NULL
	}else fit$y = Y
    fit <- c(fit, list(call = call, formula = formula, terms = mt, 
        data = data, offset = offset, control = control, method = method, 
        contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, 
            mf)))
    class(fit) <- c(fit$class, c("fbrglm", "glm", "lm"))
    fit
}

glm.fit0 <- function (x, y, weights = rep(1, nobs), start = NULL, etastart = NULL, 
    mustart = NULL, offset = rep(0, nobs), family = gaussian(), 
    control = list(), intercept = TRUE) 
{
	n='Shut up CRAN checker.'
    control <- do.call("glm.control", control)
    x <- as.matrix(x)
    xnames <- dimnames(x)[[2L]]
    ynames <- if (is.matrix(y)) 
        rownames(y)
    else names(y)
    conv <- FALSE
    nobs <- NROW(y)
    nvars <- ncol(x)
    EMPTY <- nvars == 0
    if (is.null(weights)) 
        weights <- rep.int(1, nobs)
    if (is.null(offset)) 
        offset <- rep.int(0, nobs)
    variance <- family$variance
    linkinv <- family$linkinv
    if (!is.function(variance) || !is.function(linkinv)) 
        stop("'family' argument seems not to be a valid family object", 
            call. = FALSE)
    dev.resids <- family$dev.resids
    aic <- family$aic
    mu.eta <- family$mu.eta
    unless.null <- function(x, if.null) if (is.null(x)) 
        if.null
    else x
    valideta <- unless.null(family$valideta, function(eta) TRUE)
    validmu <- unless.null(family$validmu, function(mu) TRUE)
    if (is.null(mustart)) {
        eval(family$initialize)
    }
    else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }
    if (EMPTY) {
        eta <- rep.int(0, nobs) + offset
        if (!valideta(eta)) 
            stop("invalid linear predictor values in empty model", 
                call. = FALSE)
        mu <- linkinv(eta)
        if (!validmu(mu)) 
            stop("invalid fitted means in empty model", call. = FALSE)
        dev <- sum(dev.resids(y, mu, weights))
        w <- ((weights * mu.eta(eta)^2)/variance(mu))^0.5
        residuals <- (y - mu)/mu.eta(eta)
        good <- rep(TRUE, length(residuals))
        boundary <- conv <- TRUE
        coef <- numeric()
        iter <- 0L
    }
    else {
        coefold <- NULL
        eta <- if (!is.null(etastart)) 
            etastart
        # else if (!is.null(start)) 
            # if (length(start) != nvars) 
                # stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s", 
                  # nvars, paste(deparse(xnames), collapse = ", ")), 
                  # domain = NA)
            # else {
                # coefold <- start
                # offset + as.vector(if (NCOL(x) == 1L) 
                  # x * start
                # else x %*% start)
            # }
        else family$linkfun(mustart)
        mu <- linkinv(eta)
        if (!(validmu(mu) && valideta(eta))) 
            stop("cannot find valid starting values: please specify some", 
                call. = FALSE)
        devold <- sum(dev.resids(y, mu, weights))
        boundary <- conv <- FALSE
        for (iter in seq_len(control$maxit)) {
            good <- weights > 0
            varmu <- variance(mu)[good]
            if (any(is.na(varmu))) 
                stop("NAs in V(mu)")
            if (any(varmu == 0)) 
                stop("0s in V(mu)")
            mu.eta.val <- mu.eta(eta)
            if (any(is.na(mu.eta.val[good]))) 
                stop("NAs in d(mu)/d(eta)")
            good <- (weights > 0) & (mu.eta.val != 0)
            if (all(!good)) {
                conv <- FALSE
                warning(gettextf("no observations informative at iteration %d", 
                  iter), domain = NA)
                break
            }
            z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
            w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
            ngoodobs <- as.integer(nobs - sum(!good))
			if(FALSE){ ## operational way to reproduce fit
				tmp.zw=z*w; tmp.tol=min(1e-7, control$epsilon/1000)
				fit=qr(x[good,,drop=FALSE] *w, tol=tmp.tol)
				if(FALSE){
					fit$coefficients=unname(qr.coef(fit, tmp.zw))
					fit$residuals=qr.resid(fit, tmp.zw)
					fit$effects=qr.qty(fit, tmp.zw)
				}else{ ## did not check pivoting yet
					qq=qr.Q(fit, complete=TRUE)
					seqRank=seq_len(fit$rank)
					R=qr.R(fit)
					fit$effects=drop(crossprod(qq, tmp.zw))
					fit$coefficients=backsolve(R, fit$effects[seqRank])
					fit$residuals=drop(tmp.zw-qq[,seqRank,drop=FALSE]%*%fit$effects[seqRank])
				}
				fit$pivoted=!all.equal(seq_len(ncol(x)), fit$pivot)
				fit$tol=tmp.tol
				class(fit)='list'
			}else{	### this is faster way to do it, but still tooooo slow compared to C/FORTRAN code
				lmf=lm.fit(x[good,,drop=FALSE]*w, z*w, tol=min(1e-7, control$epsilon/1000))
				fit=lmf$qr
				fit[c('coefficients','residuals','effects')]=lmf[c('coefficients','residuals','effects')]
				fit$pivoted=!all.equal(seq_len(ncol(x)), fit$pivot)
				names(fit$effects)=seq_along(fit$effects)
				names(fit$coefficients)=NULL
				class(fit)='list'
			}
           # fit <- .Call(QuasiSeq:::Cdqrls_stats, x[good, , drop = FALSE] * w, z * w, min(1e-07, control$epsilon/1000), FALSE)
            if (any(!is.finite(fit$coefficients))) { ## LQ: this shouldn't happen
                conv <- FALSE
                warning(gettextf("non-finite coefficients at iteration %d", 
                  iter), domain = NA)
                break
            }
            if (nobs < fit$rank) 
                stop(sprintf(ngettext(nobs, "X matrix has rank %d, but only %d observation", 
                  "X matrix has rank %d, but only %d observations"), 
                  fit$rank, nobs), domain = NA)
            start[fit$pivot] <- fit$coefficients
            eta <- drop(x %*% start)
            mu <- linkinv(eta <- eta + offset)
            dev <- sum(dev.resids(y, mu, weights))
            if (control$trace) 
                cat("Deviance = ", dev, " Iterations - ", iter, 
                  "\n", sep = "")
            boundary <- FALSE
				class(fit)='qr'
				hats=rowSums(qr.Q(fit)^2)
				class(fit)='list'
            if (!is.finite(dev) || !is.finite(family$adjresponse(hats,mu,weights)) ) {
                if (is.null(coefold)) 
                  stop("no valid set of coefficients has been found: please supply starting values", 
                    call. = FALSE)
                warning("step size truncated due to divergence", 
                  call. = FALSE)
                ii <- 1
                while (!is.finite(dev)  || !is.finite(family$adjresponse(hats,mu,weights)) ) {
                  if (ii > max(1e2, control$maxit) )
                    stop("inner loop 1; cannot correct step size", 
                      call. = FALSE)
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- drop(x %*% start)
                  mu <- linkinv(eta <- eta + offset)
                  dev <- sum(dev.resids(y, mu, weights))
                }
                boundary <- TRUE
                if (control$trace) 
                  cat("Step halved: new deviance = ", dev, "\n", 
                    sep = "")
            }
            if (!(valideta(eta) && validmu(mu))) {
                if (is.null(coefold)) 
                  stop("no valid set of coefficients has been found: please supply starting values", 
                    call. = FALSE)
                warning("step size truncated: out of bounds", 
                  call. = FALSE)
                ii <- 1
                while (!(valideta(eta) && validmu(mu))) {
                  if (ii > max(1e2, control$maxit) ) 
                    stop("inner loop 2; cannot correct step size", 
                      call. = FALSE)
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- drop(x %*% start)
                  mu <- linkinv(eta <- eta + offset)
                }
                boundary <- TRUE
                dev <- sum(dev.resids(y, mu, weights))
                if (control$trace) 
                  cat("Step halved: new deviance = ", dev, "\n", 
                    sep = "")
            }
            if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
                conv <- TRUE
                coef <- start
                break
            }
            else {
                devold <- dev
                coef <- coefold <- start
            }
        }
        if (!conv) 
            warning("glm.fit0: algorithm did not converge", call. = FALSE)
        if (boundary) ### LQ: this shouldn't happen
            warning("glm.fit0: algorithm stopped at boundary value", 
                call. = FALSE)
        eps <- 10 * .Machine$double.eps
        if (family$family == "binomial") {
            if (any(mu > 1 - eps) || any(mu < eps)) ### LQ: this shouldn't happen
                warning("glm.fit0: fitted probabilities numerically 0 or 1 occurred", 
                  call. = FALSE)
        }
        if (family$family == "poisson") {
            if (any(mu < eps)) ### LQ: this shouldn't happen
                warning("glm.fit0: fitted rates numerically 0 occurred", 
                  call. = FALSE)
        }
		if (substr(family$family,1,17) == "Negative Binomial") {
            if (any(mu < eps)) ### LQ: this shouldn't happen
                warning("glm.fit0: fitted means numerically 0 occurred", 
                  call. = FALSE)
		}
        if (fit$rank < nvars) 
            coef[fit$pivot][seq.int(fit$rank + 1, nvars)] <- NA
        xxnames <- xnames[fit$pivot]
        residuals <- (y - mu)/mu.eta(eta)
        fit$qr <- as.matrix(fit$qr)
        nr <- min(sum(good), nvars)
        if (nr < nvars) {
            Rmat <- diag(nvars)
            Rmat[1L:nr, 1L:nvars] <- fit$qr[1L:nr, 1L:nvars]
        }
        else Rmat <- fit$qr[1L:nvars, 1L:nvars]
        Rmat <- as.matrix(Rmat)
        Rmat[row(Rmat) > col(Rmat)] <- 0
        names(coef) <- xnames
        colnames(fit$qr) <- xxnames
        dimnames(Rmat) <- list(xxnames, xxnames)
    }
    names(residuals) <- ynames
    names(mu) <- ynames
    names(eta) <- ynames
    wt <- rep.int(0, nobs)
    wt[good] <- w^2
    names(wt) <- ynames
    names(weights) <- ynames
    names(y) <- ynames
    if (!EMPTY) 
        names(fit$effects) <- c(xxnames[seq_len(fit$rank)], rep.int("", 
            sum(good) - fit$rank))
    wtdmu <- if (intercept) 
        sum(weights * y)/sum(weights)
    else linkinv(offset)
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok - as.integer(intercept)
    rank <- if (EMPTY) 
        0
    else fit$rank
    resdf <- n.ok - rank
    aic.model <- aic(y, n, mu, weights, dev) + 2 * rank
    list(coefficients = coef, residuals = residuals, fitted.values = mu, 
        effects = if (!EMPTY) fit$effects, R = if (!EMPTY) Rmat, 
        rank = rank, qr = if (!EMPTY) structure(fit[c("qr", "rank", 
            "qraux", "pivot", "tol")], class = "qr"), family = family, 
        linear.predictors = eta, deviance = dev, aic = aic.model, 
        null.deviance = nulldev, iter = iter, weights = wt, prior.weights = weights, 
        df.residual = resdf, df.null = nulldf, y = y, converged = conv, 
        boundary = boundary)
}