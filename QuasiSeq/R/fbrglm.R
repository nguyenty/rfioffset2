###########################################
### negbin.br, glm.fit0, and fbrglm were authored by Long Qu <long.qu@wright.edu> and are used during the Firth bias correction
### for genes with low counts and extreme coefficient estimates
###########################################

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


