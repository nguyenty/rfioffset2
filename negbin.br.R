
###########################################
### negbin.br, glm.fit0, and fbrglm were authored by Long Qu <long.qu@wright.edu> and are used during the Firth bias correction
### for genes with low counts and extreme coefficient estimates
###########################################

### This is only a modified version of stats::glm.fit(). This is not to be directly used by any user.
negbin.br=function( link = "log", overdisp = stop("'overdisp' must be specified"))
{
	mu.eta=variance=linkfun=function(...){return('Shut up CRAN checker.')}
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

