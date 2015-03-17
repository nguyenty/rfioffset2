% File src/library/QuasiSeq/man/NBDev.Rd
     \name{NBDev}
     \alias{NBDev}
     \title{Fit a negative binomial GLM for given design matrix}
     \description{A function called within \code{QL.fit} to fit a negative binomial GLM of each gene for given design matrix}
     \usage{
     NBDev(counts,design,log.offset,nb.disp,print.progress=TRUE)
     }
     \arguments{
      \item{counts}{RNA-seq data matrix of integer expression counts.  Each row contains observations from a single gene. Each column contains observations from a single experimental unit.}
      \item{design}{A single element from the \code{design.list} argument given to \code{QL.fit}.}
	\item{log.offset}{A vector of log-scale, additive factors used to adjust estimated log-scale means for differences in library sizes across samples.  Commonly used offsets include,\code{log.offset=log(colSums(counts))} or \code{log.offset=log(apply(counts,2,quantile,.75))}.  The default setting  in \code{QLfit} makes no adjustment for library sizes (i.e. log.offset=0).}
	\item{nb.disp}{estimated negative binomial dispersion parameters obtained from either \code{estimateGLMTrendedDisp} or \code{estimateGLMCommonDisp} in package \code{edgeR}.  These estimates are treated as known and are used to compute deviances.}
	\item{print.progress}{logical. If TRUE, the function will provide an update on what gene (row number) is being analyzed.  Updates occur frequently to start then eventually occur every 5000 genes.}
}
	
\value{list containing:
\item{dev}{vector containing the deviance for each gene under a negative binomial model fit to design matrix specified by \code{design}. This vector is passed along within the \code{QL.fit} function.}
  \item{means}{matrix of fitted mean values for each gene}
\item{parms}{matrix of estimated coefficients for each gene. Note that these are given on the log scale. (i.e. intercept coefficient reports log(average count) and non-intercept coefficients report estimated log fold-changes.)  Genes with at least one zero count and initial absolute coefficient estimates greater than three undergo the bias correction of Kosmidis & Firth (2009) to moderate extreme coefficient estimates.}	
}

\references{
Kosmidis and Firth (2009) "Bias reduction in exponential family nonlinear models" \emph{Biometrika}, \bold{96}, 793--804.

Lund, Nettleton, McCarthy and Smyth (2012) "Detecting differential expression in RNA-sequence data using quasi-likelihood with shrunken dispersion estimates" emph{SAGMB}, \bold{11}(5).
}


\author{Steve Lund \email{lundsp@gmail.com}}

     \keyword{RNA-seq, quasi-likelihood, differential expression, negative binomial}




