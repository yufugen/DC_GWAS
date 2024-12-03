#' Return the fitted parameters for a model object
#' 
#' Note that in all cases these returns the summary statistics from the post-warmup draws.
#' 
#'  @export
#'  @param fit fitted stan model
#'  @return named list of fit parameters. 
#'  Item 1 provides full summary of pi. 
#'  Item 2 provides the median value for pi, which is used for posterior calculation

extract_fit_pi <- function(fit){
  
  fit_sum_pi <- rstan::summary(fit, pars = c("pi"), probs = c(0.025, 0.50, 0.975))$summary
  rownames(fit_sum_pi) <- sapply(0:(nrow(fit_sum_pi)-1), function(i) sprintf("pi[%i]", i))
  p <- fit_sum_pi[,c("50%")]
  names(p) <- sapply(0:(length(p)-1), function(i) sprintf("pi[%i]",i))
  return(list("summary_prop" = fit_sum_pi,
              "prop" = p))
}

extract_fit_vars <- function(fit){
	
	fitS <- rstan::summary(fit, pars = c("sigmasq"), probs = c(0.025, 0.50, 0.975))$summary
	rownames(fitS) <- sapply(1:(nrow(fitS)), function(i) sprintf("sigmasq[%i]", i))
	s <- fitS[,c("50%")]
	names(s) <- sapply(1:(length(s)), function(i) sprintf("sigmasq[%i]",i))
	return(list("summary_sigmasq" = fitS,
		    "sigmasq" = s))
}
