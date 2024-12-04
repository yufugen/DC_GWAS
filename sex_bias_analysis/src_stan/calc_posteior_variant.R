#' Calculate posterior for each variant
#'
#' @export
#' @param sum.prop, the estimated proportion of each component as a vector, posterior mean or median
#' @param sigmasq, the estimated sigmasq
#' @param dat, a list with B and SESQ elements scaled with allele frequency
#' @param alpha, a vector of four elements describe the relationship between male and female
#' @return a data frame with the posteriors for each of the variant

calc_posteriors <- function(sum.prop, sigmasq, dat, alpha){
  B.dat <- dat$B
  SESQ.dat <- dat$SESQ
  scalesq_vector <- dat$scalesq_vector
  N <- nrow(B.dat)
  p <- sum.prop
  
  posteriors <- lapply(1:N, function(i)
    .calc_posterior_variant(B.dat[i,], SESQ.dat[i,], p, sigmasq, alpha))
  posterior.df <- data.frame(do.call(rbind, posteriors))
  colnames(posterior.df) <- paste0("p",0:3)
  
  posterior.df$ID <- dat$ID
  return(posterior.df)
}

.calc_posterior_variant <- function(B, SESQ, p, sigmasq, alpha){
  zeros <- c(0,0)
  SESQ_mat <- matrix(c(SESQ[1], 0, 0, SESQ[2]), 2, 2)
  p_0 = p[1]*mnormt::dmnorm(B, zeros, SESQ_mat)
  p_1 = p[2]*mnormt::dmnorm(B, zeros, SESQ_mat +
                              matrix(c(sigmasq*alpha^2, sigmasq*alpha, sigmasq*alpha,  sigmasq),2,2))
  p_2 = p[3]*mnormt::dmnorm(B, zeros, SESQ_mat +
                              matrix(c(sigmasq, sigmasq, sigmasq, sigmasq),2,2))
  p_3 = p[4]*mnormt::dmnorm(B, zeros, SESQ_mat +
                              matrix(c(sigmasq, sigmasq*alpha, sigmasq*alpha, sigmasq*alpha^2),2,2))
  p_tot = p_0 + p_1 + p_2 + p_3
  prob_0 = exp(log(p_0) - log(p_tot))
  prob_1 = exp(log(p_1) - log(p_tot))
  prob_2 = exp(log(p_2) - log(p_tot))
  prob_3 = exp(log(p_3) - log(p_tot))
  posteriors <- data.frame(t(c(prob_0, prob_1, prob_2, prob_3)))
  #posterior.df <- data.frame(do.call(rbind, posteriors))
  #colnames(posterior.df) <- c("p0", "p1", "p2", "p3")
  return(posteriors)
}

