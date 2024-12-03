data{
  int<lower=1> K; // number of mixture components, should be 4 in this case. <lower=1> is contstrain to have a lower bond at 1.
  int<lower=1> N; // number of samples, i.e. the number of SNPs being examined
  int<lower=1> M; // dimensions, should be 2

  matrix[N, M] B; // scaled betas from observed data
  matrix[N, M] SESQ; // scaled squared SE from observed data

  //vector[K] alphas; // a vector to define the relationship between beta_m and beta_f
  real<lower=0> alpha;
}

transformed data{ // transformed data allows for preprocessing the data. This step needed?
  vector[M] zeros;
  zeros = rep_vector(0, M) ;
}

parameters{
  simplex[K] pi; // a vector with non-negative values wholse entries sum to 1
  //vector<lower=0>[3] sigmasq;
  real<lower=0> sigmasq;
}

transformed parameters{
    matrix[M, M] Sigma[K]; // for each K, define a 2x2 Sigma matrix
    vector[2] a;
    vector[2] b;
    
    a[1] = alpha;
    a[2] = 1;
    b[1] = 1;
    b[2] = alpha;

    Sigma[1] = diag_matrix(rep_vector(0,2));
    Sigma[2] = sigmasq*[[a[1]^2, a[1]*b[1]], [a[1]*b[1], b[1]^2]];
    Sigma[3] = sigmasq*[[1, 1], [1, 1]];
    Sigma[4] = sigmasq*[[a[2]^2, a[2]*b[2]], [a[2]*b[2], b[2]^2]];

}

model{
  vector[K] ps; // contribution of each
  sigmasq ~ uniform(0,1);   // priors for sigmasq
  pi ~ dirichlet(rep_vector(inv(K), K)) ;
  
  for (n in 1:N){
    for (k in 1:K) {
      ps[k] = log(pi[k]) + multi_normal_lpdf(B[n] | zeros, diag_matrix(to_vector(SESQ[n])) + Sigma[k]) ;
    }
    target += log_sum_exp(ps); 
  }
}
