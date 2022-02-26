// stan program for hetGP
// useful for when the noise depends on the input
// e.g. sin(x) + eps_1 + eps_2*x has noise that increases with |x|
  // this version assumes the input, X, has dimension K>1
functions {
  
  vector[,] diffs( matrix X1, matrix X2, int n1, int n2, int K ) {
    
    // X1 might be the cheap design
    // X2 might be the expensive design
    // or they might be the same
    // Ni is the size of design i
    
    
    vector[K] D[n1,n2];
    vector[K] tmp1;
    vector[K] tmp2;
    
    
    
    for(i in 1:n1){
      for(j in 1:n2){
        tmp1 = (X1[i,]' - X2[j,]');
        D[i,j] = tmp1 .* tmp1;
      }
    }
    
    return(D);
  }
  
  matrix exp_kern( vector[,] D, real sq_sigma, vector omega, int n1, int n2,    int K  ) {
    
    matrix[n1,n2] mat;
    vector[K] tmp1;
    real tmp2;
    
    
    int ii = 0;
    
    mat = rep_matrix(0, n1, n2);
    if(n1 != n2){
      for(i in 1:n1){
        
        //ii = ii + 1; // ii = 1:(N-1)
        
        for(j in 1:n2){ // just compute the lower triangle
          
          tmp1 = D[i,j];	
          
          mat[i,j] =  (-1)*dot_product( omega, tmp1 );
          
        }
        
        
      }
    }
    if(n1 == n2){
      for(i in 1:n1){
        
        ii = ii + 1; // ii = 1:(N-1)
        
        for(j in ii:n2){ // just compute the lower triangle
          
          tmp1 = D[i,j];	
          
          mat[i,j] =  (-1)*dot_product( omega, tmp1 );
          
        }
        
        
      }
      
      mat = mat + mat';

		}


		mat = sq_sigma*exp(mat);

		
		return(mat);
	}

}

data {
	
	int<lower = 1> m_p;		// number regression functions for the mean
	int<lower = 1> N;		// number data points
	int<lower = 2> K;		// dimension of input space
	matrix[N, K] x;			// input data (should be studentised)
	matrix[N, m_p] m_H;		// design matrix (i.e. H = h(x)) (mean)
	vector[N] y;			// code outputs (noisy) - for now assume no replication ...
	
	vector<lower = 1>[N] a;		// vector of replication level (we can have different levels of replication here, unlike pairs of homGP)
	// now describe the prior


	// first the prior for the mean GP
	vector[m_p] m_beta_m;
	matrix[m_p,m_p] m_beta_s; // mean and sd of mean function parameters for the mean


	vector<lower = 0>[K] m_a_theta;
	vector<lower = 0>[K] m_b_theta; //correlation lengthscales (will use the form 1/theta^2)

	
	real<lower = 0> m_a_sigma; //params of dist for sigma^2 (mean)
	real<lower = 0> m_b_sigma;
	
	real<lower = 0> m_nugget_a; // quantifies the noise of the noise
	real<lower = 0> m_nugget_b;
}

transformed data{
  vector[N] mu;			// the "actual" mean of the top-level GP
  matrix[N,N] var_mean;
	vector[K] Diffs[N, N];
	Diffs = diffs(x, x, N, N, K);
  
  mu = m_H * m_beta_m;
  var_mean = m_H * m_beta_s * ( m_H');
    }

parameters {
  
  real<lower = 0> m_sigma;	// mean scale param
  
  vector<lower = 0>[K] m_theta; 	// length scale parameters
  
  real<lower = 0> m_nugget;
  
}

model {
  
  matrix[N, N] m_var;		// variance matrix of the mean
  
  // first produce the variance
  
  m_var = var_mean + exp_kern(Diffs, square(m_sigma), exp(-2*log(m_theta)), N, N, K) + diag_matrix(rep_vector(m_nugget, N));
  
  y ~ multi_normal_cholesky(mu, cholesky_decompose(m_var)); // top level statement
  
  // prior beliefs
  
  //for(i in 1:m_p){
    //	m_beta[i] ~ normal(m_beta_m[i], m_beta_s[i]);
    //}
  
  
  for(k in 1:K){
    m_theta ~ gamma(m_a_theta[k], m_b_theta[k]);
  }
  
  m_sigma ~ inv_gamma(m_a_sigma, m_b_sigma);
  
  m_nugget ~ inv_gamma(m_nugget_a, m_nugget_b);
  
}


