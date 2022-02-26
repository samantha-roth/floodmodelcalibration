// HetGP but with an SML level added into it
// the cheap code will simply be modelled by a homoscedastic GP 
// because I'm not massively interested in the noise here
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
	int<lower = 1> n_c;		// number cheap points
	int<lower = 2> n_e;		// number expensive points
	int<lower = 2> K;		// dimension of input space
	matrix[n_e, K] x_e;		// expensive input variables (should be studentised)
	matrix[n_c, K] x_c;		// expensive input variables (should be studentised)
	matrix[n_c, 2*m_p] m_H;		// design matrix (i.e. H = h(x)) of the mean, should be in block form
	vector[n_c] y_c;			// code outputs (noisy) - for now assume no replication ...
	vector[n_e] y_e;
	
	// now describe the priors of GP hyperparams


	// first the prior for the mean GP
	vector[m_p] m_beta_m;
	matrix[m_p,m_p] m_beta_s; // mean and sd of mean function parameters for the mean


	vector<lower = 0>[K] m_a_theta;
	vector<lower = 0>[K] m_b_theta; //correlation lengthscales (will use the form 1/theta^2)

	
	real<lower = 0> m_a_sigma; //params of dist for sigma^2 (mean)
	real<lower = 0> m_b_sigma;

	real<lower = 0> m_nugget_a; // quantifies the noise of the noise
	real<lower = 0> m_nugget_b;
	// we might set the nugget terms to be 10^(-4)

	// prior for the cheap mean GP
	
	vector[m_p] c_beta_m;
	matrix[m_p,m_p] c_beta_s; // mean and sd of mean function parameters for the mean


	vector<lower = 0>[K] c_a_theta;
	vector<lower = 0>[K] c_b_theta; //correlation lengthscales (will use the form 1/theta^2)

	
	real<lower = 0> c_a_sigma; //params of dist for sigma^2 (mean)
	real<lower = 0> c_b_sigma;

	real<lower = 0> c_nugget_a;	// quantifies the noise of the cheap mean
	real<lower = 0> c_nugget_b;	// just assume constant noise on this for simplicity	
	
	real m_rho;
	real<lower = 0> s_rho;
}
transformed data{
		
	vector[N] y; // all the code outputs concatenated
	vector[K] D_ee[n_e,n_e]; // their pairwise differences in block form
	vector[K] D_cc[n_c,n_c];
	vector[K] D_ce[n_c,n_e];

  vector[n_c] mu_c;
  vector[n_e] mu_e;
  
  mu_c = m_H[1:n_c,1:m_p]*c_beta_m;	//all rows, first m_p columns of m_H
  mu_e = m_H[(1 + n_c - n_e):n_c ,(1+m_p):2*m_p] * m_beta_m ; 	// m_H is in block form!
	D_ee = diffs(x_e, x_e,  n_e, n_e, K);
	D_cc = diffs(x_c, x_c, n_c, n_c, K);
	D_ce = diffs(x_c, x_e, n_c, n_e, K);
	
	y = append_row(y_c, y_e);	

}


parameters {



real<lower = 0> m_sigma;	// exp mean scale param
real<lower = 0> c_sigma;	// cheap mean scale param

vector<lower = 0>[K] m_theta; 	// length scale parameters
vector<lower = 0>[K] c_theta;

real<lower = 0> m_nugget;	// variance of the log variance
real<lower = 0> c_nugget;	// variance of the cheap outputs

real rho;

}

model{

vector[N] mu;
vector[n_c] c_mu;
vector[n_e] m_mu;
matrix[N, N] Sigma_mat;
matrix[n_e, n_e] Var_Y_e;
matrix[n_c, n_c] Var_Y_c;
matrix[n_e, n_c] Cov_c_e;
matrix[N, N] var_cholesky;
matrix[N, 2*m_p] m_H2;

// m_H is in block form!
// append 1/2 the columns of m_H to a n_c x m_p matrix of zeros,
// then stack this on top of the last n_e rows of m_H
m_H2 = append_row(append_col(m_H[,1:m_p], to_matrix(rep_array(0,n_c, m_p))) , m_H[(1+n_c-n_e):n_c,] );
// m_H[,1:m_p]: design matrix for cheap runs [1,n_fp,n_ch]
// to_matrix(rep_array(0,n_c, m_p)): matrix of zeros of same size as m_H[,1:m_p]
// m_H[(1+n_c-n_e):n_c,] is just the bottom portion of m_H containing the rows 
// that contain the parameter settings for the cheap runs 

print(dims(m_H2));
m_H2[(1+n_c):N,1:m_p] = rho*m_H2[(1+n_c):N,1:m_p]; 
// in the bottom of mH_2, multiply the first m_p columns by rho

c_mu = m_H[1:n_c,1:m_p]*c_beta_m;	//m_c(X_c)
//the first n_c rows and m_p columns of m_H are weighted by the mean of the c betas

m_mu = m_H[(1+n_c-n_e):n_c ,(1+m_p):2*m_p] * (m_beta_m + rho*c_beta_m); 	
//corresponds to rho*m_c(X_e) + m_e(X_e)
// m_H is in block form!
//the last n_e rows and the last m_p columns of m_mu 
//are weighted by the sum of the mean of m betas and rho times the mean of the c betas

mu = append_row(c_mu, m_mu);
//corresponds to rbind(m_c(X_c), rho*m_c(X_e) + m_e(X_e))
//see page 10 of Kennedy et al 2020

// construct the variance matrix
Var_Y_c = exp_kern(D_cc, square(c_sigma), exp(-2*log(c_theta)), n_c, n_c, K);
// C(x_c, x_c)

// I can be clever and order the x_c such that the intersection with x_e is at end
Var_Y_e = (rho * rho * Var_Y_c[(1+n_c-n_e):n_c,(1+n_c-n_e):n_c] ) + exp_kern(D_ee, square(m_sigma), exp(-2*log(m_theta)), n_e, n_e, K) +  diag_matrix(rep_vector(m_nugget, n_e)); 
//Var(Y_e)= rho^2 * part of Var_Y_c corresponding to x_e + C(x_e, x_e) + diag(lambda)

Cov_c_e = rho * exp_kern(D_ce, square(c_sigma), exp(-2*log(c_theta)), n_c, n_e, K)' ;
//covariance between y_c and y_e

Var_Y_c = Var_Y_c + diag_matrix(rep_vector(c_nugget, n_c));
//Var(Y_c)= C(x_c, x_c) + diag(c_nugget)

Sigma_mat = m_H2 * diag_matrix(append_row(diagonal(c_beta_s), diagonal(m_beta_s))) * (m_H2)' +  append_row( append_col(Var_Y_c, Cov_c_e') , append_col(Cov_c_e, Var_Y_e) ); 
//part to account for betas not being modeled separated + Var(Y)

var_cholesky = cholesky_decompose(Sigma_mat);

y ~ multi_normal_cholesky(mu, var_cholesky); // statement about the observed code outputs
// the above is <kind of> two GPs

// prior beliefs

//for(i in 1:m_p){
//	m_beta[i] ~ normal(m_beta_m[i], m_beta_s[i]);
//}

//for(i in 1:m_p){
//	c_beta[i] ~ normal(c_beta_m[i], c_beta_s[i]);
//}

m_theta ~ gamma(m_a_theta, m_b_theta);
c_theta ~ gamma(c_a_theta, c_b_theta);

m_sigma ~ inv_gamma(m_a_sigma, m_b_sigma);
c_sigma ~ inv_gamma(c_a_sigma, c_b_sigma);

c_nugget ~ inv_gamma(c_nugget_a, c_nugget_b);
m_nugget ~ inv_gamma(m_nugget_a, m_nugget_b);

rho ~ normal(m_rho, s_rho); // a priori we expect the codes to be giving a similar output


}

