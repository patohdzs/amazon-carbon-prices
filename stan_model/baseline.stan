data {
  // Gamma regression
  int<lower=1> N_gamma; // total number of observations
  int<lower=1> K_gamma; // number of predictors
  int<lower=1> M_gamma; // number of groups
  matrix[N_gamma, K_gamma] X_gamma; // design matrix
  vector[N_gamma] y_gamma; // outcome variable
  array[N_gamma] int m_gamma; // group map

  // Theta regression
  int<lower=1> N_theta; // total number of observations
  int<lower=1> K_theta; // number of predictors
  int<lower=1> M_theta; // number of groups
  matrix[N_theta, K_theta] X_theta; // design matrix
  vector[N_theta] y_theta; // outcome variable
  array[N_theta] int m_theta; // group map
  vector[N_theta] W_theta; // weights

  // Gamma projection
  int<lower=1> num_sites; // number of samples
  matrix[num_sites, K_gamma] X_gamma_fit; // design matrix
  array[num_sites] int m_gamma_fit; // group map

  // Theta projection
  int<lower=1> C_theta_fit; // Number of municipalities
  matrix[C_theta_fit, K_theta] X_theta_fit; // design matrix
  array[C_theta_fit] int m_theta_fit; // group map

  // Sparse G_theta
  int<lower=0> G_nnz_theta; // Number of non-zero elements
  vector[G_nnz_theta] G_w_theta; // Non-zero values
  array[G_nnz_theta] int G_v_theta; // Column indices (1-based in Stan)
  array[num_sites + 1] int G_u_theta; // Row pointers (1-based in Stan)
  real pa_2017; // price of cattle in 2017
}
transformed data {
  vector[N_theta] y_theta_w = W_theta .* y_theta;
  matrix[N_theta, K_theta] X_theta_w;
  for (n in 1 : N_theta)
    X_theta_w[n] = W_theta[n] * X_theta[n];
}
parameters {
  // Gamma regression parameters
  vector[K_gamma] beta_gamma;
  vector[M_gamma] nu_gamma;
  real log_precision_u_gamma;
  real log_precision_v_gamma;

  // Theta regression parameters
  vector[K_theta] beta_theta;
  vector[M_theta] nu_theta;
  real log_precision_u_theta;
  real log_precision_v_theta;
}
transformed parameters {
  // Pre-multiply theta FE's by weights
  vector[N_theta] nu_theta_w = W_theta .* nu_theta[m_theta];
  vector[C_theta_fit] nu_theta_fit = nu_theta[m_theta_fit];
  vector[N_gamma] nu_gamma_sort = nu_gamma[m_gamma];
  vector[num_sites] nu_gamma_fit = nu_gamma[m_gamma_fit];

  real sigma_u_gamma = exp(-0.5 * log_precision_u_gamma);
  real sigma_v_gamma = exp(-0.5 * log_precision_v_gamma);
  real sigma_u_theta = exp(-0.5 * log_precision_u_theta);
  real sigma_v_theta = exp(-0.5 * log_precision_v_theta);

  // Projection
  vector<lower=0>[num_sites] gamma = exp(X_gamma_fit * beta_gamma
                                         + nu_gamma_fit);

  vector[C_theta_fit] exp_log_theta = exp(X_theta_fit * beta_theta
                                          + nu_theta_fit);

  vector<lower=0>[num_sites] theta = csr_matrix_times_vector(num_sites,
                                                             C_theta_fit,
                                                             G_w_theta,
                                                             G_v_theta,
                                                             G_u_theta,
                                                             exp_log_theta)
                                     / pa_2017;

  // for (n in 1 : N_nonzero_G_theta) {
  //   theta[row_G_theta[n]] += val_G_theta[n] * exp_log_theta[col_G_theta[n]]
  //                            / pa_2017;
  // }
}
model {
  // Priors
  nu_gamma ~ normal(0, sigma_v_gamma);
  nu_theta ~ normal(0, sigma_v_theta);

  y_gamma ~ normal(X_gamma * beta_gamma + nu_gamma_sort, sigma_u_gamma);
  y_theta_w ~ normal(X_theta_w * beta_theta + nu_theta_w, sigma_u_theta);
}
