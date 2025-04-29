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
  matrix[N_theta, N_theta] W_theta; // weights

  // Gamma projection
  int<lower=1> num_sites; // number of samples
  matrix[num_sites, K_gamma] X_gamma_fit; // design matrix
  array[num_sites] int m_gamma_fit; // group map

  // Theta projection
  int<lower=1> C_theta; // Number of municipality-site combinations
  matrix[C_theta, K_theta] X_theta_fit; // design matrix
  array[C_theta] int m_theta_fit; // group map
  matrix[num_sites, C_theta] G_theta_fit; // municipality projection
  real pa_2017; // price of cattle in 2017
}

transformed data {
  vector[N_theta] y_theta_weighted;
  matrix[N_theta, K_theta] X_theta_weighted;

  y_theta_weighted = W_theta * y_theta;
  X_theta_weighted = W_theta * X_theta;
}

parameters {
  // Gamma regression parameters
  vector[K_gamma] beta_gamma;
  vector[M_gamma] nu_gamma;
  real<lower=0> sigma_u_gamma;
  real<lower=0> sigma_v_gamma;

  // Theta regression parameters
  vector[K_theta] beta_theta;
  vector[M_theta] nu_theta;
  real<lower=0> sigma_u_theta;
  real<lower=0> sigma_v_theta;
}

transformed parameters {
  vector<lower=0>[num_sites] gamma = exp(X_gamma_fit * beta_gamma
                                         + nu_gamma[m_gamma_fit]);
  vector<lower=0>[num_sites] theta = (G_theta_fit
                                      * exp(X_theta_fit * beta_theta
                                            + nu_theta[m_theta_fit]))
                                     / pa_2017;
}
model {
  // Priors
  beta_gamma ~ normal(0, sigma_u_gamma);
  beta_theta ~ normal(0, sigma_u_theta);
  nu_gamma ~ normal(0, sigma_v_gamma);
  nu_theta ~ normal(0, sigma_v_theta);

  vector[N_theta] nu_theta_weighted=W_theta*nu_theta[m_theta];

  y_gamma ~ normal(X_gamma * beta_gamma + nu_gamma[m_gamma], sigma_u_gamma);
  y_theta_weighted ~ normal(X_theta_weighted * beta_theta + nu_theta[m_theta], sigma_u_theta);
}
