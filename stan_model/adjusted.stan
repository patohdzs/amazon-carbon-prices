functions {
  real log_value(vector gamma, vector theta, int T, int S, matrix Z,
                 vector forest_area_2017, vector alpha_p_Adym, matrix Bdym,
                 vector ds_vect, real xi, real pa, real pe,
                 matrix carbon_stock) {
    // Compute stock of carbon (X)
    row_vector[T] omega = gamma' * carbon_stock;
    vector[T + 1] X;
    X[1] = gamma' * forest_area_2017;
    X[2 : (T + 1)] = alpha_p_Adym * X[1] + Bdym * omega';

    // Compute aggregate X dot
    vector[T] Xdot_agg = (X[2 : (T + 1)] - X[1 : T]);
    real term_1 = -pe * sum(ds_vect .* (-Xdot_agg));

    // Value of agricultural output
    vector[T + 1] agri_output = pa * (theta' * Z)';
    real term_2 = sum(ds_vect .* agri_output[2 : T + 1]);

    // Overall objective value
    real obj_val = term_1 + term_2;
    real log_density_val = -1.0 / xi * obj_val;

    return log_density_val;
  }
}
data {
  // Planner problem
  int<lower=0> T; // Time horizon
  int<lower=0> S; // Number of sites
  matrix[S, T + 1] Z; // Agricultural area state
  matrix[S, T] U; // U control
  matrix[S, T] V; // V control
  vector[S] zbar_2017; // z_bar in 2017
  vector[S] forest_area_2017; // forest area in 2017
  vector[T] alpha_p_Adym;
  matrix[T, T] Bdym;
  vector[T] ds_vect; // Time discounting vector
  matrix[S, T] carbon_stock;
  real<lower=0> alpha; // Mean-reversion coefficient
  real zeta_u; // Penalty on adjustment costs
  real zeta_v; // Penalty on adjustment costs
  real xi; // Penalty on prior-posterior KL div
  real kappa; // Effect of cattle farming on emissions
  real<lower=0> pa; // Price of cattle output
  real<lower=0> pe; // Price of carbon emission transfers

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
  vector[M_gamma] nu_gamma_transform;
  real log_precision_u_gamma;
  real log_precision_v_gamma;

  // Theta regression parameters
  vector[K_theta] beta_theta;
  vector[M_theta] nu_theta_transform;
  real log_precision_u_theta;
  real log_precision_v_theta;
  // real<lower=log(1e-5), upper=log(1e5)> log_precision_v_theta;
}
transformed parameters {
  real sigma_u_gamma = exp(-0.5 * log_precision_u_gamma);
  real sigma_v_gamma = exp(-0.5 * log_precision_v_gamma);
  real sigma_u_theta = exp(-0.5 * log_precision_u_theta);
  real sigma_v_theta = exp(-0.5 * log_precision_v_theta);

  vector[M_gamma] nu_gamma = sigma_v_gamma * nu_gamma_transform;
  vector[M_theta] nu_theta = sigma_v_theta * nu_theta_transform;

  // Pre-multiply theta FE's by weights
  vector[N_theta] nu_theta_w = W_theta .* nu_theta[m_theta];
  vector[C_theta_fit] nu_theta_fit = nu_theta[m_theta_fit];
  vector[N_gamma] nu_gamma_sort = nu_gamma[m_gamma];
  vector[num_sites] nu_gamma_fit = nu_gamma[m_gamma_fit];

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
}
model {

  nu_gamma_transform ~ normal(0, 1);
  nu_theta_transform ~ normal(0, 1);

  y_gamma ~ normal(X_gamma * beta_gamma + nu_gamma_sort, sigma_u_gamma);
  y_theta_w ~ normal(X_theta_w * beta_theta + nu_theta_w, sigma_u_theta);

  target += log_value(gamma, theta, T, S, Z, forest_area_2017, alpha_p_Adym,
                      Bdym, ds_vect, xi, pa, pe, carbon_stock);
}
