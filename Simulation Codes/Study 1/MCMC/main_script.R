# rm(list=ls(all=TRUE)) # Discouraged in shared scripts; clear environment manually if needed

library("fda")        
library("fields")
library("MASS")     
library("ggplot2")
library("dplyr")
library("tidyr")
library("posterior")
library("bayesplot")
library("ggrepel")
library("coda")
library("cmdstanr")
library("snowfall")

set.seed(12345)     # Seed for reproducibility

# Set to relative paths for GitHub compatibility
# Assuming the script runs from the root of your project directory
setwd("./Data")
load("Y_obs.RData")
load("d.RData")
load("times.RData")
load("C1.RData")
load("C2.RData")

int  = 250
cpus = 4
m    = 25           # Number of spatial locations
n    = 200 
D    = d

# Functional bases
G = 10              # Number of bases for intercept β0(t)
R = 10              # Number of bases for each βp(t), p= 1,2
V = 10              # Number of bases for spatial effect w_s(t)

###############################################################################
# Construction of B-spline bases for β0(t), β1(t), β2(t) and w_s(t)
###############################################################################

# 3.1. Intercept β0(t)
basis_beta0 = create.bspline.basis(rangeval = c(0,1), nbasis = G, norder = 4)
Phi0 = eval.basis(times, basis_beta0)    # n x G matrix

# 3.2. Coefficients β1(t) and β2(t)
basis_beta1 = create.bspline.basis(rangeval = c(0,1), nbasis = R, norder = 4)
basis_beta2 = create.bspline.basis(rangeval = c(0,1), nbasis = R, norder = 4)

Phi1 = eval.basis(times, basis_beta1)    # n x R matrix
Phi2 = eval.basis(times, basis_beta2)    # n x R matrix

# 3.3. Spatial effect w_s(t)
basis_w = create.bspline.basis(rangeval = c(0,1), nbasis = V, norder = 4)
Phi_w = eval.basis(times, basis_w)       # n x V matrix

#-------------------------------------------------------------------------------
model_spatial_functional = "
functions {
  // Matérn 3/2 Function
  matrix matern32_cov(matrix D, real kappa2, real varphi) {
    int M = dims(D)[1];            // number of sites (rows of D)
    matrix[M, M] cov;
    real sqrt3 = sqrt(3.0);
    
    for (i in 1:M) {
      for (j in i:M) {
        if (i == j) {
          // Diagonal: add a jitter to ensure positive-definiteness
          cov[i, j] = kappa2 + 1e-10;
        } else {
          real h = D[i, j];
          cov[i, j] = kappa2 * (1 + sqrt3 * h / varphi) * exp(-sqrt3 * h / varphi);
          cov[j, i] = cov[i, j];   // ensure symmetry
        }
      }
    }
    return cov;
  }
}

data {
  int<lower=1> m;            // number of locations 
  int<lower=1> n;            // number of time points 
  int<lower=1> G;            // number of B-spline bases for the q term 
  int<lower=1> R;            // number of bases for beta1(t) and beta2(t)
  int<lower=1> V;            // number of B-spline bases for the spatial term
  
  matrix[n, G] Phi0;         // intercept basis matrix
  matrix[n, R] Phi1;         // basis matrix for Beta 1
  matrix[n, R] Phi2;         // basis matrix for Beta 2 
  matrix[n, V] Phi_w;        // temporal basis matrix for spatial effect
  
  vector[m] C1;              // first covariate (e.g., longitude)
  vector[m] C2;              // second covariate (e.g., latitude)
  matrix[m, m] D;            // Matrix of (Euclidean) distances between all pairs of sites [m x m]
  
  matrix[m, n] Y;            // simulated (observed) data
}

parameters {
  vector[G] q;               // global coefficients
  vector[R] z1;              // coefficients for spatial covariate p = 1
  vector[R] z2;              // coefficients for spatial covariate p = 2 
  matrix[m, V] z_E;          // each column ~ Normal(0, I_m)
  
  real<lower=0> kappa2;      // scale (variance) for the spatial effect
  real<lower=0> varphi;      // range of the Matern
  real<lower=0> sigma2_q;    // variance for q
  real<lower=0> sigma2_z;    // variance for z1 and z2
  real<lower=0> tau2;        // noise variance
}
  
transformed parameters {
  matrix[m, V] E;
  matrix[m, m] Cov_w;                          // spatial covariance matrix
  Cov_w = matern32_cov(D, kappa2, varphi);
  matrix[m, m] L = cholesky_decompose(Cov_w);  // Cholesky factor of Sigma
  E = L * z_E;
}
  
model {
  // ========== Priors ==========
  
  sigma2_q ~ inv_gamma(2,1.2);       // q ~ Normal(0, sigma2_q)
  sigma2_z ~ inv_gamma(2,0.5);       // z1, z2 ~ Normal(0, sigma2_z)
  tau2     ~ inv_gamma(2,0.5);   
  kappa2   ~ inv_gamma(2,1);
  varphi   ~ inv_gamma(2,1);
  
  q        ~ normal(0, sqrt(sigma2_q));
  z1       ~ normal(0, sqrt(sigma2_z));
  z2       ~ normal(0, sqrt(sigma2_z));
  
  to_vector(z_E) ~ normal(0, 1);  // Non-centered prior for z_E

  // ========== Likelihood ==========
 
  for (j in 1:m) {
    vector[n] mu_j;
    mu_j = Phi0  * q
         + C1[j] * (Phi1 * z1)
         + C2[j] * (Phi2 * z2)
         + Phi_w * E[j]';         // now E[j]' has dimension V x 1
    Y[j] ~ normal(mu_j, sqrt(tau2)); // broadcasting in n dimensions
  }
}
"
#-------------------------------------------------------------------------------

Model = write_stan_file(model_spatial_functional)
mod = cmdstan_model(Model)

counter = 0

MonteCarlo = function(step){
  
  # Define the name of the unique log file (Using relative paths for GitHub)
  log_file = "../Output/MonteCarlo_log.txt"  
  
  log_msg <- function(msg) {
    cat(msg, file = log_file, append = TRUE, sep = "\n")  # Writes the message to the file
    flush.console()                                       # Flushes the console and file immediately
  }
  
  ##
  Y = Y_obs %>%
      filter(id == step) %>%
      dplyr::select(dplyr::starts_with("V")) %>%
      as.matrix() %>%
      t()
  
  data_list <- list(
    m     = m,
    n     = n,
    G     = G,
    R     = R,
    V     = V,
    Phi0  = Phi0,
    Phi1  = Phi1,
    Phi2  = Phi2,
    Phi_w = Phi_w,
    C1    = C1,
    C2    = C2,
    D     = D,
    Y     = Y
  )
  
  # 5) Compile and run the Stan model
  fit_model = mod$sample(
    data            = data_list,
    seed            = 123,
    chains          = 2,
    parallel_chains = 2,
    iter_warmup     = 20000,
    iter_sampling   = 10000,
    thin            = 10
  )
  ##
  
  draws = fit_model$draws(c("q", "z1", "z2", "z_E", "kappa2" ,"varphi", "sigma2_q", "sigma2_z", "tau2"))
  variable = as_draws_matrix(draws)
  diagnostic = summarise_draws(variable, "mean","median", "sd", default_convergence_measures())
  id = rep(step, nrow(diagnostic))
  diagnostic = cbind(diagnostic, id)
  
  ###########
  ### Hpd ###
  ###########
  
  Hpd_mcmc = as.mcmc(variable)
  Hpd = HPDinterval(Hpd_mcmc, prob=0.95)
  id2 = rep(step, nrow(Hpd))
  Hpd = cbind(Hpd, id2)
  
  write.table(diagnostic, file="../Output/diagnostic.txt",
              sep = ",",
              col.names = FALSE,
              row.names = FALSE,
              append = TRUE)
  
  write.table(Hpd, file="../Output/Hpd.txt",
              sep = ",",
              col.names = FALSE,
              row.names = FALSE,
              append = TRUE)
  
  ###
  # Increments the global counter to track how many iterations were completed
  counter <<- counter + 1
  
  # Logs progress in the log file
  log_msg(paste("Finishing iteration:", step, "- Total completed iterations:", counter, "\n"))
}

sfInit(parallel=TRUE, cpus=cpus)
sfLibrary(dplyr)
sfLibrary(wavethresh)
sfLibrary(cmdstanr)
sfLibrary(posterior)
sfLibrary(bayesplot)
sfLibrary(coda)
sfExportAll()
sfLapply(1:int, fun=MonteCarlo) # Function that I want to compute multiple times using sfLapply
sfStop()
