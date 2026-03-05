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

d_lower = unname( quantile(D, probs = 0.05) )
d_upper = unname( quantile(D, probs = 0.95) )

#-------------------------------------------------------------------------------
model_spatial_functional = "
functions {
  // Matern 3/2 Covariance Function
   matrix matern32_cov(matrix Full_D, real kappa2, real varphi) {
    int M = dims(Full_D)[1];            // Number of sites (rows of D)
    real c = sqrt(3.0) / varphi;
    
    for (i in 1:M) {
      for (j in i:M) {
        if (i == j) {
          // To ensure positive definiteness, we add a small amount of diagonal jitter.
          cov[i, j] = kappa2 + 1e-10;
        } else {
          real h = Full_D[i, j] * c;
          cov[i, j] = kappa2 * (1 + h) * exp(-h);
          cov[j, i] = cov[i, j]; // Ensure symmetry
        }
      }
    }
    return cov;
  }
}

  data {
    int<lower=1> m_obs;      
    int<lower=1> m_new;       
    int<lower=1> m;            
    int<lower=1> n;           
    int<lower=1> A;           
    int<lower=1> G;           
    int<lower=1> R;           
    
    matrix[n, A] Phi0;        
    matrix[n, G] Phi1;        
    matrix[n, G] Phi2;        
    matrix[n, R] Phi_w;       
    
    vector[m_obs] C1_obs;     
    vector[m_obs] C2_obs;     
    vector[m_new] C1_new;     
    vector[m_new] C2_new;     
    matrix[m, m] D;     
    real d_lower;
    real d_upper;
    
    matrix[m_obs, n] Y;                       // simulated data (observed)      
    
  }

    parameters {
    vector[A] q_raw;                           // Global coefficients (raw)
    vector[G] z1_raw;                          // Coefficients for spatial covariate p=1 (raw)
    vector[G] z2_raw;                          // Coefficients for spatial covariate p=2 (raw) 
    matrix[m, R] z_E;                          // Latent variable ~ Normal(0, I_m)
  
    real<lower=1e-6> kappa2;                   // Scale parameter (variance) for the spatial effect
    real<lower=d_lower, upper=d_upper> varphi; // Range parameter for Matern
    real<lower=1e-6> sigma2_q;                 // Variance for q
    real<lower=1e-6> sigma2_z1;                // Variance for z1 
    real<lower=1e-6> sigma2_z2;                // Variance for z2
    real mu_q;                                 // Mean for q
    real mu_z1;                                // Mean for z1 
    real mu_z2;                                // Mean for z2 
    real<lower=1e-6> tau2;                     // Noise variance (nugget)
}
  
    transformed parameters {
    vector[A] q;                                  
    vector[G] z1;                             
    vector[G] z2;
    matrix[m, R] E;
  
    matrix[m, m] Cov_w;                           // Spatial covariance matrix
    Cov_w = matern32_cov(D, kappa2, varphi);
    matrix[m, m] L = cholesky_decompose(Cov_w);   // Cholesky factor of Sigma
    E = L * z_E;

    // Non-centered parameterization reconstruction 
    q  = rep_vector(mu_q,A) + sqrt(sigma2_q) * q_raw;
    z1 = rep_vector(mu_z1,G)+ sqrt(sigma2_z1)* z1_raw;
    z2 = rep_vector(mu_z2,G)+ sqrt(sigma2_z2)* z2_raw;
    E  = L * z_E;

    // Mean function construction for observed locations
    matrix[m_obs, n] mu;
    mu = rep_matrix((Phi0 * q)', m_obs)
     + C1_obs * (Phi1 * z1)'
     + C2_obs * (Phi2 * z2)'
     + E[1:m_obs,] * Phi_w';

    // Mean function construction for new (unobserved) locations  
    matrix[m_new, n] mu_new;
    mu_new = rep_matrix((Phi0 * q)', m_new)
     + C1_new * (Phi1 * z1)'
     + C2_new * (Phi2 * z2)'
     + E[(m_obs+1):m,] * Phi_w';
  
  }
  
    model {
    // ========== Priors ==========
  
    sigma2_q  ~ inv_gamma(2,1.2);         // q  ~ Normal(0, sigma2_q)
    sigma2_z1 ~ inv_gamma(2,0.5);         // z1 ~ Normal(0, sigma2_z1)
    sigma2_z2 ~ inv_gamma(2,0.5);         // z2 ~ Normal(0, sigma2_z2)
    mu_q      ~  normal(0,100);
    mu_z1     ~  normal(0,100);
    mu_z2     ~  normal(0,100);
  
    kappa2    ~ inv_gamma(2,1);
    varphi    ~ uniform(d_lower,d_upper);
    tau2      ~ inv_gamma(2,0.5); 

    q_raw  ~ std_normal();
    z1_raw ~ std_normal();
    z2_raw ~ std_normal();
  
    to_vector(z_E) ~ std_normal();  // Priori “não-centralizada” para z_E

    // ========== Verossimilhança ==========
 
    to_vector(Y) ~ normal(to_vector(mu), sqrt(tau2));
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
    iter_sampling   = 20000,
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


