# Load Required Packages
library("fda")        
library("MASS")     
library("ggplot2")
library("dplyr")
library("tidyr")
library("posterior")
library("bayesplot")
library("coda")
library("cmdstanr")
library("snowfall")

load("../Data/coords.RData")
load("../Data/t.RData")
load("../Data/Yobs.RData")

cpus = 20
int  = 250

n = 100          # The number of observations in each curve.
m_new = 3        # The number of unobserved spatial locations.
m_obs = 22       # The number of observed spatial locations.
m = m_obs+m_new  # Total number of spatial locations.


# Splitting coordinates into observed and unobserved sets
# Note: Indices 6, 18, and 23 correspond to the unobserved locations.
Coord_obs = coords[-c(6, 18, 23),]                                         
Coord_new = coords[c(6, 18, 23),]    

Coord_obs_new = rbind(Coord_obs,Coord_new);rownames(Coord_obs_new)=NULL    # full coordinates.
Full_D = as.matrix(dist(Coord_obs_new))                                    # Matrix of Euclidean distances between all pairs of sites [m × m].

# Standardization of covariates (Mean = 0, Variance = 1) 
coords_scaled = scale(coords)                                             .
coords_scaled_obs = coords_scaled[-c(6, 18, 23),]
coords_scaled_new = coords_scaled[c(6, 18, 23),]

# Extracting Lat/Long covariates
C1_obs = as.numeric(coords_scaled_obs[, "latitude"])   # covariate p=1.
C2_obs = as.numeric(coords_scaled_obs[, "longitude"])  # covariate p=2.

C1_new = as.numeric(coords_scaled_new[, "latitude"])   # covariate p=1.
C2_new = as.numeric(coords_scaled_new[, "longitude"])  # covariate p=2.

# Functional basis
A = 7                              # Number of bases for the interception of β0(t).
G = 7                              # Number of bases for each βₚ(t), where p = 1, 2.
R = 7                              # The number of bases for the spatial effect w_s(t).

d_lower = unname( quantile(Full_D, probs = 0.05) )
d_upper = unname( quantile(Full_D, probs = 0.95) )

###############################################################################
# Construction of B-spline basis functions for β₀(t), β₁(t), β₂(t), and wₛ(t).
###############################################################################

# Interception β0(t)
basis_beta0 = create.bspline.basis(rangeval = c(0,1), nbasis = A, norder = 4)
Phi0 = eval.basis(t, basis_beta0)    # matriz n x G

# Coefficients β1(t) e β2(t)
basis_beta1 = create.bspline.basis(rangeval = c(0,1), nbasis = G, norder = 4)
basis_beta2 = create.bspline.basis(rangeval = c(0,1), nbasis = G, norder = 4)

Phi1 = eval.basis(t, basis_beta1)    # matrix n x R
Phi2 = eval.basis(t, basis_beta2)    # matrix n x R

# Spatial effect w_s(t)
basis_w = create.bspline.basis(rangeval = c(0,1), nbasis = R, norder = 4)
Phi_w = eval.basis(t, basis_w)    

################################################################################
# Stan Model Definition 
# Defined as a character string to be compiled by cmdstanr
################################################################################

model_spatial_functional="

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
    matrix[m, m]  Full_D;     
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
    Cov_w = matern32_cov(Full_D, kappa2, varphi);
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

    generated quantities {
    matrix[m_new, n] Y_new;        // Predictive distribution for new locations
   
    for(i in 1:m_new){
    for(l in 1:n){
    Y_new[i,l]=normal_rng(mu_new[i,l], sqrt(tau2));
   }
 }
}
"
# Compile the model
Model = write_stan_file(model_spatial_functional)
mod = cmdstan_model(Model)

# Monte Carlo Simulation Function 
counter = 0
MonteCarlo = function(step){
# Log file setup - Using current directory for portability
log_file = "./MonteCarlo_log.txt"  
  
  log_msg = function(msg) {
    cat(msg, file = log_file, append = TRUE, sep = "\n")  # Grava a mensagem no arquivo
    flush.console()  # Atualiza o console e o arquivo imediatamente
  }
  
  # Data filtering for the current step
  Y = Yobs %>%
      filter(id==step) %>%
      dplyr::select(dplyr::starts_with("S") & !c("S6", "S18", "S23")) %>%
      as.matrix() %>%
      t()
  
  # List of data to pass to Stan
  data_list <- list(
    m_obs                 = m_obs,
    m_new                 = m_new,
    m                     = m,
    n                     = n,
    A                     = A,
    G                     = G,
    R                     = R,
    d_lower               = d_lower,
    d_upper               = d_upper,
    Phi0                  = Phi0,
    Phi1                  = Phi1,
    Phi2                  = Phi2,
    Phi_w                 = Phi_w,
    C1_obs                = C1_obs,
    C2_obs                = C2_obs,
    C1_new                = C1_new,
    C2_new                = C2_new,
    Full_D                = Full_D,
    Y                     = Y
  )
  
  # Sampling
  fit_model = mod$sample(
    data            = data_list,
    seed            = 123,
    chains          = 2,
    parallel_chains = 2,
    iter_warmup     = 20000,
    iter_sampling   = 20000,
    thin            = 20,
  )

  # --- Diagnostics & HPD for Y_new ---
  draws = fit_model$draws(c("Y_new"))
  variable = as_draws_matrix(draws)
  diagnostic = summarise_draws(variable, "mean","median", "sd", default_convergence_measures())
  id = rep(step,nrow(diagnostic))
  diagnostic = cbind(diagnostic,id)
  
  ###########
  ### Hpd ###
  ###########
  
  Hpd_mcmc = as.mcmc(variable)
  Hpd = HPDinterval(Hpd_mcmc,prob=0.95)
  id2 = rep(step,nrow(Hpd))
  Hpd = cbind(Hpd,id2)
  
  ##############
  ###parameters#
  ##############
  
  draws_2      = fit_model$draws(c("q", "z1", "z2", "E", "kappa2" ,"varphi", "mu_q", "mu_z1", "mu_z2" ,"sigma2_q", "sigma2_z1", "sigma2_z2", "tau2"))
  variable_2   = as_draws_matrix(draws_2)
  diagnostic_2 = summarise_draws(variable_2, "mean","median", "sd", default_convergence_measures())
  id3 = rep(step,nrow(diagnostic_2))
  diagnostic_2 = cbind(diagnostic_2,id3)

  # --- Writing to Files ---
  # Using 'append = TRUE' allows accumulating results
  write.table(diagnostic, file ="./diagnostic.txt",
              sep = ",",
              col.names = FALSE,
              row.names = FALSE ,
              append = TRUE)
  
  write.table(Hpd, file ="./Simu/Hpd.txt",
              sep = ",",
              col.names = FALSE,
              row.names = FALSE ,
              append = TRUE)
  
  write.table(diagnostic_2, file ="./diagnostic_2.txt",
              sep = ",",
              col.names = FALSE,
              row.names = FALSE ,
              append = TRUE)
  
  counter <<- counter + 1
  
  log_msg(paste("Finalizando iteração:", step, "- Total de iterações concluídas:", counter, "\n"))
}
# Parallel Execution
sfInit(parallel=TRUE, cpus=cpus)
sfLibrary(dplyr)
sfLibrary(cmdstanr)
sfLibrary(posterior)
sfLibrary(bayesplot)
sfLibrary(coda)
sfExportAll()
sfLapply(1:int, fun=MonteCarlo) # Function that I want to compute multiple times using sfLapply:
sfStop()


