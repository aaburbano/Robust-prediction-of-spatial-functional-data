rm(list=ls(all=TRUE))
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
library("sf")

setwd("")

load("coord_full.RData")
load("coordpred.RData")
load("TMP2023_Day.RData")

Coord_obs = coord_full %>%
            select(Lon, Lat)

Coord_new = coordpred[1:2,] %>%
            select(Lon, Lat)

days = 1:365

#-------------------------------------------------------------------------------
# 1. Convert the data frame to a spatial 'sf' object.
# - 'coords' specifies which columns are the coordinates (longitude, latitude).
# - 'crs' defines the original Coordinate Reference System (CRS).
#   EPSG:4326 is the code for the WGS 84 system (most common for Lon/Lat).
obs_sf = st_as_sf(Coord_obs, coords = c("Lon", "Lat"), crs = 4326)

# 2. Transform (reproject) coordinates to the UTM system.
# - The EPSG for WGS 84 / UTM Zone 14N is 32614.
#   The st_transform() function handles all the projection mathematics.
obs_utm = st_transform(obs_sf, crs = 32614)

# 3. Extract UTM coordinates to a new data frame (optional, but useful).
# The st_coordinates() function extracts the X (Easting) and Y (Northing) values
# from the geometry column of the 'sf' object.
coordinates_obs_utm = as.data.frame(st_coordinates(obs_utm))
names(coordinates_obs_utm) = c("Easting", "Northing") # Renaming columns

coordinates_obs_utm_km = coordinates_obs_utm
coordinates_obs_utm_km$Easting  = coordinates_obs_utm$Easting/1000
coordinates_obs_utm_km$Northing = coordinates_obs_utm$Northing/1000

##
##

pred_sf = st_as_sf(Coord_new, coords = c("Lon", "Lat"), crs = 4326)
pred_utm = st_transform(pred_sf, crs = 32614)

coordinates_pred_utm = as.data.frame(st_coordinates(pred_utm))
names(coordinates_pred_utm) = c("Easting", "Northing") # Renaming columns

coordinates_pred_utm_km = coordinates_pred_utm
coordinates_pred_utm_km$Easting  = coordinates_pred_utm$Easting/1000
coordinates_pred_utm_km$Northing = coordinates_pred_utm$Northing/1000

coordinatesfull_km = rbind(coordinates_obs_utm_km, coordinates_pred_utm_km)
rownames(coordinatesfull_km)=NULL    # Complete coordinates


#-------------------------------------------------------------------------------


n=365
m_new=2
m_obs=23
m=m_obs+m_new      # Number of spatial locations

Full_D = as.matrix(dist(coordinatesfull_km))


coords_scaled_obs = scale(coordinates_obs_utm_km)                   # Optional: standardizes covariates to mean 0, variance 1
coords_scaled_new = scale(coordinates_pred_utm_km)

C1_obs = as.numeric(coords_scaled_obs[, "Easting"])  # Covariate p=1
C2_obs = as.numeric(coords_scaled_obs[, "Northing"]) # Covariate p=2

C1_new = as.numeric(coords_scaled_new[, "Easting"])  # Covariate p=1
C2_new = as.numeric(coords_scaled_new[, "Northing"]) # Covariate p=2

# Functional bases
A = 147                              # Number of bases for intercept β0(t)
G = 147                              # Number of bases for each βp(t), p= 1,2
R = 147                              # Number of bases for spatial effect w_s(t)

d_lower = unname( quantile(Full_D, probs = 0.05) )
d_upper = unname( quantile(Full_D, probs = 0.95) )

###############################################################################
# Construction of B-spline bases for β0(t), β1(t), β2(t) and w_s(t)
###############################################################################

# 3.1. Intercept β0(t)
basis_beta0 = create.bspline.basis(rangeval = c(0,365), nbasis = A, norder = 4)
Phi0 = eval.basis(days, basis_beta0)    # n x A matrix

# 3.2. Coefficients β1(t) and β2(t)
basis_beta1 = create.bspline.basis(rangeval = c(0,365), nbasis = G, norder = 4)
basis_beta2 = create.bspline.basis(rangeval = c(0,365), nbasis = G, norder = 4)

Phi1 = eval.basis(days, basis_beta1)    # n x G matrix
Phi2 = eval.basis(days, basis_beta2)    # n x G matrix

# 3.3. Spatial effect w_s(t)
basis_w = create.bspline.basis(rangeval = c(0,365), nbasis = R, norder = 4)
Phi_w = eval.basis(days, basis_w)       # n x R matrix 

Y = TMP2023_Day %>%
    select(- Day) %>%
    as.matrix() %>%
    t()

#-------------------------------------------------------------------------------
model_spatial_functional="

functions {
  // Matérn 3/2 Function
   matrix matern32_cov(matrix Full_D, real kappa2, real varphi) {
    int M = dims(Full_D)[1];            // Number of sites (rows of D)
    matrix[M, M] cov;
    real c = sqrt(3.0) / varphi;
    
    for (i in 1:M) {
      for (j in i:M) {
        if (i == j) {
          // Diagonal: add a jitter to ensure positive-definiteness
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
    int<lower=1> m_obs;       // Number of observed locations (stations)
    int<lower=1> m_new;       // Number of new locations
    int<lower=1> m;           // Total number of locations 
    int<lower=1> n;           // Number of time points 
    int<lower=1> A;           // Number of B-spline bases for term q 
    int<lower=1> G;           // Number of bases for beta1(t) and beta2(t)
    int<lower=1> R;           // Number of B-spline bases for the spatial term
    
    matrix[n, A] Phi0;        // Intercept basis matrix
    matrix[n, G] Phi1;        // Beta 1 basis matrix
    matrix[n, G] Phi2;        // Beta 2 basis matrix 
    matrix[n, R] Phi_w;       // Temporal basis matrix for spatial effect
    
    vector[m_obs] C1_obs;     // First standardized covariate (e.g., longitude)
    vector[m_obs] C2_obs;     // Second standardized covariate (e.g., latitude)
    vector[m_new] C1_new;     // First standardized covariate (e.g., longitude)
    vector[m_new] C2_new;     // Second standardized covariate (e.g., latitude)
    matrix[m, m]  Full_D;     // Matrix of (Euclidean) distances between all pairs of sites [m × m]
    real d_lower;
    real d_upper;
    
    matrix[m_obs, n] Y;       // Observed data
    
  }

   parameters {
    vector[A] q_raw;                                 // Global coefficients
    vector[G] z1_raw;                                // Coefficients for spatial covariate p = 1
    vector[G] z2_raw;                                // Coefficients for spatial covariate p = 2 
    matrix[m, R] z_E;                                // Each column ~ Normal(0, I_m)
  
  real<lower=1e-6> kappa2;                           // Scale (variance) for the spatial effect
  real<lower= d_lower, upper=d_upper> varphi;        // Matérn range
  real<lower=1e-6> sigma2_q;                         // Variance for q
  real<lower=1e-6> sigma2_z1;                        // Variance for z1 
  real<lower=1e-6> sigma2_z2;                        // Variance for z2
  real mu_q;                                         // Mean of q
  real mu_z1;                                        // Mean for z1 
  real mu_z2;                                        // Mean for z2 
  real<lower=1e-6> tau2;                             // Noise variance
  }
  
  transformed parameters {
  vector[A] q;                                 
  vector[G] z1;                              
  vector[G] z2;
  matrix[m, R] E;
  
  matrix[m, m] Cov_w;                          // Spatial covariance matrix
  Cov_w = matern32_cov(Full_D, kappa2, varphi);
  matrix[m, m] L = cholesky_decompose(Cov_w);  // Cholesky factor of Sigma
  E = L * z_E;
  
 q  = rep_vector(mu_q,A) + sqrt(sigma2_q) * q_raw;
 z1 = rep_vector(mu_z1,G)+ sqrt(sigma2_z1)* z1_raw;
 z2 = rep_vector(mu_z2,G)+ sqrt(sigma2_z2)* z2_raw;
 E  = L * z_E;
  
  matrix[m_obs, n] mu;
  mu = rep_matrix((Phi0 * q)', m_obs)
     + C1_obs * (Phi1 * z1)'
     + C2_obs * (Phi2 * z2)'
     + E[1:m_obs,] * Phi_w';
     
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
  
  to_vector(z_E) ~ std_normal();  // Non-centered prior for z_E

  // ========== Likelihood ==========
  
  to_vector(Y) ~ normal(to_vector(mu), sqrt(tau2));
  }

generated quantities {
   matrix[m_new, n] Y_new;        // Predictive for new locations
   matrix[m_obs, n] log_lik;      // Matrix for pointwise log-likelihood
   
  for(i in 1:m_new){
    for(l in 1:n){
    Y_new[i,l]=normal_rng(mu_new[i,l], sqrt(tau2));
   }
  }
 
  for(i in 1:m_obs){
    for(j in 1:n){
      log_lik[i,j] = normal_lpdf(Y[i,j] | mu[i,j], sqrt(tau2));
    }
  }
}
"

Model = write_stan_file(model_spatial_functional)
mod   = cmdstan_model(Model)

data_list = list(
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

# Compile and run the Stan model

fit_model = mod$sample(
  data            = data_list,
  seed            = 123,
  chains          = 2,
  parallel_chains = 2,
  iter_warmup     = 20000,
  iter_sampling   = 20000,
  thin            = 20
)
