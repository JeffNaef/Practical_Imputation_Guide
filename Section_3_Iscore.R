
library(mice)
library(miceDRF)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("impute")

library(impute)
library(missForest)
library(miceDRF)


# Simulate from FGM Copula
# Method: Sample U1, then sample U2 | U1 from conditional distribution

# Function to sample U2 given U1 using inverse CDF method
sample_u2_given_u1 <- function(u1, alpha) {
  # The conditional CDF is: F(u2|u1) = u2 + alpha*(2*u1-1)*(u2^2 - u2)
  # We need to invert this to sample u2
  # This is a quadratic equation: alpha*(2*u1-1)*u2^2 + [1 - alpha*(2*u1-1)]*u2 - v = 0
  # where v ~ Uniform(0,1)
  
  v <- runif(1)  # Sample from uniform
  
  a <- alpha * (2*u1 - 1)
  b <- 1 - alpha * (2*u1 - 1)
  c <- -v
  
  # Use quadratic formula: u2 = (-b + sqrt(b^2 - 4ac)) / (2a)
  # We need the root in [0,1]
  
  if (abs(a) < 1e-10) {
    # Linear case (when u1 ≈ 0.5)
    u2 <- v / b
  } else {
    discriminant <- b^2 - 4*a*c
    u2 <- (-b + sqrt(discriminant)) / (2*a)
  }
  
  return(u2)
}

# Simulation function
simulate_fgm <- function(n, alpha) {
  # Step 1: Sample U1 from Uniform(0,1)
  u1 <- runif(n)
  
  # Step 2: Sample U2 | U1 for each u1
  u2 <- sapply(u1, function(x) sample_u2_given_u1(x, alpha))
  
  return(data.frame(U1 = u1, U2 = u2))
}



### Quantile Estimation


n<-5000
d<-5


  set.seed(seeds[1])
  
  # independent uniform
  #X<-matrix(runif(n=d*n), nrow=n, ncol=d)
  # uniform with Gaussian copula
  # X <- gaussian_copula_uniform_sim(n = n, d = 2)$uniform_data
  X<-simulate_fgm(n=n, alpha=1)
  X<-cbind(X,matrix(runif( (d-2)*n ), nrow=n, ncol=d-2 ))
  
  vectors <- matrix(c(
    rep(0, d),
    0, 1, rep(0,d-2),
    1, rep(0,d-1)
  ), nrow = 3, byrow = TRUE)
  
  
  # Generate random draws
  # sample() will generate indices, which we use to select rows from the matrix
  M <- vectors[apply(X,1, function(x) sample(1:3, size = 1, prob=c((x[1]+x[2])/3, (2-x[1])/3, (1-x[2])/3), replace = TRUE)), ]
  
  X.NA<-X
  X.NA[M==1]<-NA
  
  
  colnames(X)<-NULL
  colnames(X)<-paste0("X",1:d)
  colnames(X.NA)<-paste0("X",1:d)
  
  n<-nrow(X)
  
  ################################## imputations #########################################
  ########################################################################################
  
  ## Add your favorite imputations here
  imputations<-list()
  
  #knn
  impute_knn <- function(X) { return(impute.knn(as.matrix(X), colmax=0.99)$data) }
  
  #missForest
  impute_missForest <- function(X) { return(missForest(X)$ximp) }
  
  # mice_cart
  impute_mice_cart <-  miceDRF::create_mice_imputation("cart")

  # mice_rf
  impute_mice_rf <- miceDRF::create_mice_imputation("rf")
  
  
  # mice_drf
  impute_mice_drf <- miceDRF::create_mice_imputation("DRF")
  

  imputation_list<-list(knn=impute_knn,
                        missForest=impute_missForest,
                        mice_cart=impute_mice_cart,
                        mice_rf=impute_mice_rf,
                        mice_drf=impute_mice_drf)
  
  
scores<-Iscores_compare(X.NA, imputation_list, N = 10)
  
scores
  



