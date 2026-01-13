# Load required packages
library(mice)
library(MASS)  # For multivariate normal distribution
library(miceDRF)
library(ggplot2)
# Set seed for reproducibility
set.seed(123)

library(drf)


# Simulation parameters
B <- 200       # Number of simulations
n <- 1000        # Sample size per simulation
p_missing <- 0.4  # Probability of setting X1 to NA when X2 > 0
method<-"cart"




get_results <- function(data, method="DRF"){  # Introduce missingness in X1 when X2 > 0 with probability p_missing
  missing_idx <- which(data$X2 > 0)
set_to_na<-missing_idx[ rbinom(n=length(missing_idx), prob=p_missing, size=1)==1 ]
data$X1[set_to_na] <- NA

# Impute missing values using mice with norm.nob method
#imp <- mice(data, method = "norm.nob", m = 5, printFlag = FALSE)
imp <- mice(data, method = method, m = 1, printFlag = FALSE)
#imp <- mice(data, method = "DRF", m = 5, printFlag = FALSE)

# Analyze imputed datasets
fit <- with(imp, lm(X2 ~ X1))
pooled <- pool(fit)

# Extract results
return( list(estimate=summary(pooled)$estimate[2], se= summary(pooled)$std.error[2] ) )
}


newmiceimputation<-function(train_set, test_set, m=1){
  
  imp <- mice(train_set, method = method, m = 1, printFlag = FALSE)
  train_set<-complete(imp)
  
  drf<-list()
  d<-ncol(train_set)
  #Learn method, so far only with DRF!!
  for (j in 1:d){
    drf[[j]] <- drf(Y=train_set[,j, drop=F], X=train_set[,-j, drop=F], num.trees = 10, min.node.size = 1, compute.oob.predictions = F)
    
  }  
  
  imp <- mice(test_set, method = "sample", m = 1, printFlag = FALSE)
  test_setimp<-complete(imp)
  
  cols_with_na <- which(colSums(is.na(test_set)) > 0)
  
  for (s in 1:5){
    
    for (j in cols_with_na){
      indexj<-is.na(test_set[,j])
      xmis<-test_setimp[indexj,-j, drop=F]
      #yobs<-test_setimp[!indexj,j, drop=F]
      yobs<-train_set[,j, drop=F] ##Not so nice, but need observed Ys from other dataset here
      ##Continue here!!
      
      DRFw <- predict(drf[[j]], newdata = xmis)$weights
      test_setimp[indexj,j] <- vapply(1:nrow(xmis), function(s) {
        yobs[sample(1:nrow(yobs), size = 1, replace = T, prob = DRFw[s,]), ]
      }, numeric(1))
      
    }  
    
  }
  
  return(test_setimp)
  
}

get_results_sample_split <- function(data, method="DRF"){  # Introduce missingness in X1 when X2 > 0 with probability p_missing

  d<-ncol(data)
  
  train_indices <- sample(1:nrow(data), 0.5 * nrow(data))
  
  
  train_set <- data[train_indices, ]
  test_set <- data[-train_indices, ]
  
  test_set_imp1<-newmiceimputation(train_set,test_set, m=1)

  estimate1<-unname(coefficients(lm(X2 ~ X1, data=test_set_imp1))[2])
  
  
  #two folds
  
  train_set <- data[-train_indices, ]
  test_set <- data[train_indices, ]
  
  test_set_imp2<-newmiceimputation(train_set,test_set, m=1)
  
  estimate2<-unname(coefficients(lm(X2 ~ X1, data=test_set_imp2))[2])
  
  
  
  estimate<-1/2*(estimate1 + estimate2)
  #estimate<-estimate1
  
  ##Continue here, why so negative??
  
  # Extract results
  return( list(estimate=estimate ) )
}

introducemissing <- function(data){
  
  # Introduce missingness in X1 when X2 > 0 with probability p_missing
  missing_idx <- which(data$X2 > 0)
  set_to_na<-missing_idx[ rbinom(n=length(missing_idx), prob=p_missing, size=1)==1 ]
  data$X1[set_to_na] <- NA
  
  
  return(data)
  
}



L<-20
resultsboot <- data.frame(
  beta_imp = numeric(B),
  se_imp = numeric(B),
  ci_lower = numeric(B),
  ci_upper = numeric(B),
  included = logical(B)
)

# Run simulations
for (b in 1:B) {
  cat(b)
  # Generate bivariate normal data with zero correlation
  mu <- c(0, 0)
  Sigma <- matrix(c(1, 0, 0, 1), ncol = 2)
  data_complete <- as.data.frame(mvrnorm(n, mu = mu, Sigma = Sigma))
  names(data_complete) <- c("X1", "X2")
  
  # Introduce Missingness
  data <- introducemissing(data_complete)
  
  
  res<-matrix(NA, nrow=L, ncol=1)
  for (l in 1:L){
    res[l]<-get_results_sample_split(data = data[sample(1:nrow(data), size=n, replace=T),], method=method)$estimate
  }
  
  resultsboot$beta_imp[b]<- get_results_sample_split(data = data, method=method)$estimate #mean(c(res))
  resultsboot$se_imp[b]<-sqrt(mean( (c(res)-resultsboot$beta_imp[b])^2 ))
  
  # Calculate confidence intervals (95%)
  resultsboot$ci_lower[b] <- resultsboot$beta_imp[b] - 1.96 * resultsboot$se_imp[b]
  resultsboot$ci_upper[b] <- resultsboot$beta_imp[b] + 1.96 * resultsboot$se_imp[b]
  
  # Check if true value (0) is in the confidence interval
  resultsboot$included[b] <- resultsboot$ci_lower[b] <= 0 && 0 <= resultsboot$ci_upper[b]
  
  # Optional: print progress
  if (b %% 100 == 0) cat("Completed", b, "simulations\n")
}

# Summarize results
coverage_probability_Bootstrap <- mean(resultsboot$included)
cat("Coverage probability:", coverage_probability_Bootstrap, "\n")


resultsboot <- resultsboot[order(resultsboot$beta_imp), ]
resultsboot$sim_id <- 1:nrow(resultsboot)

# Create the plot
ggplot(resultsboot, aes(x = sim_id, y = beta_imp)) +
  # Color CIs by whether they cover the true value
  geom_linerange(aes(ymin = ci_lower, ymax = ci_upper, 
                     color = included), 
                 alpha = 0.6, linewidth = 0.5) +
  geom_point(aes(color = included), size = 0.8) +
  # True parameter value
  geom_hline(yintercept = 0, linetype = "dashed", 
             color = "black", linewidth = 0.8) +
  scale_color_manual(values = c("TRUE" = "gray60", "FALSE" = "red"),
                     labels = c("TRUE" = "Covers true value", 
                                "FALSE" = "Misses true value")) +
  labs(title = "Confidence Interval Coverage Bootstrap",
       subtitle = sprintf("Coverage rate: %.1f%% (expected: 95%%)", 
                          mean(resultsboot$included) * 100),
       x = "Simulation run (ordered by estimate)",
       y = expression(hat(beta)),
       color = "") +
  theme_minimal() +
  theme(legend.position = "bottom")





##Rubins Rules

resultsRub <- data.frame(
  beta_imp = numeric(B),
  se_imp = numeric(B),
  ci_lower = numeric(B),
  ci_upper = numeric(B),
  included = logical(B)
)


# Run simulations
for (b in 1:B) {
  # Generate bivariate normal data with zero correlation
  mu <- c(0, 0)
  Sigma <- matrix(c(1, 0, 0, 1), ncol = 2)
  data_complete <- as.data.frame(mvrnorm(n, mu = mu, Sigma = Sigma))
  names(data_complete) <- c("X1", "X2")
  
  # Introduce Missingness
  data <- introducemissing(data_complete)
  
  res<-get_resultsRub(data, method=method)
  resultsRub$beta_imp[b] <-res$estimate
  resultsRub$se_imp[b] <- res$se
  
  # Calculate confidence intervals (95%)
  resultsRub$ci_lower[b] <- resultsRub$beta_imp[b] - 1.96 * resultsRub$se_imp[b]
  resultsRub$ci_upper[b] <- resultsRub$beta_imp[b] + 1.96 * resultsRub$se_imp[b]
  
  # Check if true value (0) is in the confidence interval
  resultsRub$included[b] <- resultsRub$ci_lower[b] <= 0 && 0 <= resultsRub$ci_upper[b]
  
  # Optional: print progress
  if (b %% 100 == 0) cat("Completed", b, "simulations\n")
}

# Summarize results
coverage_probability_RubinsRules <- mean(resultsRub$included)
cat("Coverage probability with Rubins Rules:", coverage_probability_RubinsRules, "\n")

## It seems even for n=1000, there is still an undercoverage with mice-drf!

resultsRub <- resultsRub[order(resultsRub$beta_imp), ]
resultsRub$sim_id <- 1:nrow(resultsRub)

# Create the plot
ggplot(resultsRub, aes(x = sim_id, y = beta_imp)) +
  # Color CIs by whether they cover the true value
  geom_linerange(aes(ymin = ci_lower, ymax = ci_upper, 
                     color = included), 
                 alpha = 0.6, linewidth = 0.5) +
  geom_point(aes(color = included), size = 0.8) +
  # True parameter value
  geom_hline(yintercept = 0, linetype = "dashed", 
             color = "black", linewidth = 0.8) +
  scale_color_manual(values = c("TRUE" = "gray60", "FALSE" = "red"),
                     labels = c("TRUE" = "Covers true value", 
                                "FALSE" = "Misses true value")) +
  labs(title = "Confidence Interval Coverage Bootstrap",
       subtitle = sprintf("Coverage rate: %.1f%% (expected: 95%%)", 
                          mean(resultsRub$included) * 100),
       x = "Simulation run (ordered by estimate)",
       y = expression(hat(beta)),
       color = "") +
  theme_minimal() +
  theme(legend.position = "bottom")



