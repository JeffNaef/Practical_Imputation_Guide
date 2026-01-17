# Load required packages
library(mice)
library(MASS)  # For multivariate normal distribution
library(miceDRF)
library(missForest)
library(ggplot2)
library(patchwork)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("impute")
library(impute)
# Set seed for reproducibility
set.seed(123)

library(drf)


my_mice_rf <- function(data, m = 5, maxit = 5, ntree = 10, 
                       rfPackage = c("ranger", "randomForest", "literanger"),
                       seed = NULL, saved_forests = NULL, ...) {
  
  # Set seed if provided
  if (!is.null(seed)) set.seed(seed)
  
  # Match rfPackage argument
  rfPackage <- match.arg(rfPackage)
  
  # Find variables with missing values
  vars_to_impute <- names(data)[colSums(is.na(data)) > 0]
  
  # Initialize or use saved forests
  if (is.null(saved_forests)) {
    forests <- vector("list", length(vars_to_impute))
    names(forests) <- vars_to_impute
    fit_forests <- TRUE
  } else {
    forests <- saved_forests
    fit_forests <- FALSE
  }
  
  # Create m imputed datasets
  imputed_list <- vector("list", m)
  
  for (imp in 1:m) {
    
    # Start with original data
    dat <- data
    
    # Iterate maxit times
    for (iter in 1:maxit) {
      
      # Impute each variable
      for (var in vars_to_impute) {
        
        y <- dat[[var]]
        ry <- !is.na(data[[var]])  # Use original missingness pattern
        x <- as.matrix(dat[, setdiff(names(dat), var)])
        
        wy <- !ry
        xobs <- x[ry, , drop = FALSE]
        xmis <- x[wy, , drop = FALSE]
        yobs <- y[ry]
        
        # Either fit new forest or use saved one
        if (fit_forests && iter == maxit && imp == m) {
          # Fit and save the forest
          if (rfPackage == "ranger") {
            forests[[var]] <- ranger::ranger(
              y = yobs,
              x = data.frame(xobs),  # Convert to data.frame
              num.trees = ntree,
              ...
            )
            # Manually store training data
            forests[[var]]$xobs <- xobs
            forests[[var]]$yobs <- yobs
          } else if (rfPackage == "randomForest") {
            data_obs <- data.frame(y = yobs, xobs)
            forests[[var]] <- randomForest::randomForest(
              y ~ .,
              data = data_obs,
              ntree = ntree,
              ...
            )
          }
        }
        
        # Get donors
        if (!is.null(saved_forests) && !is.null(saved_forests[[var]])) {
          # Use saved forest to get donors
          donor_matrix <- get_donors_from_forest(
            forest = saved_forests[[var]],
            xmis = xmis,
            rfPackage = rfPackage
          )
        } else {
          # Use mice internal functions (standard fitting)
          f <- switch(rfPackage, 
                      randomForest = mice:::.randomForest.donors, 
                      ranger = mice:::.ranger.donors, 
                      literanger = mice:::.literanger.donor)
          
          donor_matrix <- f(xobs, xmis, yobs, ntree, ...)
        }
        
        # Sample imputations
        nmis <- sum(wy)
        if (nmis == 1) 
          donor_matrix <- array(donor_matrix, dim = c(1, ntree))
        
        imputations <- apply(donor_matrix, MARGIN = 1, FUN = function(s) {
          sample(unlist(s), 1)
        })
        
        # Fill in the imputations
        dat[[var]][!ry] <- imputations
      }
    }
    
    imputed_list[[imp]] <- dat
  }
  
  return(list(imputed_data = imputed_list, forests = forests))
}

get_donors_from_forest <- function(forest, xmis, rfPackage = c("ranger", "randomForest")) {
  
  rfPackage <- match.arg(rfPackage)
  
  if (rfPackage == "ranger") {
    # Use manually stored training data
    xobs <- forest$xobs
    yobs <- forest$yobs
    
    # IMPORTANT: Ensure column names match
    if (!is.null(colnames(xobs)) && !is.null(colnames(xmis))) {
      if (!all(colnames(xmis) == colnames(xobs))) {
        # Reorder xmis columns to match xobs
        xmis <- xmis[, colnames(xobs), drop = FALSE]
      }
    }
    
    # Predict terminal nodes
    pred_obs <- predict(forest, data = data.frame(xobs), type = "terminalNodes")
    pred_mis <- predict(forest, data = data.frame(xmis), type = "terminalNodes")
    
    nodes_obs <- pred_obs$predictions
    nodes_mis <- pred_mis$predictions
    ntree <- forest$num.trees
    
  } else if (rfPackage == "randomForest") {
    xobs <- forest$x
    yobs <- forest$y
    
    # IMPORTANT: Ensure column names match
    if (!is.null(colnames(xobs)) && !is.null(colnames(xmis))) {
      if (!all(colnames(xmis) == colnames(xobs))) {
        xmis <- xmis[, colnames(xobs), drop = FALSE]
      }
    }
    
    nodes_obs <- attr(predict(forest, newdata = xobs, nodes = TRUE), "nodes")
    nodes_mis <- attr(predict(forest, newdata = xmis, nodes = TRUE), "nodes")
    ntree <- forest$ntree
  }
  
  # Find donors
  nmis <- nrow(xmis)
  donor_matrix <- matrix(list(), nrow = nmis, ncol = ntree)
  
  for (i in 1:nmis) {
    for (j in 1:ntree) {
      in_same_node <- nodes_obs[, j] == nodes_mis[i, j]
      donor_matrix[i, j] <- list(yobs[in_same_node])
    }
  }
  
  return(donor_matrix)
}
# Helper for NULL coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x

get_results_sample_split <- function(data, method="DRF"){  # Introduce missingness in X1 when X2 > 0 with probability p_missing
  
  d<-ncol(data)
  
  train_indices <- sample(1:nrow(data), 0.5 * nrow(data))
  
  
  train_set <- data[train_indices, ]
  test_set <- data[-train_indices, ]
  
  result_train <- my_mice_rf(
    data = train_set,
    m = 1,
    maxit = 5,
    ntree = 10,
    rfPackage = "ranger",
    seed = 123
  )
  
  saved_forests <- result_train$forests
  
  # Step 2: Impute on different data using the SAME forests
  result_test <- my_mice_rf(
    data = test_set,
    m = 1,
    maxit = 5,
    rfPackage = "ranger",
    saved_forests = saved_forests,  # Use pre-fitted forests!
    seed = 456
  )
  
  test_set_imp1<-result_test$imputed_data[[1]]
  
  estimate1<- quantile(test_set_imp1[,1], probs=0.1) #unname(coefficients(lm(X2 ~ X1, data=test_set_imp1))[2])
  
  
  #two folds
  
  train_set <- data[-train_indices, ]
  test_set <- data[train_indices, ]
  
  result_train <- my_mice_rf(
    data = train_set,
    m = 1,
    maxit = 5,
    ntree = 10,
    rfPackage = "ranger",
    seed = 123
  )
  
  saved_forests <- result_train$forests
  
  # Step 2: Impute on different data using the SAME forests
  result_test <- my_mice_rf(
    data = test_set,
    m = 1,
    maxit = 5,
    rfPackage = "ranger",
    saved_forests = saved_forests,  # Use pre-fitted forests!
    seed = 456
  )
  
  test_set_imp2<-result_test$imputed_data[[1]]
  
  estimate2<-quantile(test_set_imp2[,1], probs=0.1) #unname(coefficients(lm(X2 ~ X1, data=test_set_imp2))[2])
  
  
  
  estimate<-1/2*(estimate1 + estimate2)
  #estimate<-estimate1
  
  
  
  # Extract results
  return( list(estimate=estimate ) )
}



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



### Define Imputation Functions
#knn
impute_knn <- function(X) { return(impute.knn(as.matrix(X), colmax=0.99)$data) }

#missForest
impute_missForest <- function(X) { return(missForest(X)$ximp) }

# mice_cart
impute_mice_cart <-  miceDRF::create_mice_imputation("cart")

# mice_rf
impute_mice_rf_m <- function(data, m=5){
  dat_imp_norm <- mice::mice(data, m = m, method = "rf", printFlag = FALSE)

  
  # Method 5: More explicit loop
  imputed_list <- vector("list", m)
  res<-c()
  for (i in 1:m) {
    imputed_list[[i]] <- complete(dat_imp_norm, i)
    
    res[i]<-quantile(imputed_list[[i]][,1], probs=0.1)
    
  }
  
  return(mean(res))
  
  
}

impute_mice_rf <- miceDRF::create_mice_imputation("rf")


# mice_drf
impute_mice_drf <- miceDRF::create_mice_imputation("DRF")


#methods<-c("knn", "missForest", "mice_cart", "mice_rf","mice_drf") 

methods<-c("mice_cart","sample_split_rf")

n<-2000
d<-5
alpha<-0.1
B<-200


L<-20
resultsboot<-list()
for (method in methods){
  
  resultsboot[[method]] <- data.frame(
    beta_imp = numeric(B),
    se_imp = numeric(B),
    ci_lower = numeric(B),
    ci_upper = numeric(B),
    included = logical(B)
  )
  
}  


# Run simulations
for (b in 1:B) {
  cat(b)
  
  
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
  
  #Bootstrapping L times
  res<-list()
  for (l in 1:L){
    ## Each time impute + calculate quantile!
    Xl<-X.NA[sample(1:nrow(X.NA), size=n, replace=T),]
    #res[["knn"]][l]<-quantile(impute_knn(Xl)[,1], probs=0.1)
    #res[["missForest"]][l]<-quantile(impute_missForest(Xl)[,1], probs=0.1)
    res[["mice_cart"]][l]<-quantile(impute_mice_cart(Xl)[,1], probs=0.1)
    #res[["mice_rf"]][l]<-quantile(impute_mice_rf(Xl)[,1], probs=0.1)
    res[["sample_split_rf"]][l]<-unlist(get_results_sample_split(Xl))
    #res[["multiple_rf"]][l]<-unlist(impute_mice_rf_m(Xl))
    #res[["mice_drf"]][l]<-quantile(impute_mice_drf(Xl)[,1], probs=0.1)
  }
  
  res0<-list()
  #res0[["knn"]][l]<-quantile(impute_knn(Xl)[,1], probs=0.1)
  #res0[["missForest"]][l]<-quantile(impute_missForest(Xl)[,1], probs=0.1)
  res0[["mice_cart"]]<-quantile(impute_mice_cart(X.NA)[,1], probs=0.1)
  #res0[["mice_rf"]][l]<-quantile(impute_mice_rf(Xl)[,1], probs=0.1)
  res0[["sample_split_rf"]]<-unlist(get_results_sample_split(X.NA))
  #res0[["multiple_rf"]][l]<-unlist(impute_mice_rf_m(Xl))
  #res0[["mice_drf"]][l]<-quantile(impute_mice_drf(Xl)[,1], probs=0.1)
  
  
  for (method in methods){
    
    
    resultsboot[[method]]$beta_imp[b]<- res0[[method]]
    resultsboot[[method]]$se_imp[b]<-sqrt(mean( (c(res[[method]])-resultsboot[[method]]$beta_imp[b])^2 ))
    
    # Calculate confidence intervals (95%)
    resultsboot[[method]]$ci_lower[b] <- resultsboot[[method]]$beta_imp[b] - 1.96 * resultsboot[[method]]$se_imp[b] # 2*resultsboot[[method]]$beta_imp[b] - quantile(res[[method]], probs=1-0.05/2)  #resultsboot[[method]]$beta_imp[b] - 1.96 * resultsboot[[method]]$se_imp[b]
    resultsboot[[method]]$ci_upper[b] <- resultsboot[[method]]$beta_imp[b] + 1.96 * resultsboot[[method]]$se_imp[b]#2*resultsboot[[method]]$beta_imp[b] - quantile(res[[method]], probs=0.05/2)   #resultsboot[[method]]$beta_imp[b] + 1.96 * resultsboot[[method]]$se_imp[b]
    
    # Check if true value (0) is in the confidence interval
    resultsboot[[method]]$included[b] <- resultsboot[[method]]$ci_lower[b] <= alpha && alpha <= resultsboot[[method]]$ci_upper[b]
  }
  
  
  # Optional: print progress
  if (b %% 10 == 0) 
  {
    cat("Completed", b, "simulations\n")
    
    # Create list to store plots
    plot_list <- list()
    
    for (method in methods){
      # Summarize results
      
      mat<-resultsboot[[method]][1:b,]
      
      coverage_probability_Bootstrap <- mean(mat$included)
      cat(sprintf("%s - Coverage probability: %.3f\n", method, coverage_probability_Bootstrap))
      
      # Order and prepare data
      mat$sim_id <- 1:b
      
      # Create the plot
      plot_list[[method]] <- ggplot(mat, aes(x = sim_id, y = beta_imp)) +
        geom_linerange(aes(ymin = ci_lower, ymax = ci_upper, color = included), 
                       alpha = 0.6, linewidth = 0.5) +
        geom_point(aes(color = included), size = 0.8) +
        geom_hline(yintercept = alpha, linetype = "dashed", 
                   color = "black", linewidth = 0.8) +
        scale_color_manual(values = c("TRUE" = "gray60", "FALSE" = "red"),
                           labels = c("TRUE" = "Covers", "FALSE" = "Misses"),
                           guide = "none") +  # Remove individual legends
        labs(title = method,
             subtitle = sprintf("Coverage: %.1f%%", 
                               mat$included) * 100),
             x = "Simulation (ordered by estimate)",
             y = expression(hat(beta))) +
        theme_minimal() +
        theme(plot.title = element_text(size = 11, face = "bold"),
              plot.subtitle = element_text(size = 9))
    }
    
    # Combine plots in 2x2 grid
    combined_plot <- wrap_plots(plot_list, ncol = 2) +
      plot_annotation(
        title = "Confidence Interval Coverage: Bootstrap Methods",
        subtitle = "Expected coverage: 95%",
        theme = theme(plot.title = element_text(size = 14, face = "bold"))
      )
    
    print(combined_plot)
    
  }
}









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



