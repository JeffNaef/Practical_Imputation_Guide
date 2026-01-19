# Load required packages
library(mice)
library(MASS)  # For multivariate normal distribution
library(ggplot2)
library(patchwork)



# Set seed for reproducibility
set.seed(123)





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





n<-500
d<-5
alpha<-0.1
B<-200

##Rubins Rules

resultsRub <- data.frame(
  beta_imp = numeric(B),
  se_imp = numeric(B),
  ci_lower = numeric(B),
  ci_upper = numeric(B),
  included = logical(B)
)

# Uniform quantile SE (no density estimation needed!)
quantile_se_uniform <- function(x, prob) {
  n <- length(x)
  range_x <- diff(range(x))
  sqrt(prob * (1 - prob) / n) * range_x
}

# Simple wrapper for pool()
quantile_fit <- function(X1, prob = 0.1) {
  q <- quantile(X1, probs = prob)
  se <- quantile_se_uniform(X1, prob = prob)
  
  result <- list(
    coefficients = q,
    vcov = matrix(se^2)
  )
  class(result) <- "quantile_fit"
  return(result)
}

# Add tidy method for the quantile_fit class
tidy.quantile_fit <- function(x, ...) {
  data.frame(
    term = "quantile",
    estimate = x$coefficients,
    std.error = sqrt(x$vcov[1,1])
  )
}


# Run simulations
for (b in 1:B) {
  
  
  ## Generate data ##
  
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

  # Analyze imputed datasets
  imp <- mice(X.NA, method = "cart", m = 30, printFlag = FALSE)
  fit <- with(imp, quantile_fit(X1, prob = 0.1))
  pooled <- pool(fit)
  
  ##Manual Implementation##
  # Extract results from each imputation
  # m <- 30
  # results <- data.frame(
  #   Q = sapply(fit$analyses, function(x) x$coefficients),
  #   SE = sapply(fit$analyses, function(x) sqrt(x$vcov[1,1]))
  # )
  # 
  # # Apply Rubin's rules manually
  # Q_bar <- mean(results$Q)           # Pooled estimate
  # W <- mean(results$SE^2)            # Within-imputation variance
  # B <- var(results$Q)                # Between-imputation variance
  # T <- W + (1 + 1/m) * B             # Total variance
  # SE_pooled <- sqrt(T)               # Pooled SE
  #####################
  
  
  tmp<-summary(pooled)
  
  resultsRub$beta_imp[b] <-tmp$estimate
  resultsRub$se_imp[b] <- tmp$std.error
  
  # Calculate confidence intervals (95%)
  resultsRub$ci_lower[b] <- resultsRub$beta_imp[b] - 1.96 * resultsRub$se_imp[b]
  resultsRub$ci_upper[b] <- resultsRub$beta_imp[b] + 1.96 * resultsRub$se_imp[b]
  
  # Check if true value (0) is in the confidence interval
  resultsRub$included[b] <- resultsRub$ci_lower[b] <= alpha && alpha <= resultsRub$ci_upper[b]
  
  
  # Optional: print progress
  if (b %% 10 == 0) 
  {
    cat("Completed", b, "simulations\n")
    
    
      # Summarize results
      
      mat<-resultsRub[1:b,]
      
      coverage_probability_Bootstrap <- mean(mat$included)
      cat(sprintf("%s - Coverage probability: %.3f\n", "cart", coverage_probability_Bootstrap))
      
      # Order and prepare data
      mat$sim_id <- 1:b
      
      # Create the plot
      ggplot(mat, aes(x = sim_id, y = beta_imp)) +
        geom_linerange(aes(ymin = ci_lower, ymax = ci_upper, color = included), 
                       alpha = 0.6, linewidth = 0.5) +
        geom_point(aes(color = included), size = 0.8) +
        geom_hline(yintercept = alpha, linetype = "dashed", 
                   color = "black", linewidth = 0.8) +
        scale_color_manual(values = c("TRUE" = "gray60", "FALSE" = "red"),
                           labels = c("TRUE" = "Covers", "FALSE" = "Misses"),
                           guide = "none") +  # Remove individual legends
        labs(title = "cart",
             subtitle = sprintf("Coverage: %.1f%%", 
                                mean(mat$included) * 100),
             x = "Simulation (ordered by estimate)",
             y = expression(hat(q))) +
        theme_minimal() +
        theme(plot.title = element_text(size = 11, face = "bold"),
              plot.subtitle = element_text(size = 9))
    }
  
    
  
}

save(resultsRub, file="RubinResults_n500.Rdata")


png(filename = "Rubin_CIs_n500.png", 
    width = 1700,    # Width in pixels
    height = 800,    # Height in pixels
    res = 120)       # Resolution in dpi


par(mfrow=c(1,1))



  # Summarize results
  
  mat<-resultsRub
  
  coverage_probability_Bootstrap <- mean(mat$included)
  cat(sprintf("%s - Coverage probability: %.3f\n", "cart", coverage_probability_Bootstrap))
  
  # Order and prepare data
  mat$sim_id <- 1:b
  
  # Create the plot
  ggplot(mat, aes(x = sim_id, y = beta_imp)) +
    geom_linerange(aes(ymin = ci_lower, ymax = ci_upper, color = included), 
                   alpha = 0.6, linewidth = 0.5) +
    geom_point(aes(color = included), size = 0.8) +
    geom_hline(yintercept = alpha, linetype = "dashed", 
               color = "black", linewidth = 0.8) +
    scale_color_manual(values = c("TRUE" = "gray60", "FALSE" = "red"),
                       labels = c("TRUE" = "Covers", "FALSE" = "Misses"),
                       guide = "none") +  # Remove individual legends
    labs(title = "cart",
         subtitle = sprintf("Coverage: %.1f%%", 
                            mean(mat$included) * 100),
         x = "Simulation (ordered by estimate)",
         y = expression(hat(q))) +
    theme_minimal() +
    theme(plot.title = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 9))




# Close the PNG device
dev.off()

