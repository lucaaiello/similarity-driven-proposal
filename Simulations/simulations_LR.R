rm(list = ls())
library(MASS)

set.seed(1234)

n <- 200
p <- 5  
aux_p <- 495 

rho <- 0.9
Sigma <- outer(1:(p+aux_p), 1:(p+aux_p), function(i, j) rho^abs(i - j))

Xstd <- scale(mvrnorm(n, mu = rep(0, (p+aux_p)), Sigma = Sigma))
X <- cbind(1, Xstd)

active_idx <- sample(1:(p + aux_p), p)
gamma_true <- rep(0, p + aux_p)
gamma_true[active_idx] <- 1

beta <- rnorm(p + 1, 0, 1)
beta
y <- rnorm(n, X[, c(TRUE, gamma_true == 1)] %*% beta, 1)

library(Rcpp)
library(RcppArmadillo)

library(huge)

# Run Graphical Lasso path directly on X
out.huge <- huge(Xstd, method = "glasso")

# Model selection using EBIC
out.select <- huge.select(out.huge, criterion = "ebic")

# Extract adjacency matrix
adj <- out.select$refit

# Get neighbors for each covariate
neighbors <- lapply(1:ncol(Xstd), function(i) which(adj[i, ] == 1))
names(neighbors) <- colnames(Xstd)

gamma_init <- rep(0, dim(X)[2] - 1)
phi_init <- 0.5

# with lambda at at the exponent

sourceCpp('sampler_LR.cpp')

library(mcmcse)

lambda <- seq(0.01, 1.5, length.out = 150)

acc_rate_bd <- vector(length = length(lambda))

for (l in 1:length(lambda)) {
  set.seed(l)
  cat("### lambda = ", lambda[l], "\n")
  cat("\n")
  
  key <- paste0("lambda_", lambda[l])
  
  results <- var_sel_distance_cpp(y, X, neighbors, 
                                  gamma_init, phi_init, 
                                  a = 2, d = 1, tau = 1, a_phi = 1, b_phi = 3,
                                  runs = 20000, burn = 10000,
                                  lambda_bd = lambda[l], adapt = FALSE, 
                                  move_step = FALSE)
  
  acc_rate_bd[l] <- results[["acc_rate_bd"]]
}

acc_rates_df <- data.frame(lambda = lambda,
                           acc_rate_bd = acc_rate_bd)

library(ggplot2)

ggplot(acc_rates_df, aes(x = lambda, y = acc_rate_bd, size = acc_rate_bd, color = acc_rate_bd)) +
  geom_point() +
  scale_x_continuous(breaks = lambda[seq(0, length(lambda), by = 10)]) +  # Exact lambda values
  scale_size_continuous(range = c(1, 5), guide = "none") +  # Adjust range
  scale_color_gradient(low = "blue", high = "red") +
  labs(x = expression(lambda), y = "Acc. rate", color = "Acc. rate") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Move step experiments --------------------------------------------------------

n_runs <- 100
hamming_list <- vector("list", n_runs)

for (r in 1:n_runs) {
  
  set.seed(r)
  
  out <- var_sel_distance_cpp(
    y, X, neighbors,
    gamma_init, phi_init,
    a = 2, d = 1, tau = 1, a_phi = 1, b_phi = 3,
    runs = 20000, burn = 10000,
    lambda_bd = 0.75, lambda_move = 1.25, adapt = TRUE,
    move_step = TRUE
  )
  
  gamma_chain <- out$gamma
  n_iter <- nrow(gamma_chain)
  
  # Compute Hamming distances for this run
  hamming <- numeric(n_iter - 1)
  for (t in 2:n_iter) {
    hamming[t - 1] <- sum(abs(gamma_chain[t, ] - gamma_chain[t - 1, ]))
  }
  
  hamming_list[[r]] <- hamming
}


hamming_all <- unlist(hamming_list)

ggplot(data.frame(hamming_all = hamming_all), aes(x = hamming_all)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 7,
                 fill = "skyblue",
                 color = "black") +
  labs(
    x = "Hamming distance",
    y = "Probability"
  ) +
  theme_minimal()


# Adaptation experiments2 ------------------------------------------------------

set.seed(123)

tic()
out_dis_LR <- var_sel_distance_cpp(y, X, neighbors, 
                                       gamma_init, phi_init, 
                                       a = 1, d = 1, tau = 1, a_phi = 1, b_phi = 1,
                                       runs = 20000, burn = 10000, 
                                       start_adapt = 100, end_adapt = 5000,
                                       lambda_bd = 0.7, lambda_move = 1.25, adapt = TRUE,
                                       move_step = FALSE)
toc()

gamma_chain <- out_dis_LR$gamma

gamma_post_mean <- apply(gamma_chain,2,mean)

plot(gamma_post_mean,type = "h")
points(active_idx, gamma_post_mean[active_idx], col = "red", pch = 16)

plot(1:dim(gamma_chain)[1], apply(gamma_chain, 1, sum))
acf(apply(gamma_chain, 1, sum))
table(apply(gamma_chain, 1, sum))/dim(gamma_chain)[1]

library(ggplot2)

# Plot 1: Posterior means with highlighted active indices
df1 <- data.frame(
  index = 1:length(gamma_post_mean),
  value = gamma_post_mean,
  active = 1:length(gamma_post_mean) %in% active_idx
)

p1 <- ggplot(df1, aes(x = index, y = value)) +
  geom_segment(aes(xend = index, yend = 0), linewidth = 0.33) +
  geom_point(data = subset(df1, active), color = "red", size = 2) +
  labs(x = "Taxa-Covariate index",
       y = "PIPs") +
  theme_minimal()

# Plot 2: Sum of gamma across iterations
gamma_sums <- apply(gamma_chain, 1, sum)
df2 <- data.frame(
  iteration = 1:length(gamma_sums),
  sum = gamma_sums
)

p2 <- ggplot(df2, aes(x = iteration, y = sum)) +
  geom_point() +
  labs(x = "Iteration",
       y = "Number of active variables") +
  scale_y_continuous(breaks = seq(floor(min(gamma_sums)), ceiling(max(gamma_sums)), by = 1)) + 
  theme_minimal()

# Plot 3: ACF of gamma sums
acf_result <- acf(gamma_sums, plot = FALSE)
df3 <- data.frame(
  lag = acf_result$lag,
  acf = acf_result$acf
)

p3 <- ggplot(df3, aes(x = lag, y = acf)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_segment(aes(xend = lag, yend = 0)) +
  geom_hline(yintercept = c(-1.96/sqrt(length(gamma_sums)), 
                            1.96/sqrt(length(gamma_sums))),
             linetype = "dashed", color = "blue") +
  labs(x = "Lag",
       y = "ACF") +
  theme_minimal()

print(p1)
print(p2)
print(p3)

#

set.seed(123)

tic()
out_dis_LR <- var_sel_distance_cpp(y, X, neighbors, 
                                       gamma_init, phi_init, 
                                       a = 1, d = 1, tau = 1, a_phi = 1, b_phi = 1,
                                       runs = 100000, burn = 100, 
                                       start_adapt = 100, end_adapt = 75000,
                                       lambda_bd = 0.7, lambda_move = 1.25, adapt = TRUE,
                                       move_step = FALSE)
toc()

acc_rates_df <- data.frame(lambda = out_dis_LR$lambda_chain[-(1:100),],
                           acc_rate_bd = out_dis_LR$acc_rate_bd_chain[-(1:100),])

library(ggplot2)
library(dplyr)

# Create segments dataframe
segments_df <- acc_rates_df %>%
  mutate(
    x = 1:n(),
    xend = lead(x),
    yend = lead(lambda)
  ) %>%
  filter(!is.na(xend))  # Remove last row (no next point)

ggplot(segments_df, aes(x = x, y = lambda, xend = xend, yend = yend, color = lambda)) +
  geom_segment() +
  scale_color_gradient(low = "blue", high = "red",
                       limits = c(1, max(segments_df$lambda)),
                       oob = scales::squish, guide = "none") +
  labs(x = "Iteration", y = expression(lambda), color = expression(lambda)) +
  theme_minimal()

# Create segments dataframe
segments_df <- acc_rates_df %>%
  mutate(
    x = 1:n(),
    xend = lead(x),
    yend = lead(acc_rate_bd)
  ) %>%
  filter(!is.na(xend))  # Remove last row (no next point)

ggplot(segments_df, aes(x = x, y = acc_rate_bd, xend = xend, yend = yend, color = acc_rate_bd)) +
  geom_segment() +
  scale_color_gradient(low = "blue", high = "red",
                       limits = c(0, max(segments_df$acc_rate_bd)),
                       oob = scales::squish, guide = "none") +
  labs(x = "Iteration", y = "Acceptance rate") +
  theme_minimal()