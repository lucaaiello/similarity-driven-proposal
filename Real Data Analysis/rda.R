rm(list = ls())

library(readr)

# combo_count_tab: contains COMBO OTU Counts at all the five levels.
df <- as.matrix(read.table("combo_count_tab.txt"))
lens <- sapply(strsplit(colnames(df), "\\."), length)

# Get the genera - We select genera with prevalence > 5%
Y <- df[, (lens == 6) & apply(df, 2, function(x) sum(x!=0)) >= round(nrow(df) * 0.25)] 

# udiet.adj.scale: FFQ nutrient data, standardized.  The nutrient intakes are
# highly correlated so we group them and select one as the representative of
# the group.

X <- as.matrix(read.table("udiet.adj.scale.txt"))
hclust.tre <- hclust(as.dist(1-cor(X)), method="average")
membership <- cutree(hclust.tre, h=0.1)
ind <- NULL
for (i in 1:max(membership)) {
  ind.i <- which(membership == i)
  ind <- c(ind, ind.i[1])
}
X <- X[, ind]

J <- dim(Y)[2]
P <- dim(X)[2]

labels <- read.table("ffq.labels2.txt")

library(Rcpp)
library(RcppArmadillo)
library(tictoc)

sourceCpp("DM_RJMCMC_VS_sampler.cpp")

set.seed(1234)

tic()

results_inf <- DM_informed_VS_RJMCMC(
  X, Y,
  s2 = rep(10, J), r2 = rep(10, J),
  a = 1, b = 9,
  g_initial = 1, c = 1,
  runs = 10000, burn = 5000,
  informed = TRUE,
  adapt_g = TRUE,
  adaptation_start = 100,
  adaptation_end = 5000,
  neigh = neighbors,
  use_move = FALSE,
  n_move_attempts = 1,
  lambda_move = 1.0
)

toc()

round(colMeans(results_inf$alpha), 2)

xi_post_mean <- apply(results_inf$xi, c(1, 2), mean)

plot(as.vector(xi_post_mean),type = "h")
sum(as.vector(xi_post_mean)> 0.5)

plot(apply(results_inf$xi, 3, sum))
acf(apply(results_inf$xi, 3, sum))

# FDR ---------------------------------------------------------------------

threshold <- seq(0.01, 1, by = 0.01)
C <- length(threshold)
FDR <- vector(length = C)

for (c in 1:C) {
  num <- 0
  den <- 0
  for (p in 1:P) {
    for (j in 1:J) {
      D_pj <- (xi_post_mean[p,j] > threshold[c])
      num <- num + (1 - xi_post_mean[p,j]) * D_pj 
      den <- den + D_pj
    }
  }
  FDR[c] <- num/den
}

threshold[min(which(FDR < 0.05))]

# Cool plots --------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)

# Posterior mean of xi over the 3rd dimension
xi_post_mean <- apply(results_inf$xi, c(1, 2), mean)

# Flatten for plotting
df_xi_mean <- data.frame(
  index = 1:length(as.vector(xi_post_mean)),
  value = as.vector(xi_post_mean)
)

# Sum of xi over first two dims for each iteration
xi_sum_iter <- apply(results_inf$xi, 3, sum)
df_xi_sum <- data.frame(
  iteration = 1:length(xi_sum_iter),
  value = xi_sum_iter
)

# ACF data
acf_xi <- acf(xi_sum_iter, plot = FALSE)
df_acf <- data.frame(
  lag = acf_xi$lag,
  acf = acf_xi$acf
)

ggplot(df_xi_mean, aes(x = index, y = value)) +
  geom_linerange(aes(ymin = 0, ymax = value), color = "gray1") +
  geom_point(data = df_xi_mean %>% filter(value > threshold[min(which(FDR < 0.05))]),
             aes(x = index, y = value),
             color = "red", size = 2) +
  scale_x_continuous(breaks = seq(0, 3500, by = 500)) +
  labs(x = "Taxa-covariate parameter index", y = "PIPs") +
  theme_minimal()

ggplot(df_xi_sum, aes(x = iteration, y = value)) +
  geom_point(color = "gray1") +
  labs(x = "Iteration", y = "Number of variables included") +
  theme_minimal()

ggplot(df_acf, aes(x = lag, y = acf)) +
  geom_col(fill = "gray1", width = 0.1) +
  geom_hline(yintercept = 0, color = "black") +
  labs(x = "Lag", y = "ACF") +
  theme_minimal()

which(xi_post_mean > threshold[min(which(FDR < 0.05))], arr.ind = TRUE)

# For illustration:
indices <- which(xi_post_mean > threshold[min(which(FDR < 0.05))], arr.ind = TRUE)

# Get covariate names
cov_names <- colnames(X)[indices[, "row"]]

# Get category names
category_names <- colnames(Y)[indices[, "col"]]#as.character((1:50)[indices[, "col"]])

# Build a data frame
df_assoc <- data.frame(
  Genera = category_names,
  Nutrients = cov_names
)

# Summarise in a table: list all covariates per category
library(dplyr)

summary_table <- df_assoc %>%
  group_by(Genera) %>%
  summarise(Nutrients = paste(Nutrients, collapse = ", "), .groups = "drop")

summary_table

# Replace acronyms with full names
df_assoc$NutrientsFull <- labels[df_assoc$Nutrients, 1]

# Build the summary table with full names
library(dplyr)

summary_table_full <- df_assoc %>%
  group_by(Genera) %>%
  summarise(Nutrients = paste(NutrientsFull, collapse = ", "), .groups = "drop")

summary_table_full

################################################################################
library(dplyr)
library(tidyr)

summary_table_full_sep <- summary_table_full %>%
  separate(
    Genera,
    into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
    sep = "\\.",
    fill = "right"
  ) %>%
  mutate(Species = NA) %>%  # Since species are not indicated
  select(Kingdom, Phylum, Class, Order, Family, Genus, Species, Nutrients)

# View the result
summary_table_full_sep

library(dplyr)
library(tidyr)

summary_table_long <- summary_table_full_sep %>%
  separate_rows(Nutrients, sep = ",\\s*")  # Split on comma and optional space


p1 <- ggplot(df_xi_mean, aes(x = index, y = value)) +
  geom_linerange(aes(ymin = 0, ymax = value), color = "gray1") +
  geom_point(data = df_xi_mean %>% filter(value > threshold[min(which(FDR < 0.05))]),
             aes(x = index, y = value),
             color = "red", size = 2) +
  scale_x_continuous(breaks = seq(0, 3500, by = 500)) +
  labs(x = "Taxa-covariate parameter index", y = "PIPs") +
  theme_minimal()

p2 <- ggplot(df_xi_sum, aes(x = iteration, y = value)) +
  geom_point(color = "gray1") +
  labs(x = "Iteration", y = "Number of variables included") +
  theme_minimal()

p3 <- ggplot(df_acf, aes(x = lag, y = acf)) +
  geom_col(fill = "gray1", width = 0.1) +
  geom_hline(yintercept = 0, color = "black") +
  labs(x = "Lag", y = "ACF") +
  theme_minimal()

print(p1)
print(p2)
print(p3)

#################################################################################

library(huge)

# Estimate a predictor graph and turn it into the neighbor list required
# by the local move step.
out.huge <- huge(X, method = "glasso")
out.select <- huge.select(out.huge, criterion = "ebic", ebic.gamma = 0.25)
adj <- out.select$refit
neighbors <- lapply(seq_len(ncol(X)), function(i) which(adj[i, ] == 1))

# If the selected graph is still empty, fall back to correlation-based
# nearest neighbors so the move step remains available.
if (all(lengths(neighbors) == 0)) {
  message("EBIC-selected graph is empty; falling back to correlation-based neighbors.")
  cor_abs <- abs(cor(X))
  diag(cor_abs) <- 0
  k_fallback <- min(5L, ncol(X) - 1L)
  neighbors <- lapply(seq_len(ncol(X)), function(i) {
    order(cor_abs[i, ], decreasing = TRUE)[seq_len(k_fallback)]
  })
}

message(
  "Neighbor counts min/median/max: ",
  min(lengths(neighbors)), "/",
  median(lengths(neighbors)), "/",
  max(lengths(neighbors))
)

names(neighbors) <- if (is.null(colnames(X))) {
  paste0("X", seq_len(ncol(X)))
} else {
  colnames(X)
}

sourceCpp("DM_RJMCMC_VS_sampler.cpp")

set.seed(1234)

results_inf_move <- DM_informed_VS_RJMCMC(
  X, Y,
  s2 = rep(10, J), r2 = rep(10, J),
  a = 1, b = 9,
  g_initial = 1, c = 1,
  runs = 10000, burn = 5000,
  informed = TRUE,
  adapt_g = TRUE,
  adaptation_start = 100,
  adaptation_end = 5000,
  neigh = neighbors,
  use_move = TRUE,
  n_move_attempts = 1,
  lambda_move = 1.0
)

round(colMeans(results_inf_move$alpha), 2)

xi_post_mean <- apply(results_inf_move$xi, c(1, 2), mean)

plot(as.vector(xi_post_mean),type = "h")
sum(as.vector(xi_post_mean)> 0.5)

plot(apply(results_inf_move$xi, 3, sum))
acf(apply(results_inf_move$xi, 3, sum))

# FDR ---------------------------------------------------------------------

threshold <- seq(0.01, 1, by = 0.01)
C <- length(threshold)
FDR <- vector(length = C)

for (c in 1:C) {
  num <- 0
  den <- 0
  for (p in 1:P) {
    for (j in 1:J) {
      D_pj <- (xi_post_mean[p,j] > threshold[c])
      num <- num + (1 - xi_post_mean[p,j]) * D_pj 
      den <- den + D_pj
    }
  }
  FDR[c] <- num/den
}

threshold[min(which(FDR < 0.05))]

# Cool plots --------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(tidyr)

# Posterior mean of xi over the 3rd dimension
xi_post_mean <- apply(results_inf_move$xi, c(1, 2), mean)

# Flatten for plotting
df_xi_mean <- data.frame(
  index = 1:length(as.vector(xi_post_mean)),
  value = as.vector(xi_post_mean)
)

# Sum of xi over first two dims for each iteration
xi_sum_iter <- apply(results_inf_move$xi, 3, sum)
df_xi_sum <- data.frame(
  iteration = 1:length(xi_sum_iter),
  value = xi_sum_iter
)

# ACF data
acf_xi <- acf(xi_sum_iter, plot = FALSE)
df_acf <- data.frame(
  lag = acf_xi$lag,
  acf = acf_xi$acf
)

ggplot(df_xi_mean, aes(x = index, y = value)) +
  geom_linerange(aes(ymin = 0, ymax = value), color = "gray1") +
  geom_point(data = df_xi_mean %>% filter(value > threshold[min(which(FDR < 0.05))]),
             aes(x = index, y = value),
             color = "red", size = 2) +
  scale_x_continuous(breaks = seq(0, 3500, by = 500)) +
  labs(x = "Taxa-covariate parameter index", y = "PIPs") +
  theme_minimal()

ggplot(df_xi_sum, aes(x = iteration, y = value)) +
  geom_point(color = "gray1") +
  labs(x = "Iteration", y = "Number of variables included") +
  theme_minimal()

ggplot(df_acf, aes(x = lag, y = acf)) +
  geom_col(fill = "gray1", width = 0.1) +
  geom_hline(yintercept = 0, color = "black") +
  labs(x = "Lag", y = "ACF") +
  theme_minimal()

which(xi_post_mean > threshold[min(which(FDR < 0.05))], arr.ind = TRUE)

# For illustration:
indices <- which(xi_post_mean > threshold[min(which(FDR < 0.05))], arr.ind = TRUE)

# Get covariate names
cov_names <- colnames(X)[indices[, "row"]]

# Get category names
category_names <- colnames(Y)[indices[, "col"]]#as.character((1:50)[indices[, "col"]])

# Build a data frame
df_assoc <- data.frame(
  Genera = category_names,
  Nutrients = cov_names
)

# Summarise in a table: list all covariates per category
library(dplyr)

summary_table <- df_assoc %>%
  group_by(Genera) %>%
  summarise(Nutrients = paste(Nutrients, collapse = ", "), .groups = "drop")

summary_table

# Replace acronyms with full names
df_assoc$NutrientsFull <- labels[df_assoc$Nutrients, 1]

# Build the summary table with full names
library(dplyr)

summary_table_full <- df_assoc %>%
  group_by(Genera) %>%
  summarise(Nutrients = paste(NutrientsFull, collapse = ", "), .groups = "drop")

summary_table_full

################################################################################
library(dplyr)
library(tidyr)

summary_table_full_sep <- summary_table_full %>%
  separate(
    Genera,
    into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
    sep = "\\.",
    fill = "right"
  ) %>%
  mutate(Species = NA) %>%  # Since species are not indicated
  select(Kingdom, Phylum, Class, Order, Family, Genus, Species, Nutrients)

# View the result
summary_table_full_sep

library(dplyr)
library(tidyr)

summary_table_long <- summary_table_full_sep %>%
  separate_rows(Nutrients, sep = ",\\s*")  # Split on comma and optional space

p1 <- ggplot(df_xi_mean, aes(x = index, y = value)) +
  geom_linerange(aes(ymin = 0, ymax = value), color = "gray1") +
  geom_point(data = df_xi_mean %>% filter(value > threshold[min(which(FDR < 0.05))]),
             aes(x = index, y = value),
             color = "red", size = 2) +
  scale_x_continuous(breaks = seq(0, 3500, by = 500)) +
  labs(x = "Taxa-covariate parameter index", y = "PIPs") +
  theme_minimal()

p2 <- ggplot(df_xi_sum, aes(x = iteration, y = value)) +
  geom_point(color = "gray1") +
  labs(x = "Iteration", y = "Number of variables included") +
  theme_minimal()

p3 <- ggplot(df_acf, aes(x = lag, y = acf)) +
  geom_col(fill = "gray1", width = 0.1) +
  geom_hline(yintercept = 0, color = "black") +
  labs(x = "Lag", y = "ACF") +
  theme_minimal()

print(p1)
print(p2)
print(p3)