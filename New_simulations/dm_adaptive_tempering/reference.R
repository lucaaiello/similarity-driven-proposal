# Independent pure-R integer-count target and curvature for short audit datasets.
# No C++ helper, optimizer, or MCMC is called by these reference functions.
jb_rising <- function(la, count) {
  if (count == 0L) return(c(correction = 0, first = 0, second = 0))
  d <- log(0:(unname(count) - 1L)) - unname(la); p <- plogis(-d)
  c(correction = sum(pmax(d, 0) + log1p(exp(-abs(d)))), first = sum(p), second = sum(p * plogis(d)))
}
jb_likelihood <- function(Y, eta) sum(vapply(seq_len(nrow(Y)), function(i) {
  m <- max(eta[i, ]); ls <- m + log(sum(exp(eta[i, ] - m)))
  sum(vapply(seq_len(ncol(Y)), function(j) Y[i, j] * (eta[i, j] - ls) +
    jb_rising(eta[i, j], Y[i, j])["correction"], 0)) - jb_rising(ls, sum(Y[i, ]))["correction"]
}, 0))
jb_target <- function(alpha, B, xi, X, Y, s2, r2, a = 1, b = 9) {
  eta <- sweep(X %*% B, 2, alpha, "+")
  priors <- -.5 * sum(alpha^2 / s2)
  for (j in seq_len(ncol(Y))) priors <- priors + sum(dnorm(B[xi[, j] == 1L, j], 0, sqrt(r2[j]), log = TRUE))
  jb_likelihood(Y, eta) + priors + sum(lgamma(xi + a) + lgamma(1 - xi + b))
}
jb_raw_precision <- function(design, Y, alpha, xb, j, beta, s2, r2) {
  eta <- sweep(xb, 2, alpha, "+"); Z <- cbind(1, design)
  weights <- vapply(seq_len(nrow(Y)), function(i) {
    top <- max(eta[i, ]); ls <- top + log(sum(exp(eta[i, ] - top)))
    own <- jb_rising(eta[i, j], Y[i, j]); total <- jb_rising(ls, sum(Y[i, ]))
    p <- exp(eta[i, j] - ls)
    -(own["second"] - p * (1 - p) * total["first"] - p^2 * total["second"])
  }, 0)
  crossprod(Z, Z * weights) + diag(c(1 / s2, rep(1 / r2, length(beta))), length(beta) + 1L)
}
jb_project <- function(Q) {
  if (!inherits(try(chol(Q), silent = TRUE), "try-error")) return(Q)
  e <- eigen(Q, symmetric = TRUE); e$vectors %*% (pmax(e$values, .1) * t(e$vectors))
}
jb_logq <- function(x, mean, Q) {
  -.5 * (length(x) * log(2 * pi) - 2 * sum(log(diag(chol(Q)))) + drop(crossprod(x - mean, Q %*% (x - mean))))
}
# Rank-normalized folded split Rhat, cross-checked against posterior::rhat by audit.
jb_rhat <- function(x) {
  if (nrow(x) < 20L || any(!is.finite(x)) || any(apply(x, 2, function(v) length(unique(v)) == 1L))) return(NA_real_)
  h <- floor(nrow(x) / 2); s <- cbind(x[seq_len(h), , drop = FALSE], tail(x, h))
  score <- function(v) matrix(qnorm((rank(v, ties.method = "average") - .375) / (length(v) + .25)), nrow = h)
  basic <- function(z) { W <- mean(apply(z, 2, var)); B <- h * var(colMeans(z));
    if (W == 0) return(NA_real_); sqrt(((h - 1) * W / h + B / h) / W) }
  max(basic(score(s)), basic(score(abs(s - median(s)))))
}
