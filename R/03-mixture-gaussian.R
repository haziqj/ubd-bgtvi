# Adapted from https://rpubs.com/cakapourani/variational-bayes-gmm

suppressPackageStartupMessages(require(matrixcalc))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(mvtnorm))
suppressPackageStartupMessages(library(Matrix))

# Define ggplot2 theme
gg_theme <- function(){
  p <- theme(
    plot.title = element_text(size = 20,face = 'bold',
                              margin = margin(0,0,3,0), hjust = 0.5),
    axis.text = element_text(size = rel(1.05), color = 'black'),
    axis.title = element_text(size = rel(1.45), color = 'black'),
    axis.title.y = element_text(margin = margin(0,10,0,0)),
    axis.title.x = element_text(margin = margin(10,0,0,0)),
    axis.ticks.x = element_line(colour = "black", size = rel(0.8)),
    axis.ticks.y = element_blank(),
    legend.position = "right",
    legend.key.size = unit(1.4, 'lines'),
    legend.title = element_text(size = 12, face = 'bold'),
    legend.text = element_text(size = 12),
    panel.border = element_blank(),
    panel.grid.major = element_line(colour = "gainsboro"),
    panel.background = element_blank()
  )
  return(p)
}

# Mixture density using approximate predictive gaussian distribution
mixture_pdf_gaussian <- function(model, data){
  mixture <- vector(mode = "numeric", length = NROW(data))
  for (k in 1:length(model$nu)) {
    tau_k <- model$W[,,k] * model$nu[k] # TODO: Is this right?
    mu_k  <- model$m[, k]
    mixture <- mixture + model$pi_k[k] *
      dmvnorm(cbind(data$x,data$y), mean = mu_k, sigma = solve(tau_k))
  }
  return(mixture)
}
# Mixture density using predictive t-distribution
mixture_pdf_t <- function(model, data){
  mixture <- vector(mode = "numeric", length = NROW(data))
  for (k in 1:length(model$nu)) {
    L_k <- solve((((model$nu[k] + 1 - NROW(model$m)) * model$beta[k]) /
                    (1 + model$beta[k])) * model$W[,,k])
    mixture <- mixture + (model$alpha[k]/sum(model$alpha)) *
      dmvt(x = cbind(data$x,data$y), delta = model$m[, k],
           sigma = L_k, df = model$nu[k] + 1 - NROW(model$m),
           log = FALSE, type = "shifted")
  }
  return(mixture)
}
# Use the log sum exp trick for having numeric stability
log_sum_exp <- function(x) {
  # Computes log(sum(exp(x))
  offset <- max(x)
  s <- log(sum(exp(x - offset))) + offset
  i <- which(!is.finite(s))
  if (length(i) > 0) { s[i] <- offset }
  return(s)
}

# Fit VBLR model
vb_gmm <- function(X, K = 3, alpha_0 = 1/K, m_0 = c(colMeans(X)),
                   beta_0 = 1, nu_0 = NCOL(X) + 50,
                   W_0 = diag(100, NCOL(X)), max_iter = 500,
                   epsilon_conv = 1e-4, is_animation = FALSE,
                   is_verbose = FALSE){

  res.mu <- res.Lambda <- list()

  # Compute logB function
  logB <- function(W, nu){
    D <- NCOL(W)
    return(-0.5*nu*log(det(W)) - (0.5*nu*D*log(2) + 0.25*D*(D - 1) *
                                    log(pi) + sum(lgamma(0.5 * (nu + 1 - 1:NCOL(W)))) ))  #log of B.79
  }
  X <- as.matrix(X)
  D <- NCOL(X)              # Number of features
  N <- NROW(X)              # Number of observations
  W_0_inv <- solve(W_0)     # Compute W^{-1}
  L <- rep(-Inf, max_iter)  # Store the lower bounds
  r_nk = log_r_nk = log_rho_nk <- matrix(0, nrow = N, ncol = K)
  x_bar_k <- matrix(0, nrow = D, ncol = K)
  S_k = W_k <- array(0, c(D, D, K ) )
  log_pi = log_Lambda <- rep(0, K)

  if (is_animation) {
    # Variables needed for plotting
    dt <- data.table(expand.grid(x = seq(from = min(X[,1]) - 2,
                                         to = max(X[,1]) + 2, length.out = 80),
                                 y = seq(from = min(X[,2]) - 8,
                                         to = max(X[,2]) + 2, length.out = 80)))
  }
  dt_all <- data.table(x = numeric(), y = numeric(), z = numeric(),
                       iter = numeric())

  m_k     <- t(kmeans(X, K, nstart = 25)$centers)  # Mean of Gaussian
  # m_k <- rbind(
  #   sample(X[, 1], size = K, replace = FALSE),
  #   sample(X[, 2], size = K, replace = FALSE)
  # )
  beta_k  <- rep(beta_0, K)                # Scale of precision matrix
  nu_k    <- rep(nu_0, K)                  # Degrees of freedom
  alpha   <- rep(alpha_0, K)               # Dirichlet parameter
  log_pi  <- digamma(alpha) - digamma(sum(alpha))
  for (k in 1:K) {
    W_k[,,k] <-  W_0  # Scale matrix for Wishart
    log_Lambda[k] <- sum(digamma((nu_k[k] + 1 - c(1:D))/2)) +
      D*log(2) + log(det(W_k[,,k]))
  }

  res.mu[[1]] <- m_k
  Lambda.tmp <- W_k
  for (k in seq_len(K)) {
    v1 <- runif(1, min = 0.005, max = 0.2)
    v2 <- runif(1, min = 185 / 10, max = 185 / 10)
    Lambda.tmp[, , k] <- diag(1 / c(v1, v2))
  }
  res.Lambda[[1]] <- Lambda.tmp

  if (is_animation) { # Create animation for initial assignments
    my_z = mixture_pdf_t(model = list(m = m_k, W = W_k, beta = beta_k,
                                      nu = nu_k, alpha = rep(1/K, K)), data = dt)
    dt_all <- rbind(dt_all, dt[, z := my_z] %>% .[, iter := 0])
  }

  # Iterate to find optimal parameters
  for (i in 2:max_iter) {
    ##-------------------------------
    # Variational E-Step
    ##-------------------------------
    for (k in 1:K) {
      diff <- sweep(X, MARGIN = 2, STATS = m_k[, k], FUN = "-")
      log_rho_nk[, k] <- log_pi[k] + 0.5*log_Lambda[k] -
        0.5*(D/beta_k[k]) - 0.5*nu_k[k] * diag(diff %*%
                                                 W_k[,,k] %*% t(diff)) # log of 10.67
    }
    # Responsibilities using the logSumExp for numerical stability
    Z        <- apply(log_rho_nk, 1, log_sum_exp)
    log_r_nk <- log_rho_nk - Z              # log of 10.49
    r_nk     <- apply(log_r_nk, 2, exp)     # 10.49

    ##-------------------------------
    # Variational M-Step
    ##-------------------------------
    N_k <- colSums(r_nk) + 1e-10  # 10.51
    for (k in 1:K) {
      x_bar_k[, k] <- (r_nk[ ,k] %*% X) / N_k[k]   # 10.52
      x_cen        <- sweep(X,MARGIN = 2,STATS = x_bar_k[, k],FUN = "-")
      S_k[, , k]   <- t(x_cen) %*% (x_cen * r_nk[, k]) / N_k[k]  # 10.53
    }
    # Update Dirichlet parameter
    alpha <- alpha_0 + N_k  # 10.58
    # # Compute expected value of mixing proportions
    # pi_k <- (alpha + N_k) / (K * alpha_0 + N)
    pi_k <- alpha / sum(alpha)
    # Update parameters for Gaussia-nWishart distribution
    beta_k <- beta_0 + N_k    # 10.60
    nu_k   <- nu_0 + N_k   # 10.63
    for (k in 1:K) {
      # 10.61
      m_k[, k]   <- (1/beta_k[k]) * (beta_0*m_0 + N_k[k]*x_bar_k[, k])
      # 10.62
      W_k_inv <- W_0_inv + N_k[k] * S_k[,,k] +
        ((beta_0*N_k[k])/(beta_0 + N_k[k])) * tcrossprod((x_bar_k[, k] - m_0))
      W_k[, , k] <- solve(W_k_inv)
    }
    # Update expectations over \pi and \Lambda
    # 10.66
    log_pi <- digamma(alpha) - digamma(sum(alpha))
    for (k in 1:K) { # 10.65
      log_Lambda[k] <- sum(digamma((nu_k[k] + 1 - 1:D)/2)) +
        D*log(2) + log(det(W_k[,,k]))
    }

    res.mu[[i]] <- m_k
    res.Lambda[[i]] <- array(rep(nu_k, each = K), dim = c(D, D, K)) * W_k

    ##-------------------------------
    # Variational lower bound
    ##-------------------------------
    lb_px = lb_pml = lb_pml2 = lb_qml <- 0
    for (k in 1:K) {
      # 10.71
      lb_px <- lb_px + N_k[k] * (log_Lambda[k] - D/beta_k[k] - nu_k[k] *
                                   matrix.trace(S_k[,,k] %*% W_k[,,k]) - nu_k[k]*t(x_bar_k[,k] -
                                                                                     m_k[,k]) %*% W_k[,,k] %*% (x_bar_k[,k] - m_k[,k]) - D*log(2*pi) )
      # 10.74
      lb_pml <- lb_pml + D*log(beta_0/(2*pi)) + log_Lambda[k] -
        (D*beta_0)/beta_k[k] - beta_0*nu_k[k]*t(m_k[,k] - m_0) %*%
        W_k[,,k] %*% (m_k[,k] - m_0)
      # 10.74
      lb_pml2 <- lb_pml2 + nu_k[k] * matrix.trace(W_0_inv %*% W_k[,,k])
      # 10.77
      lb_qml <- lb_qml + 0.5*log_Lambda[k] + 0.5*D*log(beta_k[k]/(2*pi)) -
        0.5*D - logB(W = W_k[,,k], nu = nu_k[k]) -
        0.5*(nu_k[k] - D - 1)*log_Lambda[k] + 0.5*nu_k[k]*D
    }
    lb_px  <- 0.5 * lb_px             # 10.71
    lb_pml <- 0.5*lb_pml + K*logB(W = W_0,nu = nu_0) + 0.5*(nu_0 - D - 1) *
      sum(log_Lambda) - 0.5*lb_pml2 # 10.74
    lb_pz  <- sum(r_nk %*% log_pi)    # 10.72
    lb_qz  <- sum(r_nk * log_r_nk)    # 10.75
    lb_pp  <- sum((alpha_0 - 1)*log_pi) + lgamma(sum(K*alpha_0)) -
      K*sum(lgamma(alpha_0))        # 10.73
    lb_qp  <- sum((alpha - 1)*log_pi) + lgamma(sum(alpha)) -
      sum(lgamma(alpha)) # 10.76
    # Sum all parts to compute lower bound
    L[i] <- lb_px + lb_pz + lb_pp + lb_pml - lb_qz - lb_qp - lb_qml

    ##-------------------------------
    # Evaluate mixture density for plotting
    ##-------------------------------
    if (is_animation) {
      if ( (i - 1) %% 5 == 0 | i < 10) {
        my_z = mixture_pdf_t(model = list(m = m_k, W = W_k, beta = beta_k,
                                          nu = nu_k, alpha = alpha), data = dt)
        dt_all <- rbind(dt_all, dt[, z := my_z] %>% .[, iter := i - 1])
      }
    }
    # Show VB difference
    if (is_verbose) { cat("It:\t",i,"\tLB:\t",L[i],
                          "\tLB_diff:\t",L[i] - L[i - 1],"\n")}
    # Check if lower bound decreases
    if (L[i] < L[i - 1]) { message("Warning: Lower bound decreases!\n"); }
    # Check for convergence
    if (abs(L[i] - L[i - 1]) < epsilon_conv) { break }
    # Check if VB converged in the given maximum iterations
    if (i == max_iter) {warning("VB did not converge!\n")}
  }
  obj <- structure(list(X = X, K = K, N = N, D = D, pi_k = pi_k,
                        alpha = alpha, r_nk = r_nk,  m = m_k, W = W_k,
                        beta = beta_k, nu = nu_k, L = L[2:i],
                        dt_all = dt_all, res.mu = res.mu, res.Lambda = res.Lambda),
                   class = "vb_gmm")
  cat("Iter = ", i, "\n")
  return(obj)
}

set.seed(123)  # For reproducibility
X <- as.matrix(faithful)
K <- 6     # Number of clusters
# Run vb-gmm model model
res <- vb_gmm_model <- vb_gmm(X = X, K = K, alpha_0 = 1e-5, max_iter = 300,
                       is_animation = TRUE)
                       # W_0 = diag(c(10, 0.05)) / 100)
# which(diff(diff(res$L)) > 0.01)  #1, 12, 149

# Data plot
ggplot(faithful, aes(eruptions, waiting)) +
  geom_point() +
  labs(x = "Eruption length (minutes)", y = "Time to next eruption (minutes)") +
  theme_bw() -> p
ggExtra::ggMarginal(p, type = "histogram", fill = "grey75", col = "grey75") -> p

ggsave("../figure/faithful_scatter.pdf", p, width = 7, height = 7 * 4 / 7)

# Plot of ELBO
elbows <- c(3, 22, 28, 155)
ggplot() +
  geom_line(aes(x = seq_along(res$L), y = res$L)) +
  geom_point(aes(x = elbows, y = res$L[elbows]), shape = 1, size = 5,
             col = iprior::gg_col_hue(1)) +
  labs(x = "Iteration", y = "Lower bound") +
  theme_bw()

# Plot of 1-sd contours
data_fun <- function(mean.vec = apply(X, 2, mean),
                     prec.mat = diag(1 / apply(X, 2, var))) {
  set.seed(123)
  res <- as.data.frame(rmvnorm(1000, mean = mean.vec, sigma = solve(prec.mat)))
  colnames(res) <- c("x1", "x2")
  res
}

plot_faithful_iter <- function(iter = 1) {
  mu <- res$res.mu[[iter]]
  Psi <- res$res.Lambda[[iter]]

  ggplot(faithful, aes(eruptions, waiting)) +
    geom_point() +
    labs(x = "Eruption length (minutes)", y = "Time to next eruption (minutes)") +
    theme_bw() -> p

  K <- ncol(mu)
  cols <- iprior::gg_col_hue(K)
  for (k in 1:K) {
    dat <- data_fun(mu[, k], (Psi[, , k] + t(Psi[, , k])) / 2)
    p <- p + stat_ellipse(data = dat, aes(x1, x2), geom = "polygon",
                          alpha = 0.3, fill = cols[k], level = 0.6827)
  }

  the.text <- paste0("Iteration ", iter - 1)
  return(p + annotate(geom = "text", label = the.text, x = -Inf, y = Inf, hjust = -0.1,
               vjust = 1.5))
}

ggsave("../figure/faithful_iter_1.pdf", plot_faithful_iter(1),
       width = 7, height = 7 * 4.5 / 7)
ggsave("../figure/faithful_iter_2.pdf", plot_faithful_iter(20),
       width = 7, height = 7 * 4.5 / 7)
ggsave("../figure/faithful_iter_3.pdf", plot_faithful_iter(30),
       width = 7, height = 7 * 4.5 / 7)
ggsave("../figure/faithful_iter_4.pdf", plot_faithful_iter(158),
       width = 7, height = 7 * 4.5 / 7)

# Save GIF
makeplot <- function() {
  for (i in c(1:31, seq(41, 158, by = 10), rep(158, 10))) {
    p <- plot_faithful_iter(i) + coord_cartesian(xlim = c(1.5, 5.2), ylim = c(43, 97))
    print(p)
  }
  animation::ani.pause()
}
animation::saveGIF(makeplot(), interval = 0.25, ani.width = 550, ani.height = 550 / 1.5)
