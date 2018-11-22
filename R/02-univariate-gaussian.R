source("01-prelim.R")

## ---- variational.example ----
# Simulate data from N(10, 4)
set.seed(123)
n <- 30
mu <- 0
psi <- 1
ydat <- rnorm(n, mean = 0, sd = sqrt(1 / psi))
s <- var(ydat) * (n - 1) / n

# pdf of normal gamma
dnormgamma <- function(x, mu = 0, lambda = 1, alpha = 1, beta = 1) {
  const <- alpha * log(beta) + 0.5 * log(lambda) - lgamma(alpha) - 0.5 * log(2 * pi)
  res <- const + (alpha - 0.5) * log(x[2]) - beta * x[2] -
    lambda * x[2] * (x[1] - mu) ^ 2 / 2
  exp(res)
}

# pdf of variational approximation
dvarapprox <- function(x, mean.mu, var.mu, shape.psi, rate.psi) {
  dnorm(x[1], mean = mean.mu, sd = sqrt(var.mu)) *
    dgamma(x[2], shape = shape.psi, rate = rate.psi)
}

# Lower bound
varLowerBound <- function(an, bn, mun, kappan) {
  L1 <- -0.5 * log(kappan) + 0.5 * (1 + log(2 * pi))
  L2 <- lgamma(an) - (an - 1) * digamma(an) - log(bn) + an
  L3 <- (-n / 2) * log(2 * pi) + (n / 2) * (digamma(an) - log(bn)) -
    (n / 2 * an / bn) * (s + mean(ydat) ^ 2 - 2 * mean(ydat) * mun + mun ^ 2 +
                           1 / kappan)
  L4 <- (alpha0 - 1) * (digamma(an) - log(bn)) - beta0 * an / bn +
    alpha0 * log(beta0) - lgamma(alpha0)
  L5 <- 0.5 * log(lambda0 / (2*pi)) + 0.5 * (digamma(an) - log(bn)) -
    lambda0 / 2 * an / bn * (1 / kappan + (mun - mu0) ^ 2)
  L1 + L2 + L3 + L4 + L5
}

# Exact parameters of the posterior
lambda0 <- 10
mu0 <- 0
alpha0 <- 5
beta0 <- 4
mu.post <- (lambda0 * mu0 + n * mean(ydat)) / (lambda0 + n)
lambda.post <- lambda0 + n
alpha.post <- alpha0 + n / 2
beta.post <- beta0 + 0.5 * (n * s + lambda0 * n * (mean(ydat) - mu0) ^ 2 / (lambda0 + n))
(mean.post <- mu.post)  # estimate of mean
(prec.post <- alpha.post / beta.post)  # estimate of precision
marginal.log.lik <- sum(dnorm(ydat, mean = mu.post, sd = sqrt(1 / prec.post),
                              log = TRUE))

# Data for contour plot
plotxlim <- c(-0.4, 0.5)
plotylim <- c(0.6, 2)
x <- seq(plotxlim[1], plotxlim[2], length = 200)  # mean
y <- seq(plotylim[1], plotylim[2], length = 200)  # variance
xy <- expand.grid(x, y); names(xy) <- c("x", "y")
z <- apply(xy, 1, dnormgamma, mu = mu.post, lambda = lambda.post,
           alpha = alpha.post, beta = beta.post)
contour.df <- data.frame(xy, z, iteration = "Truth", logLik = NA, alpha = 1)

# Function for variational inference contour data
varContourDat <- function(niter = 3, alpha = 0.5) {
  res <- NULL; LB <- rep(NA, 1 + niter * 2)
  varLabel <- function(x) {
    xx <- rep(1:50, each = 2)[x]
    ifelse(!isEven(x),
           paste0("Iteration ", xx, " (mu update)"),
           paste0("Iteration ", xx, " (psi update)"))
  }

  # Initialise
  q.psi.c <- 60; q.psi.d <- 40
  q.mu.a <- 0.4; q.mu.b <- 400
  q.mu.var <- q.psi.d / (q.psi.c * q.mu.b)
  res[[1]] <- apply(xy, 1, dvarapprox, mean.mu = q.mu.a, var.mu = q.mu.var,
                  shape.psi = q.psi.c, rate.psi = q.psi.d)
  LB[1] <- varLowerBound(q.psi.c, q.psi.d, q.mu.a, q.mu.b)  # LB init

  # Variational updates
  for (i in 1:niter) {
    # Update mu
    q.mu.a <- (lambda0 * mu0 + n * mean(ydat)) / (lambda0 + n)
    q.mu.b <- (lambda0 + n) * q.psi.c / q.psi.d
    q.mu.var <- 1 / q.mu.b
    res[[2 * i]] <-
      apply(xy, 1, dvarapprox, mean.mu = q.mu.a, var.mu = q.mu.var,
            shape.psi = q.psi.c, rate.psi = q.psi.d)
    LB[2 * i] <- varLowerBound(q.psi.c, q.psi.d, q.mu.a, q.mu.b)  # LB (mu)

    # Update psi
    q.psi.c <- alpha0 + n / 2
    q.psi.d <- beta0 + 0.5 * (sum(ydat ^ 2) - 2 * mean(ydat) * q.mu.a +
                                q.mu.var + q.mu.a ^ 2)
    res[[2 * i + 1]] <-
      apply(xy, 1, dvarapprox, mean.mu = q.mu.a, var.mu = q.mu.var,
            shape.psi = q.psi.c, rate.psi = q.psi.d)
    LB[2 * i + 1] <- varLowerBound(q.psi.c, q.psi.d, q.mu.a, q.mu.b)  # LB (psi)
  }

  tmp <- lapply(res, function(w) data.frame(xy, z = w))
  df <- data.frame(tmp[[1]], iteration = "Iteration 0", logLik = LB[1],
                   alpha = alpha)
  for (i in 2:(niter * 2 + 1)) {
    df <- rbind(df, data.frame(tmp[[i]], iteration = varLabel(i - 1),
                               logLik = LB[i], alpha = alpha))
  }
  df
}
contour.df <- rbind(contour.df, varContourDat())

# Contour plot
varContourPlot <- function(iter = NULL, animate = FALSE) {
  var.lb <- unique(contour.df[, 5])[-1]
  min.var.lb <- min(var.lb)
  var.iter.numbers <- as.numeric(unique(contour.df[, 4]))
  my.red <- gg_colour_hue(2)[1]
  my.green <- gg_colour_hue(2)[2]

  if (is.null(iter)) {
    df <- subset(contour.df, as.numeric(iteration) == 1)
    plot.title <- " "
    current.lb <- var.lb[1]
    legend.lab <- c("", "log p(y)")
  } else {
    df <- subset(contour.df, as.numeric(iteration) == c(1, iter + 2))
    iter.plot.to.show <- var.iter.numbers[iter + 1] - 1
    if (iter.plot.to.show == 0) {
      plot.title <- bquote(Iteration~0~(initialisation))
    } else if (!isEven(iter.plot.to.show)) {
      iter.plot.to.show <- (iter.plot.to.show + 1) / 2
      plot.title <- bquote(Iteration~.(iter.plot.to.show)~(mu~update))
    } else {
      iter.plot.to.show <- iter.plot.to.show / 2
      plot.title <- bquote(Iteration~.(iter.plot.to.show)~(psi~update))
    }
    current.lb <- var.lb[iter + 1]
    legend.lab <- c(expression(italic(L)(q)), "log p(y)")
  }

  ggplot(data = df, aes(x = x, y = y, z = z, group = iteration, col = logLik,
                        linetype = iteration, size = iteration)) +
    geom_contour() +
    # geom_point(aes(x = 0, y = 1), size = 3, col = "grey10") +
    coord_cartesian(xlim = plotxlim, ylim = plotylim) +
    labs(x = expression(Mean~(mu)), y = expression(Precision~(psi)),
         title = plot.title) +
    scale_colour_gradient(name = NULL,
                          breaks = c(current.lb, marginal.log.lik),
                          labels = legend.lab,
                          high = my.green, low = my.red,
                          limits = c(min.var.lb - 0.4, marginal.log.lik),
                          na.value = alpha(my.green, 0.25)) +
    scale_linetype_manual(values = c(1, 13), guide = FALSE) +
    scale_size_manual(values = c(1.9, 0.7), guide = FALSE) +
    # geom_dl(aes(label = iteration), method = "first.bumpup") +
    theme_bw() +
    theme(legend.text = element_text(size = 10), legend.text.align = 0,
          legend.justification = "center") +
    guides(col = guide_colorbar(barwidth = 1, barheight = 16))
}

## ---- variational.example.save ----
ggsave("../figure/vbupdate.pdf", varContourPlot(), width = 6.5, height = 6.5 / 1.5)
plot.names <- c(1:6, 0)
for (i in 1:length(plot.names)) {
  ggsave(paste0("../figure/vbupdate_", i,".pdf"), varContourPlot(plot.names[i]),
         width = 6.5, height = 6.5 / 1.5)
}

## ---- variational.example.gganimate ----
makeplot <- function() {
  plot.names <- 0:6
  for (i in 1:length(plot.names)) {
    p <- varContourPlot(plot.names[i])
    print(p)
  }
  ani.pause()
}
saveGIF(makeplot(), interval = 1, ani.width = 550, ani.height = 550 / 1.5)
# Move animate.gif to ../figure and also ./R/figure
file.copy("animation.gif", file.path("../figure", "animation.gif"))
file.rename("animation.gif", "figure/animation.gif")
