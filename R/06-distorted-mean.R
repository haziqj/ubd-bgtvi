source("01-prelim.R")

x <- seq(-5, 5, length = 250)
x <- expand.grid(x1 = x, x2 = x)

Sigma <- matrix(c(1, 0.98, 0.98, 1), 2, 2)
Psi <- solve(Sigma)

px <- function(x) mvtnorm::dmvnorm(x, mean = rep(0, 2), sigma = Sigma)
p.df <- cbind(x, z = apply(as.matrix(x), 1, px))
qx <- function(x) mvtnorm::dmvnorm(x, mean = rep(0, 2),
                                   sigma = diag(1 / diag(Psi)))
q.df <- cbind(x, z = apply(as.matrix(x), 1, qx))

plot.df <- rbind(
  cbind(p.df, type = "p(z)"),
  cbind(q.df, type = "q(z)")
)

cols <- iprior::gg_col_hue(2)

# PLOT 1
ggplot(p.df) +
  geom_contour(aes(x1, x2, z = z), col = cols[1], alpha = 0.8) +
  labs(x = expression(z[1]), y = expression(z[2])) +
  annotate("text", x = 2, y = 1.5, label = "p(z)", col = cols[1], size = 4) +
  theme_bw() +
  theme(legend.position = "none") -> p1; p1

# PLOT 2
ggplot(plot.df) +
  geom_contour(aes(x1, x2, z = z, col = type), alpha = 0.8, bins = 6) +
  labs(x = expression(z[1]), y = expression(z[2])) +
  annotate("text", x = 2, y = 1.5, label = "p(z)", col = cols[1], size = 4) +
  annotate("text", x = 0.42, y = -0.42, label = "q(z)", col = cols[2], size = 4) +
  theme_bw() +
  theme(legend.position = "none") -> p2; p2

width <- 6
ggsave("../figure/distorted_mean_1.pdf", p1, width = width, height = width * 3 / 7)
ggsave("../figure/distorted_mean_2.pdf", p2, width = width, height = width * 3 / 7)
