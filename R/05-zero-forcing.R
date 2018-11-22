source("01-prelim.R")

px <- function(x) 0.65 * dnorm(x, -4, 1.4) + 0.35 * dnorm(x, 5, 1.85)
x <- seq(-10, 12, length.out = 1000)

KL_forward <- function(theta) {
  mu <- theta[1]
  sigma <- theta[2]
  res <- integrate(
    function(x) {
      (log(px(x)) - dnorm(x, mean = mu, sd = sigma, log = TRUE)) * px(x)
    }, lower = -30, upper = 30
  )
  res$value
}

KL_reverse <- function(theta) {
  mu <- theta[1]
  sigma <- theta[2]
  res <- integrate(
    function(x) {
      (dnorm(x, mean = mu, sd = sigma, log = TRUE) - log(px(x))) * dnorm(x, mean = mu, sd = sigma)
    }, lower = -50, upper = 50
  )
  res$value
}

# Zero-forcing function (reverse KL)
res <- optim(c(-4, 1), KL_reverse, method = "L-BFGS-B", lower = c(-Inf, 1e-9))
theta1 <- res$par  # left mode
res <- optim(c(5, 2), KL_reverse, method = "L-BFGS-B", lower = c(-Inf, 1e-9))
theta2 <- res$par  # right mode

# Zero-avoiding function (forward KL)
res <- optim(c(0, 1), KL_forward, method = "L-BFGS-B", lower = c(-Inf, 1e-9))
theta3 <- res$par

# Construct data frame for plotting
plot.df <- tibble(x = x, y = px(x),
                  zf1 = dnorm(x = x, mean = theta1[1], sd = theta1[2]),
                  zf2 = dnorm(x = x, mean = theta2[1], sd = theta2[2]),
                  za  = dnorm(x = x, mean = theta3[1], sd = theta3[2]))

tmp <- iprior::gg_col_hue(2)  # ggplot colours

# PLOT 1
plot.df %>%
  select(x, y) %>%
  reshape2::melt(id = "x") %>%
  subset(value > 0.0001) %>%
  ggplot(aes(x, value, col = variable)) +
  geom_line(size = 1) +
  scale_colour_manual(name = "", labels = c("p(z)"),
                      values = c("black")) +
  theme_void() +
  coord_cartesian(ylim = c(0, 0.28), xlim = c(-10, 12)) +
  theme(legend.position = c(1, 0.95), legend.justification = c("right", "top")) -> p1

# PLOT 2 (zero forcing left mode)
plot.df %>%
  select(x, y, zf1) %>%
  reshape2::melt(id = "x") %>%
  subset(value > 0.0001) %>%
  ggplot(aes(x, value, col = variable)) +
  geom_line(size = 1) +
  annotate("text", x = -0.1, y = 0.2, label = "KL(q || p) is small", size = 5,
           col = tmp[1]) +
  scale_colour_manual(name = "", labels = c("p(z)", "q(z)"),
                      values = c("black", tmp[1])) +
  theme_void() +
  coord_cartesian(ylim = c(0, 0.28), xlim = c(-10, 12)) +
  theme(legend.position = c(1, 0.95), legend.justification = c("right", "top")) -> p2

# PLOT 3 (zero forcing right mode)
plot.df %>%
  select(x, y, zf2) %>%
  reshape2::melt(id = "x") %>%
  subset(value > 0.0001) %>%
  ggplot(aes(x, value, col = variable)) +
  geom_line(size = 1) +
  annotate("text", x = 8.9, y = 0.16, label = "KL(q || p) is small", size = 5,
           col = tmp[1]) +
  scale_colour_manual(name = "", labels = c("p(z)", "q(z)"),
                      values = c("black", tmp[1])) +
  theme_void() +
  coord_cartesian(ylim = c(0, 0.28), xlim = c(-10, 12)) +
  theme(legend.position = c(1, 0.95), legend.justification = c("right", "top")) -> p3

# PLOT 4 (forward KL/zero avoiding)
plot.df %>%
  select(x, y, za) %>%
  reshape2::melt(id = "x") %>%
  subset(value > 0.0001) %>%
  ggplot(aes(x, value, col = variable)) +
  geom_line(size = 1) +
  annotate("text", x = 2, y = 0.0962, label = "KL(p || q) is small", size = 5,
           col = tmp[2]) +
  scale_colour_manual(name = "", labels = c("p(z)", "q(z)"),
                      values = c("black", tmp[2])) +
  theme_void() +
  coord_cartesian(ylim = c(0, 0.28), xlim = c(-10, 12)) +
  theme(legend.position = c(1, 0.95), legend.justification = c("right", "top")) -> p4

ggsave("../figure/zero_force_1.pdf", p1, width = 7.5, height = 5)
ggsave("../figure/zero_force_2.pdf", p2, width = 7.5, height = 5)
ggsave("../figure/zero_force_3.pdf", p3, width = 7.5, height = 5)
ggsave("../figure/zero_force_4.pdf", p4, width = 7.5, height = 5)


