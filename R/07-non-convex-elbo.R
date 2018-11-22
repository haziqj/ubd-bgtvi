source("01-prelim.R")

dat <- gen_mixture(m = 3, n = 100)
lb <- NA
for (i in 1:10) {
  mod <- iprobit(y ~ ., dat, kernel = "fbm", common.RKHS.scale = TRUE,
                 control = list(maxit = 100))
  tmp <- mod$lower.bound
  if (length(tmp) < 100) tmp <- c(tmp, rep(max(tmp), 100 - length(tmp)))
  lb <- rbind(lb, cbind(x = 1:100, lb = tmp, i = i))
}

lb <- lb[-1, ]

p1 <- ggplot(subset(as.data.frame(lb), i != 3)) +
  geom_line(aes(x, lb, group = i, col = i), size = 0.3) +
  labs(x = "Iteration", y = "Lower bound") +
  theme_bw() +
  theme(legend.position = "none")

p2 <- ggplot(subset(as.data.frame(lb), i != 3)) +
  geom_line(aes(x, lb, group = i, col = i), size = 0.3) +
  labs(x = "Iteration", y = "Lower bound") +
  ggforce::facet_zoom(y = lb > -48) +
  theme_bw() +
  theme(legend.position = "none")

ggsave("../figure/elbo_local_1.pdf", p1, width = 6, height = 6 * 3.7 / 7)
ggsave("../figure/elbo_local_2.pdf", p2, width = 6, height = 6 * 3.7 / 7)


