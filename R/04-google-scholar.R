source("01-prelim.R")

## ---- extract.results ----
system2("python", "extract_occurrences.py 'variational inference' 1950 2017")
vi.res <- read_csv("out.csv")
system2("python", "extract_occurrences.py 'EM algorithm' 1950 2017")
em.res <- read_csv("out.csv")

## ---- hist.plots ----
# as.tibble(rbind(
#   cbind(vi.res, type = "Variational inference"),
#   cbind(em.res, type = "EM algorithm")
# )) -> tmp
# ggplot(subset(tmp, year >= 1980)) +
#   geom_freqpoly(aes(x = year, y = results, col = type),  stat = "identity") +
#   # facet_grid(. ~ type, space = "free", scales = "free") +
#   geom_dl(aes(year, results, label = type, col = type), method = "lines2") +
#   # scale_y_log10() +
#   coord_trans(y = "log10") +
#   labs(x = "Year", y = "Hits on Google Scholar") +
#   theme_bw() +
#   theme(legend.position = "none")
ggplot(subset(vi.res, year >= 1980)) +
  geom_bar(aes(year, results), stat = "identity") +
  labs(x = "Year", y = "Hits", title = "Google Scholar results for 'variational inference'") +
  theme_light() +
  theme(legend.position = "none") -> p

## ---- save.plots.for.presentation ----
ggsave("../figure/google_scholar.pdf", p, width = 7, height = 7 / 3)
