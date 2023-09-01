# ------------------------------------------------------------------------------------------------ #
# Correlating fd with Fst
#
# Correlate fd against Fst. Negative relationship suggests that gene flow is homogenising the genome
# while positive would indicate the opposite.

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

# ------------------------------------------------------------------------------------------------ #
# Import data
df.fd <- fs::dir_ls(
  path = here("results", "between-species", "genomics_general"),
  glob = "*.csv"
) |>
  read_csv(
    col_names = TRUE,
    col_types = cols(),
    id = "P2_P3"
  ) |>
  mutate(
    P2_P3 = sub("_mgr-fd.csv", "", basename(P2_P3)),
    P2_P3 = toupper(sub("^.{3}_", "", P2_P3)),
    comparison = case_when(
      P2_P3 == "HEL_HMA" ~ "HEL-HMA",
      P2_P3 == "HMA_HCU" ~ "HCU-HMA",
      P2_P3 == "HEL_HCU" ~ "HCU-HEL"
    ),
    fd = ifelse(D < 0, 0, fd)
  ) |>
  rename(chromosome = scaffold) |>
  select(comparison, chromosome, start, end, fd) |>
  mutate(
    midpoint = floor((start + end)/2),
    comparison = factor(comparison, c("HCU-HEL", "HCU-HMA", "HEL-HMA"))
  )

div.stats <- read_csv(
  file = here("results", "between-species", "fst-zfst-dxy.csv"),
  col_names = TRUE,
  col_types = cols()
) |>
  mutate(comparison = factor(comparison, c("HCU-HEL", "HCU-HMA", "HEL-HMA")))

# Fst against fd
plot.fst_fd <- df.fd |>
  left_join(div.stats, join_by(chromosome, start, end, comparison)) |>
  filter(!if_any(everything(), is.na)) |>
  ggplot(
    aes(
      x = fst,
      y = fd
    )
  ) +
  labs(
    x = bquote(F[st]),
    y = bquote("\U0192"[d])
  ) +
  geom_point(alpha = 0.3, colour = "grey80") +
  geom_smooth(method = "lm") +
  ggpubr::stat_cor(
    digits = 3,
    label.x.npc = "left",
    label.y.npc = "top",
    p.accuracy = 0.001,
    show.legend = FALSE,
    cor.coef.name = "r",
    size = 5
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  scale_y_continuous(limits = c(0, 0.8)) +
  facet_wrap(facets = vars(comparison)) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    strip.text = element_text(size = 18, face = "bold")
  )

ragg::agg_png(
  filename = here("manuscript", "figures", "fst_fd-correlation.png"),
  width = 3000,
  height = 1500,
  units = "px",
  res = 300
)
plot.fst_fd
invisible(dev.off())

# ------------------------------------------------------------------------------------------------ #
# fd boxplot
df.fd |>
  ggplot(
    aes(
      x = comparison,
      y = fd,
      fill = comparison
    )
  ) +
  geom_boxplot() +
  scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb")) +
  scale_y_continuous(trans = "log2") +
  # scale_y_sqrt() +
  theme_bw()
