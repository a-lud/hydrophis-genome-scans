# ------------------------------------------------------------------------------------------------ #
# Compare statistics in outlier/profiled regions to genomic background
#
#

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(patchwork)
})

# ------------------------------------------------------------------------------------------------ #
# Read in outlier windows: ZFst & Profiled
zfst <- read_csv(
  file = here("results", "between-species", "zfst-outliers.csv"),
  col_names = TRUE,
  col_types = cols()
)

div.prof <- read_tsv(
  here("results", "tables", "classified-windows.tsv"),
  col_names = TRUE,
  col_types = cols()
) |>
  select(comparison, chromosome, start, end, fst, zfst, dxy, profile) |>
  filter(profile != "Barrier")

# 1. Get unique windows and their profile (elevated/reduced Fst)
df.outliers <- bind_rows(zfst, div.prof) |>
  mutate(
    window_type = case_when(
      classification == -3 ~ "High diversity",
      classification == 3 ~ "Low diversity",
      .default = NA_character_
    ),
    window_type = case_when(
      profile == "Balancing" ~ "Low diversity",
      profile %in% c("Allopatric", "Recurrent") ~ "High diversity",
      .default = window_type
    )
  ) |>
  select(comparison, chromosome, start, end, window_type) |>
  distinct()

# 2. get only the unique windows
tmp <- df.outliers |>
  select(chromosome, start, end, window_type) |>
  distinct()

# ------------------------------------------------------------------------------------------------ #
# Get genomic background data
df.div.stats <- read_csv(
  file = here("results/between-species/fst-zfst-dxy.csv"),
  col_names = TRUE,
  col_types = cols()
)

# Nucleotide diversity (pi)
df.pi <- fs::dir_ls(path = here("results", "within-species"), glob = "*pi.txt") %>%
  magrittr::extract(str_detect(., "50kb")) %>%
  magrittr::extract(str_detect(., "HMA|HEL|HCU")) |>
  map(read_tsv, col_names = TRUE, col_types = cols()) |>
  list_rbind() |>
  filter(no_sites >= 5000) |>
  rename(comparison = pop, pi = avg_pi) |>
  rename(start = window_pos_1, end = window_pos_2) |>
  select(-c(no_sites, count_diffs, count_comparisons, count_missing))

# Tajima's D
df.tajD <- fs::dir_ls(path = here("results", "within-species"), glob = "*TajimaD.tsv") %>%
  magrittr::extract(str_detect(., "50kb")) %>%
  magrittr::extract(str_detect(., "HMA|HEL|HCU")) |>
  read_tsv(col_names = TRUE, col_types = cols(), id = "comparison") |>
  mutate(
    comparison = sub("\\..*", "", basename(comparison)),
    end = BIN_START + 5e4,
    start = BIN_START + 1
  ) |>
  select(comparison, chromosome = CHROM, start, end, TajimaD)

df.pi_tajD <- inner_join(df.pi, df.tajD)

# ------------------------------------------------------------------------------------------------ #
# Classify windows based on their high/low div or background
df.div.stats <- df.div.stats |>
  full_join(tmp) |>
  mutate(window_type = ifelse(is.na(window_type), "Background", window_type)) |>
  pivot_longer(names_to = "statistics", values_to = "values", c(fst, dxy))

df.pi_tajD <- df.pi_tajD |>
  full_join(tmp) |>
  mutate(window_type = ifelse(is.na(window_type), "Background", window_type)) |>
  pivot_longer(names_to = "statistics", values_to = "values", c(pi, TajimaD))

# ------------------------------------------------------------------------------------------------ #
# Explore distribution of datasets
df.div.stats |>
  ggplot(
    aes(
      x = values,
      fill = statistics
    )
  ) +
  geom_histogram() +
  facet_wrap(facets = vars(statistics), scales = "free")

df.pi_tajD |>
  ggplot(
    aes(
      x = values,
      fill = statistics
    )
  ) +
  geom_histogram() +
  facet_wrap(facets = vars(statistics), scales = "free")

# ------------------------------------------------------------------------------------------------ #
# plot the statistics
labeller_between <- as_labeller(
  c(dxy = "D[XY]",fst = "F[ST]"),
  default = label_parsed
)

labeller_within <- as_labeller(
  c(pi = "\U03C0",TajimaD = "Tajimas~D"),
  default = label_parsed
)

plot.between <- df.div.stats |>
  ggplot(
    aes(
      x = comparison,
      y = values,
      fill = window_type
    )
  ) +
  labs(
    fill = "Window type"
  ) +
  geom_boxplot() +
  scale_fill_manual(values = c("#edae49", "#d1495b", "#00798c")) +
  facet_wrap(
    facets = vars(statistics),
    scales = "free_y",
    labeller = labeller(statistics = labeller_between, value = label_value)
  ) +
  guides(fill = guide_legend(title.position = "bottom", title.hjust = 0.5)) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 16)
  )

plot.within <- df.pi_tajD |>
  ggplot(
    aes(
      x = comparison,
      y = values,
      fill = window_type
    )
  ) +
  geom_boxplot() +
  labs(
    fill = "Window type"
  ) +
  scale_fill_manual(values = c("#edae49", "#d1495b", "#00798c")) +
  facet_wrap(facets = vars(statistics), scales = "free_y", labeller = labeller(statistics = labeller_within, value = label_value)) +
  guides(fill = guide_legend(title.position = "bottom", title.hjust = 0.5)) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 16),
    strip.text = element_text(size = 16)
  )


ragg::agg_png(
  filename = here("manuscript", "figures", "figure-x-div_statitics_comp-high_low_background.png"),
  width = 3000,
  height = 3000,
  units = "px",
  res = 300
)
cairo_pdf(
  filename = here("manuscript", "figures", "figure-x-div_statitics_comp-high_low_background.pdf"),
  width = 14,
  height= 14
)
(plot.between / plot.within) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )
invisible(dev.off())

# # ------------------------------------------------------------------------------------------------ #
# # Statitical test
# df.statTest.fst <- df.div.stats |> filter(statistics == "fst")
# df.statTest.dxy <- df.div.stats |> filter(statistics == "dxy")
# df.statTest.pi <- df.pi_tajD |> filter(statistics == "pi")
# df.statTest.tajD <- df.pi_tajD |> filter(statistics == "TajimaD")
#
# # Shapiro Wilk test for normality - significant indicates distributions significantly different to bell
# shapiro.test(sample(df.statTest.fst$values, size = 5000, replace = FALSE))
# shapiro.test(sample(df.statTest.dxy$values, size = 5000, replace = FALSE))
# shapiro.test(sample(df.statTest.pi$values, size = 5000, replace = FALSE))
# shapiro.test(sample(df.statTest.tajD$values, size = 5000, replace = FALSE))
