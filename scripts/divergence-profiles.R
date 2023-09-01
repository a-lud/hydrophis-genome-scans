# ------------------------------------------------------------------------------------------------ #
# Classification of windows
#
# This script is responsible for classifying genomic windows based on Fst and Dxy values. We use
# divergence profiles from Hahn et al. 2017 and Irwin et al. 2018 (and Sheng et al. 2023).

# Libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

# ------------------------------------------------------------------------------------------------ #
# Load data

# Pairwise divergence statistics
df.fst_dxy <- read_csv(
  file = here("results", "between-species", "fst-zfst-dxy.csv"),
  col_names = TRUE,
  col_types = cols()
) |>
  mutate(chromosome = factor(chromosome, paste0(rep("chr", 15), 1:15)))

# ------------------------------------------------------------------------------------------------ #
# Classify windows
df.selective.regime <- df.fst_dxy |>
  group_by(comparison) |>
  mutate(
    top1_dxy = quantile(dxy, 0.99, na.rm = TRUE),
    bottom1_dxy = quantile(dxy, 0.01, na.rm = TRUE),
    top1_fst = quantile(fst, 0.99, na.rm = TRUE),
    bottom1_fst = quantile(fst, 0.01, na.rm = TRUE),
    # Classify rows by adaptation process
    profile = case_when(
      (dxy >= top1_dxy & fst >= top1_fst) ~ "Barrier",
      (dxy >= bottom1_dxy & dxy <=top1_dxy & fst >= top1_fst) ~ "Allopatric",
      (dxy <= bottom1_dxy & fst >= top1_fst) ~ "Recurrent",
      (dxy >= top1_dxy & fst <= bottom1_fst) ~ "Balancing",
      .default = NA_character_
    ),
    profile = factor(profile, levels = c("Allopatric", "Recurrent", "Balancing", "Barrier"))
  ) |>
  ungroup() |>
  filter(!is.na(profile))

# Total number of windows classified
# comparison     count
# 1 HCU-HEL      410
# 2 HCU-HMA      383
# 3 HEL-HMA      429
df.selective.regime |>
  group_by(comparison) |>
  summarise(count = n()) |>
  write_csv(
    file = here("results", "tables", "classified-windows-within.csv"),
    col_names = TRUE
  )

# comparison Allopatric Recurrent Balancing Barrier
# HCU-HEL    307        20        75        8
# HCU-HMA    298        27        50        8
# HEL-HMA    310        22        94        3
df.selective.regime |>
  group_by(comparison, profile) |>
  summarise(count = n()) |>
  pivot_wider(names_from = profile, values_from = count) |>
  select(comparison, Allopatric, Recurrent, Balancing, Barrier) |>
  write_csv(
    file = here("results", "tables", "classified-windows-by-profile-within.csv"),
    col_names = TRUE
  )

# ------------------------------------------------------------------------------------------------ #
# Write results to table
df.selective.regime |>
  write_tsv(
    file = here("results", "tables", "classified-windows.tsv"),
    col_names = TRUE
  )

# ------------------------------------------------------------------------------------------------ #
# Visualise classified windows via scatter plot
plot.scatter_profile <- df.fst_dxy |>
  ggplot(
    aes(
      x = fst,
      y = dxy
    )
  ) +
  labs(x = bquote(bold(F[st])), y = bquote(bold(D[xy]))) +
  # All points
  geom_point(size = 0.2, alpha = 0.1) +
  # Highlight points in each selective regime
  geom_point(
    data = df.selective.regime,
    aes(x = fst, y = dxy, colour = profile),
    alpha = 0.7,
    size = 1
  ) +
  scale_colour_manual(values = c("#66c2a5", "#ffd92f", "#8da0cb", "#fc8d62")) +
  scale_x_continuous(
    breaks = c(0, 0.5, 1),
    labels = c("0", "0.5", "1")
  ) +
  facet_wrap(facets = vars(comparison), nrow = 1) +
  guides(
    colour = guide_legend(
      title.position = "bottom",
      title.hjust = 0.5,
      override.aes = list(size = 10)
    )
  ) +
  geom_hline(data = df.selective.regime, aes(yintercept = top1_dxy), linetype = "dashed", alpha = 0.6) +
  geom_hline(data = df.selective.regime, aes(yintercept = bottom1_dxy), alpha = 0.6) +
  geom_vline(data = df.selective.regime, aes(xintercept = top1_fst), linetype = "dashed", alpha = 0.6) +
  geom_vline(data = df.selective.regime, aes(xintercept = bottom1_fst), alpha = 0.6) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),

    strip.text = element_text(size = 16, face = "bold"),

    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1, "cm")
  )

ragg::agg_png(
  filename = here("results", "figures", "figure-x-scatter-profile.png"),
  width = 1000,
  height = 500,
  units = "px"
)
plot.scatter_profile
invisible(dev.off())

# ------------------------------------------------------------------------------------------------ #
# Collapse windows within comparisons
gr.profile <- df.selective.regime |>
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# Iterate over each comparison and 'reduce' windows that are adjacent
df.profile.reduced <- split(gr.profile, gr.profile$comparison) |>
  map(\(gr) {
    # Collapse adjacent windows
    gr.tmp <- gr |>
      GenomicRanges::reduce(min.gapwidth = 1L, with.revmap = TRUE)

    # Aggregate metadata columns
    GenomicRanges::mcols(gr.tmp) <- aggregate(
      gr,
      GenomicRanges::mcols(gr.tmp)$revmap,
      comparison = comparison, profile = profile,
      drop = FALSE
    )

    # Turn back to tibble and unnest on comparison column
    gr.tmp |>
      as_tibble() |>
      select(-grouping) |>
      unnest(col = comparison) |>
      distinct()
  }) |>
  list_rbind() |>
  select(comparison, chromosome = seqnames, start, end, width, profile)

# ------------------------------------------------------------------------------------------------ #
# Total number of regions after collapsing adjacent windows
# comparison     count
# 1 HCU-HEL      287
# 2 HCU-HMA      233
# 3 HEL-HMA      247
df.profile.reduced |>
  group_by(comparison) |>
  summarise(count = n()) |>
  write_csv(
    file = here("results", "tables", "classified-windows-reduced-within.csv"),
    col_names = TRUE
  )

# Number of classified windows broken down by profile
# comparison    Allopatric `Allopatric|Recurrent` `Allopatric|Barrier` Recurrent Balancing Barrier
# HCU-HEL           206                      6                   NA        11        56       8
# HCU-HMA           165                     13                   NA         9        38       8
# HEL-HMA           158                     12                    1         9        65       2
df.profile.reduced |>
  unnest(profile) |>
  distinct() |>
  group_by(comparison, chromosome, start, end, width) |>
  mutate(profile = paste0(sort(unique(profile)), collapse = "|")) |>
  distinct() |>
  group_by(comparison, profile) |>
  summarise(count = n()) |>
  pivot_wider(names_from = profile, values_from = count) |>
  select(1, 2, 3, 7, 6, 4, 5) |>
  write_csv(
    file = here("results", "tables", "classified-windows-reduced-by-profile-within.csv"),
    col_names = TRUE,
    na = "0"
  )

# ------------------------------------------------------------------------------------------------ #
# Genome wide view of classified windows
df.profile.plot <- df.selective.regime |>
  mutate(
    chromosome = as.numeric(str_remove(chromosome, "chr")),
    chromosome = factor(chromosome, 1:15)
  )

plot.zfst_classified <- df.fst_dxy |>
  mutate(
    chromosome = as.numeric(str_remove(chromosome, "chr")),
    chromosome = factor(chromosome, 1:15)
  ) |>
  ggplot(
    aes(
      x = (start + end)/2,
      y = zfst
    )
  ) +
  labs(
    x = "\nChromosome (Mbp)",
    y = bquote(ZF[st])
  ) +
  geom_point(size = 0.5, alpha = 0.1, colour = "grey80") +
  geom_point(
    data = df.profile.plot,
    aes(x = (start + end)/2, y = zfst, colour = profile)
  ) +
  # Plotting again to make these points super clear
  geom_point(
    data = df.profile.plot |> filter(profile == "Recurrent"),
    aes(x = (start + end)/2, y = zfst),
    colour = "#ffd92f",
    size = 0.8
  ) +
  scale_colour_manual(values = c("#66c2a5", "#ffd92f", "#8da0cb", "#fc8d62")) +
  scale_x_continuous(
    labels = scales::label_number(scale = 1e-6),
    expand = c(0.03, 0)
  ) +
  guides(
    colour = guide_legend(
      title.hjust = 0.5,
      override.aes = list(size = 10)
    )
  ) +
  facet_grid(
    rows = vars(comparison),
    cols = vars(chromosome),
    scales = "free_x",
    space = "free",
    switch = "x"
  ) +
  theme_classic() +
  theme(
    # Axis
    axis.text = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 14),

    # Legend
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 14, angle = 90),
    strip.placement = "outside",
  )

ragg::agg_png(
  filename = here("results", "figures", "figure-x-genome_wide-zfst_classified.png"),
  width = 4000,
  height = 2000,
  units = "px",
  res = 300
)
plot.zfst_classified
invisible(dev.off())

# ------------------------------------------------------------------------------------------------ #
# Shared loci between comparisons - i.e. genomic windows that are classified in multiple comparisons

# Total number shared - In total there are 745 unique windows, so total should sum to that
# comparison                count
# HCU-HEL                   147
# HCU-HMA                   109
# HEL-HMA                   171
# HCU-HEL|HCU-HMA            60
# HCU-HEL|HEL-HMA            44
# HCU-HMA|HEL-HMA            55
# HCU-HEL|HCU-HMA|HEL-HMA   159
df.selective.regime |>
  group_by(chromosome, start, end) |>
  summarise(
    comparison = paste0(sort(unique(comparison)), collapse = "|"),
    profile = paste0(sort(unique(profile)), collapse = "|")
  ) |>
  ungroup() |>
  group_by(comparison) |>
  summarise(count = n()) |> pull(count) |> sum()
  mutate(
    comparison = factor(comparison, c(
      "HCU-HEL", "HCU-HMA", "HEL-HMA",
      "HCU-HEL|HCU-HMA", "HCU-HEL|HEL-HMA", "HCU-HMA|HEL-HMA",
      "HCU-HEL|HCU-HMA|HEL-HMA"
    ))
  ) |>
  arrange(comparison) |>
  write_csv(
    file = here("results", "tables", "classified-windows-between.csv"),
    col_names = TRUE,
    na = "0"
  )

# Shared loci stratified by profile
# comparison                  Allopatric Recurrent Balancing Barrier `Allopatric|Recurrent` `Allopatric|Barrier`
# HCU-HEL                        107         1        32       7                     NA                   NA
# HCU-HMA                         82         4        16       7                     NA                   NA
# HEL-HMA                        117         1        52       1                     NA                   NA
# HCU-HEL|HCU-HMA                 48         2         7       1                      2                   NA
# HCU-HEL|HEL-HMA                 25        NA        15      NA                      3                    1
# HCU-HMA|HEL-HMA                 44         1         6      NA                      3                    1
# HCU-HEL|HCU-HMA|HEL-HMA        107         7        21      NA                     24                   NA
df.selective.regime |>
  group_by(chromosome, start, end) |>
  summarise(
    comparison = paste0(sort(unique(comparison)), collapse = "|"),
    profile = paste0(sort(unique(profile)), collapse = "|")
  ) |>
  ungroup() |>
  group_by(comparison, profile) |>
  summarise(count = n()) |>
  mutate(
    comparison = factor(comparison, c(
      "HCU-HEL", "HCU-HMA", "HEL-HMA",
      "HCU-HEL|HCU-HMA", "HCU-HEL|HEL-HMA", "HCU-HMA|HEL-HMA",
      "HCU-HEL|HCU-HMA|HEL-HMA"
    )),
    profile = factor(profile, levels = c(
      "Allopatric", "Recurrent", "Balancing", "Barrier",
      "Allopatric|Recurrent", "Allopatric|Barrier"
    ))
  ) |>
  arrange(profile, comparison) |>
  pivot_wider(names_from = profile, values_from = count) |>
  write_csv(
    file = here("results", "tables", "classified-windows-by-profile-between.csv"),
    col_names = TRUE,
    na = "0"
  )

# ------------------------------------------------------------------------------------------------ #
# Merge adjacent windows into single large "regions"

# Collapse adjacent windows
gr.tmp <- gr.profile |>
  GenomicRanges::reduce(min.gapwidth = 1L, with.revmap = TRUE)

# Aggregate metadata columns
GenomicRanges::mcols(gr.tmp) <- aggregate(
  gr.profile,
  GenomicRanges::mcols(gr.tmp)$revmap,
  comparison = comparison, profile = profile,
  drop = FALSE
)

# Merge
gr.tmp |>
  as_tibble() |>
  select(-grouping, -strand) |>
  group_by(seqnames, start, width) |>
  mutate(
    profile = paste0(sort(unique(unlist(profile))), collapse = "|"),
    comparison = paste0(sort(unique(unlist(comparison))), collapse = "|"),
  ) |>
  distinct() |>
  group_by(comparison) |>
  summarise(sum = n()) |>
  mutate(
    comparison = factor(comparison, c(
      "HCU-HEL", "HCU-HMA", "HEL-HMA",
      "HCU-HEL|HCU-HMA", "HCU-HEL|HEL-HMA", "HCU-HMA|HEL-HMA",
      "HCU-HEL|HCU-HMA|HEL-HMA"
    ))
  ) |>
  arrange(comparison) |>
  write_csv(
    file = here("results", "tables", "classified-windows-by-profile-reduced-between.csv"),
    col_names = TRUE,
    na = "0"
  )

# ------------------------------------------------------------------------------------------------ #
# Pair wise Wilcoxon test - shared genomic regions larger than comparison-specific genomic regions?
dat <- gr.tmp |>
  as_tibble() |>
  select(-grouping) |>
  group_by(seqnames, start, width) |>
  mutate(
    profile = paste0(sort(unique(unlist(profile))), collapse = "|"),
    comparison = paste0(sort(unique(unlist(comparison))), collapse = "|"),
  ) |>
  distinct() |>
  mutate(
    comparison = factor(comparison, c(
      "HCU-HEL", "HCU-HMA", "HEL-HMA",
      "HCU-HEL|HCU-HMA", "HCU-HEL|HEL-HMA", "HCU-HMA|HEL-HMA",
      "HCU-HEL|HCU-HMA|HEL-HMA"
    ))
  )

wilcox.out <- pairwise.wilcox.test(
  x = dat$width,
  g = dat$comparison,
  p.adjust.method = "bonferroni",
  exact = FALSE
)

wilcox.out$p.value |>
  as_tibble(rownames = "comparison") |>
  write_csv(
    file = here("results", "tables", "classified-windows-size-comparison-wilcox-results.csv"),
    col_names = TRUE,
    na = ""
  )

# ------------------------------------------------------------------------------------------------ #
# Plot of "unique" loci - checking if they're genuinely alone or nearby to other peaks
df.profile.unique <- df.selective.regime |>
  group_by(chromosome, start, end) |>
  summarise(
    comparison = paste0(sort(unique(comparison)), collapse = "|"),
    profile = paste0(sort(unique(profile)), collapse = "|")
  ) |>
  filter(comparison %in% c("HCU-HEL", "HCU-HMA", "HEL-HMA")) |>
  left_join(df.selective.regime) |>
  mutate(
    chromosome = as.numeric(str_remove(chromosome, "chr")),
    chromosome = factor(chromosome, 1:15)
  )

plot.shared_unique <- df.fst_dxy |>
  mutate(
    chromosome = as.numeric(str_remove(chromosome, "chr")),
    chromosome = factor(chromosome, 1:15)
  ) |>
  ggplot(
    aes(
      x = (start + end)/2,
      y = zfst
    )
  ) +
  labs(x = "\nChromosome (Mbp)", y = bquote(ZF[st])) +
  geom_point(size = 0.5, alpha = 0.1, colour = "grey80") +
  geom_point(data = df.profile.plot, aes(x = (start + end)/2, y = zfst, colour = profile)) +
  # Plotting again to make these points stand out
  geom_point(
    data = df.profile.plot |> filter(profile == "Recurrent"),
    aes(x = (start + end)/2, y = zfst),
    colour = "#ffd92f"
  ) +
  # Points that are "unique"
  geom_point(
    data = df.profile.unique,
    mapping = aes(x = (start + end)/2, y = zfst),
    colour = "red",
    alpha = 0.5
  ) +
  scale_colour_manual(values = c("#66c2a5", "#ffd92f", "#8da0cb", "#fc8d62")) +
  scale_x_continuous(labels = scales::label_number(scale = 1e-6),expand = c(0.03, 0)) +
  guides(
    colour = guide_legend(
      title.hjust = 0.5,
      override.aes = list(size = 10)
    )
  ) +
  facet_grid(
    rows = vars(comparison),
    cols = vars(chromosome),
    scales = "free_x",
    space = "free",
    switch = "x"
  ) +
  theme_classic() +
  theme(
    # Axis
    axis.text = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 14),

    # Legend
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 14, angle = 90),
    strip.placement = "outside"
  )

ragg::agg_png(
  filename = here("results", "figures", "figure-x-genome_wide-zfst_classified-shared-unique.png"),
  width = 4000,
  height = 2000,
  units = "px",
  res = 300
)
plot.shared_unique
invisible(dev.off())
