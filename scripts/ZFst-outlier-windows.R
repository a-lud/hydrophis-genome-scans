# ------------------------------------------------------------------------------------------------ #
# Outlier windows
#
# Here we use ZFst to identify outlier windows in each species comparison. An outlier window is where
# the Z-score (i.e. Z-transformed Fst valuye) is +/-3. We extract these windows, merge windows that
# are next to each other, then check the following:
#
# - Number of 50kbp outlier windows per comparison (before/after merging)
# - Size of outlier windows

# Libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

# ------------------------------------------------------------------------------------------------ #
# Load divergence data
dxy <- fs::dir_ls(path = here("results", "between-species"), glob = "*dxy.txt") %>%
  magrittr::extract(str_detect(., "50kb")) |>
  map(read_tsv, col_names = TRUE, col_types = cols()) |>
  list_rbind() |>
  filter(no_sites >= 5000) |>
  unite(col = "comparison", sep = "-", pop1, pop2) |>
  select(comparison, chromosome, start = window_pos_1, end = window_pos_2, dxy=avg_dxy)

# Between species - Fst (weighted average per window)
fst <- fs::dir_ls(path = here("results", "between-species"), glob = "*fst.txt") %>%
  magrittr::extract(str_detect(., "50kb")) |>
  map(read_tsv, col_names = TRUE, col_types = cols()) |>
  list_rbind() |>
  unite(col = "comparison", sep = "-", pop1, pop2) |>
  select(comparison, chromosome, start = window_pos_1, end = window_pos_2, fst = avg_wc_fst)

# Nucleotide diversity (pi)
nuc.div <- fs::dir_ls(path = here("results", "within-species"), glob = "*pi.txt") %>%
  magrittr::extract(str_detect(., "50kb")) %>%
  magrittr::extract(str_detect(., "HMA|HEL|HCU")) |>
  map(read_tsv, col_names = TRUE, col_types = cols()) |>
  list_rbind() |>
  filter(no_sites >= 5000) |>
  rename(comparison = pop, pi = avg_pi) |>
  rename(start = window_pos_1, end = window_pos_2) |>
  select(-c(no_sites, count_diffs, count_comparisons, count_missing))

# Tajima's D
tajD <- fs::dir_ls(path = here("results", "within-species"), glob = "*TajimaD.tsv") %>%
  magrittr::extract(str_detect(., "50kb")) %>%
  magrittr::extract(str_detect(., "HMA|HEL|HCU")) |>
  read_tsv(col_names = TRUE, col_types = cols(), id = "comparison") |>
  mutate(
    comparison = sub("\\..*", "", basename(comparison)),
    end = BIN_START + 5e4,
    start = BIN_START + 1
  ) |>
  select(comparison, chromosome = CHROM, start, end, TajimaD)

fst.dxy.avg <- read_csv(
  here("results", "tables", "genome-wide-dxy_fst.csv"),
  col_names = TRUE,
  col_types = cols()
) |>
  select(comparison = Comparison, avg_fst = `Fst (weighted)`)

# Z-normalised Fst - approach taken from Han et al. 2017
# www.genome.org/cgi/doi/10.1101/gr.212522.116

# Calculate SD - use same formula as sd() function = sqrt( (sum(x - u)^2)/N-1 )
# - x = avg Fst for each window (rows in dataframe)
# - u = weighted genome-wide avg. Fst - obtained from VCFtools (reported at the end of a run)
# - N-1 = number of windows for each comparison - 1
fst.sd <- left_join(fst, fst.dxy.avg) |>
  mutate(sq_deviation = (fst  - avg_fst)^2) |>
  group_by(comparison, avg_fst) |>
  summarise(Fst_sd = sqrt(x = sum(sq_deviation, na.rm = TRUE)/(n() - 1)))

fst <- reduce(list(fst, fst.sd), left_join) |>
  # Calculate Z(Fst) using equation from Han. et al. 2017 above
  mutate(zfst = (fst - avg_fst)/Fst_sd) |>
  select(-c(avg_fst, Fst_sd))

# 100,107 windows
df.fst_dxy <- fst |>
  inner_join(dxy) |>
  mutate(chromosome = factor(chromosome, paste0(rep("chr", 15), 1:16)))

# Write object to file to save remaking this thing constantly
df.fst_dxy |>
  write_csv(
    file = here("results", "between-species","fst-zfst-dxy.csv"),
    col_names = TRUE
  )

# ------------------------------------------------------------------------------------------------ #
# Outlier loci based on Zfst +/- 3 - number of windows before merging adjacent

# comparison  ZFst type          count
# HCU-HEL     +3                 8
# HCU-HEL     -3                 6
# HCU-HMA     +3                 2
# HCU-HMA     -3                26
# HEL-HMA     -3                96

df.fst_dxy |>
  filter(zfst >= 3 | zfst <= -3) |>
  mutate(classification = ifelse(zfst >= 3, "+3", "-3")) |>
  group_by(comparison, classification) |>
  summarise(count = n())

# Write to file
df.fst_dxy |>
  filter(zfst >= 3 | zfst <= -3) |>
  mutate(classification = ifelse(zfst >= 3, "+3", "-3")) |>
  write_csv(
    file = here("results", "between-species", "zfst-outliers.csv"),
    col_names = TRUE
  )

# Write as BED file for intersecting
df.fst_dxy |>
  filter(zfst >= 3 | zfst <= -3) |>
  mutate(classification = ifelse(zfst >= 3, "+3", "-3")) |>
  relocate(chromosome, start, end, comparison, everything()) |>
  write_tsv(
    file = here("results", "between-species", "zfst-outliers.bed"),
    col_names = FALSE
  )

## Low ZFst overwhelmingly driven by the start of chromsome 2
# Number of Zfst <= -3: 128
zfst.lowlier <- df.fst_dxy |>
  filter(zfst <= -3)

# 94 windows in the first 70Mb of chromosome 2
zfst.lowlier |>
  filter(chromosome == "chr2", start <= 70e6)

zfst.lowlier |>
  filter(chromosome == "chr2" & start <= 70e6) |>
  group_by(comparison) |>
  summarise(count = n())

# ------------------------------------------------------------------------------------------------ #
# Collapse windows within 50kbp of each other - within each species
# Helpful: https://bioconductor.org/packages/release/bioc/vignettes/HelloRanges/inst/doc/tutorial.pdf

# Turn dataframe into GRanges
gr.fst_dxy <- df.fst_dxy |>
  filter(zfst >= 3 | zfst <= -3) |>
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# Iterate over each comparison and 'reduce' windows that are adjacent
df.fst_dxy.reduce <- split(gr.fst_dxy, gr.fst_dxy$comparison) |>
  map(\(gr) {
    # Collapse adjacent windows
    gr.tmp <- gr |>
      GenomicRanges::reduce(min.gapwidth = 2L, with.revmap = TRUE)

    # Aggregate metadata columns
    GenomicRanges::mcols(gr.tmp) <- aggregate(
      gr,
      GenomicRanges::mcols(gr.tmp)$revmap,
      zfst = zfst, comparison = comparison,
      drop = FALSE
    )

    # Turn back to tibble and unnest on comparison column - leave ZFst values as vector (see below)
    gr.tmp |>
      as_tibble() |>
      select(-grouping) |>
      unnest(col = comparison) |>
      distinct() |>
      group_by(start)
  }) |>
  list_rbind() |>
  select(comparison, chromosome = seqnames, start, end, width, zfst)

# Summarise number of windows with different widths
# comparison  `Window type`       `50000` `100000` `150000` `200000` `250000`
# 1 HCU-HEL    +3                  8       NA       NA       NA       NA
# 2 HCU-HEL    -3                  3       NA        1       NA       NA
# 3 HCU-HMA    -3                 15        2        1        1       NA
# 4 HCU-HMA    +3                 NA        1       NA       NA       NA
# 5 HEL-HMA    -3                 53       11        2       NA        3
df.fst_dxy.reduce.summary <- df.fst_dxy.reduce |>
  mutate(`Window type` = ifelse(all(unlist(zfst) < 0), "-3", "+3")) |>
  group_by(comparison, width, `Window type`) |>
  summarise(count = n()) |>
  pivot_wider(names_from = width, values_from = count) |>
  select(comparison, `Window type`, `50000`, `100000`, `150000`, `200000`, `250000`)

df.fst_dxy.reduce.summary |>
  write_csv(
    file = here("results", "tables", "zfst-outlier_windows-reduced-within.csv"),
    col_names = TRUE,
    na = ""
  )

# ------------------------------------------------------------------------------------------------ #
# Find shared outlier windows across species comparisons
gr.fst_dxy.reduced <- gr.fst_dxy |>
  GenomicRanges::reduce(min.gapwidth = 2, with.revmap = TRUE)

GenomicRanges::mcols(gr.fst_dxy.reduced) <- aggregate(
  gr.fst_dxy,
  GenomicRanges::mcols(gr.fst_dxy.reduced)$revmap,
  zfst = zfst, comparison = comparison,
  drop = FALSE
)

# comparison                `Window type`     `50000` `100000` `150000` `250000` `200000`
# HCU-HEL                       +3             7       NA       NA       NA       NA
# HCU-HMA                       -3             5        1       NA       NA       NA
# HEL-HMA                       -3             44        9        1        1       NA
# HCU-HEL|HCU-HMA               +3             NA        1       NA       NA       NA
# HCU-HMA|HEL-HMA               -3             6        1       NA        1        1
# HCU-HEL|HCU-HMA|HEL-HMA       -3             1       NA        2        1       NA
df.shared.outlier.between <- gr.fst_dxy.reduced |>
  as_tibble() |>
  select(-grouping) |>
  unnest(col = comparison) |>
  group_by(seqnames, start, width) |>
  mutate(
    comparison = paste0(sort(unique(comparison)), collapse = "|"),
    `Window type` = ifelse(all(unlist(zfst) < 0), "-3", "+3")
  ) |>
  ungroup() |>
  distinct() |>
  # select(comparison, chromosome = seqnames, width, `Window type`) |>
  group_by(comparison, width, `Window type`) |>
  summarise(count = n()) |>
  pivot_wider(names_from = width, values_from = count) |>
  mutate(
    comparison = factor(
      comparison, c(
        "HCU-HEL", "HCU-HMA", "HEL-HMA",
        "HCU-HEL|HCU-HMA", "HCU-HMA|HEL-HMA",
        "HCU-HEL|HCU-HMA|HEL-HMA"
      )
    )
  ) |>
  arrange(comparison)

df.shared.outlier.between |>
  write_csv(
    file = here("results/tables/zfst-outlier_windows-between.csv"),
    col_names = TRUE,
    na = ""
  )

gr.fst_dxy.reduced |>
  plyranges::filter(start >= 214550001)

# Getting range of ZFst for windows in shared outlier region
df.fst_dxy |>
  filter(chromosome == "chr1" & start >= 214550001 & end <= 214650001)

# ------------------------------------------------------------------------------------------------ #
# Outliers defined by 3 (+/-) Zfst
df.plt <- df.fst_dxy |>
  mutate(
    chromosome = str_remove(chromosome, "chr"),
    chromosome = as.numeric(chromosome),
    grouping = ifelse(chromosome %% 2 != 0, "even", "odd"),
    chromosome = factor(chromosome, levels = 1:15)
  )

plot.zfst <- df.plt |>
  ggplot(
    aes(
      x = start,
      y = zfst,
      colour = grouping
    )
  ) +
  labs(x = "\nChromosome (Mbp)", y = bquote(bold(ZF[st]))) +
  geom_hline(yintercept = c(-3, 3), linetype = "dashed", colour = "grey40", alpha = 0.5) +
  geom_point(size = 0.5, alpha = 0.1) +
  geom_point(data = df.plt |> filter(zfst >= 3), mapping = aes(x = start, y = zfst), colour = "#F24C3D", size = 0.8) +
  geom_point(data = df.plt |> filter(zfst <= -3), mapping = aes(x = start, y = zfst), colour = "#22A699", size = 0.8) +
  scale_colour_manual(values = c("grey60", "grey80")) +
  scale_y_continuous(
    limits = c(-5, 4),
    breaks = seq(-4, 4, 2),
    labels = c("-4", "-2", "0", "2", "4")
  ) +
  scale_x_continuous(
    labels = scales::label_number(scale = 1e-6),
    expand = c(0.03, 0)
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
    # Axis text/labels
    axis.text = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 16),

    # Strip text and orientation
    strip.text = element_text(size = 16),
    strip.text.x = element_text(size = 14, hjust = 0.9, vjust = 0.8),
    # strip.text.y = element_text(angle = 0),
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.spacing = unit(0.05, "cm"),

    # Legend
    legend.position = "none"
  )

ragg::agg_png(
  filename = here("results", "figures", "figure-x-genome_wide-zfst.png"),
  width = 1400,
  height = 600,
  units = "px"
)
plot.zfst
invisible(dev.off())

# ------------------------------------------------------------------------------------------------ #
# Compare Fst, Dxy, Pi and TajD in outlier windows relative to genomic Background
strip_labeller <- as_labeller(
  c(
    pi = "Avg.~\U03C0",
    dxy = "D[XY]",
    zfst = "ZF[ST]"
  ),
  default = label_parsed
)

avg.pi <- nuc.div |>
  pivot_wider(names_from = comparison, values_from = pi) |>
  mutate(
    `HCU-HEL` = (HCU + HEL)/2,
    `HCU-HMA` = (HCU + HMA)/2,
    `HEL-HMA` = (HEL + HMA)/2
  ) |>
  select(-c(HCU, HMA, HEL)) |>
  pivot_longer(values_to = "pi", names_to = "comparison", 4:6)

plot.outlier.statistics <- df.fst_dxy |>
  inner_join(avg.pi) |>
  mutate(grouping = case_when(
    zfst <= -3 ~ "ZFst (-3)",
    zfst >= 3 ~ "ZFst (+3)",
    .default = "Genomic\nbackground"
  )) |>
  select(comparison, chromosome, start, end, zfst, dxy, pi, grouping) |>
  pivot_longer(names_to = "statistic", values_to = "values", c(zfst, dxy, pi)) |>
  ggplot(
    aes(
      x = grouping,
      y = values,
      fill = comparison
    )
  ) +
  geom_boxplot() +
  facet_wrap(
    facets = vars(statistic),
    scales = "free_y",
    labeller = labeller(statistic = strip_labeller, value = label_value)
  ) +
  scale_fill_manual(values = c("#edae49", "#d1495b", "#00798c")) +
  guides(fill = guide_legend(label.position = "bottom")) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = 14),

    strip.text = element_text(size = 14),

    # legend.position = "bottom",
    # legend.key.size = unit(1, "cm"),
    legend.title = element_blank(),
    legend.text = element_text(size = 14)
  )

ragg::agg_png(
  filename = here("results", "figures", "figure-x-zfst-outliers-compare-statistics.png"),
  width = 800,
  height = 300,
  units = "px"
)
plot.outlier.statistics
invisible(dev.off())

