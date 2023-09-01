# ------------------------------------------------------------------------------------------------ #
# Profiled regions: Allopatric regions
#
# Here I manually visualise "Allopatric" loci. These are the main, obvious peaks of divergence
# across the genome, and consequently appear to be of significant interest. I manually plot each
# region and assess which loci appear interesting.

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(patchwork)
})

source(here("scripts", "accessory-plotting-scripts.R"))

# ------------------------------------------------------------------------------------------------ #
# Load data

df.div.stats <- read_csv(
  file = here("results", "between-species", "fst-zfst-dxy.csv"),
  col_names = TRUE,
  col_types = cols()
) |>
  pivot_longer(names_to = "statistics", values_to = "values", cols = c(zfst, dxy)) |>
  mutate(midpoint = floor((start + end)/2))

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

df.pi_tajD <- inner_join(df.pi, df.tajD) |>
  pivot_longer(names_to = "statistics", values_to = "values", cols = c(pi, TajimaD)) |>
  mutate(
    midpoint = floor((start + end)/2),
    statistics = factor(statistics, levels = c("pi", "TajimaD")),
    comparison = factor(comparison, levels = c("HCU", "HMA", "HEL"))
  )

# Introgression statistics (sliding window fd)
#   -  Filter sites wehere D < 0.
df.fd <- fs::dir_ls(
  path = here("results/between-species/genomics_general"),
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
  filter(sites >= 100) |>
  select(comparison, chromosome, start, end, values = fd) |>
  mutate(statistics = "fd", midpoint = floor((start + end)/2))

# Join fd with div.stats
df.div.stats <- list_rbind(list(df.div.stats, df.fd)) |>
  select(-fst) |>
  mutate(
    statistics = factor(statistics, levels = c("fd", "zfst", "dxy")),
    comparison = factor(comparison, c("HCU-HEL", "HCU-HMA", "HEL-HMA"))
  )

# Join fd with classified data
df.classified <- read_tsv(
  file = here("results", "tables", "classified-windows.tsv"),
  col_names = TRUE,
  col_types = cols()
) |>
  pivot_longer(names_to = "statistics", values_to = "values", cols = c(zfst, dxy)) |>
  select(-starts_with("top"), -starts_with("bottom"), -fst) |>
  mutate(midpoint = floor((start + end)/2))

df.tmp <- df.fd |>
  inner_join(
    df.classified |>
      select(comparison, chromosome, start, end, profile) |>
      distinct()
  )

# long form table with all statistics in it
tmp <- bind_rows(df.classified, df.tmp) |>
  arrange(chromosome, start, comparison) |>
  mutate(
    statistics = factor(statistics, levels = c("fd", "zfst", "dxy")),
    comparison = factor(comparison, c("HCU-HEL", "HCU-HMA", "HEL-HMA"))
  )

# ------------------------------------------------------------------------------------------------ #
# Feature data
gff <- read_tsv(
  file = here("data", "hydrophis_major.gff3"),
  col_names = c(
    "chromosome", "method", "feature",
    "start", "end", "score",
    "strand", "phase", "attributes"
  ),
  col_types = cols(),
  comment = "#"
) |>
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

repeats <- read_tsv(
  file = here("data", "hydrophis_major.repeats.gff3"),
  col_names = c(
    "chromosome", "method", "feature",
    "start", "end", "score",
    "strand", "phase", "attributes"
  ),
  col_types = cols(),
  comment = "#"
) |>
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# Make genome GRange
gr.genome <- read_tsv(
  file = here("data", "hydrophis_major.fa.fai"),
  col_names = c("chromosome", "end"),
  col_types = cols(),
  col_select = c("chromosome", "end")
) |>
  mutate(start = 1) |>
  filter(str_detect(chromosome, "chr") & chromosome != "chrZ") |>
  GenomicRanges::makeGRangesFromDataFrame()

# Tile genome
gr.genome.seqinfo <- GenomeInfoDb::Seqinfo(seqnames = levels(gr.genome@seqnames), seqlengths = gr.genome@ranges@width)
gr.genome.tiles <- unlist(GenomicRanges::tileGenome(seqlengths = gr.genome.seqinfo, tilewidth = 5e4))

# Count number of features overlapping each 50kb window
gr.genome.tiles$count <- GenomicRanges::countOverlaps(query = gr.genome.tiles, repeats)

# Using the repeat counts per window for plotting rather than plotting actual repetitive elements
list.features <- list("Genes" = gff, "Repeats" = gr.genome.tiles)

# ------------------------------------------------------------------------------------------------ #
# Allopatric loci

# Already accounted for in ZFst
ragg::agg_png(
  filename = here("results/figures/high-div/figure-x-div_profile-chr1-212_214.png"),
  width = 3000,
  height = 3000,
  units = "px",
  res = 300
)
plotRegions(
  df_between = df.div.stats,
  df_within = df.pi_tajD,
  df_outliers_between = tmp |> filter(chromosome == "chr1" & start >= 211e6 & end <= 216e6),
  chr = "chr1",
  s = 211000000 - 10e6,
  e = 215700000 + 10e6,
  by = 2e6,
  pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
  gr_features = list.features,
  span = 0.5
)
invisible(dev.off())

# All statistics are elevated...
ragg::agg_png(
  filename = here("results/figures/high-div/figure-x-div_profile-chr1-400mb.png"),
  width = 3000,
  height = 3000,
  units = "px",
  res = 300
)
plotRegions(
  df_between = df.div.stats,
  df_within = df.pi_tajD,
  df_outliers_between = tmp |> filter(chromosome == "chr1" & start >= 398150001 & end <= 401850000),
  chr = "chr1",
  s = 398150001 - 5e6,
  e = 401850000 + 5e6,
  by = 1e6,
  pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
  gr_features = list.features,
  span = 0.4
)
invisible(dev.off())

# Same as ZFst region
ragg::agg_png(
  filename = here("results/figures/high-div/figure-x-div_profile-chr2-240mb.png"),
  width = 3000,
  height = 3000,
  units = "px",
  res = 300
)
plotRegions(
  df_between = df.div.stats,
  df_within = df.pi_tajD,
  df_outliers_between = tmp |> filter(chromosome == "chr2" & start >= 237300001 & end <= 242650001),
  chr = "chr2",
  s = 237300001 - 5e6,
  e = 242650001 + 5e6,
  by = 2e6,
  pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
  gr_features = list.features,
  span = 0.4
)
invisible(dev.off())

# Chromosome 3: already in ZFst output - has a lowlier right before the peak
ragg::agg_png(
  filename = here("results/figures/high-div/figure-x-div_profile-chr3-118_127mb.png"),
  width = 3000,
  height = 3000,
  units = "px",
  res = 300
)
plotRegions(
  df_between = df.div.stats,
  df_within = df.pi_tajD,
  df_outliers_between = tmp |> filter(chromosome == "chr3" & start >= 119650001 & end <= 126150000),
  chr = "chr3",
  s = 119650001 - 5e6,
  e = 126150000 + 5e6,
  by = 2e6,
  pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
  gr_features = list.features,
  span = 0.4
)
invisible(dev.off())

# This one is interesting - has the drop ZFst right in the middle of a peak in HCU-HEL
ragg::agg_png(
  filename = here("results/figures/high-div/figure-x-div_profile-chr3-137_147mb.png"),
  width = 3000,
  height = 3000,
  units = "px",
  res = 300
)
plotRegions(
  df_between = df.div.stats,
  df_within = df.pi_tajD,
  df_outliers_between = tmp |> filter(chromosome == "chr3" & start >= 141850001 & end <= 146250000),
  chr = "chr3",
  s = 141850001 - 10e6,
  e = 146250000 + 10e6,
  by = 2e6,
  pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
  gr_features = list.features,
  span = 0.4
)
invisible(dev.off())

# Chromosome 4 - Massively increased TajD in H. elegans in this comparison
ragg::agg_png(
  filename = here("results/figures/high-div/figure-x-div_profile-chr4-112_123mb.png"),
  width = 3000,
  height = 3000,
  units = "px",
  res = 300
)
plotRegions(
  df_between = df.div.stats,
  df_within = df.pi_tajD,
  df_outliers_between = tmp |> filter(chromosome == "chr4" & start >= 112350001 & end <= 123650000),
  chr = "chr4",
  s = 112350001 - 5e6,
  e = 123650000 + 5e6,
  by = 2e6,
  pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
  gr_features = list.features,
  span = 0.4
)
invisible(dev.off())

# Similar story as chromosome 3 - Large region of elevated ZFst but a small window of decreased diversity
# at the centre
ragg::agg_png(
  filename = here("results/figures/high-div/figure-x-div_profile-chr5-60_74mb.png"),
  width = 3000,
  height = 3000,
  units = "px",
  res = 300
)
plotRegions(
  df_between = df.div.stats,
  df_within = df.pi_tajD,
  df_outliers_between = tmp |> filter(chromosome == "chr5" & start >= 61350001 & end <= 73200000),
  chr = "chr5",
  s = 61350001 - 10e6,
  e = 73200000 + 10e6,
  by = 2e6,
  pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
  gr_features = list.features,
  span = 0.3
)
invisible(dev.off())

# Region already reported in ZFst approach
ragg::agg_png(
  filename = here("results/figures/high-div/figure-x-div_profile-chr6-92mb.png"),
  width = 3000,
  height = 3000,
  units = "px",
  res = 300
)
plotRegions(
  df_between = df.div.stats,
  df_within = df.pi_tajD,
  df_outliers_between = tmp |> filter(chromosome == "chr6" & start >= 91150001 & end <= 94750000),
  chr = "chr6",
  s = 91150001 - 4e6,
  e = 94750000 + 4e6,
  by = 1e6,
  pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
  gr_features = list.features,
  span = 0.5
)
invisible(dev.off())

# ------------------------------------------------------------------------------------------------ #
# Balanced loci (NOT in chromosome 2 jumbled area)
tmp |> filter(profile == "Balancing") |> View()

# Chromosome 1: Found in ZFst scan
ragg::agg_png(
  filename = "results/figures/low-div/figure-x-div_profile-chr1-1_2mb.png",
  width = 3000,
  height = 3000,
  units = "px",
  res = 300
)
plotRegions(
  df_between = df.div.stats,
  df_within = df.pi_tajD,
  df_outliers_between = tmp |> filter(chromosome == "chr1" & start >= 950001 & end <= 1950000),
  chr = "chr1",
  s = 1,
  e = 1950000 + 1e6,
  by = 5e5,
  pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
  gr_features = list.features,
  span = 0.5
)
invisible(dev.off())

# Chromosome 1: Found in ZFst scan
ragg::agg_png(
  filename = here("results/figures/low-div/figure-x-div_profile-chr1-134_135mb.png"),
  width = 3000,
  height = 3000,
  units = "px",
  res = 300
)
plotRegions(
  df_between = df.div.stats,
  df_within = df.pi_tajD,
  df_outliers_between = tmp |> filter(chromosome == "chr1" & start >= 134700001 & end <= 135050000),
  chr = "chr1",
  s = 134700001 - 5e6,
  e = 135050000 + 5e6,
  by = 2e6,
  pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
  gr_features = list.features,
  span = 0.5
)
invisible(dev.off())

# Chromosome 1: Centre of "Allopatric" peak reported above
ragg::agg_png(
  filename = here("results/figures/low-div/figure-x-div_profile-chr1-399mb.png"),
  width = 3000,
  height = 3000,
  units = "px",
  res = 300
)
plotRegions(
  df_between = df.div.stats,
  df_within = df.pi_tajD,
  df_outliers_between = tmp |> filter(chromosome == "chr1" & start >= 399750001 & end <= 399800000 & profile == "Balancing"),
  chr = "chr1",
  s = 399750001 - 5e6,
  e = 399800000 + 5e6,
  by = 2e6,
  pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
  gr_features = list.features,
  span = 0.5
)
invisible(dev.off())

# Chromosome 2: Found in ZFst
ragg::agg_png(
  filename = here("results/figures/low-div/figure-x-div_profile-chr2-225mb.png"),
  width = 3000,
  height = 3000,
  units = "px",
  res = 300
)
plotRegions(
  df_between = df.div.stats,
  df_within = df.pi_tajD,
  df_outliers_between = tmp |> filter(chromosome == "chr2" & start >= 225100001 & end <= 225300000 & profile == "Balancing"),
  chr = "chr2",
  s = 225100001 - 5e6,
  e = 225300000 + 5e6,
  by = 2e6,
  pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
  gr_features = list.features,
  span = 0.5
)
invisible(dev.off())

# Chromosome 3: Identified in "Allopatric" peak above
ragg::agg_png(
  filename = here("results/figures/low-div/figure-x-div_profile-chr3-145mb.png"),
  width = 3000,
  height = 3000,
  units = "px",
  res = 300
)
plotRegions(
  df_between = df.div.stats,
  df_within = df.pi_tajD,
  df_outliers_between = tmp |> filter(chromosome == "chr3" & start >= 145300001 & end <= 145650000),
  chr = "chr3",
  s = 145300001 - 5e6,
  e = 145650000 + 5e6,
  by = 2e6,
  pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
  gr_features = list.features,
  span = 0.5
)
invisible(dev.off())
