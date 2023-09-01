# ------------------------------------------------------------------------------------------------ #
# ZFst outlier regions: Zoom in and look
#
# This script uses a custom function to explore outlier loci in more detail. I plot the genomic
# region of interest (+/- a few bases), and overlay the divergence statistics and introgression
# values. I also show the genes + repeat sequences over the region.

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

# ------------------------------------------------------------------------------------------------ #
# Genomic outliers
zfst.outliers <- read_csv(
  file = here("results", "between-species", "zfst-outliers.csv"),
  col_names = TRUE,
  col_types = cols()
) |>
  select(comparison, chromosome, start, end, classification) |>
  inner_join(df.div.stats, multiple = "all") |>
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
# Genomic islands of divergence - ZFst > 3

zfst.outliers |>
  filter(classification == 3, statistics == "zfst") |>
  arrange(chromosome, start)

# Chromosome 1: width = 2,349,999 bp, large region with elevated Fst across all three species
# Allopatric: High Fst, average Dxy, dip in pi
ragg::agg_png(
  filename = here("results", "figures", "high-div", "figure-x-zfst-outlier-chr1-212_214mb.png"),
  width = 3000,
  height = 3000,
  units = "px",
  res = 300
)
# cairo_pdf(
#   filename = here("results", "figures", "high-div", "figure-x-zfst-outlier-chr1-212_214mb.pdf"),
#   width = 14,
#   height = 14
# )
plotRegions(
  df_between = df.div.stats,
  df_within = df.pi_tajD,
  df_outliers_between = zfst.outliers |> filter(chromosome == "chr1", start %in% c(212300001, 214550001, 214600001)),
  chr = "chr1",
  s = 212300001 - 10e6,
  e = 214650000 + 10e6,
  by = 2e6,
  pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
  gr_features = list.features,
  span = 0.4
)
invisible(dev.off())

# Chromsome 2: width = 50kbp with potential introgression HEL-HMA
ragg::agg_png(
  filename = here("results/figures/figure-x-zfst-outlier-chr2-240mb.png"),
  width = 3000,
  height = 3000,
  units = "px",
  res = 300
)
# cairo_pdf(
#   filename = here("results", "figures", "high-div", "figure-x-zfst-outlier-chr2-240mb.pdf"),
#   width = 14,
#   height = 10,
# )
plotRegions(
  df_between = df.div.stats,
  df_within = df.pi_tajD,
  df_outliers_between = zfst.outliers |> filter(chromosome == "chr2", start == 240250001),
  chr = "chr2",
  s = 240250001 - 6e6,
  e = 240300000 + 6e6,
  by = 2e6,
  pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
  gr_features = list.features,
  span = 0.4
)
invisible(dev.off())

# Chromosome 6: width = 1949999 - potential introgression event + flanks a repeat spike
ragg::agg_png(
  filename = here("results/figures/figure-x-zfst-outlier-chr6-92mb.png"),
  width = 3000,
  height = 3000,
  units = "px",
  res = 300
)
# cairo_pdf(
#   filename = here("results", "figures", "high-div", "figure-x-zfst-outlier-chr6-92mb-introgressed_recurrent.pdf"),
#   width = 14,
#   height = 12
# )
plotRegions(
  df_between = df.div.stats,
  df_within = df.pi_tajD,
  df_outliers_between = zfst.outliers |> filter(chromosome == "chr6", start %in% c(92450001, 92550001, 93950001, 94350001)),
  chr = "chr6",
  s = 92450001 - 10e6,
  e = 94400000 + 10e6,
  by = 2e6,
  pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
  gr_features = list.features,
  span = 0.4
)
invisible(dev.off())

# Chromosome 11: width = 50kbp (probably just noise)
# plotRegions(
#   df_between = df.div.stats,
#   df_within = df.pi_tajD,
#   df_outliers_between = zfst.outliers |> filter(chromosome == "chr11", start == 27850001),
#   chr = "chr11",
#   s = 27850001 - 1e6,
#   e = 27900000 + 5e5,
#   by = 2e6,
#   pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
#   gr_features = list.features
# )

# ------------------------------------------------------------------------------------------------ #
# Genomic islands of conservation (low ZFst)
zfst.outliers |>
  filter(classification == -3, statistics == "fd") |>
  mutate(chromosome = factor(chromosome, paste0(rep("chr", 15), 1:15))) |>
  arrange(chromosome, start) |>
  View()

# Chromosome 1: width = 950kbp - Looks genuine but is at the start of the chromosome...
ragg::agg_png(
  filename = here("results/figures/low-div/figure-x-zfst-lowlier-chr1-1_2mb.png"),
  width = 3000,
  height = 3000,
  units = "px",
  res = 300
)
# cairo_pdf(
#   filename = here("results/figures/low-div/figure-x-zfst-lowlier-chr1-1_2mb-balancing.pdf"),
#   width = 14,
#   height = 12
# )
plotRegions(
  df_between = df.div.stats,
  df_within = df.pi_tajD,
  df_outliers_between = zfst.outliers |> filter(chromosome == "chr1", start %in% c(1000001, 1200001, 1900001)),
  chr = "chr1",
  s = 1,
  e = 1950000 + 1e6,
  by = 2e5,
  pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
  gr_features = list.features,
  span = 0.3
)
invisible(dev.off())

# Chromosome 1: width = 350kb - many sites all within a peak - inrogression from elegans into the others
ragg::agg_png(
  filename = here("results/figures/low-div/figure-x-zfst-lowlier-chr1-134_135mb.png"),
  width = 3000,
  height = 3000,
  units = "px",
  res = 300
)
# cairo_pdf(
#   filename = here("results/figures/low-div/figure-x-zfst-lowlier-chr1-134_135mb-balancing.pdf"),
#   width = 14,
#   height = 12
# )
plotRegions(
  df_between = df.div.stats,
  df_within = df.pi_tajD,
  df_outliers_between = zfst.outliers |> filter(chromosome == "chr1", start %in% c(134700001, 134900001, 134950001, 135000001)),
  chr = "chr1",
  s = 134700001 - 2e6,
  e = 135050000 + 2e6,
  by = 5e5,
  pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
  gr_features = list.features,
  span = 0.3
)
invisible(dev.off())

# Chromosome 2: width = 50kb
ragg::agg_png(
  filename = here("results/figures/low-div/figure-x-zfst-lowlier-chr2-225mb.png"),
  width = 3000,
  height = 3000,
  units = "px",
  res = 300
)
# cairo_pdf(
#   filename = here("results/figures/low-div/figure-x-zfst-lowlier-chr2-225mb-balancing.pdf"),
#   width = 14,
#   height = 12
# )
plotRegions(
  df_between = df.div.stats,
  df_within = df.pi_tajD,
  df_outliers_between = zfst.outliers |> filter(chromosome == "chr2", start %in% c(225100001)),
  chr = "chr2",
  s = 225100001 - 2e6,
  e = 225150000 + 2e6,
  by = 2e6,
  pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
  gr_features = list.features,
  span = 0.2
)
invisible(dev.off())

# Chromosome 2: width = ~66Mb (whole region is a mess)
ragg::agg_png(
  filename = here("results/figures/low-div/figure-x-zfst-lowlier-chr2-1_70mb.png"),
  width = 3000,
  height = 3000,
  units = "px",
  res = 300
)
# cairo_pdf(
#   filename = here("results/figures/low-div/figure-x-zfst-lowlier-chr2-1_70mb.pdf"),
#   width = 14,
#   height = 12
# )
plotRegions(
  df_between = df.div.stats,
  df_within = df.pi_tajD,
  df_outliers_between = zfst.outliers |> filter(chromosome == "chr2", start <= 65850001),
  chr = "chr2",
  s = 1,
  e = 65850001 + 30e6,
  by = 10e6,
  pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
  gr_features = list.features,
  span = 0.7
)
invisible(dev.off())

# Chromosome 3: width = 50kb - Lowpoint is right before a high Fst (low pi) region
#   - Might be an interesting one to explore
ragg::agg_png(
  filename = here("results/figures/low-div/figure-x-zfst-lowlier-chr3-118mb.png"),
  width = 3000,
  height = 3000,
  units = "px",
  res = 300
)
# cairo_pdf(
#   filename = here("results/figures/low-div/figure-x-zfst-lowlier-chr3-118mb-recurrent_repeat_spike.pdf"),
#   width = 14,
#   height = 12
# )
plotRegions(
  df_between = df.div.stats,
  df_within = df.pi_tajD,
  df_outliers_between = zfst.outliers |> filter(chromosome == "chr3", start == 118800001),
  chr = "chr3",
  s = 118800001 - 5e6,
  e = 118850000 + 20e6,
  by = 2e6,
  pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
  gr_features = list.features,
  span = 0.3
)
invisible(dev.off())

# ------------------------------------------------------------------------------------------------ #
# Unsure outliers

# Chromosome 1: width = Has a big gap 3' so difficult to say
# plotRegions(
#   df_between = df.div.stats,
#   df_within = df.pi_tajD,
#   df_outliers_between = zfst.outliers |> filter(chromosome == "chr1", start == 3400001),
#   chr = "chr1",
#   s = 3400001 - 2e6,
#   e = 3450000 + 2e6,
#   by = 5e5,
#   pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
#   gr_features = list.features
# )

# Chromosome 1: width = 50kbp single loci that is in a pretty unremarkable area
# plotRegions(
#   df_between = df.div.stats,
#   df_within = df.pi_tajD,
#   df_outliers_between = zfst.outliers |> filter(chromosome == "chr1", start == 9450001),
#   chr = "chr1",
#   s = 9450001 - 2e6,
#   e = 9500000 + 2e6,
#   by = 2e5,
#   pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
#   gr_features = list.features
# )

# Chromosome 2: 50Kb window that appears to have a bit of a dip in H. elegans/major comparison
# plotRegions(
#   df_between = df.div.stats,
#   df_within = df.pi_tajD,
#   df_outliers_between = zfst.outliers |> filter(chromosome == "chr2", start %in% c(227450001)),
#   chr = "chr2",
#   s = 227450001 - 1e6,
#   e = 227500000 + 1e6,
#   by = 1e6,
#   pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
#   gr_features = list.features,
#   span = 1
# )

# Chromosome 2: 50kb window that has a bit of a dip but seems to be driven by a single window
# plotRegions(
#   df_between = df.div.stats,
#   df_within = df.pi_tajD,
#   df_outliers_between = zfst.outliers |> filter(chromosome == "chr2", start == 349800001),
#   chr = "chr2",
#   s = 349800001 - 1e6,
#   e = 349850000 + 1e6,
#   by = 2e6,
#   pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
#   gr_features = list.features
# )

# Chromosome 3:50kb window. Looks like a drop but hard to say if genuine
# plotRegions(
#   df_between = df.div.stats,
#   df_within = df.pi_tajD,
#   df_outliers_between = zfst.outliers |> filter(chromosome == "chr3", start == 5800001),
#   chr = "chr3",
#   s = 5800001 - 2e6,
#   e = 5850000 + 2e6,
#   by = 2e5,
#   pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
#   gr_features = list.features
# )

# Chromosome 3 - width = 50kbp Hard to know, as not far from chromosome end/only one window
# plotRegions(
#   df_between = df.div.stats,
#   df_within = df.pi_tajD,
#   df_outliers_between = zfst.outliers |> filter(chromosome == "chr3", start == 234750001),
#   chr = "chr3",
#   s = 234750001 - 2e6,
#   e = 234800000 + 2e6,
#   by = 5e5,
#   pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
#   gr_features = list.features
# )

# Chromosome 4 - dips but is shallow
# plotRegions(
#   df_between = df.div.stats,
#   df_within = df.pi_tajD,
#   df_outliers_between = zfst.outliers |> filter(chromosome == "chr4", start == 87350001),
#   chr = "chr4",
#   s = 87350001 - 1e6,
#   e = 87400000 + 1e6,
#   by = 2e6,
#   pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
#   gr_features = list.features
# )

# Chromosome 5 - 50kbp There's a dip, but not really convincing
# plotRegions(
#   df_between = df.div.stats,
#   df_within = df.pi_tajD,
#   df_outliers_between = zfst.outliers |> filter(chromosome == "chr5", start == 42650001),
#   chr = "chr5",
#   s = 42650001 - 1e6,
#   e = 42700000 + 1e6,
#   by = 5e5,
#   pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
#   gr_features = list.features
# )

# Chromosome 5 - width 50kbp Kind of a dip, but only a single window
# plotRegions(
#   df_between = df.div.stats,
#   df_within = df.pi_tajD,
#   df_outliers_between = zfst.outliers |> filter(chromosome == "chr5", start %in% c(132000001)),
#   chr = "chr5",
#   s = 132000001 - 1e6,
#   e = 132050000 + 1e6,
#   by = 5e5,
#   pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
#   gr_features = list.features
# )

# Chromosome 5 - width 50kb little dip but also right near the end of chromosome
# plotRegions(
#   df_between = df.div.stats,
#   df_within = df.pi_tajD,
#   df_outliers_between = zfst.outliers |> filter(chromosome == "chr5", start %in% c(139950001)),
#   chr = "chr5",
#   s = 139950001 - 1e6,
#   e = 140000000 + 1e6,
#   by = 5e5,
#   pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
#   gr_features = list.features
# )

# Chromosome 6 - width 50kb low point, but not super dramatic and also right near start of chromosome
# plotRegions(
#   df_between = df.div.stats,
#   df_within = df.pi_tajD,
#   df_outliers_between = zfst.outliers |> filter(chromosome == "chr6", start == 6100001),
#   chr = "chr6",
#   s = 1,
#   e = 6100001 + 10e6,
#   by = 5e5,
#   pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
#   gr_features = list.features
# )

# Chromosome 9 - single point that is very small (almost -4 ZFst) but seems like a single outlier
# plotRegions(
#   df_between = df.div.stats,
#   df_within = df.pi_tajD,
#   df_outliers_between = zfst.outliers |> filter(chromosome == "chr9", start == 8500001),
#   chr = "chr9",
#   s = 8500001 - 2e6,
#   e = 8550000 + 2e6,
#   by = 5e5,
#   pal = c("#66c2a5", "#fc8d62", "#8da0cb", "#e41a1c", "#377eb8", "#4daf4a"),
#   gr_features = list.features
# )

# 15 unknown
