# ------------------------------------------------------------------------------------------------ #
# Overlap betwen ZFst and "Barriers to gene flow"

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

# ------------------------------------------------------------------------------------------------ #
# Load data
zfst.outliers <- read_csv(
  file = here("results/between-species/zfst-outliers.csv"),
  col_names = TRUE,
  col_types = cols()
) |>
  select(comparison, chromosome, start, end, classification)

# Outliers and lowliers
zfst.high <- zfst.outliers |> filter(classification == 3)
zfst.low <- zfst.outliers |> filter(classification == -3)

df.profile <- read_tsv(
  file = here("results", "tables", "classified-windows.tsv"),
  col_names = TRUE,
  col_types = cols()
) |>
  select(comparison, chromosome, start, end, profile)

df.profile |>
  inner_join(zfst.high)

df.profile |>
  inner_join(zfst.low)
