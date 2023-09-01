# ------------------------------------------------------------------------------------------------ #
# ZFst outliers/lowliers: Gene information
#
# This script has the code (with some manual steps) to annotate and get functional information about
# the genes within the ZFst outliers/lowliers. Specifically:
#
# 1. Extract genes within manually curated regions (see zfst-zoom-outlier-regions.R)
# 2. Export transcript IDs (even for genes with annotations)
# 3. Extract protein sequences (outside this script)
# 4. emapper (eggNOG) annotation (outside this script)
# 5. Collate gene symbols + descriptions

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

# ------------------------------------------------------------------------------------------------ #
# Read in outliers + filter for genuine candidates
gff <- read_tsv(
  file = here("data", "hydrophis_major.gff3"),
  col_names = c("chromosome", "method", "type", "start", "end", "score", "strand", "phase", "attribute"),
  col_types = cols(),
  comment = "#"
) |>
  filter(type == "gene")

# Genes in outlier regions
transcripts.outliers <- gff |>
  filter(
    # chromosome == "chr1" & start >= 210e6 & end <= 217e6 |
      chromosome == "chr2" & start >= 238e6 & end <= 242.3e6
      # chromosome == "chr3" & start >= 119.5e6 & end <= 127e6|
      # chromosome == "chr6" & start >= 91.25e6 & end <= 94.75e6
  ) |>
  pull(attribute) |>
  str_remove(";.*") |>
  str_remove("ID=")

# Genes in lowlier regions
transcripts.lowliers <- gff |>
  filter(
    chromosome == "chr1" & start >= 134.7e6 & end <= 135.05e6 |
      chromosome == "chr2" & start >= 225e6 & end <= 225.3e6 |
      chromosome == "chr3" & start >= 119.5e6 & 127e6
  ) |>
  pull(attribute) |>
  str_remove(";.*") |>
  str_remove("ID=")

# ------------------------------------------------------------------------------------------------ #
# Import emapper results (against eggNOG) and BLASTP results
emapper <- read_tsv(
  file = here("results","outlier-windows","all-proteins.emapper.annotations"),
  col_names = TRUE,
  col_types = cols(),
  comment = "##"
) |>
  select(transcript = `#query`, gene_eggnog = Preferred_name, description_eggnog = Description) |>
  mutate(transcript = str_remove(transcript, "-T.*")) |>
  filter(gene_eggnog != "-") |>
  distinct()

# ------------------------------------------------------------------------------------------------ #
# BLASTP results
blastp <- read_tsv(
  file = here("data", "hydrophis_major.outfmt6.gz"),
  col_names = c(
    'qaccver', 'saccver', 'qlen', 'slen', 'length',
    'qcovs', 'pident', 'mismatch', 'gapopen', 'qstart',
    'qend', 'sstart', 'send', 'evalue', 'bitscore'
  ),
  col_types = cols()
) |>
  mutate(qaccver = str_remove(qaccver, "-T.*")) |>
  select(transcript = qaccver, saccver, pident, qcovs) |>
  group_by(transcript) |>
  slice(1) |>
  ungroup()

# Uniprot id mappings (premade)
df.unimap <- read_tsv(
  file = here("data/idmapping-uniprot-all-genes-hmaj.tsv.gz"),
  col_names = TRUE,
  col_types = cols()
) |>
  select(saccver = From, description_uniprot = `Protein names`, gene_uniprot = `Gene Names`, Pathway, `Tissue specificity`) |>
  mutate(gene_uniprot = toupper(str_split_i(gene_uniprot, " ", 1)))

# ------------------------------------------------------------------------------------------------ #
# All annotations together
df.homology.annotations <- full_join(emapper, blastp) |>
  full_join(df.unimap, multiple = "all") |>
  select(transcript, saccver, gene_eggnog, gene_uniprot, description_uniprot, Pathway, `Tissue specificity`)

# ------------------------------------------------------------------------------------------------ #
# Subset for genes overlapping regions of interest
df.homology.annotations[df.homology.annotations$transcript %in% transcripts.outliers,]
df.homology.annotations[df.homology.annotations$transcript %in% transcripts.lowliers,]
