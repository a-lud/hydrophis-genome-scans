# ------------------------------------------------------------------------------------------------ #
# ZFst outlier's overlapping genes
#
# Using BEDtools intersect, we identified which genes Zfst outlier's intersected (by at least 5kb
# of window width). The genes are summarised below.

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

# ------------------------------------------------------------------------------------------------ #
# Read in intersected data
df.intersections <- read_tsv(
  file = here("results", "outlier-windows", "zfst-outlier-intersect-gff.tsv"),
  col_names = c("chromosome", "start", "end", "comparison", "fst", "zfst", "dxy", "window_type", "chr", "method", "type", "start_gff", "end_gff", "score", "strand", "phase", "attributes"),
  col_types = cols(),
  col_select = c("comparison", "chromosome", "start", "end", "fst", "zfst", "dxy", "window_type", "chr", "type", "start_gff", "end_gff", "attributes")
) |>
  filter(type == "gene") |>
  mutate(
    gene = ifelse(str_detect(attributes, "Name="), sub(".*Name=(.*);", "\\1", attributes), NA_character_),
    gene = str_remove(gene, "_.*"),
    transcript = sub("ID=(.*);.*", "\\1", attributes),
    transcript = str_remove(transcript, ";Name=.*")
  ) |>
  select(-attributes)

# Get unique transcript IDs
interserctions.transcripts <- df.intersections |>
  pull(transcript) |>
  unique()

# ------------------------------------------------------------------------------------------------ #
# Export transcript IDs for annotation with emapper + eggNOG DB - only needs to be run once
interserctions.transcripts |>
  paste0("-T1") |>
  sort() |>
  write_lines(
    file = here("results/outlier-windows/zfst-transcript-ids-to-annotate.txt")
  )

# ------------------------------------------------------------------------------------------------ #
# UniProt (SwissProt) BLASTP - UniProt accessions for H. major transcripts
# Export uniprot accessions for transcripts and annotate via UniProts IDmapping tool
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
  ungroup() |>
  filter(transcript %in% interserctions.transcripts)

# write UniProt identifiers to file
blastp |>
  pull(saccver) |>
  unique() |>
  write_lines(file = here("results", "outlier-windows", "zfst-outlier-uniprot-ids.txt"))

# ------------------------------------------------------------------------------------------------ #
# Import annotations
#   - table from emapper + eggNOG
#   - table from IDmapping (uniprot)
eggNOG.annotations <- read_tsv(
  file = here("results","outlier-windows","zfst-outlier-genes.emapper.annotations"),
  col_names = TRUE,
  col_types = cols(),
  comment = "##"
) |>
  select(query = `#query`, symbol = Preferred_name, Description) |>
  mutate(query = str_remove(query, "-T.*")) |>
  filter(symbol != "-")

uniprot.annotations <- read_tsv(
  file = here("results", "outlier-windows", "zfst-outlier-unirprot-mapped-ids.gz"),
  col_names = TRUE,
  col_types = cols()
) |>
  mutate(
    `Gene Names` = map(.x = str_split(`Gene Names`, " "), .f = first) |> unlist(),
    `Gene Names` = toupper(`Gene Names`)
  ) |>
  select(saccver = From, `Protein names`, `Gene Names`)

# ------------------------------------------------------------------------------------------------ #
# Gene annotation table after merging all annotation sources
df.annotated.intersections <- df.intersections |>
  full_join(blastp) |>
  full_join(uniprot.annotations) |>
  select(
    comparison, chromosome, start, end, window_type,`Gene Names`,
    `Protein names`, transcript, pident, qcovs
  )

# Write to file
df.annotated.intersections |>
  write_csv(
    file = here("results/tables/zfst-outlier-genes.csv"),
    col_names = TRUE
  )

# ------------------------------------------------------------------------------------------------ #
# Summary tables

# Number of genes intersecting each class of outlier (ZFst >= 3 or ZFst <= -3)
# window_type `HCU-HEL` `HCU-HMA` `HEL-HMA`
# -3           10        28        61
#  3           16         5        NA
df.annotated.intersections |>
  group_by(comparison, window_type) |>
  summarise(count = n())|>
  pivot_wider(names_from = comparison, values_from = count) |>
  write_csv(
    file = here("results/tables/zfst-count-genes-intersecting-outlier.csv"),
    col_names = TRUE
  )

# Formatted by shared genes
df.annotated.intersections |>
  filter(!is.na(`Gene Names`)) |>
  group_by(`Gene Names`) |>
  mutate(
    comparison = paste0(sort(unique(comparison)), collapse = "|"),
    comparison = factor(comparison, levels = c(
      "HCU-HEL", "HCU-HMA", "HEL-HMA",
      "HCU-HEL|HCU-HMA", "HCU-HMA|HEL-HMA",
      "HCU-HEL|HCU-HMA|HEL-HMA"
    ))
  ) |>
  select(comparison, `Gene Names`, transcript, window_type, `Protein names`) |>
  arrange(comparison) |>
  distinct(`Gene Names`, .keep_all = TRUE) |>
  write_csv(
    file = here("results/tables/zfst-outlier-genes-grouped.csv"),
    col_names = TRUE
  )
