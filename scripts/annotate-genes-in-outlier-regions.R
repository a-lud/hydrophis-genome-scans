# ------------------------------------------------------------------------------------------------ #
# Annotate genes within ZFst/Profiled regions
#
# We manually curated ZFst high/low divergence windows, along with the Allopatric/recurrent windows
# identified via the divergence profiling appraoch. These windows represent some of the strongest
# signals of divergence between these three species, so we decided they were the best targets to
# investigate.
#
# Genes within candidate regions were extracted using the script 'extract-genes-in-outlier.sh'
#
# This script uses annotation information from previous work + new annotations generated via emapper
# + EggNOG to assign gene symbols and descriptions to each gene within a region.

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
})

# ------------------------------------------------------------------------------------------------ #
# Read in recurrent table
zfst.genes <- read_tsv(
  file = here("results", "outlier-windows", "zfst-genes-in-high_low-divergence.tsv"),
  col_names = c("classification", "chromosome", "region_start", "region_end", "transcript"),
  col_types = cols()
)

prof.genes <- read_tsv(
  file = here("results", "outlier-windows", "profile-genes-in-allopatric_recurrent.tsv"),
  col_names = c("classification", "chromosome", "region_start", "region_end", "transcript"),
  col_types = cols()
)

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
  rename(saccver = From, description_uniprot = `Protein names`, gene_uniprot = `Gene Names`) |>
  mutate(gene_uniprot = toupper(str_split_i(gene_uniprot, " ", 1)))

# ------------------------------------------------------------------------------------------------ #
# All annotations together
df.homology.annotations <- full_join(emapper, blastp) |>
  full_join(df.unimap, multiple = "all") |>
  relocate(transcript, gene_eggnog, gene_uniprot, description_uniprot, description_eggnog, saccver)

write_csv(x = df.homology.annotations, file = here("data", "annotation-table.csv"), col_names = TRUE)

# ------------------------------------------------------------------------------------------------ #
# Annotate genes (some exist already - obviously, but assign gene names to those that are not)
zfst.genes <- zfst.genes |>
  left_join(df.homology.annotations, multiple = "all") |>
  distinct(transcript, .keep_all = TRUE) # eggnog description for one transcript causes strife

prof.genes <- prof.genes |>
  left_join(df.homology.annotations, multiple = "all") |>
  distinct(transcript, .keep_all = TRUE)

# ------------------------------------------------------------------------------------------------ #
# Summary statistics

# Number of genes that could be assigned a gene symbol - N=270 (total=426) = 63.4%
zfst.genes |>
  # One of the eggnogg descriptions has two different texts - choose first
  distinct(transcript, .keep_all = TRUE) |>
  filter(!if_all(c(gene_eggnog, gene_uniprot), is.na)) |>
  nrow()
length(unique(zfst.genes$transcript))

# N=164 (total=316) = 52.5%
prof.genes |>
  # One of the eggnogg descriptions has two different texts - choose first
  distinct(transcript, .keep_all = TRUE) |>
  filter(!if_all(c(gene_eggnog, gene_uniprot), is.na)) |>
  nrow()
length(unique(prof.genes$transcript))

# Number of genes within each region
zfst.n_genes <- zfst.genes |>
  mutate(width = region_end - region_start) |>
  group_by(classification, chromosome, region_start, region_end, width) |>
  summarise(n_genes = n())

prof.n_genes <- prof.genes |>
  mutate(width = region_end - region_start) |>
  group_by(classification, chromosome, region_start, region_end, width) |>
  summarise(n_genes = n())

bind_rows(zfst.n_genes, prof.n_genes) |>
  arrange(classification, width) |>
  write_csv(file = "manuscript/tables/number-genes-in-curated-regions.csv", col_names = TRUE)

# Number regions by chromosome
bind_rows(zfst.n_genes, prof.n_genes) |>
  group_by(classification, chromosome) |>
  summarise(n_chromosome = n())

# Correlation between region size and gene count
ragg::agg_png(
  filename = here("results/figures/figure-x-region_width-vs-gene_count.png"),
  width = 2000,
  height = 2000,
  units = "px",
  res = 300
)
bind_rows(zfst.n_genes, prof.n_genes) |>
  arrange(classification, width) |>
  mutate(classification = ifelse(classification == "high_div", "High divergence", "Low divergence")) |>
  ggplot(
    aes(
      x = width,
      y = n_genes,
      colour = classification
    )
  ) +
  geom_point(size = 3) +
  geom_smooth(se = FALSE, method = "lm") +
  scale_colour_manual(values = c("black", "red")) +
  scale_x_continuous(
    # breaks = seq(0, 12.5e6, 1e6),
    labels = scales::label_number(scale = 1e-6)
  ) +
  # scale_y_continuous(breaks = seq(0,140, 20)) +
  labs(
    colour = "Region type",
    x = "Width (Mb)",
    y = "Number of Genes"
  ) +
  ggpubr::stat_cor(
    digits = 3,
    label.x.npc = "center",
    label.y.npc = "bottom",
    p.accuracy = 0.001,
    show.legend = FALSE,
    cor.coef.name = "r",
    size = 3
  ) +
  facet_wrap(facets = vars(classification), scales = "free") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = c(0.7, 0.85),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
invisible(dev.off())

# ------------------------------------------------------------------------------------------------ #
# Wilcoxon test
tmp <- bind_rows(zfst.n_genes, prof.n_genes) |>
  arrange(classification, width)

# Are sizes of low-div regions different to high div
wilcox.test(width ~ classification, data = tmp, exact = FALSE)

# ------------------------------------------------------------------------------------------------ #
# Write to file
zfst.genes |>
  select(classification, chromosome, region_start, region_end, transcript, gene_uniprot, gene_eggnog, description_uniprot, saccver, `Tissue specificity`) |>
  mutate(chromosome = factor(chromosome, paste0(rep("chr", 15), 1:15))) |>
  arrange(chromosome, region_start) |>
  write_csv(
    file = here("results", "outlier-windows", "zfst-genes-in-high_low-divergence-annotation.csv"),
    col_names = TRUE,
    na = ""
  )

prof.genes |>
  select(classification, chromosome, region_start, region_end, transcript, gene_uniprot, gene_eggnog, description_uniprot, saccver, `Tissue specificity`) |>
  mutate(chromosome = factor(chromosome, paste0(rep("chr", 15), 1:15))) |>
  arrange(chromosome, region_start) |>
  write_csv(
    file = here("results", "outlier-windows", "profile-genes-in-allopatric_recurrent-annotation.csv"),
    col_names = TRUE,
    na = ""
  )
