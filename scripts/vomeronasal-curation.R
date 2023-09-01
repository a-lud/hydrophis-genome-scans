# ------------------------------------------------------------------------------------------------ #
# Annotate V2R genes
#

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(gggenes)
})

# ------------------------------------------------------------------------------------------------ #
# Read in BLAST table
df.blast <- read_tsv(
  file = here('results/annotate-vomeronasal/v2r-blast.tsv'),
  col_names = c(
    'qseqid','qlen','sseqid','sacc','sgi',
    'slen','qstart','qend','sstart','send',
    'evalue','bitscore','length','pident',
    'mismatch','gaps','qcovs'
  ),
  col_types = cols()
) |>
  mutate(qseqid = str_remove(qseqid, "-T1.*"))

df.annotation <- read_tsv(
  "data/uniprotkb_Vomeronasal_2023_08_30.tsv.gz",
  col_names = TRUE,
  col_types = cols()
)

# ------------------------------------------------------------------------------------------------ #
# get best hits
# FUN_000096 = nicotinamide riboside kinase 2
# FUN_000070 = hypothetical
# FUN_000072 = No hit
# FUN_000083 = No Hit

df.annotated <- df.blast |>
  group_by(qseqid) |>
  slice(1) |>
  left_join(df.annotation, join_by(sacc == Entry)) |>
  select(transcript = qseqid, sacc, pident, qcovs, `Gene Names`, `Protein names`) |>
  filter(! transcript %in% c("FUN_000096", "FUN_000070", "FUN_000072", "FUN_000083"))

# ------------------------------------------------------------------------------------------------ #
# Gene coordinates
gff <- read_tsv(
  file = here("data", "hydrophis_major.gff3"),
  col_names = c("chromosome", "method", "type", "start", "end", "score", "strand", "phase", "attribute"),
  col_types = cols(),
  comment = "#"
) |>
  filter(type == "gene") |>
  separate(col = attribute, into = c("transcript", "symbol"), sep = ";") |>
  mutate(
    symbol = str_remove(symbol, "Name="),
    transcript = str_remove(transcript, "ID=")
  )

# Annotated genes
ann <- read_csv(
  here("results/outlier-windows/zfst-genes-in-high_low-divergence-annotation.csv"),
  col_names = TRUE,
  col_types = cols()
) |>
  filter(chromosome == "chr1" &  region_start >= 8e5 & region_end <= 2e6) |>
  inner_join(df.annotated)

df.test <- gff |>
  inner_join(ann) |>
  select(start, end, chromosome, strand, gene_uniprot) |>
  mutate(
    orientation = ifelse(strand == "+", 0, 1),
    strand = ifelse(strand == "+", "forward", "reverse")
  )

cairo_pdf(
  filename = here("manuscript", "figures", "vmn-gene-casette.pdf"),
  width = 10,
  height = 4
)
ggplot(df.test, aes(xmin = start, xmax = end, y = chromosome, forward = orientation)) +
  geom_gene_arrow() +
  theme_genes() +
  facet_wrap(facets = vars(chromosome), scales = "free", nrow = 1) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 14),
    legend.position = "bottom"
  )
invisible(dev.off())
