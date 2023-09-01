plotRegions <- function(
    df_between, df_within, df_outliers_between, chr, s, e, by,
    pal, scale_factor = 1e-6,
    gr_features = NULL,
    span = 0.2
) {
  # Filter input dataframe by upper/lower bounds
  df.subset_between <- df_between |>
    filter(chromosome == chr, midpoint >= s & midpoint <= e)

  df.subset_within <- df_within |>
    filter(chromosome == chr, midpoint >= s & midpoint <= e)

  df.test <- list(df.subset_between, df.subset_within) |> list_rbind()

  # Extract outlier windows for individual samples (for within species data)
  #   - i.e. only want to highlight points for samples that actually report an outlier
  df_outliers_within <- df_outliers_between |>
    separate(col = comparison, into = c("A", "B"), sep = "-") |>
    select(A, B, chromosome, start, end) |>
    pivot_longer(names_to = "temp", values_to = "comparison", cols = c(A, B)) |>
    select(-temp) |>
    distinct() |>
    inner_join(df.subset_within, multiple = "all")

  df.test.outlier <- list(df_outliers_between, df_outliers_within) |>
    list_rbind()

  # Plot all regions
  plot.region <- df.test |>
    # Between species (ZFst/Dxy/fd)
    ggplot(
      aes(
        x = midpoint,
        y = values,
        colour = comparison
      )
    ) +
    geom_point(alpha = 0.2) +
    # geom_smooth(linewidth = 2, alpha = 0.1, span = 0.2) +
    geom_line(stat = "smooth", method = "loess", span = span, alpha = 0.5, linewidth = 1) +

    # Highligh "outlier" windows in red
    geom_point(
      data = df.test.outlier,
      mapping = aes(x = midpoint,y = values),
      size = 2
    ) +
    scale_x_continuous(
      expand = c(0,0),
      breaks = seq(s, e, by),
      limits = c(s, e),
      labels = scales::label_number(scale = 1e-6)
    ) +
    scale_y_continuous(position = "right") +
    scale_colour_manual(values = pal) +
    labs(
      x = paste0("\nPosition ", chr, " (Mb)"),
      colour = "Comparison/sample"
    ) +
    facet_grid(
      rows = vars(statistics),
      labeller = labeller(
        statistics = as_labeller(c(zfst = "ZF[st]", dxy = "D[xy]", fd = "\U0192[d]", pi = "\U03C0", TajimaD = "Tajimas~D"), default = label_parsed),
        value = label_value
      ),
      scales = "free_y",
      switch = "y"
    ) +
    guides(colour = guide_legend(title.position = "bottom", title.hjust = 0.5, nrow = 1)) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      axis.title.y = element_blank(),

      # Facet strip
      strip.background = element_blank(),
      strip.text.y = element_text(size = 16),
      strip.placement = "outside",

      #
      legend.position = "bottom"
    )

  # Plot feature tracks if provided
  plot.features <- NULL
  if (!is.null(gr_features)) {
    if("Genes" %in% names(gr_features)) {
      df.gene.features <- getOverlappingFeatures(
        df_windows = df_outliers_between,
        gr_gff = gr_features[["Genes"]],
        chr = chr,
        s = s,
        e = e
      ) |>
        mutate(Feature = "Genes") |>
        mutate(Overlap = factor(Overlap, levels = c("True", "False")))

      plot.gene_features <- df.gene.features |>
        ggplot() +
        geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = Overlap), alpha = 0.5, show.legend = FALSE) +
        scale_fill_manual(values = if(length(unique(df.gene.features$Overlap)) == 2) c("red","grey80") else "grey80") +
        scale_x_continuous(expand = c(0,0),limits = c(s, e)) +
        scale_y_continuous(expand = c(0,0)) +
        facet_grid(rows = vars(Feature), switch = "y") +
        guides(fill = "none") +
        theme_bw() +
        theme(
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          strip.background = element_blank(),
          strip.text.y = element_text(size = 14),
          legend.position = "none",
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )
      if(is.null(plot.features)) {
        plot.features <- plot.gene_features
      } else {
        plot.features <- plot.gene_features / plot.features
      }
    }

    if("Repeats" %in% names(gr_features)) {
      df.repeat.feature.counts <- gr_features[["Repeats"]] |>
        as_tibble() |>
        rename(chromosome = seqnames) |>
        mutate(
          midpoint = (start + end)/2,
          Feature = "Repeats"
        ) |>
        filter(chromosome == chr, start >= s & end <= e)

      plot.repeat_density <- df.repeat.feature.counts |>
        ggplot(
          aes(
            x = (start + end)/2,
            y = count
          )
        ) +
        geom_line() +
        scale_x_continuous(expand = c(0,0),limits = c(s, e)) +
        scale_y_continuous(expand = c(0,0), position = "right") +
        facet_grid(rows = vars(Feature), switch = "y") +
        theme_bw() +
        theme(
          axis.text = element_blank(),
          axis.text.y = element_text(size = 14),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          strip.background = element_blank(),
          strip.text.y = element_text(size = 14),
          legend.position = "none",
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        )

      if(is.null(plot.features)) {
        plot.features <- plot.repeat_density
      } else {
        plot.features <- plot.features / plot.repeat_density
      }
    }
  } else {
    plot.features <- NULL
  }

  # patchwork the plots together
  if(!is.null(plot.features)) {
    feature.height <- if(length(gr_features) > 1) c(0.1, 0.1) else 0.1
    return(
      plot.features + plot.region |>
        wrap_plots() +
        plot_layout(height = c(feature.height, 1), guides = "collect") &
        theme(
          legend.position = "bottom",
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 14)
        )
    )
  } else {
    return(plot.region)
  }
}

plotRegions10kb <- function(
    df_between, df_within, chr, s, e, by,
    pal, scale_factor = 1e-6,
    gff, os, oe,
    span = 0.2
) {
  # Filter input dataframe by upper/lower bounds
  df.subset_between <- df_between |>
    filter(chromosome == chr, midpoint >= s & midpoint <= e)

  df.subset_within <- df_within |>
    filter(chromosome == chr, midpoint >= s & midpoint <= e)

  df.test <- list(df.subset_between, df.subset_within) |> list_rbind()

  # Plot all regions
  plot.region <- df.test |>
    # Between species (Fst/Dxy/fd)
    ggplot(
      aes(
        x = midpoint,
        y = values,
        colour = comparison
      )
    ) +
    geom_point(alpha = 0.2) +
    geom_line(stat = "smooth", method = "loess", span = span, alpha = 0.7, linewidth = 1) +
    scale_x_continuous(
      expand = c(0,0),
      breaks = seq(s, e, by),
      limits = c(s, e),
      labels = scales::label_number(scale = 1e-6)
    ) +
    scale_y_continuous(position = "right") +
    scale_colour_manual(values = pal) +
    labs(
      x = paste0("\nPosition ", chr, " (Mb)"),
      colour = "Comparison/sample"
    ) +
    facet_grid(
      rows = vars(statistics),
      labeller = labeller(
        statistics = as_labeller(c(fst = "F[st]", dxy = "D[xy]", pi = "\U03C0"), default = label_parsed),
        value = label_value
      ),
      scales = "free_y",
      switch = "y"
    ) +
    guides(colour = guide_legend(title.position = "bottom", title.hjust = 0.5, nrow = 1)) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      axis.title.y = element_blank(),

      # Facet strip
      strip.background = element_blank(),
      strip.text.y = element_text(size = 16),
      strip.placement = "outside",

      #
      legend.position = "bottom"
    )

  # Get overlapping gene features
  df.genes <- gff |>
    filter(
      chromosome == chr, start >= s  & end <= e
    ) |>
    mutate(
      overlap = case_when(
        start >= os & end <= oe ~ "TRUE",
        .default = "FALSE"
      )
    )

  plot.gene_features <- df.genes |>
    ggplot() +
    geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = overlap), alpha = 0.5, show.legend = FALSE) +
    scale_fill_manual(values = if(length(unique(df.genes$overlap)) == 2) c("grey", "red") else "grey80") +
    scale_x_continuous(expand = c(0,0),limits = c(s, e)) +
    scale_y_continuous(expand = c(0,0)) +
    # facet_grid(rows = vars(Feature), switch = "y") +
    guides(fill = "none") +
    theme_bw() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      # strip.background = element_blank(),
      # strip.text.y = element_text(size = 14),
      legend.position = "none",
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  # patchwork the plots together
  return(
    plot.gene_features / plot.region +
      plot_layout(height = c(0.1, 1), guides = "collect") &
      theme(
        legend.position = "bottom",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)
      )
  )
}

# Extract genes that overlap with windows
getOverlappingFeatures <- function(df_windows, gr_gff, chr, s, e) {
  # gr.windows <- GenomicRanges::makeGRangesFromDataFrame(df_windows, keep.extra.columns = TRUE) |>
  #   plyranges::reduce_ranges()

  # Make a GRange that covers the entire window of interest - from left most outlier to rightmost
  # Opting for this as, even though all windows might not intersect with a specific feature, the
  # features within the bounds of the outlier windows will be interesting
  gr.window <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges = IRanges::IRanges(
      start = min(df_windows$start),
      end = max(df_windows$end))
  )

  # Subset for genes if feature file is a gene GFF3
  if("gene" %in% gr_gff$feature) {
    gr_gff <- gr_gff[gr_gff$feature == "gene"]
  }

  # Subset feature file for elements within region of interest
  df.features.region <- gr_gff[
    GenomicRanges::start(gr_gff) >= s
    & GenomicRanges::end(gr_gff) <= e
    & GenomicRanges::seqnames(gr_gff) == chr
  ] |>
    as_tibble() |>
    select(chromosome = seqnames, start, end)

  # Return genes that overlap with windows
  df.tmp <- plyranges::find_overlaps(
    gr_gff,
    gr.window
  ) |>
    as_tibble() |>
    select(chromosome = seqnames, start, end) |>
    mutate(Overlap = "True")

  # Join with full dataset
  df.features.region |>
    select(chromosome, start, end) |>
    left_join(df.tmp) |>
    mutate(Overlap = ifelse(is.na(Overlap), "False", Overlap))
}
