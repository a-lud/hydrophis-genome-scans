# Helper functions for 03-diversity-statistics

plotRegion <- function(
    between, within, centromere = NULL,
    chr, start, end, sstart, send,
    dxy_ymax = 0.01, dxy_y_by = 0.002,
    pi_ymax = 0.01, pi_y_by = 0.002,
    x_by = 2e6,
    tag = NULL,
    palette_between,
    palette_within,
    return_legends = FALSE
) {

  # -------------------------------- #
  # plotting variables

  # custom scales in facet
  scales_facet <- list(
    # ZFst
    scale_y_continuous(
      limits = c(-5, 3),
      breaks = c(-3, 0, 3),
      labels = as.character(c(-3, 0, 3))
    ),
    # Dxy
    scale_y_continuous(
      limits = c(0, dxy_ymax),
      breaks = seq(0, dxy_ymax, dxy_y_by),
      labels = as.character(seq(0, dxy_ymax, dxy_y_by))
    )
  )

  # Custom facet Labels
  between_labeller <- as_labeller(
    c(
      avg_dxy = "d[XY]",
      zfst = "ZF[ST]"
    ),
    default = label_parsed
  )

  # -------------------------------- #
  # Plot: Between div.
  p_between <- between |>
    filter(chromosome == chr, window_pos_1 >= start, window_pos_2 <= end) |>
    ggplot(
      aes(
        x = floor((window_pos_1 + window_pos_2)/2),
        y = values,
        colour = comp_label,
        linetype = comp_label
      )
    ) +
    geom_line(linewidth = 0.8, alpha = 0.3) +
    geom_smooth(span = 0.2, se = FALSE) +
    geom_hline(
      yintercept = c(-3, 3),
      col = "black",
      linewidth = 0.5,
      linetype = "dotted"
    ) +
    annotate(
      "rect",
      xmin = sstart, xmax = send,
      ymin = -Inf, ymax = Inf,
      alpha = .4,
      fill = "grey60"
    ) +
    scale_x_continuous(
      name = glue::glue("\nChromosome {chr} (Mb)"),
      limits = c(start, end),
      labels = scales::label_number(scale = 1e-6, drop0trailing = TRUE),
      breaks = seq(start, end, x_by),
      expand = c(0, 0)
    ) +
    labs(y = NULL) +
    scale_colour_manual(values = palette_between) +
    guides(
      colour = guide_legend(
        ncol = 1,
        override.aes = list(linewidth = 3)
      )
    ) +
    facet_grid(
      vars(statistic),
      scales = "free_y",
      labeller = labeller(
        statistic = between_labeller,
        value = label_value
      )
    ) +
    ggh4x::facetted_pos_scales(y = scales_facet) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(face = "italic"),
      plot.margin=grid::unit(c(0,0,0,0), "mm"),,
      strip.text = element_text(size = 16)
    )

  # -------------------------------- #
  # Plot: Within div
  p_within <- within |>
    filter(chromosome == chr, window_pos_1 >= start, window_pos_2 <= end) |>
    ggplot(
      aes(
        x = floor((window_pos_1 + window_pos_2)/2),
        y = avg_pi,
        colour = species,
        linetype = species
      )
    ) +
    geom_line(linewidth = 0.8, alpha = 0.3) +
    geom_smooth(span = 0.2, se = FALSE) +
    annotate(
      "rect",
      xmin = sstart, xmax = send,
      ymin = -Inf, ymax = Inf,
      alpha = .4,
      fill = "grey60"
    ) +
    scale_x_continuous(
      name = glue::glue("\nChromosome {chr} (Mb)"),
      limits = c(start, end),
      labels = scales::label_number(scale = 1e-6, drop0trailing = TRUE),
      breaks = seq(start, end, x_by),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      name = NULL,
      limits = c(0, pi_ymax),
      breaks = seq(0, pi_ymax, pi_y_by),
      labels = as.character(seq(0, pi_ymax, pi_y_by))
    ) +
    scale_colour_manual(values = palette_within) +
    guides(
      colour = guide_legend(
        ncol = 1,
        override.aes = list(linewidth = 3)
      )
    ) +
    facet_grid(vars(statistic)) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(face = "italic"),
      plot.margin=grid::unit(c(0,0,0,0), "mm"),
      strip.text = element_text(size = 16)
    )

  # -------------------------------- #
  # Return just the figure legends and exit
  if(return_legends) {
    legend_between <- ggpubr::get_legend(p_between)
    legend_within <- ggpubr::get_legend(p_within)
    return(list(
      ggpubr::as_ggplot(legend_between),
      ggpubr::as_ggplot(legend_within)
    ))
  } else {
    # -------------------------------- #
    # Plot: Centromere
    if(is.null(centromere)) {
      ( p_between / plot_spacer() / p_within) +
        plot_layout(
          axes = "collect_x",
          guides = "collect",
          ncol = 1,
          heights = c( 2, 0.05, 1)
        ) &
        theme(
          axis.title.x = element_text(size = 16),
          axis.text = element_text(size = 14),
          strip.text = element_text(size = 16),
          legend.title = element_blank(),
          strip.background = element_blank(),
          legend.position = "none",
          plot.margin = margin(0,0,0,0,"pt")
        )
    } else {
      p_centro <- centromere |>
        filter(chromosome == chr, window_pos_1 >= start, window_pos_2 <= end) |>
        ggplot( aes( x = floor((window_pos_1 + window_pos_2)/2), y = TR ) ) +
        geom_line() +
        annotate(
          "rect",
          xmin = sstart, xmax = send,
          ymin = -Inf, ymax = Inf,
          alpha = .4,
          fill = "grey60"
        ) +
        scale_x_continuous(
          name = glue::glue("\nChromosome {chr} (Mb)"),
          limits = c(start, end),
          labels = scales::label_number(scale = 1e-6),
          breaks = seq(start, end, x_by),
          expand = c(0, 0)
        ) +
        scale_y_continuous(
          name = NULL,
          labels = scales::number_format(scale = 1e-3),
          breaks = seq(0, 1e5, 5e4),
          limits = c(0, 1e5)
        ) +
        facet_grid(vars(statistic)) +
        theme_classic() +
        theme(
          legend.position = "bottom",
          plot.margin=grid::unit(c(0,0,0,0), "mm"),
          strip.text = element_blank()
        )

      if(!is.null(tag)) {
        p_centro <- p_centro +
          labs(tag = tag) +
          theme(plot.tag = element_text(size = 16, face = "bold"))
      }

      wrap_elements(
        (p_centro  / plot_spacer() / p_between / plot_spacer() / p_within) +
          plot_layout(
            axes = "collect_x",
            guides = "collect",
            ncol = 1,
            heights = c(0.6, 0.08, 2, 0.05, 1)
          ) &
          theme(
            axis.title.x = element_text(size = 16),
            axis.text = element_text(size = 14),
            legend.text = element_text(size = 14),
            legend.title = element_blank(),
            strip.background = element_blank(),
            legend.position = "none",
            plot.margin = margin(0,0,0,0,"pt")
          )
      )
    }
  }
}
