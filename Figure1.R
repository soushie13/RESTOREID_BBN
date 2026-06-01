# =========================================================
# COMPARATIVE FIGURE — FOUR DISEASE SYSTEMS
# Requires: ggplot2, patchwork, scales
install.packages("patchwork")
# =========================================================

library(ggplot2)
library(patchwork)
library(scales)


scen_levels <- c("none", "low", "moderate", "high")
scen_labels <- c("None", "Low", "Moderate", "High")

# ── Shared theme ──────────────────────────────────────────
theme_fig <- function(){
  theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor    = element_blank(),
      panel.grid.major.x  = element_blank(),
      axis.text           = element_text(size = 8, colour = "grey40"),
      axis.title.y        = element_text(size = 8, colour = "grey40"),
      axis.title.x        = element_blank(),
      legend.position     = "none",
      plot.tag            = element_text(size = 9, face = "bold",
                                         colour = "grey30"),
      plot.title          = element_text(size = 9, face = "bold"),
      plot.subtitle       = element_text(size = 7.5, colour = "grey50",
                                         face = "italic",
                                         margin = margin(b = 4))
    )
}

# ── Panel builder ─────────────────────────────────────────
make_panel <- function(df, pal, baseline_scen = "none",
                       title, subtitle,
                       ymax = 1.0, show_y = TRUE){
  
  baseline_val <- df$mean_disease[df$scenario == baseline_scen]
  
  ggplot(df,
         aes(x     = factor(scenario, levels = scen_levels,
                            labels = scen_labels),
             y     = mean_disease,
             group = 1)) +
    
    # 95% credible interval ribbon
    geom_ribbon(
      aes(ymin = lower_disease, ymax = upper_disease),
      fill  = pal[4], alpha = 0.13
    ) +
    
    # No-restoration baseline
    geom_hline(
      yintercept = baseline_val,
      linetype   = "dashed",
      colour     = "grey60",
      linewidth  = 0.45
    ) +
    
    # Mean line
    geom_line(linewidth = 1.2, colour = pal[4]) +
    
    # Scenario points
    geom_point(
      aes(fill = scenario),
      shape  = 21, size = 3.2,
      colour = "white", stroke = 1.2
    ) +
    
    scale_fill_manual(
      values = setNames(pal, scen_levels)
    ) +
    scale_x_discrete(labels = scen_labels) +
    scale_y_continuous(
      limits = c(0, ymax),
      breaks = seq(0, ymax, 0.2),
      labels = label_number(accuracy = 0.1),
      expand = expansion(mult = c(0.02, 0.06))
    ) +
    labs(
      title    = title,
      subtitle = subtitle,
      y        = if(show_y) "Disease risk probability" else NULL
    ) +
    theme_fig()
}

# ── Build four panels ─────────────────────────────────────

# A – Mosquito-borne (strong hump)
pal_mosq <- c(
  none     = "#aaaaaa",
  low      = "#bdd7e7",
  moderate = "#6baed6",
  high     = "#2171b5"
)

p_mosq <- make_panel(
  df       = mosq_summary,
  pal      = pal_mosq,
  title    = "Mosquito-borne",
  subtitle = "Standing water · vector habitat · ecological regulation",
  show_y   = TRUE
) + labs(tag = "A")

# B – Tick-borne (small hump)
pal_tick <- c("#aaaaaa","#bae4b3","#74c476","#238b45")

p_tick <- make_panel(
  df       = tick_summary,
  pal      = pal_tick,
  title    = "Tick-borne",
  subtitle = "Host connectivity before predator recovery · trophic lag",
  show_y   = FALSE
) + labs(tag = "B")

# C – Rodent-borne (steep decline)
pal_rod <- c("#aaaaaa","#fcae91","#fb6a4a","#cb181d")

p_rod <- make_panel(
  df       = rodent_summary,
  pal      = pal_rod,
  title    = "Rodent-borne",
  subtitle = "Habitat quality · biodiversity · rodent suppression",
  show_y   = TRUE
) + labs(tag = "C")

# D – Direct-contact zoonoses 
pal_zoon <- c(
  none     = "#aaaaaa",
  low      = "#cbc9e2",
  moderate = "#9e9ac8",
  high     = "#6a51a3"
)

p_zoon <- make_panel(
  df       = zoonotic_summary,
  pal      = pal_zoon,
  title    = "Direct-contact zoonoses",
  subtitle = "Fragmentation · human–wildlife interfaces · host diversity",
  show_y   = FALSE
) + labs(tag = "D")

# ── Compose with patchwork ────────────────────────────────

fig <- (p_mosq | p_tick) / (p_rod | p_zoon) +
  plot_annotation(
    caption = paste0(
      "Shaded bands: 95% Monte Carlo credible intervals (n = 1,000 iterations). ",
      "Dashed lines: no-restoration baseline risk. ",
      "All panels conditioned on high-risk region context."
    ),
    theme = theme(
      plot.caption = element_text(size = 7, colour = "grey50",
                                  hjust = 0,
                                  margin = margin(t = 8))
    )
  )

# ── Save ──────────────────────────────────────────────────

ggsave(
  "fig1_comparative_disease_systems.pdf",
  plot   = fig,
  width  = 174,     # mm — standard 2-column journal width
  height = 130,
  units  = "mm",
  device = cairo_pdf
)

ggsave(
  "fig1_comparative_disease_systems.png",
  plot   = fig,
  width  = 174,
  height = 130,
  units  = "mm",
  dpi    = 600
)
