# =========================================================
# MOSQUITO-BORNE DISEASE MODEL
# Bayesian Network with Monte Carlo uncertainty
#
# Expected pattern:
#   Risk        : strong hump (low–moderate restoration peaks,
#                 high restoration ultimately falls BELOW no-restoration)
#   Stabilisation: moderately delayed
#   Uncertainty : moderate
# =========================================================

library(bnlearn)
library(gRain)
library(dplyr)
library(ggplot2)
library(tidyr)

set.seed(123)

# =========================================================
# STATES
# =========================================================

states <- list(
  
  Region_Context =
    c("Low","Moderate","High"),
  
  Restoration_Intensity =
    c("none","low","moderate","high"),
  
  Hydrological_Recovery =
    c("low","moderate","high"),
  
  Standing_Water =
    c("low","moderate","high"),
  
  Mosquito_Disease_Risk =
    c("absent","present"),
  
  Ecological_Regulation =
    c("weak","moderate","strong"),
  
  Time_to_Stable_State =
    c("1_3_years","4_7_years","8_12_years","13_15_years")
)

# =========================================================
# DAG  (unchanged)
# =========================================================
dag <- model2network(
  paste0(
    "[Region_Context]",
    "[Restoration_Intensity]",
    "[Hydrological_Recovery|Region_Context:Restoration_Intensity]",
    "[Standing_Water|Hydrological_Recovery:Restoration_Intensity]",
    "[Ecological_Regulation|Restoration_Intensity]",
    "[Mosquito_Disease_Risk|Standing_Water:Ecological_Regulation]",
    "[Time_to_Stable_State|Ecological_Regulation]"
  ))

# =========================================================
# HELPER
# =========================================================
sample_prob <- function(mn, mx) runif(1, mn, mx)

# =========================================================
# MONTE CARLO
# =========================================================
n_iter  <- 1000
results <- data.frame()

for(i in 1:n_iter){
  
  # -------------------------------------------------------
  # REGION (flat prior)
  # -------------------------------------------------------
  cpt_region <- array(
    c(0.33, 0.34, 0.33),
    dim      = c(3),
    dimnames = list(Region_Context = states$Region_Context)
  )
  
  # -------------------------------------------------------
  # RESTORATION (flat prior)
  # -------------------------------------------------------
  cpt_restoration <- array(
    c(0.25, 0.25, 0.25, 0.25),
    dim      = c(4),
    dimnames = list(Restoration_Intensity = states$Restoration_Intensity)
  )
  
  # -------------------------------------------------------
  # HYDROLOGICAL RECOVERY
  # Higher restoration → higher probability of hydrological
  # recovery, regardless of region context.
  # -------------------------------------------------------
  hydro_vals <- c()
  
  for(region in 1:3){
    for(rest in 1:4){
      
      probs <- switch(rest,
                      
                      # none – mostly low recovery
                      `1` = c(sample_prob(0.60,0.80),
                              sample_prob(0.15,0.30),
                              sample_prob(0.05,0.15)),
                      
                      # low – mixed
                      `2` = c(sample_prob(0.20,0.35),
                              sample_prob(0.35,0.50),
                              sample_prob(0.25,0.40)),
                      
                      # moderate – leaning high
                      `3` = c(sample_prob(0.05,0.20),
                              sample_prob(0.25,0.40),
                              sample_prob(0.45,0.65)),
                      
                      # high – strongly high recovery
                      `4` = c(sample_prob(0.05,0.15),
                              sample_prob(0.20,0.35),
                              sample_prob(0.55,0.75))
      )
      
      probs      <- probs / sum(probs)
      hydro_vals <- c(hydro_vals, probs)
    }
  }
  
  cpt_hydro <- array(
    hydro_vals,
    dim      = c(3,4,3),
    dimnames = list(
      Hydrological_Recovery = states$Hydrological_Recovery,
      Restoration_Intensity = states$Restoration_Intensity,
      Region_Context        = states$Region_Context
    )
  )
  
  # -------------------------------------------------------
  # STANDING WATER
  #
  # KEY CHANGE – encode the hump mechanism explicitly:
  #
  #   none     → moderate standing water (baseline, unmanaged)
  #   low      → standing water rises (early drainage disruption,
  #               partial re-wetting without vegetation control)
  #   moderate → PEAK standing water (maximum re-wetting before
  #               mature vegetation / drainage establishes)
  #   high     → standing water returns toward baseline or lower
  #               (hydrological equilibrium restored, mature
  #               wetland absorbs/routes water efficiently)
  #
  # Hydrological recovery modifies the shape but does not
  # override the restoration-driven hump.
  # -------------------------------------------------------
  # Array dim order: [Standing_Water, Hydrological_Recovery, Restoration_Intensity]
  # R fills with dim 1 (Standing_Water) fastest, dim 3 (Restoration_Intensity) slowest.
  # Loop order must be: Restoration_Intensity outer (dim 3), Hydrological_Recovery inner (dim 2).
  
  water_vals <- c()
  
  for(rest in 1:4){           # dim 3 - Restoration_Intensity (outer = slowest)
    for(hydro in 1:3){        # dim 2 - Hydrological_Recovery
      
      probs <- switch(rest,
                      
                      # none  - moderate baseline standing water
                      `1` = c(sample_prob(0.25,0.40),   # low
                              sample_prob(0.40,0.55),   # moderate  <- dominant
                              sample_prob(0.10,0.20)),  # high
                      
                      # low  - re-wetting begins, standing water increases
                      `2` = c(sample_prob(0.10,0.20),   # low
                              sample_prob(0.35,0.50),   # moderate
                              sample_prob(0.35,0.50)),  # high       <- elevated
                      
                      # moderate - PEAK: maximum standing water
                      `3` = c(sample_prob(0.02,0.08),   # low
                              sample_prob(0.15,0.28),   # moderate
                              sample_prob(0.65,0.82)),  # high       <- peak
                      
                      # high - hydrological equilibrium, mature vegetation
                      # HIGH standing water suppressed BELOW 'none' baseline
                      `4` = c(sample_prob(0.35,0.55),   # low        <- dominant
                              sample_prob(0.30,0.45),   # moderate
                              sample_prob(0.05,0.15))   # high       <- suppressed
      )
      
      # Hydrological recovery: small secondary modifier
      if(hydro == 1) probs[1] <- probs[1] + 0.10   # low hydro -> more low water
      if(hydro == 3) probs[3] <- probs[3] + 0.08   # high hydro -> slightly more high water
      
      probs      <- probs / sum(probs)
      water_vals <- c(water_vals, probs)
    }
  }
  
  cpt_water <- array(
    water_vals,
    dim      = c(3,3,4),
    dimnames = list(
      Standing_Water        = states$Standing_Water,
      Hydrological_Recovery = states$Hydrological_Recovery,
      Restoration_Intensity = states$Restoration_Intensity
    )
  )
  
  # -------------------------------------------------------
  # ECOLOGICAL REGULATION
  # (unchanged from original – high restoration → strong regulation)
  # -------------------------------------------------------
  reg_vals <- c(
    0.80, 0.15, 0.05,   # none    → mostly weak
    0.60, 0.30, 0.10,   # low     → mostly weak
    0.25, 0.45, 0.30,   # moderate→ moderate
    0.05, 0.20, 0.75    # high    → mostly strong
  )
  
  reg_vals <- reg_vals + runif(length(reg_vals), -0.03, 0.03)
  reg_vals[reg_vals < 0.01] <- 0.01
  
  reg_vals <- unlist(lapply(
    split(reg_vals, ceiling(seq_along(reg_vals)/3)),
    function(x) x / sum(x)
  ))
  
  cpt_regulation <- array(
    reg_vals,
    dim      = c(3,4),
    dimnames = list(
      Ecological_Regulation = states$Ecological_Regulation,
      Restoration_Intensity = states$Restoration_Intensity
    )
  )
  
  # -------------------------------------------------------
  # DISEASE RISK
  #
  # KEY CHANGE – the long-run benefit of high restoration
  # must push disease risk BELOW the no-restoration baseline.
  #
  # Mechanism: high restoration → strong ecological regulation
  #            (predators, competitors, mature wetland structure)
  #            → disease risk suppressed even when some standing
  #            water is present.
  #
  # Rule of thumb for the 9 cells (water × regulation):
  #
  #   high water  + weak reg   → very high risk  (0.88–0.95 present)
  #   high water  + mod  reg   → high risk        (0.60–0.75 present)
  #   high water  + strong reg → LOW risk          (0.08–0.18 present)
  #   mod  water  + weak reg   → moderate-high    (0.55–0.70 present)
  #   mod  water  + mod  reg   → moderate          (0.35–0.55 present)
  #   mod  water  + strong reg → low               (0.10–0.22 present)
  #   low  water  + weak reg   → moderate          (0.35–0.50 present)
  #   low  water  + mod  reg   → low-moderate      (0.18–0.30 present)
  #   low  water  + strong reg → very low          (0.04–0.12 present)
  #
  # The "below baseline" effect emerges from the combination of
  # low standing water AND strong ecological regulation under high
  # restoration, both of which are individually low-probability
  # under no-restoration.
  # -------------------------------------------------------
  # Array dim order: [Mosquito_Disease_Risk, Standing_Water, Ecological_Regulation]
  # R fills arrays with the FIRST dimension varying fastest.
  # So the loop filling disease_vals must step through dimensions
  # from LAST (slowest) to FIRST (fastest):
  #   outer loop  → Ecological_Regulation  (dim 3, slowest)
  #   middle loop → Standing_Water         (dim 2)
  #   inner        → disease pair          (dim 1, fastest — just c(absent, present))
  # This ensures cell [disease, water, reg] receives the correct probabilities.
  
  disease_vals <- c()
  
  for(reg in 1:3){          # dim 3 – Ecological_Regulation (outer = slowest)
    for(water in 1:3){      # dim 2 – Standing_Water
      
      probs <- if        (water == 3 & reg == 1){
        c(sample_prob(0.05, 0.12), sample_prob(0.88, 0.95))   # high water + weak reg   → near-certain risk
        
      } else if(water == 3 & reg == 2){
        c(sample_prob(0.25, 0.40), sample_prob(0.60, 0.75))   # high water + mod reg    → high risk
        
      } else if(water == 3 & reg == 3){
        # Strong regulation suppresses risk even with high standing water
        c(sample_prob(0.82, 0.92), sample_prob(0.08, 0.18))   # high water + strong reg → LOW risk
        
      } else if(water == 2 & reg == 1){
        c(sample_prob(0.30, 0.45), sample_prob(0.55, 0.70))   # mod water  + weak reg   → moderate-high
        
      } else if(water == 2 & reg == 2){
        c(sample_prob(0.45, 0.65), sample_prob(0.35, 0.55))   # mod water  + mod reg    → moderate
        
      } else if(water == 2 & reg == 3){
        c(sample_prob(0.78, 0.90), sample_prob(0.10, 0.22))   # mod water  + strong reg → low
        
      } else if(water == 1 & reg == 1){
        c(sample_prob(0.50, 0.65), sample_prob(0.35, 0.50))   # low water  + weak reg   → moderate
        
      } else if(water == 1 & reg == 2){
        c(sample_prob(0.70, 0.82), sample_prob(0.18, 0.30))   # low water  + mod reg    → low-moderate
        
      } else {
        # low water + strong regulation → absolute floor of risk
        c(sample_prob(0.88, 0.96), sample_prob(0.04, 0.12))   # low water  + strong reg → very low
      }
      
      probs        <- probs / sum(probs)
      disease_vals <- c(disease_vals, probs)
    }
  }
  
  cpt_disease <- array(
    disease_vals,
    dim      = c(2,3,3),
    dimnames = list(
      Mosquito_Disease_Risk = states$Mosquito_Disease_Risk,
      Standing_Water        = states$Standing_Water,
      Ecological_Regulation = states$Ecological_Regulation
    )
  )
  
  # -------------------------------------------------------
  # TIME TO STABLE STATE
  #
  # KEY CHANGE – link temporal delay to restoration intensity,
  # not just ecological regulation.
  #
  # The hump implies a DELAY before risk declines; this should
  # be most pronounced at moderate restoration (peak risk,
  # slow stabilisation) and fastest at high restoration.
  #
  # Weak regulation (none/low restoration)  → quick plateau at
  #   HIGH risk (not a delay toward safety, just early lock-in)
  # Moderate regulation (moderate rest.)    → delayed stabilisation,
  #   risk lingers in 4–12 yr window
  # Strong regulation (high restoration)   → slower initial
  #   trajectory but reaches stable LOW risk by 8–15 yr
  # -------------------------------------------------------
  time_vals <- c()
  
  for(reg in 1:3){
    
    probs <- switch(reg,
                    
                    # weak regulation (none/low restoration)
                    # → stabilises quickly, but AT high risk
                    `1` = c(sample_prob(0.40, 0.55),   # 1-3 yr   ← dominant
                            sample_prob(0.25, 0.35),   # 4-7 yr
                            sample_prob(0.10, 0.18),   # 8-12 yr
                            sample_prob(0.04, 0.10)),  # 13-15 yr
                    
                    # moderate regulation (moderate restoration)
                    # → delayed; peak of stabilisation in middle windows
                    `2` = c(sample_prob(0.10, 0.20),   # 1-3 yr
                            sample_prob(0.25, 0.35),   # 4-7 yr   ← elevated
                            sample_prob(0.30, 0.42),   # 8-12 yr  ← elevated
                            sample_prob(0.15, 0.28)),  # 13-15 yr
                    
                    # strong regulation (high restoration)
                    # → takes longest but achieves true stable low-risk state
                    `3` = c(sample_prob(0.04, 0.10),   # 1-3 yr
                            sample_prob(0.12, 0.22),   # 4-7 yr
                            sample_prob(0.32, 0.45),   # 8-12 yr
                            sample_prob(0.30, 0.48))   # 13-15 yr ← dominant
    )
    
    probs     <- probs / sum(probs)
    time_vals <- c(time_vals, probs)
  }
  
  cpt_time <- array(
    time_vals,
    dim      = c(4,3),
    dimnames = list(
      Time_to_Stable_State  = states$Time_to_Stable_State,
      Ecological_Regulation = states$Ecological_Regulation
    )
  )
  
  # -------------------------------------------------------
  # FIT NETWORK
  # -------------------------------------------------------
  fit <- custom.fit(
    dag,
    dist = list(
      Region_Context        = cpt_region,
      Restoration_Intensity = cpt_restoration,
      Hydrological_Recovery = cpt_hydro,
      Standing_Water        = cpt_water,
      Mosquito_Disease_Risk = cpt_disease,
      Ecological_Regulation = cpt_regulation,
      Time_to_Stable_State  = cpt_time
    )
  )
  
  bn <- as.grain(fit)
  
  # -------------------------------------------------------
  # SCENARIO QUERIES  (conditioned on High-risk region)
  # -------------------------------------------------------
  for(rest in states$Restoration_Intensity){
    
    bn_temp <- setEvidence(
      bn,
      nodes  = c("Region_Context","Restoration_Intensity"),
      states = c("High", rest)
    )
    
    risk       <- querygrain(bn_temp, nodes = "Mosquito_Disease_Risk")$Mosquito_Disease_Risk
    regulation <- querygrain(bn_temp, nodes = "Ecological_Regulation")$Ecological_Regulation
    time       <- querygrain(bn_temp, nodes = "Time_to_Stable_State")$Time_to_Stable_State
    
    results <- bind_rows(results, data.frame(
      iteration       = i,
      scenario        = rest,
      disease_risk    = risk["present"],
      regulation_strong = regulation["strong"],
      years_1_3       = time["1_3_years"],
      years_4_7       = time["4_7_years"],
      years_8_12      = time["8_12_years"],
      years_13_15     = time["13_15_years"]
    ))
  }
  
}  # END MONTE CARLO

# =========================================================
# SUMMARY
# =========================================================
results$long_term_stability <-
  results$years_8_12 + results$years_13_15

summary_df <- results %>%
  group_by(scenario) %>%
  summarise(
    mean_disease      = mean(disease_risk),
    lower_disease     = quantile(disease_risk, 0.025),
    upper_disease     = quantile(disease_risk, 0.975),
    mean_regulation   = mean(regulation_strong),
    lower_regulation  = quantile(regulation_strong, 0.025),
    upper_regulation  = quantile(regulation_strong, 0.975),
    mean_stability    = mean(long_term_stability),
    lower_stability   = quantile(long_term_stability, 0.025),
    upper_stability   = quantile(long_term_stability, 0.975),
    .groups = "drop"
  )

print(summary_df)

# =========================================================
# COLOUR PALETTE
# =========================================================
mosquito_cols <- c(
  none     = "#aaaaaa",
  low      = "#bdd7e7",
  moderate = "#6baed6",
  high     = "#2171b5"
)

scen_levels <- c("none","low","moderate","high")

# =========================================================
# PLOT 1 – Disease Risk Hump
# =========================================================
ggplot(
  summary_df,
  aes(x = factor(scenario, levels = scen_levels),
      y = mean_disease, group = 1)
)+
  geom_line(linewidth = 1.4)+
  geom_point(aes(fill = scenario), shape = 21, size = 5)+
  geom_errorbar(aes(ymin = lower_disease, ymax = upper_disease),
                width = 0.15)+
  scale_fill_manual(values = mosquito_cols)+
  theme_minimal(base_size = 14)+
  labs(
    title = "Mosquito-Borne Disease Risk Across Restoration Intensities",
    x = "Restoration Intensity",
    y = "Disease Risk Probability"
  )

# =========================================================
# PLOT 2 – Ecological Regulation Recovery
# =========================================================
ggplot(
  summary_df,
  aes(x = factor(scenario, levels = scen_levels),
      y = mean_regulation, group = 1)
)+
  geom_line(linewidth = 1.4)+
  geom_point(aes(fill = scenario), shape = 21, size = 5)+
  geom_errorbar(aes(ymin = lower_regulation, ymax = upper_regulation),
                width = 0.15)+
  scale_fill_manual(values = mosquito_cols)+
  theme_minimal(base_size = 14)+
  labs(
    title = "Ecological Regulation Recovery",
    x = "Restoration Intensity",
    y = "Probability of Strong Regulation"
  )

# =========================================================
# PLOT 3 – Long-Term Ecological Stabilisation
# =========================================================
ggplot(
  summary_df,
  aes(x = factor(scenario, levels = scen_levels),
      y = mean_stability, group = 1)
)+
  geom_line(linewidth = 1.4)+
  geom_point(aes(fill = scenario), shape = 21, size = 5)+
  geom_errorbar(aes(ymin = lower_stability, ymax = upper_stability),
                width = 0.15)+
  scale_fill_manual(values = mosquito_cols)+
  theme_minimal(base_size = 14)+
  labs(
    title = "Long-Term Ecological Stabilisation (8–15 yr)",
    x = "Restoration Intensity",
    y = "Probability of Stabilisation in Long Window"
  )

# =========================================================
# PLOT 4 – Disease Risk vs Long-Term Stability (scatter)
# =========================================================
ggplot(
  results,
  aes(x = disease_risk, y = long_term_stability, colour = scenario)
)+
  geom_point(alpha = 0.20)+
  geom_smooth(method = "loess", se = TRUE)+
  scale_colour_manual(values = mosquito_cols)+
  theme_minimal(base_size = 14)+
  labs(
    title = "Mosquito Disease Risk vs Long-Term Stability",
    x = "Disease Risk",
    y = "Long-Term Stability"
  )

# =========================================================
# PLOT 5 – Stabilisation Trajectories over Time
# =========================================================
traj_long <- results %>%
  pivot_longer(
    cols      = c(years_1_3, years_4_7, years_8_12, years_13_15),
    names_to  = "time_period",
    values_to = "probability"
  )

traj_long$time_period <- factor(
  traj_long$time_period,
  levels = c("years_1_3","years_4_7","years_8_12","years_13_15"),
  labels = c("1–3 yr","4–7 yr","8–12 yr","13–15 yr")
)

corridor <- traj_long %>%
  group_by(scenario, time_period) %>%
  summarise(
    mean_prob = mean(probability),
    lower     = quantile(probability, 0.025),
    upper     = quantile(probability, 0.975),
    .groups   = "drop"
  )

ggplot(
  corridor,
  aes(x = time_period, y = mean_prob,
      group = scenario, colour = scenario, fill = scenario)
)+
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.20, colour = NA)+
  geom_line(linewidth = 1.4)+
  geom_point(size = 3)+
  scale_colour_manual(values = mosquito_cols)+
  scale_fill_manual(values = mosquito_cols)+
  theme_minimal(base_size = 14)+
  labs(
    title    = "Stabilisation Trajectories by Restoration Intensity",
    x = "Time to Stabilisation",
    y = "Probability"
  )

# =========================================================
# HUMP PLOT – Disease risk with uncertainty ribbon,
# baseline reference, and peak / below-baseline annotations
# =========================================================


mosquito_cols <- c(
  none     = "#aaaaaa",
  low      = "#bdd7e7",
  moderate = "#6baed6",
  high     = "#2171b5"
)

scen_levels <- c("none", "low", "moderate", "high")

# ── baseline value (no-restoration mean risk) ─────────────
baseline_val <- summary_df$mean_disease[summary_df$scenario == "none"]

# ── peak annotation position ──────────────────────────────
peak_row <- summary_df[summary_df$scenario == "moderate", ]

# ── below-baseline annotation position ───────────────────
high_row  <- summary_df[summary_df$scenario == "high", ]

ggplot(
  summary_df,
  aes(
    x     = factor(scenario, levels = scen_levels),
    y     = mean_disease,
    group = 1
  )
) +
  
  # shaded 95% credible interval ribbon
  geom_ribbon(
    aes(ymin = lower_disease, ymax = upper_disease),
    fill  = "#2171b5",
    alpha = 0.12
  ) +
  
  # dashed baseline reference line
  geom_hline(
    yintercept = baseline_val,
    linetype   = "dashed",
    colour     = "grey50",
    linewidth  = 0.6
  ) +
  
  # annotate the dashed line
  annotate(
    "text",
    x     = 0.55,
    y     = baseline_val + 0.012,
    label = "No-restoration baseline",
    hjust = 0,
    size  = 3.2,
    colour = "grey45"
  ) +
  
  # main line
  geom_line(
    linewidth = 1.5,
    colour    = "#2171b5"
  ) +
  
  # coloured dots per scenario
  geom_point(
    aes(fill = scenario),
    shape = 21,
    size  = 5,
    colour = "white",
    stroke = 1.8
  ) +
  
  # ── peak annotation ───────────────────────────────────────
  annotate(
    "segment",
    x    = peak_row$scenario,  xend = peak_row$scenario,
    y    = peak_row$upper_disease + 0.005,
    yend = peak_row$upper_disease + 0.025,
    colour    = "#6baed6",
    linewidth = 0.5
  ) +
  annotate(
    "label",
    x         = peak_row$scenario,
    y         = peak_row$upper_disease + 0.030,
    label     = paste0("Peak risk\n", round(peak_row$mean_disease, 2)),
    fill      = "#bdd7e7",
    colour    = "#6baed6",
    size      = 3.2,
    label.size = 0,
    fontface  = "bold"
  ) +
  
  # ── below-baseline annotation ─────────────────────────────
  annotate(
    "segment",
    x    = high_row$scenario,  xend = high_row$scenario,
    y    = high_row$lower_disease - 0.005,
    yend = high_row$lower_disease - 0.025,
    colour    = "#2171b5",
    linewidth = 0.5
  ) +
  annotate(
    "label",
    x         = high_row$scenario,
    y         = high_row$lower_disease - 0.030,
    label     = paste0("Below baseline\n", round(high_row$mean_disease, 2)),
    fill      = "#eeedfe",
    colour    = "#2171b5",
    size      = 3.2,
    label.size = 0,
    fontface  = "bold"
  ) +
  
  scale_fill_manual(values = mosquito_cols) +
  scale_y_continuous(
    limits = c(0, 0.60),
    breaks = seq(0, 0.6, 0.1),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position  = "none",
    panel.grid.minor = element_blank(),
    plot.subtitle    = element_text(colour = "grey50", size = 11)
  ) +
  labs(
    title    = "Mosquito-borne disease risk across restoration intensities",
    x        = "Restoration intensity",
    y        = "Disease risk probability"
  )
