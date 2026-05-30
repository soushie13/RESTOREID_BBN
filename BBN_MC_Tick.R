# =========================================================
# TICK-BORNE DISEASE SYSTEM
# MONTE CARLO BAYESIAN BELIEF NETWORK
#
# Mechanism:
#   Restoration → increased vegetation complexity
#               → increased host abundance
#               → increased host connectivity  ← hump driver
#               BEFORE predator recovery and
#               trophic regulation kick in
#
# Expected pattern:
#   Risk          : small hump at low restoration,
#                   then decline as predators recover
#   Stabilisation : very delayed (trophic cascades are slow)
#   Uncertainty   : moderate
#
# =========================================================


library(bnlearn)
library(gRain)
library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(123)

# =========================================================
# HELPER
# =========================================================

sample_prob <- function(mn, mx) runif(1, mn, mx)

# =========================================================
# SETTINGS
# =========================================================

n_iter  <- 1000
results <- data.frame()

# =========================================================
# STATES
# =========================================================

states <- list(
  
  Region_Context = c(
    "Low",
    "Moderate",
    "High"
  ),
  
  Restoration_Intensity = c(
    "none",
    "low",
    "moderate",
    "high"
  ),
  
  Host_Connectivity = c(
    "low",
    "moderate",
    "high"
  ),
  
  Predator_Recovery = c(
    "weak",
    "moderate",
    "strong"
  ),
  
  Tick_Disease_Risk = c(
    "absent",
    "present"
  ),
  
  Ecological_Regulation = c(
    "weak",
    "moderate",
    "strong"
  ),
  
  Time_to_Stable_State = c(
    "1_3_years",
    "4_7_years",
    "8_12_years",
    "13_15_years"
  )
)

# =========================================================
# DAG  (unchanged)
# =========================================================

dag <- model2network(
  paste0(
    "[Region_Context]",
    "[Restoration_Intensity]",
    "[Host_Connectivity|Restoration_Intensity]",
    "[Predator_Recovery|Restoration_Intensity]",
    "[Tick_Disease_Risk|Host_Connectivity:Predator_Recovery]",
    "[Ecological_Regulation|Predator_Recovery:Tick_Disease_Risk]",
    "[Time_to_Stable_State|Ecological_Regulation]"
  )
)

# =========================================================
# MONTE CARLO LOOP
# =========================================================

for(i in 1:n_iter){
  
  # -------------------------------------------------------
  # PRIORS (flat)
  # -------------------------------------------------------
  
  cpt_region <- array(
    c(0.33, 0.34, 0.33),
    dim      = c(3),
    dimnames = list(Region_Context = states$Region_Context)
  )
  
  cpt_restoration <- array(
    c(0.25, 0.25, 0.25, 0.25),
    dim      = c(4),
    dimnames = list(Restoration_Intensity = states$Restoration_Intensity)
  )
  
  # -------------------------------------------------------
  # HOST CONNECTIVITY
  # dim = c(3,4) = [Connectivity, Rest]
  # Single parent — flat vector fills correctly (Rest slowest).
  #
  # HUMP SHAPE — the core mechanism:
  #
  #   none     → LOW connectivity
  #              Fragmented landscape, low host density,
  #              few suitable tick habitats (open degraded land)
  #
  #   low      → PEAK connectivity  ← hump maximum
  #              Vegetation regenerates, host abundance rises
  #              (deer, rodents colonise new scrub-edge habitat),
  #              patch connectivity increases — but predators
  #              have NOT yet recovered, so no top-down control.
  #              This is the critical window of elevated risk.
  #
  #   moderate → still high connectivity but predators beginning
  #              to establish; host density starts to be checked
  #
  #   high     → moderate connectivity
  #              Mature vegetation reduces tick-suitable edge
  #              habitat; trophic regulation suppresses host
  #              density at landscape scale; structural complexity
  #              routes hosts away from human interfaces
  # -------------------------------------------------------
  
  connectivity_vals <- c(
    # none – mostly low
    sample_prob(0.55, 0.72), sample_prob(0.18, 0.30), sample_prob(0.05, 0.15),
    # low  – PEAK: strongly skewed to high connectivity
    sample_prob(0.05, 0.12), sample_prob(0.18, 0.30), sample_prob(0.60, 0.78),
    # moderate – still high but slightly pulled back
    sample_prob(0.08, 0.18), sample_prob(0.22, 0.35), sample_prob(0.50, 0.68),
    # high – back to moderate; trophic regulation reducing hosts
    sample_prob(0.15, 0.30), sample_prob(0.38, 0.52), sample_prob(0.25, 0.42)
  )
  
  # Normalise each restoration column (3 values each)
  connectivity_vals <- unlist(lapply(
    split(connectivity_vals, ceiling(seq_along(connectivity_vals)/3)),
    function(x) x / sum(x)
  ))
  
  cpt_connectivity <- array(
    connectivity_vals,
    dim      = c(3,4),
    dimnames = list(
      Host_Connectivity     = states$Host_Connectivity,
      Restoration_Intensity = states$Restoration_Intensity
    )
  )
  
  # -------------------------------------------------------
  # PREDATOR RECOVERY
  # dim = c(3,4) = [Predator, Rest]
  # Single parent — flat vector fills correctly.
  #
  # KEY: predators lag behind vegetation and host recovery.
  # At low/moderate restoration the predator signal is still
  # mostly weak — this is what allows the hump to persist.
  # Only at high restoration does strong predator recovery
  # become probable, and even then uncertainty is moderate.
  # -------------------------------------------------------
  
  predator_vals <- c(
    # none     – almost no predator recovery
    0.82, 0.14, 0.04,
    # low      – still mostly weak; predators haven't tracked yet
    0.72, 0.22, 0.06,
    # moderate – beginning to shift; weak still dominant
    0.52, 0.33, 0.15,
    # high     – strong recovery likely but uncertain (slow cascade)
    0.12, 0.28, 0.60
  )
  
  predator_vals <- predator_vals +
    runif(length(predator_vals), -0.04, 0.04)
  predator_vals[predator_vals < 0.01] <- 0.01
  
  predator_vals <- unlist(lapply(
    split(predator_vals, ceiling(seq_along(predator_vals)/3)),
    function(x) x / sum(x)
  ))
  
  cpt_predators <- array(
    predator_vals,
    dim      = c(3,4),
    dimnames = list(
      Predator_Recovery     = states$Predator_Recovery,
      Restoration_Intensity = states$Restoration_Intensity
    )
  )
  
  # -------------------------------------------------------
  # TICK DISEASE RISK
  # dim = c(2,3,3) = [Disease, Connectivity, Predator]
  # Loop: Predator outer (dim3, slowest), Connectivity inner (dim2)
  #
  # Risk is jointly determined by connectivity (host exposure)
  # and predator suppression (top-down control).
  # The hump emerges because low restoration → high connectivity
  # + weak predators → peak risk.
  # High restoration → moderate connectivity + strong predators
  # → risk suppressed below the no-restoration baseline.
  #
  # Cell design:
  #   high connectivity + weak predator   → PEAK risk (hump top)
  #   high connectivity + strong predator → LOW risk  (recovery)
  #   low  connectivity + weak predator   → moderate  (baseline)
  #   low  connectivity + strong predator → very low  (floor)
  #   all middle combinations             → broad uncertainty
  # -------------------------------------------------------
  
  risk_vals <- c()
  
  for(predator in 1:3){         # dim 3 – Predator (outer)
    for(connectivity in 1:3){   # dim 2 – Connectivity (inner)
      
      probs <- if(connectivity == 3 & predator == 1){
        # PEAK: high connectivity + weak predator
        # This is the low-restoration hump maximum
        c(sample_prob(0.05, 0.14),
          sample_prob(0.86, 0.95))
        
      } else if(connectivity == 3 & predator == 3){
        # Trophic suppression: high connectivity but strong predators
        # Risk falls well below no-restoration baseline
        c(sample_prob(0.68, 0.82),
          sample_prob(0.18, 0.32))
        
      } else if(connectivity == 1 & predator == 1){
        # No restoration baseline: low connectivity, weak predators
        # Moderate risk — few hosts but no control either
        c(sample_prob(0.42, 0.60),
          sample_prob(0.40, 0.58))
        
      } else if(connectivity == 1 & predator == 3){
        # Low connectivity + strong predators: very low risk floor
        c(sample_prob(0.75, 0.90),
          sample_prob(0.10, 0.25))
        
      } else if(connectivity == 2 & predator == 2){
        # Mid-transition: maximum uncertainty
        c(sample_prob(0.30, 0.70),
          sample_prob(0.30, 0.70))
        
      } else {
        # All remaining combinations: moderate, broad uncertainty
        c(sample_prob(0.28, 0.68),
          sample_prob(0.32, 0.72))
      }
      
      probs     <- probs / sum(probs)
      risk_vals <- c(risk_vals, probs)
    }
  }
  
  cpt_disease <- array(
    risk_vals,
    dim      = c(2,3,3),
    dimnames = list(
      Tick_Disease_Risk = states$Tick_Disease_Risk,
      Host_Connectivity = states$Host_Connectivity,
      Predator_Recovery = states$Predator_Recovery
    )
  )
  
  # -------------------------------------------------------
  # ECOLOGICAL REGULATION
  # dim = c(3,3,2) = [Regulation, Predator, Disease]
  # Loop: Disease outer (dim3, slowest), Predator inner (dim2)
  #
  # Strong predators + absent disease → strong regulation
  # Weak predators   + present disease → weak regulation
  # All other combinations: broad uncertainty (moderate range)
  # -------------------------------------------------------
  
  regulation_vals <- c()
  
  for(disease in 1:2){          # dim 3 – Disease (outer)
    for(predator in 1:3){       # dim 2 – Predator (inner)
      
      probs <- if(predator == 3 & disease == 1){
        # Strong predators + disease absent: full trophic regulation
        c(sample_prob(0.04, 0.10),
          sample_prob(0.15, 0.26),
          sample_prob(0.66, 0.80))
        
      } else if(predator == 1 & disease == 2){
        # Weak predators + disease present: regulatory failure
        c(sample_prob(0.65, 0.80),
          sample_prob(0.15, 0.26),
          sample_prob(0.04, 0.10))
        
      } else {
        # Transitional: wide, overlapping distributions
        c(sample_prob(0.18, 0.48),
          sample_prob(0.22, 0.48),
          sample_prob(0.18, 0.48))
      }
      
      probs           <- probs / sum(probs)
      regulation_vals <- c(regulation_vals, probs)
    }
  }
  
  cpt_regulation <- array(
    regulation_vals,
    dim      = c(3,3,2),
    dimnames = list(
      Ecological_Regulation = states$Ecological_Regulation,
      Predator_Recovery     = states$Predator_Recovery,
      Tick_Disease_Risk     = states$Tick_Disease_Risk
    )
  )
  
  # -------------------------------------------------------
  # TIME TO STABLE STATE
  # dim = c(4,3) = [Time, Regulation]
  # Single parent — flat vector, Regulation is slowest (correct).
  #
  # VERY DELAYED stabilisation is the key signature:
  #
  # Weak regulation  → fast plateau but at degraded state
  #                    (mass in 1-3 yr window)
  #
  # Moderate         → broad, flat distribution; no clear peak
  #                    (trophic cascade still developing)
  #
  # Strong regulation → very delayed resilient equilibrium
  #                    Mass concentrated in 8-15 yr windows
  #                    This is much more delayed than the
  #                    mosquito system (which had 8-12 peak)
  #                    — trophic regulation takes longer to
  #                    stabilise than wetland hydrology
  # -------------------------------------------------------
  
  time_vals <- c()
  
  for(reg in 1:3){
    
    probs <- switch(reg,
                    
                    # weak – rapid degraded equilibrium
                    `1` = c(sample_prob(0.38, 0.52),   # 1-3 yr  ← dominant
                            sample_prob(0.28, 0.38),   # 4-7 yr
                            sample_prob(0.08, 0.18),   # 8-12 yr
                            sample_prob(0.03, 0.10)),  # 13-15 yr
                    
                    # moderate – flat, weakly predictable
                    `2` = c(sample_prob(0.15, 0.28),   # 1-3 yr
                            sample_prob(0.22, 0.32),   # 4-7 yr
                            sample_prob(0.24, 0.34),   # 8-12 yr
                            sample_prob(0.20, 0.32)),  # 13-15 yr
                    
                    # strong – very delayed; mass in long windows
                    # Trophic cascades stabilise over decades not years
                    `3` = c(sample_prob(0.02, 0.06),   # 1-3 yr
                            sample_prob(0.05, 0.12),   # 4-7 yr
                            sample_prob(0.22, 0.34),   # 8-12 yr
                            sample_prob(0.52, 0.68))   # 13-15 yr ← dominant
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
      Host_Connectivity     = cpt_connectivity,
      Predator_Recovery     = cpt_predators,
      Tick_Disease_Risk     = cpt_disease,
      Ecological_Regulation = cpt_regulation,
      Time_to_Stable_State  = cpt_time
    )
  )
  
  bn <- as.grain(fit)
  
  # -------------------------------------------------------
  # SCENARIO QUERIES  (High-risk region)
  # -------------------------------------------------------
  
  for(rest in states$Restoration_Intensity){
    
    bn_temp <- setEvidence(
      bn,
      nodes  = c("Region_Context", "Restoration_Intensity"),
      states = c("High", rest)
    )
    
    risk       <- querygrain(bn_temp, nodes = "Tick_Disease_Risk")$Tick_Disease_Risk
    regulation <- querygrain(bn_temp, nodes = "Ecological_Regulation")$Ecological_Regulation
    time       <- querygrain(bn_temp, nodes = "Time_to_Stable_State")$Time_to_Stable_State
    
    results <- bind_rows(results, data.frame(
      iteration        = i,
      scenario         = rest,
      disease_risk     = risk["present"],
      regulation_strong = regulation["strong"],
      years_1_3        = time["1_3_years"],
      years_4_7        = time["4_7_years"],
      years_8_12       = time["8_12_years"],
      years_13_15      = time["13_15_years"]
    ))
  }
  
}  # END MONTE CARLO

# =========================================================
# POST-PROCESSING
# =========================================================

results$long_term_stability <-
  results$years_8_12 + results$years_13_15

scen_levels <- c("none", "low", "moderate", "high")

tick_cols <- c(
  none     = "#aaaaaa",
  low      = "#bae4b3",
  moderate = "#74c476",
  high     = "#238b45"
)

# =========================================================
# SUMMARY
# =========================================================

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
# PLOT 1 – Hump plot with ribbon and annotations
# =========================================================

baseline_val <- summary_df$mean_disease[summary_df$scenario == "none"]
peak_row     <- summary_df[summary_df$scenario == "low", ]       # hump peak at LOW
high_row     <- summary_df[summary_df$scenario == "high", ]

ggplot(
  summary_df,
  aes(x = factor(scenario, levels = scen_levels),
      y = mean_disease, group = 1)
) +
  geom_ribbon(
    aes(ymin = lower_disease, ymax = upper_disease),
    fill  = "#238b45",
    alpha = 0.12
  ) +
  geom_hline(
    yintercept = baseline_val,
    linetype   = "dashed",
    colour     = "grey50",
    linewidth  = 0.6
  ) +
  annotate(
    "text",
    x = 0.55, y = baseline_val + 0.012,
    label = "No-restoration baseline",
    hjust = 0, size = 3.2, colour = "grey45"
  ) +
  geom_line(linewidth = 1.5, colour = "#238b45") +
  geom_point(
    aes(fill = scenario),
    shape = 21, size = 5,
    colour = "white", stroke = 1.8
  ) +
  # peak at low restoration
  annotate(
    "segment",
    x    = peak_row$scenario, xend = peak_row$scenario,
    y    = peak_row$upper_disease + 0.005,
    yend = peak_row$upper_disease + 0.025,
    colour = "#74c476", linewidth = 0.5
  ) +
  annotate(
    "label",
    x = peak_row$scenario,
    y = peak_row$upper_disease + 0.030,
    label     = paste0("Peak\n", round(peak_row$mean_disease, 2)),
    fill      = "#e5f5e0",
    colour    = "#238b45",
    size      = 3.2, label.size = 0, fontface = "bold"
  ) +
  # endpoint below baseline at high restoration
  annotate(
    "segment",
    x    = high_row$scenario, xend = high_row$scenario,
    y    = high_row$lower_disease - 0.005,
    yend = high_row$lower_disease - 0.025,
    colour = "#238b45", linewidth = 0.5
  ) +
  annotate(
    "label",
    x = high_row$scenario,
    y = high_row$lower_disease - 0.030,
    label     = paste0("Below baseline\n", round(high_row$mean_disease, 2)),
    fill      = "#e5f5e0",
    colour    = "#238b45",
    size      = 3.2, label.size = 0, fontface = "bold"
  ) +
  scale_fill_manual(values = tick_cols) +
  scale_y_continuous(
    limits = c(0, 0.80),
    breaks = seq(0, 0.8, 0.1),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position  = "none",
    panel.grid.minor = element_blank()
  ) +
  labs(
    title    = "Tick-borne disease risk across restoration intensities",
    x        = "Restoration intensity",
    y        = "Disease risk probability"
  )

# =========================================================
# PLOT 2 – Ecological regulation recovery
# =========================================================

ggplot(
  summary_df,
  aes(x = factor(scenario, levels = scen_levels),
      y = mean_regulation, group = 1)
) +
  geom_ribbon(
    aes(ymin = lower_regulation, ymax = upper_regulation),
    fill = "#238b45", alpha = 0.12
  ) +
  geom_line(linewidth = 1.4, colour = "#238b45") +
  geom_point(aes(fill = scenario), shape = 21, size = 5,
             colour = "white", stroke = 1.8) +
  scale_fill_manual(values = tick_cols) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none", panel.grid.minor = element_blank()) +
  labs(
    title    = "Ecological regulation recovery — tick-borne system",
    x        = "Restoration intensity",
    y        = "Probability of strong regulation"
  )

# =========================================================
# PLOT 3 – Disease risk vs long-term stability
# =========================================================

ggplot(
  results,
  aes(x = disease_risk, y = long_term_stability,
      colour = scenario)
) +
  geom_point(alpha = 0.18, size = 1.8) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1.4) +
  scale_colour_manual(values = tick_cols) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank()) +
  labs(
    title  = "Tick disease risk vs long-term stability",
    x      = "Disease risk",
    y      = "Long-term stability (8–15 yr)",
    colour = "Restoration intensity"
  )

# =========================================================
# PLOT 4 – Stabilisation trajectories (very delayed signal)
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
) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.18, colour = NA) +
  geom_line(linewidth = 1.4) +
  geom_point(size = 3) +
  scale_colour_manual(values = tick_cols) +
  scale_fill_manual(values   = tick_cols) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank()) +
  labs(
    title    = "Tick-borne disease stabilisation trajectories",
    x        = "Time window",
    y        = "Probability",
    colour   = "Restoration intensity",
    fill     = "Restoration intensity"
  )
