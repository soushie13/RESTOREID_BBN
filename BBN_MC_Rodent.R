# =========================================================
# RODENT-BORNE DISEASE MODEL
# MONTE CARLO BAYESIAN BELIEF NETWORK
#
# Mechanism:
#   Restoration → improved habitat quality
#               → higher biodiversity
#               → reduced rodent dominance
#               → lower reservoir abundance
#               → lower disease risk
#
# Expected pattern:
#   Disease risk   : steady decline, steeper than multi-host
#                    (direct pathway: habitat → rodent → risk)
#   Stabilisation  : slow but monotonic increase
#   Uncertainty    : moderate
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
  
  Habitat_Quality = c(
    "low",
    "moderate",
    "high"
  ),
  
  Rodent_Abundance = c(
    "low",
    "moderate",
    "high"
  ),
  
  Rodent_Disease_Risk = c(
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
# DAG 
# =========================================================

dag <- model2network(
  paste0(
    "[Region_Context]",
    "[Restoration_Intensity]",
    "[Habitat_Quality|Region_Context:Restoration_Intensity]",
    "[Rodent_Abundance|Habitat_Quality]",
    "[Rodent_Disease_Risk|Rodent_Abundance]",
    "[Ecological_Regulation|Restoration_Intensity:Rodent_Disease_Risk]",
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
  # HABITAT QUALITY
  # dim = c(3,4,3) = [Habitat, Rest, Region]
  # Loop: Region outer (dim3), Rest inner (dim2)  ← was correct
  #
  # Mechanism: restoration directly improves habitat quality
  # through vegetation recovery, reduced soil disturbance,
  # and return of native plant communities.
  # High-risk region adds a penalty to low habitat quality.
  #
  # Pattern: steep monotonic improvement with restoration
  # (steeper than multi-host system because the pathway
  # habitat → rodent → risk is more direct)
  # -------------------------------------------------------
  
  habitat_vals <- c()
  
  for(region in 1:3){         # dim 3 – Region (outer)
    for(rest in 1:4){         # dim 2 – Restoration (inner)
      
      probs <- switch(rest,
                      
                      # none – degraded habitat, low quality dominant
                      `1` = c(sample_prob(0.62, 0.80),   # low    ← dominant
                              sample_prob(0.14, 0.28),   # moderate
                              sample_prob(0.02, 0.10)),  # high
                      
                      # low – marginal recovery, still mostly poor
                      `2` = c(sample_prob(0.38, 0.58),   # low
                              sample_prob(0.26, 0.40),   # moderate
                              sample_prob(0.08, 0.20)),  # high
                      
                      # moderate – quality clearly shifting
                      `3` = c(sample_prob(0.18, 0.35),   # low
                              sample_prob(0.32, 0.46),   # moderate
                              sample_prob(0.25, 0.42)),  # high
                      
                      # high – mostly high quality; mature diverse habitat
                      `4` = c(sample_prob(0.04, 0.14),   # low
                              sample_prob(0.18, 0.34),   # moderate
                              sample_prob(0.55, 0.75))   # high   ← dominant
      )
      
      # High-risk region: additional low-quality penalty
      if(region == 3) probs[1] <- probs[1] + 0.10
      
      probs        <- probs / sum(probs)
      habitat_vals <- c(habitat_vals, probs)
    }
  }
  
  cpt_habitat <- array(
    habitat_vals,
    dim      = c(3,4,3),
    dimnames = list(
      Habitat_Quality       = states$Habitat_Quality,
      Restoration_Intensity = states$Restoration_Intensity,
      Region_Context        = states$Region_Context
    )
  )
  
  # -------------------------------------------------------
  # RODENT ABUNDANCE
  # dim = c(3,3) = [Abundance, Habitat]
  # Single parent — flat vector fills correctly (Habitat slowest).
  #
  # Mechanism: high habitat quality supports diverse competitor
  # communities (raptors, mustelids, diverse small mammals)
  # that suppress rodent dominance. Low habitat quality
  # favours r-selected rodent irruptions.
  #
  # The rodent suppression gradient is steep:
  # low habitat  → high rodent abundance (0.72–0.82 present)
  # high habitat → low rodent abundance  (0.68–0.80 absent)
  # This steep gradient is what produces the steep risk decline.
  # -------------------------------------------------------
  
  rodent_vals <- c(
    # low habitat quality    → high rodent dominance
    0.05, 0.18, 0.77,
    # moderate habitat       → mixed, moderate uncertainty
    0.28, 0.48, 0.24,
    # high habitat quality   → rodent suppression
    0.72, 0.22, 0.06
  )
  
  rodent_vals <- rodent_vals +
    runif(length(rodent_vals), -0.03, 0.03)
  rodent_vals[rodent_vals < 0.01] <- 0.01
  
  rodent_vals <- unlist(lapply(
    split(rodent_vals, ceiling(seq_along(rodent_vals)/3)),
    function(x) x / sum(x)
  ))
  
  # FIX 1: array() was missing — cpt_rodents was never created
  cpt_rodents <- array(
    rodent_vals,
    dim      = c(3,3),
    dimnames = list(
      Rodent_Abundance = states$Rodent_Abundance,
      Habitat_Quality  = states$Habitat_Quality
    )
  )
  
  # -------------------------------------------------------
  # DISEASE RISK
  # dim = c(2,3) = [Disease, Abundance]
  # Single parent — flat vector fills correctly.
  #
  # Tight, steep relationship: rodent abundance is the
  # proximate driver of transmission risk (direct contact,
  # contaminated excreta, infected ectoparasites).
  # Small noise to reflect spillover stochasticity.
  # -------------------------------------------------------
  
  disease_vals <- c(
    # low abundance    → mostly absent
    0.92, 0.08,
    # moderate         → uncertain, near 50/50
    0.55, 0.45,
    # high abundance   → mostly present
    0.10, 0.90
  )
  
  disease_vals <- disease_vals +
    runif(length(disease_vals), -0.03, 0.03)
  disease_vals[disease_vals < 0.01] <- 0.01
  
  disease_vals <- unlist(lapply(
    split(disease_vals, ceiling(seq_along(disease_vals)/2)),
    function(x) x / sum(x)
  ))
  
  cpt_disease <- array(
    disease_vals,
    dim      = c(2,3),
    dimnames = list(
      Rodent_Disease_Risk = states$Rodent_Disease_Risk,
      Rodent_Abundance    = states$Rodent_Abundance
    )
  )
  
  # -------------------------------------------------------
  # ECOLOGICAL REGULATION
  # dim = c(3,4,2) = [Regulation, Rest, Disease]
  # Loop: Disease outer (dim3, slowest), Rest inner (dim2)
  #
  # FIX 2: loop order was inverted (Rest outer, Disease inner)
  # → regulation probabilities were assigned to wrong cells
  #
  # Restoration drives regulation through two pathways:
  #   (a) direct — biodiversity recovery, predator return
  #   (b) indirect — rodent suppression removes positive
  #       feedback between high abundance and weak regulation
  #
  # Disease present adds a small penalty to regulation
  # (ongoing outbreak pressure delays recovery).
  # -------------------------------------------------------
  
  reg_vals <- c()
  
  for(disease in 1:2){        # dim 3 – Disease (outer)
    for(rest in 1:4){         # dim 2 – Restoration (inner)
      
      probs <- switch(rest,
                      
                      # none – weak regulation dominant
                      `1` = c(sample_prob(0.65, 0.82),   # weak   ← dominant
                              sample_prob(0.12, 0.24),   # moderate
                              sample_prob(0.02, 0.10)),  # strong
                      
                      # low – weak still dominant, small shift
                      `2` = c(sample_prob(0.45, 0.62),   # weak
                              sample_prob(0.24, 0.38),   # moderate
                              sample_prob(0.08, 0.22)),  # strong
                      
                      # moderate – regulation building monotonically
                      `3` = c(sample_prob(0.22, 0.38),   # weak
                              sample_prob(0.30, 0.44),   # moderate
                              sample_prob(0.25, 0.42)),  # strong
                      
                      # high – strong regulation dominant
                      `4` = c(sample_prob(0.04, 0.12),   # weak
                              sample_prob(0.12, 0.26),   # moderate
                              sample_prob(0.65, 0.82))   # strong ← dominant
      )
      
      # Disease present: small penalty on regulation
      # (outbreak pressure delays biodiversity recovery)
      if(disease == 2) probs[1] <- probs[1] + 0.05
      
      probs    <- probs / sum(probs)
      reg_vals <- c(reg_vals, probs)
    }
  }
  
  # FIX 2 (continued): array() was missing — cpt_regulation never created
  cpt_regulation <- array(
    reg_vals,
    dim      = c(3,4,2),
    dimnames = list(
      Ecological_Regulation = states$Ecological_Regulation,
      Restoration_Intensity = states$Restoration_Intensity,
      Rodent_Disease_Risk   = states$Rodent_Disease_Risk
    )
  )
  
  # -------------------------------------------------------
  # TIME TO STABLE STATE
  # dim = c(4,3) = [Time, Regulation]
  # Single parent — Regulation is dim 2 (slowest). ✓
  #
  # FIX 3: array() was missing — cpt_time was never created
  #
  # Pattern: slow but monotonic stabilisation increase
  #   Weak regulation  → fast degraded plateau (1-3 yr mass)
  #   Moderate         → gradual shift; peak at 4-7 yr
  #   Strong           → slow, delayed; peak at 8-12 yr
  #                      (slower than tick system but less
  #                       extreme — no trophic cascade lag)
  # The key signature is MONOTONIC — each step clearly shifts
  # the distribution later, with moderate uncertainty.
  # -------------------------------------------------------
  
  time_vals <- c()
  
  for(reg in 1:3){
    
    probs <- switch(reg,
                    
                    # weak – rapid degraded equilibrium
                    `1` = c(sample_prob(0.48, 0.62),   # 1-3 yr  ← dominant
                            sample_prob(0.24, 0.34),   # 4-7 yr
                            sample_prob(0.06, 0.14),   # 8-12 yr
                            sample_prob(0.02, 0.08)),  # 13-15 yr
                    
                    # moderate – peak shifting to middle windows
                    `2` = c(sample_prob(0.14, 0.24),   # 1-3 yr
                            sample_prob(0.32, 0.44),   # 4-7 yr  ← dominant
                            sample_prob(0.22, 0.34),   # 8-12 yr
                            sample_prob(0.10, 0.22)),  # 13-15 yr
                    
                    # strong – delayed, concentrated in 8-12 yr
                    # (slower than no-restoration but not as extreme as tick)
                    `3` = c(sample_prob(0.03, 0.09),   # 1-3 yr
                            sample_prob(0.10, 0.20),   # 4-7 yr
                            sample_prob(0.42, 0.56),   # 8-12 yr  ← dominant
                            sample_prob(0.24, 0.38))   # 13-15 yr
    )
    
    probs     <- probs / sum(probs)
    time_vals <- c(time_vals, probs)
  }
  
  # FIX 3 (continued): array() was missing
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
      Habitat_Quality       = cpt_habitat,
      Rodent_Abundance      = cpt_rodents,
      Rodent_Disease_Risk   = cpt_disease,
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
    
    risk       <- querygrain(bn_temp, nodes = "Rodent_Disease_Risk")$Rodent_Disease_Risk
    regulation <- querygrain(bn_temp, nodes = "Ecological_Regulation")$Ecological_Regulation
    time       <- querygrain(bn_temp, nodes = "Time_to_Stable_State")$Time_to_Stable_State
    
    results <- bind_rows(results, data.frame(
      iteration             = i,
      scenario              = rest,
      disease_risk          = risk["present"],
      ecological_regulation = regulation["strong"],
      years_1_3             = time["1_3_years"],
      years_4_7             = time["4_7_years"],
      years_8_12            = time["8_12_years"],
      years_13_15           = time["13_15_years"]
    ))
  }
  
}  # END MONTE CARLO

# =========================================================
# POST-PROCESSING
# =========================================================

results$long_term_stability <-
  results$years_8_12 + results$years_13_15

scen_levels <- c("none", "low", "moderate", "high")

rodent_cols <- c(
  none     = "#aaaaaa",
  low      = "#fcae91",
  moderate = "#fb6a4a",
  high     = "#cb181d"
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
    mean_regulation   = mean(ecological_regulation),
    lower_regulation  = quantile(ecological_regulation, 0.025),
    upper_regulation  = quantile(ecological_regulation, 0.975),
    mean_stability    = mean(long_term_stability),
    lower_stability   = quantile(long_term_stability, 0.025),
    upper_stability   = quantile(long_term_stability, 0.975),
    .groups = "drop"
  )

print(summary_df)

# =========================================================
# PLOT 1 – Steady decline with ribbon and annotations
# =========================================================

baseline_val <- summary_df$mean_disease[summary_df$scenario == "none"]
high_row     <- summary_df[summary_df$scenario == "high", ]

ggplot(
  summary_df,
  aes(x = factor(scenario, levels = scen_levels),
      y = mean_disease, group = 1)
) +
  geom_ribbon(
    aes(ymin = lower_disease, ymax = upper_disease),
    fill  = "#cb181d",
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
    x = 0.55, y = baseline_val + 0.013,
    label = "No-restoration baseline",
    hjust = 0, size = 3.2, colour = "grey45"
  ) +
  geom_line(linewidth = 1.5, colour = "#cb181d") +
  geom_point(
    aes(fill = scenario),
    shape = 21, size = 5,
    colour = "white", stroke = 1.8
  ) +
  annotate(
    "segment",
    x    = high_row$scenario, xend = high_row$scenario,
    y    = high_row$lower_disease - 0.005,
    yend = high_row$lower_disease - 0.025,
    colour = "#cb181d", linewidth = 0.5
  ) +
  annotate(
    "label",
    x = high_row$scenario,
    y = high_row$lower_disease - 0.030,
    label     = paste0("Steady decline endpoint\n",
                       round(high_row$mean_disease, 2)),
    fill      = "#fff5f0",
    colour    = "#cb181d",
    size      = 3.2, label.size = 0, fontface = "bold"
  ) +
  scale_fill_manual(values = rodent_cols) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position  = "none",
    panel.grid.minor = element_blank()
  ) +
  labs(
    title    = "Rodent-borne disease risk across restoration intensities",
    x        = "Restoration intensity",
    y        = "Disease risk probability"
  )

# =========================================================
# PLOT 2 – Ecological regulation recovery 
# =========================================================

reg_summary <- results %>%
  group_by(scenario) %>%
  summarise(
    mean_reg = mean(ecological_regulation, na.rm = TRUE),
    lower    = quantile(ecological_regulation, 0.025, na.rm = TRUE),
    upper    = quantile(ecological_regulation, 0.975, na.rm = TRUE),
    .groups  = "drop"
  )

ggplot(
  reg_summary,
  aes(x = factor(scenario, levels = scen_levels),
      y = mean_reg, group = 1)
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    fill = "#cb181d", alpha = 0.12
  ) +
  geom_line(linewidth = 1.4, colour = "#cb181d") +
  geom_point(
    aes(fill = scenario),
    shape = 21, size = 5,
    colour = "white", stroke = 1.8
  ) +
  scale_fill_manual(values = rodent_cols) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.1),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position  = "none",
    panel.grid.minor = element_blank()
  ) +
  labs(
    title    = "Ecological regulation recovery — rodent-borne system",
    x        = "Restoration intensity",
    y        = "Probability of strong ecological regulation"
  )

# =========================================================
# PLOT 3 – Disease risk vs long-term stability
# =========================================================

ggplot(
  results,
  aes(disease_risk, long_term_stability, colour = scenario)
) +
  geom_point(alpha = 0.15, size = 1.8) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1.4) +
  scale_colour_manual(values = rodent_cols) +
  scale_fill_manual(values   = rodent_cols) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank()) +
  labs(
    title  = "Rodent-borne diseases: disease risk vs long-term stability",
    x      = "Disease risk",
    y      = "Long-term stability (8–15 yr)",
    colour = "Restoration intensity"
  )

# =========================================================
# PLOT 4 – Stabilisation trajectories
# =========================================================

traj <- results %>%
  pivot_longer(
    starts_with("years_"),
    names_to  = "time_period",
    values_to = "probability"
  )

corridor <- traj %>%
  group_by(scenario, time_period) %>%
  summarise(
    mean_prob = mean(probability),
    lower     = quantile(probability, 0.025),
    upper     = quantile(probability, 0.975),
    .groups   = "drop"
  )

corridor$time_period <- factor(
  corridor$time_period,
  levels = c("years_1_3","years_4_7","years_8_12","years_13_15"),
  labels = c("1–3 yr","4–7 yr","8–12 yr","13–15 yr")
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
  scale_colour_manual(values = rodent_cols) +
  scale_fill_manual(values   = rodent_cols) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank()) +
  labs(
    title    = "Rodent-borne disease stabilisation trajectories",
    x        = "Time window",
    y        = "Probability",
    colour   = "Restoration intensity",
    fill     = "Restoration intensity"
  )