# =========================================================
# DIRECT-CONTACT ZOONOTIC DISEASE MODEL
# MULTI-HOST SYSTEMS — NIPAH / HENIPA / FILO-LIKE
# MONTE CARLO BAYESIAN BELIEF NETWORK
#
# Mechanism:
#   Restoration → reduced fragmentation
#               → reduced edge effects
#               → reduced wildlife-human interfaces
#   Host diversity effects: weaker and less predictable
#     (dilution effect present but not dominant)
#
# Expected pattern:
#   Risk          : gradual decline with restoration
#   Stabilisation : weakly predictable
#   Uncertainty   : very broad
#
# Mature ecosystems after long-term restoration
# decrease human-domestic-host/reservoir interfaces
# → disease risk declines, but slowly and with wide uncertainty



library(bnlearn)
library(gRain)
library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(123)

# =========================================================
# STATES
# =========================================================

states <- list(
  
  Region_Context = c(
    "Low_Risk_Context",
    "Moderate_Risk_Context",
    "High_Risk_Context"
  ),
  
  Restoration_Intensity = c(
    "none",
    "low",
    "moderate",
    "high"
  ),
  
  # 2 states only — must match dimnames exactly throughout
  Resilient_Ecological_Regulation = c(
    "low",
    "high"
  ),
  
  Habitat_Disturbance = c(
    "low",
    "moderate",
    "high"
  ),
  
  Host_Diversity = c(
    "low",
    "moderate",
    "high"
  ),
  
  Human_Wildlife_Contact = c(
    "low",
    "moderate",
    "high"
  ),
  
  MultiHost_Amplification = c(
    "low",
    "high"
  ),
  
  Direct_Transmission_Risk = c(
    "absent",
    "present"
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
    "[Habitat_Disturbance|Region_Context:Restoration_Intensity]",
    "[Host_Diversity|Habitat_Disturbance]",
    "[Human_Wildlife_Contact|Habitat_Disturbance:Restoration_Intensity]",
    "[MultiHost_Amplification|Host_Diversity:Human_Wildlife_Contact]",
    "[Direct_Transmission_Risk|MultiHost_Amplification]",
    "[Resilient_Ecological_Regulation|Restoration_Intensity:Host_Diversity]",
    "[Time_to_Stable_State|",
    "Restoration_Intensity:",
    "Direct_Transmission_Risk:",
    "Resilient_Ecological_Regulation",
    "]"
  )
)

# =========================================================
# HELPER
# =========================================================

sample_prob <- function(mn, mx) runif(1, mn, mx)

# =========================================================
# MONTE CARLO LOOP
# =========================================================

n_iter  <- 1000
results <- data.frame()

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
  # HABITAT DISTURBANCE
  # dim = c(3,4,3) = [Disturbance, Rest, Region]
  # Loop: Region outer (dim3, slowest), Rest inner (dim2)
  #
  # Mechanism: restoration reduces fragmentation and edge
  # effects → disturbance decreases monotonically.
  # High-risk region adds a penalty to high disturbance.
  # -------------------------------------------------------
  
  disturbance_vals <- c()
  
  for(region in 1:3){           # dim 3 – Region (outer)
    for(rest in 1:4){           # dim 2 – Restoration (inner)
      
      probs <- switch(rest,
                      
                      # none – ongoing fragmentation, high disturbance
                      `1` = c(sample_prob(0.05, 0.15),   # low dist
                              sample_prob(0.20, 0.30),   # mod dist
                              sample_prob(0.55, 0.75)),  # high dist ← dominant
                      
                      # low – some edge reduction but fragmentation persists
                      `2` = c(sample_prob(0.10, 0.25),
                              sample_prob(0.35, 0.50),
                              sample_prob(0.30, 0.45)),
                      
                      # moderate – fragmentation meaningfully reduced
                      `3` = c(sample_prob(0.20, 0.40),
                              sample_prob(0.35, 0.45),
                              sample_prob(0.15, 0.30)),
                      
                      # high – landscape connectivity restored, low edge effects
                      `4` = c(sample_prob(0.50, 0.70),   # low dist ← dominant
                              sample_prob(0.20, 0.35),
                              sample_prob(0.05, 0.15))
      )
      
      # High-risk region: elevated background disturbance
      if(region == 3) probs[3] <- probs[3] + 0.10
      
      probs            <- probs / sum(probs)
      disturbance_vals <- c(disturbance_vals, probs)
    }
  }
  
  cpt_disturbance <- array(
    disturbance_vals,
    dim      = c(3,4,3),
    dimnames = list(
      Habitat_Disturbance   = states$Habitat_Disturbance,
      Restoration_Intensity = states$Restoration_Intensity,
      Region_Context        = states$Region_Context
    )
  )
  
  # -------------------------------------------------------
  # HOST DIVERSITY
  # dim = c(3,3) = [Diversity, Disturbance]
  # Single parent — loop not needed; flat vector fills correctly.
  #
  # Low disturbance → high diversity (intact community)
  # High disturbance → low diversity (disrupted community)
  # Small noise added to reflect ecological unpredictability.
  # -------------------------------------------------------
  
  diversity_vals <- c(
    # low disturbance  → mostly high diversity
    0.10, 0.25, 0.65,
    # mod disturbance  → mixed
    0.30, 0.45, 0.25,
    # high disturbance → mostly low diversity
    0.70, 0.20, 0.10
  )
  
  diversity_vals <- diversity_vals +
    runif(length(diversity_vals), -0.05, 0.05)
  diversity_vals[diversity_vals < 0.01] <- 0.01
  
  diversity_vals <- unlist(lapply(
    split(diversity_vals, ceiling(seq_along(diversity_vals)/3)),
    function(x) x / sum(x)
  ))
  
  cpt_diversity <- array(
    diversity_vals,
    dim      = c(3,3),
    dimnames = list(
      Host_Diversity      = states$Host_Diversity,
      Habitat_Disturbance = states$Habitat_Disturbance
    )
  )
  
  # -------------------------------------------------------
  # HUMAN-WILDLIFE CONTACT
  # dim = c(3,3,4) = [Contact, Disturbance, Rest]
  # Loop: Rest outer (dim3, slowest), Disturbance inner (dim2)
  #
  # Core mechanism for this disease system:
  # Restoration reduces fragmentation → fewer edge habitats
  # → fewer human-reservoir interfaces → lower contact rates.
  # This is the PRIMARY driver of risk decline.
  # High disturbance adds a contact penalty.
  # -------------------------------------------------------
  
  contact_vals <- c()
  
  for(rest in 1:4){             # dim 3 – Restoration (outer)
    for(disturbance in 1:3){   # dim 2 – Disturbance (inner)
      
      probs <- switch(rest,
                      
                      # none – fragmented landscape, high interface zones
                      `1` = c(sample_prob(0.10, 0.20),   # low contact
                              sample_prob(0.25, 0.35),   # mod contact
                              sample_prob(0.50, 0.65)),  # high contact ← dominant
                      
                      # low – some interface reduction, still transitional
                      `2` = c(sample_prob(0.15, 0.30),
                              sample_prob(0.35, 0.45),
                              sample_prob(0.30, 0.45)),
                      
                      # moderate – edge effects reduced, interfaces declining
                      `3` = c(sample_prob(0.30, 0.50),
                              sample_prob(0.30, 0.42),
                              sample_prob(0.15, 0.28)),
                      
                      # high – landscape connectivity restored, interfaces low
                      # Mature ecosystem → wildlife retreats from human zones
                      `4` = c(sample_prob(0.55, 0.75),   # low contact ← dominant
                              sample_prob(0.18, 0.30),
                              sample_prob(0.05, 0.15))
      )
      
      # High disturbance forces more contact regardless of restoration
      if(disturbance == 3) probs[3] <- probs[3] + 0.10
      
      probs        <- probs / sum(probs)
      contact_vals <- c(contact_vals, probs)
    }
  }
  
  cpt_contact <- array(
    contact_vals,
    dim      = c(3,3,4),
    dimnames = list(
      Human_Wildlife_Contact = states$Human_Wildlife_Contact,
      Habitat_Disturbance    = states$Habitat_Disturbance,
      Restoration_Intensity  = states$Restoration_Intensity
    )
  )
  
  # -------------------------------------------------------
  # MULTI-HOST AMPLIFICATION
  # dim = c(2,3,3) = [Amplification, Diversity, Contact]
  # Loop: Contact outer (dim3, slowest), Diversity inner (dim2)
  #
  # HOST DIVERSITY EFFECT — WEAK AND UNPREDICTABLE:
  # Unlike vector-borne systems, dilution is less reliable here.
  # High diversity can EITHER suppress transmission (dilution)
  # OR amplify it (more competent hosts in community).
  # → Most cells have wide, overlapping probability ranges.
  # → Only the extreme corner cases (low diversity + high
  #   contact, high diversity + low contact) have clearer
  #   signal; all mid-range combinations are near-random.
  # -------------------------------------------------------
  
  amplification_vals <- c()
  
  for(contact in 1:3){          # dim 3 – Contact (outer)
    for(diversity in 1:3){      # dim 2 – Diversity (inner)
      
      probs <- if(contact == 3 & diversity == 1){
        # High contact + low diversity → near-certain amplification
        # Few hosts but frequent human exposure
        c(sample_prob(0.05, 0.15),
          sample_prob(0.85, 0.95))
        
      } else if(contact == 1 & diversity == 3){
        # Low contact + high diversity → dilution suppresses risk
        # Rare exposure, transmission diluted across many species
        c(sample_prob(0.70, 0.88),
          sample_prob(0.12, 0.30))
        
      } else if(contact == 3 & diversity == 3){
        # High contact + high diversity → AMBIGUOUS
        # More competent hosts but also more dilution hosts
        # → unpredictable; wide range centred near 50/50
        c(sample_prob(0.30, 0.55),
          sample_prob(0.45, 0.70))
        
      } else if(contact == 2 & diversity == 2){
        # Moderate everything → maximum uncertainty
        c(sample_prob(0.35, 0.65),
          sample_prob(0.35, 0.65))
        
      } else {
        # All other combinations: weak signal, broad uncertainty
        c(sample_prob(0.30, 0.60),
          sample_prob(0.40, 0.70))
      }
      
      probs              <- probs / sum(probs)
      amplification_vals <- c(amplification_vals, probs)
    }
  }
  
  cpt_amplification <- array(
    amplification_vals,
    dim      = c(2,3,3),
    dimnames = list(
      MultiHost_Amplification = states$MultiHost_Amplification,
      Host_Diversity          = states$Host_Diversity,
      Human_Wildlife_Contact  = states$Human_Wildlife_Contact
    )
  )
  
  # -------------------------------------------------------
  # DIRECT TRANSMISSION RISK
  # dim = c(2,2) = [Disease, Amplification]
  # Single parent — flat vector, no loop needed.
  #
  # Wide noise on both states to reflect spillover
  # unpredictability in filovirus/henipavirus systems.
  # -------------------------------------------------------
  
  disease_vals <- c(
    # low amplification  → mostly absent
    0.80, 0.20,
    # high amplification → mostly present
    0.15, 0.85
  )
  
  disease_vals <- disease_vals +
    runif(length(disease_vals), -0.06, 0.06)
  disease_vals[disease_vals < 0.01] <- 0.01
  
  disease_vals <- unlist(lapply(
    split(disease_vals, ceiling(seq_along(disease_vals)/2)),
    function(x) x / sum(x)
  ))
  
  cpt_disease <- array(
    disease_vals,
    dim      = c(2,2),
    dimnames = list(
      Direct_Transmission_Risk = states$Direct_Transmission_Risk,
      MultiHost_Amplification  = states$MultiHost_Amplification
    )
  )
  
  # -------------------------------------------------------
  # RESILIENT ECOLOGICAL REGULATION
  # dim = c(2,4,3) = [Regulation, Rest, Diversity]
  # Loop: Diversity outer (dim3, slowest), Rest inner (dim2)
  #
  # 2 states only: "low" / "high"
  # High restoration + high diversity → strong regulation
  # but diversity effect is weak and noisy (broad ranges).
  # Under no restoration: weak regulation regardless of diversity.
  # -------------------------------------------------------
  
  regulation_vals <- c()
  
  for(diversity in 1:3){        # dim 3 – Diversity (outer)
    for(rest in 1:4){           # dim 2 – Restoration (inner)
      
      probs <- switch(rest,
                      
                      # none – poor structural complexity, weak regulation
                      `1` = c(sample_prob(0.68, 0.85),   # low reg ← dominant
                              sample_prob(0.15, 0.32)),  # high reg
                      
                      # low – marginal structural recovery
                      `2` = c(sample_prob(0.55, 0.72),
                              sample_prob(0.28, 0.45)),
                      
                      # moderate – regulation building but highly variable
                      `3` = c(sample_prob(0.38, 0.58),
                              sample_prob(0.42, 0.62)),
                      
                      # high – mature ecosystem, strong regulation likely
                      # but still uncertain — some systems never fully recover
                      `4` = c(sample_prob(0.15, 0.35),
                              sample_prob(0.65, 0.85))
      )
      
      # Biodiversity bonus: weak and noisy signal
      # High diversity increases regulation probability slightly,
      # but the effect is deliberately small and overlapping
      if(diversity == 3) probs[2] <- probs[2] + sample_prob(0.03, 0.12)
      if(diversity == 1) probs[1] <- probs[1] + sample_prob(0.02, 0.08)
      
      probs           <- probs / sum(probs)
      regulation_vals <- c(regulation_vals, probs)
    }
  }
  
  cpt_regulation <- array(
    regulation_vals,
    dim      = c(2,4,3),
    dimnames = list(
      Resilient_Ecological_Regulation = states$Resilient_Ecological_Regulation,
      Restoration_Intensity           = states$Restoration_Intensity,
      Host_Diversity                  = states$Host_Diversity
    )
  )
  
  # -------------------------------------------------------
  # TIME TO STABLE STATE
  # dim = c(4,2,4,2) = [Time, Regulation, Rest, Disease]
  # Loop: Disease outer (dim4), Rest next (dim3),
  #       Regulation inner (dim2)
  #
  # Stabilisation is weakly predictable — wide ranges throughout.
  # No restoration: rapid but degraded plateau (early time windows)
  # High restoration: delayed, broad distribution across all windows
  # -------------------------------------------------------
  
  time_vals <- c()
  
  for(disease in 1:2){          # dim 4 – Disease (outer)
    for(rest in 1:4){           # dim 3 – Restoration
      for(regulation in 1:2){   # dim 2 – Regulation (inner)
        
        probs <- if(regulation == 2 & rest == 4){
          # High regulation + high restoration
          # → delayed but genuinely resilient stabilisation
          c(sample_prob(0.03, 0.10),   # 1-3 yr
            sample_prob(0.10, 0.20),   # 4-7 yr
            sample_prob(0.25, 0.38),   # 8-12 yr
            sample_prob(0.38, 0.55))   # 13-15 yr ← dominant
          
        } else if(rest == 1){
          # No restoration: fast plateau at degraded state
          c(sample_prob(0.38, 0.52),   # 1-3 yr ← dominant
            sample_prob(0.25, 0.35),   # 4-7 yr
            sample_prob(0.10, 0.20),   # 8-12 yr
            sample_prob(0.04, 0.12))   # 13-15 yr
          
        } else {
          # Transitional: broad, flat, weakly predictable
          # This is the key signature — very wide uncertainty
          c(sample_prob(0.12, 0.32),
            sample_prob(0.18, 0.35),
            sample_prob(0.20, 0.35),
            sample_prob(0.18, 0.35))
        }
        
        probs     <- probs / sum(probs)
        time_vals <- c(time_vals, probs)
      }
    }
  }
  
  cpt_time <- array(
    time_vals,
    dim      = c(4,2,4,2),
    dimnames = list(
      Time_to_Stable_State             = states$Time_to_Stable_State,
      Resilient_Ecological_Regulation  = states$Resilient_Ecological_Regulation,
      Restoration_Intensity            = states$Restoration_Intensity,
      Direct_Transmission_Risk         = states$Direct_Transmission_Risk
    )
  )
  
  # -------------------------------------------------------
  # FIT NETWORK
  # NOTE: Resilient_Ecological_Regulation listed only ONCE
  # -------------------------------------------------------
  
  fitted_bn <- custom.fit(
    dag,
    dist = list(
      Region_Context                   = cpt_region,
      Restoration_Intensity            = cpt_restoration,
      Habitat_Disturbance              = cpt_disturbance,
      Host_Diversity                   = cpt_diversity,
      Human_Wildlife_Contact           = cpt_contact,
      MultiHost_Amplification          = cpt_amplification,
      Direct_Transmission_Risk         = cpt_disease,
      Resilient_Ecological_Regulation  = cpt_regulation,
      Time_to_Stable_State             = cpt_time
    )
  )
  
  bn_grain <- as.grain(fitted_bn)
  
  # -------------------------------------------------------
  # SCENARIO QUERIES  (High-risk region)
  # -------------------------------------------------------
  
  for(restoration_state in states$Restoration_Intensity){
    
    bn_temp <- setEvidence(
      bn_grain,
      nodes  = c("Region_Context", "Restoration_Intensity"),
      states = c("High_Risk_Context", restoration_state)
    )
    
    disease_q    <- querygrain(bn_temp, nodes = "Direct_Transmission_Risk")$Direct_Transmission_Risk
    time_q       <- querygrain(bn_temp, nodes = "Time_to_Stable_State")$Time_to_Stable_State
    regulation_q <- querygrain(bn_temp, nodes = "Resilient_Ecological_Regulation")$Resilient_Ecological_Regulation
    
    results <- bind_rows(results, data.frame(
      iteration            = i,
      scenario             = restoration_state,
      disease_risk         = disease_q["present"],
      resilient_regulation = regulation_q["high"],
      years_1_3            = time_q["1_3_years"],
      years_4_7            = time_q["4_7_years"],
      years_8_12           = time_q["8_12_years"],
      years_13_15          = time_q["13_15_years"]
    ))
  }
  
}  # END MONTE CARLO

# =========================================================
# POST-PROCESSING
# =========================================================

results$long_term_stability <-
  results$years_8_12 + results$years_13_15

scen_levels <- c("none", "low", "moderate", "high")

zoonotic_cols <- c(
  none     = "#aaaaaa",
  low      = "#cbc9e2",
  moderate = "#9e9ac8",
  high     = "#6a51a3"
)

# =========================================================
# SUMMARY TABLE
# =========================================================

summary_df <- results %>%
  group_by(scenario) %>%
  summarise(
    mean_disease      = mean(disease_risk),
    lower_disease     = quantile(disease_risk, 0.025),
    upper_disease     = quantile(disease_risk, 0.975),
    mean_regulation   = mean(resilient_regulation),
    lower_regulation  = quantile(resilient_regulation, 0.025),
    upper_regulation  = quantile(resilient_regulation, 0.975),
    mean_stability    = mean(long_term_stability),
    lower_stability   = quantile(long_term_stability, 0.025),
    upper_stability   = quantile(long_term_stability, 0.975),
    .groups = "drop"
  )

print(summary_df)

# =========================================================
# PLOT 1 – Gradual risk decline with very wide uncertainty
# =========================================================

baseline_val <- summary_df$mean_disease[summary_df$scenario == "none"]

ggplot(
  summary_df,
  aes(x = factor(scenario, levels = scen_levels),
      y = mean_disease, group = 1)
) +
  geom_ribbon(
    aes(ymin = lower_disease, ymax = upper_disease),
    fill  = "#9e9ac8",
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
    label  = "No-restoration baseline",
    hjust  = 0, size = 3.2, colour = "grey45"
  ) +
  geom_line(linewidth = 1.5, colour = "#6a51a3") +
  geom_point(
    aes(fill = scenario),
    shape = 21, size = 5,
    colour = "white", stroke = 1.8
  ) +
  scale_fill_manual(values = zoonotic_cols) +
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
    title    = "Direct-contact zoonotic disease risk across restoration intensities",
    x        = "Restoration intensity",
    y        = "Disease risk probability"
  )

# =========================================================
# PLOT 2 – Disease risk distributions (density)
# Wide, overlapping distributions = broad uncertainty signal
# =========================================================

ggplot(
  results,
  aes(x = disease_risk, fill = scenario, colour = scenario)
) +
  geom_density(alpha = 0.25, linewidth = 0.8) +
  scale_fill_manual(values   = zoonotic_cols) +
  scale_colour_manual(values = zoonotic_cols) +
  scale_x_continuous(limits = c(0, 1)) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank()) +
  labs(
    title    = "Disease risk distributions — direct-contact zoonoses",
    x        = "Disease risk probability",
    y        = "Density",
    fill     = "Restoration intensity",
    colour   = "Restoration intensity"
  )

# =========================================================
# PLOT 3 – Trade-off: disease risk vs long-term stability
# =========================================================

ggplot(
  results,
  aes(x = disease_risk, y = long_term_stability,
      colour = scenario, fill = scenario)
) +
  geom_point(alpha = 0.15, size = 1.8) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1.4) +
  scale_colour_manual(values = zoonotic_cols) +
  scale_fill_manual(values   = zoonotic_cols) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank()) +
  labs(
    title    = "Direct-contact zoonoses: disease risk vs long-term stabilisation",
    x        = "Disease risk probability",
    y        = "Long-term ecological stabilisation (8–15 yr)",
    colour   = "Restoration intensity",
    fill     = "Restoration intensity"
  )

# =========================================================
# PLOT 4 – Stabilisation trajectories with uncertainty corridors
# =========================================================

traj_long <- results %>%
  pivot_longer(
    cols      = starts_with("years_"),
    names_to  = "time_period",
    values_to = "probability"
  )

traj_long$time_period <- factor(
  traj_long$time_period,
  levels = c("years_1_3","years_4_7","years_8_12","years_13_15"),
  labels = c("1–3 yr","4–7 yr","8–12 yr","13–15 yr")
)

corridor_df <- traj_long %>%
  group_by(scenario, time_period) %>%
  summarise(
    mean_prob = mean(probability, na.rm = TRUE),
    lower     = quantile(probability, 0.025, na.rm = TRUE),
    upper     = quantile(probability, 0.975, na.rm = TRUE),
    .groups   = "drop"
  )

ggplot(
  corridor_df,
  aes(x = time_period, y = mean_prob,
      group = scenario, colour = scenario, fill = scenario)
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.18, colour = NA
  ) +
  geom_line(linewidth = 1.4) +
  geom_point(size = 3) +
  scale_colour_manual(values = zoonotic_cols) +
  scale_fill_manual(values   = zoonotic_cols) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank()) +
  labs(
    title    = "Stabilisation trajectories — direct-contact zoonoses",
    x        = "Time window",
    y        = "Probability",
    colour   = "Restoration intensity",
    fill     = "Restoration intensity"
  )
