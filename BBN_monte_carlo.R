# =========================================================
# BAYESIAN BELIEF NETWORK WITH MONTE CARLO UNCERTAINTY
# =========================================================

library(bnlearn)
library(gRain)
library(dplyr)
library(tidyr)

set.seed(123)

# =========================================================
# MONTE CARLO CPT GENERATION
# FOR:
# Region → Rodent Disease → Mitigation → Time
# =========================================================

library(bnlearn)
library(gRain)
library(dplyr)

set.seed(123)

# =========================================================
# NODE STATES
# =========================================================

node_states <- list(
  
  Region_Context = c(
    "Low_Risk_Context",
    "Moderate_Risk_Context",
    "High_Risk_Context"
  ),
  
  Rodent_Borne_Diseases = c(
    "absent",
    "present"
  ),
  
  Habitat_Design = c(
    "not_implemented",
    "implemented"
  ),
  
  Public_Health_Awareness = c(
    "not_implemented",
    "implemented"
  ),
  
  Time_to_Stable_State = c(
    "1_3_years",
    "4_7_years",
    "8_12_years",
    "13_15_years"
  )
)

# =========================================================
# NETWORK STRUCTURE
# =========================================================

all_nodes <- c(
  "Region_Context",
  "Rodent_Borne_Diseases",
  "Habitat_Design",
  "Public_Health_Awareness",
  "Time_to_Stable_State"
)

bn <- empty.graph(all_nodes)

# =========================================================
# DEFINE ARCS
# =========================================================

arcs_df <- data.frame(
  
  from = c(
    "Region_Context",
    "Rodent_Borne_Diseases",
    "Rodent_Borne_Diseases",
    "Region_Context",
    "Habitat_Design",
    "Public_Health_Awareness"
  ),
  
  to = c(
    "Rodent_Borne_Diseases",
    "Habitat_Design",
    "Public_Health_Awareness",
    "Time_to_Stable_State",
    "Time_to_Stable_State",
    "Time_to_Stable_State"
  )
)

arcs(bn) <- as.matrix(arcs_df)

# =========================================================
# HELPER FUNCTION
# RANDOM PROBABILITY WITHIN RANGE
# =========================================================

sample_prob <- function(min, max){
  
  runif(1, min, max)
  
}

# =========================================================
# GENERATE RANDOM CPTS
# =========================================================

generate_random_cpts <- function(){
  
  # =====================================================
  # REGION PRIOR
  # =====================================================
  
  cpt_region <- cptable(
    
    ~Region_Context,
    
    values = c(
      0.33,
      0.34,
      0.33
    ),
    
    levels = node_states$Region_Context
  )
  
  # =====================================================
  # RODENT DISEASE CPT
  # REGION-SPECIFIC UNCERTAINTY
  # =====================================================
  
  low_present  <- sample_prob(0.10, 0.30)
  mod_present  <- sample_prob(0.30, 0.50)
  high_present <- sample_prob(0.50, 0.80)
  
  cpt_rodent <- cptable(
    
    ~Rodent_Borne_Diseases | Region_Context,
    
    values = c(
      
      # LOW RISK
      1-low_present,
      low_present,
      
      # MODERATE RISK
      1-mod_present,
      mod_present,
      
      # HIGH RISK
      1-high_present,
      high_present
    ),
    
    levels = node_states$Rodent_Borne_Diseases
  )
  
  # =====================================================
  # HABITAT DESIGN CPT
  # =====================================================
  
  habitat_impl_absent <- sample_prob(0.15, 0.35)
  habitat_impl_present <- sample_prob(0.60, 0.90)
  
  cpt_habitat <- cptable(
    
    ~Habitat_Design | Rodent_Borne_Diseases,
    
    values = c(
      
      # disease absent
      1-habitat_impl_absent,
      habitat_impl_absent,
      
      # disease present
      1-habitat_impl_present,
      habitat_impl_present
    ),
    
    levels = node_states$Habitat_Design
  )
  
  # =====================================================
  # PUBLIC HEALTH AWARENESS CPT
  # =====================================================
  
  awareness_absent <- sample_prob(0.10, 0.30)
  awareness_present <- sample_prob(0.65, 0.95)
  
  cpt_awareness <- cptable(
    
    ~Public_Health_Awareness |
      Rodent_Borne_Diseases,
    
    values = c(
      
      # disease absent
      1-awareness_absent,
      awareness_absent,
      
      # disease present
      1-awareness_present,
      awareness_present
    ),
    
    levels = node_states$Public_Health_Awareness
  )
  
  # =====================================================
  # TEMPORAL CPT
  #
  # Parents:
  # Region_Context = 3
  # Habitat_Design = 2
  # Public_Health_Awareness = 2
  #
  # TOTAL COMBINATIONS:
  # 3 × 2 × 2 = 12
  #
  # 12 × 4 outcome probabilities = 48
  # =====================================================
  
  temporal_values <- numeric(0)
  
  n_rows <- 12
  
  for(i in 1:n_rows){
    
    # ---------------------------------------------
    # Generate stabilization probabilities
    # ---------------------------------------------
    
    early <- sample_prob(0.01, 0.50)
    medium <- sample_prob(0.05, 0.40)
    delayed <- sample_prob(0.10, 0.45)
    
    final <- 1 - (early + medium + delayed)
    
    # Prevent negative values
    
    if(final < 0){
      
      probs <- c(early, medium, delayed)
      
      probs <- probs / sum(probs)
      
      early <- probs[1] * 0.85
      medium <- probs[2] * 0.85
      delayed <- probs[3] * 0.85
      
      final <- 0.15
    }
    
    row_probs <- c(
      early,
      medium,
      delayed,
      final
    )
    
    # Normalize
    
    row_probs <- row_probs / sum(row_probs)
    
    temporal_values <- c(
      temporal_values,
      row_probs
    )
  }
  
  # SAFETY CHECK
  
  print(length(temporal_values))
  
  # MUST = 48
  
  cpt_time <- cptable(
    
    ~Time_to_Stable_State |
      Region_Context :
      Habitat_Design :
      Public_Health_Awareness,
    
    values = temporal_values,
    
    levels = node_states$Time_to_Stable_State
  )
  
  # =====================================================
  # COMPILE CPTS
  # =====================================================
  
  plist <- compileCPT(list(
    
    cpt_region,
    cpt_rodent,
    cpt_habitat,
    cpt_awareness,
    cpt_time
  ))
  
  grain(plist)
}

# =========================================================
# MONTE CARLO SIMULATION
# =========================================================

n_iter <- 500

results <- data.frame()

for(i in 1:n_iter){
  
  bn_mc <- generate_random_cpts()
  
  # =====================================================
  # HIGH-RISK EVIDENCE
  # =====================================================
  
  bn_high <- setEvidence(
    
    bn_mc,
    
    nodes = "Region_Context",
    
    states = "High_Risk_Context"
  )
  
  q <- querygrain(
    
    bn_high,
    
    nodes = "Time_to_Stable_State"
    
  )$Time_to_Stable_State
  
  temp_df <- data.frame(
    
    iteration = i,
    
    years_1_3 = q["1_3_years"],
    years_4_7 = q["4_7_years"],
    years_8_12 = q["8_12_years"],
    years_13_15 = q["13_15_years"]
  )
  
  results <- bind_rows(
    results,
    temp_df
  )
}

# =========================================================
# SUMMARY TABLE
# =========================================================

summary_df <- results %>%
  
  summarise(
    
    mean_1_3 = mean(years_1_3),
    lower_1_3 = quantile(years_1_3, 0.025),
    upper_1_3 = quantile(years_1_3, 0.975),
    
    mean_4_7 = mean(years_4_7),
    lower_4_7 = quantile(years_4_7, 0.025),
    upper_4_7 = quantile(years_4_7, 0.975),
    
    mean_8_12 = mean(years_8_12),
    lower_8_12 = quantile(years_8_12, 0.025),
    upper_8_12 = quantile(years_8_12, 0.975),
    
    mean_13_15 = mean(years_13_15),
    lower_13_15 = quantile(years_13_15, 0.025),
    upper_13_15 = quantile(years_13_15, 0.975)
  )

print(summary_df)

# =========================================================
# EXAMPLE QUERY
# =========================================================

bn_example <- generate_random_cpts()

bn_highrisk <- setEvidence(
  
  bn_example,
  
  nodes = "Region_Context",
  
  states = "High_Risk_Context"
)

querygrain(
  
  bn_highrisk,
  
  nodes = c(
    "Rodent_Borne_Diseases",
    "Time_to_Stable_State"
  )
)
print(summary_df)

# =========================================================
# LONG FORMAT FOR PLOTTING
# =========================================================

plot_df <- results %>%
  
  pivot_longer(
    
    cols = starts_with("years"),
    
    names_to = "Time_Class",
    
    values_to = "Probability"
  )

# =========================================================
# UNCERTAINTY ENVELOPE PLOT
# =========================================================

library(ggplot2)

summary_plot <- plot_df %>%
  
  group_by(Time_Class) %>%
  
  summarise(
    
    mean_prob = mean(Probability),
    
    lower = quantile(Probability, 0.025),
    
    upper = quantile(Probability, 0.975)
  )

ggplot(
  summary_plot,
  aes(
    x = Time_Class,
    y = mean_prob,
    group = 1
  )
) +
  
  geom_line(size = 1.2) +
  
  geom_point(size = 3) +
  
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    alpha = 0.3
  ) +
  
  labs(
    
    title = "Monte Carlo Uncertainty Envelope",
    subtitle = "Probability of Stabilization Time Classes",
    x = "Temporal Stabilization Class",
    y = "Probability"
  ) +
  
  theme_minimal(base_size = 14)

# =========================================================
# EXPORT RESULTS
# =========================================================

write.csv(
  results,
  "MonteCarlo_BBN_Iterations.csv",
  row.names = FALSE
)

write.csv(
  summary_df,
  "MonteCarlo_BBN_Summary.csv",
  row.names = FALSE
)
