library(bnlearn)
library(gRain)
library(dplyr)

set.seed(123)

# =========================================================
# HELPER FUNCTION
# =========================================================

sample_prob <- function(min,max){
  
  runif(1,min,max)
  
}

# =========================================================
# NODE STATES
# =========================================================

node_states <- list(
  
  Region_Context = c(
    "Low_Risk_Context",
    "Moderate_Risk_Context",
    "High_Risk_Context"
  ),
  
  Habitat_Design = c(
    "not_restored",
    "restored"
  ),
  
  Biodiversity_Recovery = c(
    "low",
    "moderate",
    "high"
  ),
  
  Rodent_Abundance = c(
    "low",
    "high"
  ),
  
  Rodent_Borne_Diseases = c(
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
# DAG STRUCTURE
# =========================================================

dag <- model2network(
  
  paste0(
    
    "[Region_Context]",
    
    "[Habitat_Design]",
    
    "[Biodiversity_Recovery|Habitat_Design:Region_Context]",
    
    "[Rodent_Abundance|Biodiversity_Recovery]",
    
    "[Rodent_Borne_Diseases|Rodent_Abundance]",
    
    "[Time_to_Stable_State|Region_Context:Habitat_Design:Rodent_Borne_Diseases]"
    
  ))

# =========================================================
# MONTE CARLO SETTINGS
# =========================================================

n_iter <- 1000

results <- data.frame()

# =========================================================
# MONTE CARLO LOOP
# =========================================================

for(i in 1:n_iter){
  
  # ======================================================
  # REGION PRIOR
  # ======================================================
  
  cpt_region <- array(
    
    c(0.33,0.34,0.33),
    
    dim = c(3),
    
    dimnames = list(
      
      Region_Context =
        node_states$Region_Context
    )
  )
  
  # ======================================================
  # HABITAT PRIOR
  # ======================================================
  
  cpt_habitat <- array(
    
    c(0.5,0.5),
    
    dim = c(2),
    
    dimnames = list(
      
      Habitat_Design =
        node_states$Habitat_Design
    )
  )
  
  # ======================================================
  # BIODIVERSITY CPT
  # ======================================================
  
  biodiversity_vals <- c()
  
  for(region in 1:3){
    
    for(habitat in 1:2){
      
      # restored
      if(habitat == 2){
        
        probs <- c(
          
          sample_prob(0.10,0.25),
          sample_prob(0.25,0.40),
          sample_prob(0.45,0.65)
          
        )
        
      } else {
        
        probs <- c(
          
          sample_prob(0.55,0.75),
          sample_prob(0.15,0.30),
          sample_prob(0.03,0.12)
          
        )
      }
      
      # regional degradation penalty
      if(region == 3){
        
        probs[1] <- probs[1] + 0.10
        
      }
      
      probs <- probs / sum(probs)
      
      biodiversity_vals <- c(
        biodiversity_vals,
        probs
      )
    }
  }
  
  cpt_biodiversity <- array(
    
    biodiversity_vals,
    
    dim = c(3,2,3),
    
    dimnames = list(
      
      Biodiversity_Recovery =
        node_states$Biodiversity_Recovery,
      
      Habitat_Design =
        node_states$Habitat_Design,
      
      Region_Context =
        node_states$Region_Context
    )
  )
  
  # ======================================================
  # RODENT ABUNDANCE CPT
  # ======================================================
  
  rodent_vals <- c(
    
    # low biodiversity
    0.20,0.80,
    
    # moderate biodiversity
    0.50,0.50,
    
    # high biodiversity
    0.80,0.20
  )
  
  cpt_rodents <- array(
    
    rodent_vals,
    
    dim = c(2,3),
    
    dimnames = list(
      
      Rodent_Abundance =
        node_states$Rodent_Abundance,
      
      Biodiversity_Recovery =
        node_states$Biodiversity_Recovery
    )
  )
  
  # ======================================================
  # DISEASE CPT
  # ======================================================
  
  disease_vals <- c(
    
    # rodents low
    0.80,0.20,
    
    # rodents high
    0.25,0.75
  )
  
  cpt_disease <- array(
    
    disease_vals,
    
    dim = c(2,2),
    
    dimnames = list(
      
      Rodent_Borne_Diseases =
        node_states$Rodent_Borne_Diseases,
      
      Rodent_Abundance =
        node_states$Rodent_Abundance
    )
  )
  
  # ======================================================
  # TEMPORAL CPT
  # ======================================================
  
  time_vals <- c()
  
  for(region in 1:3){
    
    for(habitat in 1:2){
      
      for(disease in 1:2){
        
        # ================================================
        # BEST CONDITIONS
        # ================================================
        
        if(habitat == 2 & disease == 1){
          
          probs <- c(
            
            sample_prob(0.45,0.70),
            sample_prob(0.15,0.30),
            sample_prob(0.05,0.15),
            sample_prob(0.01,0.08)
          )
          
          # ================================================
          # WORST CONDITIONS
          # ================================================
          
        } else if(habitat == 1 & disease == 2){
          
          probs <- c(
            
            sample_prob(0.01,0.08),
            sample_prob(0.05,0.18),
            sample_prob(0.20,0.35),
            sample_prob(0.45,0.70)
          )
          
          # ================================================
          # INTERMEDIATE
          # ================================================
          
        } else {
          
          probs <- c(
            
            sample_prob(0.15,0.35),
            sample_prob(0.20,0.35),
            sample_prob(0.20,0.35),
            sample_prob(0.15,0.35)
          )
        }
        
        # high-risk regional penalty
        if(region == 3){
          
          probs[4] <- probs[4] + 0.15
          probs[1] <- probs[1] * 0.6
          
        }
        
        probs <- probs / sum(probs)
        
        time_vals <- c(
          time_vals,
          probs
        )
      }
    }
  }
  
  cpt_time <- array(
    
    time_vals,
    
    dim = c(4,3,2,2),
    
    dimnames = list(
      
      Time_to_Stable_State =
        node_states$Time_to_Stable_State,
      
      Region_Context =
        node_states$Region_Context,
      
      Habitat_Design =
        node_states$Habitat_Design,
      
      Rodent_Borne_Diseases =
        node_states$Rodent_Borne_Diseases
    )
  )
  
  # ======================================================
  # FIT NETWORK
  # ======================================================
  
  fitted_bn <- custom.fit(
    
    dag,
    
    dist = list(
      
      Region_Context = cpt_region,
      
      Habitat_Design = cpt_habitat,
      
      Biodiversity_Recovery = cpt_biodiversity,
      
      Rodent_Abundance = cpt_rodents,
      
      Rodent_Borne_Diseases = cpt_disease,
      
      Time_to_Stable_State = cpt_time
    )
  )
  
  # ======================================================
  # CONVERT TO GRAIN
  # ======================================================
  
  bn_grain <- as.grain(fitted_bn)
  
  # ======================================================
  # RESTORATION SCENARIO
  # ======================================================
  
  bn_restore <- setEvidence(
    
    bn_grain,
    
    nodes = c(
      "Region_Context",
      "Habitat_Design"
    ),
    
    states = c(
      "High_Risk_Context",
      "restored"
    )
  )
  
  q_restore <- querygrain(
    bn_restore,
    nodes = "Time_to_Stable_State"
  )$Time_to_Stable_State
 
   q_restore_time <- querygrain(
    bn_restore,
    nodes = "Time_to_Stable_State"
  )$Time_to_Stable_State
  
  q_restore_disease <- querygrain(
    bn_restore,
    nodes = "Rodent_Borne_Diseases"
  )$Rodent_Borne_Diseases
  # ======================================================
  # NO RESTORATION SCENARIO
  # ======================================================
  
  bn_no_restore <- setEvidence(
    
    bn_grain,
    
    nodes = c(
      "Region_Context",
      "Habitat_Design"
    ),
    
    states = c(
      "High_Risk_Context",
      "not_restored"
    )
  )
  
  q_no_restore <- querygrain(
    bn_no_restore,
    nodes = "Time_to_Stable_State"
  )$Time_to_Stable_State
  
  q_no_time <- querygrain(
    bn_no_restore,
    nodes = "Time_to_Stable_State"
  )$Time_to_Stable_State
  
  q_no_disease <- querygrain(
    bn_no_restore,
    nodes = "Rodent_Borne_Diseases"
  )$Rodent_Borne_Diseases
  # ======================================================
  # STORE RESULTS
  # ======================================================
  
  temp_df <- data.frame(
    
    iteration = c(i,i),
    
    scenario = c(
      "Restoration",
      "No_Restoration"
    ),
    
    disease_risk = c(
      q_restore_disease["present"],
      q_no_disease["present"]
    ),
    
    years_1_3 = c(
      q_restore_time["1_3_years"],
      q_no_time["1_3_years"]
    ),
    
    years_4_7 = c(
      q_restore_time["4_7_years"],
      q_no_time["4_7_years"]
    ),
    
    years_8_12 = c(
      q_restore_time["8_12_years"],
      q_no_time["8_12_years"]
    ),
    
    years_13_15 = c(
      q_restore_time["13_15_years"],
      q_no_time["13_15_years"]
    )
  )
  
  results <- bind_rows(
    results,
    temp_df
  )
}

# =========================================================
# SUMMARY
# =========================================================

summary_df <- results %>%
  
  group_by(scenario) %>%
  
  summarise(
    
    mean_1_3 = mean(years_1_3),
    lower_1_3 = quantile(years_1_3,0.025),
    upper_1_3 = quantile(years_1_3,0.975),
    
    mean_4_7 = mean(years_4_7),
    lower_4_7 = quantile(years_4_7,0.025),
    upper_4_7 = quantile(years_4_7,0.975),
    
    mean_8_12 = mean(years_8_12),
    lower_8_12 = quantile(years_8_12,0.025),
    upper_8_12 = quantile(years_8_12,0.975),
    
    mean_13_15 = mean(years_13_15),
    lower_13_15 = quantile(years_13_15,0.025),
    upper_13_15 = quantile(years_13_15,0.975)
  )

print(summary_df)

# =========================================================
# TRADEOFF CURVE PLOT
# =========================================================
results$long_term_instability <-
  results$years_8_12 +
  results$years_13_15

library(ggplot2)

ggplot(
  
  results,
  
  aes(
    x = disease_risk,
    y = long_term_instability,
    color = scenario
  )
  
) +
  
  geom_point(
    alpha = 0.2
  ) +
  
  geom_smooth(
    method = "loess",
    se = TRUE
  ) +
  
  labs(
    
    x = "Disease Risk Probability",
    y = "Long-Term Stabilization Probability",
    title = "Disease–Restoration Tradeoff Curves"
    
  ) +
  
  theme_minimal()

# =========================================================
# UNCERTAINTY CORRIDOR PLOTS
# =========================================================

trajectory_summary <- results %>%
  
  group_by(scenario) %>%
  
  summarise(
    
    mean_1_3 = mean(years_1_3),
    low_1_3 = quantile(years_1_3,0.025),
    high_1_3 = quantile(years_1_3,0.975),
    
    mean_4_7 = mean(years_4_7),
    low_4_7 = quantile(years_4_7,0.025),
    high_4_7 = quantile(years_4_7,0.975),
    
    mean_8_12 = mean(years_8_12),
    low_8_12 = quantile(years_8_12,0.025),
    high_8_12 = quantile(years_8_12,0.975),
    
    mean_13_15 = mean(years_13_15),
    low_13_15 = quantile(years_13_15,0.025),
    high_13_15 = quantile(years_13_15,0.975)
  )

library(tidyr)

plot_df <- data.frame(
  
  scenario = rep(
    trajectory_summary$scenario,
    each = 4
  ),
  
  time_period = rep(
    c("1_3","4_7","8_12","13_15"),
    times = 2
  ),
  
  mean = c(
    trajectory_summary$mean_1_3,
    trajectory_summary$mean_4_7,
    trajectory_summary$mean_8_12,
    trajectory_summary$mean_13_15
  ),
  
  lower = c(
    trajectory_summary$low_1_3,
    trajectory_summary$low_4_7,
    trajectory_summary$low_8_12,
    trajectory_summary$low_13_15
  ),
  
  upper = c(
    trajectory_summary$high_1_3,
    trajectory_summary$high_4_7,
    trajectory_summary$high_8_12,
    trajectory_summary$high_13_15
  )
)

plot_df$time_period <- factor(
  
  plot_df$time_period,
  
  levels = c(
    "1_3",
    "4_7",
    "8_12",
    "13_15"
  )
)

ggplot(
  
  plot_df,
  
  aes(
    x = time_period,
    y = mean,
    group = scenario,
    color = scenario,
    fill = scenario
  )
  
) +
  
  geom_line(size = 1.2) +
  
  geom_ribbon(
    
    aes(
      ymin = lower,
      ymax = upper
    ),
    
    alpha = 0.2,
    color = NA
  ) +
  
  labs(
    
    x = "Stabilization Horizon",
    y = "Probability",
    title = "Restoration Stabilization Trajectories with Uncertainty Corridors"
    
  ) +
  
  theme_minimal()
