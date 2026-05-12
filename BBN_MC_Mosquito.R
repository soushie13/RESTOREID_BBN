# =========================================================
# MOSQUITO-BORNE DISEASE MODEL
# =========================================================

library(bnlearn)
library(gRain)
library(dplyr)
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
  
  Restoration = c(
    "no_restoration",
    "restoration"
  ),
  
  Wetland_Persistence = c(
    "low",
    "moderate",
    "high"
  ),
  
  Mosquito_Abundance = c(
    "low",
    "high"
  ),
  
  Mosquito_Borne_Diseases = c(
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
    
    "[Restoration]",
    
    "[Wetland_Persistence|Restoration:Region_Context]",
    
    "[Mosquito_Abundance|Wetland_Persistence]",
    
    "[Mosquito_Borne_Diseases|Mosquito_Abundance]",
    
    "[Time_to_Stable_State|Region_Context:Restoration:Mosquito_Borne_Diseases]"
    
  ))

# =========================================================
# HELPER FUNCTION
# =========================================================

sample_prob <- function(min,max){
  
  runif(1,min,max)
  
}

# =========================================================
# MONTE CARLO
# =========================================================

n_iter <- 1000

results <- data.frame()

for(i in 1:n_iter){
  
  # =======================================================
  # REGION
  # =======================================================
  
  cpt_region <- array(
    
    c(0.33,0.34,0.33),
    
    dim = 3,
    
    dimnames = list(
      Region_Context = states$Region_Context
    )
  )
  
  # =======================================================
  # RESTORATION
  # =======================================================
  
  cpt_restoration <- array(
    
    c(0.5,0.5),
    
    dim = 2,
    
    dimnames = list(
      Restoration = states$Restoration
    )
  )
  
  # =======================================================
  # WETLAND PERSISTENCE
  # =======================================================
  
  wetland_vals <- c()
  
  for(region in 1:3){
    
    for(rest in 1:2){
      
      if(rest == 2){
        
        probs <- c(
          sample_prob(0.05,0.20),
          sample_prob(0.25,0.40),
          sample_prob(0.45,0.70)
        )
        
      } else {
        
        probs <- c(
          sample_prob(0.45,0.70),
          sample_prob(0.20,0.35),
          sample_prob(0.05,0.15)
        )
      }
      
      probs <- probs / sum(probs)
      
      wetland_vals <- c(
        wetland_vals,
        probs
      )
    }
  }
  
  cpt_wetland <- array(
    
    wetland_vals,
    
    dim = c(3,2,3),
    
    dimnames = list(
      
      Wetland_Persistence =
        states$Wetland_Persistence,
      
      Restoration =
        states$Restoration,
      
      Region_Context =
        states$Region_Context
    )
  )
  
  # =======================================================
  # MOSQUITO ABUNDANCE
  # Wetlands may increase mosquitoes
  # =======================================================
  
  mosquito_vals <- c(
    
    0.80,0.20,
    0.50,0.50,
    0.20,0.80
    
  )
  
  cpt_mosquito <- array(
    
    mosquito_vals,
    
    dim = c(2,3),
    
    dimnames = list(
      
      Mosquito_Abundance =
        states$Mosquito_Abundance,
      
      Wetland_Persistence =
        states$Wetland_Persistence
    )
  )
  
  # =======================================================
  # DISEASE CPT
  # =======================================================
  
  disease_vals <- c(
    
    0.85,0.15,
    0.30,0.70
    
  )
  
  cpt_disease <- array(
    
    disease_vals,
    
    dim = c(2,2),
    
    dimnames = list(
      
      Mosquito_Borne_Diseases =
        states$Mosquito_Borne_Diseases,
      
      Mosquito_Abundance =
        states$Mosquito_Abundance
    )
  )
  
  # =======================================================
  # TEMPORAL CPT
  # =======================================================
  
  time_vals <- c()
  
  for(region in 1:3){
    
    for(rest in 1:2){
      
      for(disease in 1:2){
        
        # restoration but mosquito amplification
        
        if(rest == 2 & disease == 2){
          
          probs <- c(
            sample_prob(0.05,0.15),
            sample_prob(0.15,0.30),
            sample_prob(0.25,0.40),
            sample_prob(0.30,0.55)
          )
          
          # restoration + low disease
          
        } else if(rest == 2 & disease == 1){
          
          probs <- c(
            sample_prob(0.45,0.65),
            sample_prob(0.20,0.30),
            sample_prob(0.05,0.15),
            sample_prob(0.01,0.08)
          )
          
        } else {
          
          probs <- c(
            sample_prob(0.15,0.35),
            sample_prob(0.20,0.35),
            sample_prob(0.20,0.35),
            sample_prob(0.15,0.35)
          )
        }
        
        if(region == 3){
          
          probs[4] <- probs[4] + 0.15
          probs[1] <- probs[1] * 0.6
        }
        
        probs <- probs / sum(probs)
        
        time_vals <- c(time_vals, probs)
      }
    }
  }
  
  cpt_time <- array(
    
    time_vals,
    
    dim = c(4,3,2,2),
    
    dimnames = list(
      
      Time_to_Stable_State =
        states$Time_to_Stable_State,
      
      Region_Context =
        states$Region_Context,
      
      Restoration =
        states$Restoration,
      
      Mosquito_Borne_Diseases =
        states$Mosquito_Borne_Diseases
    )
  )
  
  # =======================================================
  # FIT NETWORK
  # =======================================================
  
  fitted_bn <- custom.fit(
    
    dag,
    
    dist = list(
      
      Region_Context =
        cpt_region,
      
      Restoration =
        cpt_restoration,
      
      Wetland_Persistence =
        cpt_wetland,
      
      Mosquito_Abundance =
        cpt_mosquito,
      
      Mosquito_Borne_Diseases =
        cpt_disease,
      
      Time_to_Stable_State =
        cpt_time
    )
  )
  
  bn_grain <- as.grain(fitted_bn)
  
  # =======================================================
  # RESTORATION SCENARIO
  # =======================================================
  
  bn_restore <- setEvidence(
    
    bn_grain,
    
    nodes = c(
      "Region_Context",
      "Restoration"
    ),
    
    states = c(
      "High_Risk_Context",
      "restoration"
    )
  )
  
  q_restore_time <- querygrain(
    bn_restore,
    nodes = "Time_to_Stable_State"
  )$Time_to_Stable_State
  
  q_restore_disease <- querygrain(
    bn_restore,
    nodes = "Mosquito_Borne_Diseases"
  )$Mosquito_Borne_Diseases
  
  # =======================================================
  # NO RESTORATION
  # =======================================================
  
  bn_no <- setEvidence(
    
    bn_grain,
    
    nodes = c(
      "Region_Context",
      "Restoration"
    ),
    
    states = c(
      "High_Risk_Context",
      "no_restoration"
    )
  )
  
  q_no_time <- querygrain(
    bn_no,
    nodes = "Time_to_Stable_State"
  )$Time_to_Stable_State
  
  q_no_disease <- querygrain(
    bn_no,
    nodes = "Mosquito_Borne_Diseases"
  )$Mosquito_Borne_Diseases
  
  # =======================================================
  # STORE RESULTS
  # =======================================================
  
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
