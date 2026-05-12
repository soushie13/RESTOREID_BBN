# =========================================================
# AVIAN DISEASE MODEL
# Broad uncertainty envelopes
# Strong regional sensitivity
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
  
  Region_Context = c(
    "Low_Risk_Context",
    "Moderate_Risk_Context",
    "High_Risk_Context"
  ),
  
  Restoration = c(
    "no_restoration",
    "restoration"
  ),
  
  Bird_Aggregation = c(
    "low",
    "moderate",
    "high"
  ),
  
  Migratory_Connectivity = c(
    "low",
    "high"
  ),
  
  Avian_Diseases = c(
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
    
    "[Bird_Aggregation|Restoration:Region_Context]",
    
    "[Migratory_Connectivity|Bird_Aggregation]",
    
    "[Avian_Diseases|Migratory_Connectivity]",
    
    "[Time_to_Stable_State|Region_Context:Restoration:Avian_Diseases]"
    
  )
  
)

# =========================================================
# HELPER FUNCTION
# =========================================================

sample_prob <- function(min,max){
  
  runif(1,min,max)
  
}

# =========================================================
# MONTE CARLO SETTINGS
# =========================================================

n_iter <- 1000

results <- data.frame()

# =========================================================
# MONTE CARLO LOOP
# =========================================================

for(i in 1:n_iter){
  
  # =======================================================
  # REGION PRIOR
  # =======================================================
  
  cpt_region <- array(
    
    c(0.33,0.34,0.33),
    
    dim = 3,
    
    dimnames = list(
      Region_Context = states$Region_Context
    )
  )
  
  # =======================================================
  # RESTORATION PRIOR
  # =======================================================
  
  cpt_restoration <- array(
    
    c(0.5,0.5),
    
    dim = 2,
    
    dimnames = list(
      Restoration = states$Restoration
    )
  )
  
  # =======================================================
  # BIRD AGGREGATION CPT
  # Restoration can increase bird congregation
  # especially in high-risk migratory regions
  # =======================================================
  
  aggregation_vals <- c()
  
  for(region in 1:3){
    
    for(rest in 1:2){
      
      # restoration
      
      if(rest == 2){
        
        probs <- c(
          
          sample_prob(0.05,0.20),
          sample_prob(0.20,0.35),
          sample_prob(0.45,0.70)
          
        )
        
        # no restoration
        
      } else {
        
        probs <- c(
          
          sample_prob(0.45,0.70),
          sample_prob(0.20,0.35),
          sample_prob(0.05,0.15)
          
        )
      }
      
      # ===================================================
      # STRONG REGIONAL SENSITIVITY
      # High-risk regions amplify aggregation uncertainty
      # ===================================================
      
      if(region == 3){
        
        probs[3] <- probs[3] + sample_prob(0.10,0.25)
        
      }
      
      probs <- probs / sum(probs)
      
      aggregation_vals <- c(
        aggregation_vals,
        probs
      )
    }
  }
  
  cpt_aggregation <- array(
    
    aggregation_vals,
    
    dim = c(3,2,3),
    
    dimnames = list(
      
      Bird_Aggregation =
        states$Bird_Aggregation,
      
      Restoration =
        states$Restoration,
      
      Region_Context =
        states$Region_Context
    )
  )
  
  # =======================================================
  # MIGRATORY CONNECTIVITY CPT
  # Aggregation increases migratory mixing
  # =======================================================
  
  connectivity_vals <- c(
    
    # low aggregation
    0.80,0.20,
    
    # moderate aggregation
    0.45,0.55,
    
    # high aggregation
    0.10,0.90
  )
  
  # extra uncertainty
  
  connectivity_vals <- connectivity_vals +
    runif(length(connectivity_vals), -0.12, 0.12)
  
  connectivity_vals[connectivity_vals < 0.01] <- 0.01
  
  connectivity_vals <- unlist(
    
    lapply(
      
      split(
        connectivity_vals,
        ceiling(seq_along(connectivity_vals)/2)
      ),
      
      function(x){
        
        x / sum(x)
        
      }
    )
  )
  
  cpt_connectivity <- array(
    
    connectivity_vals,
    
    dim = c(2,3),
    
    dimnames = list(
      
      Migratory_Connectivity =
        states$Migratory_Connectivity,
      
      Bird_Aggregation =
        states$Bird_Aggregation
    )
  )
  
  # =======================================================
  # AVIAN DISEASE CPT
  # Strong uncertainty envelopes
  # =======================================================
  
  disease_vals <- c(
    
    # low connectivity
    0.75,0.25,
    
    # high connectivity
    0.30,0.70
  )
  
  # LARGE uncertainty perturbation
  
  disease_vals <- disease_vals +
    runif(length(disease_vals), -0.18, 0.18)
  
  disease_vals[disease_vals < 0.01] <- 0.01
  
  disease_vals <- unlist(
    
    lapply(
      
      split(
        disease_vals,
        ceiling(seq_along(disease_vals)/2)
      ),
      
      function(x){
        
        x / sum(x)
        
      }
    )
  )
  
  cpt_disease <- array(
    
    disease_vals,
    
    dim = c(2,2),
    
    dimnames = list(
      
      Avian_Diseases =
        states$Avian_Diseases,
      
      Migratory_Connectivity =
        states$Migratory_Connectivity
    )
  )
  
  # =======================================================
  # TEMPORAL CPT
  # Broad uncertainty envelopes
  # =======================================================
  
  time_vals <- c()
  
  for(region in 1:3){
    
    for(rest in 1:2){
      
      for(disease in 1:2){
        
        # =================================================
        # Restoration + disease present
        # Large uncertainty spread
        # =================================================
        
        if(rest == 2 & disease == 2){
          
          probs <- c(
            
            sample_prob(0.05,0.30),
            sample_prob(0.10,0.35),
            sample_prob(0.15,0.40),
            sample_prob(0.20,0.50)
          )
          
          # =================================================
          # Restoration + disease absent
          # =================================================
          
        } else if(rest == 2 & disease == 1){
          
          probs <- c(
            
            sample_prob(0.20,0.50),
            sample_prob(0.15,0.35),
            sample_prob(0.05,0.25),
            sample_prob(0.01,0.15)
          )
          
          # =================================================
          # No restoration + disease present
          # =================================================
          
        } else if(rest == 1 & disease == 2){
          
          probs <- c(
            
            sample_prob(0.05,0.20),
            sample_prob(0.10,0.30),
            sample_prob(0.20,0.40),
            sample_prob(0.25,0.55)
          )
          
          # =================================================
          # Intermediate uncertainty
          # =================================================
          
        } else {
          
          probs <- c(
            
            sample_prob(0.10,0.35),
            sample_prob(0.15,0.35),
            sample_prob(0.15,0.35),
            sample_prob(0.10,0.35)
          )
        }
        
        # =================================================
        # STRONG REGIONAL SENSITIVITY
        # =================================================
        
        if(region == 3){
          
          probs[4] <- probs[4] + sample_prob(0.10,0.25)
          
          probs[1] <- probs[1] * sample_prob(0.40,0.75)
          
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
        states$Time_to_Stable_State,
      
      Region_Context =
        states$Region_Context,
      
      Restoration =
        states$Restoration,
      
      Avian_Diseases =
        states$Avian_Diseases
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
      
      Bird_Aggregation =
        cpt_aggregation,
      
      Migratory_Connectivity =
        cpt_connectivity,
      
      Avian_Diseases =
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
    nodes = "Avian_Diseases"
  )$Avian_Diseases
  
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
    nodes = "Avian_Diseases"
  )$Avian_Diseases
  
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
# LONG-TERM INSTABILITY
# =========================================================

results$long_term_instability <-
  results$years_8_12 +
  results$years_13_15

# =========================================================
# TRADEOFF CURVES
# =========================================================

ggplot(
  
  results,
  
  aes(
    x = disease_risk,
    y = long_term_instability,
    color = scenario
  )
  
) +
  
  geom_point(alpha = 0.15) +
  
  geom_smooth(
    method = "loess",
    se = TRUE
  ) +
  
  labs(
    
    title =
      "Avian Disease Risk vs Stabilization Tradeoff",
    
    x =
      "Avian Disease Risk Probability",
    
    y =
      "Long-Term Stabilization Probability"
    
  ) +
  
  theme_minimal()

# =========================================================
# UNCERTAINTY CORRIDOR PLOT
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

plot_df <- data.frame(
  
  scenario = rep(summary_df$scenario, each = 4),
  
  time_period = rep(
    c("1_3","4_7","8_12","13_15"),
    times = 2
  ),
  
  mean = c(
    summary_df$mean_1_3,
    summary_df$mean_4_7,
    summary_df$mean_8_12,
    summary_df$mean_13_15
  ),
  
  lower = c(
    summary_df$lower_1_3,
    summary_df$lower_4_7,
    summary_df$lower_8_12,
    summary_df$lower_13_15
  ),
  
  upper = c(
    summary_df$upper_1_3,
    summary_df$upper_4_7,
    summary_df$upper_8_12,
    summary_df$upper_13_15
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
    color = scenario,
    fill = scenario,
    group = scenario
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
    
    title =
      "Avian Disease Stabilization Trajectories",
    
    x =
      "Stabilization Horizon",
    
    y =
      "Probability"
    
  ) +
  
  theme_minimal()
