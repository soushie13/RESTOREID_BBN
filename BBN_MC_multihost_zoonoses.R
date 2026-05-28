# =========================================================
# DIRECT-CONTACT ZOONOTIC DISEASE MODEL
# MULTI-HOST SYSTEMS
# MONTE CARLO BAYESIAN BELIEF NETWORK
#
# Ecological assumptions:
# - Disease transmission driven by disturbance + contact rates
# - Multi-host amplification systems
# - Weak stabilization predictability
# - Broad uncertainty envelopes
# - Nonlinear restoration trajectories
#
# Suitable for:
# Nipah-like systems
# Hantavirus systems
# Bat-mammal#2-human spillover systems
# Multi-host generalist pathogens
# =========================================================

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
  
  Restoration = c(
    "No_Restoration",
    "Restoration"
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
    
    "[Restoration]",
    
    "[Habitat_Disturbance|Restoration:Region_Context]",
    
    "[Host_Diversity|Habitat_Disturbance]",
    
    "[Human_Wildlife_Contact|Habitat_Disturbance]",
    
    "[MultiHost_Amplification|Host_Diversity:Human_Wildlife_Contact]",
    
    "[Direct_Transmission_Risk|MultiHost_Amplification]",
    
    "[Time_to_Stable_State|Region_Context:Restoration:Direct_Transmission_Risk]"
    
  )
)

# =========================================================
# HELPER FUNCTION
# =========================================================

sample_prob <- function(min_val,max_val){
  
  runif(1,min_val,max_val)
  
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
  # HABITAT DISTURBANCE
  #
  # Restoration reduces disturbance
  # but uncertainty remains broad
  # =======================================================
  
  disturbance_vals <- c()
  
  for(region in 1:3){
    
    for(rest in 1:2){
      
      # -----------------------------------------------
      # RESTORATION
      # -----------------------------------------------
      
      if(rest == 2){
        
        probs <- c(
          
          sample_prob(0.35,0.65),
          sample_prob(0.20,0.40),
          sample_prob(0.05,0.25)
        )
        
      } else {
        
        probs <- c(
          
          sample_prob(0.05,0.20),
          sample_prob(0.20,0.35),
          sample_prob(0.45,0.75)
        )
      }
      
      # high-risk regional penalty
      
      if(region == 3){
        
        probs[3] <- probs[3] + 0.15
      }
      
      probs <- probs / sum(probs)
      
      disturbance_vals <- c(
        disturbance_vals,
        probs
      )
    }
  }
  
  cpt_disturbance <- array(
    
    disturbance_vals,
    
    dim = c(3,2,3),
    
    dimnames = list(
      
      Habitat_Disturbance =
        states$Habitat_Disturbance,
      
      Restoration =
        states$Restoration,
      
      Region_Context =
        states$Region_Context
    )
  )
  
  # =======================================================
  # HOST DIVERSITY
  #
  # Disturbance reduces diversity
  # =======================================================
  
  diversity_vals <- c(
    
    # low disturbance
    0.10,0.25,0.65,
    
    # moderate disturbance
    0.30,0.45,0.25,
    
    # high disturbance
    0.70,0.20,0.10
  )
  
  diversity_vals <- diversity_vals +
    runif(length(diversity_vals),-0.08,0.08)
  
  diversity_vals[diversity_vals < 0.01] <- 0.01
  
  diversity_vals <- unlist(
    
    lapply(
      
      split(
        diversity_vals,
        ceiling(seq_along(diversity_vals)/3)
      ),
      
      function(x) x/sum(x)
    )
  )
  
  cpt_diversity <- array(
    
    diversity_vals,
    
    dim = c(3,3),
    
    dimnames = list(
      
      Host_Diversity =
        states$Host_Diversity,
      
      Habitat_Disturbance =
        states$Habitat_Disturbance
    )
  )
  
  # =======================================================
  # HUMAN-WILDLIFE CONTACT
  #
  # Disturbance increases spillover interfaces
  # =======================================================
  
  contact_vals <- c(
    
    # low disturbance
    0.60,0.30,0.10,
    
    # moderate disturbance
    0.25,0.50,0.25,
    
    # high disturbance
    0.10,0.30,0.60
  )
  
  contact_vals <- contact_vals +
    runif(length(contact_vals),-0.10,0.10)
  
  contact_vals[contact_vals < 0.01] <- 0.01
  
  contact_vals <- unlist(
    
    lapply(
      
      split(
        contact_vals,
        ceiling(seq_along(contact_vals)/3)
      ),
      
      function(x) x/sum(x)
    )
  )
  
  cpt_contact <- array(
    
    contact_vals,
    
    dim = c(3,3),
    
    dimnames = list(
      
      Human_Wildlife_Contact =
        states$Human_Wildlife_Contact,
      
      Habitat_Disturbance =
        states$Habitat_Disturbance
    )
  )
  
  # =======================================================
  # MULTI-HOST AMPLIFICATION
  #
  # High contact + low diversity
  # strongly amplifies spillover
  #
  # Broad uncertainty intentionally included
  # =======================================================
  
  amplification_vals <- c()
  
  for(diversity in 1:3){
    
    for(contact in 1:3){
      
      # worst ecological conditions
      
      if(diversity == 1 & contact == 3){
        
        probs <- c(
          sample_prob(0.05,0.20),
          sample_prob(0.80,0.95)
        )
        
        # best ecological conditions
        
      } else if(diversity == 3 & contact == 1){
        
        probs <- c(
          sample_prob(0.70,0.90),
          sample_prob(0.10,0.30)
        )
        
        # transitional systems
        
      } else {
        
        probs <- c(
          sample_prob(0.25,0.70),
          sample_prob(0.30,0.75)
        )
      }
      
      probs <- probs / sum(probs)
      
      amplification_vals <- c(
        amplification_vals,
        probs
      )
    }
  }
  
  cpt_amplification <- array(
    
    amplification_vals,
    
    dim = c(2,3,3),
    
    dimnames = list(
      
      MultiHost_Amplification =
        states$MultiHost_Amplification,
      
      Host_Diversity =
        states$Host_Diversity,
      
      Human_Wildlife_Contact =
        states$Human_Wildlife_Contact
    )
  )
  
  # =======================================================
  # DIRECT TRANSMISSION RISK
  # =======================================================
  
  disease_vals <- c(
    
    # low amplification
    0.80,0.20,
    
    # high amplification
    0.20,0.80
  )
  
  disease_vals <- disease_vals +
    runif(length(disease_vals),-0.10,0.10)
  
  disease_vals[disease_vals < 0.01] <- 0.01
  
  disease_vals <- unlist(
    
    lapply(
      
      split(
        disease_vals,
        ceiling(seq_along(disease_vals)/2)
      ),
      
      function(x) x/sum(x)
    )
  )
  
  cpt_disease <- array(
    
    disease_vals,
    
    dim = c(2,2),
    
    dimnames = list(
      
      Direct_Transmission_Risk =
        states$Direct_Transmission_Risk,
      
      MultiHost_Amplification =
        states$MultiHost_Amplification
    )
  )
  
  # =======================================================
  # TEMPORAL TRAJECTORIES
  #
  # Weak stabilization predictability
  # Broad uncertainty corridors
  # =======================================================
  
  time_vals <- c()
  
  for(region in 1:3){
    
    for(rest in 1:2){
      
      for(disease in 1:2){
        
        # -----------------------------------------------
        # restoration + low disease
        # -----------------------------------------------
        
        if(rest == 2 & disease == 1){
          
          probs <- c(
            
            sample_prob(0.10,0.35),
            sample_prob(0.15,0.35),
            sample_prob(0.20,0.40),
            sample_prob(0.20,0.45)
          )
          
          # -----------------------------------------------
          # no restoration + high disease
          # -----------------------------------------------
          
        } else if(rest == 1 & disease == 2){
          
          probs <- c(
            
            sample_prob(0.15,0.40),
            sample_prob(0.15,0.35),
            sample_prob(0.15,0.35),
            sample_prob(0.15,0.40)
          )
          
          # -----------------------------------------------
          # highly uncertain transitions
          # -----------------------------------------------
          
        } else {
          
          probs <- c(
            
            sample_prob(0.10,0.40),
            sample_prob(0.10,0.40),
            sample_prob(0.10,0.40),
            sample_prob(0.10,0.40)
          )
        }
        
        # regional instability
        
        if(region == 3){
          
          probs[4] <- probs[4] + 0.10
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
      
      Direct_Transmission_Risk =
        states$Direct_Transmission_Risk
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
      
      Habitat_Disturbance =
        cpt_disturbance,
      
      Host_Diversity =
        cpt_diversity,
      
      Human_Wildlife_Contact =
        cpt_contact,
      
      MultiHost_Amplification =
        cpt_amplification,
      
      Direct_Transmission_Risk =
        cpt_disease,
      
      Time_to_Stable_State =
        cpt_time
    )
  )
  
  bn_grain <- as.grain(fitted_bn)
  
  # =======================================================
  # RESTORATION SCENARIOS
  # =======================================================
  
  for(restoration_state in states$Restoration){
    
    bn_temp <- setEvidence(
      
      bn_grain,
      
      nodes = c(
        "Region_Context",
        "Restoration"
      ),
      
      states = c(
        "High_Risk_Context",
        restoration_state
      )
    )
    
    disease_q <- querygrain(
      bn_temp,
      nodes = "Direct_Transmission_Risk"
    )$Direct_Transmission_Risk
    
    time_q <- querygrain(
      bn_temp,
      nodes = "Time_to_Stable_State"
    )$Time_to_Stable_State
    
    temp_df <- data.frame(
      
      iteration = i,
      
      scenario = restoration_state,
      
      disease_risk =
        disease_q["present"],
      
      years_1_3 =
        time_q["1_3_years"],
      
      years_4_7 =
        time_q["4_7_years"],
      
      years_8_12 =
        time_q["8_12_years"],
      
      years_13_15 =
        time_q["13_15_years"]
    )
    
    results <- bind_rows(
      results,
      temp_df
    )
  }
}

# =========================================================
# SUMMARY
# =========================================================

summary_df <- results %>%
  
  group_by(scenario) %>%
  
  summarise(
    
    mean_disease =
      mean(disease_risk),
    
    lower_disease =
      quantile(disease_risk,0.025),
    
    upper_disease =
      quantile(disease_risk,0.975),
    
    mean_13_15 =
      mean(years_13_15),
    
    lower_13_15 =
      quantile(years_13_15,0.025),
    
    upper_13_15 =
      quantile(years_13_15,0.975)
  )

print(summary_df)

# =========================================================
# TRADEOFF PLOT
# =========================================================

results$long_term_stability <-
  results$years_8_12 +
  results$years_13_15

ggplot(
  
  results,
  
  aes(
    x = disease_risk,
    y = long_term_stability,
    color = scenario
  )
  
) +
  
  geom_point(alpha = 0.2) +
  
  geom_smooth(
    method = "loess",
    se = TRUE
  ) +
  
  labs(
    
    title =
      "Direct-Transmission Zoonoses",
    
    x =
      "Disease Risk Probability",
    
    y =
      "Long-Term Stabilization Probability"
  ) +
  
  theme_minimal(base_size = 14)

traj_long <- results %>%
  
  pivot_longer(
    
    cols = starts_with("years_"),
    
    names_to = "time_period",
    
    values_to = "probability"
  )

# =========================================================
# ORDER TEMPORAL AXIS
# =========================================================

traj_long$time_period <- factor(
  
  traj_long$time_period,
  
  levels = c(
    "years_1_3",
    "years_4_7",
    "years_8_12",
    "years_13_15"
  ),
  
  labels = c(
    "1-3 years",
    "4-7 years",
    "8-12 years",
    "13-15 years"
  )
)

# =========================================================
# SUMMARISE UNCERTAINTY ENVELOPES
# =========================================================

corridor_df <- traj_long %>%
  
  group_by(
    scenario,
    time_period
  ) %>%
  
  summarise(
    
    mean_prob =
      mean(probability, na.rm = TRUE),
    
    lower =
      quantile(
        probability,
        0.025,
        na.rm = TRUE
      ),
    
    upper =
      quantile(
        probability,
        0.975,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )

# =========================================================
# UNCERTAINTY CORRIDOR PLOT
# =========================================================

ggplot(
  
  corridor_df,
  
  aes(
    x = time_period,
    y = mean_prob,
    group = scenario,
    color = scenario,
    fill = scenario
  )
) +
  
  geom_ribbon(
    
    aes(
      ymin = lower,
      ymax = upper
    ),
    
    alpha = 0.20,
    
    color = NA
  ) +
  
  geom_line(
    
    linewidth = 1.5
  ) +
  
  geom_point(
    
    size = 3
  ) +
  
  labs(
    
    title =
      "Direct-Contact Zoonoses:\nStabilization Trajectories with Uncertainty Corridors",
    
    x =
      "Time to Ecological Stabilization",
    
    y =
      "Probability",
    
    color =
      "Scenario",
    
    fill =
      "Scenario"
  ) +
  
  theme_minimal(base_size = 14)

# =========================================================
# DISEASE RISK DISTRIBUTION PLOT
# =========================================================

ggplot(
  
  results,
  
  aes(
    x = disease_risk,
    fill = scenario
  )
) +
  
  geom_density(
    
    alpha = 0.35
  ) +
  
  labs(
    
    title =
      "Disease Risk Uncertainty Distribution\nDirect-Contact Zoonoses",
    
    x =
      "Disease Risk Probability",
    
    y =
      "Density",
    
    fill =
      "Scenario"
  ) +
  
  theme_minimal(base_size = 14)
