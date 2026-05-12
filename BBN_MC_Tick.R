# =========================================================
# TICK-BORNE DISEASE MODEL
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
  
  Host_Connectivity = c(
    "low",
    "moderate",
    "high"
  ),
  
  Tick_Abundance = c(
    "low",
    "high"
  ),
  
  Tick_Borne_Diseases = c(
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
    
    "[Host_Connectivity|Restoration:Region_Context]",
    
    "[Tick_Abundance|Host_Connectivity]",
    
    "[Tick_Borne_Diseases|Tick_Abundance]",
    
    "[Time_to_Stable_State|Region_Context:Restoration:Tick_Borne_Diseases]"
    
  ))

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
  # HOST CONNECTIVITY CPT
  # Restoration increases connectivity
  # but with uncertainty
  # =======================================================
  
  connectivity_vals <- c()
  
  for(region in 1:3){
    
    for(rest in 1:2){
      
      # restoration
      
      if(rest == 2){
        
        probs <- c(
          
          sample_prob(0.05,0.20),
          sample_prob(0.20,0.40),
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
      
      # high-risk region amplifies connectivity
      
      if(region == 3){
        
        probs[3] <- probs[3] + 0.15
        
      }
      
      probs <- probs / sum(probs)
      
      connectivity_vals <- c(
        connectivity_vals,
        probs
      )
    }
  }
  
  cpt_connectivity <- array(
    
    connectivity_vals,
    
    dim = c(3,2,3),
    
    dimnames = list(
      
      Host_Connectivity =
        states$Host_Connectivity,
      
      Restoration =
        states$Restoration,
      
      Region_Context =
        states$Region_Context
    )
  )
  
  # =======================================================
  # TICK ABUNDANCE CPT
  # Nonlinear response
  # Moderate connectivity can still produce ticks
  # =======================================================
  
  tick_vals <- c(
    
    # LOW connectivity
    0.80,0.20,
    
    # MODERATE connectivity
    0.45,0.55,
    
    # HIGH connectivity
    0.15,0.85
  )
  
  tick_vals <- tick_vals +
    runif(length(tick_vals), -0.08, 0.08)
  
  tick_vals[tick_vals < 0.01] <- 0.01
  
  tick_vals <- unlist(
    lapply(
      split(
        tick_vals,
        ceiling(seq_along(tick_vals)/2)
      ),
      function(x) x/sum(x)
    )
  )
  
  cpt_ticks <- array(
    
    tick_vals,
    
    dim = c(2,3),
    
    dimnames = list(
      
      Tick_Abundance =
        states$Tick_Abundance,
      
      Host_Connectivity =
        states$Host_Connectivity
    )
  )
  
  # =======================================================
  # DISEASE CPT
  # Tick abundance drives disease
  # =======================================================
  
  disease_vals <- c(
    
    # low ticks
    0.85,0.15,
    
    # high ticks
    0.20,0.80
  )
  
  disease_vals <- disease_vals +
    runif(length(disease_vals), -0.08, 0.08)
  
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
      
      Tick_Borne_Diseases =
        states$Tick_Borne_Diseases,
      
      Tick_Abundance =
        states$Tick_Abundance
    )
  )
  
  # =======================================================
  # TEMPORAL CPT
  # Delayed equilibrium dynamics
  # =======================================================
  
  time_vals <- c()
  
  for(region in 1:3){
    
    for(rest in 1:2){
      
      for(disease in 1:2){
        
        # =================================================
        # Restoration but disease persists
        # Delayed ecological equilibrium
        # =================================================
        
        if(rest == 2 & disease == 2){
          
          probs <- c(
            
            sample_prob(0.03,0.10),
            sample_prob(0.10,0.22),
            sample_prob(0.25,0.40),
            sample_prob(0.35,0.60)
          )
          
          # =================================================
          # Restoration + disease absent
          # =================================================
          
        } else if(rest == 2 & disease == 1){
          
          probs <- c(
            
            sample_prob(0.25,0.45),
            sample_prob(0.25,0.35),
            sample_prob(0.15,0.30),
            sample_prob(0.05,0.15)
          )
          
          # =================================================
          # No restoration + disease present
          # =================================================
          
        } else if(rest == 1 & disease == 2){
          
          probs <- c(
            
            sample_prob(0.05,0.15),
            sample_prob(0.15,0.25),
            sample_prob(0.25,0.35),
            sample_prob(0.30,0.55)
          )
          
          # =================================================
          # Intermediate uncertainty
          # =================================================
          
        } else {
          
          probs <- c(
            
            sample_prob(0.15,0.30),
            sample_prob(0.20,0.35),
            sample_prob(0.20,0.35),
            sample_prob(0.15,0.30)
          )
        }
        
        # =================================================
        # High-risk region penalty
        # =================================================
        
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
        states$Time_to_Stable_State,
      
      Region_Context =
        states$Region_Context,
      
      Restoration =
        states$Restoration,
      
      Tick_Borne_Diseases =
        states$Tick_Borne_Diseases
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
      
      Host_Connectivity =
        cpt_connectivity,
      
      Tick_Abundance =
        cpt_ticks,
      
      Tick_Borne_Diseases =
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
    nodes = "Tick_Borne_Diseases"
  )$Tick_Borne_Diseases
  
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
    nodes = "Tick_Borne_Diseases"
  )$Tick_Borne_Diseases
  
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
# TRADEOFF CURVE
# =========================================================

ggplot(
  
  results,
  
  aes(
    x = disease_risk,
    y = long_term_instability,
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
      "Tick-Borne Disease Tradeoff Curves",
    
    x =
      "Disease Risk Probability",
    
    y =
      "Long-Term Stabilization Probability"
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
# =========================================================
# TICK-BORNE DISEASES
# RESTORATION INTENSITY GRADIENT
# MONTE CARLO BBN
# =========================================================

library(bnlearn)
library(gRain)
library(dplyr)
library(tidyr)
library(ggplot2)

set.seed(123)

# =========================================================
# NODE STATES
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
  
  Host_Connectivity = c(
    "poor",
    "moderate",
    "optimal"
  ),
  
  Tick_Abundance = c(
    "low",
    "moderate",
    "high"
  ),
  
  Tick_Borne_Diseases = c(
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
# NETWORK STRUCTURE
# =========================================================

dag <- model2network(
  
  paste0(
    
    "[Region_Context]",
    
    "[Restoration_Intensity]",
    
    "[Host_Connectivity|Region_Context:Restoration_Intensity]",
    
    "[Tick_Abundance|Host_Connectivity]",
    
    "[Tick_Borne_Diseases|Tick_Abundance]",
    
    "[Time_to_Stable_State|Region_Context:Restoration_Intensity:Tick_Borne_Diseases]"
    
  ))

# =========================================================
# HELPER FUNCTION
# =========================================================

sample_prob <- function(min_val, max_val){
  
  runif(1, min_val, max_val)
  
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
    
    c(0.25,0.25,0.25,0.25),
    
    dim = 4,
    
    dimnames = list(
      Restoration_Intensity =
        states$Restoration_Intensity
    )
  )
  
  # =======================================================
  # HOST connectivity CPT
  # =======================================================
  
  connect_vals <- c()
  
  for(region in 1:3){
    
    for(restoration in 1:4){
      
      # ---------------------------------------------
      # NONE
      # ---------------------------------------------
      
      if(restoration == 1){
        
        probs <- c(
          
          sample_prob(0.60,0.80),
          sample_prob(0.15,0.30),
          sample_prob(0.01,0.10)
        )
        
        # ---------------------------------------------
        # LOW
        # ---------------------------------------------
        
      } else if(restoration == 2){
        
        probs <- c(
          
          sample_prob(0.35,0.55),
          sample_prob(0.25,0.40),
          sample_prob(0.10,0.25)
        )
        
        # ---------------------------------------------
        # MODERATE
        # ---------------------------------------------
        
      } else if(restoration == 3){
        
        probs <- c(
          
          sample_prob(0.15,0.35),
          sample_prob(0.30,0.45),
          sample_prob(0.25,0.45)
        )
        
        # ---------------------------------------------
        # HIGH
        # ---------------------------------------------
        
      } else {
        
        probs <- c(
          
          sample_prob(0.05,0.20),
          sample_prob(0.20,0.35),
          sample_prob(0.45,0.70)
        )
      }
      
      # ---------------------------------------------
      # HIGH-RISK REGIONAL INSTABILITY
      # ---------------------------------------------
      
      if(region == 3){
        
        probs[1] <- probs[1] + 0.10
      }
      
      probs <- probs / sum(probs)
      
      connect_vals <- c(connect_vals, probs)
    }
  }
  
  cpt_connectivity <- array(
    
    connect_vals,
    
    dim = c(3,3,4),
    
    dimnames = list(
      
      Host_Connectivity =
        states$Host_Connectivity,
      
      Region_Context =
        states$Region_Context,
      
      Restoration_Intensity =
        states$Restoration_Intensity
    )
  )
  
  # =======================================================
  # TICK abundance CPT
  # recovery can initially
  # increase tick abundance
  # =======================================================
  
  tick_vals <- c(
    
    # poor connectivity
    0.10,0.25,0.65,
    
    # moderate connectivity
    0.25,0.45,0.30,
    
    # optimal connectivity
    0.40,0.20,0.05
  )
  
  tick_vals <- tick_vals +
    runif(length(tick_vals), -0.05, 0.05)
  
  tick_vals[tick_vals < 0.01] <- 0.01
  
  tick_vals <- unlist(
    lapply(
      split(
        tick_vals,
        ceiling(seq_along(tick_vals)/3)
      ),
      function(x) x/sum(x)
    )
  )
  
  cpt_tick_abundance <- array(
    
    tick_vals,
    
    dim = c(3,3),
    
    dimnames = list(
      
      Tick_Abundance =
        states$Tick_Abundance,
      
      Host_Connectivity =
        states$Host_Connectivity
    )
  )
  
  # =======================================================
  # DISEASE CPT
  # =======================================================
  
  disease_vals <- c(
    
    # low abundance
    0.70,0.05,
    
    # moderate abundance
    0.50,0.40,
    
    # high abundance
    0.20,0.80
  )
  
  disease_vals <- disease_vals +
    runif(length(disease_vals), -0.06, 0.06)
  
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
    
    dim = c(2,3),
    
    dimnames = list(
      
      Tick_Borne_Diseases =
        states$Tick_Borne_Diseases,
      
      Tick_Abundance =
        states$Tick_Abundance
    )
  )
  
  # =======================================================
  # TEMPORAL TRAJECTORIES
  # =======================================================
  
  time_vals <- c()
  
  for(region in 1:3){
    
    for(restoration in 1:4){
      
      for(disease in 1:2){
        
        # =================================================
        # HIGH RESTORATION + LOW DISEASE
        # =================================================
        
        if(restoration == 4 & disease == 1){
          
          probs <- c(
            
            sample_prob(0.45,0.65),
            sample_prob(0.20,0.30),
            sample_prob(0.05,0.15),
            sample_prob(0.01,0.08)
          )
          
          # =================================================
          # NO RESTORATION + HIGH DISEASE
          # =================================================
          
        } else if(restoration == 1 & disease == 2){
          
          probs <- c(
            
            sample_prob(0.01,0.08),
            sample_prob(0.05,0.15),
            sample_prob(0.20,0.35),
            sample_prob(0.45,0.70)
          )
          
          # =================================================
          # MODERATE TRANSITION STATES
          # =================================================
          
        } else {
          
          probs <- c(
            
            sample_prob(0.15,0.35),
            sample_prob(0.20,0.35),
            sample_prob(0.20,0.35),
            sample_prob(0.10,0.30)
          )
        }
        
        # =================================================
        # HIGH-RISK REGIONAL DELAY
        # =================================================
        
        if(region == 3){
          
          probs[4] <- probs[4] + 0.12
          probs[1] <- probs[1] * 0.7
        }
        
        probs <- probs / sum(probs)
        
        time_vals <- c(time_vals, probs)
      }
    }
  }
  
  cpt_time <- array(
    
    time_vals,
    
    dim = c(4,3,4,2),
    
    dimnames = list(
      
      Time_to_Stable_State =
        states$Time_to_Stable_State,
      
      Region_Context =
        states$Region_Context,
      
      Restoration_Intensity =
        states$Restoration_Intensity,
      
      Tick_Borne_Diseases =
        states$Tick_Borne_Diseases
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
      
      Restoration_Intensity =
        cpt_restoration,
      
      Host_Connectivity =
        cpt_connectivity,
      
      Tick_Abundance =
        cpt_tick_abundance,
      
      Tick_Borne_Diseases =
        cpt_disease,
      
      Time_to_Stable_State =
        cpt_time
    )
  )
  
  bn_grain <- as.grain(fitted_bn)
  
  # =======================================================
  # QUERY EACH RESTORATION LEVEL
  # =======================================================
  
  for(restoration_level in states$Restoration_Intensity){
    
    bn_temp <- setEvidence(
      
      bn_grain,
      
      nodes = c(
        "Region_Context",
        "Restoration_Intensity"
      ),
      
      states = c(
        "High_Risk_Context",
        restoration_level
      )
    )
    
    # ---------------------------------------------
    # DISEASE QUERY
    # ---------------------------------------------
    
    disease_q <- querygrain(
      bn_temp,
      nodes = "Tick_Borne_Diseases"
    )$Tick_Borne_Diseases
    
    disease_risk <- disease_q["present"]
    
    # ---------------------------------------------
    # STABILIZATION QUERY
    # ---------------------------------------------
    
    stab_q <- querygrain(
      bn_temp,
      nodes = "Time_to_Stable_State"
    )$Time_to_Stable_State
    
    rapid_stabilization <-
      stab_q["1_3_years"] +
      stab_q["4_7_years"]
    
    delayed_stabilization <-
      stab_q["8_12_years"] +
      stab_q["13_15_years"]
    
    # ---------------------------------------------
    # STORE
    # ---------------------------------------------
    
    temp_df <- data.frame(
      
      iteration = i,
      
      restoration = restoration_level,
      
      disease_risk = disease_risk,
      
      rapid_stabilization = rapid_stabilization,
      
      delayed_stabilization = delayed_stabilization
    )
    
    results <- bind_rows(
      results,
      temp_df
    )
  }
}
# =========================================================
# SUMMARISE UNCERTAINTY
# =========================================================

summary_df <- results %>%
  
  group_by(restoration) %>%
  
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
results <- results %>%
  
  filter(
    !is.na(disease_risk),
    !is.na(rapid_stabilization)
  )
tradeoff_summary <- results %>%
  
  group_by(restoration) %>%
  
  summarise(
    
    mean_disease =
      mean(disease_risk),
    
    lower_disease =
      quantile(disease_risk,0.025),
    
    upper_disease =
      quantile(disease_risk,0.975),
    
    mean_stability =
      mean(rapid_stabilization),
    
    lower_stability =
      quantile(rapid_stabilization,0.025),
    
    upper_stability =
      quantile(rapid_stabilization,0.975)
  )

print(tradeoff_summary)

ggplot(
  
  tradeoff_summary,
  
  aes(
    x = mean_disease,
    y = mean_stability,
    color = restoration
  )
  
) +
  
  geom_point(size = 5) +
  
  geom_line(aes(group = 1),
            linewidth = 1.2) +
  
  # horizontal uncertainty
  geom_errorbar(
    
    aes(
      ymin = lower_stability,
      ymax = upper_stability
    ),
    
    width = 0
  ) +
  
  # vertical uncertainty
  geom_errorbarh(
    
    aes(
      xmin = lower_disease,
      xmax = upper_disease
    ),
    
    height = 0
  ) +
  
  labs(
    
    title =
      "Restoration–Disease Trade-off Curve",
    
    subtitle =
      "Tick-borne disease risk versus ecosystem stabilization",
    
    x =
      "Probability of tick-borne disease emergence",
    
    y =
      "Probability of rapid stabilization (1–7 years)"
  ) +
  
  theme_minimal(base_size = 14)

# =========================================================
# LONG FORMAT
# =========================================================

plot_df <- results %>%
  
  pivot_longer(
    
    cols = starts_with("years"),
    
    names_to = "trajectory",
    
    values_to = "probability"
  )

plot_df$trajectory <- factor(
  
  plot_df$trajectory,
  
  levels = c(
    "years_1_3",
    "years_4_7",
    "years_8_12",
    "years_13_15"
  )
)

# =========================================================
# UNCERTAINTY CORRIDOR PLOT
# =========================================================

ggplot(
  
  plot_df,
  
  aes(
    x = trajectory,
    y = probability,
    color = restoration,
    fill = restoration,
    group = restoration
  )
  
) +
  
  stat_summary(
    fun = mean,
    geom = "line",
    linewidth = 1.3
  ) +
  
  stat_summary(
    fun.data = mean_cl_normal,
    geom = "ribbon",
    alpha = 0.18,
    color = NA
  ) +
  
  labs(
    
    title =
      "Tick Disease Stabilization Trajectories",
    
    subtitle =
      "Restoration intensity gradient with Monte Carlo uncertainty",
    
    x = "Time to stabilization",
    
    y = "Probability"
  ) +
  
  theme_minimal(base_size = 14)
