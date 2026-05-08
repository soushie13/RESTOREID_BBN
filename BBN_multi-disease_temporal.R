# ============================================================
# BAYESIAN BELIEF NETWORK (BBN)
# MULTI-DISEASE + TEMPORAL SCENARIOS
# ============================================================

# ============================================================
# 1. LOAD PACKAGES
# ============================================================

install.packages(c(
  "bnlearn",
  "gRain",
  "gRbase",
  "dplyr",
  "tidyr",
  "readr"
))

library(bnlearn)
library(gRain)
library(gRbase)
library(dplyr)
library(tidyr)
library(readr)

# ============================================================
# 2. DEFINE NODE STATES
# ============================================================

node_states <- list(
  
  Region_Context = c(
    "Low_Risk_Context",
    "Moderate_Risk_Context",
    "High_Risk_Context"
  ),
  
  Habitat_Design = c(
    "not_implemented",
    "implemented"
  ),
  
  Vector_Management_Strategies = c(
    "not_implemented",
    "implemented"
  ),
  
  Public_Health_Awareness = c(
    "not_implemented",
    "implemented"
  ),
  
  Connectivity_Existing_Forests = c(
    "not_connected",
    "connected"
  ),
  
  Human_Activity_Restrictions = c(
    "not_implemented",
    "implemented"
  ),
  
  Wetland_Design = c(
    "not_implemented",
    "implemented"
  ),
  
  Culling_of_Birds = c(
    "not_implemented",
    "implemented"
  ),
  
  Rodent_Borne_Diseases = c(
    "absent",
    "present"
  ),
  
  Tick_Borne_Diseases = c(
    "absent",
    "present"
  ),
  
  Waterborne_Diseases = c(
    "absent",
    "present"
  ),
  
  Mosquito_Borne_Diseases = c(
    "absent",
    "present"
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

# ============================================================
# 3. DEFINE ALL NODES
# ============================================================

all_nodes <- names(node_states)

# ============================================================
# 4. DEFINE NETWORK STRUCTURE (ARCS)
# ============================================================

arcs_df <- data.frame(

  from = c(

    # Restoration
    "Wetland_Restoration",
    "Grassland_Restoration",
    "Urban_Green_Spaces",
    "Afforestation",

    # Ecological interactions
    "Host_Vector_Abundance",
    "Rodent_Abundance",
    "Tick_Population",
    "Human_Presence",
    "Wildlife_Migration",

    # Disease → management
    "Rodent_Borne_Diseases",
    "Tick_Borne_Diseases",
    "Waterborne_Diseases",
    "Mosquito_Borne_Diseases",
    "Avian_Diseases"

  ),

  to = c(

    # Restoration impacts
    "Host_Vector_Abundance",
    "Host_Vector_Abundance",
    "Rodent_Abundance",
    "Tick_Population",

    # Ecological → disease
    "Mosquito_Borne_Diseases",
    "Rodent_Borne_Diseases",
    "Tick_Borne_Diseases",
    "Waterborne_Diseases",
    "Avian_Diseases",

    # Disease → management
    "Habitat_Design",
    "Connectivity_Existing_Forests",
    "Public_Health_Awareness",
    "Vector_Management_Strategies",
    "Culling_of_Birds"
  )
)

# ============================================================
# 5. CREATE NETWORK STRUCTURE
# ============================================================

bn_structure <- empty.graph(all_nodes)

arcs(bn_structure) <- as.matrix(arcs_df)

# ============================================================
# 6. DEFINE CPTS
# ============================================================

cpts_best <- list()
cpts_worst <- list()

# ============================================================
# 7. REGION PRIOR
# ============================================================

cpts_best[["Region_Context"]] <- cptable(
  ~Region_Context,
  values = c(0.33, 0.34, 0.33),
  levels = node_states$Region_Context
)

cpts_worst[["Region_Context"]] <- cptable(
  ~Region_Context,
  values = c(0.33, 0.34, 0.33),
  levels = node_states$Region_Context
)

# ============================================================
# 8. DISEASE CPTS
# ============================================================

disease_nodes <- c(
  "Rodent_Borne_Diseases",
  "Tick_Borne_Diseases",
  "Waterborne_Diseases",
  "Mosquito_Borne_Diseases",
  "Avian_Diseases"
)

for(disease in disease_nodes){
  
  # BEST CASE
  
  cpts_best[[disease]] <- cptable(
    
    as.formula(
      paste0("~", disease, "| Region_Context")
    ),
    
    values = c(
      
      0.80, 0.20,
      0.60, 0.40,
      0.40, 0.60
      
    ),
    
    levels = node_states[[disease]]
  )
  
  # WORST CASE
  
  cpts_worst[[disease]] <- cptable(
    
    as.formula(
      paste0("~", disease, "| Region_Context")
    ),
    
    values = c(
      
      0.60, 0.40,
      0.40, 0.60,
      0.20, 0.80
      
    ),
    
    levels = node_states[[disease]]
  )
}

# ============================================================
# 9. MANAGEMENT CPTS
# ============================================================

management_nodes <- c(
  "Habitat_Design",
  "Vector_Management_Strategies",
  "Public_Health_Awareness",
  "Connectivity_Existing_Forests",
  "Human_Activity_Restrictions",
  "Wetland_Design",
  "Culling_of_Birds"
)

for(node in management_nodes){
  
  levels_here <- node_states[[node]]
  
  cpts_best[[node]] <- cptable(
    
    as.formula(
      paste0("~", node)
    ),
    
    values = c(0.30, 0.70),
    
    levels = levels_here
  )
  
  cpts_worst[[node]] <- cptable(
    
    as.formula(
      paste0("~", node)
    ),
    
    values = c(0.70, 0.30),
    
    levels = levels_here
  )
}

# ============================================================
# 10. TEMPORAL CPTS
# ============================================================

cpts_best[["Time_to_Stable_State"]] <- cptable(
  
  ~ Time_to_Stable_State |
    Region_Context,
  
  values = c(
    
    # LOW RISK
    0.55, 0.30, 0.10, 0.05,
    
    # MODERATE RISK
    0.30, 0.35, 0.25, 0.10,
    
    # HIGH RISK
    0.10, 0.25, 0.35, 0.30
    
  ),
  
  levels = node_states$Time_to_Stable_State
)

cpts_worst[["Time_to_Stable_State"]] <- cptable(
  
  ~ Time_to_Stable_State |
    Region_Context,
  
  values = c(
    
    # LOW RISK
    0.30, 0.30, 0.25, 0.15,
    
    # MODERATE RISK
    0.15, 0.25, 0.35, 0.25,
    
    # HIGH RISK
    0.05, 0.10, 0.35, 0.50
    
  ),
  
  levels = node_states$Time_to_Stable_State
)

# ============================================================
# 11. COMPILE NETWORKS
# ============================================================

plist_best <- compileCPT(cpts_best)

plist_worst <- compileCPT(cpts_worst)

bn_best <- grain(plist_best)

bn_worst <- grain(plist_worst)

# ============================================================
# 12. CREATE REGIONAL SCENARIOS
# ============================================================

bn_low <- setEvidence(
  bn_best,
  nodes = "Region_Context",
  states = "Low_Risk_Context"
)

bn_moderate <- setEvidence(
  bn_best,
  nodes = "Region_Context",
  states = "Moderate_Risk_Context"
)

bn_highrisk <- setEvidence(
  bn_best,
  nodes = "Region_Context",
  states = "High_Risk_Context"
)

# ============================================================
# 13. TEMPORAL QUERY FUNCTION
# ============================================================

extract_time_results <- function(network, scenario_name){
  
  result <- querygrain(
    network,
    nodes = "Time_to_Stable_State"
  )$Time_to_Stable_State
  
  data.frame(
    
    Scenario = scenario_name,
    
    Time = names(result),
    
    Probability = as.numeric(result)
  )
}

# ============================================================
# 14. RUN ALL TEMPORAL QUERIES
# ============================================================

time_results <- bind_rows(
  
  extract_time_results(
    bn_best,
    "Best_Case"
  ),
  
  extract_time_results(
    bn_worst,
    "Worst_Case"
  ),
  
  extract_time_results(
    bn_low,
    "Low_Risk_Region"
  ),
  
  extract_time_results(
    bn_moderate,
    "Moderate_Risk_Region"
  ),
  
  extract_time_results(
    bn_highrisk,
    "High_Risk_Region"
  )
)

print(time_results)

# ============================================================
# 15. DISEASE QUERY FUNCTION
# ============================================================

extract_disease_results <- function(
    network,
    scenario_name
){
  
  query <- querygrain(
    network,
    nodes = disease_nodes
  )
  
  data.frame(
    
    Scenario = scenario_name,
    
    Disease = disease_nodes,
    
    Present_Probability = c(
      
      query$Rodent_Borne_Diseases["present"],
      
      query$Tick_Borne_Diseases["present"],
      
      query$Waterborne_Diseases["present"],
      
      query$Mosquito_Borne_Diseases["present"],
      
      query$Avian_Diseases["present"]
    )
  )
}

# ============================================================
# 16. RUN ALL DISEASE QUERIES
# ============================================================

disease_results <- bind_rows(
  
  extract_disease_results(
    bn_best,
    "Best_Case"
  ),
  
  extract_disease_results(
    bn_worst,
    "Worst_Case"
  ),
  
  extract_disease_results(
    bn_low,
    "Low_Risk_Region"
  ),
  
  extract_disease_results(
    bn_moderate,
    "Moderate_Risk_Region"
  ),
  
  extract_disease_results(
    bn_highrisk,
    "High_Risk_Region"
  )
)

print(disease_results)

# ============================================================
# 17. LONG-TERM RISK
# ============================================================

calculate_long_term_risk <- function(
    network,
    scenario_name
){
  
  result <- querygrain(
    network,
    nodes = "Time_to_Stable_State"
  )$Time_to_Stable_State
  
  long_term <- result["8_12_years"] +
    result["13_15_years"]
  
  data.frame(
    
    Scenario = scenario_name,
    
    Long_Term_Risk = as.numeric(long_term)
  )
}

long_term_results <- bind_rows(
  
  calculate_long_term_risk(
    bn_best,
    "Best_Case"
  ),
  
  calculate_long_term_risk(
    bn_worst,
    "Worst_Case"
  ),
  
  calculate_long_term_risk(
    bn_low,
    "Low_Risk_Region"
  ),
  
  calculate_long_term_risk(
    bn_moderate,
    "Moderate_Risk_Region"
  ),
  
  calculate_long_term_risk(
    bn_highrisk,
    "High_Risk_Region"
  )
)

print(long_term_results)

# ============================================================
# 18. EXPORT RESULTS
# ============================================================

write.csv(
  time_results,
  "BBN_Time_Results.csv",
  row.names = FALSE
)

write.csv(
  disease_results,
  "BBN_Disease_Results.csv",
  row.names = FALSE
)

write.csv(
  long_term_results,
  "BBN_Long_Term_Risk.csv",
  row.names = FALSE
)

# ============================================================
# END
# ============================================================