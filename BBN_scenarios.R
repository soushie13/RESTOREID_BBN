# =========================================================
# LOAD LIBRARIES
# =========================================================

library(bnlearn)
library(gRain)
library(gRbase)
library(dplyr)
library(ggplot2)
library(tidyr)

# =====================================================
# LOAD PACKAGES
# =====================================================

library(bnlearn)
library(gRain)
library(gRbase)

# =====================================================
# DEFINE NODE STATES
# =====================================================

node_states <- list()

# Region
node_states[["Region_Context"]] <- c(
  "Low_Risk_Context",
  "Moderate_Risk_Context",
  "High_Risk_Context"
)

# Disease nodes
node_states[["Rodent_Borne_Diseases"]] <- c(
  "absent",
  "present"
)

# Mitigation nodes
node_states[["Habitat_Design"]] <- c(
  "not_implemented",
  "implemented"
)

node_states[["Vector_Management_Strategies"]] <- c(
  "not_implemented",
  "implemented"
)

node_states[["Public_Health_Awareness"]] <- c(
  "not_implemented",
  "implemented"
)

# Temporal node
node_states[["Time_to_Stable_State"]] <- c(
  "1_3_years",
  "4_7_years",
  "8_12_years",
  "13_15_years"
)

# =====================================================
# CREATE NETWORK STRUCTURE
# =====================================================

all_nodes <- c(
  "Region_Context",
  "Rodent_Borne_Diseases",
  "Habitat_Design",
  "Vector_Management_Strategies",
  "Public_Health_Awareness",
  "Time_to_Stable_State"
)

bn <- empty.graph(all_nodes)

# =====================================================
# DEFINE ARCS
# =====================================================

arcs_df <- data.frame(
  
  from = c(
    "Region_Context",
    "Rodent_Borne_Diseases",
    "Rodent_Borne_Diseases",
    "Region_Context",
    "Habitat_Design",
    "Vector_Management_Strategies"
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

# Add arcs
arcs(bn) <- as.matrix(arcs_df)

# =====================================================
# REGION PRIOR CPT
# =====================================================

cpt_region <- cptable(
  ~Region_Context,
  
  values = c(
    0.33,
    0.34,
    0.33
  ),
  
  levels = node_states[["Region_Context"]]
)

# =====================================================
# RODENT DISEASE CPT
# =====================================================

cpt_rodent_best <- cptable(
  
  ~Rodent_Borne_Diseases | Region_Context,
  
  values = c(
    
    # LOW RISK
    0.80, 0.20,
    
    # MODERATE RISK
    0.60, 0.40,
    
    # HIGH RISK
    0.40, 0.60
  ),
  
  levels = node_states[["Rodent_Borne_Diseases"]]
)

# =====================================================
# HABITAT DESIGN CPT
# =====================================================

cpt_habitat <- cptable(
  
  ~Habitat_Design | Rodent_Borne_Diseases,
  
  values = c(
    
    # disease absent
    0.70, 0.30,
    
    # disease present
    0.20, 0.80
  ),
  
  levels = node_states[["Habitat_Design"]]
)

# =====================================================
# PUBLIC HEALTH AWARENESS CPT
# =====================================================

cpt_awareness <- cptable(
  
  ~Public_Health_Awareness | Rodent_Borne_Diseases,
  
  values = c(
    
    # disease absent
    0.75, 0.25,
    
    # disease present
    0.20, 0.80
  ),
  
  levels = node_states[["Public_Health_Awareness"]]
)

# =====================================================
# VECTOR MANAGEMENT PRIOR
# =====================================================

cpt_vector <- cptable(
  
  ~Vector_Management_Strategies,
  
  values = c(
    0.50,
    0.50
  ),
  
  levels = node_states[["Vector_Management_Strategies"]]
)

# =====================================================
# TEMPORAL CPT
# =====================================================

cpt_time_best <- cptable(
  
  ~Time_to_Stable_State |
    Region_Context :
    Habitat_Design :
    Vector_Management_Strategies,
  
  values = c(
    
    # =================================================
    # LOW RISK REGION
    # =================================================
    
    # Habitat implemented + Vector implemented
    0.55, 0.30, 0.10, 0.05,
    
    # Habitat implemented + Vector NOT implemented
    0.40, 0.30, 0.20, 0.10,
    
    # Habitat NOT implemented + Vector implemented
    0.35, 0.30, 0.22, 0.13,
    
    # Habitat NOT implemented + Vector NOT implemented
    0.22, 0.28, 0.28, 0.22,
    
    # =================================================
    # MODERATE RISK REGION
    # =================================================
    
    0.30, 0.35, 0.25, 0.10,
    0.20, 0.30, 0.30, 0.20,
    0.18, 0.27, 0.30, 0.25,
    0.05, 0.15, 0.35, 0.45,
    
    # =================================================
    # HIGH RISK REGION
    # =================================================
    
    0.10, 0.25, 0.35, 0.30,
    0.06, 0.14, 0.30, 0.50,
    0.05, 0.12, 0.28, 0.55,
    0.02, 0.08, 0.30, 0.60
  ),
  
  levels = node_states[["Time_to_Stable_State"]]
)
# =====================================================
# WORST-CASE RODENT DISEASE CPT
# =====================================================

cpt_rodent_worst <- cptable(
  
  ~Rodent_Borne_Diseases | Region_Context,
  
  values = c(
    
    # LOW RISK
    0.60, 0.40,
    
    # MODERATE RISK
    0.40, 0.60,
    
    # HIGH RISK
    0.20, 0.80
  ),
  
  levels = node_states[["Rodent_Borne_Diseases"]]
)

# =====================================================
# WORST-CASE TEMPORAL CPT
# =====================================================

cpt_time_worst <- cptable(
  
  ~Time_to_Stable_State |
    Region_Context :
    Habitat_Design :
    Vector_Management_Strategies,
  
  values = c(
    
    # =================================================
    # LOW RISK REGION
    # =================================================
    
    0.30, 0.30, 0.25, 0.15,
    0.20, 0.25, 0.30, 0.25,
    0.15, 0.25, 0.30, 0.30,
    0.10, 0.20, 0.30, 0.40,
    
    # =================================================
    # MODERATE RISK REGION
    # =================================================
    
    0.15, 0.20, 0.35, 0.30,
    0.10, 0.15, 0.35, 0.40,
    0.08, 0.12, 0.35, 0.45,
    0.05, 0.10, 0.30, 0.55,
    
    # =================================================
    # HIGH RISK REGION
    # =================================================
    
    0.05, 0.10, 0.35, 0.50,
    0.03, 0.07, 0.30, 0.60,
    0.02, 0.06, 0.27, 0.65,
    0.01, 0.04, 0.25, 0.70
  ),
  
  levels = node_states[["Time_to_Stable_State"]]
)
# =====================================================
# COMPILE CPTS
# =====================================================

cpts_best <- compileCPT(list(
  cpt_region,
  cpt_rodent_best,
  cpt_habitat,
  cpt_awareness,
  cpt_vector,
  cpt_time_best
))
cpts_worst <- compileCPT(list(
  
  cpt_region,
  cpt_rodent_worst,
  cpt_habitat,
  cpt_awareness,
  cpt_vector,
  cpt_time_worst
))
# =====================================================
# BUILD NETWORK
# =====================================================

bn_best <- grain(cpts_best)
bn_worst <- grain(cpts_worst)

# =====================================================
# SET EVIDENCE
# =====================================================

bn_highrisk <- setEvidence(
  bn_best,
  nodes = "Region_Context",
  states = "High_Risk_Context"
)
bn_worst_highrisk <- setEvidence(
  bn_worst,
  nodes = "Region_Context",
  states = "High_Risk_Context"
)

# =====================================================
# QUERY TEMPORAL DISTRIBUTION
# =====================================================

querygrain(
  bn_highrisk,
  nodes = "Time_to_Stable_State"
)

querygrain(
  bn_best,
  nodes = "Time_to_Stable_State"
)

querygrain(
  bn_worst,
  nodes = "Time_to_Stable_State"
)

bn_low <- setEvidence(
  bn_best,
  nodes = "Region_Context",
  states = "Low_Risk_Context"
)

querygrain(
  bn_low,
  nodes = "Time_to_Stable_State"
)

bn_moderate <- setEvidence(
  bn_best,
  nodes = "Region_Context",
  states = "Moderate_Risk_Context"
)

querygrain(
  bn_moderate,
  nodes = "Time_to_Stable_State"
)

bn_habitat <- setEvidence(
  bn_best,
  nodes = "Habitat_Design",
  states = "implemented"
)

querygrain(
  bn_habitat,
  nodes = "Time_to_Stable_State"
)

bn_no_habitat <- setEvidence(
  bn_best,
  nodes = "Habitat_Design",
  states = "not_implemented"
)

querygrain(
  bn_no_habitat,
  nodes = "Time_to_Stable_State"
)

bn_management <- setEvidence(
  bn_best,
  nodes = c(
    "Habitat_Design",
    "Vector_Management_Strategies"
  ),
  states = c(
    "implemented",
    "implemented"
  )
)

querygrain(
  bn_management,
  nodes = "Time_to_Stable_State"
)

bn_region <- setEvidence(
  bn_best,
  nodes = "Region_Context",
  states = "High_Risk_Context"
)

querygrain(
  bn_region,
  nodes = "Rodent_Borne_Diseases"
)

bn_sensitivity <- setEvidence(
  bn_best,
  nodes = c(
    "Region_Context",
    "Habitat_Design"
  ),
  states = c(
    "High_Risk_Context",
    "implemented"
  )
)

querygrain(
  bn_sensitivity,
  nodes = "Time_to_Stable_State"
)

result <- querygrain(
  bn_highrisk,
  nodes = "Time_to_Stable_State"
)$Time_to_Stable_State

long_term_risk <- result["8_12_years"] +
  result["13_15_years"]

print(long_term_risk)

# =====================================================
# DEFINE NODE STATES for AVIAN TICK MOSQUITO WATER BORNE 
# =====================================================

node_states[["Tick_Borne_Diseases"]] <- c(
  "absent",
  "present"
)

node_states[["Waterborne_Diseases"]] <- c(
  "absent",
  "present"
)

node_states[["Mosquito_Borne_Diseases"]] <- c(
  "absent",
  "present"
)

node_states[["Avian_Diseases"]] <- c(
  "absent",
  "present"
)

node_states[["Connectivity_Existing_Forests"]] <- c(
  "not_implemented",
  "implemented"
)

node_states[["Human_Activity_Restrictions"]] <- c(
  "not_implemented",
  "implemented"
)

node_states[["Wetland_Design"]] <- c(
  "not_implemented",
  "implemented"
)

node_states[["Vector_Management"]] <- c(
  "not_implemented",
  "implemented"
)

node_states[["Culling_of_Birds"]] <- c(
  "not_implemented",
  "implemented"
)

# =====================================================
# BEST-CASE CPTs
# =====================================================

# ---------------------------------
# Tick-Borne Diseases → Habitat Design
# ---------------------------------

cpt_tick_habitat_best <- cptable(
  
  ~Habitat_Design | Tick_Borne_Diseases,
  
  values = c(
    
    # disease absent
    0.70, 0.30,
    
    # disease present
    0.20, 0.80
  ),
  
  levels = node_states[["Habitat_Design"]]
)

# ---------------------------------
# Tick-Borne Diseases → Connectivity
# ---------------------------------

cpt_tick_connectivity_best <- cptable(
  
  ~Connectivity_Existing_Forests | Tick_Borne_Diseases,
  
  values = c(
    
    # disease absent
    0.60, 0.40,
    
    # disease present
    0.15, 0.85
  ),
  
  levels = node_states[["Connectivity_Existing_Forests"]]
)

# ---------------------------------
# Waterborne Diseases → Awareness
# ---------------------------------

cpt_water_awareness_best <- cptable(
  
  ~Public_Health_Awareness | Waterborne_Diseases,
  
  values = c(
    
    # disease absent
    0.75, 0.25,
    
    # disease present
    0.20, 0.80
  ),
  
  levels = node_states[["Public_Health_Awareness"]]
)

# ---------------------------------
# Waterborne Diseases → Restrictions
# ---------------------------------

cpt_water_restrictions_best <- cptable(
  
  ~Human_Activity_Restrictions | Waterborne_Diseases,
  
  values = c(
    
    # disease absent
    0.70, 0.30,
    
    # disease present
    0.10, 0.90
  ),
  
  levels = node_states[["Human_Activity_Restrictions"]]
)

# ---------------------------------
# Mosquito-Borne Diseases → Wetland Design
# ---------------------------------

cpt_mosquito_wetland_best <- cptable(
  
  ~Wetland_Design | Mosquito_Borne_Diseases,
  
  values = c(
    
    # disease absent
    0.65, 0.35,
    
    # disease present
    0.15, 0.85
  ),
  
  levels = node_states[["Wetland_Design"]]
)

# ---------------------------------
# Mosquito-Borne Diseases → Vector Management
# ---------------------------------

cpt_mosquito_vector_best <- cptable(
  
  ~Vector_Management | Mosquito_Borne_Diseases,
  
  values = c(
    
    # disease absent
    0.60, 0.40,
    
    # disease present
    0.10, 0.90
  ),
  
  levels = node_states[["Vector_Management"]]
)

# ---------------------------------
# Avian Diseases → Public Awareness
# ---------------------------------

cpt_avian_awareness_best <- cptable(
  
  ~Public_Health_Awareness | Avian_Diseases,
  
  values = c(
    
    # disease absent
    0.75, 0.25,
    
    # disease present
    0.20, 0.80
  ),
  
  levels = node_states[["Public_Health_Awareness"]]
)

# ---------------------------------
# Avian Diseases → Culling
# ---------------------------------

cpt_avian_culling_best <- cptable(
  
  ~Culling_of_Birds | Avian_Diseases,
  
  values = c(
    
    # disease absent
    0.80, 0.20,
    
    # disease present
    0.10, 0.90
  ),
  
  levels = node_states[["Culling_of_Birds"]]
)

# =====================================================
# WORST-CASE CPTs
# =====================================================

# Tick-Borne Diseases → Habitat Design

cpt_tick_habitat_worst <- cptable(
  
  ~Habitat_Design | Tick_Borne_Diseases,
  
  values = c(
    
    0.80, 0.20,
    0.45, 0.55
  ),
  
  levels = node_states[["Habitat_Design"]]
)

# Tick-Borne Diseases → Connectivity

cpt_tick_connectivity_worst <- cptable(
  
  ~Connectivity_Existing_Forests | Tick_Borne_Diseases,
  
  values = c(
    
    0.75, 0.25,
    0.40, 0.60
  ),
  
  levels = node_states[["Connectivity_Existing_Forests"]]
)

# Waterborne Diseases → Awareness

cpt_water_awareness_worst <- cptable(
  
  ~Public_Health_Awareness | Waterborne_Diseases,
  
  values = c(
    
    0.80, 0.20,
    0.50, 0.50
  ),
  
  levels = node_states[["Public_Health_Awareness"]]
)

# Waterborne Diseases → Restrictions

cpt_water_restrictions_worst <- cptable(
  
  ~Human_Activity_Restrictions | Waterborne_Diseases,
  
  values = c(
    
    0.80, 0.20,
    0.40, 0.60
  ),
  
  levels = node_states[["Human_Activity_Restrictions"]]
)

# Mosquito-Borne Diseases → Wetland Design

cpt_mosquito_wetland_worst <- cptable(
  
  ~Wetland_Design | Mosquito_Borne_Diseases,
  
  values = c(
    
    0.75, 0.25,
    0.40, 0.60
  ),
  
  levels = node_states[["Wetland_Design"]]
)

# Mosquito-Borne Diseases → Vector Management

cpt_mosquito_vector_worst <- cptable(
  
  ~Vector_Management | Mosquito_Borne_Diseases,
  
  values = c(
    
    0.75, 0.25,
    0.35, 0.65
  ),
  
  levels = node_states[["Vector_Management"]]
)

# Avian Diseases → Awareness

cpt_avian_awareness_worst <- cptable(
  
  ~Public_Health_Awareness | Avian_Diseases,
  
  values = c(
    
    0.80, 0.20,
    0.45, 0.55
  ),
  
  levels = node_states[["Public_Health_Awareness"]]
)

# Avian Diseases → Culling

cpt_avian_culling_worst <- cptable(
  
  ~Culling_of_Birds | Avian_Diseases,
  
  values = c(
    
    0.85, 0.15,
    0.45, 0.55
  ),
  
  levels = node_states[["Culling_of_Birds"]]
)
# =====================================================
# TEMPORAL + REGION-SPECIFIC CPTs
# FOR ALL DISEASE TYPES
# =====================================================

# =====================================================
# DEFINE STATES
# =====================================================

node_states <- list(
  
  Region_Context = c(
    "Low_Risk_Context",
    "Moderate_Risk_Context",
    "High_Risk_Context"
  ),
  
  Time_to_Stable_State = c(
    "1_3_years",
    "4_7_years",
    "8_12_years",
    "13_15_years"
  ),
  
  Habitat_Design = c(
    "not_implemented",
    "implemented"
  ),
  
  Vector_Management = c(
    "not_implemented",
    "implemented"
  ),
  
  Wetland_Design = c(
    "not_implemented",
    "implemented"
  ),
  
  Connectivity_Existing_Forests = c(
    "not_implemented",
    "implemented"
  ),
  
  Public_Health_Awareness = c(
    "not_implemented",
    "implemented"
  ),
  
  Human_Activity_Restrictions = c(
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
  )
)

# =====================================================
# REGION PRIOR
# =====================================================

cpt_region <- cptable(
  ~Region_Context,
  values = c(0.33, 0.34, 0.33),
  levels = node_states$Region_Context
)

# =====================================================
# DISEASE PRIORS
# =====================================================

# -----------------------------
# BEST CASE
# -----------------------------

cpt_rodent_best <- cptable(
  ~Rodent_Borne_Diseases | Region_Context,
  values = c(
    0.80,0.20,
    0.60,0.40,
    0.40,0.60
  ),
  levels = node_states$Rodent_Borne_Diseases
)

cpt_tick_best <- cptable(
  ~Tick_Borne_Diseases | Region_Context,
  values = c(
    0.75,0.25,
    0.55,0.45,
    0.35,0.65
  ),
  levels = node_states$Tick_Borne_Diseases
)

cpt_water_best <- cptable(
  ~Waterborne_Diseases | Region_Context,
  values = c(
    0.85,0.15,
    0.65,0.35,
    0.45,0.55
  ),
  levels = node_states$Waterborne_Diseases
)

cpt_mosquito_best <- cptable(
  ~Mosquito_Borne_Diseases | Region_Context,
  values = c(
    0.70,0.30,
    0.50,0.50,
    0.30,0.70
  ),
  levels = node_states$Mosquito_Borne_Diseases
)

cpt_avian_best <- cptable(
  ~Avian_Diseases | Region_Context,
  values = c(
    0.80,0.20,
    0.60,0.40,
    0.40,0.60
  ),
  levels = node_states$Avian_Diseases
)

# -----------------------------
# WORST CASE
# -----------------------------

cpt_rodent_worst <- cptable(
  ~Rodent_Borne_Diseases | Region_Context,
  values = c(
    0.60,0.40,
    0.40,0.60,
    0.20,0.80
  ),
  levels = node_states$Rodent_Borne_Diseases
)

cpt_tick_worst <- cptable(
  ~Tick_Borne_Diseases | Region_Context,
  values = c(
    0.55,0.45,
    0.35,0.65,
    0.15,0.85
  ),
  levels = node_states$Tick_Borne_Diseases
)

cpt_water_worst <- cptable(
  ~Waterborne_Diseases | Region_Context,
  values = c(
    0.65,0.35,
    0.45,0.55,
    0.20,0.80
  ),
  levels = node_states$Waterborne_Diseases
)

cpt_mosquito_worst <- cptable(
  ~Mosquito_Borne_Diseases | Region_Context,
  values = c(
    0.50,0.50,
    0.30,0.70,
    0.10,0.90
  ),
  levels = node_states$Mosquito_Borne_Diseases
)

cpt_avian_worst <- cptable(
  ~Avian_Diseases | Region_Context,
  values = c(
    0.60,0.40,
    0.40,0.60,
    0.20,0.80
  ),
  levels = node_states$Avian_Diseases
)

# =====================================================
# TEMPORAL CPT
# =====================================================

cpt_time_best <- cptable(
  
  ~Time_to_Stable_State |
    Region_Context :
    Habitat_Design :
    Vector_Management,
  
  values = c(
    
    # LOW RISK
    0.55,0.30,0.10,0.05,
    0.40,0.30,0.20,0.10,
    0.35,0.30,0.22,0.13,
    0.22,0.28,0.28,0.22,
    
    # MODERATE RISK
    0.30,0.35,0.25,0.10,
    0.20,0.30,0.30,0.20,
    0.18,0.27,0.30,0.25,
    0.05,0.15,0.35,0.45,
    
    # HIGH RISK
    0.10,0.25,0.35,0.30,
    0.06,0.14,0.30,0.50,
    0.05,0.12,0.28,0.55,
    0.02,0.08,0.30,0.60
  ),
  
  levels = node_states$Time_to_Stable_State
)

# =====================================================
# WORST CASE TEMPORAL CPT
# =====================================================

cpt_time_worst <- cptable(
  
  ~Time_to_Stable_State |
    Region_Context :
    Habitat_Design :
    Vector_Management,
  
  values = c(
    
    # LOW RISK
    0.30,0.30,0.25,0.15,
    0.20,0.25,0.30,0.25,
    0.15,0.25,0.30,0.30,
    0.10,0.20,0.30,0.40,
    
    # MODERATE RISK
    0.15,0.20,0.35,0.30,
    0.10,0.15,0.35,0.40,
    0.08,0.12,0.35,0.45,
    0.05,0.10,0.30,0.55,
    
    # HIGH RISK
    0.05,0.10,0.35,0.50,
    0.03,0.07,0.30,0.60,
    0.02,0.06,0.27,0.65,
    0.01,0.04,0.25,0.70
  ),
  
  levels = node_states$Time_to_Stable_State
)

# =====================================================
# COMPILE NETWORKS
# =====================================================

cpts_best <- compileCPT(list(
  cpt_region,
  cpt_rodent_best,
  cpt_tick_best,
  cpt_water_best,
  cpt_mosquito_best,
  cpt_avian_best,
  cpt_time_best
))

cpts_worst <- compileCPT(list(
  cpt_region,
  cpt_rodent_worst,
  cpt_tick_worst,
  cpt_water_worst,
  cpt_mosquito_worst,
  cpt_avian_worst,
  cpt_time_worst
))

bn_best <- grain(cpts_best)
bn_worst <- grain(cpts_worst)

# =====================================================
# QUERY FUNCTIONS
# =====================================================

# -----------------------------
# HIGH-RISK BEST CASE
# -----------------------------

bn_best_high <- setEvidence(
  bn_best,
  nodes = "Region_Context",
  states = "High_Risk_Context"
)

querygrain(
  bn_best_high,
  nodes = c(
    "Rodent_Borne_Diseases",
    "Tick_Borne_Diseases",
    "Waterborne_Diseases",
    "Mosquito_Borne_Diseases",
    "Avian_Diseases",
    "Time_to_Stable_State"
  )
)

# -----------------------------
# HIGH-RISK WORST CASE
# -----------------------------

bn_worst_high <- setEvidence(
  bn_worst,
  nodes = "Region_Context",
  states = "High_Risk_Context"
)

querygrain(
  bn_worst_high,
  nodes = c(
    "Rodent_Borne_Diseases",
    "Tick_Borne_Diseases",
    "Waterborne_Diseases",
    "Mosquito_Borne_Diseases",
    "Avian_Diseases",
    "Time_to_Stable_State"
  )
)

# =====================================================
# LOW-RISK REGION
# =====================================================

bn_low <- setEvidence(
  bn_best,
  nodes = "Region_Context",
  states = "Low_Risk_Context"
)

querygrain(
  bn_low,
  nodes = c(
    "Time_to_Stable_State",
    "Rodent_Borne_Diseases",
    "Tick_Borne_Diseases",
    "Waterborne_Diseases"
  )
)

# =====================================================
# MODERATE-RISK REGION
# =====================================================

bn_moderate <- setEvidence(
  bn_best,
  nodes = "Region_Context",
  states = "Moderate_Risk_Context"
)

querygrain(
  bn_moderate,
  nodes = c(
    "Time_to_Stable_State",
    "Mosquito_Borne_Diseases",
    "Avian_Diseases"
  )
)

# =====================================================
# MANAGEMENT INTERVENTION SCENARIO
# =====================================================

bn_management <- setEvidence(
  
  bn_best,
  
  nodes = c(
    "Habitat_Design",
    "Vector_Management"
  ),
  
  states = c(
    "implemented",
    "implemented"
  )
)

querygrain(
  bn_management,
  nodes = "Time_to_Stable_State"
)

# =====================================================
# SENSITIVITY ANALYSIS
# =====================================================

bn_sensitivity <- setEvidence(
  
  bn_best,
  
  nodes = c(
    "Region_Context",
    "Habitat_Design"
  ),
  
  states = c(
    "High_Risk_Context",
    "implemented"
  )
)

querygrain(
  bn_sensitivity,
  nodes = c(
    "Time_to_Stable_State",
    "Rodent_Borne_Diseases",
    "Tick_Borne_Diseases"
  )
)
