install.packages(c("bnlearn", "gRain", "igraph", "readr"))
library(bnlearn)
library(gRain)
library(igraph)
library(readr)
##library(Rgraphviz)## unavailable
install.packages("igraph")
library(igraph)
library(dplyr)

# Define the nodes of the Bayesian Belief Network
nodes <- c("Habitat_Fragmentation", "Deforestation", "Agricultural_Expansion", "Biodiversity_Loss", 
           "Urbanization", "Urban_Green_Spaces", "Reforestation", "Grassland_Restoration", 
           "Wetland_Restoration", "Afforestation", "Wildlife_Abundance", "Rodent_Abundance", 
           "Tick_Population", "Human_Presence", "Wildlife_Migration", "Rodent_Borne_Diseases", 
           "Tick_Borne_Diseases", "Waterborne_Diseases", "Mosquito_Borne_Diseases", 
           "Avian_Diseases", "Disease_Spread", "Habitat_Design", "Public_Health_Campaigns", 
           "Connectivity_Existing_Forests", "Public_Health_Awareness", "Human_Activity_Restrictions", 
           "Wetland_Design", "Vector_Management", "Culling_of_Birds")

# Create an empty Bayesian network
bbn_structure <- empty.graph(nodes)

# Define arcs (directed dependencies between nodes)
arc.set <- data.frame(
  from = c(
    # Landscape Changes → Restoration Types
    "Habitat_Fragmentation", "Habitat_Fragmentation", "Deforestation", "Deforestation", 
    "Agricultural_Expansion", "Agricultural_Expansion", "Biodiversity_Loss", "Biodiversity_Loss",
    "Biodiversity_Loss", "Biodiversity_Loss", "Biodiversity_Loss", "Urbanization",
    
    # Restoration Types → Wildlife-Human Interactions
    "Wetland_Restoration", "Wetland_Restoration", "Reforestation", "Reforestation",
    "Urban_Green_Spaces", "Urban_Green_Spaces", "Grassland_Restoration", "Grassland_Restoration",
    "Afforestation", "Afforestation",
    
    # Wildlife-Human Interactions → Zoonotic Disease Risks
    "Wildlife_Abundance", "Wildlife_Abundance", "Human_Presence", "Human_Presence", 
    "Wildlife_Migration", "Wildlife_Migration",
    
    # Zoonotic Disease Risks → Mitigation Strategies
    "Rodent_Borne_Diseases", "Tick_Borne_Diseases", "Tick_Borne_Diseases", "Waterborne_Diseases",
    "Waterborne_Diseases", "Mosquito_Borne_Diseases", "Mosquito_Borne_Diseases",
    "Avian_Diseases", "Avian_Diseases"
  ),
  to = c(
    # Landscape Changes → Restoration Types
    "Urban_Green_Spaces", "Reforestation", "Reforestation", "Grassland_Restoration",
    "Grassland_Restoration", "Wetland_Restoration", "Urban_Green_Spaces", "Reforestation",
    "Grassland_Restoration", "Wetland_Restoration", "Afforestation", "Urban_Green_Spaces",
    
    # Restoration Types → Wildlife-Human Interactions
    "Wildlife_Abundance", "Human_Presence", "Wildlife_Migration", "Human_Presence",
    "Rodent_Abundance", "Wildlife_Migration", "Tick_Population", "Human_Presence",
    "Tick_Population", "Human_Presence",
    
    # Wildlife-Human Interactions → Zoonotic Disease Risks
    "Rodent_Borne_Diseases", "Tick_Borne_Diseases", "Waterborne_Diseases", "Mosquito_Borne_Diseases",
    "Disease_Spread", "Avian_Diseases",
    
    # Zoonotic Disease Risks → Mitigation Strategies
    "Habitat_Design", "Habitat_Design", "Connectivity_Existing_Forests", "Public_Health_Awareness",
    "Human_Activity_Restrictions", "Wetland_Design", "Vector_Management",
    "Public_Health_Awareness", "Culling_of_Birds"
  )
)

# Assign arcs to the Bayesian Network
arcs(bbn_structure) <- as.matrix(arc.set)

# Convert bnlearn arcs to igraph object
edge_list <- arcs(bbn_structure)

g <- graph_from_data_frame(edge_list, directed = TRUE)
# Plot
plot(g,
     vertex.size = 45,
     vertex.label.cex = 0.8,
     vertex.color = "lightblue",
     edge.arrow.size = 0.4,
     layout = layout_as_tree(g))




# Install if not already installed
install.packages("readxl")

# Load the package
library(readxl)

# Read the Excel file 
df <- read_excel("bbn_cpt.xlsx")

# Get unique nodes
nodes <- unique(c(df$parent, df$child))

# Create empty graph
bn <- empty.graph(nodes)

# Add arcs (edges)
arc.set <- data.frame(from = df$parent, to = df$child)
arcs(bn) <- as.matrix(arc.set)

# Define levels for nodes
levels <- c("low", "high")

# Combine all node names and states
nodes_from_parent <- unique(df[, c("parent", "parent_state")])
nodes_from_child  <- unique(df[, c("child", "child_state")])

# Combine and clean
combined_nodes <- rbind(
  setNames(nodes_from_parent, c("node", "state")),
  setNames(nodes_from_child,  c("node", "state"))
)

# Build node_states list with unique levels per node
node_states <- split(combined_nodes$state, combined_nodes$node)
node_states <- lapply(node_states, function(x) unique(na.omit(x)))

# View the result
print(node_states)
# Generate CPTs using midpoint of p_low and p_high
df <- df %>%
  mutate(
    p_low = as.numeric(gsub(",", ".", trimws(p_low))),
    p_high = as.numeric(gsub(",", ".", trimws(p_high)))
  )
str(df)
summary(df$p_low)
summary(df$p_high)
df <- df %>%
  mutate(prob = (p_low + p_high)/2)

# Get all unique nodes
all_nodes <- unique(c(df$parent, df$child))

cpts <- list()

for (child_node in unique(df$child)) {
  
  parents <- df %>%
    filter(child == child_node) %>%
    pull(parent) %>%
    unique()
  
  # PRIOR NODE (no parents)

  if (length(parents) == 0) {
    
    cpts[[child_node]] <- cptable(
      as.formula(paste0("~", child_node)),
      values = c(0.5, 0.5),
      levels = c("absent", "present")
    )
    
  } else {
    
    # Build formula
    formula_str <- paste0("~", child_node, "|", paste(parents, collapse = ":"))
    
    # Number of parent state combinations
    n_combinations <- 2^length(parents)
    
    # Pull probabilities from Excel
    probs <- df %>%
      filter(child == child_node) %>%
      pull(prob)
    
    # If Excel rows not exhaustive, fill safely
    expected_length <- n_combinations * 2
    
    if(length(probs) < expected_length){
      probs <- rep(mean(probs, na.rm = TRUE), expected_length)
    }
    
    probs <- probs[1:expected_length]
    
    # normalize each pair
    final_probs <- c()
    for(i in seq(1, length(probs), by = 2)){
      pair <- probs[i:(i+1)]
      if(sum(pair)==0) pair <- c(0.5,0.5)
      pair <- pair/sum(pair)
      final_probs <- c(final_probs, pair)
    }
    
    cpts[[child_node]] <- cptable(
      as.formula(formula_str),
      values = final_probs,
      levels = c("absent", "present")
    )
  }
}
all_nodes <- nodes(bn_structure)

for(n in all_nodes){
  if(is.null(cpts[[n]])){
    
    cpts[[n]] <- cptable(
      as.formula(paste0("~", n)),
      values = c(0.5,0.5),
      levels = c("absent","present")
    )
  }
}
# Compile network
plist <- compileCPT(cpts)
# Create an empty Bayesian network structure
bn_structure <- empty.graph(nodes = all_nodes)

# Fit the custom CPTs to this structure
fitted_bn <- custom.fit(bn_structure, plist)


# Create scenario tables
df$p_low <- as.numeric(df$p_low)
low_scenario <- xtabs(p_low ~ child + parent, data = df)
df$p_high <- as.numeric(df$p_high)
high_scenario <- xtabs(p_high ~ child + parent, data = df)

# Convert to CPTs
cpt_low <- as.table(low_scenario)
cpt_high <- as.table(high_scenario)

# Load packages
library(gRain)
library(gRbase)



# Extract all node names and their states
nodes_from_parent <- unique(df[, c("parent", "parent_state")])
nodes_from_child  <- unique(df[, c("child", "child_state")])

# Combine and rename
combined_nodes <- rbind(
  setNames(nodes_from_parent, c("node", "state")),
  setNames(nodes_from_child,  c("node", "state"))
)

# Step 3: Create the node_states list
node_states <- split(combined_nodes$state, combined_nodes$node)
node_states <- lapply(node_states, function(x) unique(na.omit(x)))

# Step 4: Create a list of all nodes
all_nodes <- names(node_states)

# Step 5: Create dummy CPTs using uniform distributions (replace later with actual probs)
cpts <- list()

for (node in all_nodes) {
  # Find parent(s) from df
  parents <- unique(df$parent[df$child == node])
  parents <- intersect(parents, all_nodes)  # ensure they're in node_states
  
  if (length(parents) == 0) {
    # No parents — prior probability
    values <- rep(1 / length(node_states[[node]]), length(node_states[[node]]))
    cpts[[node]] <- cptable(as.formula(paste0("~", node)), values = values, levels = node_states[[node]])
  } else {
    # Has parents — define CPT shape
    formula_str <- paste("~", paste(c(node, parents), collapse = "+"))
    
    n_combinations <- prod(sapply(parents, function(p) length(node_states[[p]])))
    values <- rep(1 / length(node_states[[node]]), length(node_states[[node]]) * n_combinations)
    
    cpts[[node]] <- cptable(as.formula(formula_str), values = values, levels = c(node_states[[node]], unlist(node_states[parents])))
  }
}

# Step 6: Compile the CPTs into a grain object
plist <- compileCPT(cpts)
fitted_bn <- grain(plist)

# Plot the structure
plot(fitted_bn)
library(igraph)
# -----------------------------
##for pretty plot
# -----------------------------
# Build igraph object from bn structure
edge_list <- arcs(bn)
g <- graph_from_data_frame(edge_list, directed = TRUE)

# Define manual layers

layer1 <- c("Habitat_Fragmentation", "Deforestation", "Agricultural_Expansion",
            "Biodiversity_Loss", "Urbanization", "Region_Context")

layer2 <- c("Urban_Green_Spaces", "Reforestation", "Grassland_Restoration",
            "Wetland_Restoration", "Afforestation")

layer3 <- c("Wildlife_Abundance", "Rodent_Abundance", "Tick_Population",
            "Human_Presence", "Wildlife_Movement")

layer4 <- c("Rodent_Borne_Diseases", "Tick_Borne_Diseases", "Waterborne_Diseases",
            "Mosquito_Borne_Diseases", "Avian_Diseases", "Disease_Spread")

layer5 <- c("Habitat_Design", "Public_Health_Campaigns", "Connectivity_Existing_Forests",
            "Public_Health_Awareness", "Human_Activity_Restrictions", "Wetland_Design",
            "Vector_Management", "Culling_of_Birds", "Time_to_Stable_State")

# Create coordinate matrix
all_nodes <- V(g)$name
coords <- matrix(NA, nrow = length(all_nodes), ncol = 2)
rownames(coords) <- all_nodes

place_layer <- function(nodes_layer, y_val){
  
  nodes_layer <- nodes_layer[nodes_layer %in% rownames(coords)]
  
  x_vals <- seq(-1, 1, length.out = length(nodes_layer))
  
  for(i in seq_along(nodes_layer)){
    coords[nodes_layer[i], ] <<- c(x_vals[i], y_val)
  }
}

place_layer(layer1, 5)
place_layer(layer2, 4)
place_layer(layer3, 3)
place_layer(layer4, 2)
place_layer(layer5, 1)

# Plot

plot(g,
     layout = coords,
     vertex.size = 22,
     vertex.color = "lightblue",
     vertex.frame.color = "grey40",
     vertex.label.cex = 0.3,
     vertex.label.family = "sans",
     vertex.label.dist = 0.2,
     edge.arrow.size = 0.3,
     edge.width = 1.2,
     edge.curved = 0.05,
     margin = 0.15)


##add region specific prior
c("Low_Risk_Context", "Moderate_Risk_Context", "High_Risk_Context")

# Add region node
node_states[["Region_Context"]] <-
  c("Low_Risk_Context", "Moderate_Risk_Context", "High_Risk_Context")

nodes <- c(nodes, "Region_Context")

# Update structure
bn_structure <- empty.graph(nodes)
arcs(bn_structure) <- as.matrix(arcs_df)

# Region influences baseline disease pressure
region_arcs <- matrix(c(
  "Region_Context", "Rodentborne_Diseases",
  "Region_Context", "Mosquitoborne_Diseases",
  "Region_Context", "Tickborne_Diseases"
), byrow = TRUE, ncol = 2)

arcs(bn_structure) <- rbind(arcs(bn_structure), region_arcs)

cpt_rodent_best <- cptable(
  ~ Rodentborne_Diseases | Region_Context,
  values = c(
    0.8, 0.2,   # Low-risk region
    0.6, 0.4,   # Moderate
    0.4, 0.6    # High-risk
  ),
  levels = c("absent", "present")
)

cpt_rodent_worst <- cptable(
  ~ Rodentborne_Diseases | Region_Context,
  values = c(
    0.6, 0.4,
    0.4, 0.6,
    0.2, 0.8
  ),
  levels = c("absent", "present")
)


cpt_time_best <- cptable(
  ~ Time_to_Stable_State |
    Region_Context :
    Habitat_Design :
    Vector_Management_Strategies,
  values = c(
    # Low-risk region
    0.55, 0.30, 0.10, 0.05,
    # Moderate-risk
    0.30, 0.35, 0.25, 0.10,
    # High-risk
    0.10, 0.25, 0.35, 0.30
  ),
  levels = node_states[["Time_to_Stable_State"]]
)

bn_best <- grain(compileCPT(cpts_best))
bn_worst <- grain(compileCPT(cpts_worst))

# Condition on region
bn_best_high <- setEvidence(bn_best, nodes = "Region_Context", states = "High_Risk_Context")
bn_worst_high <- setEvidence(bn_worst, nodes = "Region_Context", states = "High_Risk_Context")

best_dist <- querygrain(bn_best_high, nodes = "Time_to_Stable_State")$Time_to_Stable_State
worst_dist <- querygrain(bn_worst_high, nodes = "Time_to_Stable_State")$Time_to_Stable_State

plot_df <- data.frame(
  Time = factor(names(best_dist), levels = names(best_dist)),
  Best_Case = as.numeric(best_dist),
  Worst_Case = as.numeric(worst_dist)
)

library(ggplot2)
library(tidyr)

plot_df_long <- pivot_longer(
  plot_df,
  cols = c("Best_Case", "Worst_Case"),
  names_to = "Scenario",
  values_to = "Probability"
)

ggplot(plot_df_long, aes(x = Time, y = Probability, fill = Scenario)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("#4CAF50", "#D32F2F")) +
  labs(
    title = "Time to Reach a Stable Low Zoonotic Risk State",
    subtitle = "Comparison of Best- and Worst-Case Restoration Scenarios (High-Risk Regions)",
    x = "Years After Restoration",
    y = "Probability",
    fill = "Scenario"
  ) +
  theme_minimal(base_size = 13)
