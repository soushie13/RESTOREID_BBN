# =========================================================
# FUNDING SCENARIO SIMULATIONS
# SEQUENTIAL RESTORATION TRANSITIONS ACROSS FOUR DISEASE SYSTEMS
#
# Simulates two real-world funding scenarios:
#
# SCENARIO A — "Ramp up"
#   Phase 1 (Years 1-5 or 1-10): low or moderate restoration
#   Phase 2 (remaining years to Year 15): high restoration
#   Reflects: project starts under-resourced, funding secured later
#
# SCENARIO B — "Funding collapse"
#   Phase 1 (Years 1-5 or 1-10): high restoration
#   Phase 2 (remaining years to Year 15): low restoration
#   Reflects: well-funded project loses funding mid-way
#
# MODELLING APPROACH — Ecological legacy weighting
#
# At the transition point, the ecosystem is NOT in a pristine
# state — it carries the ecological legacy of Phase 1.
# This is modelled by computing a weighted mixture of:
#   (a) the disease risk trajectory under Phase 1 restoration
#       continued to Year 15 (no transition counterfactual)
#   (b) the disease risk trajectory under Phase 2 restoration
#       applied from Year 1 (ideal counterfactual)
#
# The mixture weight is determined by how far along the
# Time_to_Stable_State distribution Phase 1 has progressed
# by the transition year — ecosystems that have stabilised
# more under Phase 1 carry more ecological legacy into Phase 2.
#
# Formally:
#   risk_transition(t) =
#     legacy_weight(t) * risk_phase1 +
#     (1 - legacy_weight(t)) * risk_phase2
#
# where legacy_weight(t) = cumulative P(stabilised by year t)
# under Phase 1 restoration, derived from Time_to_Stable_State.
#
# This produces ecologically realistic trajectories:
# - Ecosystems that stabilised early under Phase 1 are harder
#   to redirect in Phase 2 (high legacy weight)
# - Ecosystems still in flux at the transition respond more
#   readily to the Phase 2 restoration intensity change
#
# =========================================================

library(bnlearn)
library(gRain)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

set.seed(123)

# =========================================================
# SHARED HELPERS
# =========================================================

sample_prob <- function(mn, mx) runif(1, mn, mx)
norm        <- function(x){ x[x < 0.01] <- 0.01; x / sum(x) }
norm_groups <- function(x, n){
  unlist(lapply(split(x, ceiling(seq_along(x)/n)), norm))
}

# =========================================================
# TRANSITION YEAR DEFINITIONS
# =========================================================

transition_years <- c(5, 10)   # simulate both
total_years      <- 15

# Cumulative probability of stabilisation by transition year
# derived from Time_to_Stable_State distributions.
# Windows: 1-3 yr, 4-7 yr, 8-12 yr, 13-15 yr
# Assign midpoints: 2, 5.5, 10, 14
window_midpoints <- c(2, 5.5, 10, 14)

get_legacy_weight <- function(time_probs, transition_year){
  # Proportion of ecosystem expected to have stabilised
  # by the transition year under Phase 1 restoration
  # = sum of P(window) for windows whose midpoint <= transition year
  sum(time_probs[window_midpoints <= transition_year])
}

# =========================================================
# CORE FUNCTION: run BBN for one system across scenarios
# Returns a data frame of disease risk per scenario
# =========================================================

run_funding_scenarios <- function(system_name,
                                  build_bn_fn,
                                  query_node,
                                  query_state,
                                  ev_nodes,
                                  ev_states_template,
                                  restoration_states,
                                  n_iter = 1000){
  
  # Scenarios to run:
  # - Pure scenarios (each intensity for full 15 years)
  # - Ramp-up: low→high at year 5, low→high at year 10
  # - Ramp-up: moderate→high at year 5, moderate→high at year 10
  # - Collapse: high→low at year 5, high→low at year 10
  # - Collapse: high→moderate at year 5, high→moderate at year 10
  
  scenarios <- list(
    # Pure baselines
    list(name="None (15 yr)",         p1="none",     p2=NA,    t=NA),
    list(name="Low (15 yr)",          p1="low",      p2=NA,    t=NA),
    list(name="Moderate (15 yr)",     p1="moderate", p2=NA,    t=NA),
    list(name="High (15 yr)",         p1="high",     p2=NA,    t=NA),
    
    # Ramp-up scenarios
    list(name="Low→High (yr 5)",      p1="low",      p2="high",     t=5),
    list(name="Low→High (yr 10)",     p1="low",      p2="high",     t=10),
    list(name="Moderate→High (yr 5)", p1="moderate", p2="high",     t=5),
    list(name="Moderate→High (yr 10)",p1="moderate", p2="high",     t=10),
    
    # Funding collapse scenarios
    list(name="High→Low (yr 5)",      p1="high",     p2="low",      t=5),
    list(name="High→Low (yr 10)",     p1="high",     p2="low",      t=10),
    list(name="High→Moderate (yr 5)", p1="high",     p2="moderate", t=5),
    list(name="High→Moderate (yr 10)",p1="high",     p2="moderate", t=10)
  )
  
  all_results <- data.frame()
  
  for(i in seq_len(n_iter)){
    
    bn <- build_bn_fn()
    
    # Query disease risk and time distribution for each
    # pure restoration intensity
    pure_risk <- list()
    pure_time <- list()
    
    for(rest in c("none","low","moderate","high")){
      states_i <- ev_states_template
      states_i[ev_nodes == "Restoration_Intensity"] <- rest
      
      bn_ev <- setEvidence(bn,
                           nodes  = ev_nodes,
                           states = states_i)
      
      pure_risk[[rest]] <- querygrain(
        bn_ev, nodes = query_node)[[query_node]][query_state]
      
      pure_time[[rest]] <- querygrain(
        bn_ev, nodes = "Time_to_Stable_State")$Time_to_Stable_State
    }
    
    # Compute scenario risk using ecological legacy weighting
    for(sc in scenarios){
      
      if(is.na(sc$p2)){
        # Pure scenario — no transition
        risk_i <- pure_risk[[sc$p1]]
        type_i <- "Pure baseline"
        
      } else {
        
        # Legacy weight = how locked-in Phase 1 state is
        # at the transition year
        lw <- get_legacy_weight(pure_time[[sc$p1]], sc$t)
        
        # Phase 2 risk under new intensity
        risk_p2 <- pure_risk[[sc$p2]]
        
        # Weighted mixture: legacy from Phase 1, new trajectory Phase 2
        # Remaining years after transition = total_years - sc$t
        # Phase 2 has less time to take effect if transition is late
        time_discount <- (total_years - sc$t) / total_years
        
        risk_i <- lw * pure_risk[[sc$p1]] +
          (1 - lw) * (time_discount * risk_p2 +
                        (1 - time_discount) * pure_risk[[sc$p1]])
        
        type_i <- if(grepl("→High", sc$name)) "Ramp-up" else "Funding collapse"
      }
      
      all_results <- rbind(all_results, data.frame(
        iteration   = i,
        system      = system_name,
        scenario    = sc$name,
        phase1      = sc$p1,
        phase2      = if(is.na(sc$p2)) sc$p1 else sc$p2,
        transition  = if(is.na(sc$t)) 0 else sc$t,
        scenario_type = type_i,
        disease_risk  = risk_i,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  all_results
}

# =========================================================
# SYSTEM-SPECIFIC BBN BUILDERS
# Each returns a compiled gRain object ready for querying
# =========================================================

# ── MOSQUITO ──────────────────────────────────────────────

build_mosquito_bn <- function(){
  states <- list(
    Region_Context        = c("Low","Moderate","High"),
    Restoration_Intensity = c("none","low","moderate","high"),
    Hydrological_Recovery = c("low","moderate","high"),
    Standing_Water        = c("low","moderate","high"),
    Mosquito_Disease_Risk = c("absent","present"),
    Ecological_Regulation = c("weak","moderate","strong"),
    Time_to_Stable_State  = c("1_3_years","4_7_years","8_12_years","13_15_years")
  )
  dag <- model2network(paste0(
    "[Region_Context][Restoration_Intensity]",
    "[Hydrological_Recovery|Region_Context:Restoration_Intensity]",
    "[Standing_Water|Hydrological_Recovery:Restoration_Intensity]",
    "[Ecological_Regulation|Restoration_Intensity]",
    "[Mosquito_Disease_Risk|Standing_Water:Ecological_Regulation]",
    "[Time_to_Stable_State|Ecological_Regulation]"
  ))
  
  hyd_vals <- c()
  for(region in 1:3) for(rest in 1:4){
    probs <- switch(rest,
                    `1`=c(sample_prob(.60,.80),sample_prob(.15,.30),sample_prob(.05,.15)),
                    `2`=c(sample_prob(.20,.35),sample_prob(.35,.50),sample_prob(.25,.40)),
                    `3`=c(sample_prob(.05,.20),sample_prob(.25,.40),sample_prob(.45,.65)),
                    `4`=c(sample_prob(.05,.15),sample_prob(.20,.35),sample_prob(.55,.75)))
    hyd_vals <- c(hyd_vals, norm(probs))
  }
  
  sw_vals <- c()
  for(rest in 1:4) for(hydro in 1:3){
    probs <- switch(rest,
                    `1`=c(sample_prob(.25,.40),sample_prob(.40,.55),sample_prob(.10,.20)),
                    `2`=c(sample_prob(.10,.20),sample_prob(.35,.50),sample_prob(.35,.50)),
                    `3`=c(sample_prob(.02,.08),sample_prob(.15,.28),sample_prob(.65,.82)),
                    `4`=c(sample_prob(.35,.55),sample_prob(.30,.45),sample_prob(.05,.15)))
    if(hydro==1) probs[1] <- probs[1]+.10
    if(hydro==3) probs[3] <- probs[3]+.08
    sw_vals <- c(sw_vals, norm(probs))
  }
  
  reg_vals <- norm_groups(c(
    sample_prob(.78,.84),sample_prob(.13,.15),sample_prob(.03,.05),
    sample_prob(.58,.64),sample_prob(.24,.28),sample_prob(.09,.11),
    sample_prob(.25,.27),sample_prob(.42,.44),sample_prob(.29,.31),
    sample_prob(.05,.07),sample_prob(.15,.17),sample_prob(.72,.78)
  ),3)
  
  dis_vals <- c()
  for(reg in 1:3) for(sw in 1:3){
    probs <- if(sw==3&reg==1) c(sample_prob(.05,.12),sample_prob(.88,.95))
    else if(sw==3&reg==2)     c(sample_prob(.25,.40),sample_prob(.60,.75))
    else if(sw==3&reg==3)     c(sample_prob(.82,.92),sample_prob(.08,.18))
    else if(sw==2&reg==1)     c(sample_prob(.30,.45),sample_prob(.55,.70))
    else if(sw==2&reg==2)     c(sample_prob(.45,.65),sample_prob(.35,.55))
    else if(sw==2&reg==3)     c(sample_prob(.78,.90),sample_prob(.10,.22))
    else if(sw==1&reg==1)     c(sample_prob(.50,.65),sample_prob(.35,.50))
    else if(sw==1&reg==2)     c(sample_prob(.70,.82),sample_prob(.18,.30))
    else                      c(sample_prob(.88,.96),sample_prob(.04,.12))
    dis_vals <- c(dis_vals, norm(probs))
  }
  
  time_vals <- c()
  for(reg in 1:3){
    probs <- switch(reg,
                    `1`=c(sample_prob(.38,.52),sample_prob(.28,.38),sample_prob(.08,.18),sample_prob(.03,.10)),
                    `2`=c(sample_prob(.14,.24),sample_prob(.32,.44),sample_prob(.22,.34),sample_prob(.10,.22)),
                    `3`=c(sample_prob(.03,.09),sample_prob(.10,.20),sample_prob(.42,.56),sample_prob(.24,.38)))
    time_vals <- c(time_vals, norm(probs))
  }
  
  fit <- custom.fit(dag, dist=list(
    Region_Context        = array(c(.33,.34,.33),dim=c(3),dimnames=list(Region_Context=states$Region_Context)),
    Restoration_Intensity = array(rep(.25,4),dim=c(4),dimnames=list(Restoration_Intensity=states$Restoration_Intensity)),
    Hydrological_Recovery = array(hyd_vals,dim=c(3,4,3),dimnames=list(Hydrological_Recovery=states$Hydrological_Recovery,Restoration_Intensity=states$Restoration_Intensity,Region_Context=states$Region_Context)),
    Standing_Water        = array(sw_vals,dim=c(3,3,4),dimnames=list(Standing_Water=states$Standing_Water,Hydrological_Recovery=states$Hydrological_Recovery,Restoration_Intensity=states$Restoration_Intensity)),
    Ecological_Regulation = array(reg_vals,dim=c(3,4),dimnames=list(Ecological_Regulation=states$Ecological_Regulation,Restoration_Intensity=states$Restoration_Intensity)),
    Mosquito_Disease_Risk = array(dis_vals,dim=c(2,3,3),dimnames=list(Mosquito_Disease_Risk=states$Mosquito_Disease_Risk,Standing_Water=states$Standing_Water,Ecological_Regulation=states$Ecological_Regulation)),
    Time_to_Stable_State  = array(time_vals,dim=c(4,3),dimnames=list(Time_to_Stable_State=states$Time_to_Stable_State,Ecological_Regulation=states$Ecological_Regulation))
  ))
  as.grain(fit)
}

# ── TICK ──────────────────────────────────────────────────

build_tick_bn <- function(){
  states <- list(
    Region_Context        = c("Low","Moderate","High"),
    Restoration_Intensity = c("none","low","moderate","high"),
    Host_Connectivity     = c("low","moderate","high"),
    Predator_Recovery     = c("weak","moderate","strong"),
    Tick_Disease_Risk     = c("absent","present"),
    Ecological_Regulation = c("weak","moderate","strong"),
    Time_to_Stable_State  = c("1_3_years","4_7_years","8_12_years","13_15_years")
  )
  dag <- model2network(paste0(
    "[Region_Context][Restoration_Intensity]",
    "[Host_Connectivity|Restoration_Intensity]",
    "[Predator_Recovery|Restoration_Intensity]",
    "[Tick_Disease_Risk|Host_Connectivity:Predator_Recovery]",
    "[Ecological_Regulation|Predator_Recovery:Tick_Disease_Risk]",
    "[Time_to_Stable_State|Ecological_Regulation]"
  ))
  
  con_vals <- norm_groups(c(
    sample_prob(.55,.72),sample_prob(.18,.30),sample_prob(.05,.15),
    sample_prob(.05,.12),sample_prob(.18,.30),sample_prob(.60,.78),
    sample_prob(.08,.18),sample_prob(.22,.35),sample_prob(.50,.68),
    sample_prob(.15,.30),sample_prob(.38,.52),sample_prob(.25,.42)
  ),3)
  
  pred_vals <- norm_groups(c(
    sample_prob(.78,.86),sample_prob(.10,.18),sample_prob(.02,.06),
    sample_prob(.68,.76),sample_prob(.18,.26),sample_prob(.04,.08),
    sample_prob(.48,.56),sample_prob(.29,.37),sample_prob(.11,.19),
    sample_prob(.08,.16),sample_prob(.24,.32),sample_prob(.56,.64)
  ),3)
  
  risk_vals <- c()
  for(pred in 1:3) for(con in 1:3){
    probs <- if(con==3&pred==1) c(sample_prob(.05,.14),sample_prob(.86,.95))
    else if(con==3&pred==3)     c(sample_prob(.68,.82),sample_prob(.18,.32))
    else if(con==1&pred==1)     c(sample_prob(.42,.60),sample_prob(.40,.58))
    else if(con==1&pred==3)     c(sample_prob(.75,.90),sample_prob(.10,.25))
    else if(con==2&pred==2)     c(sample_prob(.30,.70),sample_prob(.30,.70))
    else                        c(sample_prob(.28,.68),sample_prob(.32,.72))
    risk_vals <- c(risk_vals, norm(probs))
  }
  
  reg_vals <- c()
  for(disease in 1:2) for(pred in 1:3){
    probs <- if(pred==3&disease==1) c(sample_prob(.04,.10),sample_prob(.15,.26),sample_prob(.66,.80))
    else if(pred==1&disease==2)     c(sample_prob(.65,.80),sample_prob(.15,.26),sample_prob(.04,.10))
    else                            c(sample_prob(.18,.48),sample_prob(.22,.48),sample_prob(.18,.48))
    reg_vals <- c(reg_vals, norm(probs))
  }
  
  time_vals <- c()
  for(reg in 1:3){
    probs <- switch(reg,
                    `1`=c(sample_prob(.38,.52),sample_prob(.28,.38),sample_prob(.08,.18),sample_prob(.03,.10)),
                    `2`=c(sample_prob(.15,.28),sample_prob(.22,.32),sample_prob(.24,.34),sample_prob(.20,.32)),
                    `3`=c(sample_prob(.02,.06),sample_prob(.05,.12),sample_prob(.22,.34),sample_prob(.52,.68)))
    time_vals <- c(time_vals, norm(probs))
  }
  
  fit <- custom.fit(dag, dist=list(
    Region_Context        = array(c(.33,.34,.33),dim=c(3),dimnames=list(Region_Context=states$Region_Context)),
    Restoration_Intensity = array(rep(.25,4),dim=c(4),dimnames=list(Restoration_Intensity=states$Restoration_Intensity)),
    Host_Connectivity     = array(con_vals,dim=c(3,4),dimnames=list(Host_Connectivity=states$Host_Connectivity,Restoration_Intensity=states$Restoration_Intensity)),
    Predator_Recovery     = array(pred_vals,dim=c(3,4),dimnames=list(Predator_Recovery=states$Predator_Recovery,Restoration_Intensity=states$Restoration_Intensity)),
    Tick_Disease_Risk     = array(risk_vals,dim=c(2,3,3),dimnames=list(Tick_Disease_Risk=states$Tick_Disease_Risk,Host_Connectivity=states$Host_Connectivity,Predator_Recovery=states$Predator_Recovery)),
    Ecological_Regulation = array(reg_vals,dim=c(3,3,2),dimnames=list(Ecological_Regulation=states$Ecological_Regulation,Predator_Recovery=states$Predator_Recovery,Tick_Disease_Risk=states$Tick_Disease_Risk)),
    Time_to_Stable_State  = array(time_vals,dim=c(4,3),dimnames=list(Time_to_Stable_State=states$Time_to_Stable_State,Ecological_Regulation=states$Ecological_Regulation))
  ))
  as.grain(fit)
}

# ── RODENT ────────────────────────────────────────────────

build_rodent_bn <- function(){
  states <- list(
    Region_Context        = c("Low","Moderate","High"),
    Restoration_Intensity = c("none","low","moderate","high"),
    Habitat_Quality       = c("low","moderate","high"),
    Rodent_Abundance      = c("low","moderate","high"),
    Rodent_Disease_Risk   = c("absent","present"),
    Ecological_Regulation = c("weak","moderate","strong"),
    Time_to_Stable_State  = c("1_3_years","4_7_years","8_12_years","13_15_years")
  )
  dag <- model2network(paste0(
    "[Region_Context][Restoration_Intensity]",
    "[Habitat_Quality|Region_Context:Restoration_Intensity]",
    "[Rodent_Abundance|Habitat_Quality]",
    "[Rodent_Disease_Risk|Rodent_Abundance]",
    "[Ecological_Regulation|Restoration_Intensity:Rodent_Disease_Risk]",
    "[Time_to_Stable_State|Ecological_Regulation]"
  ))
  
  hab_vals <- c()
  for(region in 1:3) for(rest in 1:4){
    probs <- switch(rest,
                    `1`=c(sample_prob(.62,.80),sample_prob(.14,.28),sample_prob(.02,.10)),
                    `2`=c(sample_prob(.38,.58),sample_prob(.26,.40),sample_prob(.08,.20)),
                    `3`=c(sample_prob(.18,.35),sample_prob(.32,.46),sample_prob(.25,.42)),
                    `4`=c(sample_prob(.04,.14),sample_prob(.18,.34),sample_prob(.55,.75)))
    if(region==3) probs[1] <- probs[1]+.10
    hab_vals <- c(hab_vals, norm(probs))
  }
  
  rod_vals <- norm_groups(c(
    .05,.18,.77, .28,.48,.24, .72,.22,.06
  ),3)
  
  dis_vals <- norm_groups(c(.92,.08, .55,.45, .10,.90),2)
  
  reg_vals <- c()
  for(disease in 1:2) for(rest in 1:4){
    probs <- switch(rest,
                    `1`=c(sample_prob(.65,.82),sample_prob(.12,.24),sample_prob(.02,.10)),
                    `2`=c(sample_prob(.45,.62),sample_prob(.24,.38),sample_prob(.08,.22)),
                    `3`=c(sample_prob(.22,.38),sample_prob(.30,.44),sample_prob(.25,.42)),
                    `4`=c(sample_prob(.04,.12),sample_prob(.12,.26),sample_prob(.65,.82)))
    if(disease==2) probs[1] <- probs[1]+.05
    reg_vals <- c(reg_vals, norm(probs))
  }
  
  time_vals <- c()
  for(reg in 1:3){
    probs <- switch(reg,
                    `1`=c(sample_prob(.48,.62),sample_prob(.24,.34),sample_prob(.06,.14),sample_prob(.02,.08)),
                    `2`=c(sample_prob(.14,.24),sample_prob(.32,.44),sample_prob(.22,.34),sample_prob(.10,.22)),
                    `3`=c(sample_prob(.03,.09),sample_prob(.10,.20),sample_prob(.42,.56),sample_prob(.24,.38)))
    time_vals <- c(time_vals, norm(probs))
  }
  
  fit <- custom.fit(dag, dist=list(
    Region_Context        = array(c(.33,.34,.33),dim=c(3),dimnames=list(Region_Context=states$Region_Context)),
    Restoration_Intensity = array(rep(.25,4),dim=c(4),dimnames=list(Restoration_Intensity=states$Restoration_Intensity)),
    Habitat_Quality       = array(hab_vals,dim=c(3,4,3),dimnames=list(Habitat_Quality=states$Habitat_Quality,Restoration_Intensity=states$Restoration_Intensity,Region_Context=states$Region_Context)),
    Rodent_Abundance      = array(rod_vals,dim=c(3,3),dimnames=list(Rodent_Abundance=states$Rodent_Abundance,Habitat_Quality=states$Habitat_Quality)),
    Rodent_Disease_Risk   = array(dis_vals,dim=c(2,3),dimnames=list(Rodent_Disease_Risk=states$Rodent_Disease_Risk,Rodent_Abundance=states$Rodent_Abundance)),
    Ecological_Regulation = array(reg_vals,dim=c(3,4,2),dimnames=list(Ecological_Regulation=states$Ecological_Regulation,Restoration_Intensity=states$Restoration_Intensity,Rodent_Disease_Risk=states$Rodent_Disease_Risk)),
    Time_to_Stable_State  = array(time_vals,dim=c(4,3),dimnames=list(Time_to_Stable_State=states$Time_to_Stable_State,Ecological_Regulation=states$Ecological_Regulation))
  ))
  as.grain(fit)
}

# ── DIRECT-CONTACT ZOONOSES ───────────────────────────────

build_zoonotic_bn <- function(){
  states <- list(
    Region_Context                  = c("Low_Risk_Context","Moderate_Risk_Context","High_Risk_Context"),
    Restoration_Intensity           = c("none","low","moderate","high"),
    Resilient_Ecological_Regulation = c("low","high"),
    Habitat_Disturbance             = c("low","moderate","high"),
    Host_Diversity                  = c("low","moderate","high"),
    Human_Wildlife_Contact          = c("low","moderate","high"),
    MultiHost_Amplification         = c("low","high"),
    Direct_Transmission_Risk        = c("absent","present"),
    Time_to_Stable_State            = c("1_3_years","4_7_years","8_12_years","13_15_years")
  )
  dag <- model2network(paste0(
    "[Region_Context][Restoration_Intensity]",
    "[Habitat_Disturbance|Region_Context:Restoration_Intensity]",
    "[Host_Diversity|Habitat_Disturbance]",
    "[Human_Wildlife_Contact|Habitat_Disturbance:Restoration_Intensity]",
    "[MultiHost_Amplification|Host_Diversity:Human_Wildlife_Contact]",
    "[Direct_Transmission_Risk|MultiHost_Amplification]",
    "[Resilient_Ecological_Regulation|Restoration_Intensity:Host_Diversity]",
    "[Time_to_Stable_State|Restoration_Intensity:Direct_Transmission_Risk:Resilient_Ecological_Regulation]"
  ))
  
  dist_vals <- c()
  for(region in 1:3) for(rest in 1:4){
    probs <- switch(rest,
                    `1`=c(sample_prob(.05,.15),sample_prob(.20,.30),sample_prob(.55,.75)),
                    `2`=c(sample_prob(.10,.25),sample_prob(.35,.50),sample_prob(.30,.45)),
                    `3`=c(sample_prob(.20,.40),sample_prob(.35,.45),sample_prob(.15,.30)),
                    `4`=c(sample_prob(.50,.70),sample_prob(.20,.35),sample_prob(.05,.15)))
    if(region==3) probs[3] <- probs[3]+.10
    dist_vals <- c(dist_vals, norm(probs))
  }
  
  div_vals <- norm_groups(c(.10,.25,.65,.30,.45,.25,.70,.20,.10),3)
  
  con_vals <- c()
  for(rest in 1:4) for(disturbance in 1:3){
    probs <- switch(rest,
                    `1`=c(sample_prob(.10,.20),sample_prob(.25,.35),sample_prob(.50,.65)),
                    `2`=c(sample_prob(.15,.30),sample_prob(.35,.45),sample_prob(.30,.45)),
                    `3`=c(sample_prob(.30,.50),sample_prob(.30,.42),sample_prob(.15,.28)),
                    `4`=c(sample_prob(.55,.75),sample_prob(.18,.30),sample_prob(.05,.15)))
    if(disturbance==3) probs[3] <- probs[3]+.10
    con_vals <- c(con_vals, norm(probs))
  }
  
  amp_vals <- c()
  for(contact in 1:3) for(diversity in 1:3){
    probs <- if(contact==3&diversity==1) c(sample_prob(.05,.15),sample_prob(.85,.95))
    else if(contact==1&diversity==3)     c(sample_prob(.70,.88),sample_prob(.12,.30))
    else if(contact==3&diversity==3)     c(sample_prob(.30,.55),sample_prob(.45,.70))
    else if(contact==2&diversity==2)     c(sample_prob(.35,.65),sample_prob(.35,.65))
    else                                 c(sample_prob(.30,.60),sample_prob(.40,.70))
    amp_vals <- c(amp_vals, norm(probs))
  }
  
  disease_vals <- norm_groups(c(.80,.20,.15,.85),2)
  
  reg_vals <- c()
  for(diversity in 1:3) for(rest in 1:4){
    probs <- switch(rest,
                    `1`=c(sample_prob(.68,.85),sample_prob(.15,.32)),
                    `2`=c(sample_prob(.55,.72),sample_prob(.28,.45)),
                    `3`=c(sample_prob(.38,.58),sample_prob(.42,.62)),
                    `4`=c(sample_prob(.15,.35),sample_prob(.65,.85)))
    if(diversity==3) probs[2] <- probs[2]+sample_prob(.03,.12)
    if(diversity==1) probs[1] <- probs[1]+sample_prob(.02,.08)
    reg_vals <- c(reg_vals, norm(probs))
  }
  
  time_vals <- c()
  for(disease in 1:2) for(rest in 1:4) for(regulation in 1:2){
    probs <- if(regulation==2&rest==4) norm(c(sample_prob(.03,.10),sample_prob(.10,.20),sample_prob(.25,.38),sample_prob(.38,.55)))
    else if(rest==1)                   norm(c(sample_prob(.38,.52),sample_prob(.25,.35),sample_prob(.10,.20),sample_prob(.04,.12)))
    else                               norm(c(sample_prob(.12,.32),sample_prob(.18,.35),sample_prob(.20,.35),sample_prob(.18,.35)))
    time_vals <- c(time_vals, probs)
  }
  
  fit <- custom.fit(dag, dist=list(
    Region_Context                  = array(c(.33,.34,.33),dim=c(3),dimnames=list(Region_Context=states$Region_Context)),
    Restoration_Intensity           = array(rep(.25,4),dim=c(4),dimnames=list(Restoration_Intensity=states$Restoration_Intensity)),
    Habitat_Disturbance             = array(dist_vals,dim=c(3,4,3),dimnames=list(Habitat_Disturbance=states$Habitat_Disturbance,Restoration_Intensity=states$Restoration_Intensity,Region_Context=states$Region_Context)),
    Host_Diversity                  = array(div_vals,dim=c(3,3),dimnames=list(Host_Diversity=states$Host_Diversity,Habitat_Disturbance=states$Habitat_Disturbance)),
    Human_Wildlife_Contact          = array(con_vals,dim=c(3,3,4),dimnames=list(Human_Wildlife_Contact=states$Human_Wildlife_Contact,Habitat_Disturbance=states$Habitat_Disturbance,Restoration_Intensity=states$Restoration_Intensity)),
    MultiHost_Amplification         = array(amp_vals,dim=c(2,3,3),dimnames=list(MultiHost_Amplification=states$MultiHost_Amplification,Host_Diversity=states$Host_Diversity,Human_Wildlife_Contact=states$Human_Wildlife_Contact)),
    Direct_Transmission_Risk        = array(disease_vals,dim=c(2,2),dimnames=list(Direct_Transmission_Risk=states$Direct_Transmission_Risk,MultiHost_Amplification=states$MultiHost_Amplification)),
    Resilient_Ecological_Regulation = array(reg_vals,dim=c(2,4,3),dimnames=list(Resilient_Ecological_Regulation=states$Resilient_Ecological_Regulation,Restoration_Intensity=states$Restoration_Intensity,Host_Diversity=states$Host_Diversity)),
    Time_to_Stable_State            = array(time_vals,dim=c(4,2,4,2),dimnames=list(Time_to_Stable_State=states$Time_to_Stable_State,Resilient_Ecological_Regulation=states$Resilient_Ecological_Regulation,Restoration_Intensity=states$Restoration_Intensity,Direct_Transmission_Risk=states$Direct_Transmission_Risk))
  ))
  as.grain(fit)
}

# =========================================================
# RUN ALL FOUR SYSTEMS
# =========================================================

mosq_funding <- run_funding_scenarios(
  system_name       = "Mosquito-borne",
  build_bn_fn       = build_mosquito_bn,
  query_node        = "Mosquito_Disease_Risk",
  query_state       = "present",
  ev_nodes          = c("Region_Context","Restoration_Intensity"),
  ev_states_template = c("High","moderate"),
  restoration_states = c("none","low","moderate","high")
)

tick_funding <- run_funding_scenarios(
  system_name       = "Tick-borne",
  build_bn_fn       = build_tick_bn,
  query_node        = "Tick_Disease_Risk",
  query_state       = "present",
  ev_nodes          = c("Region_Context","Restoration_Intensity"),
  ev_states_template = c("High","moderate"),
  restoration_states = c("none","low","moderate","high")
)

rodent_funding <- run_funding_scenarios(
  system_name       = "Rodent-borne",
  build_bn_fn       = build_rodent_bn,
  query_node        = "Rodent_Disease_Risk",
  query_state       = "present",
  ev_nodes          = c("Region_Context","Restoration_Intensity"),
  ev_states_template = c("High","moderate"),
  restoration_states = c("none","low","moderate","high")
)

zoonotic_funding <- run_funding_scenarios(
  system_name       = "Direct-contact zoonoses",
  build_bn_fn       = build_zoonotic_bn,
  query_node        = "Direct_Transmission_Risk",
  query_state       = "present",
  ev_nodes          = c("Region_Context","Restoration_Intensity"),
  ev_states_template = c("High_Risk_Context","moderate"),
  restoration_states = c("none","low","moderate","high")
)

all_funding <- bind_rows(mosq_funding, tick_funding, rodent_funding, zoonotic_funding)

# =========================================================
# SUMMARY
# =========================================================

funding_summary <- all_funding %>%
  group_by(system, scenario, scenario_type, phase1, phase2, transition) %>%
  summarise(
    mean_risk  = mean(disease_risk),
    lower_risk = quantile(disease_risk, 0.025),
    upper_risk = quantile(disease_risk, 0.975),
    .groups    = "drop"
  )

print(funding_summary %>% select(system, scenario, mean_risk, lower_risk, upper_risk),
      n = 60)

# =========================================================
# COLOUR AND ORDERING
# =========================================================

# Separate colour palettes per system (consistent with main figures)
system_palettes <- list(
  "Mosquito-borne"         = c("#aaaaaa","#bdd7e7","#6baed6","#2171b5","#e7298a","#ce1256","#feb24c","#f46d43"),
  "Tick-borne"             = c("#aaaaaa","#bae4b3","#74c476","#238b45","#41ab5d","#006d2c","#feb24c","#f46d43"),
  "Rodent-borne"           = c("#aaaaaa","#fcae91","#fb6a4a","#cb181d","#ef3b2c","#99000d","#feb24c","#f46d43"),
  "Direct-contact zoonoses"= c("#aaaaaa","#cbc9e2","#9e9ac8","#6a51a3","#b30000","#4d0000","#feb24c","#f46d43")
)

# Scenario display order and line types
scenario_order <- c(
  "None (15 yr)", "Low (15 yr)", "Moderate (15 yr)", "High (15 yr)",
  "Low→High (yr 5)","Low→High (yr 10)",
  "Moderate→High (yr 5)","Moderate→High (yr 10)",
  "High→Low (yr 5)","High→Low (yr 10)",
  "High→Moderate (yr 5)","High→Moderate (yr 10)"
)

linetypes <- c(
  "None (15 yr)"          = "solid",
  "Low (15 yr)"           = "solid",
  "Moderate (15 yr)"      = "solid",
  "High (15 yr)"          = "solid",
  "Low→High (yr 5)"       = "dashed",
  "Low→High (yr 10)"      = "dotdash",
  "Moderate→High (yr 5)"  = "dashed",
  "Moderate→High (yr 10)" = "dotdash",
  "High→Low (yr 5)"       = "dashed",
  "High→Low (yr 10)"      = "dotdash",
  "High→Moderate (yr 5)"  = "dashed",
  "High→Moderate (yr 10)" = "dotdash"
)

# =========================================================
# PLOT FUNCTION
# =========================================================

make_funding_plot <- function(df_summary, system_name, pal, tag_label){
  
  d <- df_summary %>%
    filter(system == system_name) %>%
    mutate(scenario = factor(scenario, levels = scenario_order))
  
  # Separate pure baselines from transition scenarios
  d_pure  <- d %>% filter(scenario_type == "Pure baseline")
  d_trans <- d %>% filter(scenario_type != "Pure baseline")
  
  x_labels <- c("None","Low","Moderate","High",
                "Low→High\n(yr5)","Low→High\n(yr10)",
                "Mod→High\n(yr5)","Mod→High\n(yr10)",
                "High→Low\n(yr5)","High→Low\n(yr10)",
                "High→Mod\n(yr5)","High→Mod\n(yr10)")
  
  # Colour map: pure baselines get system colours, transitions get type colours
  col_map <- c(
    "None (15 yr)"          = pal[1],
    "Low (15 yr)"           = pal[2],
    "Moderate (15 yr)"      = pal[3],
    "High (15 yr)"          = pal[4],
    "Low→High (yr 5)"       = "#fde0dd",
    "Low→High (yr 10)"      = "#fa9fb5",
    "Moderate→High (yr 5)"  = "#f768a1",
    "Moderate→High (yr 10)" = "#c51b8a",
    "High→Low (yr 5)"       = "#ffffb2",
    "High→Low (yr 10)"      = "#fecc5c",
    "High→Moderate (yr 5)"  = "#fe9929",
    "High→Moderate (yr 10)" = "#d95f0e"
  )
  
  # Background shading for scenario groups
  shade_df <- data.frame(
    xmin = c(0.5, 4.5, 8.5),
    xmax = c(4.5, 8.5, 12.5),
    fill = c("grey97","#e8f4fd","#fff0f0")
  )
  
  ggplot(d, aes(x = scenario, y = mean_risk,
                colour = scenario, group = scenario)) +
    
    # Background group shading
    annotate("rect", xmin=0.5, xmax=4.5, ymin=-Inf, ymax=Inf,
             fill="grey96", alpha=0.5) +
    annotate("rect", xmin=4.5, xmax=8.5, ymin=-Inf, ymax=Inf,
             fill="#edf5fb", alpha=0.5) +
    annotate("rect", xmin=8.5, xmax=12.5, ymin=-Inf, ymax=Inf,
             fill="#fff0f0", alpha=0.5) +
    
    # Group labels
    annotate("text", x=2.5,  y=Inf, label="Pure baselines",
             vjust=1.5, size=2.8, colour="grey50", fontface="italic") +
    annotate("text", x=6.5,  y=Inf, label="Ramp-up scenarios",
             vjust=1.5, size=2.8, colour="#c51b8a", fontface="italic") +
    annotate("text", x=10.5, y=Inf, label="Funding collapse",
             vjust=1.5, size=2.8, colour="#feb24c", fontface="italic") +
    
    # Error bars
    geom_errorbar(
      aes(ymin=lower_risk, ymax=upper_risk, colour=scenario),
      width=0.3, linewidth=0.6
    ) +
    
    # Points
    geom_point(aes(fill=scenario), shape=21, size=3.5,
               colour="white", stroke=1.4) +
    
    scale_colour_manual(values=col_map, guide="none") +
    scale_fill_manual(values=col_map, guide="none") +
    
    scale_x_discrete(labels=x_labels) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      labels = scales::label_number(accuracy=0.1)
    ) +
    
    theme_minimal(base_size=10) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.x        = element_text(size=7.5, angle=30, hjust=1, colour="grey30"),
      axis.text.y        = element_text(size=8),
      axis.title         = element_text(size=8, colour="grey40"),
      plot.title         = element_text(size=10, face="bold"),
      plot.subtitle      = element_text(size=8,  colour="grey50", face="italic"),
      plot.tag           = element_text(size=9,  face="bold", colour="grey30")
    ) +
    labs(
      title    = system_name,
      subtitle = "Points = mean · bars = 95% UI · high-risk region context",
      x        = NULL,
      y        = "Disease risk probability",
      tag      = tag_label
    )
}

# =========================================================
# BUILD FOUR PANELS
# =========================================================

p1 <- make_funding_plot(funding_summary, "Mosquito-borne",
                        system_palettes[["Mosquito-borne"]], "A")
p2 <- make_funding_plot(funding_summary, "Tick-borne",
                        system_palettes[["Tick-borne"]], "B")
p3 <- make_funding_plot(funding_summary, "Rodent-borne",
                        system_palettes[["Rodent-borne"]], "C")
p4 <- make_funding_plot(funding_summary, "Direct-contact zoonoses",
                        system_palettes[["Direct-contact zoonoses"]], "D")

# =========================================================
# COMPOSE FIGURE
# =========================================================

fig_funding <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title   = "Supplementary Figure 3 — Funding scenario simulations: disease risk under sequential restoration transitions",
    subtitle = paste0(
      "Ramp-up: project begins at low or moderate intensity, transitions to high restoration at year 5 or 10. ",
      "Funding collapse: project begins at high intensity, drops to low or moderate at year 5 or 10. ",
      "Disease risk estimated using ecological legacy weighting at the transition point. ",
      "Pure baselines (grey background) shown for reference."
    ),
    theme = theme(
      plot.title    = element_text(size=11, face="bold"),
      plot.subtitle = element_text(size=8,  colour="grey50", margin=margin(b=10))
    )
  )


ggsave("supplementaryfig_funding_scenarios.png",
       fig_funding, width=260, height=200, units="mm", dpi=600)

# =========================================================
# SUMMARY TABLE — key comparisons for results text
# =========================================================

key_comparisons <- funding_summary %>%
  filter(scenario %in% c(
    "High (15 yr)",
    "Low→High (yr 5)","Low→High (yr 10)",
    "Moderate→High (yr 5)","Moderate→High (yr 10)",
    "High→Low (yr 5)","High→Low (yr 10)"
  )) %>%
  mutate(across(c(mean_risk, lower_risk, upper_risk), ~round(., 3))) %>%
  arrange(system, scenario)

print(key_comparisons, n=40)
write.csv(key_comparisons, "table_funding_scenario_key_comparisons.csv", row.names=FALSE)

# =========================================================
# FUNDING SCENARIO FIGURE — THREE SYSTEMS
# Mosquito, Tick, Rodent
#
# Requires: funding_summary from funding_scenario_simulations.R
#
# Visual design:
#   - Connected line + coloured points for pure baselines
#   - Dashed horizontal line = no-restoration reference
#   - Points with error bars for transition scenarios
#   - Shaded ribbon = 95% UI envelope across scenario group
#   - Filled circle = year 5 transition
#   - Open circle  = year 10 transition
#   - purple tones   = ramp-up scenarios
#   - grey tones    = funding collapse scenarios
# =========================================================

# =========================================================
# COLOUR PALETTES
# =========================================================

mosquito_cols <- c(
  none     = "#aaaaaa",
  low      = "#bdd7e7",
  moderate = "#6baed6",
  high     = "#2171b5"
)

tick_cols <- c(
  none     = "#aaaaaa",
  low      = "#bae4b3",
  moderate = "#74c476",
  high     = "#238b45"
)

rodent_cols <- c(
  none     = "#aaaaaa",
  low      = "#fcae91",
  moderate = "#fb6a4a",
  high     = "#cb181d"
)

# =========================================================
# X-AXIS POSITION MAP
# Pure baselines sit at integer positions 0-3.
# Transition scenario points are offset slightly so they
# cluster visually around the restoration level they relate
# to, without overlapping the baseline points.
#
# Ramp-up scenarios (Low→High, Moderate→High):
#   Low→High   yr5  = 1.15,  yr10 = 1.35
#   Mod→High   yr5  = 2.15,  yr10 = 2.35
#
# Collapse scenarios (High→Low, High→Moderate):
#   High→Low   yr5  = 2.65,  yr10 = 2.85
#   High→Mod   yr5  = 1.65,  yr10 = 1.85
# =========================================================

x_positions <- c(
  "none"     = 0,
  "low"      = 1,
  "moderate" = 2,
  "high"     = 3
)

transition_x <- c(
  "Low→High (yr 5)"       = 1.15,
  "Low→High (yr 10)"      = 1.35,
  "Moderate→High (yr 5)"  = 2.15,
  "Moderate→High (yr 10)" = 2.35,
  "High→Low (yr 5)"       = 2.65,
  "High→Low (yr 10)"      = 2.85,
  "High→Moderate (yr 5)"  = 1.65,
  "High→Moderate (yr 10)" = 1.85
)

# =========================================================
# SHARED THEME
# =========================================================

theme_funding <- function(){
  theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x        = element_text(size = 9,  colour = "grey35"),
      axis.text.y        = element_text(size = 9,  colour = "grey35"),
      axis.title.x       = element_text(size = 9,  colour = "grey40",
                                        margin = margin(t = 6)),
      axis.title.y       = element_text(size = 9,  colour = "grey40",
                                        margin = margin(r = 6)),
      plot.title         = element_text(size = 11, face = "bold"),
      plot.subtitle      = element_text(size = 8,  colour = "grey50",
                                        face = "italic",
                                        margin = margin(b = 6)),
      plot.tag           = element_text(size = 10, face = "bold",
                                        colour = "grey25"),
      legend.position    = "none"
    )
}

# =========================================================
# PANEL BUILDER FUNCTION
# =========================================================

make_funding_panel <- function(sys_name,
                               sys_title,
                               intensity_cols,
                               tag_label,
                               subtitle = NULL){
  
  # ── Split funding_summary into baseline and transition ──
  
  df_base <- funding_summary %>%
    filter(
      system        == sys_name,
      scenario_type == "Pure baseline"
    ) %>%
    mutate(
      x     = x_positions[phase1],
      scen  = factor(phase1, levels = names(x_positions))
    ) %>%
    arrange(x)
  
  df_trans <- funding_summary %>%
    filter(
      system        == sys_name,
      scenario_type != "Pure baseline"
    ) %>%
    mutate(
      x            = transition_x[scenario],
      transition_yr = factor(transition),
      col_group     = scenario_type   # "Ramp-up" or "Funding collapse"
    )
  
  # No-restoration baseline value (dashed reference line)
  no_rest_val <- df_base$mean_risk[df_base$phase1 == "none"]
  
  # ── Ribbon envelope: min lower / max upper per scenario type ──
  ribbon_df <- df_trans %>%
    group_by(col_group) %>%
    summarise(
      x_min = min(x) - 0.08,
      x_max = max(x) + 0.08,
      y_min = min(lower_risk),
      y_max = max(upper_risk),
      .groups = "drop"
    )
  
  # ── Transition point colours and shapes ──
  rampup_col   <- "#c51b8a"
  collapse_col <- "#feb24c"
  
  df_trans <- df_trans %>%
    mutate(
      pt_colour = case_when(
        col_group == "Ramp-up"          ~ rampup_col,
        col_group == "Funding collapse"  ~ collapse_col
      ),
      pt_fill = case_when(
        col_group == "Ramp-up"   & transition == 5  ~ rampup_col,
        col_group == "Ramp-up"   & transition == 10 ~ "white",
        col_group == "Funding collapse" & transition == 5  ~ collapse_col,
        col_group == "Funding collapse" & transition == 10 ~ "white"
      )
    )
  
  # ── Build plot ────────────────────────────────────────────
  
  ggplot() +
    
    # No-restoration dashed reference line
    geom_hline(
      yintercept = no_rest_val,
      linetype   = "dashed",
      colour     = "grey55",
      linewidth  = 0.55
    ) +
    
    # Shaded ribbon: 95% UI envelope for each scenario group
    geom_rect(
      data = ribbon_df,
      aes(xmin = x_min, xmax = x_max,
          ymin = y_min, ymax = y_max,
          fill = col_group),
      alpha = 0.10
    ) +
    scale_fill_manual(
      values = c(
        "Ramp-up"          = rampup_col,
        "Funding collapse"  = collapse_col
      )
    ) +
    
    # Pure baseline connected line
    geom_line(
      data    = df_base,
      aes(x = x, y = mean_risk),
      colour  = "grey30",
      linewidth = 1.4
    ) +
    
    # Pure baseline 95% UI ribbon
    geom_ribbon(
      data  = df_base,
      aes(x = x, ymin = lower_risk, ymax = upper_risk),
      fill  = "grey70",
      alpha = 0.18
    ) +
    
    # Pure baseline points — coloured by restoration intensity
    geom_point(
      data  = df_base,
      aes(x = x, y = mean_risk, fill = phase1),
      shape = 21,
      size  = 4.5,
      colour = "white",
      stroke = 1.5
    ) +
    scale_fill_manual(
      values = c(intensity_cols,
                 "Ramp-up"          = rampup_col,
                 "Funding collapse"  = collapse_col)
    ) +
    
    # Transition scenario error bars
    geom_errorbar(
      data = df_trans,
      aes(x    = x,
          ymin = lower_risk,
          ymax = upper_risk,
          colour = col_group),
      width     = 0.04,
      linewidth = 0.75
    ) +
    scale_colour_manual(
      values = c(
        "Ramp-up"          = rampup_col,
        "Funding collapse"  = collapse_col
      )
    ) +
    
    # Transition scenario points
    # Filled circle = year 5 transition
    # Open circle   = year 10 transition
    geom_point(
      data  = df_trans,
      aes(x = x, y = mean_risk),
      shape  = 21,
      size   = 3.8,
      colour = df_trans$pt_colour,
      fill   = df_trans$pt_fill,
      stroke = 1.8
    ) +
    
    # Axis scales
    scale_x_continuous(
      breaks = 0:3,
      labels = c("None","Low","Moderate","High"),
      limits = c(-0.35, 3.35)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      labels = scales::label_number(accuracy = 0.1),
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    
    theme_funding() +
    
    labs(
      title    = sys_title,
      subtitle = subtitle,
      x        = "Restoration intensity",
      y        = "Disease risk probability",
      tag      = tag_label
    )
}

# =========================================================
# BUILD THREE PANELS
# =========================================================

p_mosq <- make_funding_panel(
  sys_name      = "Mosquito-borne",
  sys_title     = "Mosquito-borne",
  intensity_cols = mosquito_cols,
  tag_label     = "A",
  subtitle      = "Peak risk at moderate restoration · Hydrological mechanism"
)

p_tick <- make_funding_panel(
  sys_name      = "Tick-borne",
  sys_title     = "Tick-borne",
  intensity_cols = tick_cols,
  tag_label     = "B",
  subtitle      = "Risk hump at low restoration · Trophic lag mechanism"
)

p_rod <- make_funding_panel(
  sys_name      = "Rodent-borne",
  sys_title     = "Rodent-borne",
  intensity_cols = rodent_cols,
  tag_label     = "C",
  subtitle      = "Monotonic decline · Habitat-mediated rodent suppression"
)

# =========================================================
# SHARED LEGEND (built manually as a plot)
# =========================================================

legend_data <- data.frame(
  x     = c(1,2,3,4,5,6,7,8,9,10),
  y     = rep(1, 10),
  group = c("none","low","moderate","high",
            "ru5","ru10","cl5","cl10",
            "baseline_line","no_rest_ref")
)

legend_plot <- ggplot() +
  
  # Baseline line segment
  annotate("segment", x=0.5, xend=1.5, y=3, yend=3,
           colour="grey30", linewidth=1.2) +
  annotate("point",   x=1.0, y=3, shape=21, size=3.5,
           colour="white", fill="#2171b5", stroke=1.4) +
  annotate("text",    x=1.7, y=3, label="Sustained baseline (line + coloured points)",
           hjust=0, size=3.2, colour="grey30") +
  
  # Dashed reference
  annotate("segment", x=0.5, xend=1.5, y=2.3, yend=2.3,
           colour="grey55", linewidth=0.8, linetype="dashed") +
  annotate("text",    x=1.7, y=2.3, label="No-restoration reference (dashed)",
           hjust=0, size=3.2, colour="grey30") +
  
  # Ramp-up yr5
  annotate("point", x=0.9, y=1.55, shape=21, size=3.5,
           colour="#c51b8a", fill="#c51b8a", stroke=1.8) +
  annotate("text",  x=1.7, y=1.7, label="Ramp-up scenario (filled = yr 5)",
           hjust=0, size=3.2, colour="grey30") +
  
  # Ramp-up yr10
  annotate("point", x=0.9, y=1.1, shape=21, size=3.5,
           colour="#c51b8a", fill="white", stroke=1.8) +
  annotate("text",  x=1.7, y=1.1, label="Ramp-up scenario (open = yr 10)",
           hjust=0, size=3.2, colour="grey30") +
  
  # Collapse yr5
  annotate("point", x=4.0, y=1.7, shape=21, size=3.5,
           colour="#feb24c", fill="#feb24c", stroke=1.8) +
  annotate("text",  x=4.2, y=1.7, label="Funding collapse (filled = yr 5)",
           hjust=0, size=3.2, colour="grey30") +
  
  # Collapse yr10
  annotate("point", x=4.0, y=1.1, shape=21, size=3.5,
           colour="#feb24c", fill="white", stroke=1.8) +
  annotate("text",  x=4.2, y=1.1, label="Funding collapse (open = yr 10)",
           hjust=0, size=3.2, colour="grey30") +
  
  # Shaded ribbon swatches
  annotate("rect", xmin=0.5, xmax=1.5, ymin=0.55, ymax=0.85,
           fill="#c51b8a", alpha=0.15) +
  annotate("text", x=1.7, y=0.70,
           label="Shaded band = 95% UI envelope across scenario group",
           hjust=0, size=3.2, colour="grey30") +
  
  xlim(0, 9) + ylim(0, 3.6) +
  theme_void()

# =========================================================
# COMPOSE FIGURE
# =========================================================

fig_funding <- (p_mosq | p_tick | p_rod) /
  legend_plot +
  plot_layout(heights = c(5, 1)) +
  plot_annotation(
    caption = paste0(
      "Lines = mean disease risk for sustained restoration scenarios (pure baselines). ",
      "Points = mean risk under funding transition scenarios. ",
      "Error bars = 95% Monte Carlo uncertainty interval (n = 1,000 iterations). ",
      "All scenarios conditioned on high-risk region context. ",
      "Transition points positioned horizontally relative to the restoration ",
      "intensity level of the originating phase."
    ),
    theme = theme(
      plot.caption = element_text(size = 7.5, colour = "grey50",
                                  hjust = 0, margin = margin(t = 8))
    )
  )

# =========================================================
# SAVE
# =========================================================

ggsave(
  "fig_funding_scenarios_three_systems.pdf",
  plot   = fig_funding,
  width  = 240,
  height = 160,
  units  = "mm",
  device = cairo_pdf
)

ggsave(
  "fig_funding_scenarios_three_systems.png",
  plot   = fig_funding,
  width  = 300,
  height = 200,
  units  = "mm",
  dpi    = 600
)

