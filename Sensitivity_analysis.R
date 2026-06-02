# =========================================================
# SENSITIVITY ANALYSIS + TORNADO PLOTS  v2
# FOUR ZOONOTIC DISEASE SYSTEMS
#
# Method: One-At-a-Time (OAT) — FULL DISTRIBUTIONAL SHIFT
#
#
# PARAMETER DESIGN:
#   Each param is a named list with:
#     $name    : label for plot
#     $node    : which CPT to modify ("sw", "dis", etc.)
#     $col_idx : which column (scenario) within the CPT matrix
#     $lo_probs: full probability vector at lower-bound assumption
#     $hi_probs: full probability vector at upper-bound assumption
#
# =========================================================

library(bnlearn)
library(gRain)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(stringr)

set.seed(42)

# =========================================================
# SHARED HELPERS
# =========================================================

mid        <- function(mn, mx) (mn + mx) / 2
norm       <- function(x){ x[x < 0.01] <- 0.01; x / sum(x) }
norm_cols  <- function(mat){
  # normalise each column of a matrix
  apply(mat, 2, norm)
}

query_risk <- function(bn_grain, risk_node, risk_state,
                       ev_nodes, ev_states){
  bn_ev <- setEvidence(bn_grain,
                       nodes  = ev_nodes,
                       states = ev_states)
  querygrain(bn_ev, nodes = risk_node)[[risk_node]][risk_state]
}

# Replace one column of a matrix with a new probability vector
# (used to inject perturbations into CPT matrices)
replace_col <- function(mat, col_idx, new_probs){
  mat[, col_idx] <- norm(new_probs)
  mat
}

# =========================================================
# SYSTEM 1: MOSQUITO-BORNE
# =========================================================

run_mosquito_sensitivity <- function(){
  
  states <- list(
    Region_Context        = c("Low","Moderate","High"),
    Restoration_Intensity = c("none","low","moderate","high"),
    Hydrological_Recovery = c("low","moderate","high"),
    Standing_Water        = c("low","moderate","high"),
    Mosquito_Disease_Risk = c("absent","present"),
    Ecological_Regulation = c("weak","moderate","strong"),
    Time_to_Stable_State  = c("1_3_years","4_7_years",
                              "8_12_years","13_15_years")
  )
  
  dag <- model2network(paste0(
    "[Region_Context]",
    "[Restoration_Intensity]",
    "[Hydrological_Recovery|Region_Context:Restoration_Intensity]",
    "[Standing_Water|Hydrological_Recovery:Restoration_Intensity]",
    "[Ecological_Regulation|Restoration_Intensity]",
    "[Mosquito_Disease_Risk|Standing_Water:Ecological_Regulation]",
    "[Time_to_Stable_State|Ecological_Regulation]"
  ))
  
  # --------------------------------------------------------
  # BASELINE CPT MATRICES (columns = conditioning scenarios)
  # Each column sums to 1 after norm().
  # --------------------------------------------------------
  
  # Standing water: dim=c(3,12) — 3 SW states x (4 rest * 3 hydro)
  # Column order (rest outer, hydro inner):
  #   col 1-3:  rest=none,  hydro=low/mod/high
  #   col 4-6:  rest=low,   hydro=low/mod/high
  #   col 7-9:  rest=mod,   hydro=low/mod/high
  #   col 10-12:rest=high,  hydro=low/mod/high
  sw_base <- matrix(c(
    # rest=none
    mid(.25,.40)+.10, mid(.40,.55), mid(.10,.20),   # hydro=low  (+.10 low-water adj)
    mid(.25,.40),     mid(.40,.55), mid(.10,.20),   # hydro=mod
    mid(.25,.40),     mid(.40,.55), mid(.10,.20)+.08,# hydro=high (+.08 high-water adj)
    # rest=low
    mid(.10,.20)+.10, mid(.35,.50), mid(.35,.50),
    mid(.10,.20),     mid(.35,.50), mid(.35,.50),
    mid(.10,.20),     mid(.35,.50), mid(.35,.50)+.08,
    # rest=moderate (PEAK standing water)
    mid(.02,.08)+.10, mid(.15,.28), mid(.65,.82),
    mid(.02,.08),     mid(.15,.28), mid(.65,.82),
    mid(.02,.08),     mid(.15,.28), mid(.65,.82)+.08,
    # rest=high (suppressed)
    mid(.35,.55)+.10, mid(.30,.45), mid(.05,.15),
    mid(.35,.55),     mid(.30,.45), mid(.05,.15),
    mid(.35,.55),     mid(.30,.45), mid(.05,.15)+.08
  ), nrow=3, ncol=12)
  sw_base <- norm_cols(sw_base)
  
  # Ecological regulation: dim=c(3,4) — cols = rest scenarios
  reg_base <- matrix(c(
    mid(.78,.84), mid(.13,.15), mid(.03,.05),  # none
    mid(.58,.64), mid(.24,.28), mid(.09,.11),  # low
    mid(.25,.27), mid(.42,.44), mid(.29,.31),  # moderate
    mid(.05,.07), mid(.15,.17), mid(.72,.78)   # high
  ), nrow=3, ncol=4)
  reg_base <- norm_cols(reg_base)
  
  # Disease: dim=c(2,9) — 9 cols = (3 SW) x (3 Reg)
  # Col order: reg outer, SW inner
  # reg=weak:   SW=low(1), SW=mod(2), SW=high(3)
  # reg=mod:    SW=low(4), SW=mod(5), SW=high(6)
  # reg=strong: SW=low(7), SW=mod(8), SW=high(9)
  dis_base <- matrix(c(
    mid(.50,.65), mid(.35,.50),   # reg=weak,   SW=low
    mid(.30,.45), mid(.55,.70),   # reg=weak,   SW=mod
    mid(.05,.12), mid(.88,.95),   # reg=weak,   SW=high  ← PEAK
    mid(.70,.82), mid(.18,.30),   # reg=mod,    SW=low
    mid(.45,.65), mid(.35,.55),   # reg=mod,    SW=mod
    mid(.25,.40), mid(.60,.75),   # reg=mod,    SW=high
    mid(.88,.96), mid(.04,.12),   # reg=strong, SW=low
    mid(.78,.90), mid(.10,.22),   # reg=strong, SW=mod
    mid(.82,.92), mid(.08,.18)    # reg=strong, SW=high  ← KEY: strong reg suppresses
  ), nrow=2, ncol=9, byrow=TRUE)
  dis_base <- norm_cols(dis_base)
  
  # Time: dim=c(4,3) — cols=regulation states
  time_base <- matrix(c(
    mid(.38,.52), mid(.28,.38), mid(.08,.18), mid(.03,.10),  # weak reg
    mid(.14,.24), mid(.32,.44), mid(.22,.34), mid(.10,.22),  # mod reg
    mid(.03,.09), mid(.10,.20), mid(.42,.56), mid(.24,.38)   # strong reg
  ), nrow=4, ncol=3)
  time_base <- norm_cols(time_base)
  
  # Hydrological recovery: dim=c(3,12) — same layout as SW
  # (averaged across regions for sensitivity; region effect is secondary)
  hyd_base <- matrix(c(
    mid(.60,.80), mid(.15,.30), mid(.05,.15),  # none,  hydro states per rest
    mid(.60,.80), mid(.15,.30), mid(.05,.15),
    mid(.60,.80), mid(.15,.30), mid(.05,.15),
    mid(.20,.35), mid(.35,.50), mid(.25,.40),  # low
    mid(.20,.35), mid(.35,.50), mid(.25,.40),
    mid(.20,.35), mid(.35,.50), mid(.25,.40),
    mid(.05,.20), mid(.25,.40), mid(.45,.65),  # moderate
    mid(.05,.20), mid(.25,.40), mid(.45,.65),
    mid(.05,.20), mid(.25,.40), mid(.45,.65),
    mid(.05,.15), mid(.20,.35), mid(.55,.75),  # high
    mid(.05,.15), mid(.20,.35), mid(.55,.75),
    mid(.05,.15), mid(.20,.35), mid(.55,.75)
  ), nrow=3, ncol=12)
  hyd_base <- norm_cols(hyd_base)
  
  # --------------------------------------------------------
  # SENSITIVITY PARAMETERS
  # Each entry: which CPT column to perturb, and the full
  # probability vector at lower-bound vs upper-bound assumption.
  # --------------------------------------------------------
  
  # Helper: make a 3-element prob vector from pessimistic/optimistic
  # framing of each scenario
  params <- list(
    
    # --- STANDING WATER (the hump mechanism) ---
    
    list(name = "Standing water: none rest (baseline level)",
         node = "sw", col = 2,  # rest=none, hydro=mod (representative)
         lo   = norm(c(.45, .38, .08)),   # optimistic: less standing water at baseline
         hi   = norm(c(.18, .42, .30))),  # pessimistic: more standing water at baseline
    
    list(name = "Standing water: low rest → peak SW",
         node = "sw", col = 5,  # rest=low, hydro=mod
         lo   = norm(c(.25, .42, .25)),   # hump smaller: low rest doesn't raise SW much
         hi   = norm(c(.04, .22, .72))),  # hump larger: low rest sharply raises SW
    
    list(name = "Standing water: moderate rest → peak SW",
         node = "sw", col = 8,  # rest=mod, hydro=mod
         lo   = norm(c(.20, .35, .38)),   # hump smaller at moderate rest
         hi   = norm(c(.01, .10, .85))),  # hump peak higher at moderate rest
    
    list(name = "Standing water: high rest → low SW (recovery)",
         node = "sw", col = 11, # rest=high, hydro=mod
         lo   = norm(c(.52, .36, .06)),   # strong recovery: SW well suppressed
         hi   = norm(c(.22, .42, .28))),  # weak recovery: SW not well suppressed
    
    # --- ECOLOGICAL REGULATION ---
    
    list(name = "Regulation: none rest → weak reg",
         node = "reg", col = 1,
         lo   = norm(c(.68, .22, .08)),   # slightly better baseline regulation
         hi   = norm(c(.88, .09, .02))),  # near-total regulatory failure
    
    list(name = "Regulation: high rest → strong reg",
         node = "reg", col = 4,
         lo   = norm(c(.12, .22, .58)),   # weaker regulatory recovery
         hi   = norm(c(.02, .10, .86))),  # stronger regulatory recovery
    
    list(name = "Regulation: moderate rest (transition)",
         node = "reg", col = 3,
         lo   = norm(c(.18, .42, .36)),   # faster regulation at moderate rest
         hi   = norm(c(.36, .44, .18))),  # slower regulation at moderate rest
    
    # --- DISEASE CPT ---
    
    list(name = "Disease: high SW + weak reg → risk (hump peak)",
         node = "dis", col = 3,  # reg=weak, SW=high
         lo   = norm(c(.18, .78)),   # lower peak risk
         hi   = norm(c(.03, .97))),  # higher peak risk
    
    list(name = "Disease: high SW + strong reg → risk (recovery)",
         node = "dis", col = 9,  # reg=strong, SW=high
         lo   = norm(c(.92, .05)),   # strong suppression by regulation
         hi   = norm(c(.65, .32))),  # weak suppression — risk persists
    
    list(name = "Disease: mod SW + mod reg → risk (transition)",
         node = "dis", col = 5,  # reg=mod, SW=mod
         lo   = norm(c(.72, .26)),   # moderate reg effective
         hi   = norm(c(.32, .66))),  # moderate reg ineffective
    
    # --- HYDROLOGICAL RECOVERY ---
    
    list(name = "Hydro recovery: none rest → low recovery",
         node = "hyd", col = 2,  # rest=none, hydro=mod (representative)
         lo   = norm(c(.52, .30, .14)),   # some natural hydro recovery
         hi   = norm(c(.78, .16, .04))),  # near-complete hydro degradation
    
    list(name = "Hydro recovery: high rest → full recovery",
         node = "hyd", col = 11, # rest=high, hydro=mod
         lo   = norm(c(.08, .26, .62)),   # slower hydro recovery
         hi   = norm(c(.02, .14, .82)))   # faster hydro recovery
  )
  
  # --------------------------------------------------------
  # BUILD BN FROM CPT MATRICES
  # --------------------------------------------------------
  
  build_bn <- function(sw_mat  = sw_base,
                       reg_mat = reg_base,
                       dis_mat = dis_base,
                       hyd_mat = hyd_base){
    
    cpt_region <- array(c(0.33,0.34,0.33), dim=c(3),
                        dimnames=list(Region_Context=states$Region_Context))
    cpt_rest   <- array(rep(0.25,4), dim=c(4),
                        dimnames=list(Restoration_Intensity=states$Restoration_Intensity))
    
    # Hydrological recovery: dim=c(3,4,3)=[Hydro,Rest,Region]
    # Replicate hyd_mat (12 cols) across 3 regions
    hyd_vec <- as.vector(hyd_mat)
    hyd_full <- rep(hyd_vec, 3)
    cpt_hyd  <- array(hyd_full, dim=c(3,4,3),
                      dimnames=list(
                        Hydrological_Recovery = states$Hydrological_Recovery,
                        Restoration_Intensity = states$Restoration_Intensity,
                        Region_Context        = states$Region_Context))
    
    # Standing water: dim=c(3,3,4)=[SW,Hydro,Rest]
    cpt_sw <- array(as.vector(sw_mat), dim=c(3,3,4),
                    dimnames=list(
                      Standing_Water        = states$Standing_Water,
                      Hydrological_Recovery = states$Hydrological_Recovery,
                      Restoration_Intensity = states$Restoration_Intensity))
    
    # Ecological regulation: dim=c(3,4)=[Reg,Rest]
    cpt_reg <- array(as.vector(reg_mat), dim=c(3,4),
                     dimnames=list(
                       Ecological_Regulation = states$Ecological_Regulation,
                       Restoration_Intensity = states$Restoration_Intensity))
    
    # Disease: dim=c(2,3,3)=[Disease,SW,Reg]
    cpt_dis <- array(as.vector(dis_mat), dim=c(2,3,3),
                     dimnames=list(
                       Mosquito_Disease_Risk  = states$Mosquito_Disease_Risk,
                       Standing_Water        = states$Standing_Water,
                       Ecological_Regulation = states$Ecological_Regulation))
    
    # Time: dim=c(4,3)=[Time,Reg]
    cpt_time <- array(as.vector(time_base), dim=c(4,3),
                      dimnames=list(
                        Time_to_Stable_State  = states$Time_to_Stable_State,
                        Ecological_Regulation = states$Ecological_Regulation))
    
    fit <- custom.fit(dag, dist=list(
      Region_Context        = cpt_region,
      Restoration_Intensity = cpt_rest,
      Hydrological_Recovery = cpt_hyd,
      Standing_Water        = cpt_sw,
      Ecological_Regulation = cpt_reg,
      Mosquito_Disease_Risk = cpt_dis,
      Time_to_Stable_State  = cpt_time
    ))
    as.grain(fit)
  }
  
  bn_base  <- build_bn()
  baseline <- query_risk(bn_base, "Mosquito_Disease_Risk", "present",
                         c("Region_Context","Restoration_Intensity"),
                         c("High","moderate"))
  
  # --------------------------------------------------------
  # OAT LOOP — perturb full column, one at a time
  # --------------------------------------------------------
  
  results <- bind_rows(lapply(params, function(p){
    
    # Identify which matrix to modify
    mat_lo <- mat_hi <- switch(p$node,
                               "sw"  = sw_base,
                               "reg" = reg_base,
                               "dis" = dis_base,
                               "hyd" = hyd_base
    )
    mat_lo[, p$col] <- p$lo
    mat_hi[, p$col] <- p$hi
    
    args_lo <- setNames(list(mat_lo), paste0(p$node, "_mat"))
    args_hi <- setNames(list(mat_hi), paste0(p$node, "_mat"))
    
    bn_lo <- do.call(build_bn, args_lo)
    bn_hi <- do.call(build_bn, args_hi)
    
    risk_lo <- query_risk(bn_lo, "Mosquito_Disease_Risk","present",
                          c("Region_Context","Restoration_Intensity"),
                          c("High","moderate"))
    risk_hi <- query_risk(bn_hi, "Mosquito_Disease_Risk","present",
                          c("Region_Context","Restoration_Intensity"),
                          c("High","moderate"))
    
    tibble(parameter = p$name,
           risk_low  = risk_lo,
           risk_high = risk_hi,
           swing     = risk_hi - risk_lo,
           abs_swing = abs(risk_hi - risk_lo),
           baseline  = baseline,
           system    = "Mosquito-borne")
  }))
  
  results
}

# =========================================================
# SYSTEM 2: TICK-BORNE
# =========================================================

run_tick_sensitivity <- function(){
  
  states <- list(
    Region_Context        = c("Low","Moderate","High"),
    Restoration_Intensity = c("none","low","moderate","high"),
    Host_Connectivity     = c("low","moderate","high"),
    Predator_Recovery     = c("weak","moderate","strong"),
    Tick_Disease_Risk     = c("absent","present"),
    Ecological_Regulation = c("weak","moderate","strong"),
    Time_to_Stable_State  = c("1_3_years","4_7_years",
                              "8_12_years","13_15_years")
  )
  
  dag <- model2network(paste0(
    "[Region_Context]",
    "[Restoration_Intensity]",
    "[Host_Connectivity|Restoration_Intensity]",
    "[Predator_Recovery|Restoration_Intensity]",
    "[Tick_Disease_Risk|Host_Connectivity:Predator_Recovery]",
    "[Ecological_Regulation|Predator_Recovery:Tick_Disease_Risk]",
    "[Time_to_Stable_State|Ecological_Regulation]"
  ))
  
  # Host connectivity: dim=c(3,4)=[HC,Rest]  cols=rest scenarios
  con_base <- matrix(c(
    mid(.55,.72), mid(.18,.30), mid(.05,.15),  # none
    mid(.05,.12), mid(.18,.30), mid(.60,.78),  # low  ← PEAK
    mid(.08,.18), mid(.22,.35), mid(.50,.68),  # moderate
    mid(.15,.30), mid(.38,.52), mid(.25,.42)   # high
  ), nrow=3, ncol=4)
  con_base <- norm_cols(con_base)
  
  # Predator recovery: dim=c(3,4)=[Pred,Rest]  cols=rest scenarios
  pred_base <- matrix(c(
    mid(.78,.86), mid(.10,.18), mid(.02,.06),  # none
    mid(.68,.76), mid(.18,.26), mid(.04,.08),  # low  ← lag persists
    mid(.48,.56), mid(.29,.37), mid(.11,.19),  # moderate
    mid(.08,.16), mid(.24,.32), mid(.56,.64)   # high
  ), nrow=3, ncol=4)
  pred_base <- norm_cols(pred_base)
  
  # Disease: dim=c(2,9)=[Disease,9 combos]
  # Col order: pred outer (weak/mod/strong), HC inner (low/mod/high)
  dis_base <- matrix(c(
    mid(.50,.65), mid(.35,.50),  # weak pred, low HC
    mid(.30,.55), mid(.45,.70),  # weak pred, mod HC
    mid(.05,.14), mid(.86,.95),  # weak pred, high HC  ← PEAK
    mid(.28,.68), mid(.32,.72),  # mod pred,  low HC
    mid(.30,.70), mid(.30,.70),  # mod pred,  mod HC
    mid(.28,.68), mid(.32,.72),  # mod pred,  high HC
    mid(.75,.90), mid(.10,.25),  # strong pred, low HC
    mid(.28,.68), mid(.32,.72),  # strong pred, mod HC
    mid(.68,.82), mid(.18,.32)   # strong pred, high HC ← suppressed
  ), nrow=2, ncol=9, byrow=TRUE)
  dis_base <- norm_cols(dis_base)
  
  # Regulation: dim=c(3,6)=[Reg, pred x disease combos]
  # Col order: disease outer, pred inner
  reg_base <- matrix(c(
    mid(.65,.80), mid(.15,.26), mid(.04,.10),  # dis=absent, pred=weak
    mid(.18,.48), mid(.22,.48), mid(.18,.48),  # dis=absent, pred=mod
    mid(.04,.10), mid(.15,.26), mid(.66,.80),  # dis=absent, pred=strong ← KEY
    mid(.65,.80), mid(.15,.26), mid(.04,.10),  # dis=present, pred=weak
    mid(.18,.48), mid(.22,.48), mid(.18,.48),  # dis=present, pred=mod
    mid(.18,.48), mid(.22,.48), mid(.18,.48)   # dis=present, pred=strong
  ), nrow=3, ncol=6, byrow=TRUE)
  reg_base <- norm_cols(reg_base)
  
  time_base <- matrix(c(
    mid(.38,.52), mid(.28,.38), mid(.08,.18), mid(.03,.10),
    mid(.15,.28), mid(.22,.32), mid(.24,.34), mid(.20,.32),
    mid(.02,.06), mid(.05,.12), mid(.22,.34), mid(.52,.68)
  ), nrow=4, ncol=3)
  time_base <- norm_cols(time_base)
  
  params <- list(
    
    list(name = "Host connectivity: low rest → peak connectivity",
         node = "con", col = 2,
         lo   = norm(c(.18, .35, .38)),   # small hump: low rest ≈ baseline
         hi   = norm(c(.02, .12, .84))),  # large hump: low rest drives very high HC
    
    list(name = "Host connectivity: none rest (baseline)",
         node = "con", col = 1,
         lo   = norm(c(.70, .20, .06)),   # very low baseline connectivity
         hi   = norm(c(.35, .35, .26))),  # moderate baseline connectivity
    
    list(name = "Host connectivity: high rest → suppressed",
         node = "con", col = 4,
         lo   = norm(c(.42, .42, .12)),   # poor recovery of connectivity control
         hi   = norm(c(.12, .32, .50))),  # wait — high rest should REDUCE HC via trophic
    
    list(name = "Predator recovery: low rest → weak (lag strength)",
         node = "pred", col = 2,
         lo   = norm(c(.55, .32, .10)),   # predators recover earlier at low rest
         hi   = norm(c(.86, .11, .02))),  # strong lag: predators almost absent
    
    list(name = "Predator recovery: none rest → weak",
         node = "pred", col = 1,
         lo   = norm(c(.65, .25, .08)),   # some background predator presence
         hi   = norm(c(.92, .06, .01))),  # near-total predator absence
    
    list(name = "Predator recovery: high rest → strong",
         node = "pred", col = 4,
         lo   = norm(c(.22, .36, .38)),   # slow trophic recovery
         hi   = norm(c(.04, .16, .78))),  # fast trophic recovery
    
    list(name = "Disease: high HC + weak pred → present (hump peak)",
         node = "dis", col = 3,
         lo   = norm(c(.20, .76)),        # lower peak risk
         hi   = norm(c(.02, .97))),       # higher peak risk
    
    list(name = "Disease: high HC + strong pred → present (recovery)",
         node = "dis", col = 9,
         lo   = norm(c(.88, .10)),        # strong predator suppression
         hi   = norm(c(.52, .46))),       # weak predator suppression
    
    list(name = "Regulation: strong pred + absent disease → strong",
         node = "reg", col = 3,
         lo   = norm(c(.02, .12, .84)),   # efficient trophic regulation
         hi   = norm(c(.18, .32, .46))),  # weak trophic regulation
    
    list(name = "Regulation: weak pred + present disease → weak",
         node = "reg", col = 4,
         lo   = norm(c(.52, .28, .16)),   # moderate regulatory failure
         hi   = norm(c(.84, .12, .03))),  # near-complete regulatory failure
    
    list(name = "Predator recovery: moderate rest (transition speed)",
         node = "pred", col = 3,
         lo   = norm(c(.32, .42, .22)),   # faster predator recovery at mod rest
         hi   = norm(c(.62, .28, .08))),  # slower predator recovery at mod rest
    
    list(name = "Disease: low HC + weak pred (no-rest baseline risk)",
         node = "dis", col = 1,
         lo   = norm(c(.68, .30)),        # lower no-restoration baseline
         hi   = norm(c(.38, .60)))        # higher no-restoration baseline
  )
  
  build_bn <- function(con_mat  = con_base,
                       pred_mat = pred_base,
                       dis_mat  = dis_base,
                       reg_mat  = reg_base){
    
    cpt_region <- array(c(0.33,0.34,0.33), dim=c(3),
                        dimnames=list(Region_Context=states$Region_Context))
    cpt_rest   <- array(rep(0.25,4), dim=c(4),
                        dimnames=list(Restoration_Intensity=states$Restoration_Intensity))
    
    cpt_con <- array(as.vector(con_mat), dim=c(3,4),
                     dimnames=list(Host_Connectivity=states$Host_Connectivity,
                                   Restoration_Intensity=states$Restoration_Intensity))
    cpt_pred <- array(as.vector(pred_mat), dim=c(3,4),
                      dimnames=list(Predator_Recovery=states$Predator_Recovery,
                                    Restoration_Intensity=states$Restoration_Intensity))
    cpt_dis <- array(as.vector(dis_mat), dim=c(2,3,3),
                     dimnames=list(Tick_Disease_Risk=states$Tick_Disease_Risk,
                                   Host_Connectivity=states$Host_Connectivity,
                                   Predator_Recovery=states$Predator_Recovery))
    cpt_reg <- array(as.vector(reg_mat), dim=c(3,3,2),
                     dimnames=list(Ecological_Regulation=states$Ecological_Regulation,
                                   Predator_Recovery=states$Predator_Recovery,
                                   Tick_Disease_Risk=states$Tick_Disease_Risk))
    cpt_time <- array(as.vector(time_base), dim=c(4,3),
                      dimnames=list(Time_to_Stable_State=states$Time_to_Stable_State,
                                    Ecological_Regulation=states$Ecological_Regulation))
    
    fit <- custom.fit(dag, dist=list(
      Region_Context        = cpt_region,
      Restoration_Intensity = cpt_rest,
      Host_Connectivity     = cpt_con,
      Predator_Recovery     = cpt_pred,
      Tick_Disease_Risk     = cpt_dis,
      Ecological_Regulation = cpt_reg,
      Time_to_Stable_State  = cpt_time
    ))
    as.grain(fit)
  }
  
  bn_base  <- build_bn()
  baseline <- query_risk(bn_base, "Tick_Disease_Risk","present",
                         c("Region_Context","Restoration_Intensity"),c("High","low"))
  
  results <- bind_rows(lapply(params, function(p){
    mat_lo <- mat_hi <- switch(p$node,
                               "con"  = con_base, "pred" = pred_base,
                               "dis"  = dis_base, "reg"  = reg_base)
    mat_lo[, p$col] <- p$lo
    mat_hi[, p$col] <- p$hi
    bn_lo <- do.call(build_bn, setNames(list(mat_lo), paste0(p$node,"_mat")))
    bn_hi <- do.call(build_bn, setNames(list(mat_hi), paste0(p$node,"_mat")))
    risk_lo <- query_risk(bn_lo,"Tick_Disease_Risk","present",
                          c("Region_Context","Restoration_Intensity"),c("High","low"))
    risk_hi <- query_risk(bn_hi,"Tick_Disease_Risk","present",
                          c("Region_Context","Restoration_Intensity"),c("High","low"))
    tibble(parameter=p$name, risk_low=risk_lo, risk_high=risk_hi,
           swing=risk_hi-risk_lo, abs_swing=abs(risk_hi-risk_lo),
           baseline=baseline, system="Tick-borne")
  }))
  results
}

# =========================================================
# SYSTEM 3: RODENT-BORNE
# =========================================================

run_rodent_sensitivity <- function(){
  
  states <- list(
    Region_Context        = c("Low","Moderate","High"),
    Restoration_Intensity = c("none","low","moderate","high"),
    Habitat_Quality       = c("low","moderate","high"),
    Rodent_Abundance      = c("low","moderate","high"),
    Rodent_Disease_Risk   = c("absent","present"),
    Ecological_Regulation = c("weak","moderate","strong"),
    Time_to_Stable_State  = c("1_3_years","4_7_years",
                              "8_12_years","13_15_years")
  )
  
  dag <- model2network(paste0(
    "[Region_Context]",
    "[Restoration_Intensity]",
    "[Habitat_Quality|Region_Context:Restoration_Intensity]",
    "[Rodent_Abundance|Habitat_Quality]",
    "[Rodent_Disease_Risk|Rodent_Abundance]",
    "[Ecological_Regulation|Restoration_Intensity:Rodent_Disease_Risk]",
    "[Time_to_Stable_State|Ecological_Regulation]"
  ))
  
  # Habitat quality: dim=c(3,4) cols=rest scenarios (mid-region)
  hab_base <- matrix(c(
    mid(.62,.80), mid(.14,.28), mid(.02,.10),  # none
    mid(.38,.58), mid(.26,.40), mid(.08,.20),  # low
    mid(.18,.35), mid(.32,.46), mid(.25,.42),  # moderate
    mid(.04,.14), mid(.18,.34), mid(.55,.75)   # high
  ), nrow=3, ncol=4)
  hab_base <- norm_cols(hab_base)
  
  # Rodent abundance: dim=c(3,3) cols=habitat quality states
  rod_base <- matrix(c(
    mid(.02,.08), mid(.15,.21), mid(.72,.80),  # low HQ
    mid(.25,.31), mid(.44,.52), mid(.21,.27),  # mod HQ
    mid(.68,.76), mid(.19,.25), mid(.03,.07)   # high HQ
  ), nrow=3, ncol=3)
  rod_base <- norm_cols(rod_base)
  
  # Disease: dim=c(2,3) cols=abundance states
  dis_base <- matrix(c(
    mid(.88,.96), mid(.04,.10),  # low abundance
    mid(.52,.58), mid(.42,.48),  # mod abundance
    mid(.06,.14), mid(.86,.94)   # high abundance
  ), nrow=2, ncol=3, byrow=TRUE)
  dis_base <- norm_cols(dis_base)
  
  # Regulation: dim=c(3,8) cols = (4 rest) x (2 disease)
  # Col order: disease outer, rest inner
  reg_base <- matrix(c(
    mid(.65,.82), mid(.12,.24), mid(.02,.10),  # dis=absent, none
    mid(.45,.62), mid(.24,.38), mid(.08,.22),  # dis=absent, low
    mid(.22,.38), mid(.30,.44), mid(.25,.42),  # dis=absent, moderate
    mid(.04,.12), mid(.12,.26), mid(.65,.82),  # dis=absent, high
    mid(.68,.85), mid(.12,.24), mid(.02,.10),  # dis=present, none
    mid(.48,.65), mid(.24,.38), mid(.08,.22),  # dis=present, low
    mid(.25,.41), mid(.30,.44), mid(.25,.42),  # dis=present, moderate
    mid(.07,.15), mid(.12,.26), mid(.65,.82)   # dis=present, high
  ), nrow=3, ncol=8, byrow=TRUE)
  reg_base <- norm_cols(reg_base)
  
  time_base <- matrix(c(
    mid(.48,.62), mid(.24,.34), mid(.06,.14), mid(.02,.08),
    mid(.14,.24), mid(.32,.44), mid(.22,.34), mid(.10,.22),
    mid(.03,.09), mid(.10,.20), mid(.42,.56), mid(.24,.38)
  ), nrow=4, ncol=3)
  time_base <- norm_cols(time_base)
  
  params <- list(
    
    list(name = "Rodent abundance: low habitat → high rodents",
         node = "rod", col = 1,
         lo   = norm(c(.05, .20, .68)),   # lower peak abundance
         hi   = norm(c(.01, .10, .87))),  # higher peak abundance
    
    list(name = "Rodent abundance: high habitat → low rodents (dilution)",
         node = "rod", col = 3,
         lo   = norm(c(.82, .14, .02)),   # strong habitat-mediated suppression
         hi   = norm(c(.48, .34, .14))),  # weak suppression
    
    list(name = "Rodent abundance: moderate habitat",
         node = "rod", col = 2,
         lo   = norm(c(.42, .46, .10)),   # moderate habitat strongly suppresses
         hi   = norm(c(.14, .44, .38))),  # moderate habitat weakly suppresses
    
    list(name = "Disease: high abundance → present",
         node = "dis", col = 3,
         lo   = norm(c(.18, .78)),        # lower spillover from high rodents
         hi   = norm(c(.02, .96))),       # near-certain spillover from high rodents
    
    list(name = "Disease: low abundance → absent",
         node = "dis", col = 1,
         lo   = norm(c(.96, .03)),        # very low background risk
         hi   = norm(c(.72, .26))),       # moderate background risk
    
    list(name = "Disease: moderate abundance (transition)",
         node = "dis", col = 2,
         lo   = norm(c(.72, .26)),        # moderate abundance mostly safe
         hi   = norm(c(.28, .70))),       # moderate abundance mostly risky
    
    list(name = "Habitat quality: none rest → low quality",
         node = "hab", col = 1,
         lo   = norm(c(.52, .32, .12)),   # degraded but not extreme
         hi   = norm(c(.88, .09, .01))),  # severely degraded
    
    list(name = "Habitat quality: high rest → high quality",
         node = "hab", col = 4,
         lo   = norm(c(.10, .28, .58)),   # slow habitat recovery
         hi   = norm(c(.01, .10, .87))),  # fast habitat recovery
    
    list(name = "Habitat quality: moderate rest (steepness)",
         node = "hab", col = 3,
         lo   = norm(c(.08, .30, .58)),   # steep improvement at mod rest
         hi   = norm(c(.38, .42, .16))),  # flat improvement at mod rest
    
    list(name = "Regulation: high rest → strong reg",
         node = "reg", col = 4,
         lo   = norm(c(.02, .10, .86)),   # fast regulation under high rest
         hi   = norm(c(.18, .30, .48))),  # slow regulation under high rest
    
    list(name = "Regulation: none rest → weak reg",
         node = "reg", col = 1,
         lo   = norm(c(.55, .28, .14)),   # moderate baseline regulation
         hi   = norm(c(.88, .09, .01))),  # complete regulatory failure
    
    list(name = "Regulation: disease penalty on recovery",
         node = "reg", col = 8,  # dis=present, high rest
         lo   = norm(c(.06, .16, .75)),   # small penalty: regulation recovers well
         hi   = norm(c(.28, .34, .34)))   # large penalty: disease delays regulation
  )
  
  build_bn <- function(hab_mat = hab_base,
                       rod_mat = rod_base,
                       dis_mat = dis_base,
                       reg_mat = reg_base){
    
    cpt_region <- array(c(0.33,0.34,0.33), dim=c(3),
                        dimnames=list(Region_Context=states$Region_Context))
    cpt_rest   <- array(rep(0.25,4), dim=c(4),
                        dimnames=list(Restoration_Intensity=states$Restoration_Intensity))
    
    # Habitat: expand to 3 regions
    hab_full <- matrix(rep(as.vector(hab_mat), 3), nrow=3)
    cpt_hab  <- array(as.vector(hab_full), dim=c(3,4,3),
                      dimnames=list(Habitat_Quality=states$Habitat_Quality,
                                    Restoration_Intensity=states$Restoration_Intensity,
                                    Region_Context=states$Region_Context))
    cpt_rod <- array(as.vector(rod_mat), dim=c(3,3),
                     dimnames=list(Rodent_Abundance=states$Rodent_Abundance,
                                   Habitat_Quality=states$Habitat_Quality))
    cpt_dis <- array(as.vector(dis_mat), dim=c(2,3),
                     dimnames=list(Rodent_Disease_Risk=states$Rodent_Disease_Risk,
                                   Rodent_Abundance=states$Rodent_Abundance))
    cpt_reg <- array(as.vector(reg_mat), dim=c(3,4,2),
                     dimnames=list(Ecological_Regulation=states$Ecological_Regulation,
                                   Restoration_Intensity=states$Restoration_Intensity,
                                   Rodent_Disease_Risk=states$Rodent_Disease_Risk))
    cpt_time <- array(as.vector(time_base), dim=c(4,3),
                      dimnames=list(Time_to_Stable_State=states$Time_to_Stable_State,
                                    Ecological_Regulation=states$Ecological_Regulation))
    
    fit <- custom.fit(dag, dist=list(
      Region_Context=cpt_region, Restoration_Intensity=cpt_rest,
      Habitat_Quality=cpt_hab, Rodent_Abundance=cpt_rod,
      Rodent_Disease_Risk=cpt_dis, Ecological_Regulation=cpt_reg,
      Time_to_Stable_State=cpt_time))
    as.grain(fit)
  }
  
  bn_base  <- build_bn()
  baseline <- query_risk(bn_base,"Rodent_Disease_Risk","present",
                         c("Region_Context","Restoration_Intensity"),c("High","moderate"))
  
  results <- bind_rows(lapply(params, function(p){
    mat_lo <- mat_hi <- switch(p$node,
                               "hab"=hab_base,"rod"=rod_base,"dis"=dis_base,"reg"=reg_base)
    mat_lo[, p$col] <- p$lo
    mat_hi[, p$col] <- p$hi
    bn_lo <- do.call(build_bn, setNames(list(mat_lo), paste0(p$node,"_mat")))
    bn_hi <- do.call(build_bn, setNames(list(mat_hi), paste0(p$node,"_mat")))
    risk_lo <- query_risk(bn_lo,"Rodent_Disease_Risk","present",
                          c("Region_Context","Restoration_Intensity"),c("High","moderate"))
    risk_hi <- query_risk(bn_hi,"Rodent_Disease_Risk","present",
                          c("Region_Context","Restoration_Intensity"),c("High","moderate"))
    tibble(parameter=p$name,risk_low=risk_lo,risk_high=risk_hi,
           swing=risk_hi-risk_lo,abs_swing=abs(risk_hi-risk_lo),
           baseline=baseline,system="Rodent-borne")
  }))
  results
}

# =========================================================
# SYSTEM 4: DIRECT-CONTACT ZOONOSES
# =========================================================

run_zoonotic_sensitivity <- function(){
  
  states <- list(
    Region_Context                  = c("Low_Risk_Context","Moderate_Risk_Context","High_Risk_Context"),
    Restoration_Intensity           = c("none","low","moderate","high"),
    Resilient_Ecological_Regulation = c("low","high"),
    Habitat_Disturbance             = c("low","moderate","high"),
    Host_Diversity                  = c("low","moderate","high"),
    Human_Wildlife_Contact          = c("low","moderate","high"),
    MultiHost_Amplification         = c("low","high"),
    Direct_Transmission_Risk        = c("absent","present"),
    Time_to_Stable_State            = c("1_3_years","4_7_years",
                                        "8_12_years","13_15_years")
  )
  
  dag <- model2network(paste0(
    "[Region_Context]",
    "[Restoration_Intensity]",
    "[Habitat_Disturbance|Region_Context:Restoration_Intensity]",
    "[Host_Diversity|Habitat_Disturbance]",
    "[Human_Wildlife_Contact|Habitat_Disturbance:Restoration_Intensity]",
    "[MultiHost_Amplification|Host_Diversity:Human_Wildlife_Contact]",
    "[Direct_Transmission_Risk|MultiHost_Amplification]",
    "[Resilient_Ecological_Regulation|Restoration_Intensity:Host_Diversity]",
    "[Time_to_Stable_State|Restoration_Intensity:Direct_Transmission_Risk:Resilient_Ecological_Regulation]"
  ))
  
  # Habitat disturbance: dim=c(3,4) cols=rest (mid-region)
  dist_base <- matrix(c(
    mid(.05,.15),mid(.20,.30),mid(.55,.75),  # none
    mid(.10,.25),mid(.35,.50),mid(.30,.45),  # low
    mid(.20,.40),mid(.35,.45),mid(.15,.30),  # moderate
    mid(.50,.70),mid(.20,.35),mid(.05,.15)   # high
  ), nrow=3, ncol=4)
  dist_base <- norm_cols(dist_base)
  
  # Host diversity: dim=c(3,3) cols=disturbance states
  div_base <- matrix(c(
    0.10,0.25,0.65,  # low disturbance
    0.30,0.45,0.25,  # mod disturbance
    0.70,0.20,0.10   # high disturbance
  ), nrow=3, ncol=3, byrow=TRUE)
  div_base <- norm_cols(div_base)
  
  # Human-wildlife contact: dim=c(3,12) cols=(4 rest)x(3 dist)
  # Col order: rest outer, dist inner
  con_base <- matrix(c(
    mid(.10,.20),mid(.25,.35),mid(.50,.65),  # none, dist=low
    mid(.15,.25),mid(.30,.40),mid(.54,.69),  # none, dist=mod  (+disturbance penalty)
    mid(.08,.18),mid(.23,.33),mid(.55,.68),  # none, dist=high
    mid(.15,.30),mid(.35,.45),mid(.30,.45),  # low, dist=low
    mid(.20,.35),mid(.38,.48),mid(.33,.48),  # low, dist=mod
    mid(.13,.28),mid(.33,.43),mid(.32,.46),  # low, dist=high
    mid(.30,.50),mid(.30,.42),mid(.15,.28),  # mod, dist=low
    mid(.35,.55),mid(.33,.45),mid(.18,.31),  # mod, dist=mod
    mid(.28,.48),mid(.28,.40),mid(.16,.29),  # mod, dist=high
    mid(.55,.75),mid(.18,.30),mid(.05,.15),  # high, dist=low
    mid(.60,.80),mid(.21,.33),mid(.08,.18),  # high, dist=mod
    mid(.53,.73),mid(.16,.28),mid(.06,.16)   # high, dist=high
  ), nrow=3, ncol=12, byrow=FALSE)
  con_base <- norm_cols(con_base)
  
  # Multi-host amplification: dim=c(2,9) cols=(3 contact)x(3 div)
  # Col order: contact outer, div inner
  amp_base <- matrix(c(
    mid(.70,.88),mid(.12,.30),  # contact=low, div=low
    mid(.30,.60),mid(.40,.70),  # contact=low, div=mod
    mid(.70,.88),mid(.12,.30),  # contact=low, div=high (dilution)
    mid(.30,.60),mid(.40,.70),  # contact=mod, div=low
    mid(.35,.65),mid(.35,.65),  # contact=mod, div=mod (max ambiguity)
    mid(.30,.60),mid(.40,.70),  # contact=mod, div=high
    mid(.05,.15),mid(.85,.95),  # contact=high, div=low  ← PEAK
    mid(.30,.55),mid(.45,.70),  # contact=high, div=mod  (ambiguous)
    mid(.30,.55),mid(.45,.70)   # contact=high, div=high (ambiguous)
  ), nrow=2, ncol=9, byrow=TRUE)
  amp_base <- norm_cols(amp_base)
  
  # Disease: dim=c(2,2) cols=amplification states
  dis_base <- matrix(c(
    mid(.74,.86),mid(.14,.26),  # low amp
    mid(.10,.20),mid(.80,.90)   # high amp
  ), nrow=2, ncol=2, byrow=TRUE)
  dis_base <- norm_cols(dis_base)
  
  # Regulation: dim=c(2,12) cols=(4 rest)x(3 div)
  reg_base <- matrix(c(
    mid(.68,.85),mid(.15,.32),  # none, div=low
    mid(.68,.85),mid(.15,.32),  # none, div=mod
    mid(.68,.85),mid(.15,.32)+mid(.03,.12),  # none, div=high (bonus)
    mid(.55,.72),mid(.28,.45),  # low,  div=low
    mid(.55,.72),mid(.28,.45),
    mid(.55,.72),mid(.28,.45)+mid(.03,.12),
    mid(.38,.58),mid(.42,.62),  # mod,  div=low
    mid(.38,.58),mid(.42,.62),
    mid(.38,.58),mid(.42,.62)+mid(.03,.12),
    mid(.15,.35),mid(.65,.85),  # high, div=low
    mid(.15,.35),mid(.65,.85),
    mid(.15,.35),mid(.65,.85)+mid(.03,.12)   # high, div=high (bonus)
  ), nrow=2, ncol=12, byrow=FALSE)
  reg_base <- norm_cols(reg_base)
  
  params <- list(
    
    list(name = "Contact: high rest → low contact (primary mechanism)",
         node = "con", col = 10,  # rest=high, dist=low
         lo   = norm(c(.72, .20, .06)),   # restoration very effective at reducing contact
         hi   = norm(c(.32, .38, .26))),  # restoration weakly reduces contact
    
    list(name = "Contact: none rest → high contact (baseline)",
         node = "con", col = 2,   # rest=none, dist=mod
         lo   = norm(c(.22, .38, .36)),   # lower baseline contact
         hi   = norm(c(.04, .18, .76))),  # very high baseline contact
    
    list(name = "Contact: disturbance × rest interaction",
         node = "con", col = 5,   # rest=low, dist=mod
         lo   = norm(c(.28, .42, .26)),   # disturbance penalty small
         hi   = norm(c(.06, .26, .64))),  # disturbance penalty large
    
    list(name = "Amplification: high contact + low diversity → high amp",
         node = "amp", col = 7,
         lo   = norm(c(.22, .74)),        # lower peak amplification
         hi   = norm(c(.02, .97))),       # near-certain amplification
    
    list(name = "Amplification: low contact + high diversity → low amp",
         node = "amp", col = 3,
         lo   = norm(c(.90, .08)),        # strong dilution effect
         hi   = norm(c(.42, .56))),       # weak/absent dilution effect
    
    list(name = "Amplification: mid conditions (ambiguity)",
         node = "amp", col = 5,
         lo   = norm(c(.72, .26)),        # moderate conditions → mostly safe
         hi   = norm(c(.18, .80))),       # moderate conditions → risky
    
    list(name = "Disease: high amplification → spillover",
         node = "dis", col = 2,
         lo   = norm(c(.28, .70)),        # lower spillover efficiency
         hi   = norm(c(.04, .95))),       # near-certain spillover
    
    list(name = "Disease: low amplification → absent",
         node = "dis", col = 1,
         lo   = norm(c(.94, .05)),        # very low background spillover
         hi   = norm(c(.58, .40))),       # moderate background spillover
    
    list(name = "Regulation: high rest → high reg",
         node = "reg", col = 10, # rest=high, div=low
         lo   = norm(c(.06, .92)),        # fast regulation recovery
         hi   = norm(c(.42, .56))),       # slow regulation recovery
    
    list(name = "Regulation: biodiversity bonus (strength)",
         node = "reg", col = 12, # rest=high, div=high
         lo   = norm(c(.04, .95)),        # strong diversity-regulation link
         hi   = norm(c(.38, .60))),       # weak diversity-regulation link
    
    list(name = "Disturbance: none rest → high disturbance",
         node = "dist", col = 1,
         lo   = norm(c(.18, .32, .44)),   # moderate baseline disturbance
         hi   = norm(c(.04, .14, .80))),  # severe baseline disturbance
    
    list(name = "Disturbance: high rest → low disturbance",
         node = "dist", col = 4,
         lo   = norm(c(.72, .20, .06)),   # strong disturbance reduction
         hi   = norm(c(.28, .42, .26)))   # weak disturbance reduction
  )
  
  build_bn <- function(dist_mat = dist_base,
                       con_mat  = con_base,
                       amp_mat  = amp_base,
                       dis_mat  = dis_base,
                       reg_mat  = reg_base){
    
    cpt_region <- array(c(0.33,0.34,0.33), dim=c(3),
                        dimnames=list(Region_Context=states$Region_Context))
    cpt_rest   <- array(rep(0.25,4), dim=c(4),
                        dimnames=list(Restoration_Intensity=states$Restoration_Intensity))
    
    dist_full <- matrix(rep(as.vector(dist_mat), 3), nrow=3)
    cpt_dist  <- array(as.vector(dist_full), dim=c(3,4,3),
                       dimnames=list(Habitat_Disturbance=states$Habitat_Disturbance,
                                     Restoration_Intensity=states$Restoration_Intensity,
                                     Region_Context=states$Region_Context))
    cpt_div <- array(as.vector(div_base), dim=c(3,3),
                     dimnames=list(Host_Diversity=states$Host_Diversity,
                                   Habitat_Disturbance=states$Habitat_Disturbance))
    cpt_con <- array(as.vector(con_mat), dim=c(3,3,4),
                     dimnames=list(Human_Wildlife_Contact=states$Human_Wildlife_Contact,
                                   Habitat_Disturbance=states$Habitat_Disturbance,
                                   Restoration_Intensity=states$Restoration_Intensity))
    cpt_amp <- array(as.vector(amp_mat), dim=c(2,3,3),
                     dimnames=list(MultiHost_Amplification=states$MultiHost_Amplification,
                                   Host_Diversity=states$Host_Diversity,
                                   Human_Wildlife_Contact=states$Human_Wildlife_Contact))
    cpt_dis <- array(as.vector(dis_mat), dim=c(2,2),
                     dimnames=list(Direct_Transmission_Risk=states$Direct_Transmission_Risk,
                                   MultiHost_Amplification=states$MultiHost_Amplification))
    
    reg_full <- matrix(rep(as.vector(reg_mat), 1), nrow=2)
    cpt_reg  <- array(as.vector(reg_mat), dim=c(2,4,3),
                      dimnames=list(
                        Resilient_Ecological_Regulation=states$Resilient_Ecological_Regulation,
                        Restoration_Intensity=states$Restoration_Intensity,
                        Host_Diversity=states$Host_Diversity))
    
    time_vals <- c()
    for(disease in 1:2) for(rest in 1:4) for(regulation in 1:2){
      probs <- if(regulation==2 & rest==4) norm(c(.06,.14,.32,.46))
      else if(rest==1)            norm(c(.44,.30,.15,.08))
      else                        norm(c(.22,.26,.28,.22))
      time_vals <- c(time_vals, probs)
    }
    cpt_time <- array(time_vals, dim=c(4,2,4,2),
                      dimnames=list(
                        Time_to_Stable_State=states$Time_to_Stable_State,
                        Resilient_Ecological_Regulation=states$Resilient_Ecological_Regulation,
                        Restoration_Intensity=states$Restoration_Intensity,
                        Direct_Transmission_Risk=states$Direct_Transmission_Risk))
    
    fit <- custom.fit(dag, dist=list(
      Region_Context=cpt_region, Restoration_Intensity=cpt_rest,
      Habitat_Disturbance=cpt_dist, Host_Diversity=cpt_div,
      Human_Wildlife_Contact=cpt_con, MultiHost_Amplification=cpt_amp,
      Direct_Transmission_Risk=cpt_dis,
      Resilient_Ecological_Regulation=cpt_reg,
      Time_to_Stable_State=cpt_time))
    as.grain(fit)
  }
  
  bn_base  <- build_bn()
  baseline <- query_risk(bn_base,"Direct_Transmission_Risk","present",
                         c("Region_Context","Restoration_Intensity"),
                         c("High_Risk_Context","moderate"))
  
  results <- bind_rows(lapply(params, function(p){
    mat_lo <- mat_hi <- switch(p$node,
                               "dist"=dist_base,"con"=con_base,"amp"=amp_base,
                               "dis"=dis_base,"reg"=reg_base)
    mat_lo[, p$col] <- p$lo
    mat_hi[, p$col] <- p$hi
    node_arg <- paste0(p$node,"_mat")
    bn_lo <- do.call(build_bn, setNames(list(mat_lo), node_arg))
    bn_hi <- do.call(build_bn, setNames(list(mat_hi), node_arg))
    risk_lo <- query_risk(bn_lo,"Direct_Transmission_Risk","present",
                          c("Region_Context","Restoration_Intensity"),
                          c("High_Risk_Context","moderate"))
    risk_hi <- query_risk(bn_hi,"Direct_Transmission_Risk","present",
                          c("Region_Context","Restoration_Intensity"),
                          c("High_Risk_Context","moderate"))
    tibble(parameter=p$name,risk_low=risk_lo,risk_high=risk_hi,
           swing=risk_hi-risk_lo,abs_swing=abs(risk_hi-risk_lo),
           baseline=baseline,system="Direct-contact zoonoses")
  }))
  results
}

# =========================================================
# RUN ALL SYSTEMS
# =========================================================

mosq_sens     <- run_mosquito_sensitivity()

tick_sens     <- run_tick_sensitivity()

rodent_sens   <- run_rodent_sensitivity()

zoonotic_sens <- run_zoonotic_sensitivity()

all_sens <- bind_rows(mosq_sens, tick_sens, rodent_sens, zoonotic_sens)

# =========================================================
# TORNADO PLOT FUNCTION
# =========================================================

make_tornado <- function(df, system_name, accent_col, top_n = 10,
                         query_scenario = NULL){
  d <- df %>%
    filter(system == system_name) %>%
    arrange(desc(abs_swing)) %>%
    slice_head(n = top_n) %>%
    mutate(
      parameter = str_wrap(parameter, width = 44),
      parameter = factor(parameter, levels = parameter[order(abs_swing)]),
      direction = if_else(swing >= 0, "Increases risk", "Decreases risk"),
      xlo = pmin(risk_low, risk_high),
      xhi = pmax(risk_low, risk_high)
    )
  
  x_pad <- diff(range(c(d$xlo, d$xhi))) * 0.18
  
  ggplot(d) +
    geom_segment(
      aes(x=xlo, xend=xhi, y=parameter, yend=parameter, colour=direction),
      linewidth=5.5, lineend="butt"
    ) +
    geom_vline(xintercept=unique(d$baseline),
               linetype="dashed", colour="grey35", linewidth=0.5) +
    geom_text(
      aes(x     = if_else(swing >= 0, xhi + x_pad*0.18, xlo - x_pad*0.18),
          y     = parameter,
          label = sprintf("%+.3f", swing),
          hjust = if_else(swing >= 0, 0, 1)),
      size=2.8, colour="grey25"
    ) +
    scale_colour_manual(
      values = c("Increases risk"="#d73027", "Decreases risk"=accent_col),
      name   = NULL
    ) +
    scale_x_continuous(
      limits = c(min(d$xlo) - x_pad, max(d$xhi) + x_pad*1.4),
      labels = scales::label_number(accuracy=0.01)
    ) +
    theme_minimal(base_size=10) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.y        = element_text(size=8, colour="grey20"),
      axis.text.x        = element_text(size=8),
      axis.title         = element_text(size=8, colour="grey40"),
      legend.position    = "bottom",
      legend.text        = element_text(size=8),
      plot.title         = element_text(size=10, face="bold"),
      plot.subtitle      = element_text(size=7.5, colour="grey50", face="italic"),
      plot.tag           = element_text(size=9, face="bold", colour="grey30")
    ) +
    labs(title    = system_name,
         subtitle = if(!is.null(query_scenario))
           paste0("Queried at: ", query_scenario) else NULL,
         x = "Disease risk probability",
         y = NULL)
}

# =========================================================
# BUILD + COMPOSE FIGURE
# =========================================================

p1 <- make_tornado(all_sens, "Mosquito-borne",    "#2171b5",
                   query_scenario="moderate restoration, high-risk region") + labs(tag="A")
p2 <- make_tornado(all_sens, "Tick-borne",         "#238b45",
                   query_scenario="low restoration (hump peak), high-risk region") + labs(tag="B")
p3 <- make_tornado(all_sens, "Rodent-borne",       "#cb181d",
                   query_scenario="moderate restoration, high-risk region") + labs(tag="C")
p4 <- make_tornado(all_sens, "Direct-contact zoonoses", "#6a51a3",
                   query_scenario="moderate restoration, high-risk region") + labs(tag="D")

fig_tornado <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title    = "Supplementary Figure S2 — One-at-a-time sensitivity analysis",
    subtitle = paste0(
      "Bars span disease risk probability from lower-bound to upper-bound assumption ",
      "for each parameter group, with all other parameters held at mid-range values. ",
      "Parameters ranked by absolute swing (top 10 per system). ",
      "Dashed line = baseline risk at mid-range values."
    ),
    theme = theme(
      plot.title    = element_text(size=11, face="bold"),
      plot.subtitle = element_text(size=8,  colour="grey50", margin=margin(b=10))
    )
  )

ggsave("supp_fig_S2_sensitivity_tornado.pdf",
       fig_tornado, width=260, height=210, units="mm", device=cairo_pdf)
ggsave("supp_fig_S2_sensitivity_tornado.png",
       fig_tornado, width=260, height=210, units="mm", dpi=600)

# =========================================================
# SUMMARY TABLE
# =========================================================

sens_table <- all_sens %>%
  arrange(system, desc(abs_swing)) %>%
  group_by(system) %>% slice_head(n=10) %>% ungroup() %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  select(system, parameter, baseline, risk_low, risk_high, swing, abs_swing)

print(sens_table, n=40)
write.csv(sens_table, "supp_table_S2_sensitivity.csv", row.names=FALSE)

																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		
																																		