# INPUT:
# θ       -- soil water content, 3 layers
# I       -- total water entering into the soil surface (mm)
# pEc     -- potential ET allocated to plants (mm)
# pEs     -- potential ET allocated to soil surface (mm)
# soilpar -- soil-related parameters
# pftpar  -- plant-related parameters
# wet     -- wetness index
# Δz      -- soil layer depths, 3 layers
# zwt     -- groundwater table depth (mm)
function swb_case1(I, pEc, pEs, s_tem, s_vod, soilpar, pftpar, fwet, state::State)
  (; θ, Δz, zwt, Ec_pot, fsm_Ec, fsm_Es) = state
  # Unsaturated depth in layer #1
  d1 = zwt
  (; Ksat, θ_sat, θ_fc, θ_wp) = soilpar
  
  # ====== Water supplement ====== #
  # Layer #1 - Unsaturated zone
  vw1 = SM_recharge!(θ, I; Δz, θ_sat)
  wa1_unsat, j = find_θ_unsat(θ, zwt; z₊ₕ, Δz, θ_sat)
  wa1, wa2, wa3 = θ

  # Layer #2 and #3 - Fully saturated

  # ====== Water consumption ====== #
  # Evapotranspiration
  pTr_partition!(pEc, fwet, soilpar, pftpar, state)
  swc_stress!(pEc, soilpar, pftpar, state)
  # Moisture constraints
  # fsm_Ec[1], fsm_Es[1] = swc_stress(wa1, pEc, soilpar, pftpar)


  # Transpiration from unsaturated and saturated zones in layer #1
  Tr_p1_u = Ec_pot[1] * (d1 * wa1_unsat) / (d1 * wa1_unsat + (Δz[1] - d1) * θ_sat)
  Tr_p1_g = Ec_pot[1] * ((Δz[1] - d1) * θ_sat) / (d1 * wa1_unsat + (Δz[1] - d1) * θ_sat)

  # Actual transpiration
  Tr1_u = clamp( fsm_Ec[1] * s_vod * s_tem * Tr_p1_u, 0, d1 * (wa1_unsat - θ_wp))
  Tr1_g = s_vod * s_tem * Tr_p1_g
  Tr1 = Tr1_u + Tr1_g
  Tr2 = s_vod * s_tem * Ec_pot[2]
  Tr3 = s_vod * s_tem * Ec_pot[3]
  Tr = Tr1 + Tr2 + Tr3

  # Actual soil evaporation
  Es = fsm_Es[1] * pEs
  Es_u = clamp(Es * (d1 * wa1_unsat) / (d1 * wa1_unsat + (Δz[1] - d1) * θ_sat), 0, d1 * wa1_unsat)

  # ====== Soil water drainage (unsaturated zone) ====== #
  f1 = soil_drainage(wa1_unsat, θ_sat, Ksat, 0.048, 4.8)
  wa1_unsat = max((wa1_unsat * d1 - f1 - Es_u - Tr1_u) / d1, 0)

  # ====== Groundwater table depth update ====== #
  F1 = f1 + vw1  # Total water recharge to groundwater
  Tr_g = Tr1_g + Tr2 + Tr3  # Total transpiration from groundwater

  # Groundwater discharge
  R_sb_max = 39  # mm day-1
  f = 1.25e-3  # mm-1
  R_sb = R_sb_max * exp(-f * zwt)

  # Variation of water stored in the saturated zone
  delta_w = F1 - Tr_g - R_sb

  # Change in groundwater table depth
  sy = θ_sat - wa1_unsat
  zwt -= delta_w / sy
  uex = 0  # Excess water to the soil surface

  # Update soil moisture and groundwater table depth
  if zwt > z₊ₕ[3]
    wa1 = (wa1_unsat * d1 + θ_fc * (Δz[1] - d1)) / Δz[1]
    wa2 = θ_fc
    wa3 = θ_fc
  elseif z₊ₕ[2] < zwt <= z₊ₕ[3]
    wa1 = (wa1_unsat * d1 + θ_fc * (Δz[1] - d1)) / Δz[1]
    wa2 = θ_fc
    wa3 = (θ_fc * (zwt - z₊ₕ[2]) + θ_sat * (z₊ₕ[3] - zwt)) / Δz[3]
  elseif z₊ₕ[1] < zwt <= z₊ₕ[2]
    wa1 = (wa1_unsat * d1 + θ_fc * (Δz[1] - d1)) / Δz[1]
    wa2 = (θ_fc * (zwt - z₊ₕ[1]) + θ_sat * (z₊ₕ[2] - zwt)) / Δz[2]
    wa3 = θ_sat
  elseif 0 < zwt <= z₊ₕ[1] # [x]
    wa1 = (wa1_unsat * zwt + θ_sat * (z₊ₕ[1] - zwt)) / Δz[1]
    wa2 = θ_sat
    wa3 = θ_sat
  elseif zwt <= 0
    wa1 = θ_sat
    wa2 = θ_sat
    wa3 = θ_sat
    uex = -zwt * θ_fc  # Excess water to soil surface
  end

  # Updated soil water content
  state.θ .= [wa1, wa2, wa3]
  state.zwt = max(0, zwt)
  return Tr, Es, uex
end
