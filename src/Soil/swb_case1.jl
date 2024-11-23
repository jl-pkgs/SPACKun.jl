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
function swb_case1(I, pEc, pEs, s_tem, s_vod, soilpar, pftpar, fwet, soil::Soil)
  (; θ, Δz, zwt, fsm_Es, Ec_sm, Ec_gw) = soil
  # Unsaturated depth in layer #1
  d1 = zwt
  (; Ksat, θ_sat, θ_fc, θ_wp) = soilpar
  
  # ====== Water supplement ====== #
  # Layer #1 - Unsaturated zone
  vw1 = SM_recharge!(θ, I; Δz, θ_sat)
  wa1_unsat = find_θ_unsat(θ, zwt; z₊ₕ, Δz, θ_sat)[1]
  wa1, wa2, wa3 = θ

  # Layer #2 and #3 - Fully saturated

  # ====== Water consumption ====== #
  # Evapotranspiration
  pTr_partition!(soil, pEc, fwet, soilpar, pftpar)
  swc_stress!(soil, pEc, soilpar, pftpar, s_tem * s_vod)

  # Transpiration from unsaturated and saturated zones in layer #1
  Tr1_u = Ec_sm[1]
  Tr = sum(Ec_sm) + sum(Ec_gw)

  # Actual soil evaporation
  Es = fsm_Es[1] * pEs
  Es_u = clamp(Es * (d1 * wa1_unsat) / (d1 * wa1_unsat + (Δz[1] - d1) * θ_sat), 0, d1 * wa1_unsat)

  # ====== Soil water drainage (unsaturated zone) ====== #
  f1 = soil_drainage(wa1_unsat, θ_sat, Ksat, 0.048, 4.8)
  wa1_unsat = max((wa1_unsat * d1 - f1 - Es_u - Tr1_u) / d1, 0)

  # ====== Groundwater table depth update ====== #
  F1 = f1 + vw1  # Total water recharge to groundwater
  delta_w = F1 - sum(Ec_gw) - GW_Rsb(zwt)

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
  soil.θ .= [wa1, wa2, wa3]
  soil.zwt = max(0, zwt)
  return Tr, Es, uex
end
