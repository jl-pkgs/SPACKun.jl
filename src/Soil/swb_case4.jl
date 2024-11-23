# INPUT:
# θ       -- soil water content, 3 layers
# I       -- total water enter into soil surface, mm
# pEc     -- potential ET allocate to plant, mm
# pEs     -- potential ET allocate to soil surface, mm
# soilpar -- soil-related parameters
# pftpar  -- plant-related parameters
# wet     -- wetness indice
# Δz      -- soil layer depth, 3 layers
# zwt     -- groundwater table depth, mm
function swb_case4(I, pEc, pEs, s_tem, s_vod, soilpar, pftpar, fwet, soil::Soil)
  (; θ_prev, θ, Δz, zwt, 
    fsm_Es, Ec_sm, Ec_gw) = soil
  # Unsaturated depth in layer #1~3
  d1 = Δz[1]
  d2 = Δz[2]
  d3 = Δz[3]

  (; Ksat, θ_sat, θ_wp) = soilpar
  
  # # ====== Water Supplement ====== #
  # θ_prev .= θ
  wa1_unsat, wa2_unsat, wa3_unsat = θ # 需要更新
  vw3 = SM_recharge!(θ, I; Δz, θ_sat)
  wa1, wa2, wa3 = θ

  pTr_partition!(soil, pEc, fwet, soilpar, pftpar)
  swc_stress!(soil, pEc, soilpar, pftpar, s_tem * s_vod)

  # Evapotranspiration
  Tr1 = Ec_sm[1] + Ec_gw[1]
  Tr2 = Ec_sm[2] + Ec_gw[2]
  Tr3 = Ec_sm[3] + Ec_gw[3]
  Tr = sum(Ec_sm) + sum(Ec_gw)

  Es = max(fsm_Es[1] * pEs, 0.0)  # Only consider about the first layer

  # ====== Soil Water Drainage (Unsaturated Zone) ====== #
  # ---------------------------------------------------------------- Layer #1
  # Update the soil moisture after ET, layer #1
  Tr1 = max(Tr1, 0)

  if wa1 > 0 && Es + Tr1 > d1 * wa1  # wilting point
    Tr1 = d1 * wa1 * Tr1 / (Tr1 + Es)
    Es = d1 * wa1 - Tr1
  end

  # Drainage from unsaturated zone, #1
  f1 = soil_drainage(wa1_unsat, θ_sat, Ksat, 0.048, 4.8)

  # Update the soil moisture after drainage, layer #1
  wa1 = (wa1 * d1 - f1 - Es - Tr1) / d1
  wa1 = clamp(wa1, 0, 1)

  # ---------------------------------------------------------------- Layer #2
  # Update the soil moisture after ET, layer #2
  Tr2 = clamp(Tr2, 0, d2 * (wa2 - θ_wp))

  # Drainage from unsaturated zone, #2
  f2 = soil_drainage(wa2_unsat, θ_sat, Ksat, 0.012, 1.2)

  # Update the soil moisture after drainage, layer #2
  wa2 = (wa2 * d2 + f1 - f2 - Tr2) / d2
  wa2 = clamp(wa2, 0, 1)  # > wwp 0

  if wa2 > θ_sat
    ff2 = max((wa2 - θ_sat) * d2, 0)  # extra water from upper layer
    wa2 = θ_sat
  else
    ff2 = 0
  end

  # ---------------------------------------------------------------- Layer #3
  # Update the soil moisture after ET, layer #3, unsat-zone
  Tr3 = clamp(Tr3, 0, d3 * (wa3 - θ_wp))

  # Drainage from unsaturated zone, #3
  f3 = soil_drainage(wa3_unsat, θ_sat, Ksat, 0.012, 1.2)

  # Update the soil moisture after drainage, layer #3
  wa3 = (wa3 * d3 + f2 + ff2 - f3 - Tr3) / d3

  wa3 = clamp(wa3, 0, 1)
  if wa3 > θ_sat
    ff3 = max((wa3 - θ_sat) * d3, 0)
    wa3_unsat = θ_sat
  else
    ff3 = 0
    wa3_unsat = wa3
  end

  # ====== The Groundwater Table Depth ====== #
  # Total water recharge to groundwater
  F1 = f3 + ff3 + vw3

  # Total transpiration from groundwater
  Tr_g = 0

  # R_sb groundwater discharge
  R_sb = GW_Rsb(zwt)
  R_sb_max = 39  # mm day-1
  f = 1.25e-3  # mm-1
  R_sb = R_sb_max * exp(-f * zwt)

  # Variation of water stored in the saturated zone
  delta_w = F1 - Tr_g - R_sb

  # Changes in groundwater table depth
  sy = 0.2 # specific yield as 0.2
  zwt = zwt - delta_w / sy
  uex = 0  # excess water to soil surface, mm

  # Update soil moisture and groundwater table depth
  if z₊ₕ[2] < zwt < z₊ₕ[3]
    wa3 = (wa3_unsat * (zwt - z₊ₕ[2]) + θ_sat * (z₊ₕ[3] - zwt)) / Δz[3]
  elseif z₊ₕ[1] < zwt < z₊ₕ[2]
    wa2 = (wa2 * (zwt - z₊ₕ[1]) + θ_sat * (z₊ₕ[2] - zwt)) / Δz[2]
    wa3 = θ_sat
  elseif 0 < zwt < z₊ₕ[1]
    wa1 = (wa1 * zwt + θ_sat * (Δz[1] - zwt)) / Δz[1]
    wa2 = θ_sat
    wa3 = θ_sat
  elseif zwt <= 0
    wa1 = θ_sat
    wa2 = θ_sat
    wa3 = θ_sat
    uex = -zwt * θ_sat  # excess water to soil surface, mm
  end

  # Updated soil water content
  soil.θ .= [wa1, wa2, wa3]
  soil.zwt = max(0, zwt)
  return Tr, Es, uex
end
