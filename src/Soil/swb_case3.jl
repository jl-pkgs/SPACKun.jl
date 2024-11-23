## 如果exceed依然>0，则补给到地下水

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
function swb_case3(I, pEc, pEs, s_tem, s_vod, soilpar, pftpar, fwet, soil)
  (; θ, Δz, zwt, Ec_sm, Ec_gw) = soil
  (; Ksat, θ_sat, θ_fc, θ_wp) = soilpar
  # Unsaturated depth in layer #1~3
  d1 = Δz[1]
  d2 = Δz[2]
  d3 = zwt - z₊ₕ[2]

  # ====== Water Supplement ====== #
  wa1_unsat = θ[1] # 需要更新
  wa2_unsat = θ[2] # 需要更新
  vw3 = SM_recharge!(θ, I; Δz, θ_sat)
  wa3_unsat = find_θ_unsat(θ, zwt; z₊ₕ, Δz, θ_sat)[1]
  wa1, wa2, wa3 = θ

  # ====== Water Consumption ====== #
  f_cons = s_tem * s_vod
  Tr, Es = Evapotranspiration!(soil, pEc, pEs, fwet, f_cons, soilpar, pftpar)

  Tr1 = Ec_sm[1] + Ec_gw[1]
  Tr2 = Ec_sm[2] + Ec_gw[2]
  Tr3_u = Ec_sm[3]
  
  # ====== Soil Water Drainage (Unsaturated Zone) ====== #
  Tr1 = max(Tr1, 0)

  if wa1 > 0 && Es + Tr1 > d1 * wa1
    Tr1 = d1 * wa1 * Tr1 / (Tr1 + Es)
    Es = d1 * wa1 - Tr1
  end

  f1 = soil_drainage(wa1_unsat, θ_sat, Ksat, 0.048, 4.8)
  wa1 = (wa1 * d1 - f1 - Es - Tr1) / d1
  wa1 = max(wa1, 0)

  Tr2 = clamp(Tr2, 0, d2 * (wa2 - θ_wp))
  f2 = soil_drainage(wa2_unsat, θ_sat, Ksat, 0.012, 1.2)
  wa2 = (wa2 * d2 + f1 - f2 - Tr2) / d2
  wa2 = clamp(wa2, 0, 1)

  if wa2 > θ_sat
    ff2 = max((wa2 - θ_sat) * d2, 0)
    wa2 = θ_sat
  else
    ff2 = 0
  end

  Tr3_u = clamp(Tr3_u, 0, d3 * (wa3_unsat - θ_wp))
  
  f3 = soil_drainage(wa3_unsat, θ_sat, Ksat, 0.012, 1.2)
  wa3_unsat = (wa3_unsat * d3 + f2 + ff2 - f3 - Tr3_u) / d3
  wa3_unsat = clamp(wa3_unsat, 0, 1)

  if wa3_unsat > θ_sat
    ff3 = max((wa3_unsat - θ_sat) * d3, 0)
    wa3_unsat = θ_sat
  else
    ff3 = 0
  end

  # ====== The Groundwater Table Depth ====== #
  F1 = f3 + ff3 + vw3
  # Tr_g = Tr3_g
  Δw = F1 - sum(Ec_gw) - GW_Rsb(zwt)

  sy = θ_sat - (wa1 + wa2 + wa3_unsat) / 3
  zwt = zwt - Δw / sy
  uex = 0

  if zwt > z₊ₕ[3]
    wa3 = (wa3_unsat * d3 + θ_fc * (Δz[3] - d3)) / Δz[3]
  elseif z₊ₕ[2] < zwt <= z₊ₕ[3] # [x]
    wa3 = (wa3_unsat * (zwt - z₊ₕ[2]) + θ_sat * (z₊ₕ[3] - zwt)) / Δz[3]
  elseif z₊ₕ[1] < zwt <= z₊ₕ[2] # 水位上升
    wa2 = (wa2 * (zwt - z₊ₕ[1]) + θ_sat * (z₊ₕ[2] - zwt)) / Δz[2]
    wa3 = θ_sat
  elseif 0 < zwt <= Δz[1]
    wa1 = (wa1 * zwt + θ_sat * (Δz[1] - zwt)) / Δz[1]
    wa2 = θ_sat
    wa3 = θ_sat
  elseif zwt <= 0
    wa1 = θ_sat
    wa2 = θ_sat
    wa3 = θ_sat
    uex = -zwt * θ_sat
  end

  soil.θ .= [wa1, wa2, wa3]
  soil.zwt = max(0, zwt)
  return Tr, Es, uex
end

