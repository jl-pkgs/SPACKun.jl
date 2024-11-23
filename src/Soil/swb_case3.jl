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

  Tr1 = max(Ec_sm[1] + Ec_gw[1], 0)
  Tr2 = Ec_sm[2] + Ec_gw[2]
  Tr3_u = Ec_sm[3]

  # 限制蒸发
  if wa1 > 0 && Es + Tr1 > d1 * wa1
    Tr1 = d1 * wa1 * Tr1 / (Tr1 + Es)
    Es = d1 * wa1 - Tr1
  end
  Tr2 = clamp(Tr2, 0, d2 * (wa2 - θ_wp))
  Tr3_u = clamp(Tr3_u, 0, d3 * (wa3_unsat - θ_wp))

  # 新版
  # # ====== Soil Water Drainage (Unsaturated Zone) ====== #  
  sink = [Tr1 + Es, Tr2, Tr3_u]
  θ_unsat = [wa1_unsat, wa2_unsat, wa3_unsat]
  exceed = SM_discharge!(soil, θ_unsat, sink, soilpar)
  wa1, wa2, wa3_unsat = θ_unsat
  # wa3 = wa3_unsat
  
  # ====== The Groundwater Table Depth ====== #
  Δw = exceed + vw3 - sum(Ec_gw) - GW_Rsb(zwt)

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

