# INPUT:
# θ       -- soil water content, 3 layers
# I       -- total water entering into soil surface (mm)
# pEc     -- potential ET allocated to plants (mm)
# pEs     -- potential ET allocated to soil surface (mm)
# soilpar -- soil-related parameters
# pftpar  -- plant-related parameters
# wet     -- wetness index
# Δz      -- soil layer depths, 3 layers
# zwt     -- groundwater table depth (mm)
function swb_case2(I, pEc, pEs, s_tem, s_vod, soilpar, pftpar, fwet, soil::Soil)
  (; θ, Δz, zwt, Ec_sm, Ec_gw) = soil
  (; Ksat, θ_sat, θ_fc, θ_wp) = soilpar

  # Unsaturated depth in layer #1~2
  d1 = Δz[1]
  d2 = zwt - d1

  # ====== Water Supplement ====== #  
  wa1_unsat = θ[1] # 需要更新
  vw2 = SM_recharge!(θ, I; Δz, θ_sat)
  wa2_unsat = find_θ_unsat(θ, zwt; z₊ₕ, Δz, θ_sat)[1]
  wa1, wa2, wa3 = θ

  # Layer #3 - Fully saturated with groundwater

  # ====== Water Consumption ====== #
  f_cons = s_tem * s_vod
  Tr, Es = Evapotranspiration!(soil, pEc, pEs, fwet, f_cons, soilpar, pftpar)

  # Transpiration from unsaturated and saturated zones in layer #2
  Tr1 = max(Ec_sm[1] + Ec_gw[1], 0.0)
  if wa1 > 0 && Es + Tr1 > d1 * wa1
    Tr1 = d1 * wa1 * Tr1 / (Tr1 + Es)
    Es = d1 * wa1 - Tr1
  end
  Tr2_u = Ec_sm[2]

  # ====== Soil Water Drainage (Unsaturated Zone) ====== #
  sink = [Tr1 + Es, Tr2_u, 0.0]
  θ_unsat = [wa1_unsat, wa2_unsat, θ[3]]
  exceed = SM_discharge!(soil, θ_unsat, sink, soilpar)
  wa1, wa2_unsat, _ = θ_unsat

  # ====== Groundwater Table Depth Update ====== #
  Δw = exceed + vw2 - sum(Ec_gw) - GW_Rsb(zwt)

  sy = θ_sat - (wa1 + wa2_unsat) / 2
  zwt -= Δw / sy
  uex = 0  # Excess water to soil surface

  # Update soil moisture and groundwater table depth
  if zwt > z₊ₕ[3]
    wa2 = (wa2_unsat * d2 + θ_fc * (Δz[2] - d2)) / Δz[2]
    wa3 = θ_fc
  elseif z₊ₕ[2] < zwt <= z₊ₕ[3]
    wa2 = (wa2_unsat * d2 + θ_fc * (Δz[2] - d2)) / Δz[2]
    wa3 = (θ_fc * (zwt - z₊ₕ[2]) + θ_sat * (z₊ₕ[3] - zwt)) / Δz[3]
  elseif z₊ₕ[1] < zwt <= z₊ₕ[2] # [x]
    wa2 = (wa2_unsat * (zwt - z₊ₕ[1]) + θ_sat * (z₊ₕ[2] - zwt)) / Δz[2]
    wa3 = θ_sat
  elseif 0 < zwt <= Δz[1]
    wa1 = (wa1 * zwt + θ_sat * (Δz[1] - zwt)) / Δz[1]
    wa2 = θ_sat
    wa3 = θ_sat
  elseif zwt <= 0
    wa1 = θ_sat
    wa2 = θ_sat
    wa3 = θ_sat
    uex = -zwt * θ_fc
  end

  soil.θ .= [wa1, wa2, wa3]
  soil.zwt = max(0, zwt)
  return Tr, Es, uex
end
