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
function swb_case2(soil::Soil, I, pEc, pEs, fwet, s_tem, s_vod, soilpar, pftpar)
  (; θ, Δz, zwt, Ec_gw, sink) = soil
  (; θ_sat, θ_fc) = soilpar

  # ====== Water Supplement ====== #  
  wa1_unsat = θ[1] # 需要更新
  vw2 = SM_recharge!(θ, I; Δz, θ_sat)
  wa2_unsat = find_θ_unsat(θ, zwt; z₊ₕ, Δz, θ_sat)[1]
  wa1, wa2, wa3 = θ

  # ====== Water Consumption ====== #
  f_cons = s_tem * s_vod
  Tr, Es = Evapotranspiration!(soil, pEc, pEs, fwet, f_cons, soilpar, pftpar)

  # ====== Soil Water Drainage (Unsaturated Zone) ====== #
  θ_unsat = [wa1_unsat, wa2_unsat, θ[3]]
  exceed = SM_discharge!(soil, θ_unsat, sink, soilpar)
  wa1, wa2_unsat, _ = θ_unsat

  # ====== Groundwater Table Depth Update ====== #
  sy = θ_sat - (wa1 + wa2_unsat) / 2
  Δw = exceed + vw2 - sum(Ec_gw) - GW_Rsb(zwt)
  zwt = zwt - Δw / sy
  
  uex = update_wa!(soil, θ_unsat, soil.zwt, zwt)
  soil.zwt = max(0, zwt)
  return Tr, Es, uex
end
