# INPUT:
# θ       -- soil water content, 3 layers
# I       -- total water enter into soil surface, mm
# pEc     -- potential ET allocate to plant, mm
# pEs     -- potential ET allocate to soil surface, mm
# soilpar -- soil-related parameters
# pftpar  -- plant-related parameters
# fwet     -- wetness indice
function swb_case3(soil, I, pEc, pEs, fwet, s_tem, s_vod, soilpar, pftpar)
  (; θ, Δz, zwt, Ec_gw, sink) = soil
  (; θ_sat, θ_fc) = soilpar

  wa1_unsat = θ[1] # 需要更新
  wa2_unsat = θ[2] # 需要更新
  vw3 = SM_recharge!(θ, I; Δz, θ_sat)
  wa3_unsat = find_θ_unsat(θ, zwt; z₊ₕ, Δz, θ_sat)[1]
  wa1, wa2, wa3 = θ

  # ====== Water Consumption ====== #
  f_cons = s_tem * s_vod
  Tr, Es = Evapotranspiration!(soil, pEc, pEs, fwet, f_cons, soilpar, pftpar)

  # ====== Soil Water Drainage (Unsaturated Zone) ====== #  
  θ_unsat = [wa1_unsat, wa2_unsat, wa3_unsat]
  exceed = SM_discharge!(soil, θ_unsat, sink, soilpar)
  wa1, wa2, wa3_unsat = θ_unsat
  
  # ====== The Groundwater Table Depth ====== #
  Δw = exceed + vw3 - sum(Ec_gw) - GW_Rsb(zwt)
  sy = θ_sat - (wa1 + wa2 + wa3_unsat) / 3
  zwt = zwt - Δw / sy

  uex = update_wa!(soil, θ_unsat, soil.zwt, zwt)
  soil.zwt = max(0, zwt)
  return Tr, Es, uex
end
