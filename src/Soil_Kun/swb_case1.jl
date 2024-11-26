# INPUT:
# θ       -- soil water content, 3 layers
# I       -- total water entering into the soil surface (mm)
# pEc     -- potential ET allocated to plants (mm)
# pEs     -- potential ET allocated to soil surface (mm)
# soilpar -- soil-related parameters
# pftpar  -- plant-related parameters
# fwet     -- wetness index
function swb_case1(soil::Soil, I, pEc, pEs, fwet, s_tem, s_vod, soilpar, pftpar)
  (; θ, Δz, zwt, Ec_gw, sink) = soil
  (; θ_sat, θ_fc) = soilpar

  # ====== Water supplement ====== #
  # Layer #1 - Unsaturated zone
  vw1 = SM_recharge!(θ, I; Δz, θ_sat)
  wa1_unsat = find_θ_unsat(θ, zwt; z₊ₕ, Δz, θ_sat)[1]

  # Layer #2 and #3 - Fully saturated
  f_cons = s_tem * s_vod
  Tr, Es = Evapotranspiration!(soil, pEc, pEs, fwet, f_cons, soilpar, pftpar)

  # ====== Soil Water Drainage (Unsaturated Zone) ====== #
  θ_unsat = [wa1_unsat, θ[2], θ[3]]
  exceed = SM_discharge!(soil, θ_unsat, sink, soilpar)
  wa1_unsat = θ_unsat[1]

  # ====== Groundwater table depth update ====== #
  sy = θ_sat - wa1_unsat
  Δw = exceed + vw1 - sum(Ec_gw) - GW_Rsb(zwt)
  zwt -= Δw / sy

  uex = update_wa!(soil, θ_unsat, soil.zwt, zwt)
  soil.zwt = max(0, zwt)
  return Tr, Es, uex
end
