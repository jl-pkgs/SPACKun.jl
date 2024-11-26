# INPUT:
# θ       -- soil water content, 3 layers
# I       -- total water entering into soil surface (mm)
# pEc     -- potential ET allocated to plants (mm)
# pEs     -- potential ET allocated to soil surface (mm)
# soilpar -- soil-related parameters
# pftpar  -- plant-related parameters
# fwet     -- wetness index
function swb_case0(soil::Soil, I, pEc, pEs, fwet, s_tem, s_vod, soilpar, pftpar)
  (; θ, Δz, zwt, Ec_gw, Ec_gw) = soil
  (; θ_sat, θ_fc) = soilpar

  f_cons = s_tem * s_vod
  Tr, Es = Evapotranspiration!(soil, pEc, pEs, fwet, f_cons, soilpar, pftpar)

  # Change in groundwater table depth
  sy = θ_sat - θ_fc
  Δw = I - sum(Ec_gw) - GW_Rsb(zwt)
  zwt -= Δw / sy

  uex = update_wa!(soil, θ, soil.zwt, zwt)
  soil.zwt = max(0, zwt)  # ensure groundwater table depth is non-negative
  return Tr, Es, uex
end
