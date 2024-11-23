# INPUT:
# θ       -- soil water content, 3 layers
# I       -- total water entering into soil surface (mm)
# pEc     -- potential ET allocated to plants (mm)
# pEs     -- potential ET allocated to soil surface (mm)
# soilpar -- soil-related parameters
# pftpar  -- plant-related parameters
# fwet     -- wetness index
function swb_case0(I, pEc, pEs, s_tem, s_vod, soilpar, pftpar, fwet, soil::Soil)
  (; θ, Δz, zwt, Ec_gw, Ec_gw) = soil
  (; θ_sat, θ_fc) = soilpar

  wa1, wa2, wa3 = θ  # Old soil water content in layer 1-3

  f_cons = s_tem * s_vod
  Tr, Es = Evapotranspiration!(soil, pEc, pEs, fwet, f_cons, soilpar, pftpar)

  # Change in groundwater table depth
  sy = θ_sat - θ_fc
  Δw = I - sum(Ec_gw) - GW_Rsb(zwt)
  zwt -= Δw / sy
  uex = 0             # excess water to soil surface

  # Update soil moisture and groundwater table depth
  if zwt > z₊ₕ[3]
    wa1 = θ_fc
    wa2 = θ_fc
    wa3 = θ_fc
  elseif z₊ₕ[2] < zwt <= z₊ₕ[3]
    wa1 = θ_fc
    wa2 = θ_fc
    wa3 = (θ_fc * (zwt - z₊ₕ[2]) + θ_sat * (z₊ₕ[3] - zwt)) / Δz[3]
  elseif z₊ₕ[1] < zwt <= z₊ₕ[2]
    wa1 = θ_fc
    wa2 = (θ_fc * (zwt - z₊ₕ[1]) + θ_sat * (z₊ₕ[2] - zwt)) / Δz[2]
    wa3 = θ_sat
  elseif 0 < zwt <= Δz[1]
    wa1 = (θ_fc * zwt + θ_sat * (Δz[1] - zwt)) / Δz[1]
    wa2 = θ_sat
    wa3 = θ_sat
  elseif zwt <= 0 # [x]
    wa1 = θ_sat
    wa2 = θ_sat
    wa3 = θ_sat
    uex = -zwt * θ_fc  # excess water to soil surface
  end

  # Updated soil water content
  soil.θ .= [wa1, wa2, wa3]
  soil.zwt = max(0, zwt)  # ensure groundwater table depth is non-negative

  return Tr, Es, uex
end
