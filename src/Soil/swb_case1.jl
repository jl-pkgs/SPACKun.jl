# INPUT:
# θ       -- soil water content, 3 layers
# I       -- total water entering into the soil surface (mm)
# pEc     -- potential ET allocated to plants (mm)
# pEs     -- potential ET allocated to soil surface (mm)
# soilpar -- soil-related parameters
# pftpar  -- plant-related parameters
# fwet     -- wetness index
function swb_case1(I, pEc, pEs, s_tem, s_vod, soilpar, pftpar, fwet, soil::Soil)
  (; θ, Δz, zwt, Ec_sm, Ec_gw, sink) = soil
  (; Ksat, θ_sat, θ_fc, θ_wp) = soilpar
  d1 = zwt  # Unsaturated depth in layer #1

  # ====== Water supplement ====== #
  # Layer #1 - Unsaturated zone
  vw1 = SM_recharge!(θ, I; Δz, θ_sat)
  wa1_unsat, frac_unsat = find_θ_unsat(θ, zwt; z₊ₕ, Δz, θ_sat)#[1]
  wa1, wa2, wa3 = θ

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
