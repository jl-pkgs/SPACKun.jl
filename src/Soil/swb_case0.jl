# INPUT:
# θ       -- soil water content, 3 layers
# I       -- total water entering into soil surface (mm)
# pEc     -- potential ET allocated to plants (mm)
# pEs     -- potential ET allocated to soil surface (mm)
# soilpar -- soil-related parameters
# pftpar  -- plant-related parameters
# wet     -- wetness index
# Δz      -- soil layer depth (3 layers)
# zwt     -- groundwater table depth (mm)
function swb_case0(I, pEc, pEs, s_tem, s_vod, soilpar, pftpar, fwet, state::State)
  (; θ, Δz, zwt) = state

  # Old soil water content in layer 1-3
  wa1, wa2, wa3 = θ
  (; θ_sat, θ_fc) = soilpar

  # ====== Water consumption ====== #
  # Evapotranspiration
  Tr_p1, Tr_p2, Tr_p3 = pTr_partition(pEc, fwet, θ, soilpar, pftpar, Δz)

  # Actual transpiration
  Tr1 = s_vod * s_tem * Tr_p1
  Tr2 = s_vod * s_tem * Tr_p2
  Tr3 = s_vod * s_tem * Tr_p3
  Tr = Tr1 + Tr2 + Tr3

  # Actual soil evaporation (only first layer)
  Es = pEs

  # Groundwater recharge and discharge calculations
  F1 = I
  Tr_g = Tr
  R_sb_max = 39  # maximum groundwater discharge (mm day-1)
  f = 1.25e-3    # discharge decay rate (mm-1)
  R_sb = R_sb_max * exp(-f * zwt)

  # Variation of water stored in the saturated zone
  delta_w = F1 - Tr_g - R_sb

  # Change in groundwater table depth
  sy = θ_sat - θ_fc
  zwt -= delta_w / sy
  uex = 0  # excess water to soil surface

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
  state.θ .= [wa1, wa2, wa3]
  state.zwt = max(0, zwt)  # ensure groundwater table depth is non-negative

  return Tr, Es, uex
end
