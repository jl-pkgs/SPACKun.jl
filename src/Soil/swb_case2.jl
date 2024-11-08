# INPUT:
# θ      -- soil water content, 3 layers
# IWS     -- total water entering into soil surface (mm)
# pEc     -- potential ET allocated to plants (mm)
# pEs     -- potential ET allocated to soil surface (mm)
# soilpar -- soil-related parameters
# pftpar  -- plant-related parameters
# wet     -- wetness index
# Δz      -- soil layer depths, 3 layers
# zwt     -- groundwater table depth (mm)
function swb_case2(θ, IWS, pEc, pEs, s_tem, s_vod, soilpar, pftpar, fwet, Δz, zwt)
  # Unsaturated depth in layer #1~2
  d1 = Δz[1]
  d2 = zwt - d1
  wa1, wa2, wa3 = θ
  (; Ksat, θ_sat, θ_fc, θ_wp) = soilpar

  # ====== Water Supplement ====== #
  # Layer #1 - Unsaturated zone
  wa1_unsat = wa1
  wc_s1 = d1 * wa1_unsat
  wc_m1 = d1 * θ_sat

  if wc_s1 + IWS >= wc_m1
    wa1 = θ_sat
    vw1 = wc_s1 + IWS - wc_m1
  else
    wa1 = wa1_unsat + IWS / d1
    vw1 = 0
  end

  # Layer #2 - Unsaturated zone
  wa2_unsat = (wa2 * Δz[2] - θ_sat * (Δz[2] - d2)) / d2
  wc_s2 = d2 * wa2_unsat
  wc_m2 = d2 * θ_sat

  if wc_s2 + vw1 >= wc_m2
    wa2_unsat = θ_sat
    vw2 = wc_s2 + vw1 - wc_m2
  else
    wa2_unsat += vw1 / d2
    wa2 = (wa2_unsat * d2 + θ_sat * (Δz[2] - d2)) / Δz[2]
    vw2 = 0
  end

  # Layer #3 - Fully saturated with groundwater

  # ====== Water Consumption ====== #
  # Evapotranspiration
  Tr_p1, Tr_p2, Tr_p3 = pTr_partition(pEc, fwet, wa1, wa2, wa3, soilpar, pftpar, Δz)

  # Transpiration from unsaturated and saturated zones in layer #2
  Tr_p2_u = Tr_p2 * (d2 * wa2_unsat) / (d2 * wa2_unsat + (Δz[2] - d2) * θ_sat)
  Tr_p2_g = Tr_p2 * ((Δz[2] - d2) * θ_sat) / (d2 * wa2_unsat + (Δz[2] - d2) * θ_sat)

  # Moisture constraints
  f_sm1, f_sm_s1 = swc_stress(wa1, pEc, soilpar, pftpar)
  f_sm2, _ = swc_stress(wa2_unsat, pEc, soilpar, pftpar)

  # Actual transpiration
  Tr1 = f_sm1 * s_vod * s_tem * Tr_p1
  Tr2_u = clamp(f_sm2 * s_vod * s_tem * Tr_p2_u, 0, d2 * (wa2_unsat - θ_wp))
  Tr2_g = s_vod * s_tem * Tr_p2_g
  Tr2 = Tr2_u + Tr2_g
  Tr3 = s_vod * s_tem * Tr_p3
  Tr = Tr1 + Tr2 + Tr3

  # Actual soil evaporation
  Es = f_sm_s1 * pEs

  # ====== Soil Water Drainage (Unsaturated Zone) ====== #
  # Layer #1
  Es = max(Es, 0)
  Tr1 = max(Tr1, 0)

  if wa1 > 0 && Es + Tr1 > d1 * wa1
    Tr1 = d1 * wa1 * Tr1 / (Tr1 + Es)
    Es = d1 * wa1 - Tr1
  end

  f1 = soil_drainage(wa1_unsat, θ_sat, Ksat, 0.048, 4.8)
  wa1 = max((wa1 * d1 - f1 - Es - Tr1) / d1, 0)

  # Layer #2
  f2 = soil_drainage(wa2_unsat, θ_sat, Ksat, 0.012, 1.2)
  wa2_unsat = clamp((wa2_unsat * d2 + f1 - f2 - Tr2_u) / d2, 0, 1)

  ff2 = max((wa2_unsat - θ_sat) * d2, 0)
  wa2_unsat = min(wa2_unsat, θ_sat)

  # Layer #3 - Fully saturated with groundwater

  # ====== Groundwater Table Depth Update ====== #
  F1 = f2 + ff2 + vw2
  Tr_g = Tr2_g + Tr3

  R_sb_max = 39  # mm day-1
  f = 1.25e-3  # mm-1
  R_sb = R_sb_max * exp(-f * zwt)

  delta_w = F1 - Tr_g - R_sb
  delta_zgw = delta_w / (θ_sat - (wa1 + wa2_unsat) / 2)
  zwt -= delta_zgw
  uex = 0  # Excess water to soil surface

  # Update soil moisture and groundwater table depth
  if zwt > z₊ₕ[3]
    wa2 = (wa2_unsat * d2 + θ_fc * (Δz[2] - d2)) / Δz[2]
    wa3 = θ_fc
  elseif zwt > z₊ₕ[2] && zwt <= z₊ₕ[3]
    wa2 = (wa2_unsat * d2 + θ_fc * (Δz[2] - d2)) / Δz[2]
    wa3 = (θ_fc * (zwt - z₊ₕ[2]) + θ_sat * (z₊ₕ[3] - zwt)) / Δz[3]
  elseif zwt > Δz[1] && zwt <= z₊ₕ[2]
    wa2 = (wa2_unsat * (zwt - z₊ₕ[1]) + θ_sat * (z₊ₕ[2] - zwt)) / Δz[2]
    wa3 = θ_sat
  elseif zwt > 0 && zwt <= Δz[1]
    wa1 = (wa1 * zwt + θ_sat * (Δz[1] - zwt)) / Δz[1]
    wa2 = θ_sat
    wa3 = θ_sat
  elseif zwt <= 0
    wa1 = θ_sat
    wa2 = θ_sat
    wa3 = θ_sat
    uex = -zwt * θ_fc
  end

  θ = [wa1, wa2, wa3]
  zwt = max(0, zwt)

  return θ, zwt, Tr, Es, uex
end
