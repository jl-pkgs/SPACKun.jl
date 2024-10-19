function swb_case4(wa, IWS, pEc, pEs, s_tem, s_vod, soilpar, pftpar, wet, zm, zgw)
  # INPUT:
  # wa      -- soil water content, 3 layers
  # IWS     -- total water entering into the soil surface (mm)
  # pEc     -- potential ET allocated to plants (mm)
  # pEs     -- potential ET allocated to soil surface (mm)
  # soilpar -- soil-related parameters
  # pftpar  -- plant-related parameters
  # wet     -- wetness index
  # zm      -- soil layer depths, 3 layers
  # zgw     -- groundwater table depth (mm)

  # Unsaturated depth in layer #1~3
  d1 = zm[1]
  d2 = zm[2]
  d3 = zm[3]

  wa1, wa2, wa3 = wa

  ks = soilpar[1]  # hydraulic conductivity
  θ_sat = soilpar[3]  # saturated soil water content
  wwp = soilpar[7]  # wilting point

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
  wa2_unsat = wa2
  wc_s2 = d2 * wa2_unsat
  wc_m2 = d2 * θ_sat

  if wc_s2 + vw1 >= wc_m2
    wa2 = θ_sat
    vw2 = wc_s2 + vw1 - wc_m2
  else
    wa2 = wa2_unsat + vw1 / d2
    vw2 = 0
  end

  # Layer #3 - Unsaturated zone
  wa3_unsat = wa3
  wc_s3 = d3 * wa3_unsat
  wc_m3 = d3 * θ_sat

  if wc_s3 + vw2 >= wc_m3
    wa3 = θ_sat
    vw3 = wc_s3 + vw2 - wc_m3
  else
    wa3 = wa3_unsat + vw2 / d3
    vw3 = 0
  end

  # ====== Water Consumption ====== #
  # Evapotranspiration
  Tr_p1, Tr_p2, Tr_p3 = pTr_partition(pEc, wa1, wa2, wa3, soilpar, pftpar, wet, zm)

  # Moisture constraints for each layer
  f_sm1, f_sm_s1 = swc_stress(wa1, soilpar, pEc, pftpar)
  f_sm2, _ = swc_stress(wa2, soilpar, pEc, pftpar)
  f_sm3, _ = swc_stress(wa3, soilpar, pEc, pftpar)

  # Actual transpiration
  Tr1 = f_sm1 * s_vod * s_tem * Tr_p1
  Tr2 = f_sm2 * s_vod * s_tem * Tr_p2
  Tr3 = f_sm3 * s_vod * s_tem * Tr_p3
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

  f1 = soil_drainage(wa1_unsat, θ_sat, ks, 0.048, 4.8)
  wa1 = clamp((wa1 * d1 - f1 - Es - Tr1) / d1, 0, θ_sat)

  # Layer #2
  Tr2 = max(Tr2, 0)
  Tr2 = min(Tr2, d2 * (wa2 - wwp))

  f2 = soil_drainage(wa2_unsat, θ_sat, ks, 0.012, 1.2)
  wa2 = clamp((wa2 * d2 + f1 - f2 - Tr2) / d2, 0, θ_sat)

  ff2 = max((wa2 - θ_sat) * d2, 0)
  wa2 = min(wa2, θ_sat)

  # Layer #3
  Tr3 = max(Tr3, 0)
  Tr3 = min(Tr3, d3 * (wa3 - wwp))

  f3 = soil_drainage(wa3_unsat, θ_sat, ks, 0.012, 1.2)
  wa3 = clamp((wa3 * d3 + f2 + ff2 - f3 - Tr3) / d3, 0, θ_sat)

  ff3 = max((wa3 - θ_sat) * d3, 0)
  wa3 = min(wa3, θ_sat)

  # ====== Groundwater Table Depth Update ====== #
  F1 = f3 + ff3 + vw3
  Tr_g = 0

  R_sb_max = 39  # mm day-1
  f = 1.25e-3  # mm-1
  R_sb = R_sb_max * exp(-f * zgw)

  delta_w = F1 - Tr_g - R_sb
  delta_zgw = delta_w / 0.2  # Specific yield = 0.2
  zgw -= delta_zgw
  uex = 0

  # Update soil moisture and groundwater table depth
  if zgw > zm[1] + zm[2] && zgw < zm[1] + zm[2] + zm[3]
    wa3 = (wa3_unsat * (zgw - zm[1] - zm[2]) + θ_sat * (zm[1] + zm[2] + zm[3] - zgw)) / zm[3]
  elseif zgw > zm[1] && zgw < zm[1] + zm[2]
    wa2 = (wa2 * (zgw - zm[1]) + θ_sat * (zm[1] + zm[2] - zgw)) / zm[2]
    wa3 = θ_sat
  elseif zgw > 0 && zgw < zm[1]
    wa1 = (wa1 * zgw + θ_sat * (zm[1] - zgw)) / zm[1]
    wa2 = θ_sat
    wa3 = θ_sat
  elseif zgw <= 0
    wa1 = θ_sat
    wa2 = θ_sat
    wa3 = θ_sat
    uex = -zgw * θ_sat
  end

  wa = [wa1, wa2, wa3]
  zgw = max(0, zgw)

  return wa, zgw, Tr, Es, uex
end
