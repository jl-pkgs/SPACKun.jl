function swb_case2(wa, IWS, pEc, pEs, s_tem, s_vod, soilpar, pftpar, wet, ZM, zgw)
  # INPUT:
  # wa      -- soil water content, 3 layers
  # IWS     -- total water entering into soil surface (mm)
  # pEc     -- potential ET allocated to plants (mm)
  # pEs     -- potential ET allocated to soil surface (mm)
  # soilpar -- soil-related parameters
  # pftpar  -- plant-related parameters
  # wet     -- wetness index
  # zm      -- soil layer depths, 3 layers
  # zgw     -- groundwater table depth (mm)

  # Unsaturated depth in layer #1~2
  d1 = ZM[1]
  d2 = zgw - d1

  wa1, wa2, wa3 = wa

  ks = soilpar[1]  # hydraulic conductivity
  theta_sat = soilpar[3]  # saturated soil water content
  theta_fc = soilpar[5]  # field capacity
  wwp = soilpar[7]  # wilting point

  # ====== Water Supplement ====== #
  # Layer #1 - Unsaturated zone
  wa1_unsat = wa1
  wc_s1 = d1 * wa1_unsat
  wc_m1 = d1 * theta_sat

  if wc_s1 + IWS >= wc_m1
    wa1 = theta_sat
    vw1 = wc_s1 + IWS - wc_m1
  else
    wa1 = wa1_unsat + IWS / d1
    vw1 = 0
  end

  # Layer #2 - Unsaturated zone
  wa2_unsat = (wa2 * ZM[2] - theta_sat * (ZM[2] - d2)) / d2
  wc_s2 = d2 * wa2_unsat
  wc_m2 = d2 * theta_sat

  if wc_s2 + vw1 >= wc_m2
    wa2_unsat = theta_sat
    vw2 = wc_s2 + vw1 - wc_m2
  else
    wa2_unsat += vw1 / d2
    wa2 = (wa2_unsat * d2 + theta_sat * (ZM[2] - d2)) / ZM[2]
    vw2 = 0
  end

  # Layer #3 - Fully saturated with groundwater

  # ====== Water Consumption ====== #
  # Evapotranspiration
  Tr_p1, Tr_p2, Tr_p3 = pTr_partition(pEc, wa1, wa2, wa3, soilpar, pftpar, wet, ZM)

  # Transpiration from unsaturated and saturated zones in layer #2
  Tr_p2_u = Tr_p2 * (d2 * wa2_unsat) / (d2 * wa2_unsat + (ZM[2] - d2) * theta_sat)
  Tr_p2_g = Tr_p2 * ((ZM[2] - d2) * theta_sat) / (d2 * wa2_unsat + (ZM[2] - d2) * theta_sat)

  # Moisture constraints
  f_sm1, f_sm_s1 = swc_stress(wa1, soilpar, pEc, pftpar)
  f_sm2, _ = swc_stress(wa2_unsat, soilpar, pEc, pftpar)

  # Actual transpiration
  Tr1 = f_sm1 * s_vod * s_tem * Tr_p1
  Tr2_u = clamp(f_sm2 * s_vod * s_tem * Tr_p2_u, 0, d2 * (wa2_unsat - wwp))
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

  f1 = soil_drainage(wa1_unsat, theta_sat, ks, 0.048, 4.8)
  wa1 = max((wa1 * d1 - f1 - Es - Tr1) / d1, 0)

  # Layer #2
  f2 = soil_drainage(wa2_unsat, theta_sat, ks, 0.012, 1.2)
  wa2_unsat = clamp((wa2_unsat * d2 + f1 - f2 - Tr2_u) / d2, 0, 1)

  ff2 = max((wa2_unsat - theta_sat) * d2, 0)
  wa2_unsat = min(wa2_unsat, theta_sat)

  # Layer #3 - Fully saturated with groundwater

  # ====== Groundwater Table Depth Update ====== #
  F1 = f2 + ff2 + vw2
  Tr_g = Tr2_g + Tr3

  R_sb_max = 39  # mm day-1
  f = 1.25e-3  # mm-1
  R_sb = R_sb_max * exp(-f * zgw)

  delta_w = F1 - Tr_g - R_sb
  delta_zgw = delta_w / (theta_sat - (wa1 + wa2_unsat) / 2)
  zgw -= delta_zgw
  uex = 0  # Excess water to soil surface

  # Update soil moisture and groundwater table depth
  if zgw > ZM[1] + ZM[2] + ZM[3]
    wa2 = (wa2_unsat * d2 + theta_fc * (ZM[2] - d2)) / ZM[2]
    wa3 = theta_fc
  elseif zgw > ZM[1] + ZM[2] && zgw <= ZM[1] + ZM[2] + ZM[3]
    wa2 = (wa2_unsat * d2 + theta_fc * (ZM[2] - d2)) / ZM[2]
    wa3 = (theta_fc * (zgw - ZM[1] - ZM[2]) + theta_sat * (ZM[1] + ZM[2] + ZM[3] - zgw)) / ZM[3]
  elseif zgw > ZM[1] && zgw <= ZM[1] + ZM[2]
    wa2 = (wa2_unsat * (zgw - ZM[1]) + theta_sat * (ZM[1] + ZM[2] - zgw)) / ZM[2]
    wa3 = theta_sat
  elseif zgw > 0 && zgw <= ZM[1]
    wa1 = (wa1 * zgw + theta_sat * (ZM[1] - zgw)) / ZM[1]
    wa2 = theta_sat
    wa3 = theta_sat
  elseif zgw <= 0
    wa1 = theta_sat
    wa2 = theta_sat
    wa3 = theta_sat
    uex = -zgw * theta_fc
  end

  wa = [wa1, wa2, wa3]
  zgw = max(0, zgw)

  return wa, zgw, Tr, Es, uex
end
