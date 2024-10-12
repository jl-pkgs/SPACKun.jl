function swb_case2(wa, IWS, pEc, pEs, s_tem, s_vod, soilpar, pftpar, wet, zm, zgw)
  # INPUT:
  # wa      -- soil water content, 3 layers
  # IWS     -- total water enter into soil surface, mm
  # pEc     -- potential ET allocate to plant, mm
  # pEs     -- potential ET allocate to soil surface, mm
  # soilpar -- soil-related parameters
  # pftpar  -- plant-related parameters
  # wet     -- wetness indice
  # zm      -- soil layer depth, 3 layers
  # zgw     -- groundwater table depth, mm

  # unsaturated depth in layer #2
  d2 = zgw - zm[1]

  wa1 = wa[1]
  wa2 = wa[2]
  wa3 = wa[3]

  ks = soilpar[1]  # hydraulic conductivity
  theta_sat = soilpar[3]  # saturated swc
  theta_fc = soilpar[5]  # field water capacity
  wwp = soilpar[7]  # wilting point

  # ====== water supplement ====== #
  # layer #2
  # existed water column in the unsaturated zone #2
  wa2_unsat = (wa2 * zm[2] - theta_sat * (zm[2] - d2)) / d2  # 未饱和区域SM, [m3 m-3]
  wc_s2 = d2 * wa2_unsat  # [mm]
  wc_m2 = d2 * theta_sat  # maximum water column in d2, [mm]

  if wc_s2 + IWS >= wc_m2
    wa2_unsat = theta_sat  # current soil water content
    vw2 = wc_s2 + IWS - wc_m2  # exceeded water
  else
    # soil water content in unsaturated zone
    wa2_unsat += IWS / d2
    # calculate the adjusted swc#2 with considering the groundwater depth
    wa2 = (wa2_unsat * d2 + theta_sat * (zm[2] - d2)) / zm[2]
    vw2 = 0
  end

  # layer #1 - saturated
  # layer #3 - saturated

  # ====== water consumption ====== #
  # Evapotranspiration
  # distributed the potential Tr to different layers
  Tr_p1, Tr_p2, Tr_p3 = pTr_partition(pEc, wa1, wa2, wa3, soilpar, pftpar, wet, zm)

  # divide Tr_p2 into unsaturated zone and saturated zone
  Tr_p2_u = Tr_p2 * (d2 * wa2_unsat) / (d2 * wa2_unsat + (zm[2] - d2) * theta_sat)
  Tr_p2_g = Tr_p2 * ((zm[2] - d2) * theta_sat) / (d2 * wa2_unsat + (zm[2] - d2) * theta_sat)

  # calculate the moisture constrains for plant and soil in unsaturated zone
  f_sm2, f_sm_s2 = swc_stress(wa2, soilpar, pEc, pftpar)

  # actual transpiration
  Tr2_u = f_sm2 * s_vod * s_tem * Tr_p2_u
  Tr2_u = clamp(Tr2_u, 0, d2 * (wa2_unsat - wwp))  # less than maximum available water

  Tr2_g = s_vod * s_tem * Tr_p2_g
  Tr2 = Tr2_u + Tr2_g
  Tr1 = s_vod * s_tem * Tr_p1
  Tr3 = s_vod * s_tem * Tr_p3

  Tr = Tr1 + Tr2 + Tr3

  # actual soil evaporation
  Es = f_sm_s2 * pEs  # only consider about the first layer
  Es_u = Es * (d2 * wa2_unsat) / (d2 * wa2_unsat + (zm[2] - d2) * theta_sat)
  Es_u = clamp(Es_u, 0, d2 * wa2_unsat)

  # soil water drainage (unsaturated zone)
  # layer #2
  # drainage from unsaturated zone, #2
  f2 = soil_drainage(wa2_unsat, theta_sat, ks, 0.012, 1.2)

  # update the soil moisture after ET & drainage, layer #2
  wa2_unsat = (wa2_unsat * d2 - f2 - Es_u - Tr2_u) / d2
  wa2_unsat = max(wa2_unsat, 0)

  # layer #1   full filled with groundwater
  # layer #3   full filled with groundwater

  # The groundwater table depth
  # total water recharge to groundwater
  F2 = f2 + vw2

  # total transpiration from groundwater
  Tr_g = Tr2_g + Tr1 + Tr3

  # R_sb groundwater discharge
  R_sb_max = 39  # mm day-1
  f = 1.25e-3  # mm-1
  R_sb = R_sb_max * exp(-f * zgw)

  # variation of water stored in the saturated zone
  delta_w = F2 - Tr_g - R_sb

  # changes of groundwater table depth
  delta_zgw = delta_w / (theta_sat - (wa1 + wa2_unsat) / 2)
  zgw -= delta_zgw
  uex = 0  # excess water to soil surface, mm

  # update soil moisture and groundwater table depth
  if zgw > zm[1] + zm[2] + zm[3]
    wa2 = (wa2_unsat * d2 + theta_fc * (zm[2] - d2)) / zm[2]
    wa3 = theta_fc
  elseif zgw > zm[1] + zm[2] && zgw < zm[1] + zm[2] + zm[3]
    wa2 = (wa2_unsat * d2 + theta_fc * (zm[2] - d2)) / zm[2]
    wa3 = (theta_fc * (zgw - zm[1] - zm[2]) + theta_sat * (zm[1] + zm[2] + zm[3] - zgw)) / zm[3]
  elseif zgw > zm[1] && zgw < zm[1] + zm[2]
    wa2 = (wa2_unsat * (zgw - zm[1]) + theta_sat * (zm[1] + zm[2] - zgw)) / zm[2]
    wa3 = theta_sat
  elseif zgw > 0 && zgw < zm[1]
    wa1 = (wa1 * zgw + theta_sat * (zm[1] - zgw)) / zm[1]
    wa2 = theta_sat
    wa3 = theta_sat
  elseif zgw <= 0
    wa1 = theta_sat
    wa2 = theta_sat
    wa3 = theta_sat
    uex = -zgw * theta_fc  # excess water to soil surface, mm
  end

  # updated soil water content
  wa = [wa1, wa2, wa3]
  zgw = max(0, zgw)

  return wa, zgw, Tr, Es, uex
end
