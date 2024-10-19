function sw_balance(IWS, pEc, pEs, Ta, Topt, s_VOD, soilpar, pftpar, wet, ZM, wa, zgw)
  # ----INPUT:
  # IWS     -- total water enter into soil surface, mm
  # pEc     -- potential ET allocate to plant, mm
  # pEs     -- potential ET allocate to soil surface, mm
  # Ta      -- air temperature, C
  # Topt    -- optimal growth temperature for plant, C
  # s_VOD   -- constrains of VOD,[0,1]
  # wa      -- previous soil water content, 3 layers
  # soilpar -- soil-related parameters
  # pftpar  -- plant-related parameters
  # wet     -- wetness fraction indice
  # ZM      -- soil layer depth, 3 layers
  # zgw     -- groundwater table depth, mm
  # ----OUTPUT:
  # Tr      -- actual plant transpiration, mm
  # Es      -- actual soil evaporation, mm
  # wa      -- updated soil water content, 3 layers, %
  # zgw     -- groundwater table depth, mm
  # uex     -- goundwater overflow soil surface, mm
  # ----------

  # Constrains of temperature
  s_tem = temp_stress(Topt, Ta)

  # Case 0: groundwater overflow
  if zgw <= 0
    wa, zgw, Tr, Es, uex = swb_case0(wa, IWS, pEc, pEs, s_tem, s_VOD, soilpar, pftpar, wet, ZM, zgw)
    # Case 1: groundwater table in layer 1
  elseif zgw > 0 && zgw <= ZM[1]
    wa, zgw, Tr, Es, uex = swb_case1(wa, IWS, pEc, pEs, s_tem, s_VOD, soilpar, pftpar, wet, ZM, zgw)
    # Case 2: groundwater table in layer 2
  elseif zgw > ZM[1] && zgw <= ZM[1] + ZM[2]
    wa, zgw, Tr, Es, uex = swb_case2(wa, IWS, pEc, pEs, s_tem, s_VOD, soilpar, pftpar, wet, ZM, zgw)
    # Case 3: groundwater table in layer 3
  elseif zgw > ZM[1] + ZM[2] && zgw < ZM[1] + ZM[2] + ZM[3]
    wa, zgw, Tr, Es, uex = swb_case3(wa, IWS, pEc, pEs, s_tem, s_VOD, soilpar, pftpar, wet, ZM, zgw)
    # Case 4: groundwater table below layer 3
  else
    wa, zgw, Tr, Es, uex = swb_case4(wa, IWS, pEc, pEs, s_tem, s_VOD, soilpar, pftpar, wet, ZM, zgw)
  end

  return wa, zgw, Tr, Es, uex
end

# Temperature Constrains for plant growing
# INPUT:
# - Topt  : Optimum temperature for plant growing
# - Ta    : Air temperature
# 
# OUTPUT:
# - St    : Temperature stress
function temp_stress(Topt, Ta)
  return exp(-((Ta - Topt) / Topt)^2)
end

export sw_balance
