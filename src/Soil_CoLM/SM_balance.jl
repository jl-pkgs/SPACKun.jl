# 按照CoLM的框架更新土壤水
function SW_balance(I, pEc, pEs, s_tem, s_vod, soilpar, pftpar, fwet, soil)
  (; θ, Δz, zwt, Ec_gw) = soil
  (; θ_sat) = soilpar

  extra = SM_recharge!(θ, I; Δz, θ_sat)
  
  f_cons = s_tem * s_vod
  Tr, Es = Evapotranspiration!(soil, pEc, pEs, fwet, f_cons, soilpar, pftpar)

  # ====== Soil Water Drainage (Unsaturated Zone) ====== #  
  exceed = SM_discharge!(soil, soilpar)

  Δw = exceed + extra - sum(Ec_gw) - GW_Rsb(zwt)
  
  ## 给水度的更新
  sy = θ_sat - (wa1 + wa2 + wa3_unsat) / 3
  zwt = zwt - Δw / sy
  uex = 0

  ## 更新zwt, θ, uex
  soil.θ .= [wa1, wa2, wa3]
  soil.zwt = max(0, zwt)
  return Tr, Es, uex
end
