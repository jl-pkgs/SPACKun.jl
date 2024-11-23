# 按照CoLM的框架更新土壤水
function SW_balance(I, pEc, pEs, s_tem, s_vod, soilpar, pftpar, fwet, soil)
  (; θ, Δz, zwt, Ec_gw, sink) = soil
  (; θ_sat, θ_fc) = soilpar
  # d3 = zwt - z₊ₕ[2]
  # wa1_unsat = θ[1] # 需要更新
  # wa2_unsat = θ[2] # 需要更新
  
  wa = 0.0
  GW_UpdateRecharge!(soil, wa, zwt, Δt, I)
  # vw3 = SM_recharge!(θ, I; Δz, θ_sat)
  # wa3_unsat = find_θ_unsat(θ, zwt; z₊ₕ, Δz, θ_sat)[1]
  # wa1, wa2, wa3 = θ

  # ====== Water Consumption ====== #
  f_cons = s_tem * s_vod
  Tr, Es = Evapotranspiration!(soil, pEc, pEs, fwet, f_cons, soilpar, pftpar)

  # ====== Soil Water Drainage (Unsaturated Zone) ====== #  
  # θ_unsat = [wa1_unsat, wa2_unsat, wa3_unsat]
  exceed = SM_discharge!(soil, θ_unsat, sink, soilpar)
  # wa1, wa2, wa3_unsat = θ_unsat

  # ====== The Groundwater Table Depth ====== #
  Δw = exceed + vw3 - sum(Ec_gw) - GW_Rsb(zwt)
  sy = θ_sat - (wa1 + wa2 + wa3_unsat) / 3
  zwt = zwt - Δw / sy
  uex = 0

  ## 更新zwt, θ, uex
  soil.θ .= [wa1, wa2, wa3]
  soil.zwt = max(0, zwt)
  return Tr, Es, uex
end
