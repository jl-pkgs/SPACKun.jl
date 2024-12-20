"""
按照CoLM的框架更新土壤水

- `补给` : SM_recharge!
- `排泄` : 
  + `非饱和带`: SM_discharge!
  + `饱和带`  : sum(Ec_gw) + GW_Rsb(zwt)
"""
function sw_balance_CoLM(soil::Soil, I::T, pEc::T, pEs::T, Ta::T, Topt::T, fwet, s_VOD::T) where {T<:Real}

  (; θ, Δz, zwt, wa, Ec_gw) = soil
  (; θ_sat) = soilpar

  s_tem = exp(-((Ta - Topt) / Topt)^2)
  f_cons = s_tem * s_VOD

  extra = SM_recharge!(θ, I; Δz, θ_sat) # 补给-> θ，不影响水位

  Tr, Es = Evapotranspiration!(soil, pEc, pEs, fwet, f_cons)

  exceed = SM_discharge!(soil, soilpar) # 排泄 -> θ, 非饱和带

  Δw = exceed + extra - sum(Ec_gw) - GW_Rsb(zwt)
  zwt, wa, uex = Update_zwt_theta!(soil, zwt, wa, Δw; θ, wa_max=5000.0) # Δw -> [θ, zwt], CoLM

  # sy = θ_sat - (wa1 + wa2 + wa3_unsat) / 3
  # zwt = zwt - Δw / sy
  # uex = 0
  return Tr, Es, uex
end
