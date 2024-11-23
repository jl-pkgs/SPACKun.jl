# s_vod * s_tem
function Evapotranspiration!(soil::Soil, pEc::T, pEs::T, fwet::T, f_cons::T,
  soilpar, pftpar) where {T<:Real}
  pTr_partition!(soil, pEc, fwet, soilpar, pftpar) # EC_POT划分

  (; θ_sat, θ_wp) = soilpar
  (; θ, Ec_pot, fsm_Ec, fsm_Es, 
    Ec_sm, Ec_gw, sink, 
    zwt, z₊ₕ, N) = soil

  j = find_jwt(z₊ₕ, zwt)

  # 非饱和
  for i = 1:max(j - 1, 0) # 全部非饱和
    fsm_Ec[i], fsm_Es[i] = swc_stress(θ[i], pEc, soilpar, pftpar)

    Ec_sm[i] = Ec_pot[i] * f_cons * fsm_Ec[i]
    # Ec_sm[j] = clamp(Ec_sm[i], 0, Δz[i] * (θ[i] - θ_wp)) # 限制最大蒸发量
    Ec_gw[i] = 0.0
  end

  # 半饱和
  if 1 <= j <= N
    i = j
    θ_unsat, frac_unsat = find_θ_unsat(soil, θ_sat) # 这里可能需要求，非饱和的比例（或深度）
    fsm_Ec[i], fsm_Es[i] = swc_stress(θ_unsat, pEc, soilpar, pftpar)

    _Ec_pot_sat = Ec_pot[i] * (1 - frac_unsat) # GW
    _Ec_pot_unsat = Ec_pot[i] * frac_unsat     # SM

    Ec_sm[i] = _Ec_pot_unsat * f_cons * fsm_Ec[i]
    Ec_sm[i] = clamp(Ec_sm[i], 0, Δz[i] * (θ_unsat - θ_wp) * frac_unsat) # 限制最大蒸发量

    Ec_gw[i] = _Ec_pot_sat * f_cons
  end

  # 全部饱和，不受水分限制
  for i = j+1:N
    fsm_Ec[i] = 1.0
    fsm_Es[i] = 1.0

    Ec_sm[i] = 0.0
    Ec_gw[i] = Ec_pot[i] * f_cons
  end

  Es = max(fsm_Es[1] * pEs, 0.0)
  Tr = sum(Ec_sm) + sum(Ec_gw)

  # 土壤蒸发中的[地下水蒸发]，处理的不妥当
  if j == 1
    d1 = zwt
    θ_unsat, frac_unsat = find_θ_unsat(soil, θ_sat) # 这里可能需要求，非饱和的比例（或深度）
    Tr1_u = Ec_sm[1]
    Es_u = clamp(Es * frac_unsat, 0, d1 * θ_unsat) # TODO: update Es?
    sink[1] = Tr1_u + Es_u
  else
    Tr1 = Ec_sm[1] + Ec_gw[1]
    # TODO: 第一层能否蒸发掉所有水？
    if θ[1] > 0 && Es + Tr1 > Δz[1] * θ[1]
      Tr1 = Δz[1] * θ[1] * Tr1 / (Tr1 + Es) # update Tr1 and Tr
      Es = Δz[1] * θ[1] - Tr1
    end
    sink[1] = Tr1 + Es
  end

  ## 对第一层蒸发进行限制和调整
  return Tr, Es
end
