# s_vod * s_tem
function Evapotranspiration!(soil::Soil, pEc::T, pEs::T, fwet::T, f_cons::T) where {T<:Real}
  (; θ, Ec_pot, fsm_Ec, fsm_Es,
    Ec_sm, Ec_gw, sink,
    zwt, z₊ₕ, N, param) = soil
  (; θ_sat, θ_fc, θ_wp, hc) = param

  PT_partition!(soil, pEc, fwet) # EC_POT划分

  j = find_jwt(z₊ₕ, zwt)
  _θ_unsat, frac_unsat = find_θ_unsat(soil, θ_sat[1]) 
  # TODO: fix `θ_sat`; 这里可能需要求，非饱和的比例（或深度）

  # 非饱和
  for i = 1:max(j - 1, 0) # 全部非饱和
    fsm_Ec[i], fsm_Es[i] = swc_stress(θ[i], pEc, θ_fc[i], θ_wp[i], hc)

    Ec_sm[i] = Ec_pot[i] * f_cons * fsm_Ec[i]
    # Ec_sm[j] = clamp(Ec_sm[i], 0, Δz[i] * (θ[i] - θ_wp)) # 限制最大蒸发量
    Ec_gw[i] = 0.0
  end

  # 半饱和
  if 1 <= j <= N
    i = j
    fsm_Ec[i], fsm_Es[i] = swc_stress(_θ_unsat, pEc, θ_fc[i], θ_wp[i], hc)

    _Ec_pot_sat = Ec_pot[i] * (1 - frac_unsat) # GW
    _Ec_pot_unsat = Ec_pot[i] * frac_unsat     # SM

    Ec_sm[i] = _Ec_pot_unsat * f_cons * fsm_Ec[i]
    Ec_sm[i] = clamp(Ec_sm[i], 0, max(Δz[i] * (_θ_unsat - θ_wp[i]) * frac_unsat, 0)) # 限制最大蒸发量

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
  
  ## 分离出土壤蒸发 -> sink 只考虑非饱和部分的土壤蒸发
  if j == 1
    d1 = zwt
    Tr1_u = Ec_sm[1]
    Es_u = clamp(Es * frac_unsat, 0, d1 * _θ_unsat) # TODO: update Es?
    sink[1] = Tr1_u + Es_u
  else
    Tr1 = Ec_sm[1] + Ec_gw[1]
    # TODO: 第一层能蒸发掉所有水？
    if θ[1] > 0 && Es + Tr1 > Δz[1] * θ[1]
      Tr1 = Δz[1] * θ[1] * Tr1 / (Tr1 + Es) # update Tr1 and Tr
      Es = Δz[1] * θ[1] - Tr1
    end
    sink[1] = Tr1 + Es
  end

  for i = 2:N
    _θ = i == j ? _θ_unsat : θ[i]
    depth = i == j ? zwt - z₊ₕ[i-1] : Δz[i]
    sink[i] = clamp(Ec_sm[i], 0, depth * (_θ - θ_wp[i]))
  end
  return Tr, Es
end
