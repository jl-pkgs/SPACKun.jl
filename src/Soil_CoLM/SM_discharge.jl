"""
    SM_discharge!(soil, soilpar; Q_in=0.0)

> 按照CoLM的方式进行修改

Kun的做法是合理的，非饱和带和饱和带应该分开，否则地下水层会出现镂空。

该函数仅处理非饱和带的土壤移动，非饱和带土壤水分的变化不会引起地下水水位的变化。
饱和带的水分移动，会引起地下水水位的变化。

# TODO: 
1. `soil_drainage`存在bug，土壤的厚度影响排泄量, 需要重新率定参数，以解决该问题；
"""
function SM_discharge!(soil::Soil, soilpar; Q_in::T=0.0) where {T<:Real}
  (; z₊ₕ, Δz, zwt, θ, Q, Dmin, Dmax, N) = soil
  (; θ_sat, θ_wp, Ksat) = soilpar

  j = find_jwt(z₊ₕ, zwt)
  d1, d2 = 0.0, 0.0

  ## 非饱和带
  for i = 1:j
    _θ_wp = i == 1 ? 0.01 : θ_wp # 土壤蒸发可超越凋萎含水量的限制
    i == N + 1 && continue
    z0 = i == 1 ? 0 : z₊ₕ[i-1]
    z1 = i == j ? zwt : z₊ₕ[i]
    depth = z1 - z0

    _θ = θ[i]
    if i == j
      d1 = zwt - z0 # 非饱和带的厚度
      d2 = z1 - zwt # 饱和带的厚度
      _θ_unsat = (θ[j] * Δz[j] - θ_sat * d2) / d1
      _θ = _θ_unsat
    end

    # Q = K (∂ψ/∂z + 1)
    Q[i] = soil_drainage(_θ, θ_sat, Ksat, Dmin[i], Dmax[i]) # 向下排泄量, [mm/d]
    ΔQ = Q_in - Q[i]

    # _sink = sink[i]
    _sink = clamp(sink[i], 0, (_θ - _θ_wp) * depth)         # 对蒸腾量进行限制
    _θ = clamp((_θ * depth + ΔQ - _sink) / depth, 0.0, 1.0) # 更新含水量

    Q_in = Q[i] # 一下次迭代
    if _θ > θ_sat
      _θ = θ_sat
      Q_in += (θ[i] - θ_sat) * depth
    elseif _θ < _θ_wp
      _θ = _θ_wp
      Q_in += (θ[i] - _θ_wp) * depth
    end

    ## θ_unsat -> θ
    i == j && (_θ = (_θ * d1 + θ_sat * d2) / Δz[j])
    θ[i] = _θ
  end
  ## 在排泄过程中，水位会下降，也会上升
  return Q_in # exceed
end
