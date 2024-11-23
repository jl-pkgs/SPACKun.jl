"""
- 未饱和: z₊ₕ[jwt] ~ zwt
- 饱和  : zwt ~ z₊ₕ[jwt+1]
"""
function find_θ_unsat(θ, zwt; z₊ₕ, Δz, θ_sat)
  N = length(θ)
  j = find_jwt(z₊ₕ, zwt)

  j == 0 && return 0.0, 0.0         # 全部饱和
  j == N + 1 && return θ[N], 1.0    # 全部非饱和

  z0 = j == 1 ? 0 : z₊ₕ[j-1]
  z1 = z₊ₕ[j]

  d_unsat = zwt - z0
  θ_unsat = (θ[j] * Δz[j] - θ_sat * (z1 - zwt)) / d_unsat

  frac = d_unsat * θ_unsat / (θ[j] * Δz[j]) # 非饱和的比例
  # frac = d_unsat / Δz[j] # 非饱和的比例
  return θ_unsat, frac
end

function find_θ_unsat(soil::Soil, θ_sat)
  (; θ, zwt, z₊ₕ, Δz) = soil
  find_θ_unsat(θ, zwt; z₊ₕ, Δz, θ_sat)
end


function SM_recharge!(θ, I; Δz, θ_sat)
  N = length(θ)
  for i = 1:N
    if θ[i] > θ_sat # 提前处理超饱和
      I += (θ[i] - θ_sat) * Δz[i]
      θ[i] = θ_sat
    end

    _layer = clamp(I, 0, (θ_sat - θ[i]) * Δz[i]) # 剩余空间
    θ[i] += (_layer / Δz[i])
    I -= _layer
    # I <= 0 && break # 超饱和现象，不能提前退出
  end
  I # 形成到向下的通量
end


# SM_discharge!(soil, θ_unsat, soilpar; Q_in=0.0)
function SM_discharge!(soil::Soil, θ_unsat::Vector{T}, sink::Vector{T},
  soilpar; Q_in::T=0.0) where {T<:Real}
  (; z₊ₕ, zwt, θ, Q, Dmin, Dmax, N) = soil
  (; θ_sat, θ_wp, Ksat) = soilpar
  j = find_jwt(z₊ₕ, zwt)

  ## 非饱和带
  for i = 1:j
    i == N + 1 && continue

    z0 = i == 1 ? 0 : z₊ₕ[i-1]
    z1 = i == j ? zwt : z₊ₕ[i]
    depth = z1 - z0

    # TODO: 这里存在bug，土壤的厚度影响排泄量
    Q[i] = soil_drainage(θ_unsat[i], θ_sat, Ksat, Dmin[i], Dmax[i]) # 向下排泄量, [mm/d]
    ΔQ = Q_in - Q[i]

    # 对蒸发量也进行限制
    _θ = i == j ? θ_unsat[i] : θ[i]
    _sink = sink[i]
    # _sink = clamp(sink[i], 0, (_θ - θ_wp) * depth)

    θ_unsat[i] = (_θ * depth + ΔQ - _sink) / depth
    θ_unsat[i] = clamp(θ_unsat[i], 0, 1)

    Q_in = Q[i] # 一下次迭代
    if θ_unsat[i] > θ_sat
      Q_in += (θ_unsat[i] - θ_sat) * depth
      θ_unsat[i] = θ_sat
    end
  end
  return Q_in # exceed
end
