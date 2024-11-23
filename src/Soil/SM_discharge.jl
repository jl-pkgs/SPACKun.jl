# SM_discharge!(soil, θ_unsat, soilpar; Q_in=0.0)
function SM_discharge!(soil::Soil, θ_unsat::Vector{T}, sink, 
  soilpar; Q_in::T=0.0) where {T<:Real}
  (; z₊ₕ, zwt, Dmin, Dmax, N) = soil
  (; θ_sat, θ_wp, Ksat) = soilpar
  j = find_jwt(z₊ₕ, zwt)

  ## 非饱和带
  for i = 1:j
    i == N + 1 && continue

    z0 = i == 1 ? 0 : z₊ₕ[i-1]
    z1 = i == j ? zwt : z₊ₕ[i]
    depth = z1 - z0

    # TODO: 这里存在bug，土壤的厚度影响排泄量
    Q_out = soil_drainage(θ_unsat[i], θ_sat, Ksat, Dmin[i], Dmax[i]) # 向下排泄量, [mm/d]
    ΔQ = Q_in - Q_out

    # 对蒸发量也进行限制
    sink[i] = clamp(sink[i], 0, (θ_unsat[i] - θ_wp) * depth)

    θ_unsat[i] = (θ_unsat[i] * depth + ΔQ - sink[i]) / depth
    θ_unsat[i] = clamp(θ_unsat[i], 0, 1)

    Q_in = Q_out # 一下次迭代
    if θ_unsat[i] > θ_sat
      Q_in += (θ_unsat[i] - θ_sat) * depth
      θ_unsat[i] = θ_sat
    end
  end

  return Q_in # exceed
end
