# drainage < 0, [m s-1]
function GW_Correctθ!(soil::Soil{T}, θ::AbstractVector{T}, zwt, wa, Δt, drainage; θ_sat=0.4) where {T<:Real}
  (; N, Δz, Sy) = soil
  θ_sat = fill(θ_sat, N)
  # (; θ_sat) = soil.param

  # jwt = find_jwt(z₊ₕ, zwt)
  zwt = clamp(zwt, 0.0, 80.0) # 地下水水位在[0, 80m]
  # 强制限制水位，不考虑水量平衡是否合适？

  ## 1. 超饱和？
  for j = N:-1:2
    exceed = max((θ[j] - θ_sat[j]) * Δz[j], 0.0)
    if exceed > 0.0
      θ[j] = θ_sat[j]
      θ[j-1] = θ[j-1] + exceed / Δz[j-1]
    end
  end

  uex = 0.0
  j = 1
  exceed = max((θ[j] - θ_sat[j]) * Δz[j], 0.0)
  if exceed > 0.0
    θ[j] = θ_sat[j] # 第一层如果超了，则产生地表径流
    uex += exceed * 1000
  end

  ## 2. 亏损？
  ∑_neg = 0.0
  for j = 1:N
    if θ[j] < 0
      ∑_neg += θ[j] * Δz[j]
      θ[j] = 0.0
    end
  end

  # 1. 少排一点水; 2. drainage扣完之后，从地下水中扣除
  drainage = drainage + ∑_neg / Δt * 1000 # [mm s-1] or [kg s-1]
  if drainage < 0
    wa = wa + drainage * Δt
    zwt = zwt - drainage * Δt / Sy[N] / 1000 # 地下水变深
    drainage = 0.0
  end
  @pack! soil = zwt, wa, uex, drainage
  (; zwt, wa, uex, drainage)
end
