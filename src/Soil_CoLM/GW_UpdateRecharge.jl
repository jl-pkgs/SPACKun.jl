# wa: 地下水量，[m]
# zwt: 地下水水位，[m]
# Δt: 时间步长，[s]
# recharge: 补给量，[mm s-1]
function GW_UpdateRecharge!(soil::Soil{T}, wa, zwt, Δt, recharge) where {T<:Real}
  (; N, z₊ₕ, Sy) = soil
  jwt = find_jwt(z₊ₕ, zwt) # 地下水所在层
  ∑ = recharge * Δt / 1000 # [mm s-1] to [m]

  wa = wa + ∑ * 1000 # [m]
  uex = 0.0   # excess water to soil surface，多余的水产生地表径流

  # 这种排列，已经考虑了`jwt=N`的情况
  if ∑ > 0
    # 水位上升，从下至上补给
    # 半饱和 + 非饱和
    for i = jwt:-1:1
      _sy = i == N + 1 ? Sy[N] : Sy[i] # unitless
      z0 = i == 1 ? 0.0 : z₊ₕ[i-1]
      _recharge = clamp(∑, 0, -_sy * (z0 - zwt)) # 补给量，正值
      
      zwt = zwt - _recharge / _sy # 如果补给增加，则水位应该上升，无需限制zwt，因为上界是_z₊ₕ
      ∑ -= _recharge
      ∑ <= 0 && break
    end
    ∑ > 0 && (uex = ∑ * 1000) # excess water to soil surface, [m]
  else
    # 水位下降
    for i = jwt:N
      _sy = Sy[i] # unitless
      _discharge = clamp(∑, -_sy * (z₊ₕ[i] - zwt), 0) # 排泄量，负值

      zwt = zwt - _discharge / _sy # 如果排泄减少，则水位应该下降，水位最小为z₊ₕ[i]
      ∑ -= _discharge
      ∑ >= 0 && break # 只可能等于0，不可能大于0
    end
    # 到最后一层，仍然有剩余排泄能力
    _sy = Sy[N]
    ∑ < 0 && (zwt = zwt - ∑ / Sy[N])
  end
  @pack! soil = zwt, wa, uex
  (; zwt, wa, uex)
end



# drainage：排出为负，同时更新土壤含水量θ
function GW_UpdateDrainage!(soil::Soil{T}, θ::AbstractVector{T}, zwt, wa, Δt, drainage) where {T<:Real}
  (; N, z₊ₕ, Δz, Sy) = soil
  jwt = find_jwt(z₊ₕ, zwt)
  ∑ = -drainage * Δt / 1000 # 负值, [mm] to [m]
  wa = wa + ∑ * 1000

  for j = jwt:N # 混合 + 饱和带
    # 补充给水度的计算方法
    
    _sy = Sy[j] # unitless
    _drainage = clamp(∑, -_sy * (z₊ₕ[j] - zwt), 0) # 排泄量，负值, [m]
    θ[j] = θ[j] + _drainage / Δz[j]

    zwt = zwt - _drainage / _sy
    ∑ -= _drainage
    ∑ >= 0 && break # 只可能等于0，不可能大于0
  end

  if ∑ < 0
    zwt = zwt - ∑ / Sy[N]
    if wa > 5000.0
      # 由于∑<0，这里应该不会生效
      θ[N] = θ[N] + (wa / 1000 - 5.0) / Δz[N] # 超出的部分，放到土壤水中，[10^3 kg m-2], [m^3 m-2]
    end
  end
  @pack! soil = zwt, wa
  (; zwt, wa)
end



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
