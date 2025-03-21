function cal_specific_yield(soil::Soil)
  (; N, z₊ₕ, Δz, zwt, θ) = soil
  (; θ_sat, θ_fc) = soil.soilpar

  jwt = find_jwt(z₊ₕ, zwt) # 地下水所在层  

  jwt == 0 && return θ_sat - θ_fc
  jwt == N + 1 && return 0.2 # Kun, 2019

  ∑ = 0.0    # 孔隙度
  ∑_d = 0.0  # depth
  
  for i = 1:jwt
    z0 = i == 1 ? 0.0 : z₊ₕ[i-1]
    z1 = z₊ₕ[i]

    depth = i == jwt ? zwt - z0 : Δz[i]
    ∑_d += depth

    _θ = i == jwt ? (θ[i] * Δz[i] - θ_sat * (z1 - zwt)) / depth : θ[i]
    ∑ += max(θ_sat - _θ, 0.0) * depth
  end
  sy = ∑ / ∑_d
  return clamp(sy, 0.0, 0.3) # sy, 对每一层的孔隙度进行加权。
end

# wa: 地下水量，[m]
# zwt: 地下水水位，[m]
# Δt: 时间步长，[s]
# recharge: 补给量，[mm s-1]
# GW_Update_zwtθ!(soil, zwt, wa, recharge; θ=nothing, wa_max=5000.0)
function Update_zwt_theta!(soil::Soil{T}, zwt, wa, recharge;
  θ::Union{AbstractVector{T},Nothing}=nothing,
  wa_max=5000.0) where {T<:Real}

  update_θ = !isnothing(θ)
  constrain_wa = wa_max > 0

  (; N, z₊ₕ, Δz, Sy, Syₙ) = soil
  (; θ_sat, θ_wp) = soil.soilpar

  jwt = find_jwt(z₊ₕ, zwt) # 地下水所在层
  ∑ = recharge / 1000 # [mm s-1] to [m]

  wa = wa + ∑ * 1000 # [m]
  uex = 0.0   # excess water to soil surface，多余的水产生地表径流

  # 这种排列，已经考虑了`jwt=N`的情况
  if ∑ > 0
    # 半饱和 + 非饱和；水位上升，从下至上补给
    for i = jwt:-1:1
      # _sy = i == N + 1 ? Sy[N] : Sy[i] # unitless
      _sy = cal_specific_yield(soil)
      z0 = i == 1 ? 0.0 : z₊ₕ[i-1]

      # 补给量，正值
      if update_θ
        left = i == N + 1 ? _sy * (zwt - z0) :
               max(min((θ_sat - θ[i]) * Δz[i], _sy * (zwt - z0)), 0.0)
        _recharge = clamp(∑, 0, left)
        i <= N && (θ[i] = θ[i] + _recharge / Δz[i])
      else
        _recharge = clamp(∑, 0, _sy * (zwt - z0))
      end
      # 如果补给增加，则水位应该上升，无需限制zwt，因为上界是_z₊ₕ
      zwt = zwt - _recharge / _sy
      ∑ -= _recharge
      ∑ <= 0 && break
    end
    ∑ > 0 && (uex = ∑ * 1000) # excess water to soil surface, [m]
  else
    # 水位下降，排泄
    for i = jwt:N
      _sy = cal_specific_yield(soil)
      _θ_wp = i == 1 ? 0.01 : θ_wp # 土壤蒸发可超越凋萎含水量的限制
      @show i, _sy

      if update_θ
        left = i == 0 ? _sy * (0 - zwt) :
               max(min((θ[i] - _θ_wp) * Δz[i], _sy * (z₊ₕ[i] - zwt)), 0.0)
        _discharge = clamp(∑, -left, 0) # 排泄量，负值
        
        i >= 1 && (θ[i] = θ[i] + _discharge / Δz[i])
      else
        _discharge = clamp(∑, -_sy * (z₊ₕ[i] - zwt), 0) # 排泄量，负值
      end
      zwt = zwt - _discharge / _sy # 如果排泄减少，则水位应该下降，水位最小为z₊ₕ[i]
      ∑ -= _discharge
      ∑ >= 0 && break # 只可能等于0，不可能大于0
    end
    # 到最后一层，仍然有剩余排泄能力
    ∑ < 0 && (zwt = zwt - ∑ / Syₙ)
  end

  if constrain_wa && wa > wa_max
    # 超出的部分，放到土壤水中，[10^3 kg m-2], [m^3 m-2]
    wa = wa_max
    update_θ && (θ[N] = θ[N] + (wa - wa_max) / 1000 / Δz[N])
  end

  @pack! soil = zwt, wa, uex
  (; zwt, wa, uex)
end


"""
drainage：

排出为负，同时更新土壤含水量θ

# Arguments
- `update_θ`: 是否更新θ
- `wa_max`: if `wa_max` = 0 or NaN, then no constraint on `wa`, [mm].
"""
function GW_UpdateDrainage!(soil::Soil{T}, zwt, wa, drainage;
  θ::Union{AbstractVector{T},Nothing}=nothing,
  wa_max=5000.0) where {T<:Real}

  update_θ = !isnothing(θ)
  constrain_wa = wa_max > 0

  (; N, z₊ₕ, Δz, Sy) = soil
  jwt = find_jwt(z₊ₕ, zwt)
  ∑ = -drainage / 1000 # 负值, [mm] to [m]
  wa = wa + ∑ * 1000

  for j = jwt:N # 混合 + 饱和带
    # 补充给水度的计算方法
    _sy = Sy[j] # unitless
    _drainage = clamp(∑, -_sy * (z₊ₕ[j] - zwt), 0) # 排泄量，负值, [m]
    update_θ && (θ[j] = θ[j] + _drainage / Δz[j])

    zwt = zwt - _drainage / _sy
    ∑ -= _drainage
    ∑ >= 0 && break # 只可能等于0，不可能大于0
  end

  if ∑ < 0
    zwt = zwt - ∑ / Sy[N]
    if constrain_wa && wa > wa_max
      # 由于∑<0，这里应该不会生效
      # 超出的部分，放到土壤水中，[10^3 kg m-2], [m^3 m-2]
      wa = wa_max
      update_θ && (θ[N] = θ[N] + (wa - wa_max) / 1000 / Δz[N])
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
