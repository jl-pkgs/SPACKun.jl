# 根据补给，更新地下水水位
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
