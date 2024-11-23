## 参考CoLM的方案，考虑地下水的补给与排泄，只考虑地下水垂直方向与土壤水的交互作用。
# 从下至上，不受地下水影响的第一层

# drainage：排出为负，同时更新土壤含水量θ
function GW_UpdateDrainage!(soil::Soil{T}, θ::AbstractVector{T}, zwt, wa, Δt, drainage) where {T<:Real}
  (; N, z₊ₕ, Δz, Sy) = soil
  jwt = find_jwt(z₊ₕ, zwt)
  ∑ = -drainage * Δt / 1000 # 负值, [mm] to [m]
  wa = wa + ∑ * 1000

  for j = jwt:N # 混合 + 饱和带
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

