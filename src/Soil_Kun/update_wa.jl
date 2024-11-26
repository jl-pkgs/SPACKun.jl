function update_wa!(soil::Soil, θ_unsat, zwt1, zwt2)
  (; θ, z₊ₕ, N) = soil
  (; θ_sat, θ_fc) = soil.soilpar
  dz = zwt2 - zwt1

  j0 = find_jwt(z₊ₕ, zwt1) # 地下水所在层
  for i = 1:N
    i != j0 && (θ[i] = θ_unsat[i])
  end

  j1 = find_jwt(z₊ₕ, zwt2) # 地下水所在层  
  uex = 0.0

  if dz <= 0 # 水位上升
    for i = min(j0, N):-1:j1
      z0 = i == 1 ? 0 : z₊ₕ[i-1]
      if i <= 0
        uex = -zwt2 * θ_fc
        break
      end

      if i == j1
        θ[i] = (θ_unsat[i] * (zwt2 - z0) + θ_sat * (z₊ₕ[i] - zwt2)) / Δz[i]
      else
        θ[i] = θ_sat
      end
    end

  else # 水位下降
    # 逻辑完整
    i = j0

    if i >= 1
      z0 = i == 1 ? 0 : z₊ₕ[i-1]
      if j0 == j1
        θ[i] = (θ_unsat[i] * (zwt2 - z0) + θ_sat * (z₊ₕ[i] - zwt2)) / Δz[i]
      else
        θ[i] = (θ_unsat[i] * (zwt1 - z0) + θ_fc * (z₊ₕ[i] - zwt1)) / Δz[i]
      end
    end

    for i = j0+1:j1
      z0 = i == 1 ? 0 : z₊ₕ[i-1]
      if i == j1
        θ[i] = (θ_fc * (zwt2 - z0) + θ_sat * (z₊ₕ[i] - zwt2)) / Δz[i]
      else
        θ[i] = θ_fc
      end
    end
  end

  return uex
end
