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
# θ_unsat, frac = find_θ_unsat(soil)
function find_θ_unsat(soil::Soil, θ_sat=NaN)
  (; θ, zwt, z₊ₕ, Δz) = soil
  isnan(θ_sat) && (θ_sat = soil.θ_sat)
  find_θ_unsat(θ, zwt; z₊ₕ, Δz, θ_sat)
end
