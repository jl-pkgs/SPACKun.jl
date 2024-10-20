const Cp = 1013  # Specific heat (J kg-1 C-1)
const ϵ = 0.622  # e (unitless) is the ratio of molecular weight of water to dry air (equal to 0.622)

"""
Potential ET partition

## INPUT:
- Ta    : air temperature, C
- Rn    : average daily net radiation, W/m^2
- Pa    : atmospheric pressure, kPa
- LAI   : Leaf area index, 1
- G     : Soil heat flux, W/m^2

## OUTPUT
- pEc   : potential Transpiration, mm/day
- pEs   : potential Soil evaporation, mm/day
"""
function potentialET(Rn::T, G::T, LAI::T, Ta::T, Pa::T) where {T<:Real}
  k = 0.6  # the empirical extinction coefficient set as 0.6
  # Radiation located into soil and canopy, separately
  Rns = exp(-k * LAI) * Rn
  Rnc = Rn - Rns

  λ = cal_lambda(Ta)                   # [J/kg]  
  Δ = cal_slope(Ta)                    # [kPa/degC]
  γ = (Cp * Pa) / (ϵ * λ)              # Psychrometric constant (kPa/degC)

  # Potential Transpiration and Soil evaporation, mm/day
  pEc = ET0_PT1972(Rnc, Δ, γ, λ)
  pEs = ET0_PT1972(Rns - G, Δ, γ, λ)
  return pEc, pEs
end


# Rn: W m-2
function ET0_PT1972(Rn::T, Δ::T, γ::T, λ::T) where {T<:Real}
  α = 1.26  # PT coefficient for water saturated surface
  ET0 = α * Rn * Δ / (Δ + γ)
  ET0 = max(ET0, 0)
  W2mm(ET0, λ)
end

# λ: latent heat of vaporization (J kg-1)
W2mm(w::T, λ::T) where {T<:Real} = w * 86400 / λ

# λ: latent heat of vaporization (J kg-1)
cal_lambda(Ta::T) where {T<:Real} = 2.501e6 - 2361 * Ta

# cal_gamma(Pa, λ) = (Cp * Pa) / (ϵ * λ)
cal_es(Ta::T) where {T<:Real} = 0.6108 * exp((17.27 * Ta) / (Ta + 237.3))

cal_slope(Ta::T) where {T<:Real} = 4098 * cal_es(Ta) / ((Ta + 237.3)^2)
