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
function potentialET(Rn, G, LAI, Ta, Pa)
  k = 0.6  # the empirical extinction coefficient set as 0.6
  alpha = 1.26  # PT coefficient for water saturated surface

  # Radiation located into soil and canopy, separately
  Rns = exp(-k * LAI) * Rn
  Rnc = Rn - Rns

  es = 0.6108 * exp((17.27 * Ta) / (Ta + 237.3)) # [kPa]
  λ = 2.501e6 - 2361 * Ta              # [J/kg]  
  Δ = (4098 * es) / ((Ta + 237.3)^2)   # [kPa/degC]
  γ = (Cp * Pa) / (ϵ * λ)              # Psychrometric constant (kPa/degC)

  # Potential Transpiration and Soil evaporation, mm/day
  pEc = ((Rnc * alpha * Δ / (Δ + γ)) / λ) * 24 * 3600
  pEs = (((Rns - G) * alpha * Δ / (Δ + γ)) / λ) * 24 * 3600

  pEc = max(pEc, 0)
  pEs = max(pEs, 0)
  return pEc, pEs
end
