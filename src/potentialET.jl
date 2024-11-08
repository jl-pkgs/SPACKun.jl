# λ: [J/kg] change to [MJ kg-1]
const Cp = 1.013 * 1e-3  # Specific heat (MJ kg-1 C-1)
const ϵ = 0.622  # e (unitless) is the ratio of molecular weight of water to dry air (equal to 0.622)

period_wheat(doy::Int) = 32 <= doy <= 166
period_corn(doy::Int) = 196 <= doy <= 288

# [wheat, corn, non-gw]
function iGrowthPeriod(doy::Int)
  i = 3
  period_wheat(doy) && (i = 1)
  period_corn(doy) && (i = 2)
  return i
end


"""
Potential ET partition

## INPUT:
- Ta  : air temperature, C
- Rn  : average daily net radiation, W/m^2
- Pa  : atmospheric pressure, kPa
- LAI : Leaf area index, 1
- G   : Soil heat flux, W/m^2

## Optional:
- k   : 消光系数，default `[0.6, 0.6, 0.6]`
- Kc  : 作物折腾转换系数，default `[0.75, 0.9, 1.0]`

## OUTPUT
- pEc   : potential Transpiration, mm/day
- pEs   : potential Soil evaporation, mm/day
"""
function potentialET(Rn::T, G::T, LAI::T, Ta::T, Pa::T, 
  VPD::T, U2::T, doy::Int;
  method="PT1972",
  Kc::Vector=[0.75, 0.9, 1.0],
  k::Vector=[0.6, 0.6, 0.6]) where {T<:Real}

  i = iGrowthPeriod(doy)
  _Kc = Kc[i]
  _k = k[i]

  # Radiation located into soil and canopy, separately
  Rns = exp(-_k * LAI) * Rn
  Rnc = Rn - Rns

  ## Potential Transpiration and Soil evaporation, mm/day
  if method == "PT1972"
    ET0 = ET0_PT1972(Rn, Ta, Pa)
    pEc = ET0_PT1972(Rnc, Ta, Pa) * _Kc
    pEs = ET0_PT1972(Rns - G, Ta, Pa)
  elseif method == "Penman48"
    ET0 = ET0_Penman48(Rn, Ta, VPD, U2)             # all
    pEc = ET0_Penman48(Rnc, Ta, VPD, U2, Pa) * _Kc  # canopy
    pEs = ET0_Penman48(Rns - G, Ta, VPD, U2, Pa)    # soil
  end
  return pEc, pEs, ET0
end


function ET0_eq(Rn::T, Tair::T, Pa::T=atm, args...) where {T<:Real}
  λ::T = cal_lambda(Tair)     # [MJ kg-1]
  Δ::T = cal_slope(Tair)      # [kPa degC-1]
  γ::T = Cp * Pa / (ϵ * λ)    # [kPa degC-1], Psychrometric constant

  Eeq::T = Δ / (Δ + γ) * Rn |> x -> W2mm(x, λ)
  Eeq = max(Eeq, 0)
  (; Eeq, λ, Δ, γ)
end

# - Rn: W m-2
# - α: PT coefficient for water saturated surface
function ET0_PT1972(Rn::T, Tair::T, Pa::T=atm; α=1.26) where {T<:Real}
  Eeq = ET0_eq(Rn, Tair, Pa)[1]
  α * Eeq # [mm d-1]
end

"""
ET0_Penman48(Rn, Tair, VPD, Uz)
"""
function ET0_Penman48(Rn::T, Tair::T, VPD::T, Uz::T, Pa::T=atm; z_wind=2) where {T<:Real}
  Eeq, λ, Δ, γ = ET0_eq(Rn, Tair, Pa)
  U2::T = cal_U2(Uz, z_wind)
  Evp::T = γ / (Δ + γ) * 6.43 * (1 + 0.536 * U2) * VPD / λ
  ET0::T = Eeq + Evp
  ET0
end


# λ: latent heat of vaporization (J kg-1)
W2mm(w::T, λ::T) where {T<:Real} = w * 0.086400 / λ

# λ: latent heat of vaporization (J kg-1)
cal_lambda(Ta::T) where {T<:Real} = 2.501 - 2.361 / 1000 * Ta

# cal_gamma(Pa, λ) = (Cp * Pa) / (ϵ * λ)
cal_es(Ta::T) where {T<:Real} = 0.6108 * exp((17.27 * Ta) / (Ta + 237.3))

cal_slope(Ta::T) where {T<:Real} = 4098 * cal_es(Ta) / ((Ta + 237.3)^2)

function cal_U2(Uz::T, z=10.0) where {T<:Real}
  z == 2 && (return Uz)
  Uz * 4.87 / log(67.8 * z - 5.42)
end


export ET0_eq, ET0_Penman48, ET0_PT1972
