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
- β   : 土壤水分限制因子, [0~1], 1充分供水

## Optional:
- kA  : 消光系数，default `[0.6, 0.6, 0.6]`
- Kc  : 作物折腾转换系数，default `[0.75, 0.9, 1.0]`

## 优化参数
- `rs`: 

## 不敏感参数
- `hc`: 指定参数，小麦和玉米，0.4m, 1.0m；
- `Kc`: 需要测试

## OUTPUT
- pEc   : potential Transpiration, mm/day
- pEs   : potential Soil evaporation, mm/day
"""
function cal_PET(Rn::FT, G::FT, LAI::FT, Ta::FT, Pa::FT, VPD::FT, U2::FT, doy::Int;
  method="PT1972", β=1.0, 
  α_soil=1.26,
  # Kc::Vector{FT}=[0.75, 0.9, 1.0],
  Hc::Vector{FT}=[0.12, 0.12, 0.12],
  rs::Vector{FT}=[70.0, 70.0, 70.0],
  kA::Vector{FT}=[0.6, 0.6, 0.6]) where {FT<:Real}

  # 考虑不同作物，rs, hc, kA的不同；这里引入了rs，因此不再需要Kc
  i = iGrowthPeriod(doy)
  _kA = kA[i]
  _rs = rs[i]
  hc = Hc[i]

  # Radiation located into soil and canopy, separately
  Rns = exp(-_kA * LAI) * Rn
  Rnc = Rn - Rns

  ## Potential Transpiration and Soil evaporation, mm/day
  # ET0 = ET0_PT1972(Rn, Ta, Pa)
  # ET0 = ET0_Penman48(Rn, Ta, VPD, U2, Pa) # all

  ## canopy
  if method == "Penman48"
    pEc = ET0_Penman48(Rnc, Ta, VPD, U2, Pa)
  elseif method == "PT1972"
    pEc = ET0_PT1972(Rnc, Ta, Pa)
  elseif method == "Monteith65"
    pEc = ET0_Monteith65(Rnc, Ta, VPD, U2, Pa, β; rs=_rs, hc)
  end
  pEs = ET0_eq(Rns - G, Ta, Pa)[1] * α_soil # PML
  return pEc, pEs
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


"""
    ET0_FAO98(Rn::T, Tair::T, VPD::T, wind::T, Pa::T=atm; z_wind=2, tall_crop=false)
# Examples
```julia
ET0_Penman48(200., 20., 2., 2.)
PET = ET0_FAO98(200.0, 20.0, 2.0, 2.0) # mm
PET = ET0_FAO98(Rn, Tair, VPD, U2, Pa; tall_crop = false) # mm
```
"""
function ET0_FAO98(Rn::T, Tair::T, VPD::T, wind::T, Pa::T=atm; z_wind=2, tall_crop=false) where {T<:Real}
  # T = eltype(Rn)
  U2 = cal_U2(wind, z_wind)

  if tall_crop
    p1 = T(1600.0)
    p2 = T(0.38)
  else
    p1 = T(900.0)
    p2 = T(0.34)
  end
  Eeq, λ, Δ, γ = ET0_eq(Rn, Tair, Pa)
  Eeq::T = Δ / (Δ + (γ * (1.0 + p2 * U2))) * Rn |> x -> W2mm(x, λ)
  Evp::T = γ * p1 / (Tair + 273.15) * U2 * VPD / (Δ + (γ * (1.0 + p2 * U2)))
  Eeq + Evp
end


# β: 土壤水分限制因子, [0~1], 1充分供水。
function ET0_Monteith65(Rn::T, Tair::T, VPD::T, wind::T, Pa::T=atm, β::T=1.0;
  z_wind=2.0, rs=70.0, hc=0.12) where {T<:Real}

  Eeq, λ, Δ, γ = ET0_eq(Rn, Tair, Pa)
  U2 = cal_U2(wind, z_wind)
  rH = aerodynamic_resistance(U2, hc)
  rs = rs / β

  rho_a = cal_rho_a(Tair, Pa) # FAO56, Eq. 3-5, kg m-3
  # Cp = 1.013 * 1e-3 # MJ kg-1 degC-1
  # rho_a * Cp * dT * gH (in MJ m-2 s-1)
  # = kg m-3 * MJ kg-1 degC-1 * degC * m s-1
  # = MJ m-2 s-1
  Eeq = Δ / (Δ + (γ * (1 + rs / rH))) * Rn |> x -> W2mm(x, λ)
  Evp = (rho_a * Cp * VPD / rH) / (Δ + (γ * (1 + rs / rH))) * 86400 / λ
  Eeq + Evp
end


# 空气动力学阻力
function aerodynamic_resistance(U2::T, h::T=0.12) where {T<:Real}
  k = 0.41
  d = 2 / 3 * h
  z_om = 0.123 * h
  z_oh = 0.1 * z_om
  z_m = z_h = 2
  log((z_m - d) / z_om) * log((z_h - d) / z_oh) / (k^2 * U2)
end

# λ: latent heat of vaporization (MJ kg-1)
W2mm(w::T, λ::T) where {T<:Real} = w * 0.086400 / λ

# λ: latent heat of vaporization (MJ kg-1)
cal_lambda(Ta::T) where {T<:Real} = 2.501 - 2.361 / 1000 * Ta

# cal_gamma(Pa, λ) = (Cp * Pa) / (ϵ * λ)
cal_es(Ta::T) where {T<:Real} = 0.6108 * exp((17.27 * Ta) / (Ta + 237.3))

cal_slope(Ta::T) where {T<:Real} = 4098 * cal_es(Ta) / ((Ta + 237.3)^2)

function cal_U2(Uz::T, z=10.0) where {T<:Real}
  z == 2 && (return Uz)
  Uz * 4.87 / log(67.8 * z - 5.42)
end

# kg m-3
cal_rho_a(Tair, Pa) = 3.486 * Pa / cal_TvK(Tair)

## 几种计算虚温的方法

# FAO56, Eq. 3-7
cal_TvK(Tair) = 1.01 * (Tair + 273)

# 这个是最精确的版本
# FAO56, Eq. 3-6
cal_TvK(Tair, ea, Pa) = (Tair + K0) * (1 + (1 - epsilon) * ea / Pa)

# https://github.com/CUG-hydro/class2022_CUG_HydroMet/blob/master/ch02_大气的基本特征.md
# q ≈ ϵ*ea/Pa
# q = ϵ*ea/(Pa - (1 - ϵ)*ea)
function cal_TvK(Tair, q)
  # ea / Pa = q/(ϵ + (1 - ϵ) * q)
  (Tair + K0) * (1 + (1 - epsilon) * q / (ϵ + (1 - ϵ) * q))
end


"""
Partition transpiration into soil layers

# INPUT:
- pEc    :  Potential Evaporation on canopy
- w      :  Initialized values for soil moisture
- pftpar :  PFT parameters

# OUTPUT:
- Tr_p :  separate potential Transpiration
"""
function PT_partition!(soil::Soil, pEc::T, fwet::T) where {T<:Real}
  (; θ, Δz, Ec_pot, param) = soil
  (; D50, D95) = param
  b = param.b[1]
  θ_sat = param.θ_sat[1]
  
  c = -2.944 / log(D95 / D50)
  r1 = (1 / (1 + (Δz[1] / D50)^c)) # Zhang 2019, Eq. 21, root depths function
  r2 = (1 / (1 + (Δz[2] / D50)^c)) - (1 / (1 + (Δz[1] / D50)^c))
  r3 = (1 / (1 + (Δz[3] / D50)^c)) - (1 / (1 + (Δz[2] / D50)^c))

  # the maximum transpiration rate of each soil layer, Tr_p
  # Get all available water contents through root distribution
  wr = r1 * (θ[1] / θ_sat)^b + r2 * (θ[2] / θ_sat)^b + r3 * (θ[3] / θ_sat)^b

  # Root distribution adjusted by soil water content
  β1 = r1 * (θ[1] / θ_sat)^b / wr
  β2 = r2 * (θ[2] / θ_sat)^b / wr
  β3 = r3 * (θ[3] / θ_sat)^b / wr

  # potential transpiration rate for different layers
  Ec_pot[1] = (1 - fwet) * β1 * pEc
  Ec_pot[2] = (1 - fwet) * β2 * pEc
  Ec_pot[3] = (1 - fwet) * β3 * pEc
end


export aerodynamic_resistance
export cal_rho_a
export ET0_eq, ET0_Penman48, ET0_Monteith65, ET0_PT1972, ET0_FAO98
