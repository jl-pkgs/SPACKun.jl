export sw_balance
export find_jwt, GW_Rsb, SM_recharge!
export soil_drainage

include("find_θ_unsat.jl")
include("SM_discharge.jl")
include("SM_recharge.jl")
include("update_wa.jl")
include("swc_stress.jl")


# 水位向下为正，地表为0
function find_jwt(z₊ₕ::AbstractVector, zwt::Real)
  N = length(z₊ₕ)
  zwt <= 0 && return 0
  zwt >= z₊ₕ[end] && return N + 1

  for j in 1:N
    zwt <= z₊ₕ[j] && return j
  end
end

function GW_Rsb(zwt::Real)
  R_sb_max = 39.0 # mm day-1
  f = 1.25e-3     # mm-1
  return R_sb_max * exp(-f * zwt)
end


"""
    soil_drainage(se, ks, Dmin, Dmax; dd = 1.5)

```julia
# layer1
Dmin = 0.048; # mm day-1
Dmax = 4.8; # mm day-1

# layer2
Dmin = 0.012; # 0.0005*24, mm day-1
Dmax = 1.2; # 0.05*24,   mm day-1
```

# Reference
1. Ducharne & Polcher, 1998
"""
function soil_drainage(se, Ksat, Dmin, Dmax; dd=1.5)
  # se = θ_unsat / θ_sat
  if se < 0.75
    perc = Dmin * se
  else
    perc = Dmin * se + (Dmax - Dmin) * se^dd
  end
  return min(Ksat, perc) # drainage from unsaturated zone
end


"""
    sw_balance(I::T, pEc::T, pEs::T, Ta::T, Topt::T, s_VOD::T,
        soilpar, pftpar, wet::::T, Δz::::T, θ, zwt::T) where {T<:Real}

# INPUT
- `I`      : total water enter into soil surface, mm
- `pEc`    : potential ET allocate to plant, mm
- `pEs`    : potential ET allocate to soil surface, mm
- `Ta`     : air temperature, C
- `Topt`   : optimal growth temperature for plant, C
- `s_VOD`  : constrains of VOD,[0,1]
- `θ`     : previous soil water content, 3 layers
- `soilpar`: soil-related parameters
- `pftpar` : plant-related parameters
- `wet`    : wetness fraction indice
- `Δz`     : soil layer depth, 3 layers
- `zwt`    : groundwater table depth, mm

# OUTPUT
- `Tr`     : actual plant transpiration, mm
- `Es`     : actual soil evaporation, mm
- `θ`     : updated soil water content, 3 layers, %
- `zwt`    : groundwater table depth, mm
- `uex`    : goundwater overflow soil surface, mm
"""
function sw_balance(soil::Soil, I::T, pEc::T, pEs::T, Ta::T, Topt::T, fwet, s_vod::T) where {T<:Real}
  θ_sat = soil.param.θ_sat[1]
  θ_fc = soil.param. θ_fc[1]

  s_tem = exp(-((Ta - Topt) / Topt)^2)

  (; θ_unsat, θ, N, Δz, zwt, Ec_gw, sink) = soil

  θ_unsat .= θ
  extra = SM_recharge!(θ, I; Δz, θ_sat)

  jwt = find_jwt(z₊ₕ, zwt)
  1 <= jwt <= N && (θ_unsat[jwt] = find_θ_unsat(θ, zwt; z₊ₕ, Δz, θ_sat)[1])

  # ====== Water Consumption ====== #
  f_cons = s_tem * s_vod
  Tr, Es = Evapotranspiration!(soil, pEc, pEs, fwet, f_cons)

  # ====== Soil Water Drainage (Unsaturated Zone) ====== #  
  exceed = SM_discharge!(soil, θ_unsat, sink)

  # ====== The Groundwater Table Depth ====== #
  Δw = exceed + extra - sum(Ec_gw) - GW_Rsb(zwt)

  1 <= jwt <= N && (sy = θ_sat - mean(θ_unsat[1:jwt]))
  jwt == 0 && (sy = θ_sat - θ_fc)
  jwt > N && (sy = 0.2)

  zwt = zwt - Δw / sy

  uex = update_wa!(soil, θ_unsat, soil.zwt, zwt)
  soil.zwt = max(0, zwt)
  return Tr, Es, uex
end
