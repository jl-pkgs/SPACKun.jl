export sw_balance


include("swc_stress.jl")
include("swb_case0.jl")
include("swb_case1.jl")
include("swb_case2.jl")
include("swb_case3.jl")
include("swb_case4.jl")



"""
    soil_drainage(wa_unsat, θ_sat, ks, Dmin, Dmax; dd = 1.5)

# layer1
Dmin = 0.048; # mm day-1
Dmax = 4.8; # mm day-1

# layer2
Dmin = 0.012; # 0.0005*24, mm day-1
Dmax = 1.2; # 0.05*24,   mm day-1

# Reference
1. Ducharne & Polcher, 1998
"""
function soil_drainage(wa_unsat, θ_sat, ks, Dmin, Dmax; dd=1.5)
  thx1 = wa_unsat / θ_sat
  if thx1 < 0.75
    Perc1 = Dmin * thx1
  else
    Perc1 = Dmin * thx1 + (Dmax - Dmin) * thx1^dd
  end
  return min(ks, Perc1) # drainage from unsaturated zone
end


"""
    sw_balance(IWS::T, pEc::T, pEs::T, Ta::T, Topt::T, s_VOD::T,
        soilpar, pftpar, wet::::T, Δz::::T, θ, zwt::T) where {T<:Real}

# INPUT
- `IWS`    : total water enter into soil surface, mm
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
function sw_balance(IWS::T, pEc::T, pEs::T, Ta::T, Topt::T, s_VOD::T,
  soilpar, pftpar, wet, Δz, state::State) where {T<:Real}
  (; θ, zwt) = state
  # s_tem = Temp_stress(Topt, Ta) # Constrains of temperature
  s_tem = exp(-((Ta - Topt) / Topt)^2)

  if zwt <= 0
    fun = swb_case0 # Case 0: groundwater overflow
  elseif 0 < zwt <= Δz[1]
    fun = swb_case1 # Case 1: groundwater table in layer 1
  elseif Δz[1] < zwt <= z₊ₕ[2]
    fun = swb_case2 # Case 2: groundwater table in layer 2
  elseif z₊ₕ[2] < zwt < z₊ₕ[3]
    fun = swb_case3 # Case 3: groundwater table in layer 3
  else
    fun = swb_case4 # Case 4: groundwater table below layer 3
  end

  state.θ, state.zwt, Tr, Es, uex = fun(θ, IWS, pEc, pEs, s_tem, s_VOD, soilpar, pftpar, wet, Δz, zwt)
  return Tr, Es, uex
end

# # Temperature Constrains for plant growing
# # INPUT:
# # - Topt  : Optimum temperature for plant growing
# # - Ta    : Air temperature
# function Temp_stress(Topt::T, Ta::T) where {T<:Real}
#   return exp(-((Ta - Topt) / Topt)^2)
# end
