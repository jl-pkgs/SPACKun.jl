export sw_balance

# include("GroundWater.jl")
# include("SM_discharge.jl")

include("swc_stress.jl")
include("swb_case0.jl")
include("swb_case1.jl")
include("swb_case2.jl")
include("swb_case3.jl")
include("swb_case4.jl")


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
function sw_balance(I::T, pEc::T, pEs::T, Ta::T, Topt::T, s_VOD::T,
  soilpar, pftpar, wet, soil::Soil) where {T<:Real}
  (; zwt) = soil
  # s_tem = Temp_stress(Topt, Ta) # Constrains of temperature
  s_tem = exp(-((Ta - Topt) / Topt)^2)

  if zwt <= 0
    fun = swb_case0 # Case 0: groundwater overflow
  elseif 0 < zwt <= z₊ₕ[1]
    fun = swb_case1 # Case 1: groundwater table in layer 1
  elseif z₊ₕ[1] < zwt <= z₊ₕ[2]
    fun = swb_case2 # Case 2: groundwater table in layer 2
  elseif z₊ₕ[2] < zwt < z₊ₕ[3]
    fun = swb_case3 # Case 3: groundwater table in layer 3
  else
    fun = swb_case4 # Case 4: groundwater table below layer 3
  end

  return fun(I, pEc, pEs, s_tem, s_VOD, soilpar, pftpar, wet, soil) # Tr, Es, uex
end

# # Temperature Constrains for plant growing
# # INPUT:
# # - Topt  : Optimum temperature for plant growing
# # - Ta    : Air temperature
# function Temp_stress(Topt::T, Ta::T) where {T<:Real}
#   return exp(-((Ta - Topt) / Topt)^2)
# end
