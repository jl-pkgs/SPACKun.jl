"""
    sw_balance(IWS::T, pEc::T, pEs::T, Ta::T, Topt::T, s_VOD::T,
        soilpar, pftpar, wet::::T, zm::::T, wa, zgw::T) where {T<:Real}

# INPUT
- `IWS`    : total water enter into soil surface, mm
- `pEc`    : potential ET allocate to plant, mm
- `pEs`    : potential ET allocate to soil surface, mm
- `Ta`     : air temperature, C
- `Topt`   : optimal growth temperature for plant, C
- `s_VOD`  : constrains of VOD,[0,1]
- `wa`     : previous soil water content, 3 layers
- `soilpar`: soil-related parameters
- `pftpar` : plant-related parameters
- `wet`    : wetness fraction indice
- `ZM`     : soil layer depth, 3 layers
- `zgw`    : groundwater table depth, mm

# OUTPUT
- `Tr`     : actual plant transpiration, mm
- `Es`     : actual soil evaporation, mm
- `wa`     : updated soil water content, 3 layers, %
- `zgw`    : groundwater table depth, mm
- `uex`    : goundwater overflow soil surface, mm
"""
function sw_balance(IWS::T, pEc::T, pEs::T, Ta::T, Topt::T, s_VOD::T,
  soilpar, pftpar, wet, zm, wa, zgw::T) where {T<:Real}
  s_tem = Temp_stress(Topt, Ta) # Constrains of temperature

  if zgw <= 0
    fun = swb_case0 # Case 0: groundwater overflow
  elseif 0 < zgw <= zm[1]
    fun = swb_case1 # Case 1: groundwater table in layer 1
  elseif zm[1] < zgw <= zm[1] + zm[2]
    fun = swb_case2 # Case 2: groundwater table in layer 2
  elseif zm[1] + zm[2] < zgw < zm[1] + zm[2] + zm[3]
    fun = swb_case3 # Case 3: groundwater table in layer 3
  else
    fun = swb_case4 # Case 4: groundwater table below layer 3
  end

  wa, zgw, Tr, Es, uex = fun(wa, IWS, pEc, pEs, s_tem, s_VOD, soilpar, pftpar, wet, zm, zgw)
  return wa, zgw, Tr, Es, uex
end

# Temperature Constrains for plant growing
# INPUT:
# - Topt  : Optimum temperature for plant growing
# - Ta    : Air temperature
function Temp_stress(Topt::T, Ta::T) where {T<:Real}
  return exp(-((Ta - Topt) / Topt)^2)
end

export sw_balance
