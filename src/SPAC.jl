module SPAC

using Parameters
using LabelledArrays

export State
export get_pftpar, get_soilpar
export potentialET, pTr_partition, interception, snp_balance, runoff_up
export ZM

# Set the soil depth for three soil layers
const ZM = [50, 1450, 3500]  # mm

include("DataType.jl")
include("Param/get_pftpar.jl")
include("Param/get_soilpar.jl")

include("pTr_partition.jl")
include("ultilize.jl")
include("potentialET.jl")
include("cal_SITH_site.jl")
include("interception.jl")
include("snp_balance.jl")
include("runoff_up.jl")

include("Soil/Soil.jl")



"""

# Model input
- Rn     : Net Radiation, W/m-2
- Ta     : Near surface air temperature, C
- Topt   : Optimal growth temperature for plant, C
- Pe     : Precipitation, mm day-1
- Pa     : Near surface air pressure, kPa
- s_VOD  : Constrains from VOD,[0,1]
- G      : Soil heat flux, W/m-2
- LAI    : Leaf area index
- soilpar: Parameters related to Soil types
- pftpar : Parameters related to PFTs
- wa     : Soil moisture (last step)
- zgw    : Groundwater table depth, mm
- snp    : Snow package (old), mm day-1

# Model output
- Et     : Total Evapotranspiration, mm day-1
- Tr     : Plant Transpiration, mm day-1
- Es     : Soil Evaporation, mm day-1
- Ei     : Intercepted Evaporation, mm day-1
- Esb    : Snow sublimation, mm day-1
- wa     : Soil moisture (three layers)
- srf    : Surface runoff, mm day-1
- zgw    : Groundwater table depth, mm
- snp    : Snow package (new), mm day-1

"""
function SiTHv2(Rn, Ta, Tas, Topt, Pe, Pa, s_VOD, G, LAI, soilpar, pftpar, state::State)
  (; sm, zg, snowpack) = state
  # Potential Evaporation allocated to canopy and soil surface
  pEc, pEs = potentialET(Rn, G, LAI, Ta, Pa)

  # Interception evaporation
  Ei, fwet, _ = interception(pEc, LAI, Pe, pftpar)

  # Snow sublimation, snow melt
  new_Pe = max(Pe - Ei, 0)
  snowpack, Esb, _, Pnet = snp_balance(new_Pe, Ta, Tas, snowpack, pEs)

  srf, IWS, Vmax = runoff_up(Pnet, zg, ZM, sm, soilpar)

  # Variables associated with soil water balance
  new_pEs = max(pEs - Esb, 0)
  sm, zg, Tr, Es, uex = sw_balance(IWS, pEc, new_pEs, Ta, Topt, s_VOD, sm, soilpar, pftpar, fwet, ZM, zg)

  # Total Evapotranspiration
  Et = Tr + Es + Ei + Esb
  srf += uex
  update_state!(state, sm, zg, snowpack)
  return Et, Tr, Es, Ei, Esb, srf, Pnet, IWS, Vmax
end


end # module SiTH
