module SPAC

using Parameters
using LabelledArrays
using UnPack

export Soil
export get_pftpar, get_soilpar
export potentialET, pTr_partition!, interception, snp_balance, runoff_up
export Δz, z₊ₕ

# Set the soil depth for three soil layers
const Δz = [50.0, 1450.0, 3500.0]  # mm
const z₊ₕ = cumsum(Δz)

const atm = 101.325  # atmospheric pressure (kPa)

include("DataType.jl")
include("get_pftpar.jl")
include("get_soilpar.jl")

include("pTr_partition.jl")
# include("ultilize.jl")
include("potentialET.jl")
include("Evapotranspiration.jl")
include("interception.jl")
include("snp_balance.jl")
include("runoff_up.jl")

include("GroundWater/GroundWater.jl")

include("Soil_Kun/sw_balance.jl")
include("SiTHv2.jl")


end # module SiTH
