module SPAC

using Parameters
using LabelledArrays
using UnPack
using Statistics: mean

export Soil
export get_pftpar, get_soilpar
export cal_PET, PT_partition!, interception, snp_balance, runoff_up
export Δz, z₊ₕ
export atm

# Set the soil depth for three soil layers
const Δz = [50.0, 1450.0, 3500.0]  # mm
const z₊ₕ = cumsum(Δz)

const atm = 101.325  # atmospheric pressure (kPa)

include("Param/Soil.jl")
include("SpacOutput.jl")
# include("get_pftpar.jl")
# include("get_soilpar.jl")
# include("Parameter.jl")
# include("ultilize.jl")
include("PET.jl")
include("Evapotranspiration.jl")

include("interception.jl")
include("snp_balance.jl")
include("runoff_up.jl")

# include("GroundWater/GroundWater.jl")
include("Soil_Kun/sw_balance.jl")
include("Soil_CoLM/CoLM.jl")

include("SiTHv2.jl")

end # module SiTH
