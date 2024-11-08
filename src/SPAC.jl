module SPAC

using Parameters
using LabelledArrays

export State
export get_pftpar, get_soilpar
export potentialET, pTr_partition, interception, snp_balance, runoff_up
export ZM

# Set the soil depth for three soil layers
const ZM = [50.0, 1450.0, 3500.0]  # mm
const atm = 101.325  # atmospheric pressure (kPa)

include("DataType.jl")
include("Param/get_pftpar.jl")
include("Param/get_soilpar.jl")

include("pTr_partition.jl")
# include("ultilize.jl")
include("potentialET.jl")
include("interception.jl")
include("snp_balance.jl")
include("runoff_up.jl")
include("Soil/Soil.jl")
include("SiTHv2.jl")


end # module SiTH
