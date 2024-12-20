export Soil, update_soil


@with_kw mutable struct Soil{FT}
  Δz::Vector{FT} = [50.0, 1450.0, 3500.0]  # mm
  z₊ₕ::Vector{FT} = cumsum(Δz)
  N::Int = length(Δz)

  jwt::Int = 0                       # index of groundwater table, 追踪地下水，所在位置的上一层
  θ::Vector{FT} = ones(N) .* 0.2
  θ_prev::Vector{FT} = fill(0.2, N)  # previous soil water content
  θ_unsat::Vector{FT} = fill(0.3, N)

  Q::Vector{FT} = fill(0.0, N)       # drainage, discharge, [mm d-1]

  Ec_pot::Vector{FT} = fill(0.0, N)  # potential transpiration
  Ec_sm::Vector{FT} = fill(0.0, N)   # actual transpiration from unsaturated zone
  Ec_gw::Vector{FT} = fill(0.0, N)   # actual transpiration from groundwater

  Es_pot::Vector{FT} = fill(0.0, N)  # potential soil evaporation
  Es_sm::Vector{FT} = fill(0.0, N)   # actual transpiration from unsaturated zone
  Es_gw::Vector{FT} = fill(0.0, N)   # actual transpiration from groundwater

  sink::Vector{FT} = fill(0.0, N)    # sink, 来自土壤的Ec + Es

  fsm_Es::Vector{FT} = fill(1.0, N)  # SM constraint for soil evaporation
  fsm_Ec::Vector{FT} = fill(1.0, N)  # SM constraint for transpiration

  zwt::FT = 0.0                      # groundwater depths [mm]
  Sy::Vector{FT} = fill(0.02, N)     # specific yield, [m3 m-3]
  Syₙ::FT = 0.2                      # specific yield of `N+1` layer, [m3 m-3] 
  wa::FT = 0.0                       # water amount in aquifer, [mm]，潜水含水层
  uex::FT = 0.0                      # 超出地表的水量, [mm], [kg m-2] 以地表径流的形式排放
  recharge::FT = 0.0                 # recharge rate, [mm/s]
  drainage::FT = 0.0                 # drainage rate, [mm/s]
  snowpack::FT = 0.0                 # snowpack depth [mm]

  ## Parameters
  Dmin::Vector{FT} = [0.048, 0.012, 0.012]  # drainage Parameters
  Dmax::Vector{FT} = [4.8, 1.2, 1.2]
  
  soiltype::Int = 2
  LC::Int = 11
  soilpar::LArray{Float64,1} = get_soilpar(soiltype)
  pftpar::LArray{Float64,1} = get_pftpar(LC)
end


"""
状态变量需要连续，传递到下一年中；模型对初始状态soil敏感。
- θ  : 采用warming-up的方式获取，warming-up period可取3年
- zwt: 采用spin-up的方式获取，spin-up period可取100年

# Argument Specification
- θ: soil water content in three layers, [m^3 m^-3]
- zwt: groundwater depth, [mm]
- snowpack: snowpack depth
"""
function update_soil!(soil, θ, zwt, snowpack)
  θ[θ.<0] .= 0.01  # set the minimum value for soil moisture
  soil.θ = θ
  soil.zwt = zwt
  soil.snowpack = snowpack
  return soil
end
