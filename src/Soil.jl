export Soil, SoilParam, update_soil
using Printf


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
  Dmin::Vector{FT} = [0.048, 0.012, 0.012]  # drainage Parameters, ! 后期移除的参数
  Dmax::Vector{FT} = [4.8, 1.2, 1.2]
  
  soiltype::Int = 2
  LC::Int = 11
  soilpar::LArray{Float64,1} = get_soilpar(soiltype)
  pftpar::LArray{Float64,1} = get_pftpar(LC)
end


@with_kw mutable struct SoilParam{FT}
  ## Parameter: 土壤水力
  N::Int = 10
  method::String = "van_Genuchten"     # "van_Genuchten" or "Campbell"
  use_m::Bool = false
  same_layer = true

  # van Genuchten
  θ_sat::Vector{FT} = fill(0.4, N)     # saturated water content, [m3 m-3]
  θ_wp::Vector{FT} = fill(0.1, N)      # residual water content, [m3 m-3]
  Ksat::Vector{FT} = fill(400.0, N)    # saturated hydraulic conductivity, [mm d-1]
  α::Vector{FT} = fill(0.01, N)        # [m-1]
  n::Vector{FT} = fill(2.0, N)         # [-]
  m::Vector{FT} = fill(0.5, N)         # [-]，优化时的可选参数

  # Campbell
  ψ_sat::Vector{FT} = fill(-10.0, N)   # [cm]
  b::Vector{FT} = fill(4.0, N)         # [-]
  
  # 其他
  θ_c::Vector{FT} = fill(0.2, N)       # [not used], critical water content, [m3 m-3]
end


function Base.show(io::IO, param::SoilParam{T}) where {T<:Real}
  (; use_m, same_layer) = param
  printstyled(io, "Parameters: \n", color=:blue, bold=true)
  # println("[use_m = $use_m, same_layer = $same_layer]")

  println(io, "-----------------------------")
  # print_var(io, param, :κ)
  # print_var(io, param, :cv; scale=1e6)
  # println(io, "-----------------------------")

  method = param.method
  subfix = same_layer ? " * 1" : " * N"
  np = use_m ? 6 : 5
  print_selected(io, "van_Genuchten ($(np)p$subfix)", method)
  print_var(io, param, :θ_sat)
  print_var(io, param, :θ_wp)
  print_var(io, param, :Ksat; scale=1e-3)
  print_var(io, param, :α)
  print_var(io, param, :n)
  use_m && print_var(io, param, :m; used=use_m)
  print_selected(io, "Campbell (4p$subfix)", method)
  printstyled(io, " - θ_sat, Ksat \n", color=:blue)

  print_var(io, param, :ψ_sat)
  print_var(io, param, :b)
  return nothing
end


function print_var(io::IO, x, var; scale=nothing, digits=3, color=:blue, used=true)
  value = getfield(x, var)
  name = @sprintf("%-5s", string(var))
  _color = used ? color : :white
  printstyled(io, " - $name: "; color=_color)
  if isnothing(scale)
    println(io, round.(value; digits))
  else
    println(io, "$(round.(value/scale; digits)) * $scale")
  end
end

function print_selected(io::IO, name::String, method::String)
  if name[1:5] == method[1:5]
    printstyled(io, "   [$name]\n", bold=true, color=:green)
  else
    printstyled(io, "   [$name]\n", bold=true)
  end
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
