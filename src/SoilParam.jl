@with_kw mutable struct SoilParam{FT}
  N::Int = 10
  
  ## 土壤水力
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
  θ_fc::Vector{FT} = fill(0.3, N)      # field capacity, [m3 m-3]
  θ_c::Vector{FT} = fill(0.2, N)       # [not used], critical water content, [m3 m-3]

  ## 植被
  β::Vector{FT} = fill(1.0, N)         # 土壤水分限制因子, [0~1], 1充分供水。
  D50::Vector{FT} = fill(0.1, N)       # [cm]
  D95::Vector{FT} = fill(0.2, N)       # [cm]
  Hc::Vector{FT} = fill(0.1, N)        # [m]
end

function SoilParam(N::Int; FT::Type=Float64, soiltype::Int=2, lc::Int=11, kw...)
  param = SoilParam{FT}(; N, kw...)
  SoilParam!(param; soiltype, lc)
end

function SoilParam!(param::SoilParam; soiltype=2, lc=11, same_layer=true)
  param.same_layer = same_layer
  soilpar = get_soilpar(soiltype)
  pftpar = get_pftpar(lc)

  # 土壤参数
  param.θ_sat .= soilpar.θ_sat
  param.θ_wp .= soilpar.θ_wp
  param.Ksat .= soilpar.Ksat
  param.ψ_sat .= soilpar.ψ_sat
  param.b .= soilpar.b
  param.θ_fc .= soilpar.θ_fc

  ## 后期启用的土壤参数
  # param.α .= soilpar.α
  # param.n .= soilpar.n
  # param.m .= soilpar.m

  # 植被参数
  param.β .= pftpar.β
  param.D50 .= pftpar.D50
  param.D95 .= pftpar.D95
  param.Hc .= pftpar.Hc
  return param
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
