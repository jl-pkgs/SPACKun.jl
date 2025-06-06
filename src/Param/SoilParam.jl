export SoilParam
export period_wheat, period_corn


period_wheat(doy::Int) = 32 <= doy <= 166
period_corn(doy::Int) = 196 <= doy <= 288

# [wheat, corn, non-gw]
function iGrowthPeriod(doy::Int)
  i = 3
  period_wheat(doy) && (i = 1)
  period_corn(doy) && (i = 2)
  return i
end


# 为多层土壤参数做好铺垫
@with_kw mutable struct SoilParam{FT}
  N::Int = 10

  ## 土壤水力
  method::String = "Campbell"     # "van_Genuchten" or "Campbell"
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
  β::FT = 1.0                          # 冠层截留系数，[-]
  D50::FT = 0.1                        # [cm]
  D95::FT = 0.2                        # [cm]

  α_soil::FT = 1.26
  hc::FT = 0.1                         # 冠层高度，[m]
  rs::FT = 70.0                        # 气孔导度阻力，[s m-1]
  kA::FT = 0.6                         # 消光系数，[m-1]
  Kc::FT = 1.0                         # ET0转换系数
  ## 植被类型变动
  _Kc::Vector{FT} = [1.0, 1.0, 1.0]    # [ungrowing, wheat, cron]
  _hc::Vector{FT} = [0.12, 0.5, 1.0]
  _rs::Vector{FT} = [70.0, 70.0, 70.0]
  _kA::Vector{FT} = [0.6, 0.6, 0.6]
end

## 引入VegParam更新机制
"update_VegParam!(soil.param, doy)"
function update_VegParam!(param::SoilParam, doy::Int)
  doy <= 0 && return
  i = iGrowthPeriod(doy)
  param.hc = param._hc[i]
  param.rs = param._rs[i]
  param.kA = param._kA[i]
  param.Kc = param._Kc[i]
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
  (; β, D50, D95, hc) = pftpar
  @pack! param = β, D50, D95, hc
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
  print_var(io, param, :Ksat) # 1e-3
  print_var(io, param, :α)
  print_var(io, param, :n)
  use_m && print_var(io, param, :m; used=use_m)

  print_selected(io, "Campbell (4p$subfix)", method)
  # printstyled(io, " - θ_sat, Ksat \n", color=:blue)
  print_var(io, param, :θ_sat)
  print_var(io, param, :Ksat) # 1e-3

  print_var(io, param, :ψ_sat)
  print_var(io, param, :b)

  println(io, "-----------------------------")
  print_var(io, param, :β)
  print_var(io, param, :D50)
  print_var(io, param, :D95)
  print_var(io, param, :hc)

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
