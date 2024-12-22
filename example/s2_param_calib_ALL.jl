using Test, SPAC
using DataFrames, RTableTools, Dates
using Plots
import HydroTools: GOF, sceua
import Ipaper: drop_missing
GOF_df(yobs, ysim; kw...) = DataFrame(; kw..., GOF(yobs, ysim)...)
include("main_pkgs.jl")

df = fread("data/Forcing_SPAC_CRO_3sp_ERA5L_FluxMet.csv")
replace_miss!(df)
sites = ["固城", "禹城", "栾城"]

function build_forcing(d)
  doy = Dates.dayofyear.(d.date)
  (; Rn, Pa, Prcp, Tavg, LAI, LAI_lu, VOD, VPD, U2) = d # ET_obs, GPP_obs
  LAI .= LAI_lu

  Tas = deepcopy(Tavg) # Effective accumulated temperature
  Tas[Tas.<0] .= 0 # Remove values less than 0
  Tas = cumsum(Tas)

  topt = 24.0
  Gi = 0.4 .* Rn .* exp.(-0.5 .* LAI) # G_soil
  s_VODi = (VOD ./ maximum(VOD)) .^ 0.5 # VOD-stress
  (; Rn, Tavg, Tas, Prcp, Pa, Gi, LAI, s_VODi, topt, VPD, U2, doy)
end

function ModelSim(forcing, par, theta=nothing; method_PET="Monteith65")
  soil = init_param(theta; par)
  SiTHv2_site(soil, forcing...; spin=false, method_PET)[1] # return ET
end

function goal(forcing, par, theta, ET_obs; method_PET="Monteith65")
  ET = ModelSim(forcing, par, theta; method_PET)
  r = GOF(ET_obs, ET)
  -r.NSE
  # -r.KGE
end


function ModelCalib(site::String)
  printstyled("[Running]: $site ... \n", color=:green, bold=true)
  d = df[df.site.==site, :]
  forcing = build_forcing(d)

  parNames = [:_rs] # :HC
  par = select_param(parNames)

  ET_Kun = ModelSim(forcing, par, nothing; method_PET="PT72")
  r_kun = GOF_df(d.ET_obs, ET_Kun; site, model="Kun")

  theta0 = unlist_field(par, :x0)
  lower = unlist_field(par, :lower)
  upper = unlist_field(par, :upper)
  f(theta) = goal(forcing, par, theta, d.ET_obs; method_PET="Monteith65")
  theta, feval, exitflag = sceua(f, theta0, lower, upper)
  @show theta

  theta = [70.0, 300.0, 600.0] # 固定参数也能取得不错的表现
  ET = ModelSim(forcing, par, theta; method_PET="Monteith65")
  r_kong = GOF_df(d.ET_obs, ET; site, model="Kong")
  vcat(r_kun, r_kong)
end

@time RES = map(ModelCalib, sites);
# @profview RES = map(ModelCalib, sites);
R = vcat(RES...)

# - 固城: 0.42
# - 禹城: 0.68
# - 栾城: 0.61

# ModelCalib(sites[1])
# ModelCalib(sites[2])
# ModelCalib(sites[3])

# p1 = plot(ET_obs, label="ET_obs")
# plot!(ET, label="ET")
# p2 = plot(ET_obs, label="ET_obs")
# plot!(ET_Kun, label="ET_Kun")

# plot(p1, p2, 
#   plot(LAI),
#   plot(VOD),
#   size=(1000, 700)) # VOD数据质量不佳
