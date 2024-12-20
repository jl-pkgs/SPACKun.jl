using Test, SPAC
using DataFrames, RTableTools
import HydroTools: GOF
using Plots


function absmax(d::AbstractDataFrame; fun=maximum)
  absmax(x) = fun(abs.(x))
  values = map(j -> absmax(d[:, j]), 1:ncol(d))
  keys = tuple(names(d)...)
  NamedTuple{Symbol.(keys)}(values)
end

function init_param(soiltype=2, lc=11)
  soil = Soil{Float64}(; zwt=0.0, snowpack=0.0, soiltype, lc)
  soil.θ .= soil.param.θ_sat
  soil
end

dir_root = "$(@__DIR__)/.."
d = fread("$dir_root/data/dat_栾城_ERA5L_2010.csv")


begin
  method_sw = "Kun"
  soil = init_param()
  soil.zwt = 0.0
  topt = 24.0

  ## 之后尝试优化模型参数
  (; Rn, Pa, Prcp, Tavg, LAI, VOD, ET_obs, GPP_obs) = d

  Tas = deepcopy(Tavg) # Effective accumulated temperature
  Tas[Tas.<0] .= 0 # Remove values less than 0
  Tas = cumsum(Tas)

  Gi = 0.4 .* Rn .* exp.(-0.5 .* LAI) # G_soil
  s_VODi = (VOD ./ maximum(VOD)) .^ 0.5 # VOD-stress

  @time ET, Tr, Es, Ei, Esb, RF, GW, SM =
    SiTHv2_site(soil, Rn, Tavg, Tas, Prcp, Pa, Gi, LAI, s_VODi, topt; spin=false, method_sw)

  SM1 = SM[:, 1]
  SM2 = SM[:, 2]
  SM3 = SM[:, 3]
  df_jl = DataFrame(; ET, Tr, Es, Ei, Esb, RF, GW, SM1, SM2, SM3)
  
  GOF(ET_obs, ET)
  # fwrite(df_jl, "$dir_root/data/OUTPUT/OUTPUT_栾城_2010_zgw=$(Int(zwt)).csv")
end
# NSE = 0.36559951360621057, R2 = 0.5191778406212801, KGE = 0.5759085145817919

## 尝试根据ET优化模型参数
begin
  gr(framestyle = :box)
  plot(ET, label="ET_sim")
  plot!(ET_obs, label="ET_obs")
end
