using Test, SPAC
using DataFrames, RTableTools, Dates
using Plots
import HydroTools: GOF, sceua
GOF_df(yobs, ysim) = DataFrame([GOF(yobs, ysim)])

# > 固城，禹城KGE只有0.4
sites = ["固城", "禹城", "栾城"]
site = sites[2]

begin
  # 装弹
  d = fread("data/calib/dat_ERA5L_$site.csv")
  (; Rn_obs, Prcp_obs, Rn, Pa, Prcp, Tavg, LAI, VOD, VPD, U2, ET_obs, GPP_obs) = d
  Rn_obs = drop_missing(Rn_obs, 0.0)
  # Rn_obs = drop_missing(Rn_obs)
  doy = Dates.dayofyear.(d.date)
  ET_obs = drop_missing(ET_obs, NaN)
  GPP_obs = drop_missing(GPP_obs, NaN)

  Tas = deepcopy(Tavg) # Effective accumulated temperature
  Tas[Tas.<0] .= 0 # Remove values less than 0
  Tas = cumsum(Tas)

  topt = 24.0
  Gi = 0.4 .* Rn_obs .* exp.(-0.5 .* LAI) # G_soil
  s_VODi = (VOD ./ maximum(VOD)) .^ 0.5 # VOD-stress
end

# 全局变量: 气象驱动, par
function ModelSim(theta; method_PET="Monteith65")
  soil = init_param(theta; par)
  # ET, Tr, Es, Ei, Esb, RF, GW, SM =
  SiTHv2_site(soil, Rn, Tavg, Tas, Prcp, Pa, Gi, LAI, s_VODi, topt, VPD, U2, doy;
    spin=false, method_PET)[1]
end

function goal(theta; method_PET="Monteith65")
  ET = ModelSim(theta; method_PET)
  r = GOF(ET_obs, ET)
  -r.NSE
  # -r.KGE
end


begin
  parNames = [:rs] # :Hc
  par = select_param(parNames)
  
  ET_Kun = ModelSim(nothing; method_PET="PT72")
  GOF_df(ET_obs, ET_Kun)
  #  Row │ NSE       R2        KGE       R         RMSE     MAE       bias       bias_perc  n_valid 
  #    1 │ 0.534438  0.701603  0.547878  0.837618  1.11244  0.770788  -0.647735   -33.3154     4018

  theta0 = unlist_field(par, :x0)
  lower = unlist_field(par, :lower)
  upper = unlist_field(par, :upper)

  theta, feval, exitflag = sceua(goal, theta0, lower, upper)
  @show theta

  theta0 = [70.0, 300.0, 600.0] # rs
  theta0 = theta
  ET = ModelSim(theta0; method_PET="Monteith65")
  GOF_df(ET_obs, ET) |> print # 仍然存在系统偏差
  # 固定参数，与优化参数，模型表现差别不大
end


## 可能是LAI不对
p1 = plot(ET_obs, label="ET_obs")
plot!(ET, label="ET")

p2 = plot(ET_obs, label="ET_obs")
plot!(ET_Kun, label="ET_Kun")

plot(p1, p2, 
  plot(LAI),
  plot(VOD),
  size=(1000, 700)) # VOD数据质量不佳
