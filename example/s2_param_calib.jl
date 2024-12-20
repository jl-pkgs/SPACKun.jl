using Test, SPAC
using DataFrames, RTableTools
import HydroTools: GOF, sceua
using Plots


function init_param(theta::Nothing=nothing; soiltype=2, lc=11)
  soil = Soil{Float64}(; zwt=0.0, snowpack=0.0, soiltype, lc)
  soil.θ .= soil.param.θ_sat
  soil
end

function init_param(theta::Vector{Float64}; soiltype=2, lc=11)
  soil = init_param(nothing; soiltype, lc)
  soil.param.rs .= theta[1:3]
  soil
end

# 气象驱动为全局变量
function ModelSim(theta; method_PET="Monteith65")
  soil = init_param(theta)
  ET, Tr, Es, Ei, Esb, RF, GW, SM =
    SiTHv2_site(soil, Rn, Tavg, Tas, Prcp, Pa, Gi, LAI, s_VODi, topt, VPD, U2;
      spin=false, method_PET)
  ET
end

function goal(theta; method_PET="Monteith65")
  ET = ModelSim(theta; method_PET)
  gof = GOF(ET_obs, ET)
  -gof.NSE
end


begin
  # 装弹
  dir_root = "$(@__DIR__)/.."
  d = fread("$dir_root/data/dat_栾城_ERA5L_2010.csv")
  (; Rn, Pa, Prcp, Tavg, LAI, VOD, VPD, U2, ET_obs, GPP_obs) = d
  Tas = deepcopy(Tavg) # Effective accumulated temperature
  Tas[Tas.<0] .= 0 # Remove values less than 0
  Tas = cumsum(Tas)

  topt = 24.0
  Gi = 0.4 .* Rn .* exp.(-0.5 .* LAI) # G_soil
  s_VODi = (VOD ./ maximum(VOD)) .^ 0.5 # VOD-stress
end



ET_Kun = ModelSim(nothing; method_PET="PT72")
GOF(ET_obs, ET_Kun)

begin
  rs0 = [70.0, 70.0, 70.0]
  lower = [10, 40, 40.0]
  upper = [1000, 1000, 1000.0]
  # theta = [70.0, 140.0, 140.0]
  theta, feval, exitflag = sceua(goal, rs0, lower, upper)
  @show theta

  # 固定参数，与优化参数，模型表现差别不大
  theta = [50, 100, 300.];
  ET = ModelSim(theta; method_PET="Monteith65")
  @show GOF(ET_obs, ET)
end

# NSE = 0.36559951360621057, R2 = 0.5191778406212801, KGE = 0.5759085145817919
# 0.36559951360621057, 0.6730262936469504
begin
  gr(framestyle = :box)
  plot(ET_obs, label="ET_obs", color=:black)
  plot!(ET_Kun, label="ET_kun")
  plot!(ET, label="ET_kong")
end
