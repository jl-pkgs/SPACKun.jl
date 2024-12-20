using Test, SPAC, RTableTools
using DataFrames


function absmax(d::AbstractDataFrame; fun=maximum)
  absmax(x) = fun(abs.(x))
  values = map(j -> absmax(d[:, j]), 1:ncol(d))
  keys = tuple(names(d)...)
  NamedTuple{Symbol.(keys)}(values)
end

function init_param(soil_type=2, lc=11)
  soilpar = get_soilpar(soil_type)
  pftpar = get_pftpar(lc)

  θ_sat = soilpar.θ_sat
  θ = ones(3) * θ_sat
  zwt = 0.0
  snowpack = 0.0
  soil = Soil(; θ, zwt, snowpack)
  soilpar, pftpar, soil
end


dir_root = "$(@__DIR__)/.."
d = fread("$dir_root/data/dat_栾城_ERA5L_2010.csv")

method = "Kun"
test_LuanCheng(; zwt=0.0, method)
@run test_LuanCheng(; zwt=0.0, method)
# @testset "SPAC.jl" 
dir_root = "$(@__DIR__)/.." |> abspath

begin
  zgws = [0, 25, 1000, 2000, 6000]
  zwt = zgws[1]
  # for zwt = zgws[1:5]
  test_LuanCheng(; zwt, method)

  f_jl = "$dir_root/data/OUTPUT/OUTPUT_栾城_2010_zgw=$(Int(zwt)).csv"
  f_mat = "$dir_root/data/OUTPUT/OUTPUT_栾城_2010_MATLAB_zgw=$(Int(zwt)).csv"

  df_jl = fread(f_jl)
  df_mat = fread(f_mat)

  diff = df_mat .- df_jl
  absmax(diff)
  # maximum(absmax(diff))
  @test maximum(absmax(diff)) <= 1e-9
  # end
end

begin
  using Plots
  function plot_var(var=:ET)
    plot(df_mat[:, var], label="MATLAB", title=string(var))
    plot!(df_jl[:, var], label="Julia")
  end

  gr(; framestyle=:box)
  plot(
    plot_var(:ET),
    plot_var(:Tr),
    plot_var(:Es),
    plot_var(:Ei),
    plot_var(:Esb),
    plot_var(:RF),
    plot_var(:GW),
    plot_var(:SM1),
    plot_var(:SM2),
    plot_var(:SM3),
    size=(1400, 800)
  )
end
