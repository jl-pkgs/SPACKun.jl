using Test, SPAC, RTableTools
using DataFrames


function absmax(d::AbstractDataFrame; fun=maximum)
  absmax(x) = fun(abs.(x))
  values = map(j -> absmax(d[:, j]), 1:ncol(d))
  keys = tuple(names(d)...)
  NamedTuple{Symbol.(keys)}(values)
end

function init_param(soil_type=2, PFTi=11)
  soilpar = get_soilpar(soil_type)
  pftpar = get_pftpar(PFTi)

  θ_sat = soilpar.θ_sat
  θ = ones(3) * θ_sat
  zwt = 0.0
  snowpack = 0.0
  soil = Soil{Float64}(; θ, zwt, snowpack)
  soilpar, pftpar, soil
end

dir_root = "$(@__DIR__)/.."
d = fread("$dir_root/data/dat_栾城_ERA5L_2010.csv")


function test_LuanCheng(; zwt=0.0, method="Kun")
  soilpar, pftpar, soil = init_param()
  soil.zwt = zwt
  topt = 24.0

  (; Rn, Pa, Prcp, Tavg, LAI, VOD) = d

  Tas = deepcopy(Tavg) # Effective accumulated temperature
  Tas[Tas.<0] .= 0 # Remove values less than 0
  Tas = cumsum(Tas)

  Gi = 0.4 .* Rn .* exp.(-0.5 .* LAI) # G_soil
  s_VODi = (VOD ./ maximum(VOD)) .^ 0.5 # VOD-stress

  @time ET, Tr, Es, Ei, Esb, RF, GW, SM =
    SiTHv2_site(Rn, Tavg, Tas, Prcp, Pa, Gi, LAI, s_VODi, topt, soilpar, pftpar, soil, false; method)

  SM1 = SM[:, 1]
  SM2 = SM[:, 2]
  SM3 = SM[:, 3]
  df_jl = DataFrame(; ET, Tr, Es, Ei, Esb, RF, GW, SM1, SM2, SM3)

  fwrite(df_jl, "$dir_root/data/OUTPUT/OUTPUT_栾城_2010_zgw=$(Int(zwt)).csv")
end



# @run test_LuanCheng(; zwt=6000.0)
@testset "SPAC.jl" begin
  zgws = [0, 25, 1000, 2000, 6000]
  for zwt = zgws[1:5]
    test_LuanCheng(; zwt)

    f_jl = "$dir_root/data/OUTPUT/OUTPUT_栾城_2010_zgw=$(Int(zwt)).csv"
    f_mat = "$dir_root/data/OUTPUT/OUTPUT_栾城_2010_MATLAB_zgw=$(Int(zwt)).csv"

    df_jl = fread(f_jl)
    df_mat = fread(f_mat)

    diff = df_mat .- df_jl
    absmax(diff)
    # maximum(absmax(diff))
    @test maximum(absmax(diff)) <= 1e-9
  end
end

# begin
#   using Plots
#   function plot_var(var=:ET)
#     plot(df_mat[:, var], label="MATLAB", title=string(var))
#     plot!(df_jl[:, var], label="Julia")
#   end

#   gr(; framestyle=:box)
#   plot(
#     plot_var(:ET),
#     plot_var(:Tr),
#     plot_var(:Es),
#     plot_var(:Ei),
#     plot_var(:Esb),
#     plot_var(:RF),
#     plot_var(:GW),
#     plot_var(:SM1),
#     plot_var(:SM2),
#     plot_var(:SM3),
#     size=(1400, 800)
#   )
# end
