using Test, SPAC, RTableTools
using DataFrames


function absmax(d::AbstractDataFrame; fun=maximum)
  absmax(x) = fun(abs.(x))
  values = map(j -> absmax(d[:, j]), 1:ncol(d))
  keys = tuple(names(d)...)
  NamedTuple{Symbol.(keys)}(values)
end

function init_param(soil_type=2, PFTi=22)
  soilpar = get_soilpar(soil_type)
  pftpar = get_pftpar(PFTi)

  θ_sat = soilpar.θ_sat
  wa = ones(3) * θ_sat
  zgw = 0.0
  snowpack = 0.0
  state = State(; wa, zgw, snowpack)
  soilpar, pftpar, state
end

dir_root = "$(@__DIR__)/.."
d = fread("$dir_root/data/dat_栾城_ERA5L_2010.csv")

function test_LuanCheng(; zgw=0.0)
  soilpar, pftpar, state = init_param()
  state.zgw = zgw
  topt = 24.0

  (; Rn, Pa, Prcp, Tavg, LAI, VOD) = d

  Tas = deepcopy(Tavg) # Effective accumulated temperature
  Tas[Tas.<0] .= 0 # Remove values less than 0
  Tas = cumsum(Tas)

  Gi = 0.4 .* Rn .* exp.(-0.5 .* LAI) # G_soil
  s_VODi = (VOD ./ maximum(VOD)) .^ 0.5 # VOD-stress

  @time ET, Tr, Es, Ei, Esb, SM, RF, GW =
    SiTHv2_site(Rn, Tavg, Tas, Prcp, Pa, Gi, LAI, s_VODi, topt, soilpar, pftpar, state, false)

  SM1 = SM[:, 1]
  SM2 = SM[:, 2]
  SM3 = SM[:, 3]
  df_jl = DataFrame(; ET, Tr, Es, Ei, Esb, RF, GW, SM1, SM2, SM3)

  fwrite(df_jl, "$dir_root/data/OUTPUT_栾城_2010_zgw=$(Int(zgw)).csv")
end

# const ZM = [50, 1450, 3500]  # mm
test_LuanCheng(; zgw=0.0)
test_LuanCheng(; zgw=25.0)
test_LuanCheng(; zgw=1000.0)
test_LuanCheng(; zgw=2000.0)
test_LuanCheng(; zgw=6000.0)

@testset "SiTHv2_site" begin
  zgw = 0.0
  f_jl = "$dir_root/data/OUTPUT_栾城_2010_zgw=$(Int(zgw)).csv"
  df_jl = fread(f_jl)
  df_mat = fread("$dir_root/data/OUTPUT_栾城_2010_MATLAB.csv")

  diff = df_mat .- df_jl
  @test maximum(absmax(diff)) <= 1e-9
end

# begin
#   using Plots
#   gr(framestyle = :box)
#   plot(ET)
# end
