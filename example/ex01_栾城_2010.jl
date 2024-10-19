using SPAC
using RTableTools, DataFrames, NaNStatistics

function init_param(soil_type=2, PFTi = 22)
  soilpar = get_soilpar(soil_type)
  pftpar = get_pftpar(PFTi)

  θ_sat = soilpar.θ_sat
  sm = ones(3) * θ_sat
  zg = 0.0
  snowpack = 0.0
  state = State(; sm, zg, snowpack)
  soilpar, pftpar, state
end


soilpar, pftpar, state = init_param()
topt = 24.0

d = fread("data/dat_栾城_ERA5L_2010.csv")
(; Rn, Pa, Prcp, Tavg, LAI, VOD) = d

Tas = deepcopy(Tavg) # Effective accumulated temperature
Tas[Tas.<0] .= 0 # Remove values less than 0
Tas = cumsum(Tas)

Gi = 0.4 .* Rn .* exp.(-0.5 .* LAI) # G_soil
s_VODi = (VOD ./ nanmaximum(VOD)) .^ 0.5 # VOD-stress

@time ET, Tr, Es, Ei, Esb, SM, RF, GW = 
  cal_SiTHv2_site(Rn, Tavg, Tas, Prcp, Pa, Gi, LAI, s_VODi, topt, soilpar, pftpar, state, false)

SM1 = SM[:, 1]
SM2 = SM[:, 2]
SM3 = SM[:, 3]
df_out = DataFrame(; ET, Tr, Es, Ei, Esb, RF, GW, SM1, SM2, SM3)
fwrite(df_out, "data/OUTPUT_栾城_2010.csv")


df_mat = fread("./data/OUTPUT_栾城_2010_MATLAB.csv")
df_jl = fread("./data/OUTPUT_栾城_2010.csv")
df_mat .- df_jl

# begin
#   using Plots
#   gr(framestyle = :box)
#   plot(ET)
# end
