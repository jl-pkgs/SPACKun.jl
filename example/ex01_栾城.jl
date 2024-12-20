# begin
#   k = 1
#   x, y = st[k, [:lon, :lat]]
#   i, j = findnear(x, y, lon, lat)
#   soil_type = Soil[i, j]
#   soilpar = get_soilpar(soil_type)

#   topt = Float64(Topt[i, j])
#   soil = init_param(;soil_type, lc=22)
# end

using SPAC, Ipaper, Dates
using RTableTools, DataFrames, NaNStatistics
# using Ipaper, Ipaper.sf, ArchGDAL
# using GLMakie, MakieLayers

function init_param(soiltype=2, lc=11)
  soil = Soil{Float64}(; zwt=0.0, snowpack=0.0, soiltype, lc)
  soil.θ .= soil.param.θ_sat
  soil
end

# Load necessary data
df = fread("data/dat_栾城_ERA5L_1982-2019.csv")
dates = df.date
soil = init_param()

inds = findall(year.(dates) .== 2000)
d = df[inds, :]
d.LAI = d.LAI |> drop_missing
d.VOD = d.VOD |> drop_missing
(; Rn, Pa, Prcp, Tavg, LAI, VOD) = d

Tas = deepcopy(Tavg) # Effective accumulated temperature
Tas[Tas.<0] .= 0 # Remove values less than 0
Tas = cumsum(Tas)

Gi = 0.4 .* Rn .* exp.(-0.5 .* LAI) # G_soil
s_VODi = (VOD ./ nanmaximum(VOD)) .^ 0.5 # VOD-stress
topt = 24.0

ET, Tr, Es, Ei, Esb, SM, RF, GW = 
  SiTHv2_site(soil, Rn, Tavg, Tas, Prcp, Pa, Gi, LAI, s_VODi, topt; spin=false)

# df_out = DataFrame(; ET, Tr, Es, Ei, Esb, SM1, SM2, SM3, RF, GW)
# fwrite(df_out, "data/Output_栾城_2010.csv")
begin
  using Plots
  gr(framestyle = :box)
  plot(ET)
end
