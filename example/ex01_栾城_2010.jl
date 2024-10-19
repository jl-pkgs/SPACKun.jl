using SITH, Ipaper, Ipaper.sf, ArchGDAL
using RTableTools, DataFrames, NaNStatistics
# using GLMakie, MakieLayers

cellsize = 0.1
b = bbox(-180, -90, 180, 90)
lon, lat = bbox2dims(b; cellsize)
Soil = read_gdal("data/param_Soil_G010.tif", 1) 
Topt = read_gdal("data/param_Topt_G010.tif", 1) |> x -> Float32.(x)
st = fread("./data/siteInfo_CRO_6sp.csv")
# serialize("data/param_GO10", (; lon, lat, Soil, Topt))

# x, y = (110.23, 20.3)
# imagesc(lon, lat, Topt)
begin
  k = 1
  x, y = st[k, [:lon, :lat]]
  i, j = findnear(x, y, lon, lat)
  soil_type = Soil[i, j]

  topt = Float64(Topt[i, j])
  PFTi = 22
end

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

# Load necessary data
begin
  df = fread("data/dat_栾城_ERA5L_1982-2019.csv")
  dates = df.date

  inds = findall(year.(dates) .== 2010)
  d = df[inds, :]
  d.LAI = d.LAI |> drop_missing
  d.VOD = d.VOD |> drop_missing
  fwrite(d, "data/dat_栾城_ERA5L_2010.csv")
end


soilpar, pftpar, state = init_param()
topt = 24.0

d = fread("data/dat_栾城_ERA5L_2010.csv")
(; Rn, Pa, Prcp, Tavg, LAI, VOD) = d

Tas = Tavg # Effective accumulated temperature
Tas[Tas.<0] .= 0 # Remove values less than 0
Tas = cumsum(Tas)

Gi = 0.4 .* Rn .* exp.(-0.5 .* LAI) # G_soil
s_VODi = (VOD ./ nanmaximum(VOD)) .^ 0.5 # VOD-stress

ET, Tr, Es, Ei, Esb, SM, RF, GW = 
  cal_SiTHv2_site(Rn, Tavg, Tas, Prcp, Pa, Gi, LAI, s_VODi, topt, soilpar, pftpar, state, false)

SM1 = SM[:, 1]
SM2 = SM[:, 2]
SM3 = SM[:, 3]
df_out = DataFrame(; ET, Tr, Es, Ei, Esb, SM1, SM2, SM3, RF, GW)
fwrite(df_out, "data/Output_栾城_2010.csv")

begin
  using Plots
  gr(framestyle = :box)
  plot(ET)
end
