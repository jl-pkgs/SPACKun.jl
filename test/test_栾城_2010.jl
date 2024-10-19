using Test, SPAC, RTableTools

function init_param(soil_type=2, PFTi = 22)
  soilpar = get_soilpar(soil_type)
  pftpar = get_pftpar(PFTi)

  θ_sat = soilpar.θ_sat
  wa = ones(3) * θ_sat
  zg = 0.0
  snowpack = 0.0
  state = State(; wa, zg, snowpack)
  soilpar, pftpar, state
end


dir_root = "$(@__DIR__)/.."
d = fread("$dir_root/data/dat_栾城_ERA5L_2010.csv")

# @testset "SiTHv2_site" 
begin
  soilpar, pftpar, state = init_param()
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
  df_out = DataFrame(; ET, Tr, Es, Ei, Esb, RF, GW, SM1, SM2, SM3)
  fwrite(df_out, "$dir_root/data/OUTPUT_栾城_2010.csv")
end

# begin
#   using Plots
#   gr(framestyle = :box)
#   plot(ET)
# end
