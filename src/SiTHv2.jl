export SiTHv2_site

"""

# Model input
- Rn     : Net Radiation, W/m-2
- Ta     : Near surface air temperature, C
- Topt   : Optimal growth temperature for plant, C
- Pe     : Precipitation, mm day-1
- Pa     : Near surface air pressure, kPa
- s_VOD  : Constrains from VOD,[0,1]
- G      : Soil heat flux, W/m-2
- LAI    : Leaf area index
- soilpar: Parameters related to Soil types
- pftpar : Parameters related to PFTs
- θ     : Soil moisture (last step)
- zwt    : Groundwater table depth, mm
- snp    : Snow package (old), mm day-1

# Model output
- Et     : Total Evapotranspiration, mm day-1
- Tr     : Plant Transpiration, mm day-1
- Es     : Soil Evaporation, mm day-1
- Ei     : Intercepted Evaporation, mm day-1
- Esb    : Snow sublimation, mm day-1
- θ     : Soil moisture (three layers)
- srf    : Surface runoff, mm day-1
- zwt    : Groundwater table depth, mm
- snp    : Snow package (new), mm day-1

"""
function SiTHv2!(output::SpacOutput{T},
  Rn::T, Ta::T, Tas::T, Topt::T, P::T, Pa::T, s_VOD::T,
  G::T, LAI::T, soilpar, pftpar, soil::Soil; 
  Kc=1.0, method="Kun") where {T<:Real}
  
  funcs = Dict(
    "Kun" => sw_balance,
    "CoLM" => sw_balance_CoLM)
  _sw_balance = funcs[method]
  
  (; θ, zwt, snowpack) = soil
  Kc = [1.0, 1.0, 1.0]
  VPD, U2, doy = 0.0, 0.0, 0
  pEc, pEs, ET0 = potentialET(Rn, G, LAI, Ta, Pa, VPD, U2, doy; Kc) # PET allocated to canopy and soil surface
  Ei, fwet, PE = interception(P, pEc, LAI, pftpar)  # Interception evaporation

  # Snow sublimation, snow melt
  soil.snowpack, Esb, _, Pnet = snp_balance(PE, Ta, Tas, snowpack, pEs)

  RS, I, Vmax = runoff_up(Pnet, θ, zwt, Δz, soilpar)

  # Variables associated with soil water balance
  _pEs = max(pEs - Esb, 0)
  Tr, Es, uex = _sw_balance(soil, I, pEc, _pEs, Ta, Topt, fwet, s_VOD, soilpar, pftpar)

  ET = Tr + Es + Ei + Esb # Total Evapotranspiration
  RS += uex

  @pack! output = ET, Tr, Es, Ei, Esb, RS
  output.GW = soil.zwt
  output.SM .= soil.θ
  return nothing
  # return ET, Tr, Es, Ei, Esb, RS, Pnet, I, Vmax
end


function _run_model!(res::SpacOutputs{FT}, Rn::T, Ta::T, Tas::T, Prcp::T, Pa::T, G::T, LAI::T, s_VOD::T,
  Top::FT, soilpar, pftpar, soil::Soil; Kc=1.0, method="Kun") where {FT<:Real,T<:AbstractVector{FT}}
  ntime = size(Rn, 1) # 1365

  output = SpacOutput{FT}()
  for i in 1:ntime
    _Kc = 180 <= i <= 220 ? Kc : 1.0
    SiTHv2!(output, Rn[i], Ta[i], Tas[i], Top, Prcp[i], Pa[i], s_VOD[i], G[i], LAI[i],
      soilpar, pftpar, soil; Kc=_Kc, method)
    res[i] = output
  end
  return res
end


"""
## Arguments
- `Top`   : optimal growth temperature for plant, degC
- `soil`: 
  + `θ`    : Soil moisture (last step)
  + `zwt`   : groundwater table depth, mm
  + `snp`   : Snow package (old), mm day-1
- `spinfg`: spin-up flag, 1 for spin-up, 0 for normal calculation. 循环重复跑100次。
"""
function SiTHv2_site(Rn, Ta, Tas, Prcp, Pa, G, LAI, s_VOD,
  Top, soilpar, pftpar, soil, spin::Bool=false; Kc=1.0, method="Kun")

  ntime = length(Rn)
  res = SpacOutputs{Float64}(; ntime)

  if spin == 1 # spin-up
    for k in 1:100 # set the spin-up time (100 years)
      _run_model!(res, Rn, Ta, Tas, Prcp, Pa, G, LAI, s_VOD, Top, soilpar, pftpar, soil; Kc, method)
    end
  else
    _run_model!(res, Rn, Ta, Tas, Prcp, Pa, G, LAI, s_VOD, Top, soilpar, pftpar, soil; Kc, method)
  end

  (; ET, Tr, Es, Ei, Esb, RS, GW, SM) = res
  return ET, Tr, Es, Ei, Esb, RS, GW, SM
end
