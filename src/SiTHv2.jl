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
- wa     : Soil moisture (last step)
- zgw    : Groundwater table depth, mm
- snp    : Snow package (old), mm day-1

# Model output
- Et     : Total Evapotranspiration, mm day-1
- Tr     : Plant Transpiration, mm day-1
- Es     : Soil Evaporation, mm day-1
- Ei     : Intercepted Evaporation, mm day-1
- Esb    : Snow sublimation, mm day-1
- wa     : Soil moisture (three layers)
- srf    : Surface runoff, mm day-1
- zgw    : Groundwater table depth, mm
- snp    : Snow package (new), mm day-1

"""
function SiTHv2(Rn, Ta, Tas, Topt, P, Pa, s_VOD, G, LAI, soilpar, pftpar, state::State; Kc=1.0)
  (; wa, zgw, snowpack) = state
  Kc = [1.0, 1.0, 1.0]
  VPD, U2, doy = 0.0, 0.0, 0
  pEc, pEs, ET0 = potentialET(Rn, G, LAI, Ta, Pa, VPD, U2, doy; Kc) # PET allocated to canopy and soil surface
  Ei, fwet, PE = interception(P, pEc, LAI, pftpar)  # Interception evaporation

  # Snow sublimation, snow melt
  state.snowpack, Esb, _, Pnet = snp_balance(PE, Ta, Tas, snowpack, pEs)

  srf, IWS, Vmax = runoff_up(Pnet, wa, zgw, ZM, soilpar)

  # Variables associated with soil water balance
  new_pEs = max(pEs - Esb, 0)
  Tr, Es, uex = sw_balance(IWS, pEc, new_pEs, Ta, Topt, s_VOD, soilpar, pftpar, fwet, ZM, state)

  Et = Tr + Es + Ei + Esb # Total Evapotranspiration
  srf += uex
  return Et, Tr, Es, Ei, Esb, srf, Pnet, IWS, Vmax
end


function _run_model(Rn::T, Ta::T, Tas::T, Prcp::T, Pa::T, G::T, LAI::T, s_VOD::T,
  Top::FT, soilpar, pftpar, state::State, output; Kc=1.0) where {FT<:Real,T<:AbstractVector{FT}}
  ntime = size(Rn, 1) # 1365

  for i in 1:ntime
    _Kc = 180 <= i <= 220 ? Kc : 1.0
    Et, Tr, Es, Ei, Esb, srf, _, _, _ = SiTHv2(Rn[i], Ta[i], Tas[i], Top, Prcp[i], Pa[i], s_VOD[i], G[i], LAI[i], soilpar, pftpar, state; Kc=_Kc)

    output[:ETs][i] = Et
    output[:Trs][i] = Tr
    output[:Ess][i] = Es
    output[:Eis][i] = Ei
    output[:Esbs][i] = Esb
    output[:SM][i, :] = state.wa
    output[:RF][i] = srf
    output[:GW][i] = state.zgw
  end
  return output
end


"""
## Arguments
- `Top`   : optimal growth temperature for plant, degC
- `state`: 
  + `wa`    : Soil moisture (last step)
  + `zgw`   : groundwater table depth, mm
  + `snp`   : Snow package (old), mm day-1
- `spinfg`: spin-up flag, 1 for spin-up, 0 for normal calculation. 循环重复跑100次。
"""
function SiTHv2_site(Rn, Ta, Tas, Prcp, Pa, G, LAI, s_VOD,
  Top, soilpar, pftpar, state, spin::Bool=false; Kc=1.0)

  ntime = length(Rn)
  output = Dict(
    :ETs => zeros(ntime),
    :Trs => zeros(ntime),
    :Ess => zeros(ntime),
    :Eis => zeros(ntime),
    :Esbs => zeros(ntime),
    :SM => zeros(ntime, 3),
    :RF => zeros(ntime),
    :GW => zeros(ntime)
  )
  if spin == 1 # spin-up
    for k in 1:100 # set the spin-up time (100 years)
      output = _run_model(Rn, Ta, Tas, Prcp, Pa, G, LAI, s_VOD, Top, soilpar, pftpar, state, output; Kc)
    end
  else
    output = _run_model(Rn, Ta, Tas, Prcp, Pa, G, LAI, s_VOD, Top, soilpar, pftpar, state, output; Kc)
  end

  ETs = output[:ETs]
  Trs = output[:Trs]
  Ess = output[:Ess]
  Eis = output[:Eis]
  Esbs = output[:Esbs]
  SM = output[:SM]
  RF = output[:RF]
  GW = output[:GW]
  return ETs, Trs, Ess, Eis, Esbs, SM, RF, GW
end
