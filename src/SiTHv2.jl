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
function SiTHv2(Rn, Ta, Tas, Topt, Pe, Pa, s_VOD, G, LAI, soilpar, pftpar, state::State)
  (; wa, zgw, snowpack) = state
  # Potential Evaporation allocated to canopy and soil surface
  pEc, pEs = potentialET(Rn, G, LAI, Ta, Pa)

  # Interception evaporation
  Ei, fwet, _ = interception(Pe, pEc, LAI, pftpar)

  # Snow sublimation, snow melt
  new_Pe = max(Pe - Ei, 0)
  snowpack, Esb, _, Pnet = snp_balance(new_Pe, Ta, Tas, snowpack, pEs)

  srf, IWS, Vmax = runoff_up(Pnet, zgw, ZM, wa, soilpar)

  # Variables associated with soil water balance
  new_pEs = max(pEs - Esb, 0)
  wa, zgw, Tr, Es, uex = sw_balance(IWS, pEc, new_pEs, Ta, Topt, s_VOD, soilpar, pftpar, fwet, ZM, wa, zgw)

  # Total Evapotranspiration
  Et = Tr + Es + Ei + Esb
  srf += uex

  update_state!(state, wa, zgw, snowpack)
  return Et, Tr, Es, Ei, Esb, srf, Pnet, IWS, Vmax
end


function _run_model(Rni::T, Tai::T, Tasi::T, Prcp::T, Pai::T, Gi::T, LAIii::T, s_VODi::T,
  Top::FT, soilpar, pftpar, state::State, output) where {FT<:Real,T<:AbstractVector{FT}}
  ntime = size(Rni, 1)

  for i in 1:ntime
    Rn = Rni[i]
    Ta = Tai[i]
    Tas = Tasi[i]
    Pe = Prcp[i]
    Pa = Pai[i]
    G = Gi[i]
    LAI = LAIii[i]
    s_VOD = s_VODi[i]

    Et, Tr, Es, Ei, Esb, srf, _, _, _ = SiTHv2(Rn, Ta, Tas, Top, Pe, Pa, s_VOD, G, LAI, soilpar, pftpar, state)

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
  Top, soilpar, pftpar, state, spin::Bool=false)

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
      output = _run_model(Rn, Ta, Tas, Prcp, Pa, G, LAI, s_VOD, Top, soilpar, pftpar, state, output)
    end
  else
    output = _run_model(Rn, Ta, Tas, Prcp, Pa, G, LAI, s_VOD, Top, soilpar, pftpar, state, output)
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
