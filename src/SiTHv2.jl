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
  G::T, LAI::T, soilpar, pftpar, state::State; Kc=1.0) where {T<:Real}

  (; θ, zwt, snowpack) = state
  Kc = [1.0, 1.0, 1.0]
  VPD, U2, doy = 0.0, 0.0, 0
  pEc, pEs, ET0 = potentialET(Rn, G, LAI, Ta, Pa, VPD, U2, doy; Kc) # PET allocated to canopy and soil surface
  Ei, fwet, PE = interception(P, pEc, LAI, pftpar)  # Interception evaporation

  # Snow sublimation, snow melt
  state.snowpack, Esb, _, Pnet = snp_balance(PE, Ta, Tas, snowpack, pEs)

  RS, IWS, Vmax = runoff_up(Pnet, θ, zwt, Δz, soilpar)

  # Variables associated with soil water balance
  new_pEs = max(pEs - Esb, 0)
  Tr, Es, uex = sw_balance(IWS, pEc, new_pEs, Ta, Topt, s_VOD, soilpar, pftpar, fwet, Δz, state)

  ET = Tr + Es + Ei + Esb # Total Evapotranspiration
  RS += uex

  @pack! output = ET, Tr, Es, Ei, Esb, RS
  output.GW = state.zwt
  output.SM .= state.θ
  return nothing
  # return ET, Tr, Es, Ei, Esb, RS, Pnet, IWS, Vmax
end


function _run_model!(res::SpacOutputs{FT}, Rn::T, Ta::T, Tas::T, Prcp::T, Pa::T, G::T, LAI::T, s_VOD::T,
  Top::FT, soilpar, pftpar, state::State; Kc=1.0) where {FT<:Real,T<:AbstractVector{FT}}
  ntime = size(Rn, 1) # 1365

  output = SpacOutput{FT}()
  for i in 1:ntime
    _Kc = 180 <= i <= 220 ? Kc : 1.0
    SiTHv2!(output, Rn[i], Ta[i], Tas[i], Top, Prcp[i], Pa[i], s_VOD[i], G[i], LAI[i],
      soilpar, pftpar, state; Kc=_Kc)
    res[i] = output
  end
  return res
end


"""
## Arguments
- `Top`   : optimal growth temperature for plant, degC
- `state`: 
  + `θ`    : Soil moisture (last step)
  + `zwt`   : groundwater table depth, mm
  + `snp`   : Snow package (old), mm day-1
- `spinfg`: spin-up flag, 1 for spin-up, 0 for normal calculation. 循环重复跑100次。
"""
function SiTHv2_site(Rn, Ta, Tas, Prcp, Pa, G, LAI, s_VOD,
  Top, soilpar, pftpar, state, spin::Bool=false; Kc=1.0)

  ntime = length(Rn)
  res = SpacOutputs{Float64}(; ntime)

  if spin == 1 # spin-up
    for k in 1:100 # set the spin-up time (100 years)
      _run_model!(res, Rn, Ta, Tas, Prcp, Pa, G, LAI, s_VOD, Top, soilpar, pftpar, state; Kc)
    end
  else
    _run_model!(res, Rn, Ta, Tas, Prcp, Pa, G, LAI, s_VOD, Top, soilpar, pftpar, state; Kc)
  end

  (; ET, Tr, Es, Ei, Esb, RS, GW, SM) = res
  return ET, Tr, Es, Ei, Esb, RS, GW, SM
end
