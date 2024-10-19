export cal_SiTHv2_site;

"""
## Arguments
- `Top`   : optimal growth temperature for plant, degC
- `state`: 
  + `wa`    : Soil moisture (last step)
  + `zgw`   : groundwater table depth, mm
  + `snp`   : Snow package (old), mm day-1
- `spinfg`: spin-up flag, 1 for spin-up, 0 for normal calculation. 循环重复跑100次。
"""
function cal_SiTHv2_site(Rn, Ta, Tas, Prcp, Pa, G, LAI, s_VOD,
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
      output = run_model(Rn, Ta, Tas, Prcp, Pa, G, LAI, s_VOD, Top, soilpar, pftpar, state, output)
    end
  else
    output = run_model(Rn, Ta, Tas, Prcp, Pa, G, LAI, s_VOD, Top, soilpar, pftpar, state, output)
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


function run_model(Rni::T, Tai::T, Tasi::T, Prcp::T, Pai::T, Gi::T, LAIii::T, s_VODi::T,
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
    output[:SM][i, :] = state.sm
    output[:RF][i] = srf
    output[:GW][i] = state.zg
  end
  return output
end
