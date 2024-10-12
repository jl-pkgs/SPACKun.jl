"""
## Arguments
- `Top`   : optimal growth temperature for plant, degC
- `state`: 
  + `wa`    : Soil moisture (last step)
  + `zgw`   : groundwater table depth, mm
  + `snp`   : Snow package (old), mm day-1
- `spinfg`: spin-up flag, 1 for spin-up, 0 for normal calculation. 循环重复跑100次。
"""
function cal_SiTHv2_site(Rni, Tai, Tasi, Precii, Pai, Gi, LAIii, Top, s_VODi, soilpar, pftpar, state, spinfg)
  ntime = size(Rni, 1)

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

  if spinfg == 1 # spin-up
    for k in 1:100 # set the spin-up time (100 years)
      state, output = run_model(Rni, Tai, Tasi, Precii, Pai, Gi, LAIii, Top, s_VODi, soilpar, pftpar, state, output)
    end
  else
    state, output = run_model(Rni, Tai, Tasi, Precii, Pai, Gi, LAIii, Top, s_VODi, soilpar, pftpar, state, output)
  end

  ETs = output[:ETs]
  Trs = output[:Trs]
  Ess = output[:Ess]
  Eis = output[:Eis]
  Esbs = output[:Esbs]
  SM = output[:SM]
  RF = output[:RF]
  GW = output[:GW]
  return ETs, Trs, Ess, Eis, Esbs, SM, RF, GW, state
end


function run_model(Rni, Tai, Tasi, Precii, Pai, Gi, LAIii, Top, s_VODi, soilpar, pftpar, state, output)
  wa, zgw, snp = mSiTH.get_state(state)
  ntime = size(Rni, 1)

  for i in 1:ntime
    Rn = Rni[i]
    Ta = Tai[i]
    Tas = Tasi[i]
    Pe = Precii[i]
    Pa = Pai[i]
    G = Gi[i]
    LAI = LAIii[i]
    s_VOD = s_VODi[i]

    Et, Tr, Es, Ei, Esb, wa, srf, zgw, snp, _, _, _ = SiTHv2(Rn, Ta, Tas, Top, Pe, Pa, s_VOD, G, LAI, soilpar, pftpar, wa, zgw, snp)

    output[:ETs][i] = Et
    output[:Trs][i] = Tr
    output[:Ess][i] = Es
    output[:Eis][i] = Ei
    output[:Esbs][i] = Esb
    output[:SM][i, :] = wa
    output[:RF][i] = srf
    output[:GW][i] = zgw
  end
  state = mSiTH.update_state(state, wa, zgw, snp)
  return state, output
end
