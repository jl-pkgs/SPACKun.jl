"""
## Model input
- Rn     : Net Radiation, W/m-2
- Ta     : Near surface air temperature, C
- Topt   : Optimal growth temperature for plant, C
- P      : Precipitation, mm day-1
- Pa     : Near surface air pressure, kPa
- s_VOD  : Constrains from VOD,[0,1]
- G      : Soil heat flux, W/m-2
- LAI    : Leaf area index

## intermediate variables in `soil`
- θ      : Soil moisture (last step)
- zwt    : Groundwater table depth, mm
- snp    : Snow package (old), mm day-1

## Model output
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
function SiTHv2!(output::SpacOutput{T}, soil::Soil,
  Rn::T, Ta::T, Tas::T, Topt::T, P::T, Pa::T,
  s_VOD::T, G::T, LAI::T,
  VPD::T=0.0, U2::T=0.0; # Monteith method
  doy::Int=1, method_SW="Kun", method_PET="PT72") where {T<:Real}

  (; θ, zwt, snowpack, param) = soil
  θ_sat = param.θ_sat[1]

  funcs = Dict(
    "Kun" => sw_balance,
    "CoLM" => sw_balance_CoLM)
  _sw_balance = funcs[method_SW]

  ## 这里还需要输入，VPD和风速数据
  pEc, pEs = cal_PET(Rn, G, LAI, Ta, Pa, VPD, U2, doy;
    param, method=method_PET)
  Ei, fwet, PE = interception(P, pEc, LAI, param.β)  # Interception evaporation

  # Snow sublimation, snow melt
  soil.snowpack, Esb, _, Pnet = snp_balance(PE, Ta, Tas, snowpack, pEs)

  RS, I, Vmax = runoff_up(Pnet, θ, zwt, Δz, θ_sat)

  # Variables associated with soil water balance
  _pEs = max(pEs - Esb, 0)
  Tr, Es, uex = _sw_balance(soil, I, pEc, _pEs, Ta, Topt, fwet, s_VOD)

  ET = Tr + Es + Ei + Esb # Total Evapotranspiration
  RS += uex

  @pack! output = ET, Tr, Es, Ei, Esb, RS
  output.GW = soil.zwt
  output.SM .= soil.θ
  return nothing
  # return ET, Tr, Es, Ei, Esb, RS, Pnet, I, Vmax
end


## WITH VPD and U2 as INPUT
function _run_model!(res::SpacOutputs{FT}, soil::Soil,
  Rn::V, Ta::V, Tas::V, Prcp::V, Pa::V,
  G::V, LAI::V, s_VOD::V, Top::FT,
  VPD::V, U2::V, doy::AbstractVector; kw...) where {FT<:Real,V<:AbstractVector{FT}}
  
  ntime = size(Rn, 1) # 1365
  output = SpacOutput{FT}()
  for i in 1:ntime
    SiTHv2!(output, soil, 
      Rn[i], Ta[i], Tas[i], Top, Prcp[i], Pa[i], s_VOD[i], G[i], LAI[i], VPD[i], U2[i]; 
      doy=doy[i], kw...)
    res[i] = output
  end
  return res
end

function _run_model!(res::SpacOutputs{FT}, soil::Soil,
  Rn::V, Ta::V, Tas::V, Prcp::V, Pa::V,
  G::V, LAI::V, s_VOD::V, Top::FT; kw...) where {FT<:Real,V<:AbstractVector{FT}}

  ntime = size(Rn, 1) # 1365
  output = SpacOutput{FT}()
  for i in 1:ntime
    SiTHv2!(output, soil,
      Rn[i], Ta[i], Tas[i], Top, Prcp[i], Pa[i], s_VOD[i], G[i], LAI[i]; kw...)
    res[i] = output
  end
  return res
end


"""
    SiTHv2_site(soil, Rn, Ta, Tas, Prcp, Pa, G, LAI, s_VOD, Top, args...; spin::Bool=false, kw...)

状态变量需要连续，传递到下一年中；模型对初始状态soil敏感。

- θ  : 采用warming-up的方式获取，warming-up period可取3年
- zwt: 采用spin-up的方式获取，spin-up period可取100年

## Arguments
- `Top`  : optimal growth temperature for plant, degC
- `soil` : 
  + `θ`     : Soil moisture (last step)
  + `zwt`   : groundwater table depth, mm
  + `snp`   : Snow package (old), mm day-1
- `spin` : spin-up flag, 1 for spin-up, 0 for normal calculation. 循环重复跑100次。
- `kw`   : other keyword arguments
  + `method_SW` : method for soil water balance, default "Kun"
  + `method_PET`: method for potential ET, default "PT72"
"""
function SiTHv2_site(soil, Rn, Ta, Tas, Prcp, Pa, G, LAI, s_VOD, Top, args...; spin::Bool=false, kw...)
  ntime = length(Rn)
  res = SpacOutputs{Float64}(; ntime)

  if spin == 1 # spin-up
    for k in 1:100 # set the spin-up time (100 years)
      _run_model!(res, soil, Rn, Ta, Tas, Prcp, Pa, G, LAI, s_VOD, Top, args...; kw...)
    end
  else
    _run_model!(res, soil, Rn, Ta, Tas, Prcp, Pa, G, LAI, s_VOD, Top, args...; kw...)
  end

  (; ET, Tr, Es, Ei, Esb, RS, GW, SM) = res
  return ET, Tr, Es, Ei, Esb, RS, GW, SM
end


export SiTHv2_site
