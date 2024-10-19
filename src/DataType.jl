export update_state

@with_kw mutable struct State{FT}
  wa::Vector{FT} = ones(3) .* 0.2
  zgw::FT = 0.0
  snowpack::FT = 0.0
end

# Update state variables
function update_state!(state, wa, zgw, snowpack)
  # Create a structure to store state variables
  #
  # 状态变量需要连续，传递到下一年中；模型对初始状态state敏感。
  # - wa: 采用warming-up的方式获取，warming-up period可取3年
  # - zgw: 采用spin-up的方式获取，spin-up period可取100年
  #
  ## Argument Specification
  # - wa: soil water content in three layers, [m^3 m^-3]
  # - zgw: groundwater depth, [mm]
  # - snowpack: snowpack depth
  wa[wa.<0] .= 0.01  # set the minimum value for soil moisture
  state.wa = wa
  state.zgw = zgw
  state.snowpack = snowpack
  return state
end
