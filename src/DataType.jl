export update_state

@with_kw mutable struct State{FT}
  wa::Vector{FT} = ones(3) .* 0.2
  zg::FT = 0.0
  snowpack::FT = 0.0
end

# Update state variables
function update_state!(state, wa, zg, snowpack)
  # Create a structure to store state variables
  #
  # 状态变量需要连续，传递到下一年中；模型对初始状态state敏感。
  # - wa: 采用warming-up的方式获取，warming-up period可取3年
  # - zg: 采用spin-up的方式获取，spin-up period可取100年
  #
  ## Argument Specification
  # - wa: soil water content in three layers, [m^3 m^-3]
  # - zg: groundwater depth, [mm]
  # - snowpack: snowpack depth
  wa[wa.<0] .= 0.01  # set the minimum value for soil moisture
  state.wa = wa
  state.zg = zg
  state.snowpack = snowpack
  return state
end

# # Get the state variables
# function get_state(state)
#   # Get the state variables
#   ## Argument Specification
#   # - State: a structure to store state variables
#   wa = state.wa
#   zg = state.zg
#   snowpack = state.snowpack
#   return wa, zg, snowpack
# end
