export find_jwt, GW_Rsb, SM_recharge!


include("find_θ_unsat.jl")
include("SM_discharge.jl")
include("SM_recharge.jl")
include("GW_update.jl")


# 水位向下为正，地表为0
function find_jwt(z₊ₕ::AbstractVector, zwt::Real)
  N = length(z₊ₕ)
  zwt <= 0 && return 0
  zwt >= z₊ₕ[end] && return N + 1

  for j in 1:N
    zwt <= z₊ₕ[j] && return j
  end
end


function GW_Rsb(zwt::Real)
  R_sb_max = 39.0 # mm day-1
  f = 1.25e-3     # mm-1
  return R_sb_max * exp(-f * zwt)
end
