export find_jwt

# 水位向下为正，地表为0
function find_jwt(z₊ₕ::AbstractVector, zwt::Real)
  N = length(z₊ₕ)
  zwt <= 0 && return 0
  zwt >= z₊ₕ[end] && return N + 1

  for j in 1:N
    zwt <= z₊ₕ[j] && return j
  end
end


"""
- 未饱和: z₊ₕ[jwt] ~ zwt
- 饱和  : zwt ~ z₊ₕ[jwt+1]
"""
function find_θ_unsat(θ, zwt; z₊ₕ, Δz, θ_sat)
  N = length(θ)
  j = find_jwt(z₊ₕ, zwt)

  j == 0 && return 0.0, 0.0         # 全部饱和
  j == N + 1 && return θ[N], 1.0    # 全部非饱和

  z0 = j == 1 ? 0 : z₊ₕ[j-1]
  z1 = z₊ₕ[j]

  d_unsat = zwt - z0
  θ_unsat = (θ[j] * Δz[j] - θ_sat * (z1 - zwt)) / d_unsat

  frac = d_unsat * θ_unsat / (θ[j] * Δz[j]) # 非饱和的比例
  # frac = d_unsat / Δz[j] # 非饱和的比例
  return θ_unsat, frac
end

function find_θ_unsat(soil::Soil, θ_sat)
  (; θ, zwt, z₊ₕ, Δz) = soil
  find_θ_unsat(θ, zwt; z₊ₕ, Δz, θ_sat)
end


function SM_recharge!(θ, I; Δz, θ_sat)
  N = length(θ)
  for i = 1:N
    if θ[i] > θ_sat # 提前处理超饱和
      I += (θ[i] - θ_sat) * Δz[i]
      θ[i] = θ_sat
    end

    _layer = clamp(I, 0, (θ_sat - θ[i]) * Δz[i]) # 剩余空间
    θ[i] += (_layer / Δz[i])
    I -= _layer
    # I <= 0 && break # 超饱和现象，不能提前退出
  end
  I # 形成到向下的通量
end


function GW_Rsb(zwt::Real)
  R_sb_max = 39.0 # mm day-1
  f = 1.25e-3     # mm-1
  return R_sb_max * exp(-f * zwt)
end


# sy = 0.02 # 给水度
function GW_update!(θ::Vector{T}, zwt::T, ∑::T, sy::T=0.2) where {T<:Real}
  uex = 0.0 # exceed water to soil surface
  if ∑ > 0 # 补给，水位上升
    ## 从地下水所在的水位开始补给
    jwt = find_jwt(z₊ₕ, zwt)
    for j in jwt+1:-1:1
      z0 = j == 1 ? 0.0 : z₊ₕ[j-1]
      _layer = clamp(∑, 0, sy * (zwt - z0)) # 剩余补给空间

      θ[j] += _layer / Δz[j]   # TODO: 土壤不能超饱和
      zwt -= _layer / sy      # 补给，水位上升
      ∑ -= _layer
      ∑ <= 0 && break
    end
    ∑ > 0 && (uex = ∑) # excess water to soil surface, [m]
  else # 排泄，水位下降
    for j = jwt+1:N
      z1 = z₊ₕ[j]
      _layer = clamp(∑, sy * (zwt - z1), 0) # 剩余排泄空间, Note: ∑ < 0

      θ[j] += _layer / Δz[j]   # 可能存在亏损
      zwt -= _layer / sy      # 需要注意单位转换
      ∑ -= _layer
      ∑ >= 0 && break
    end
    ∑ < 0 && (zwt -= ∑ / sy)
  end
end
