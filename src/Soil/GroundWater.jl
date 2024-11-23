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


function GW_Rsb(zwt::Real)
  R_sb_max = 39.0 # mm day-1
  f = 1.25e-3     # mm-1
  return R_sb_max * exp(-f * zwt)
end


# sy = 0.02 # 给水度
# TODO: jwt有待矫正
function GW_update!(θ::Vector{T}, zwt::T, ∑::T, sy::T=0.2) where {T<:Real}
  
  uex = 0.0 # exceed water to soil surface
  if ∑ > 0 # 补给，水位上升
    ## 从地下水所在的水位开始补给
    jwt = find_jwt(z₊ₕ, zwt)
    # TODO: 如果是N+1层开始补给
    for j in jwt:-1:1
      z0 = j == 1 ? 0.0 : z₊ₕ[j-1]
      # z1 = j == N+1 ? zwt : z₊ₕ[j]
      _layer = clamp(∑, 0, sy * (zwt - z0)) # 剩余补给空间
      # N+1层不计算θ[j]
      j <= N && (θ[j] += _layer / Δz[j])   # TODO: 土壤不能超饱和
      zwt -= _layer / sy      # 补给，水位上升
      ∑ -= _layer
      ∑ <= 0 && break
    end
    ∑ > 0 && (uex = ∑) # excess water to soil surface, [m]
  else # 排泄，水位下降
    for j = jwt:N
      # TODO: 这里没有考虑jwt = 0的情况
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
