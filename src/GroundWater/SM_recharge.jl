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
