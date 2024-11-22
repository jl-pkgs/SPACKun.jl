# Soil moisture Constrains
"""
    swc_stress(θ::T, pET::T, soilpar, pftpar) where {T<:Real}

# INPUT
- `θ`      : The antecedent soil water content expressed as a function of the WHC in that layer
- `soilpar` : Soil parameters according to Soil type

# OUTPUT
- `S_plant` : Soil moisture stress to plant transpiration
- `S_soil`  : Soil moisture stress to soil evaporation
 
--------
Stress function for plant transpiration and soil evaporation:
                (θ_c-θ_wp)
wc =  ------------------------------    (about 0.4; Choudhury & Digirolamo, 1998)
                (θ_fc-θ_wp)
where θ_c : the critical soil water content at which plant stress start
Stree Function (Martens et al., 2017)
                    θ_c-θ
S_plant   =   1- (-------------------)^2   =   1-(1-w/wc)^2
                    θ_c-θ_wp

                  θ_c-θ
S_soil    =   1- -----------------    =   1-(1-w/wc)=w/wc
                  θ_c-θ_r
"""
function swc_stress(θ::T, pET::T, soilpar, pftpar) where {T<:Real}
  (; θ_fc, θ_wp) = soilpar
  (; Hc) = pftpar # canopy height, Zhang 2022

  k = Hc^0.5
  k = 4 * ((k - 0.7) / 4.3) + 1 # scale [1, 25] to [1, 5], `CH_scalar`

  b = 0.1
  p = 1 / (1 + pET) - b / (1 + Hc) # Zhang 2022, Eq. 9
  θ_wpCH = θ_wp / k

  # critical soil moisture for different PFTs
  θ_c = (1 - p) * (θ_fc - θ_wpCH) + θ_wpCH
  θ_c = clamp(θ_c, θ_wpCH, θ_fc)

  if θ <= θ_wpCH
    f_sm = 0.0
  elseif θ >= θ_c
    f_sm = 1.0
  else
    f_sm = 1.0 - ((θ_c - θ) / (θ_c - θ_wpCH))^k
  end

  # constraint for soil evaporation
  θ_wp_soil = 0
  if θ <= θ_wp_soil
    f_sm_s = 0.0
  elseif θ >= θ_fc
    f_sm_s = 1.0
  else
    # f_sm_s = ((θ - θ_wp) / (θ_fc - θ_wp))^1
    f_sm_s = (θ - θ_wp_soil) / (θ_fc - θ_wp_soil)  # for soil evaporation only
  end
  return f_sm, f_sm_s
end


function swc_stress!(pET::T, soilpar, pftpar, soil::State) where {T<:Real}
  (; θ_sat) = soilpar
  (; θ, fsm_Ec, fsm_Es, zwt, z₊ₕ, N) = soil

  jwt = find_jwt(z₊ₕ, zwt)
  j = jwt + 1

  # 非饱和、半饱和、饱和
  for i = 1:jwt # 全部非饱和
    fsm_Ec[i], fsm_Es[i] = swc_stress(θ[i], pET, soilpar, pftpar)
  end

  # j是半饱和
  if j <= N
    θ_unsat = find_θ_unsat(soil, θ_sat)[1] # 这里可能需要求，非饱和的比例（或深度）
    fsm_Ec[j], fsm_Es[j] = swc_stress(θ_unsat, pET, soilpar, pftpar)
  end

  # 全部饱和的部分，不受水分限制
  for i = j+1:N
    fsm_Ec[i] = 1.0
    fsm_Es[i] = 1.0
  end
end


# --- old version
# # wc = (θ_c - θ_wp) / (θ_fc - θ_wp)
# if θ <= θ_wp
#     f_sm = 0
# elseif θ >= θ_c
#     f_sm = 1
# else
# #     f_sm = 1 - (1 - θ / wc)^2
#     f_sm = 1 - ((θ_c - θ) / (θ_c - θ_wp))^2
# end
# --- old version

# f_sm_s = θ / wc
