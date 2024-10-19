# Soil moisture Constrains
function swc_stress(wa, soilpar, pET, pftpar)
  # INPUT:
  # wa      : The antecedent soil water content expressed
  #           as a function of the WHC in that layer
  # soilpar : Soil parameters according to Soil type
  # OUTPUT:
  # S_plant : Soil moisture stress to plant transpiration
  # S_soil  : Soil moisture stress to soil evaporation
  # --------
  # Stress function for plant transpiration and soil evaporation:
  #                  (θ_c-θ_wp)
  # wc =  ------------------------------    (about 0.4; Choudhury & Digirolamo, 1998)
  #                  (θ_fc-θ_wp)
  # where θ_c : the critical soil water content at which plant stress start
  # Stree Function (Martens et al., 2017)
  #                      θ_c-θ
  # S_plant   =   1- (-------------------)^2   =   1-(1-w/wc)^2
  #                     θ_c-θ_wp
  #
  #                    θ_c-θ
  # S_soil    =   1- -----------------    =   1-(1-w/wc)=w/wc
  #                    θ_c-θ_r
  # -------------------------------------------------------------------------
  (; θ_fc, θ_wp) = soilpar
  (; Zr) = pftpar # root depth

  k = Zr^0.5
  # scale [1, 25] to [1, 5]
  k = 4 * ((k - 0.7) / 4.3) + 1

  a = 0.1
  p = 1 / (1 + pET) - a * (1 / (1 + Zr))

  θ_wpCH = θ_wp / k

  # critical soil moisture for different PFTs
  θ_c = (1 - p) * (θ_fc - θ_wpCH) + θ_wpCH
  θ_c = clamp(θ_c, θ_wpCH, θ_fc)

  if wa <= θ_wpCH
    f_sm = 0
  elseif wa >= θ_c
    f_sm = 1
  else
    f_sm = 1 - ((θ_c - wa) / (θ_c - θ_wpCH))^k
  end

  # constraint for soil evaporation
  θ_wp_soil = 0
  if wa <= θ_wp_soil
    f_sm_s = 0
  elseif wa >= θ_fc
    f_sm_s = 1
  else
    # f_sm_s = ((wa - θ_wp) / (θ_fc - θ_wp))^1
    f_sm_s = (wa - θ_wp_soil) / (θ_fc - θ_wp_soil)  # for soil evaporation only
  end

  return f_sm, f_sm_s
end

# --- old version
# # wc = (θ_c - θ_wp) / (θ_fc - θ_wp)
# if wa <= θ_wp
#     f_sm = 0
# elseif wa >= θ_c
#     f_sm = 1
# else
# #     f_sm = 1 - (1 - wa / wc)^2
#     f_sm = 1 - ((θ_c - wa) / (θ_c - θ_wp))^2
# end
# --- old version

# f_sm_s = wa / wc
