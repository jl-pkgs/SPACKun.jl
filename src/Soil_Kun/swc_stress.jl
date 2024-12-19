# Soil moisture Constrains
"""
    swc_stress(θ::T, pET::T, soilpar, pftpar) where {T<:Real}

# INPUT
- `θ`      : The antecedent soil water content expressed as a function of the
  WHC in that layer
- `soilpar`: Soil parameters according to Soil type

# OUTPUT
- `S_plant`: Soil moisture stress to plant transpiration
- `S_soil` : Soil moisture stress to soil evaporation
 
--------
Stress function for plant transpiration and soil evaporation:

wc =  \frac{(θ_c-θ_wp)}{(θ_fc-θ_wp)} % (about 0.4; Choudhury & Digirolamo, 1998)

where `θ_c` : the critical soil water content at which plant stress start Stree
Function (Martens et al., 2017)

```math
S_plant = 1 - \frac{θ_c-θ}{θ_c-θ_wp}^2 = 1-(1-w/wc)^2
S_soil  = 1 - \frac{θ_c-θ}{θ_c-θ_r}    = 1-(1-w/wc)=w/wc
```
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


function constrain_sm(θ; θ_fc, θ_wp)
  θ <= θ_wp && (fsm = 0.0)
  θ >= θ_fc && (fsm = 1.0)
  (θ - θ_wp) / (θ_fc - θ_wp)
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
