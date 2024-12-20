# 
"""
Partition PET into three soil layers

# INPUT:
- pEc    :  Potential Evaporation on canopy
- w      :  Initialized values for soil moisture
- pftpar :  PFT parameters

# OUTPUT:
- Tr_p :  separate potential Transpiration
"""
function pTr_partition!(soil::Soil, pEc::T, fwet::T) where {T<:Real}
  (; θ, Δz, Ec_pot) = soil
  
  b = soil.param.b[1]
  θ_sat = soil.param.θ_sat[1]
  D50 = soil.param.D50[1]
  D95 = soil.param.D95[1]
  
  c = -2.944 / log(D95 / D50)

  r1 = (1 / (1 + (Δz[1] / D50)^c)) # Zhang 2019, Eq. 21, root depths function
  r2 = (1 / (1 + (Δz[2] / D50)^c)) - (1 / (1 + (Δz[1] / D50)^c))
  r3 = (1 / (1 + (Δz[3] / D50)^c)) - (1 / (1 + (Δz[2] / D50)^c))

  # the maximum transpiration rate of each soil layer, Tr_p
  # Get all available water contents through root distribution
  wr = r1 * (θ[1] / θ_sat)^b + r2 * (θ[2] / θ_sat)^b + r3 * (θ[3] / θ_sat)^b

  # Root distribution adjusted by soil water content
  β1 = r1 * (θ[1] / θ_sat)^b / wr
  β2 = r2 * (θ[2] / θ_sat)^b / wr
  β3 = r3 * (θ[3] / θ_sat)^b / wr

  # potential transpiration rate for different layers
  Ec_pot[1] = (1 - fwet) * β1 * pEc
  Ec_pot[2] = (1 - fwet) * β2 * pEc
  Ec_pot[3] = (1 - fwet) * β3 * pEc
end
