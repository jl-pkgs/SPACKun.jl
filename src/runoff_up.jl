"""
# INPUT:
  Δz       : Soil depth of different layers (mm)
  θ       : The antecedent soil water content (mm) expressed as a function of the WHC in that layer
  soilpar  : Soil parameters according to Soil type
  Pnet     : Net precipitation = P-I+Snowmelt
  zwt      : Groundwater table depth
  
# OUTPUT:
  srf      : Surface Runoff (mm)
  IWS      : Water that enters into the soil surface (mm)
  Vmax     : Maximum water retention capacity (mm)
"""
function runoff_up(Pnet, θ, zwt, Δz, soilpar)
  (; θ_sat) = soilpar

  Vmax = 0.0
  if zwt <= 0
    # Exceeded groundwater on the soil surface
    Vmax = 0.0
    srf = -zwt * θ_sat
  else
    if 0 < zwt <= Δz[1]
      d1 = zwt  # Thickness of unsaturated soil (mm)
      wa1_unsat = (θ[1] * Δz[1] - θ_sat * (Δz[1] - d1)) / d1
      # soil water retention capacity (Vmax), 剩余持水能力
      Vmax = d1 * (θ_sat - wa1_unsat) 
    elseif Δz[1] < zwt <= z₊ₕ[2]
      d1 = Δz[1]
      d2 = zwt - z₊ₕ[1]
      wa1_unsat = θ[1]
      wa2_unsat = (θ[2] * Δz[2] - θ_sat * (Δz[2] - d2)) / d2
      Vmax = d1 * (θ_sat - wa1_unsat) + d2 * (θ_sat - wa2_unsat)
    elseif zwt > z₊ₕ[2]
      d1 = Δz[1]
      d2 = Δz[2]
      wa1_unsat = θ[1]
      wa2_unsat = θ[2]
      Vmax = d1 * (θ_sat - wa1_unsat) + d2 * (θ_sat - wa2_unsat)
    end

    # Calculate surface runoff (srf)
    if Pnet > 0.2 * Vmax
      srf = (Pnet - 0.2 * Vmax)^2 / (Pnet + 0.8 * Vmax) # Zhang 2019, Eq. 6
    else
      srf = 0
    end
  end

  IWS = clamp(Pnet - srf, 0, Vmax) # Actual water entering the soil surface (IWS)
  # Redundant water -> runoff (balance) under boundary conditions
  if IWS == Vmax || Vmax <= 0
    srf = Pnet - Vmax
  end
  return srf, IWS, Vmax
end
