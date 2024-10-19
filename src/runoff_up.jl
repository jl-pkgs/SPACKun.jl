"""
# INPUT:
  ZM       : Soil depth of different layers (mm)
  wa       : The antecedent soil water content (mm)
             expressed as a function of the WHC in that layer
  soilpar  : Soil parameters according to Soil type
  Pnet     : Net precipitation = P-I+Snowmelt
  zgw      : Groundwater table depth
  
# OUTPUT:
  srf      : Surface Runoff (mm)
  IWS      : Water that enters into the soil surface (mm)
  Vmax     : Maximum water retention capacity (mm)
"""
function runoff_up(Pnet, zgw, ZM, wa, soilpar)
  theta_sat = soilpar[3]  # saturated soil moisture

  Vmax = 0.0
  if zgw <= 0
    # Exceeded groundwater on the soil surface
    Vmax = 0.0
    srf = -zgw * theta_sat
    Pnet = 0.0
  else
    if zgw > 0 && zgw <= ZM[1]
      d1 = zgw  # Thickness of unsaturated soil (mm)
      # Unsaturated soil water (mm)
      wa1_unsat = (wa[1] * ZM[1] - theta_sat * (ZM[1] - d1)) / d1
      # Calculate overall soil water retention capacity (Vmax)
      Vmax = d1 * (theta_sat - wa1_unsat)
    elseif zgw > ZM[1] && zgw <= ZM[1] + ZM[2]
      d1 = ZM[1]
      d2 = zgw - ZM[1]
      # Unsaturated soil water (mm)
      wa1_unsat = wa[1]
      wa2_unsat = (wa[2] * ZM[2] - theta_sat * (ZM[2] - d2)) / d2
      # Calculate overall soil water retention capacity (Vmax)
      Vmax = d1 * (theta_sat - wa1_unsat) + d2 * (theta_sat - wa2_unsat)
    elseif zgw > ZM[1] + ZM[2]
      d1 = ZM[1]
      d2 = ZM[2]
      # Unsaturated soil water (mm)
      wa1_unsat = wa[1]
      wa2_unsat = wa[2]
      # Calculate overall soil water retention capacity (Vmax)
      Vmax = d1 * (theta_sat - wa1_unsat) + d2 * (theta_sat - wa2_unsat)
    end

    # Calculate surface runoff (srf)
    if Pnet > 0.2 * Vmax
      srf = (Pnet - 0.2 * Vmax)^2 / (Pnet + 0.8 * Vmax)
    else
      srf = 0
    end
  end

  # Actual water entering the soil surface (IWS)
  IWS = min(Vmax, Pnet - srf)
  IWS = max(0, IWS)

  # Redundant water -> runoff (balance) under boundary conditions
  if IWS == Vmax || Vmax <= 0
    srf = Pnet - Vmax
  end
  return srf, IWS, Vmax
end
