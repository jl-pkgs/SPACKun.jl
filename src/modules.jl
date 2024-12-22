"""
Interception evaporation

# INPUT:
- pEc    : potential Evaporation on canopy
- Rain   : precipitation
- LAI    : Leaf area index
- pftpar : PFT parameters

# OUTPUT:
- Ei     : Interception evaporation, mm/day
- fwet   : wetness fraction
- PE     : remaining precipitation after interception

# Reference:
Choudhury BJ, Digirolamo NE, 1998, A biological process-based estimate
of global land surface evaporation using satellite and ancillary data I.
Model description and comparison with observation. Journal of Hydrology.
"""
function interception(P, pEc, LAI, β)
  Sc = min(P, β * LAI * P)

  if pEc < 1e-3
    fwet = 1.0
    Ei = 0
  else
    fwet = min(0.7 * Sc / pEc, 1)
    Ei = pEc * fwet
  end
  PE = max(P - Ei, 0.0)
  return Ei, fwet, PE
end

# Interception (I, mm) was calculated using Horton's model:
#   I       =   min(P, aP + b)
# Where P   =   precipitation in mm
#  a, b     =   parameters with values derived from published records

"""
# INPUT:
  Δz       : Soil depth of different layers (mm)
  θ       : The antecedent soil water content (mm) expressed as a function of the WHC in that layer
  
  Pnet     : Net precipitation = P-I+Snowmelt
  zwt      : Groundwater table depth
  
# OUTPUT:
  srf      : Surface Runoff (mm)
  I      : Water that enters into the soil surface (mm)
  Vmax     : Maximum water retention capacity (mm)
"""
function runoff_up(Pnet, θ, zwt, Δz, θ_sat)
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

  I = clamp(Pnet - srf, 0, Vmax) # Actual water entering the soil surface (I)
  # Redundant water -> runoff (balance) under boundary conditions
  if I == Vmax || Vmax <= 0
    srf = Pnet - Vmax
  end
  return srf, I, Vmax
end


"""
Snowpack balance

- snowpack : available snow storage
- snowmelt : snow melt
- Esb      : snow sublimation
"""
function snowpack_balance(Prcp, Ta, Tas, snowpack, pEs)
  # Esnow_emp = 0.84 * (0.864 * (7.093 * Ta + 28.26)) / (Ta^2 - 3.593 * Ta + 5.175)
  Esnow = pEs  # Simple equivalent (Needs further development)

  # only snowfall occurs at Ta below zero
  if Ta <= 0
    # Add new snowfall, Ta <= 0
    newsnow = Prcp
    snowpack += newsnow

    snowmelt = 0  # snow melt
    Esb = clamp(snowpack, 0, Esnow)  # real snow sublimation

    Pnet = 0  # net Precipitation into soil surface
    snowpack -= Esb  # new snowpack
  else
    # real snow sublimation
    Esb = clamp(snowpack, 0, Esnow)

    snowpack -= Esb
    # snow melt, Ta > 0
    snowmelt_x = (1.5 + 0.007 * Prcp) * Tas  # Tas, accumulated Ta > 0
    snowmelt = min(snowpack, snowmelt_x)
    snowpack -= snowmelt

    Pnet = max(0, Prcp + snowmelt)  # net water into soil surface
  end
  return snowpack, Esb, snowmelt, Pnet
end
