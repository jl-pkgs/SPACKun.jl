"""
Interception evaporation

# INPUT:
- pEc    : potential Evaporation on canopy
- Rain   : precipitation
- LAI    : Leaf area index
- pftpar : PFT parameters

# OUTPUT:
- Ei     : Interception evaporation, mm/day
- wetT   : wetness fraction
- Pe     : remaining precipitation after interception

# Reference:
Choudhury BJ, Digirolamo NE, 1998, A biological process-based estimate
of global land surface evaporation using satellite and ancillary data I.
Model description and comparison with observation. Journal of Hydrology.
"""
function interception(pEc, LAI, Rain, pftpar)
  # Interception (I, mm) was calculated using Horton's model:
  #   I       =   min(P, aP + b)
  # Where P   =   precipitation in mm
  #  a, b     =   parameters with values derived from published records

  # x, parameter values, used for calculating daily rainfall interception
  inc = pftpar[1]
  Sc = min(Rain, inc * LAI * Rain)
  wetT = min(0.7 * Sc / pEc, 1)

  if pEc < 1e-3
    Ei = 0
  else
    Ei = pEc * wetT
  end

  Pe = Rain - Ei
  return Ei, wetT, Pe
end
