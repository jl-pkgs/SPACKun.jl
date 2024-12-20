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
