"""
    soil_drainage(wa_unsat, θ_sat, ks, Dmin, Dmax; dd = 1.5)

# layer1
Dmin = 0.048; # mm day-1
Dmax = 4.8; # mm day-1

# layer2
Dmin = 0.012; # 0.0005*24, mm day-1
Dmax = 1.2; # 0.05*24,   mm day-1

# Reference
1. Ducharne & Polcher, 1998
"""
function soil_drainage(θ_unsat, θ_sat, ks, Dmin, Dmax; dd=1.5)
  thx = θ_unsat / θ_sat

  if thx < 0.75
    perc = Dmin * thx
  else
    perc = Dmin * thx + (Dmax - Dmin) * thx^dd
  end
  return min(ks, perc) # drainage from unsaturated zone
end
