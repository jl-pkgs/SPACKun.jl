using SPAC, Plots
gr(framestyle=:box)


# layer2
Dmin = 0.012; # 0.0005*24, mm day-1
Dmax = 1.2; # 0.05*24,   mm day-1

# Reference
# 1. Ducharne & Polcher, 1998


begin

  function plot_figure1(; Dmax=4.8, Dmin=0.048, kw...)
    se = 0.0:0.01:1.0
    Ksat = 1.0 # [mm/d]
    drainage = soil_drainage.(se, Ksat, Dmin, Dmax; dd=1.5)
    plot(se, drainage; kw...)
  end

  p1 = plot_figure1(; Dmax=4.8, Dmin=0.048, title="(a) shadow soil", xlabel="se", ylabel="Drainage [mm/d]")
  p2 = plot_figure1(; Dmax=1.2, Dmin=0.012, title="(b) deep soil", xlabel="se", ylabel="Drainage [mm/d]")
  plot(p1, p2; size=(800, 400), margin=(3.0, :mm))
end
savefig("soil_drainage.png")
