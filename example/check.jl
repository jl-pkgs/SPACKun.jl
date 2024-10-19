using MATLAB

function mat_start(show=false)
  s = MSession()    # creates a MATLAB session
  # similarily
  # show ? mat_show(s) : mat_hide(s)
  s
end
s = mat_start();

mat"""
cd Z:/GitHub/cug-hydro/SiTHv2
install
"""
mxcall(:potentialET, 2, Rn, G, LAI, Ta, Pa)

function plot_var(var = :ET)
  plot(df_mat[:, var], label="MATLAB", title=string(var))
  plot!(df_jl[:, var], label="Julia")
end

begin
  using Plots
  gr(; framestyle = :box)
  plot(
    plot_var(:ET), 
    plot_var(:Tr), 
    plot_var(:Es), 
    plot_var(:Ei), 
    plot_var(:Esb), 
    plot_var(:RF), 
    plot_var(:GW), 
    plot_var(:SM1), 
    plot_var(:SM2), 
    plot_var(:SM3), 
    size = (1400, 800)
  )
end
