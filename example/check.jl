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

begin
  using Plots

  plot(df_mat.ET, label="MATLAB")
  plot!(df_jl.ET, label="Julia")
end


