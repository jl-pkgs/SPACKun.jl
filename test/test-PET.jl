using SPAC, Test

@testset "ET0 models" begin
  Rn = 200.0    # W/m²
  Ta = 25.0     # °C
  VPD = 1.5     # kPa
  Uz = 3.0      # m/s
  
  @test ET0_Penman48(Rn, Ta, VPD, Uz) == 7.9265544809634365
  @test ET0_PT1972(Rn, Ta) == 6.564860193238437
  @test ET0_FAO98(Rn, Ta, VPD, Uz) == 6.928646397419433

  @test aerodynamic_resistance(1.0) == 207.66407000788683 # s/m
  @test ET0_Monteith65(250.0, 25.0, 1.0, 2.0; hc=0.12) == 6.871881995698257
end

# ET0_Monteith65(250.0, 25.0, 1.0, 2.0; hc=0.7) == 7.2115556355648955
# ET0_Monteith65(250.0, 25.0, 1.0, 2.0; hc=2.0)

# using Plots
# β = 0:0.05:1.0
# ET = ET0_Monteith65.(250.0, 25.0, 1.0, 2.0, atm, β; hc=0.12) 
# plot(β, PET, label="hc=0.12", lw=2, xlabel="β", ylabel="PET (mm/d)", title="ET0_Monteith65")
