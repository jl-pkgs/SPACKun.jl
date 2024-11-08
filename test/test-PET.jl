using SPAC, Test

@testset "ET0_Penman48" begin
  Rn = 200.0    # W/m²
  Ta = 25.0     # °C
  VPD = 1.5     # kPa
  Uz = 3.0      # m/s
  
  @test ET0_Penman48(Rn, Ta, VPD, Uz) == 7.9265544809634365
  @test ET0_PT1972(Rn, Ta) == 6.564860193238437
end
