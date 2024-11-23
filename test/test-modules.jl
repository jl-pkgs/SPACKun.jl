@testset "potentialET" begin
  Rn = 100.0
  G = 5.0
  LAI = 3.0
  Ta = 20.0
  Pa = 100.0

  VPD, U2, doy = 0.0, 0.0, 0
  pEc, pEs, ET0 = potentialET(Rn, G, LAI, Ta, Pa, VPD, U2, doy)

  @test pEc ≈ 2.5389604966643824
  @test pEs ≈ 0.3507115521629298
end

@testset "runoff_up" begin
  Pnet = 20.0 # mm
  zwt = 1000.0 # mm
  θ = [0.3, 0.3, 0.3]
  soilpar = get_soilpar(2)
  srf, I, Vmax = runoff_up(Pnet, θ, zwt, Δz, soilpar)
  @test Vmax ≈ 129
  @test I ≈ 20
end

@testset "pTr_partition" begin
  pEc = 20.0
  soil = Soil()
  soil.θ = [0.3, 0.2, 0.1]

  fwet = 0.5
  soilpar = get_soilpar(2) #|> collect
  pftpar = get_pftpar(22) #|> collect
  pTr_partition!(soil, pEc, fwet, soilpar, pftpar)
  @test soil.Ec_pot == [0.004508186497043249, 9.864271368015835, 0.1312204454871218]
end
