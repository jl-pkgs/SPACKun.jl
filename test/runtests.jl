using Test, SPAC

include("test_栾城_2010.jl")

@testset "potentialET" begin
  Rn = 100.0
  G = 5.0
  LAI = 3.0
  Ta = 20.0
  Pa = 100.0
  pEc, pEs = potentialET(Rn, G, LAI, Ta, Pa)

  @test pEc ≈ 2.5389604966643824
  @test pEs ≈ 0.3507115521629298
end

@testset "runoff_up" begin
  Pnet = 20.0 # mm
  zgw = 1000.0 # mm
  ZM = [50, 1450, 3500]
  wa = [0.3, 0.3, 0.3]
  soilpar = get_soilpar(2)
  srf, IWS, Vmax = runoff_up(Pnet, zgw, ZM, wa, soilpar)
  @test Vmax ≈ 129
  @test IWS ≈ 20
end
# mxcall(:runoff_up, 3, Pnet, zgw, ZM, wa, collect(soilpar))

@testset "pTr_partition" begin
  pEc = 20.0
  wa1, wa2, wa3 = 0.3, 0.2, 0.1
  fwet = 0.5
  soilpar = get_soilpar(2) |> collect
  pftpar = get_pftpar(22) |> collect
  r = pTr_partition(pEc, wa1, wa2, wa3, soilpar, pftpar, fwet, ZM)
  @test r == (0.004508186497043249, 9.864271368015835, 0.1312204454871218)
end
# mxcall(:pTr_partition, 3, pEc, wa1, wa2, wa3, soilpar, pftpar, fwet, ZM)
