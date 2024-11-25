using SPAC, Test

begin
  N = 100
  Δz = fill(0.02, N)
  θ = fill(0.3, N)
  soil = Soil(; Δz, θ)
  z₊ₕ = soil.z₊ₕ

  wa = 4000.0 # [mm]
  zwt = 0.5
  Δt = 60 # [s]
  recharge = 1 / 3600 * Δt # [mm s-1], namely [1 mm h-1]
end

@testset "find_jwt" begin
  @test find_jwt(z₊ₕ, -0.01) == 0
  @test find_jwt(z₊ₕ, 0.01) == 1
  @test find_jwt(z₊ₕ, 0.02) == 1
  @test find_jwt(z₊ₕ, 0.03) == 2
  @test find_jwt(z₊ₕ, 100.0) == 101
end


@testset "Update_zwt_theta! with θ" begin
  θ = fill(0.3, N)
  Update_zwt_theta!(soil, 0.5, wa, recharge; θ) == # # 1 mm h-1
  (zwt=0.49916666666666665, wa=4000.016666666667, uex=0.0)
  @test θ[25] == 0.30083333333333334

  θ = fill(0.3, N)
  @test Update_zwt_theta!(soil, 0.5, wa, -recharge; θ) == # # 1 mm h-1
        (zwt=0.5008333333333334, wa=3999.983333333333, uex=0.0)
  @test θ[26] == 0.29916666666666675
end


@testset "Update_zwt_theta!" begin
  @test Update_zwt_theta!(soil, 0.5, wa, recharge) == # 1 mm h-1
        (zwt=0.49916666666666665, wa=4000.016666666667, uex=0.0)

  @test Update_zwt_theta!(soil, 0.5, wa, 1000 / 60) == # 1000 mm h-1
        (zwt=0.0, wa=4016.6666666666665, uex=6.666666666666663)

  @test Update_zwt_theta!(soil, 0.5, wa, -recharge) == # -1 mm h-1
        (zwt=0.5008333333333334, wa=3999.983333333333, uex=0.0)

  @test Update_zwt_theta!(soil, 0.5, wa, -10000 / 60) == # -10,000 mm h-1
        (zwt=8.833333333333291, wa=3833.3333333333335, uex=0.0)
end


@testset "GW_UpdateDrainage!" begin
  θ = fill(0.3, N)
  jwt = find_jwt(z₊ₕ, 0.5)
  jwt2 = find_jwt(z₊ₕ, 1.0)
  r = GW_UpdateDrainage!(soil, 0.5, 4000.0, 600 / 60; θ)

  @test r == (zwt=1.0, wa=3990.0)
  @test all(θ[jwt:jwt2] .< 0.3)

  @test GW_UpdateDrainage!(soil, 2.5, 4000.0, 1 / 60; θ) ==
        (zwt=2.5008333333333335, wa=3999.983333333333)
end


@testset "GW_Correctθ!" begin
  # 亏损
  θ = fill(0.3, N)
  @test GW_Correctθ!(soil, θ, 0.5, 4000.0, 60, 600 / 3600) == (zwt=0.5, wa = 4000.0, uex=0.0, drainage=0.16666666666666666)

  θ[2:3] .= -0.1
  @test GW_Correctθ!(soil, θ, 0.5, 4000.0, 60, 600 / 3600) == (zwt=0.5, wa=4000.0, uex=0.0, drainage=0.09999999999999999)

  θ[2:3] .= -10.0 # [m3 m-3]
  @test GW_Correctθ!(soil, θ, 0.5, 4000.0, 60, 600 / 3600) == (zwt = 20, wa=3610.0, uex=0.0, drainage=0.0)
  
  # 超饱和
  θ = fill(0.3, N)
  θ[2:3] .= 0.6
  @test GW_Correctθ!(soil, θ, 0.5, 4000.0, 60, 600 / 3600) == (zwt = 0.5, wa=4000.0, uex=5.999999999999998, drainage=0.16666666666666666)

  θ = fill(0.3, N)
  θ[100] = 12.0
  r = GW_Correctθ!(soil, θ, 0.5, 4000.0, 60, 600 / 3600) 
  @test all(θ .== 0.4)
  @test r == (zwt=0.5, wa = 4000.0, uex=33.999999999999986, drainage=0.16666666666666666)
end
