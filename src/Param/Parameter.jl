"""
  get_pftpar(LC::Int)

# Output:
- pftpar
  + 1  `β`   : the coefficient for calculation of interceptions
  + 2  `D50` : the depth above which 50% of the root mas is located, mm
  + 3  `c`   : the shape parameters of logistic dose-response root distribution model
  + 4  `Hc`  : root depth (m)

Standard numbering rule for different PFTs based on MCD12, IGBP
- `0` : 'water'
- `1` : 'evergreen needleleaf forest'
- `2` : 'evergreen broadleaf forest'
- `3` : 'deciduous needleleaf forest'
- `4` : 'deciduous broadleaf forest'
- `5` : 'mixed forests'
- `6` : 'closed shrubland'
- `7` : 'open shrublands'
- `8` : 'woody savannas'
- `9` : 'savannas'
- `10`: 'grasslands'
- `11`: 'permanent wetlands'
- `12`: 'croplands'
- `13`: 'urban and built-up'
- `14`: 'cropland/natural vegetation mosaic'
- `15`: 'snow and ice'
- `16`: 'barren or sparsely vegetated'
"""
function get_pftpar(LC::Int)
  # The PFTpar look up table
  PLUT = [
    # :β, :D50, :D95, :Hc
    0.06 1771 3998 20;   # 1  Evergreen Needleleaf Forest   ------ ENF
    0.02 2187 4316 18;   # 2  Evergreen Broadleaf  Forest   ------ EBF
    0.06 1668 3936 18;   # 3  Deciduous Needleleaf Forest   ------ DNF
    0.06 1401 3888 18;   # 4  Deciduous Broadleaf  Forest   ------ DBF
    0.04 1748 3866 15;   # 5  Mixed                Forest   ------ MF
    0.01 701 2663 1.5;   # 6  Closed Shrubland              ------ CSH
    0.01 753 2872 1.5;   # 7  Open   Shrubland              ------ OSH
    0.01 765 2253 1.2;   # 8  Woody  Savannas               ------ WSA
    0.01 559 1877 1.0;   # 9  Savannas                      ------ SAV
    0.01 598 1823 1.0;   # 10 Grassland                     ------ GRA
    0.01 896 2203 1.5;   # 11 Cropland                      ------ CRO
    0.01 752 2809 1.5    # 12 Wetland                       ------ WET
  ]
  # Adjust the PFT code to match the IGBP-LC
  # IGBP <-----> PFTpar table
  pftpar = PLUT[LC, :]
  @LArray pftpar (:β, :D50, :D95, :Hc)
end

"""
# INPUT
Soil Type index

# OUTPUT
  soilpar ::
     1. `Ksat`  : hydraulic conductivity, mm day-1
     2. `ψ_sat` : water potential at saturation condition
     3. `θ_sat` : water content at saturation condition
     4. `b`     : an empirical parameter
     5. `θ_fc`  : soil moisture at field capacity
     6. `θ_c`   : soil moisture at critical value
     7. `θ_wp`  : soil moisture at wilting point

 ----------- 3-D rasters that contain the seven soil characteristic parameters ---------------
  1     2        3       4      5      6      7
  Ksat  ψ_sat    θ_sat   b      θ_fc   θ_c    θ_wp   Soil Type
"""
function get_soilpar(SC)
  store = [
    # Ksat, ψ_sat, θ_sat, b,  θ_fc, θ_c,  θ_wp
    1.3841 -0.0232 0.373 3.39 0.151 0.109 0.035;  # 1 Sand
    0.8229 -0.0175 0.386 3.86 0.189 0.142 0.052;  # 2 Loamy sand
    0.5353 -0.0316 0.416 4.50 0.265 0.208 0.087;  # 3 Sandy loam
    0.4086 -0.0655 0.435 5.77 0.331 0.274 0.139;  # 4 Loam
    0.4331 -0.0562 0.455 5.32 0.371 0.317 0.149;  # 5 Silt
    0.4427 -0.0471 0.468 4.98 0.382 0.320 0.156;  # 6 Silty loam
    0.4991 -0.0310 0.416 7.20 0.314 0.270 0.157;  # 7 Sandy clay loam
    0.3552 -0.0599 0.449 8.32 0.387 0.339 0.212;  # 8 Clay loam
    0.3848 -0.0414 0.476 8.32 0.344 0.319 0.203;  # 9 Silty clay loam
    0.6157 -0.0269 0.423 9.59 0.349 0.322 0.207;  # 10 Sandy clay
    0.2967 -0.0453 0.481 10.4 0.390 0.384 0.274;  # 11 Silty clay
    0.2580 -0.0531 0.461 12.1 0.417 0.390 0.282   # 12 Clay
  ]
  store[:, 1] .= 1000 .* store[:, 1]  # m --> mm
  soilpar = store[SC, :]
  return @LArray(soilpar, (:Ksat, :ψ_sat, :θ_sat, :b, :θ_fc, :θ_c, :θ_wp))
end


# if Typenum in [0, 13, 15, 16]
#   pftpar = PLUT[10, :]  # set to Grassland
# elseif Typenum == 11
#   pftpar = PLUT[12, :]  # set to Wetland 
# elseif Typenum in [12, 14]
#   pftpar = PLUT[11, :]  # set to Cropland
# else
#   pftpar = PLUT[Typenum, :]
# end
# function query_LC(LC::Int)
#   keys = Dict(
#     11 => 13,
#     22 => 12,
#     33 => 10,
#     40 => 4,
#     41 => 1,
#     42 => 2,
#     43 => 3,
#     44 => 4,
#     45 => 5,
#     55 => 10,
#     66 => 10,
#     77 => 11,
#     0 => 10
#   )
#   haskey(keys, LC) ? keys[LC] : 10
# end
