function query_LC(LC::Int)
  keys = Dict(
    11 => 13,
    22 => 12,
    33 => 10,
    40 => 4,
    41 => 1,
    42 => 2,
    43 => 3,
    44 => 4,
    45 => 5,
    55 => 10,
    66 => 10,
    77 => 11,
    0 => 10
  )
  haskey(keys, LC) ? keys[LC] : 10
end

"""
    get_pftpar(LC::Int)

# Output:
- pftpar
  + 1  `beta`: the coefficient for calculation of interceptions
  + 2  `D50` : the depth above which 50% of the root mas is located, mm
  + 3  `c`   : the shape parameters of logistic dose-response root distribution model
  + 4  `Zr`  : root depth (m)

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
    # :β, :D50, :D95, :Zr
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
  Typenum = query_LC(LC)

  # IGBP <-----> PFTpar table
  if Typenum in [0, 13, 15, 16]
    pftpar = PLUT[10, :]  # set to Grassland
  elseif Typenum == 11
    pftpar = PLUT[12, :]  # set to Wetland 
  elseif Typenum in [12, 14]
    pftpar = PLUT[11, :]  # set to Cropland
  else
    pftpar = PLUT[Typenum, :]
  end
  @LArray pftpar (:β, :D50, :D95, :Zr)
end
