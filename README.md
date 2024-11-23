# SPAC.jl

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://CUG-hydro.github.io/SPACKun.jl/stable) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://CUG-hydro.github.io/SPACKun.jl/dev)
[![CI](https://github.com/CUG-hydro/SPACKun.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/CUG-hydro/SPACKun.jl/actions/workflows/CI.yml)
[![Codecov](https://codecov.io/gh/CUG-hydro/SPACKun.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/CUG-hydro/SPACKun.jl/tree/master)

## Installation

```
using Pkg
Pkg.add(url="https://github.com/CUG-hydro/SPACKun.jl")
```

## TODO

- [ ] 测试三种土壤水运动方案
  + [x] SITH
  + [ ] BEPS
  + [ ] Bonan2019


## ChangeLog

1. 修复`swb_case4`逻辑漏洞

   ```julia
   if wa3 > θ_sat
     ff3 = max((wa3 - θ_sat) * d3, 0)
     wa3 = θ_sat # fix by kongdd, 2024-11-23
     wa3_unsat = θ_sat
   end
   ```

2. 引入`SM_discharge!` 统一土壤水排泄。

   **原版：**

   ```julia
   # Drainage from unsaturated zone, #1
   f1 = soil_drainage(wa1_unsat, θ_sat, Ksat, 0.048, 4.8)
   wa1 = (wa1 * d1 - f1 - Es - Tr1) / d1
   wa1 = clamp(wa1, 0, 1)

   # Drainage from unsaturated zone, #2
   f2 = soil_drainage(wa2_unsat, θ_sat, Ksat, 0.012, 1.2)
   wa2 = (wa2 * d2 + f1 - f2 - Tr2) / d2
   wa2 = clamp(wa2, 0, 1)  # > wwp 0

   if wa2 > θ_sat
      ff2 = max((wa2 - θ_sat) * d2, 0)  # extra water from upper layer
      wa2 = θ_sat
   else
      ff2 = 0
   end

   # Drainage from unsaturated zone, #3
   f3 = soil_drainage(wa3_unsat, θ_sat, Ksat, 0.012, 1.2)
   wa3 = (wa3 * d3 + f2 + ff2 - f3 - Tr3) / d3
   wa3 = clamp(wa3, 0, 1)

   if wa3 > θ_sat
      ff3 = max((wa3 - θ_sat) * d3, 0)
      wa3 = θ_sat
      wa3_unsat = θ_sat
   else
      ff3 = 0
      wa3_unsat = wa3
   end
   exceed = f3 + ff3
   ```

   **新版：**

   ```julia
   sink = [Tr1 + Es, Tr2, Tr3]
   θ_unsat = [wa1_unsat, wa2_unsat, wa3_unsat]
   exceed = SM_discharge!(soil, θ_unsat, sink, soilpar)
   wa1, wa2, wa3_unsat = θ_unsat
   wa3 = wa3_unsat
   ```

## 存在问题

- [ ] 地下水水位更新过程汇总，$\theta$不统一，一会是原始的$\theta_{unsat}$，一会考虑补给之后的$\theta$
