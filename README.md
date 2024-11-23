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

### 地下水更新的计算方案

现在的版本是分开计算的。

- `ET_SM`: `SM_discharge!`，简化版`Richards`方程，非饱和流方程

- `ET_GW`: 直接从`Δw`中扣除
  `Δw = exceed + vw2 - sum(Ec_gw) - GW_Rsb(zwt)`

## ChangeLog

1. 修复`swb_case4`逻辑漏洞

   ```julia
   if wa3 > θ_sat
     ff3 = max((wa3 - θ_sat) * d3, 0)
     wa3 = θ_sat # fix by kongdd, 2024-11-23
     wa3_unsat = θ_sat
   end
   ```


## 存在问题

- [ ] 地下水水位更新过程汇总，$\theta$不统一，一会是原始的$\theta_{unsat}$，一会考虑补给之后的$\theta$
