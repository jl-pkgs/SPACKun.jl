## 参考CoLM的方案，考虑地下水的补给与排泄，只考虑地下水垂直方向与土壤水的交互作用。
# 从下至上，不受地下水影响的第一层
include("GW_Update.jl")
include("SM_discharge.jl")
include("SM_balance.jl")

export GW_Update_zwtθ!, GW_UpdateDrainage!, GW_Correctθ!
