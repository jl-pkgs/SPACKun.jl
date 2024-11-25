## 参考CoLM的方案，考虑地下水的补给与排泄，只考虑地下水垂直方向与土壤水的交互作用。
# 从下至上，不受地下水影响的第一层
include("Update_zwt_theta.jl")
include("SM_discharge.jl")
include("sw_balance_CoLM.jl")


export Update_zwt_theta!, GW_UpdateDrainage!, GW_Correctθ!
export sw_balance_CoLM
