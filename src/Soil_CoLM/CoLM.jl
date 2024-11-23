include("GW_UpdateRecharge.jl")
export GW_UpdateRecharge!, GW_UpdateDrainage!, GW_Correctθ!


## 参考CoLM的方案，考虑地下水的补给与排泄，只考虑地下水垂直方向与土壤水的交互作用。
# 从下至上，不受地下水影响的第一层
module CoLM

end
