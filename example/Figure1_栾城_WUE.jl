using NaNStatistics

begin
  gr(framestyle=:box)
  p_et = plot(ET_obs, label="ET_obs", color=:black)
  plot!(p_et, ET_Kun, label="ET_kun")
  plot!(p_et, ET, label="ET_kong")
end

function plot_gpp()
  plot(GPP_obs, label="GPP_obs")
  plot!(ET_obs * 2, label="ET_obs x 2")
end

## 检查驱动数据
begin
  function plot_var(var=:Rn)
    plot(d[:, var], label=string(var))
  end
  plot(
    plot_var(:Rn),
    plot_var(:Prcp),
    plot_var(:VPD),
    plot_var(:U2),
    plot_var(:LAI),
    plot_var(:VOD),
    plot_gpp(),
    p_et,
    size=(1000, 700)
  )
end

## 去掉有降水的日期
## WUE
function plot_wue(ET; label)
  trs = quantile(GPP_obs, 0.5) # 生长季阈值

  plot(GPP_obs * 0.6, label="GPP_obs")

  wue = GPP_obs ./ ET
  wue[wue.>15] .= NaN      # 去除异常值

  con_met = @. (P_obs > 1 && Rn_obs > 100)
  # plot!(wue; label="")
  wue[(GPP_obs .< trs) .|| con_met] .= NaN # 去除非生长季

  t1 = 70:170
  t2 = 190:285
  y1 = nanmean(wue[t1])
  y2 = nanmean(wue[t2])
  plot!(wue; label)
  plot!(t1, t1 * 0 .+ y1; label="") # , label="wheat"
  plot!(t2, t2 * 0 .+ y2; label="")
end

begin
  plot(
    plot_wue(ET_obs, label="WUE_obs"),
    plot_wue(ET_Kun, label="WUE_kun"),
    plot_wue(ET, label="WUE_kong"),
  )
end

## 辐射的原因，导致WUE的跳动
