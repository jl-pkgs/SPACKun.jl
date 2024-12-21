pacman::p_load(
  Ipaper, data.table, dplyr, lubridate, hydroTools
)

## 栾城站
# MJ m-2 d-1
df_obs = read_ufile("Project_SM&VPD/data/Flux/CRO_栾城_Day_Flux_200710-201809_WithGPP.csv")
d_obs = df_obs$data[, 
  .(date = ymd(date), 
    Rn_obs = MJ_2W(Rn), H = MJ_2W(Hs), 
    Prcp_obs = Prcp_M, Prcp_A = as.numeric(Prcp_A), 
    ET_obs = ET, GPP_obs = GPP)]

d_forcing = fread("data/calib/raw/dat_栾城_ERA5L.csv") %>% 
  select(-LULC, -ET_obs, -GPP_obs) 
dat = merge(d_forcing, d_obs)
fwrite(dat, "data/dat_栾城_ERA5L_V2.csv")


## 禹城
df_obs = read_ufile("./Project_SM&VPD/data/Flux/CRO_禹城_Day_FluxMet_2003-2010.csv")
d_obs = df_obs$data[, 
  .(date = make_date(year, month, day),
    Rn_obs = Rn, 
    # LE, Hs, Ta, 
    Prcp_obs = Prcp, 
    ET_obs = LE * 0.0864/ cal_lambda(Ta), 
    GPP_obs = GPP
    # RH, U2, ea
  )
]

d_obs[, sum(ET_obs), .(year(date))]

d_forcing <- fread("data/calib/raw/dat_禹城_ERA5L.csv") %>%
  select(-LULC, -ET_obs, -GPP_obs)
dat = merge(d_forcing, d_obs)
fwrite(dat, "data/calib/dat_ERA5L_禹城.csv")


## 固城
df_obs = read_ufile("./Project_SM&VPD/data/Flux/CRO_固城_Day_FluxMet_2020_2022.csv")
d_obs = df_obs$data[, 
  .(date = make_date(year, month, day),
    Rn_obs = Rn, 
    # LE, Hs, Ta, 
    Prcp_obs = Prcp, 
    ET_obs = LE * 0.0864/ cal_lambda(Ta), 
    GPP_obs = GPP
    # RH, U2, ea
  )
]

d_forcing <- fread("data/calib/raw/dat_固城_ERA5L.csv") %>%
  select(-LULC, -ET_obs, -GPP_obs)
dat = merge(d_forcing, d_obs)
fwrite(dat, "data/dat_固城_ERA5L_V2.csv")


## 分析导致不好的原因
library(ppcor)

## 这里不清楚，到底出了什么问题
df = fread("data/calib/dat_ERA5L_禹城.csv")
dat = df[, .(ET_obs, 
  Rs_toa = cal_Rsi_toa(36.95, yday(df$date)), 
  Rn, Prcp, VOD, LAI, VPD, U2, Rn_obs, Prcp_obs)] %>% na.omit()

pcor(dat, method = "pearson")
cor(dat)
