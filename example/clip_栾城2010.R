pacman::p_load(
  Ipaper, data.table, dplyr, lubridate, hydroTools
)

# MJ m-2 d-1
df_obs = read_ufile("Project_SM&VPD/data/Flux/CRO_栾城_Day_Flux_200710-201809_WithGPP.csv")
d_obs = df_obs$data[, 
  .(date = ymd(date), 
    Rn_obs = MJ_2W(Rn), H = MJ_2W(Hs), 
    P_obs = Prcp_M, Prcp_A = as.numeric(Prcp_A), 
    ET_obs = ET, GPP_obs = GPP)]

d_forcing = fread("data/calib/dat_栾城_ERA5L.csv") %>% 
  select(-LULC, -ET_obs, -GPP_obs) 

dat = merge(d_forcing, d_obs)
fwrite(dat, "data/dat_栾城_ERA5L_V2.csv")
# fwrite(d_obs, "data/dat_栾城_Flux_2010.csv")
# u10 = fread("data/temp/EAR5L_U10_01deg_1979-2023_luancheng.csv")
# vpd = fread("data/temp/EAR5L_VPD_01deg_1979-2023_luancheng.csv")
# d_met = merge(rename(u10, U10 = value), rename(vpd, VPD = value)) %>% 
#   mutate(U2 = cal_U2(U10, 10), time = ymd(time)) %>% 
#   rename(date = time)
# d_met = merge(d_met[, .(date, VPD, U2)], d_obs)
# d_era5 = fread("data/dat_栾城_ERA5L_2010.csv")[, 1:8] %>% select(-LULC)
# dat = merge(d_era5, d_met)
fwrite(dat, "data/dat_栾城_ERA5L_V2.csv")
