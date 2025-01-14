coords<- read_csv("Data/coord.csv")
coords2=read.csv("./Data/coord_lluvi.csv")


coords_mod=coords %>% mutate(SampleID=paste0(pol, Sitio, Transecto, "D")) %>% dplyr::select(
  SampleID, Longitude, Latitude)
coords_mod2=coords2 %>% mutate(SampleID=paste0(pol, Sitio, Transecto, "D")) %>% dplyr::select(
  SampleID, Longitude, Latitude)

coordi = rbind(coords_mod, coords_mod2)

library(elevatr)
examp <- sf::st_as_sf(coordi, coords = c("Longitude", "Latitude"), crs = crs_dd)
df_elev_epqs <- get_elev_point(examp, prj = crs_dd, src = "epqs")
df_elev_epqs
writexl::write_xlsx(df_elev_epqs, "Data/elev.xlsx")