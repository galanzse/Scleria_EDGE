

# ECOREGIONS WHERE SCLERIA IS PRESENT FOLLOWING OLSON ET AL. 2001


library(tidyverse)
library(terra)

load("C:/Users/user/OneDrive/PUBLICACIONES/SCLERIA/Scleria_EDGE/results/data_final.RData") # data


# ecoregions
ecoregions <- vect('C:/Users/user/Desktop/ecoregions/wwf_terr_ecos.shp')
df_ecoregions <- data_final[['occurrences']][,c('scientific_name','x','y')]


# biome_rast <- rasterize(x=wwf_biomes, y= disagg(RAST1_eqearth, 10), field='BIOME')


df_ecoregions <- data_final[['occurrences']][,c('scientific_name','x','y')] %>%
  cbind(
    terra::extract(x=ecoregions, y=data_final[['occurrences']][,c('x','y')]) %>%
      dplyr::select(ECO_NAME, REALM, BIOME, ECO_NUM, GBL_STAT, G200_REGIO)
  )


# correct levels
df_ecoregions$REALM <- as.factor(df_ecoregions$REALM)
levels(df_ecoregions$REALM) <- c("Australasia", "Afrotropics", "IndoMalay", "Nearctic", "Neotropics", "Palearctic")

df_ecoregions$BIOME <- as.factor(df_ecoregions$BIOME)
levels(df_ecoregions$BIOME) <- c("Tropical & Subtropical Moist Broadleaf Forests","Tropical & Subtropical Dry Broadleaf Forests","Tropical & Subtropical Coniferous Forests","Temperate Broadleaf & Mixed Forests","Temperate Conifer Forests","Tropical & Subtropical Grasslands, Savannas & Shrublands","Temperate Grasslands, Savannas & Shrublands","Flooded Grasslands & Savannas","Montane Grasslands & Shrublands","Deserts & Xeric Shrublands","Mangroves",NA)


# save
df_ecoregions <- merge(data_final[['traits']][,c('scientific_name', 'subgenus', 'section')], df_ecoregions, by='scientific_name')
write.table(df_ecoregions, 'results/df_ecoregions.txt')

