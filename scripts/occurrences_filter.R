
source('scripts/occurrences.R')
library(terra)


# occurrencias
str(occ_data)
colnames(occ_data)[colnames(occ_data)=='decimallatitude'] <- 'y'
colnames(occ_data)[colnames(occ_data)=='decimallongitude'] <- 'x'


# correct some common mistakes from widespread species
occ_data$species[occ_data$species=='Scleria secans' & occ_data$x > 0] <- 'Scleria boivinii' # Scleria secans es Americana
occ_data$species[occ_data$species=='Scleria hirtella' & occ_data$x > -30] <- 'Scleria distans' # Scleria hirtella es Americana
# check individual ranges
occ_data <- occ_data %>% filter(!(gbifid %in% c( occ_filtered$gbifid[occ_filtered$species=='Scleria gaertneri' & occ_filtered$x > 100],
                                                 occ_filtered$gbifid[occ_filtered$species=='Scleria distans' & occ_filtered$x > 100],
                                                 occ_filtered$gbifid[occ_filtered$species=='Scleria terrestris' & occ_filtered$x < 0],
                                                 occ_filtered$gbifid[occ_filtered$species=='Scleria racemosa' & occ_filtered$x < -30],
                                                 occ_filtered$gbifid[occ_filtered$species=='Scleria brownii' & occ_filtered$x > -100 & occ_filtered$x < 0],
                                                 occ_filtered$gbifid[occ_filtered$species=='Scleria dregeana' & occ_filtered$x < 0],
                                                 occ_filtered$gbifid[occ_filtered$species=='Scleria levis' & occ_filtered$x < 0],
                                                 occ_filtered$gbifid[occ_filtered$species=='Scleria distans' & occ_filtered$x > 100],
                                                 occ_filtered$gbifid[occ_filtered$species=='Scleria novae-hollandiae' & occ_filtered$x < 75],
                                                 occ_filtered$gbifid[occ_filtered$species=='Scleria parvula' & occ_filtered$x < 60],
                                                 occ_filtered$gbifid[occ_filtered$species=='Scleria pauciflora' & occ_filtered$y < 0],
                                                 occ_filtered$gbifid[occ_filtered$species=='Scleria rugosa' & occ_filtered$x < 0]
                                                 ) ) )


# will use same variables as in Larridon et al 2021
MAT <- rast('C:/Users/user/Desktop/wc2.1_5m_bio/wc2.1_5m_bio_1.tif') # mean annual temperature
AP <- rast('C:/Users/user/Desktop/wc2.1_5m_bio/wc2.1_5m_bio_12.tif') # annual precipitation
# 5 min is 9.28 km at the equator

# raw plot 
p_occ_data <- vect(occ_data, geom=c('x','y'), 'epsg:4326')
plot(MAT, main='raw'); points(p_occ_data)

# extract variables
occ_data$cell <- terra::extract(MAT, p_occ_data, method='simple', cells=T)[,3]
occ_data$MAT <- terra::extract(MAT, p_occ_data, method='simple')[,2]
occ_data$AP <- terra::extract(AP, p_occ_data, method='simple')[,2]


# filter
occ_filtered <- list()
for (sp in unique(occ_data$species)) {
  
  # select a species
  ss_df <- occ_data %>% filter(species == sp)
  
  # select one observation per cell
  ss_df <- ss_df %>% group_by(cell) %>% slice_sample(n=1)
  
  # climatic filter
  ss_df <- ss_df %>% filter(MAT > boxplot.stats(ss_df$MAT, coef=1.5)$stats[1] & MAT < boxplot.stats(ss_df$MAT, coef=1.5)$stats[5] &
                            AP > boxplot.stats(ss_df$AP, coef=1.5)$stats[1] & AP < boxplot.stats(ss_df$AP, coef=1.5)$stats[5] )
  
  # save
  occ_filtered[[sp]] <- ss_df
}

# list to dataframe
occ_filtered <- do.call(rbind, occ_filtered)
write.table(occ_filtered, 'results/occ_filtered.txt')


# filtered plot 
p_occ_filtered <- vect(occ_filtered, geom=c('x','y'), 'epsg:4326')
plot(MAT, main='filtered'); points(p_occ_filtered)

# check individual species
# a.spp <- unique(occ_filtered$species[order(occ_filtered$species)])
# for (sp in a.spp) {
#   ss_df <- occ_filtered %>% filter(species == sp)
#   p_ss_df <- vect(ss_df, geom=c('x','y'), 'epsg:4326')
#   plot(MAT, main=sp); points(p_ss_df)
#   Sys.sleep(2)
# }

  
  
  
  
  
  

