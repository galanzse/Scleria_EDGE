

library(tidyverse); library(readxl)
library(rnaturalearth)
library(terra)
library(ggplot2)
# library(viridis)
library(ggpubr)



# import datasets and results ####

# world map
scl_wrld_map <- ne_countries(scale=10, type="countries", continent=NULL,
                     country=NULL, geounit=NULL, sovereignty=NULL, returnclass="sf") %>%
  terra::vect() %>% terra::crop(ext(-180,180,-50,70)) %>% terra::project('+proj=eqearth')

# ecoregions
scl_ecor_map <- terra::vect('C:/Users/user/Desktop/ecoregions/wwf_terr_ecos.shp') %>%
  terra::crop(ext(-180,180,-50,70)) %>% terra::project('+proj=eqearth')


# data import
load('results/data_final.RData')


# lists
EDGE2_list <- read_excel("results/results.xlsx", sheet="EDGE2_list")
EcoDGE2_list <- read_excel("results/results.xlsx", sheet="EcoDGE2_list")


# occurrences
scleria_occ <- data_final[['occurrences']] %>%
  dplyr::select(scientific_name, x, y) %>%
  merge(data_final[['taxa']][,c('scientific_name','section','subgenus')])
scleria_occ$subgenus <- as.factor(scleria_occ$subgenus)
scleria_occ$section <- as.factor(scleria_occ$section)



# fix names
EDGE2_list$scientific_name[!(EDGE2_list$scientific_name %in% scleria_occ$scientific_name)]
EDGE2_list$scientific_name[which(EDGE2_list$scientific_name=='Scleria rubrostriata A.C.AraÃºjo & N.A.Brummitt')] <- 'Scleria rubrostriata A.C.Araújo & N.A.Brummitt'
EDGE2_list$scientific_name[which(EDGE2_list$scientific_name=='Scleria adpressohirta (KÃ¼k.) E.A.Rob.')] <- 'Scleria adpressohirta (Kük.) E.A.Rob.'

EcoDGE2_list$scientific_name[!(EcoDGE2_list$scientific_name %in% scleria_occ$scientific_name)]
EcoDGE2_list$scientific_name[which(EcoDGE2_list$scientific_name=='Scleria rubrostriata A.C.AraÃºjo & N.A.Brummitt')] <- 'Scleria rubrostriata A.C.Araújo & N.A.Brummitt'
EDGE2_list$scientific_name[which(EcoDGE2_list$scientific_name=='Scleria adpressohirta (KÃ¼k.) E.A.Rob.')] <- 'Scleria adpressohirta (Kük.) E.A.Rob.'



# points file
p_scleria_occ <- vect(scleria_occ, geom=c('x','y'), 'epsg:4326') %>% terra::project('+proj=eqearth') # points file

# plot
plot(scl_wrld_map)
points(p_scleria_occ, col='green')
lines(scl_wrld_map)



# species lists vectors
E2_list <- EDGE2_list$scientific_name[!is.na(EDGE2_list$list)]
E2_main <- EDGE2_list %>% filter(list=='main') %>% dplyr::select(scientific_name) %>% deframe()
E2_25 <- EDGE2_list[order(EDGE2_list$EDGE2, decreasing=T)[1:25],] %>% dplyr::select(scientific_name) %>% deframe()

Eco2_list <- EcoDGE2_list$scientific_name[!is.na(EcoDGE2_list$list)]
Eco2_main <- EcoDGE2_list %>% filter(list=='main') %>% dplyr::select(scientific_name) %>% deframe()
Eco2_25 <- EcoDGE2_list[order(EcoDGE2_list$EcoDGE, decreasing=T)[1:25],] %>% dplyr::select(scientific_name) %>% deframe()



# countries ####

# scl_countries <- scl_wrld_map %>% terra::extract(p_scleria_occ) %>%
#   dplyr::select(sovereignt, admin, geounit, subunit, name, name_long, pop_year, gdp_md, gdp_year, economy)
# scl_countries <- cbind(scleria_occ[,c('scientific_name', 'x', 'y', 'section', 'subgenus')],
#                        scl_countries) %>% select(scientific_name,name) %>% na.omit() %>% unique()
# colnames(scl_countries) <- c('scientific_name','country')

write.table(scl_countries, 'results/scl_countries.txt')


EDGE_countries <- matrix(nrow=length(unique(scl_countries$country)), ncol=10) %>% as.data.frame()
colnames(EDGE_countries) <- c('country','richness','EDGE2_list','EDGE2_main', 'EDGE2_top25','EDGE2_sum','EcoDGE2_list','EcoDGE2_main', 'EcoDGE2_top25','EcoDGE2_sum')
EDGE_countries$country <- unique(scl_countries$country)

EDGE_countries$country[EDGE_countries$country=='France'] <- 'French Guiana'
EDGE_countries <- EDGE_countries %>% filter(country!='Canada')

scl_wrld_map$richness <- NA
scl_wrld_map$EDGE2_sum <- NA
scl_wrld_map$EcoDGE2_sum <- NA
scl_wrld_map$EDGE2_main <- NA
scl_wrld_map$EcoDGE2_main <- NA
scl_wrld_map$EDGE2_list <- NA
scl_wrld_map$EcoDGE2_list <- NA
scl_wrld_map$EDGE2_top25 <- NA
scl_wrld_map$top25_EcoDGE2 <- NA

for (s in 1:nrow(EDGE_countries)) {
  
  # species in a given country
  sp1 <- scl_countries$scientific_name[scl_countries$country==EDGE_countries$country[s]] %>% na.omit() %>%
    unique() %>% as.vector()
  
  # richness
  EDGE_countries$richness[s] <- length(sp1)

  # lists
  EDGE_countries$EDGE2_list[s] <- length(E2_list[E2_list %in% sp1])
  EDGE_countries$EDGE2_main[s] <- length(E2_main[E2_main %in% sp1])
  EDGE_countries$EDGE2_top25[s] <- length(E2_25[E2_25 %in% sp1])
  EDGE_countries$EDGE2_sum[s] <- EDGE2_list %>% filter(scientific_name %in% sp1) %>%
    dplyr::select(EDGE2) %>% colSums()

  # EcoDGE2
  EDGE_countries$EcoDGE2_list[s] <- length(Eco2_list[Eco2_list %in% sp1])
  EDGE_countries$EcoDGE2_main[s] <- length(Eco2_main[Eco2_main %in% sp1])
  EDGE_countries$EcoDGE2_top25[s] <- length(Eco2_25[Eco2_25 %in% sp1])
  EDGE_countries$EcoDGE2_sum[s] <- EcoDGE2_list %>% filter(scientific_name %in% sp1) %>%
    dplyr::select(EcoDGE) %>% colSums()

  # save in spatvector
  scl_wrld_map$richness[scl_wrld_map$name==EDGE_countries$country[s]] <- EDGE_countries$richness[s]
  scl_wrld_map$EDGE2_sum[scl_wrld_map$name==EDGE_countries$country[s]] <- EDGE_countries$EDGE2_sum[s]
  scl_wrld_map$EcoDGE2_sum[scl_wrld_map$name==EDGE_countries$country[s]] <- EDGE_countries$EcoDGE2_sum[s]
  scl_wrld_map$EDGE2_main[scl_wrld_map$name==EDGE_countries$country[s]] <- EDGE_countries$EDGE2_main[s]
  scl_wrld_map$EcoDGE2_main[scl_wrld_map$name==EDGE_countries$country[s]] <- EDGE_countries$EcoDGE2_main[s]
  scl_wrld_map$EDGE2_list[scl_wrld_map$name==EDGE_countries$country[s]] <- EDGE_countries$EDGE2_list[s]
  scl_wrld_map$EcoDGE2_list[scl_wrld_map$name==EDGE_countries$country[s]] <- EDGE_countries$EcoDGE2_list[s]
  scl_wrld_map$EDGE2_top25[scl_wrld_map$ECO_NAME==EDGE_countries$ECO_NAME[s]] <- EDGE_countries$EDGE2_top25[s]
  scl_wrld_map$top25_EcoDGE2[scl_wrld_map$ECO_NAME==EDGE_countries$ECO_NAME[s]] <- EDGE_countries$EcoDGE2_top25[s]
  
  
  print(paste('--- ', round(s/nrow(EDGE_countries)*100, 2), '% ---', sep=''))
  
}

# save
# write.table(EDGE_countries, 'results/EDGE_countries.txt', sep=',')
# save(scl_wrld_map, file="results/scl_wrld_map.RData")

# check
terra::plot(scl_wrld_map, y='EDGE2_sum', breaks=30, col=rev(grDevices::heat.colors(30)), border=NA)
lines(scl_wrld_map, col='grey40')



# ecoregions ####

# scl_ecoregions <- scl_ecor_map %>% terra::extract(p_scleria_occ) %>%
#   dplyr::select(REALM, BIOME, ECO_NAME)
# scl_ecoregions <- cbind(scleria_occ[,c('scientific_name', 'x', 'y', 'section', 'subgenus')],
#                         scl_ecoregions) %>% select(scientific_name, REALM, BIOME, ECO_NAME) %>%
#   na.omit() %>% unique()
# colnames(scl_ecoregions) <- c('scientific_name','REALM','BIOME','ECO_NAME')

write.table(scl_ecoregions, 'results/scl_ecoregions.txt')


EDGE_scl_ecor_map <- matrix(nrow=length(unique(scl_ecoregions$ECO_NAME)), ncol=10) %>% as.data.frame()
colnames(EDGE_scl_ecor_map) <- c('ECO_NAME','richness','EDGE2_list','EDGE2_main', 'EDGE2_top25','EDGE2_sum','EcoDGE2_list','EcoDGE2_main', 'EcoDGE2_top25','EcoDGE2_sum')
EDGE_scl_ecor_map$ECO_NAME <- unique(scl_ecoregions$ECO_NAME)

scl_ecor_map$richness <- NA
scl_ecor_map$EDGE2_sum <- NA
scl_ecor_map$EcoDGE2_sum <- NA
scl_ecor_map$EDGE2_main <- NA
scl_ecor_map$EcoDGE2_main <- NA
scl_ecor_map$EDGE2_list <- NA
scl_ecor_map$EcoDGE2_list <- NA
scl_ecor_map$EDGE2_top25 <- NA
scl_ecor_map$top25_EcoDGE2 <- NA

for (s in 1:nrow(EDGE_scl_ecor_map)) {
  
  # species in a given country
  sp1 <- scl_ecoregions$scientific_name[scl_ecoregions$ECO_NAME==EDGE_scl_ecor_map$ECO_NAME[s]] %>% na.omit() %>%
    unique() %>% as.vector()
  
  # richness
  EDGE_scl_ecor_map$richness[s] <- length(sp1)
  
  # lists
  EDGE_scl_ecor_map$EDGE2_list[s] <- length(E2_list[E2_list %in% sp1])
  EDGE_scl_ecor_map$EDGE2_main[s] <- length(E2_main[E2_main %in% sp1])
  EDGE_scl_ecor_map$EDGE2_top25[s] <- length(E2_25[E2_25 %in% sp1])
  EDGE_scl_ecor_map$EDGE2_sum[s] <- EDGE2_list %>% filter(scientific_name %in% sp1) %>%
    dplyr::select(EDGE2) %>% colSums()
  
  # EcoDGE2
  EDGE_scl_ecor_map$EcoDGE2_list[s] <- length(Eco2_list[Eco2_list %in% sp1])
  EDGE_scl_ecor_map$EcoDGE2_main[s] <- length(Eco2_main[Eco2_main %in% sp1])
  EDGE_scl_ecor_map$EcoDGE2_top25[s] <- length(Eco2_25[Eco2_25 %in% sp1])
  EDGE_scl_ecor_map$EcoDGE2_sum[s] <- EcoDGE2_list %>% filter(scientific_name %in% sp1) %>%
    dplyr::select(EcoDGE) %>% colSums()
  
  # save variable in SpatVector
  scl_ecor_map$richness[scl_ecor_map$ECO_NAME==EDGE_scl_ecor_map$ECO_NAME[s]] <- EDGE_scl_ecor_map$richness[s]
  scl_ecor_map$EDGE2_sum[scl_ecor_map$ECO_NAME==EDGE_scl_ecor_map$ECO_NAME[s]] <- EDGE_scl_ecor_map$EDGE2_sum[s]
  scl_ecor_map$EcoDGE2_sum[scl_ecor_map$ECO_NAME==EDGE_scl_ecor_map$ECO_NAME[s]] <- EDGE_scl_ecor_map$EcoDGE2_sum[s]
  scl_ecor_map$EDGE2_main[scl_ecor_map$ECO_NAME==EDGE_scl_ecor_map$ECO_NAME[s]] <- EDGE_scl_ecor_map$EDGE2_main[s]
  scl_ecor_map$EcoDGE2_main[scl_ecor_map$ECO_NAME==EDGE_scl_ecor_map$ECO_NAME[s]] <- EDGE_scl_ecor_map$EcoDGE2_main[s]
  scl_ecor_map$EDGE2_list[scl_ecor_map$ECO_NAME==EDGE_scl_ecor_map$ECO_NAME[s]] <- EDGE_scl_ecor_map$EDGE2_list[s]
  scl_ecor_map$EcoDGE2_list[scl_ecor_map$ECO_NAME==EDGE_scl_ecor_map$ECO_NAME[s]] <- EDGE_scl_ecor_map$EcoDGE2_list[s]
  scl_ecor_map$EDGE2_top25[scl_ecor_map$ECO_NAME==EDGE_scl_ecor_map$ECO_NAME[s]] <- EDGE_scl_ecor_map$EDGE2_top25[s]
  scl_ecor_map$top25_EcoDGE2[scl_ecor_map$ECO_NAME==EDGE_scl_ecor_map$ECO_NAME[s]] <- EDGE_scl_ecor_map$EcoDGE2_top25[s]
  
  
  print(paste('--- ', round(s/nrow(EDGE_scl_ecor_map)*100, 2), '% ---', sep=''))
  
}

# save
# write.table(EDGE_scl_ecor_map, 'results/EDGE_scl_ecor_map.txt', sep=',')
# save(scl_ecor_map, file="results/scl_ecor_map.RData")


par(mfrow=c(2,1))

ref_rast <- rast('C:/Users/user/Desktop/worldclim/wc2.1_2.5m_elev/wc2.1_2.5m_elev.tif') %>%
  project('+proj=eqearth')

ecor_EDGE2sum <- scl_ecor_map %>% rasterize(y=ref_rast, field='EDGE2_sum')
plot(ecor_EDGE2sum, col=rev(grDevices::heat.colors(30)), main='Cumulative EDGE2 x ecorregion')
lines(scl_wrld_map, col='grey40')

ecor_EDGE2sum <- scl_ecor_map %>% rasterize(y=ref_rast, field='EcoDGE2_sum')
plot(ecor_EDGE2sum, col=rev(grDevices::heat.colors(30)), main='Cumulative EcoDGE x ecorregion')
lines(scl_wrld_map, col='grey40')



# vulnerability

scl_stat <- scl_ecor_map[,c('ECO_NAME','GBL_STAT','G200_STAT')] %>% as.data.frame() %>% unique()
scl_stat <- merge(EDGE_scl_ecor_map, scl_stat, all.x=T)

ggplot(aes(x=as.factor(GBL_STAT), y=EDGE2_sum), data=scl_stat) +
  geom_boxplot()

ggplot(aes(x=as.factor(G200_STAT), y=EDGE2_sum), data=scl_stat) +
  geom_boxplot()


