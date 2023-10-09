

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

scl_ecor_map <- terra::vect('C:/Users/user/Desktop/ecoregions/wwf_terr_ecos.shp') %>%
  terra::crop(ext(-180,180,-50,70)) %>% terra::project('+proj=eqearth')


# data import
load('results/data_final.RData')


# lists
EDGE2_list <- read_excel("results/final_lists.xlsx", sheet="EDGE2_list")
EcoDGE2_list <- read_excel("results/final_lists.xlsx", sheet="EcoDGE2_list")


# occurrences
scleria_occ <- data_final[['occurrences']] %>%
  dplyr::select(scientific_name, x, y) %>%
  merge(data_final[['assessments']][,c('scientific_name','section','subgenus')])
scleria_occ$subgenus <- as.factor(str_to_title(scleria_occ$subgenus))
scleria_occ$section <- as.factor(str_to_title(scleria_occ$section))



# Fix names
EDGE2_list$scientific_name[!(EDGE2_list$scientific_name %in% scleria_occ$scientific_name)]

EDGE2_list$scientific_name[which(EDGE2_list$scientific_name=='Scleria rubrostriata A.C.AraÃºjo & N.A.Brummitt')] <- 'Scleria rubrostriata A.C.Araujo & N.A.Brummitt'
scleria_occ$scientific_name[which(scleria_occ$scientific_name=='Scleria rubrostriata A.C.Araújo & N.A.Brummitt')] <- 'Scleria rubrostriata A.C.Araujo & N.A.Brummitt'

EDGE2_list$scientific_name[which(EDGE2_list$scientific_name=='Scleria adpressohirta (KÃ¼k.) E.A.Rob.')] <- 'Scleria adpressohirta (Kuk.) E.A.Rob.'
scleria_occ$scientific_name[which(scleria_occ$scientific_name=='Scleria adpressohirta (Kük.) E.A.Rob.')] <- 'Scleria adpressohirta (Kuk.) E.A.Rob.'


# points file
p_scleria_occ <- vect(scleria_occ, geom=c('x','y'), 'epsg:4326') %>% terra::project('+proj=eqearth') # points file


# plot
# plot(scl_wrld_map)
# points(p_scleria_occ, col='green')
# lines(scl_wrld_map)

# Check NAs
# table(is.na(terra::extract(RAST1_eqearth, p_scleria_occ, xy=T, cell=T)$cell))


# species lists vectors
E2_list <- EDGE2_list$scientific_name[!is.na(EDGE2_list$list)]
E2_main <- EDGE2_list %>% filter(list=='main') %>% dplyr::select(scientific_name) %>% deframe()
E2_25 <- EDGE2_list[order(EDGE2_list$EDGE2, decreasing=T)[1:25],] %>% dplyr::select(scientific_name) %>% deframe()

Eco2_list <- EcoDGE2_list$scientific_name[!is.na(EcoDGE2_list$list)]
Eco2_main <- EcoDGE2_list %>% filter(list=='main') %>% dplyr::select(scientific_name) %>% deframe()
Eco2_25 <- EcoDGE2_list[order(EcoDGE2_list$EcoDGE2, decreasing=T)[1:25],] %>% dplyr::select(scientific_name) %>% deframe()



# countries ####

scl_countries <- scl_wrld_map %>% terra::extract(p_scleria_occ) %>%
  dplyr::select(sovereignt, admin, geounit, subunit, name, name_long, pop_year, gdp_md, gdp_year, economy)
scl_countries <- cbind(scleria_occ[,c('scientific_name', 'x', 'y', 'section', 'subgenus')],
                       scl_countries) %>% select(scientific_name,name) %>% na.omit() %>% unique()
colnames(scl_countries) <- c('scientific_name','country')

write.table(scl_countries, 'results/scl_countries.txt')


EDGE_countries <- matrix(nrow=length(unique(scl_countries$country)), ncol=10) %>% as.data.frame()
colnames(EDGE_countries) <- c('country','richness','EDGE2_list','EDGE2_main', 'EDGE2_top25','sum_EDGE2','EcoDGE2_list','EcoDGE2_main', 'EcoDGE2_top25','sum_EcoDGE2')
EDGE_countries$country <- unique(scl_countries$country)

EDGE_countries$country[EDGE_countries$country=='France'] <- 'French Guiana'
EDGE_countries <- EDGE_countries %>% filter(country!='Canada')

scl_wrld_map$richness <- NA
scl_wrld_map$sum_EDGE2 <- NA
scl_wrld_map$sum_EcoDGE2 <- NA

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
  EDGE_countries$sum_EDGE2[s] <- EDGE2_list %>% filter(scientific_name %in% sp1) %>%
    dplyr::select(EDGE2) %>% colSums()

  # EcoDGE2
  EDGE_countries$EcoDGE2_list[s] <- length(Eco2_list[Eco2_list %in% sp1])
  EDGE_countries$EcoDGE2_main[s] <- length(Eco2_main[Eco2_main %in% sp1])
  EDGE_countries$EcoDGE2_top25[s] <- length(Eco2_25[Eco2_25 %in% sp1])
  EDGE_countries$sum_EcoDGE2[s] <- EcoDGE2_list %>% filter(scientific_name %in% sp1) %>%
    dplyr::select(EcoDGE2) %>% colSums()

  # save in spatvector
  scl_wrld_map$richness[scl_wrld_map$name==EDGE_countries$country[s]] <- EDGE_countries$richness[s]
  scl_wrld_map$sum_EDGE2[scl_wrld_map$name==EDGE_countries$country[s]] <- EDGE_countries$sum_EDGE2[s]
  scl_wrld_map$sum_EcoDGE2[scl_wrld_map$name==EDGE_countries$country[s]] <- EDGE_countries$sum_EcoDGE2[s]
  
  
  print(paste('--- ', round(s/nrow(EDGE_countries)*100, 2), '% ---', sep=''))
  
}

# save
write.table(EDGE_countries, 'results/EDGE_countries.txt')
save(scl_wrld_map, file="results/scl_wrld_map.RData")



# ecoregions ####

scl_ecor_map <- scl_ecor_map %>% terra::extract(p_scleria_occ) %>%
  dplyr::select(REALM, BIOME, ECO_NAME)
scl_ecor_map <- cbind(scleria_occ[,c('scientific_name', 'x', 'y', 'section', 'subgenus')],
                        scl_ecor_map) %>% select(scientific_name, REALM, BIOME, ECO_NAME) %>%
  na.omit() %>% unique()
colnames(scl_ecor_map) <- c('scientific_name','REALM','BIOME','ECO_NAME')

write.table(scl_ecor_map, 'results/scl_ecor_map.txt')


EDGE_scl_ecor_map <- matrix(nrow=length(unique(scl_ecor_map$ECO_NAME)), ncol=10) %>% as.data.frame()
colnames(EDGE_scl_ecor_map) <- c('ECO_NAME','richness','EDGE2_list','EDGE2_main', 'EDGE2_top25','sum_EDGE2','EcoDGE2_list','EcoDGE2_main', 'EcoDGE2_top25','sum_EcoDGE2')
EDGE_scl_ecor_map$ECO_NAME <- unique(scl_ecor_map$ECO_NAME)

scl_ecor_map$richness[s] <- NA
scl_ecor_map$sum_EDGE2[s] <- NA
scl_ecor_map$sum_EcoDGE2[s] <- NA

for (s in 1:nrow(EDGE_scl_ecor_map)) {
  
  # species in a given country
  sp1 <- scl_ecor_map$scientific_name[scl_ecor_map$ECO_NAME==EDGE_scl_ecor_map$ECO_NAME[s]] %>% na.omit() %>%
    unique() %>% as.vector()
  
  # richness
  EDGE_scl_ecor_map$richness[s] <- length(sp1)
  
  # lists
  EDGE_scl_ecor_map$EDGE2_list[s] <- length(E2_list[E2_list %in% sp1])
  EDGE_scl_ecor_map$EDGE2_main[s] <- length(E2_main[E2_main %in% sp1])
  EDGE_scl_ecor_map$EDGE2_top25[s] <- length(E2_25[E2_25 %in% sp1])
  EDGE_scl_ecor_map$sum_EDGE2[s] <- EDGE2_list %>% filter(scientific_name %in% sp1) %>%
    dplyr::select(EDGE2) %>% colSums()
  
  # EcoDGE2
  EDGE_scl_ecor_map$EcoDGE2_list[s] <- length(Eco2_list[Eco2_list %in% sp1])
  EDGE_scl_ecor_map$EcoDGE2_main[s] <- length(Eco2_main[Eco2_main %in% sp1])
  EDGE_scl_ecor_map$EcoDGE2_top25[s] <- length(Eco2_25[Eco2_25 %in% sp1])
  EDGE_scl_ecor_map$sum_EcoDGE2[s] <- EcoDGE2_list %>% filter(scientific_name %in% sp1) %>%
    dplyr::select(EcoDGE2) %>% colSums()
  
  # save variable in SpatVector
  scl_ecor_map$richness[scl_ecor_map$ECO_NAME==EDGE_scl_ecor_map$ECO_NAME[s]] <- EDGE_scl_ecor_map$richness[s]
  scl_ecor_map$sum_EDGE2[scl_ecor_map$ECO_NAME==EDGE_scl_ecor_map$ECO_NAME[s]] <- EDGE_scl_ecor_map$sum_EDGE2[s]
  scl_ecor_map$sum_EcoDGE2[scl_ecor_map$ECO_NAME==EDGE_scl_ecor_map$ECO_NAME[s]] <- EDGE_scl_ecor_map$sum_EcoDGE2[s]
  
  
  print(paste('--- ', round(s/nrow(EDGE_scl_ecor_map)*100, 2), '% ---', sep=''))
  
}

# save
write.table(EDGE_scl_ecor_map, 'results/EDGE_scl_ecor_map.txt')
save(scl_ecor_map, file="results/scl_ecor_map.RData")



# exploratory ####


# plot
scl_wrld_map$log_richness <- log(scl_wrld_map$richness)
terra::plot(scl_ecor_map, y='richness', breaks=30, col=rev(grDevices::heat.colors(30)), border=NA)
lines(world_lines)

# plot
scl_ecor_map$log_sum_EDGE2 <- log(scl_ecor_map$sum_EDGE2 + 1)
terra::plot(scl_ecor_map, y='log_sum_EDGE2', breaks=30, col=rev(grDevices::heat.colors(30)), border=NA)
lines(world_lines)

# table for manuscript
df_countries <- world_EDGE %>% dplyr::select(region, richness, EDGE2_list, EDGE2_main, EDGE2_top25, sum_EDGE2, EcoDGE2_list, EcoDGE2_main, EcoDGE2_top25, sum_EcoDGE2) %>%
  unique() %>%
  filter(!(is.na(richness)))
df_countries[is.na(df_countries)] <- 0



