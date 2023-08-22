
library(tidyverse); library(readxl)
library(rnaturalearth)
library(terra)
library(ggplot2)
library(viridis)
library(ggpubr)



# base maps ####

# world map
wrld <- ne_countries(scale = 10, type = "countries", continent = NULL,
                     country = NULL, geounit = NULL, sovereignty = NULL, returnclass = "sf") %>%
  terra::vect() %>% terra::crop(ext(-180,180,-50,70)) %>% terra::project('+proj=eqearth')


# base map
RAST1_lonlat <- rast('C:/Users/user/Desktop/wc2.1_2.5m_elev/wc2.1_2.5m_elev.tif')
RAST1_lonlat[!is.na(RAST1_lonlat)] <- 0
RAST1_eqearth <- RAST1_lonlat %>% terra::project('+proj=eqearth') %>%
  terra::crop(wrld) %>%
  terra::aggregate(fact=40, fun='min') 
terra::cellSize(RAST1_eqearth, unit='km') # 109.5 km2


# biomes map
wwf_biomes <- vect('C:/Users/user/Desktop/ecoregions/wwf_terr_ecos.shp') %>% aggregate(by='BIOME') %>%
  terra::project('+proj=eqearth') %>% terra::crop(RAST1_eqearth)


plot(RAST1_eqearth, legend=F); lines(wrld)
plot(RAST1_eqearth, legend=F); lines(wwf_biomes, col='gray28')



# data import and pixel calculation ####


load('results/data_final.RData') # data


EDGE2_list <- read_excel("results/final_lists.xlsx", sheet="all_EDGE2") # lists
EcoDGE2_list <- read_excel("results/final_lists.xlsx", sheet="all_EcoDGE2")


scleria_occ <- data_final[['occurrences']] %>% # occurrences
  dplyr::select(scientific_name, x, y) %>%
  merge(data_final[['assessments']][,c('scientific_name','section','subgenus')])
scleria_occ$subgenus <- as.factor(scleria_occ$subgenus)
scleria_occ$section <- as.factor(scleria_occ$section)


# check missing species: we need to fix S. rubrostiata name
v_miss <- EDGE2_list$scientific_name[!(EDGE2_list$scientific_name %in% scleria_occ$scientific_name)]
EDGE2_list$scientific_name[which(EDGE2_list$scientific_name==v_miss[2])] <- 'Scleria rubrostriata A.C.Araujo & N.A.Brummitt'
scleria_occ$scientific_name[which(scleria_occ$scientific_name=='Scleria rubrostriata A.C.Araújo & N.A.Brummitt')] <- 'Scleria rubrostriata A.C.Araujo & N.A.Brummitt'
EDGE2_list$scientific_name[!(EDGE2_list$scientific_name %in% scleria_occ$scientific_name)]


p_scleria_occ <- vect(scleria_occ, geom=c('x','y'), 'epsg:4326') %>% terra::project('+proj=eqearth') # points file



# extract cell values and check NAs
scleria_occ$cell <- terra::extract(RAST1_eqearth, p_scleria_occ, xy=T, cell=T)$cell

is.na(scleria_occ$cell) %>% table() # remove NA cells and redo points file
scleria_occ <- scleria_occ %>% filter(!(is.na(cell)))
p_scleria_occ <- vect(scleria_occ, geom=c('x','y'), 'epsg:4326') %>% terra::project('+proj=eqearth')



# plot observations
par(mfrow=c(1,1))
# plot(RAST1_eqearth, col='azure2', legend=F); lines(wwf_biomes, col='gray28')
# points(p_scleria_occ, col=p_scleria_occ$subgenus)
# legend( x="topright", 
#         col=c("cornflowerblue","hotpink2","green3","black"), 
#         legend=c('Trachylomia','Hypoporum','Scleria','Browniae'), lwd=0, lty=0, bty="n",
#         pch=15)



# prepare data for maps ####


# countries
# scl_countries <- wrld %>% terra::extract(p_scleria_occ) %>%
#   dplyr::select(sovereignt, admin, geounit, subunit, name, name_long, pop_year, gdp_md, gdp_year, economy)
# 
# scl_countries <- cbind(scleria_occ[,c('scientific_name', 'x', 'y', 'section', 'subgenus', 'cell', 'biome')],
#                        scl_countries)
# 
# write.table(scl_countries, 'results/scl_countries.txt')

scl_countries <- read.csv("results/scl_countries.txt", sep="")
spp_x_ctr <- scl_countries[,c('scientific_name','name')] %>% na.omit() %>% unique()
colnames(spp_x_ctr) <- c('scientific_name','country')



# species x list vectors
E2_list <- EDGE2_list$scientific_name[!is.na(EDGE2_list$list)]
E2_main <- EDGE2_list %>% filter(list=='main') %>% dplyr::select(scientific_name) %>% deframe()
E2_100 <- EDGE2_list[order(EDGE2_list$EDGE2, decreasing=T)[1:100],] %>% dplyr::select(scientific_name) %>% deframe()

Eco2_list <- EcoDGE2_list$scientific_name[!is.na(EcoDGE2_list$list)]
Eco2_main <- EcoDGE2_list %>% filter(list=='main') %>% dplyr::select(scientific_name) %>% deframe()
Eco2_100 <- EcoDGE2_list[order(EcoDGE2_list$EDGE2, decreasing=T)[1:100],] %>% dplyr::select(scientific_name) %>% deframe()



# results
EDGE_countries <- matrix(nrow=length(unique(spp_x_ctr$country)), ncol=12) %>% as.data.frame()
colnames(EDGE_countries) <- c('country','richness','EDGE2_list','EDGE2_main', 'EDGE2_top100','sum_EDGE2','mean_EDGE2','EcoDGE2_list','EcoDGE2_main', 'EcoDGE2_top100','sum_EcoDGE2','mean_EcoDGE2')
EDGE_countries$country <- unique(spp_x_ctr$country)


for (s in 1:nrow(EDGE_countries)) {
  
  # species in a given country
  sp1 <- spp_x_ctr$scientific_name[spp_x_ctr$country==EDGE_countries$country[s]]
  
  # richness
  EDGE_countries$richness[s] <- length(sp1)

  # lists
  EDGE_countries$EDGE2_list[s] <- length(E2_list[E2_list %in% sp1])
  EDGE_countries$EDGE2_main[s] <- length(E2_main[E2_main %in% sp1])
  EDGE_countries$EDGE2_top100[s] <- length(E2_100[E2_100 %in% sp1])
  EDGE_countries$sum_EDGE2[s] <- EDGE2_list %>% filter(scientific_name %in% sp1) %>%
    dplyr::select(EDGE2) %>% colSums()
  EDGE_countries$mean_EDGE2[s] <- EDGE2_list %>% filter(scientific_name %in% sp1) %>%
    dplyr::select(EDGE2) %>% colMeans()
  
  # EcoDGE2
  EDGE_countries$EcoDGE2_list[s] <- length(Eco2_list[Eco2_list %in% sp1])
  EDGE_countries$EcoDGE2_main[s] <- length(Eco2_main[Eco2_main %in% sp1])
  EDGE_countries$EcoDGE2_top100[s] <- length(Eco2_100[Eco2_100 %in% sp1])
  EDGE_countries$sum_EcoDGE2[s] <- EcoDGE2_list %>% filter(scientific_name %in% sp1) %>%
    dplyr::select(EDGE2) %>% colSums()
  EDGE_countries$mean_EcoDGE2[s] <- EcoDGE2_list %>% filter(scientific_name %in% sp1) %>%
    dplyr::select(EDGE2) %>% colMeans()
  
}

# fix some entries
EDGE_countries$country[EDGE_countries$country=='France'] <- 'French Guiana'
EDGE_countries <- EDGE_countries %>% filter(country!='Canada')

# save
# write.table(EDGE_countries, 'results/EDGE_countries.txt')



# exploratory ####

# par(mfrow=c(2,5), mar=c(5,5,5,5))
# pairs(EDGE_countries[,2:12], upper.panel=NULL)


# prepare results to merge with world dataframe
colnames(EDGE_countries)[colnames(EDGE_countries)=='country'] <- 'region'

# match names of ne_countries and map_world for reference
world <- map_data("world")

EDGE_countries$region[!(EDGE_countries$region %in% world$region)]

temp_c <- c('Democratic Republic of the Congo','Ivory Coast','USA','Republic of Congo','Equatorial Guinea','Central African Republic','Trinidad','Dominican Republic','Swaziland','South Sudan','Cuba','Cayman Islands')

EDGE_countries$region[!(EDGE_countries$region %in% world$region)] <- temp_c



# replace 0 with NA
EDGE_countries[EDGE_countries==0] <- NA


# merge with world dataframe
world_EDGE <- map_data("world")
world_EDGE$region[world_EDGE$subregion=='Alaska'] <- 'NA' # remove Alaska
world_EDGE <- merge(world_EDGE, EDGE_countries, by="region", all.x=T)
world_EDGE <- world_EDGE %>% arrange(group, order)


# richness
ggplot(aes(x=long, y=lat, group=group, fill=richness), data=world_EDGE) +
  geom_polygon() + theme_bw() +
  geom_polygon(color="white", size=0.2) +
  labs(x='', y='') + ggtitle('Richness') +
  scale_fill_viridis('', option='inferno', na.value = "gray80")


# EDGE2
g1 <- ggplot(aes(x=long, y=lat, group=group, fill=EDGE2_list), data=world_EDGE) +
  geom_polygon() + theme_bw() +
  geom_polygon(color="white", size=0.2) +
  labs(x='', y='') + ggtitle('Nº listed species EDGE2') +
  scale_fill_viridis('', option='inferno', na.value = "gray80")

g2 <- ggplot(aes(x=long, y=lat, group=group, fill=EDGE2_main), data=world_EDGE) +
  geom_polygon() + theme_bw() +
  geom_polygon(color="white", size=0.2) +
  labs(x='', y='') + ggtitle('Nº species main list EDGE2') +
  scale_fill_viridis('', option='inferno', na.value = "gray80")

g3 <- ggplot(aes(x=long, y=lat, group=group, fill=sum_EDGE2), data=world_EDGE) +
  geom_polygon() + theme_bw() +
  geom_polygon(color="white", size=0.2) +
  labs(x='', y='') + ggtitle('Sum EDGE2') +
  scale_fill_viridis('', option='inferno', na.value = "gray80")

g4 <- ggplot(aes(x=long, y=lat, group=group, fill=mean_EDGE2), data=world_EDGE) +
  geom_polygon() + theme_bw() +
  geom_polygon(color="white", size=0.2) +
  labs(x='', y='') + ggtitle('Mean EDGE2') +
  scale_fill_viridis('', option='inferno', na.value = "gray80")

ggarrange(g1, g2, g3, g4, nrow=2, ncol=2, labels=c('A','B','C','D'))


# EcoDGE2
g1 <- ggplot(aes(x=long, y=lat, group=group, fill=EcoDGE2_list), data=world_EDGE) +
  geom_polygon() + theme_bw() +
  geom_polygon(color="white", size=0.2) +
  labs(x='', y='') + ggtitle('Nº listed species EcoDGE2') +
  scale_fill_viridis('', option='inferno', na.value = "gray80")

g2 <- ggplot(aes(x=long, y=lat, group=group, fill=EcoDGE2_main), data=world_EDGE) +
  geom_polygon() + theme_bw() +
  geom_polygon(color="white", size=0.2) +
  labs(x='', y='') + ggtitle('Nº species main list EcoDGE2') +
  scale_fill_viridis('', option='inferno', na.value = "gray80")

g3 <- ggplot(aes(x=long, y=lat, group=group, fill=sum_EcoDGE2), data=world_EDGE) +
  geom_polygon() + theme_bw() +
  geom_polygon(color="white", size=0.2) +
  labs(x='', y='') + ggtitle('Sum EcoDGE2') +
  scale_fill_viridis('', option='inferno', na.value = "gray80")

g4 <- ggplot(aes(x=long, y=lat, group=group, fill=mean_EcoDGE2), data=world_EDGE) +
  geom_polygon() + theme_bw() +
  geom_polygon(color="white", size=0.2) +
  labs(x='', y='') + ggtitle('Mean EcoDGE2') +
  scale_fill_viridis('', option='inferno', na.value = "gray80")

ggarrange(g1, g2, g3, g4, nrow=2, ncol=2, labels=c('A','B','C','D'))

