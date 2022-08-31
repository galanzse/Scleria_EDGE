
library(tidyverse)
library(rnaturalearth)
library(terra)
library(ggplot2)
library(viridis)
library(ggpubr)


# world map
wrld <- ne_countries(scale = 110, type = "countries", continent = NULL, country = NULL, geounit = NULL, sovereignty = NULL, returnclass = "sf") %>% terra::vect()

# occurrences from JofBio
scleria_occ <- read.csv("data/alldata_jofbiogeo.txt", sep="") %>% dplyr::select(species,latitudecorrected,longitudecorrected)
colnames(scleria_occ)[colnames(scleria_occ)%in%c('latitudecorrected','longitudecorrected')] <- c('y','x')
scleria_occ$species <- gsub(" ", "_", scleria_occ$species)
scleria_occ$species[scleria_occ$species=='Scleria_sobolifera'] <- 'Scleria_sobolifer'

# EDGE2 results
EDGE2_list <- read.csv("results/EDGE2_list.csv")


# check missing species
table(EDGE2_list$species %in% scleria_occ$species)
EDGE2_list$species[!(EDGE2_list$species %in% scleria_occ$species)]

# add species to dataset (country level)
# scleria_occ <- rbind(scleria_occ,
#       c('Scleria_abortiva', -19.378200204815556, 46.4684760912838),
#       c("Scleria_cheekii", 5.05823925656733, 12.680990638893737),
#       c("Scleria_guineensis", 10.40328327453429, -10.064662907262854))
# str(scleria_occ)
# scleria_occ$y <- as.numeric(scleria_occ$y)
# scleria_occ$x <- as.numeric(scleria_occ$x)

# points file
p_scleria_occ <- vect(scleria_occ, geom=c('x','y'), 'epsg:4326')


# 
par(mfrow=c(1,1))
plot(wrld); points(p_scleria_occ)

# map EDGE2 to raster
EDGE_rast <- terra::rast('E:/wc2.1_2.5m_bio/wc2.1_2.5m_bio_1.tif') %>% terra::aggregate(50, fun='modal')
EDGE_rast[!is.na(EDGE_rast)] <- 0

plot(EDGE_rast); lines(wrld)
points(p_scleria_occ)

EDGE_rast_sum <- EDGE_rast
EDGE_rast_med <- EDGE_rast
EDGE_rast_list <- EDGE_rast
EDGE_rast_t25 <- EDGE_rast

scleria_occ$cell <- terra::extract(EDGE_rast, p_scleria_occ, xy=T, cell=T)$cell
E_list <- EDGE2_list$species[!is.na(EDGE2_list$list)]
top25 <- EDGE2_list$species[order(EDGE2_list$EDGE2, decreasing=T)][1:25]

for (c in unique(scleria_occ$cell)) {
  ss1 <- scleria_occ %>% filter(cell==c) %>% dplyr::select(species) %>% deframe() %>% unique() %>% as.vector()
  EDGE_rast_sum[c] <- EDGE2_list$EDGE2[EDGE2_list$species %in% ss1] %>% sum()
  EDGE_rast_med[c] <- EDGE2_list$EDGE2[EDGE2_list$species %in% ss1] %>% median()
  EDGE_rast_list[c] <- length(E_list[E_list %in% ss1])
  EDGE_rast_t25[c] <- length(top25[top25 %in% ss1])
}

par(mfrow=c(2,1))
plot(log(EDGE_rast_sum)+4); lines(wrld)
plot(log(EDGE_rast_med)+4); lines(wrld)

EDGE_rast_list[EDGE_rast_list==0] <- NA 
EDGE_rast_t25[EDGE_rast_t25==0] <- NA 
plot(EDGE_rast_list); lines(wrld) # better at the country level
plot(EDGE_rast_t25); lines(wrld) # better at the country level


# countries
scleria_occ$region <- wrld %>% terra::extract(p_scleria_occ) %>% dplyr::select(name) %>% as.vector() %>% deframe()
scleria_occ <- na.omit(scleria_occ)

# list of countries with greater EDGE2 scores
spp_x_ctr <- unique(scleria_occ[,c('species','region')])
EDGE_x_ctr <- matrix(nrow=length(unique(spp_x_ctr$region)), ncol=6) %>% as.data.frame()
colnames(EDGE_x_ctr) <- c('region','richness','sum_EDGE','median_EDGE','n_list','top25')
EDGE_x_ctr$region <- unique(spp_x_ctr$region)

for (s in 1:nrow(EDGE_x_ctr)) {
  sp1 <- spp_x_ctr$species[spp_x_ctr$region==EDGE_x_ctr$region[s]]
  EDGE_x_ctr$n_list[s] <- length(E_list[E_list %in% sp1])
  EDGE_x_ctr$top25[s] <- length(top25[top25 %in% sp1])
  
  EDGE_x_ctr$richness[s] <- length(sp1)
  ed1 <- EDGE2_list$EDGE2[EDGE2_list$species %in% sp1]
  EDGE_x_ctr$sum_EDGE[s] <- sum(ed1)
  EDGE_x_ctr$median_EDGE[s] <- median(ed1)
}

EDGE_x_ctr$region[EDGE_x_ctr$region=='France'] <- 'French Guiana'
EDGE_x_ctr <- EDGE_x_ctr %>% filter(region!='Canada')

write.table(EDGE_x_ctr, 'results/EDGE_x_ctr.txt')


# exploratory
par(mfrow=c(2,3), mar=c(5,5,5,5))
plot(sum_EDGE ~ richness, EDGE_x_ctr)
plot(median_EDGE ~ richness, EDGE_x_ctr[EDGE_x_ctr$richness>5,])
plot(median_EDGE ~ sum_EDGE, EDGE_x_ctr[EDGE_x_ctr$richness>5,])
plot(n_list ~ top25, EDGE_x_ctr)
plot(n_list ~ richness, EDGE_x_ctr)
plot(top25 ~ richness, EDGE_x_ctr)


# listed and top25 EDGE2 Scleria x country
world <- map_data("world")
world <- merge(world, EDGE_x_ctr, by="region", all.x=T)
world <- world %>% arrange(group, order)

# match names of ne_countries and map_world
EDGE_x_ctr$region[!(EDGE_x_ctr$region %in% world$region)] <- c('USA','Democratic Republic of the Congo','Republic of Congo','Ivory Coast','Dominican Republic','Central African Republic','Trinidad','Laos', 'Equatorial Guinea','Solomon Islands')

g1 <- ggplot(aes(x=long, y=lat, group=group, fill=n_list), data=world) +
  geom_polygon() + theme_bw() +
  geom_polygon(color="white", size=0.2) +
  labs(x='', y='') +
  scale_fill_viridis('', option='inferno', na.value = "gray90")

g2 <- ggplot(aes(x=long, y=lat, group=group, fill=top25), data=world) +
  geom_polygon() + theme_bw() +
  geom_polygon(color="white", size=0.2) +
  labs(x='', y='') +
  scale_fill_viridis('', option='inferno', na.value = "gray90")

ggarrange(g2, g1, nrow=2, labels=c('A','B'))


