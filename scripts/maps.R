
library(tidyverse)
library(rnaturalearth)
library(terra)
library(ggplot2)
library(viridis)
library(ggpubr)
library(picante)
source('scripts/traits.R')


# world map
wrld <- ne_countries(scale = 110, type = "countries", continent = NULL, country = NULL, geounit = NULL, sovereignty = NULL, returnclass = "sf") %>% terra::vect() %>% terra::crop(ext(-180,180,-50,70)) %>% terra::project('epsg:8857')

# base map
RAST1_lonlat <- rast('C:/Users/user/Desktop/wc2.1_5m_bio/wc2.1_5m_bio_1.tif')
RAST1_lonlat[!is.na(RAST1_lonlat)] <- 0

RAST1_8857 <- RAST1_lonlat %>% terra::project('epsg:8857') %>% # equal area, see https://epsg.io/8857
  terra::aggregate(fact=18.26152, fun='min') %>% terra::crop(wrld, mask=T)
terra::cellSize(RAST1_8857) # 100 km aprox

plot(RAST1_lonlat); lines(wrld)


# DATA IMPORT AND PIXEL CALCULATION ####

# data
EDGE2_list <- read.csv("results/EDGE2_list.csv") # data
EcoDGE2_list <- read.csv("results/EcoDGE2_list.csv")

# occurrences
scleria_occ <- read.csv("results/occ_filtered.txt", sep="") %>% dplyr::select(species,y,x)
table(EDGE2_list$species %in% scleria_occ$species) # check missing species
p_scleria_occ <- vect(scleria_occ, geom=c('x','y'), 'epsg:4326') %>% terra::project('epsg:8857') # points file

# extract cell values
scleria_occ$cell <- terra::extract(RAST1_8857, p_scleria_occ, xy=T, cell=T)$cell
scleria_occ$x <- terra::extract(RAST1_8857, p_scleria_occ, xy=T, cell=T)$x
scleria_occ$y <- terra::extract(RAST1_8857, p_scleria_occ, xy=T, cell=T)$y

head(scleria_occ)
is.na(scleria_occ$cell) %>% table()
scleria_occ <- scleria_occ %>% filter(!(is.na(cell))) # remove NA cells

# plot
par(mfrow=c(1,1))
plot(RAST1_8857, col='lightgreen'); lines(wrld); points(p_scleria_occ)



# for each cell, I compute: richness, Faith, MPD, MFD, mean_EDGE, mean_EcoDGE, sum_EDGE, sum_EcoDGE

# DATA
EDGE2_list
EcoDGE2_list
daisy.mat
imputed_trees


# results
df_scleria_map <- scleria_occ[c('cell','x','y')] %>% unique()
df_scleria_map$richness <- NA
df_scleria_map$Faith <- NA
df_scleria_map$MPD <- NA
df_scleria_map$MFD <- NA
df_scleria_map$mean_EDGE <- NA
df_scleria_map$mean_EcoDGE <- NA
df_scleria_map$sum_EDGE <- NA
df_scleria_map$sum_EcoDGE <- NA

# community data matrix for 'picante'
cdm_scleria <- scleria_occ[c('species','cell')] %>% group_by(species,cell) %>% summarise(count=n())
cdm_scleria$species <- gsub(' ', '_', cdm_scleria$species)
cdm_scleria <- cdm_scleria %>% pivot_wider(names_from=species, values_from=count) %>% as.data.frame()
rownames(cdm_scleria) <- cdm_scleria$cell
cdm_scleria$cell <- NULL
cdm_scleria <- cdm_scleria[order(rownames(cdm_scleria)),]
cdm_scleria[is.na(cdm_scleria)] <- 0 # presence/absence
cdm_scleria[cdm_scleria > 0] <- 1
table(rowSums(cdm_scleria)) # histogram nº of observations

# trait mean pairwise diss. for MFD
ss_tr <- as.matrix(daisy.mat)
ss_tr[upper.tri(ss_tr)] <- NA
rownames(ss_tr) <- gsub(' ', '_', rownames(ss_tr))
colnames(ss_tr) <- gsub(' ', '_', colnames(ss_tr))
diag(ss_tr) <- NA

# for phylogenetic indexes the values across imputed trees are equal
tree1 <- imputed_trees[[1]]
tree1$tip.label <- scleria_iucn$species[order(match(scleria_iucn$sectxspp, tree1$tip.label))] # change labels
tree1$tip.label <- gsub(' ', '_', tree1$tip.label)
cph1 <- cophenetic(tree1)

# obtain indexes per pixel
for (r in 1:nrow(df_scleria_map)) {
  
  # species
  ss_sp <- scleria_occ %>% filter(cell==df_scleria_map$cell[r]) %>% dplyr::select(species) %>% unique() %>% deframe()
  
  # richness
  df_scleria_map$richness[r] <- length(ss_sp)
  
  # EDGE2
  df_scleria_map$sum_EDGE[r] <- EDGE2_list %>% filter(species %in% ss_sp) %>% dplyr::select(EDGE2) %>% deframe() %>% sum(na.rm=T)
  df_scleria_map$mean_EDGE[r] <- EDGE2_list %>% filter(species %in% ss_sp) %>% dplyr::select(EDGE2) %>% deframe() %>% mean(na.rm=T)

  # EDGE2
  df_scleria_map$sum_EcoDGE[r] <- EcoDGE2_list %>% filter(species %in% ss_sp) %>% dplyr::select(EDGE2) %>% deframe() %>% sum(na.rm=T)
  df_scleria_map$mean_EcoDGE[r] <- EcoDGE2_list %>% filter(species %in% ss_sp) %>% dplyr::select(EDGE2) %>% deframe() %>% mean(na.rm=T)

  ss_sp <- gsub(" ", "_", ss_sp) # check names
  ss_cm <- cdm_scleria[as.character(df_scleria_map$cell[r]),] # community matrix
  
  ss_tr_r <- ss_tr[ss_sp,ss_sp] # unweighted mean pairwise distance MFD: de Bello et al. 2016 (Oecologia)
  df_scleria_map$MFD[r] <- mean(ss_tr_r, na.rm=T)
  
  df_scleria_map$Faith[r] <-  pd(samp=ss_cm, tree=tree1, include.root=T)$PD
  
  cph_n <- cph1[ss_sp,ss_sp]
  df_scleria_map$MPD[r] <-  cph_n[lower.tri(cph_n)] %>% mean()
  
  print(r)

}

head(df_scleria_map)
colSums(is.na(df_scleria_map))

df_scleria_map[is.na(df_scleria_map)] <- 0 # replace NA by 0 for plotting
write.table(df_scleria_map, 'results/df_scleria_map.txt')

ggpairs(df_scleria_map[,c("richness","Faith","MPD","mean_EDGE","sum_EDGE","MFD","mean_EcoDGE","sum_EcoDGE")])

colnames(df_scleria_map)


# create rasters
m_richness <- RAST1_8857
m_faith <- RAST1_8857
m_MPD <- RAST1_8857
m_EDGE <- RAST1_8857
m_MFD <- RAST1_8857
m_EcoDGE <- RAST1_8857

# remove diversity and values 
df_scleria_map

for (r in c(1:nrow(df_scleria_map))) {
  m_richness[df_scleria_map$cell[r]] <- df_scleria_map$richness[r]
  m_faith[df_scleria_map$cell[r]] <- df_scleria_map$Faith[r]
  m_MPD[df_scleria_map$cell[r]] <- df_scleria_map$MPD[r]
  m_EDGE[df_scleria_map$cell[r]] <- df_scleria_map$mean_EDGE[r]
  m_MFD[df_scleria_map$cell[r]] <- df_scleria_map$MFD[r]
  m_EcoDGE[df_scleria_map$cell[r]] <- df_scleria_map$mean_EcoDGE[r]
  print(r)
}


# MAPS ####


par(mfrow=c(1,1))
plot(log(m_richness), main='log Richness'); lines(wrld)
plot(log(m_faith), main='log Faith'); lines(wrld)
plot(m_MPD, main='MPD'); lines(wrld)
plot(m_MFD, main='MFD'); lines(wrld)


plot(m_EDGE, main='EDGE2'); lines(wrld)
plot(log(m_EcoDGE+1), main='EcoDGE2'); lines(wrld)



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


