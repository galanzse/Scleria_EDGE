
library(tidyverse)
library(rnaturalearth)
library(terra)
library(ggplot2)
library(viridis)
library(ggpubr)
library(picante)
source('scripts/traits.R')


# world map
wrld <- ne_countries(scale = 110, type = "countries", continent = NULL, country = NULL, geounit = NULL, sovereignty = NULL, returnclass = "sf") %>%
  terra::vect() %>% terra::crop(ext(-180,180,-50,70)) %>% terra::project('epsg:8857')

# base map
RAST1_lonlat <- rast('C:/Users/user/Desktop/wc2.1_5m_bio/wc2.1_5m_bio_1.tif')
RAST1_lonlat[!is.na(RAST1_lonlat)] <- 0

RAST1_8857 <- RAST1_lonlat %>% terra::project('epsg:8857') %>% # equal area, see https://epsg.io/8857
  terra::aggregate(fact=18.26152, fun='min') %>% terra::crop(wrld)
terra::cellSize(RAST1_8857) # 100 km aprox

# biomes map
wwf_biomes <- terra::vect('C:/Users/user/Desktop/WWF_regions/wwf_terr_ecos.shp') %>% aggregate(by='BIOME') %>%
  terra::project('epsg:8857') %>% terra::crop(RAST1_8857)

plot(RAST1_8857, legend=F); lines(wrld)
plot(RAST1_8857, legend=F); lines(wwf_biomes, col='gray28')


# DATA IMPORT AND PIXEL CALCULATION ####

# data
EDGE2_list <- read.csv("results/EDGE2_list.csv") # data
EcoDGE2_list <- read.csv("results/EcoDGE2_list.csv")

# occurrences
scleria_occ <- read.csv("results/occ_filtered.txt", sep="")
scleria_occ <- merge(scleria_occ, scleria_iucn[,c('species','section','subgenus')])
str(scleria_occ)
scleria_occ$subgenus <- as.factor(scleria_occ$subgenus)
scleria_occ$section <- as.factor(scleria_occ$section)

EDGE2_list$species[!(EDGE2_list$species %in% scleria_occ$species)] # check missing species
p_scleria_occ <- vect(scleria_occ, geom=c('x','y'), 'epsg:4326') %>% terra::project('epsg:8857') # points file

# extract cell values and check NAs
scleria_occ$cell <- terra::extract(RAST1_8857, p_scleria_occ, xy=T, cell=T)$cell
is.na(scleria_occ$cell) %>% table()

scleria_occ <- scleria_occ %>% filter(!(is.na(cell))) # remove NA cells and redo points file
p_scleria_occ <- vect(scleria_occ, geom=c('x','y'), 'epsg:4326') %>% terra::project('epsg:8857') # points file
scleria_occ$cell <- terra::extract(RAST1_8857, p_scleria_occ, xy=T, cell=T)$cell
scleria_occ$x <- terra::extract(RAST1_8857, p_scleria_occ, xy=T, cell=T)$x
scleria_occ$y <- terra::extract(RAST1_8857, p_scleria_occ, xy=T, cell=T)$y


# plot observations
par(mfrow=c(1,1))
plot(RAST1_8857, col='azure2', legend=F); lines(wwf_biomes, col='gray28'); points(p_scleria_occ, col=p_scleria_occ$subgenus)
legend( x="topright", 
        col=c("cornflowerblue","hotpink2","green3","black"), 
        legend=c('Trachylomia','Hypoporum','Scleria','Browniae'), lwd=0, lty=0, bty="n",
        pch=15)


# extract biome information
biome_rast <- rasterize(x=wwf_biomes, y= disagg(RAST1_8857, 10), field='BIOME')
biome_rast[biome_rast>14] <- NA
biome_rast <- as.factor(biome_rast)
levels(biome_rast)[[1]]['label'] <- c('Tropical & Subtropical Moist Broadleaf Forests',
                                      'Tropical & Subtropical Dry Broadleaf Forests',
                                      'Tropical & Subtropical Coniferous Forests',
                                      'Temperate Broadleaf & Mixed Forests',
                                      'Temperate Conifer Forests',
                                      'Boreal Forests/Taiga',
                                      'Tropical & Subtropical Grasslands, Savannas & Shrublands',
                                      'Temperate Grasslands, Savannas & Shrublands',
                                      'Flooded Grasslands & Savannas',
                                      'Montane Grasslands & Shrublands',
                                      'Tundra',
                                      'Mediterranean Forests, Woodlands & Scrub',
                                      'Deserts & Xeric Shrublands',
                                      'Mangroves')
plot(biome_rast)

# extract
scleria_occ$biome <- terra::extract(biome_rast, p_scleria_occ)$label
t1 <- table(scleria_occ$subgenus, scleria_occ$biome)/rowSums(table(scleria_occ$subgenus, scleria_occ$biome))
t1 <- t1 %>% as.data.frame() %>% pivot_wider(names_from=Var2, values_from=Freq)
write.table(t1, 'results/subgenus_x_biome.txt')
t1 <- table(scleria_occ$section, scleria_occ$biome)/rowSums(table(scleria_occ$section, scleria_occ$biome))
t1 <- t1 %>% as.data.frame() %>% pivot_wider(names_from=Var2, values_from=Freq)
write.table(t1, 'results/section_x_biome.txt')

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
plot(m_richness, main='Richness'); lines(wrld)
plot(m_faith, main='Faith'); lines(wrld)
plot(m_MPD, main='MPD'); lines(wrld)
plot(m_MFD, main='MFD'); lines(wrld)


# countries
scleria_occ$country <- wrld %>% terra::extract(p_scleria_occ) %>% dplyr::select(name) %>% as.vector() %>% deframe()
spp_x_ctr <- scleria_occ[,c('species','country')] %>% na.omit() %>% unique()

# add indexes to country
EDGE_countries <- matrix(nrow=length(unique(spp_x_ctr$country)), ncol=6) %>% as.data.frame()
colnames(EDGE_countries) <- c('country','richness','sum_EDGE','median_EDGE','n_list','n_main')
EDGE_countries$country <- unique(spp_x_ctr$country)
EcoDGE_countries <- EDGE_countries

# lists
E2_list <- EDGE2_list$species[!is.na(EDGE2_list$list)]
Eco2_list <- EcoDGE2_list$species[!is.na(EcoDGE2_list$list)]
E2_main <- EDGE2_list %>% filter(list=='main') %>% dplyr::select(species) %>% deframe()
Eco2_main <- EcoDGE2_list %>% filter(list=='main') %>% dplyr::select(species) %>% deframe()

for (s in 1:nrow(EDGE_countries)) {
  
  sp1 <- spp_x_ctr$species[spp_x_ctr$country==EDGE_countries$country[s]]
  EDGE_countries$n_list[s] <- length(E2_list[E2_list %in% sp1])
  EDGE_countries$n_main[s] <- length(E2_main[E2_main %in% sp1])
  EDGE_countries$richness[s] <- length(sp1)
  ed1 <- EDGE2_list$EDGE2[EDGE2_list$species %in% sp1]
  EDGE_countries$sum_EDGE[s] <- sum(ed1)
  EDGE_countries$median_EDGE[s] <- median(ed1)
  
  sp1 <- spp_x_ctr$species[spp_x_ctr$country==EcoDGE_countries$country[s]]
  EcoDGE_countries$n_list[s] <- length(Eco2_list[Eco2_list %in% sp1])
  EcoDGE_countries$n_main[s] <- length(Eco2_main[Eco2_main %in% sp1])
  EcoDGE_countries$richness[s] <- length(sp1)
  ed1 <- EcoDGE2_list$EDGE2[EcoDGE2_list$species %in% sp1]
  EcoDGE_countries$sum_EDGE[s] <- sum(ed1)
  EcoDGE_countries$median_EDGE[s] <- median(ed1)
  
}

# fix some entries
EDGE_countries$country[EDGE_countries$country=='France'] <- 'French Guiana'
EDGE_countries <- EDGE_countries %>% filter(country!='Canada')
EcoDGE_countries$country[EcoDGE_countries$country=='France'] <- 'French Guiana'
EcoDGE_countries <- EcoDGE_countries %>% filter(country!='Canada')

# save
write.table(EDGE_countries, 'results/EDGE_countries.txt')
write.table(EcoDGE_countries, 'results/EcoDGE_countries.txt')
ind_x_countries <- merge(x=EDGE_countries, y=EcoDGE_countries, by=c('country','richness'))
colnames(ind_x_countries) <- c("country","richness","sum_EDGE","median_EDGE","n_list_EDGE","n_main_EDGE","sum_EcoDGE","median_EcoDGE","n_list_EcoDGE","n_main_EcoDGE")
write.table(ind_x_countries, 'results/ind_x_countries.txt')



# exploratory
par(mfrow=c(2,5), mar=c(5,5,5,5))
plot(n_main ~ n_list, ylab='Nº species main list', xlab='Nº listed species', main='EDGE2', data=EDGE_countries)
plot(n_list ~ richness, ylab='Nº listed species', xlab='Richness', main='EDGE2', data=EDGE_countries)
plot(n_main ~ richness, ylab='Nº species main list', xlab='Richness', main='EDGE2', data=EDGE_countries)
plot(log(sum_EDGE) ~ log(sum_EcoDGE), data=ind_x_countries)
plot(log(median_EDGE) ~ log(median_EcoDGE), data=ind_x_countries)
plot(n_main ~ n_list, ylab='Nº species main list', xlab='Nº listed species', main='EcoDGE2', data=EcoDGE_countries)
plot(n_list ~ richness, ylab='Nº listed species', xlab='Richness', main='EcoDGE2', data=EcoDGE_countries)
plot(n_main ~ richness, ylab='Nº species main list', xlab='Richness', main='EcoDGE2', data=EcoDGE_countries)
plot(n_list_EDGE ~ n_list_EcoDGE, data=ind_x_countries)
plot(n_main_EDGE ~ n_main_EcoDGE, data=ind_x_countries)


# prepare results to merge with world dataframe
colnames(EDGE_countries)[colnames(EDGE_countries)=='country'] <- 'region'
colnames(EcoDGE_countries)[colnames(EcoDGE_countries)=='country'] <- 'region'

# match names of ne_countries and map_world for reference
world <- map_data("world")
temp_c <- c('Ivory Coast','Democratic Republic of the Congo','USA','Republic of Congo','Equatorial Guinea','Trinidad','Central African Republic','Dominican Republic','Laos','Solomon Islands','South Korea')
EDGE_countries$region[!(EDGE_countries$region %in% world$region)] <- temp_c
EcoDGE_countries$region[!(EcoDGE_countries$region %in% world$region)] <- temp_c

# replace 0 with NA
EDGE_countries[EDGE_countries==0] <- NA
EcoDGE_countries[EcoDGE_countries==0] <- NA

# EDGE2: merge with world dataframe
world_EDGE <- map_data("world")
world_EDGE$region[world_EDGE$subregion=='Alaska'] <- 'NA' # remove Alaska
world_EDGE <- merge(world_EDGE, EDGE_countries, by="region", all.x=T)
world_EDGE <- world_EDGE %>% arrange(group, order)

# maps
g1 <- ggplot(aes(x=long, y=lat, group=group, fill=richness), data=world_EDGE) +
  geom_polygon() + theme_bw() +
  geom_polygon(color="white", size=0.2) +
  labs(x='', y='') + ggtitle('Richness') +
  scale_fill_viridis('', option='inferno', na.value = "gray90")

g2 <- ggplot(aes(x=long, y=lat, group=group, fill=sum_EDGE), data=world_EDGE) +
  geom_polygon() + theme_bw() +
  geom_polygon(color="white", size=0.2) +
  labs(x='', y='') + ggtitle('Sum EDGE2') +
  scale_fill_viridis('', option='inferno', na.value = "gray90")

g3 <- ggplot(aes(x=long, y=lat, group=group, fill=n_list), data=world_EDGE) +
  geom_polygon() + theme_bw() +
  geom_polygon(color="white", size=0.2) +
  labs(x='', y='') + ggtitle('Nº listed species EDGE2') +
  scale_fill_viridis('', option='inferno', na.value = "gray90")

g4 <- ggplot(aes(x=long, y=lat, group=group, fill=n_main), data=world_EDGE) +
  geom_polygon() + theme_bw() +
  geom_polygon(color="white", size=0.2) +
  labs(x='', y='') + ggtitle('Nº species main list EDGE2') +
  scale_fill_viridis('', option='inferno', na.value = "gray90")

ggarrange(g1, g2, g3, g4, nrow=2, ncol=2, labels=c('A','B','C','D'))


# EcoDGE2: merge with world dataframe
world_EcoDGE <- map_data("world")
world_EcoDGE$region[world_EcoDGE$subregion=='Alaska'] <- 'NA' # remove Alaska
world_EcoDGE <- merge(world_EcoDGE, EcoDGE_countries, by="region", all.x=T)
world_EcoDGE <- world_EcoDGE %>% arrange(group, order)

# maps
g1 <- ggplot(aes(x=long, y=lat, group=group, fill=richness), data=world_EcoDGE) +
  geom_polygon() + theme_bw() +
  geom_polygon(color="white", size=0.2) +
  labs(x='', y='') + ggtitle('Richness') +
  scale_fill_viridis('', option='inferno', na.value = "gray90")

g2 <- ggplot(aes(x=long, y=lat, group=group, fill=sum_EDGE), data=world_EcoDGE) +
  geom_polygon() + theme_bw() +
  geom_polygon(color="white", size=0.2) +
  labs(x='', y='') + ggtitle('Sum EcoDGE2') +
  scale_fill_viridis('', option='inferno', na.value = "gray90")

g3 <- ggplot(aes(x=long, y=lat, group=group, fill=n_list), data=world_EcoDGE) +
  geom_polygon() + theme_bw() +
  geom_polygon(color="white", size=0.2) +
  labs(x='', y='') + ggtitle('Nº listed species EcoDGE2') +
  scale_fill_viridis('', option='inferno', na.value = "gray90")

g4 <- ggplot(aes(x=long, y=lat, group=group, fill=n_main), data=world_EcoDGE) +
  geom_polygon() + theme_bw() +
  geom_polygon(color="white", size=0.2) +
  labs(x='', y='') + ggtitle('Nº species main list EcoDGE2') +
  scale_fill_viridis('', option='inferno', na.value = "gray90")

ggarrange(g1, g2, g3, g4, nrow=2, ncol=2, labels=c('A','B','C','D'))

