

# EXPLORE THE PROPORTION OF THREATENED SPECIES PER BIOME AND SECTION 


library(tidyverse)
library(terra)
library(ggpubr)
library(readxl)


# results automated assessments
comb_AA <- read_excel("results/results.xlsx", sheet = "combined_assessments")
comb_AA <- comb_AA[,!(colnames(comb_AA)%in%c('Notes','Rcat'))]
str(comb_AA)
colnames(comb_AA)[1] <- 'scientific_name'

comb_AA$RF[is.na(comb_AA$RF) & comb_AA$IUCN%in%c('LC','NT')] <- 'NT'
comb_AA$RF[is.na(comb_AA$RF) & comb_AA$IUCN%in%c('CR','VU','EN')] <- 'T'

comb_AA$RF[comb_AA$RF=='NT'] <- 'nonthreatened'
comb_AA$RF[comb_AA$RF=='T'] <- 'threatened'



# occurrences
load("results/data_final.RData") # data
pts_scl <- data_final[['occurrences']][,c('scientific_name','x','y')] %>% vect(geom=c('x','y'), 'epsg:4326')


# ecoregions
ecoregions <- vect('C:/Users/user/Desktop/ecoregions/wwf_terr_ecos.shp')
df_ecoregions <- ecoregions %>%
  terra::intersect(y=pts_scl) %>%
  as.data.frame() %>%
  dplyr::select(scientific_name, REALM, BIOME)


# merge datasets
comb_AA <- merge(comb_AA, df_ecoregions)
comb_AA <- merge(comb_AA, data_final[['taxa']][,c('scientific_name','subgenus','section')])
comb_AA$IUCN <- NULL


# change factor levels
comb_AA$REALM <- as.factor(comb_AA$REALM)
levels(comb_AA$REALM) <- c("Australasia", "Afrotropics", "IndoMalay", "Nearctic", "Neotropics", "Oceania", "Palearctic")

comb_AA$BIOME <- as.factor(comb_AA$BIOME)
levels(comb_AA$BIOME) <- c("Tropical & Subtropical Moist Broadleaf Forests", "Tropical & Subtropical Dry Broadleaf Forests", "Tropical & Subtropical Coniferous Forests", "Temperate Broadleaf & Mixed Forests", "Temperate Conifer Forests", "Tropical & Subtropical Grasslands, Savannas & Shrublands", "Temperate Grasslands, Savannas & Shrublands", "Flooded Grasslands & Savannas", "Montane Grasslands & Shrublands", "Deserts & Xeric Shrublands", "Mangroves", NA)



# plots
t1 <- comb_AA[,c('scientific_name', 'REALM', 'RF')] %>% unique() %>%
  group_by(REALM, RF) %>% summarise(n_spp=n()) %>% na.omit()
t1 <- rbind(t1, data.frame(REALM=c('Nearctic','Oceania','Palearctic'),
                           RF=c('threatened','threatened','threatened'),
                           n_spp=c(0,0,0)))

g1 <- ggplot(aes(x=REALM, y=n_spp, fill=RF), data=t1) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values=c("forestgreen","red")) +
  theme_bw() +
  theme(legend.position="top", legend.title=element_blank()) +
  ylab('Number of species') + xlab('Realm')



t2 <- comb_AA[,c('scientific_name', 'section', 'RF')] %>% unique() %>%
  group_by(section, RF) %>% summarise(n_spp=n()) %>% na.omit()
t2 <- rbind(t2, data.frame(section=c('Lithospermae','Margaleia','Melanomphalae','Trachylomia'),
                           RF=c('threatened','threatened','threatened','threatened'),
                           n_spp=c(0,0,0,0)))

g2 <- ggplot(aes(x=section, y=n_spp, fill=RF), data=t2) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values=c("forestgreen","red")) +
  theme_bw() +
  theme(legend.position="none",
        axis.text.x=element_text(angle = -45, vjust=1, hjust=0)) +
  ylab('Number of species') + xlab('Section')
  


ggarrange(g1, g2, nrow=2, labels=c('a)','b)'))


