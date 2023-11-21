

# EXPLORE THE PROPORTION OF THREATENED SPECIES PER BIOME AND SECTION 


library(tidyverse)
library(terra)
library(ggpubr)


# results automated assessments
df_AA <- read.csv("results/AA_all.txt", sep="")


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
df_AA <- merge(df_AA, df_ecoregions) 
df_AA2 <- merge(df_AA, data_final[['taxa']][,c('scientific_name','subgenus','section')])

df_AA2$Category_CriteriaB <- NULL
df_AA2$rf_threatened[is.na(df_AA2$rf_threatened) & df_AA2$IUCN_category%in%c('LC','NT')] <- 'nonthreatened'
df_AA2$rf_threatened[is.na(df_AA2$rf_threatened) & df_AA2$IUCN_category%in%c('CR','VU','EN')] <- 'threatened'
df_AA2$IUCN_category <- NULL


# change factor levels
df_AA2$REALM <- as.factor(df_AA2$REALM)
levels(df_AA2$REALM) <- c("Australasia", "Afrotropics", "IndoMalay", "Nearctic", "Neotropics", "Oceania", "Palearctic")

df_AA2$BIOME <- as.factor(df_AA2$BIOME)
levels(df_AA2$BIOME) <- c("Tropical & Subtropical Moist Broadleaf Forests", "Tropical & Subtropical Dry Broadleaf Forests", "Tropical & Subtropical Coniferous Forests", "Temperate Broadleaf & Mixed Forests", "Temperate Conifer Forests", "Tropical & Subtropical Grasslands, Savannas & Shrublands", "Temperate Grasslands, Savannas & Shrublands", "Flooded Grasslands & Savannas", "Montane Grasslands & Shrublands", "Deserts & Xeric Shrublands", "Mangroves", NA)



# plots
t1 <- df_AA2[,c('scientific_name', 'REALM', 'rf_threatened')] %>% unique() %>%
  group_by(REALM, rf_threatened) %>% summarise(n_spp=n()) %>% na.omit()
t1 <- rbind(t1, data.frame(REALM=c('Nearctic','Oceania','Palearctic'),
                           rf_threatened=c('threatened','threatened','threatened'),
                           n_spp=c(0,0,0)))

g1 <- ggplot(aes(x=REALM, y=n_spp, fill=rf_threatened), data=t1) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values=c("forestgreen","red")) +
  theme_bw() +
  theme(legend.position="top", legend.title=element_blank()) +
  ylab('Number of species') + xlab('Realm')



t2 <- df_AA2[,c('scientific_name', 'section', 'rf_threatened')] %>% unique() %>%
  group_by(section, rf_threatened) %>% summarise(n_spp=n()) %>% na.omit()
t2 <- rbind(t2, data.frame(section=c('Lithospermae','Margaleia','Melanomphalae'),
                           rf_threatened=c('threatened','threatened','threatened'),
                           n_spp=c(0,0,0)))

g2 <- ggplot(aes(x=section, y=n_spp, fill=rf_threatened), data=t2) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values=c("forestgreen","red")) +
  theme_bw() +
  theme(legend.position="none",
        axis.text.x=element_text(angle = -45, vjust=1, hjust=0)) +
  ylab('Number of species') + xlab('Section')
  


ggarrange(g1, g2, nrow=2, labels=c('a)','b)'))


