

# EXPLORE CONSERVATION STATUS OF ECOREGIONS WHERE SCLERIA IS PRESENT FOLLOWING DINERSTEIN ET AL. 2017


library(tidyverse)
library(terra)


load("results/data_final.RData") # data


# ecoregions
ecoregions <- vect('C:/Users/user/Desktop/Ecoregions2017/Ecoregions2017.shp')
rst_ecoregions <- rast('C:/Users/user/Desktop/worldclim/wc2.1_30s_bio/wc2.1_30s_bio_1.tif')

# rasterize
rst_ecoregions <- rasterize(x=ecoregions, y=rst_ecoregions, field='NNH')
writeRaster(rst_ecoregions, 'C:/Users/user/Desktop/Ecoregions2017/rst_ecoregions.tif')

# points
pts_scleria <- data_final[['occurrences']][,c('scientific_name','x','y')] %>% vect(geom=c('x','y'), 'epsg:4326')


# extract values
df_ecoregions <- pts_scleria %>% terra::extract(x=rst_ecoregions)
df_ecoregions$scientific_name <- pts_scleria$scientific_name
df_ecoregions$ID <- NULL
head(df_ecoregions)


# boxplot
par(mfrow=c(1,1),mar=c(4,4,2,2))
bb <- barplot(table(df_ecoregions$NNH), ylab='No. of observations', ylim=c(0,10000), cex.names=0.8, cex.axis=0.8,
              col=c("forestgreen","darkolivegreen1","orange2","red"))
a <- round(as.numeric(table(df_ecoregions$NNH)/sum(table(df_ecoregions$NNH))),2)
text(bb,a,labels=a, pos=3, cex=0.8)


uq_df_ecoregions <- unique(df_ecoregions) %>% na.omit()
uq_df_ecoregions$values <- 1
NNH_scleria <- pivot_wider(uq_df_ecoregions, names_from='NNH', values_from='values')
write.table(NNH_scleria, 'results/NNH_scleria.txt')


