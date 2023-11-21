

# PREPARE A STACK OF PREDICTORS FOR AUTOMATED ASSESSMENTS
# work at 30s, it is not necesary to reproject but we will resample

library(tidyverse); library(terra)



# climatic variables: 30s
AnnualMeanTemperature <- rast('C:/Users/user/Desktop/worldclim/wc2.1_30s_bio/wc2.1_30s_bio_1.tif')
names(AnnualMeanTemperature) <- 'AnnualMeanTemperature'
MinTemperatureColdestMonth <- rast('C:/Users/user/Desktop/worldclim/wc2.1_30s_bio/wc2.1_30s_bio_6.tif')
names(MinTemperatureColdestMonth) <- 'MinTemperatureColdestMonth'
TemperatureAnnualRange <- rast('C:/Users/user/Desktop/worldclim/wc2.1_30s_bio/wc2.1_30s_bio_7.tif')
names(TemperatureAnnualRange) <- 'TemperatureAnnualRange'
AnnualPrecipitation <- rast('C:/Users/user/Desktop/worldclim/wc2.1_30s_bio/wc2.1_30s_bio_12.tif')
names(AnnualPrecipitation) <- 'AnnualPrecipitation'
PrecipitationDriesttMonth <- rast('C:/Users/user/Desktop/worldclim/wc2.1_30s_bio/wc2.1_30s_bio_14.tif')
names(PrecipitationDriesttMonth) <- 'PrecipitationDriesttMonth'
PrecipitationSeasonality <- rast('C:/Users/user/Desktop/worldclim/wc2.1_30s_bio/wc2.1_30s_bio_15.tif')
names(PrecipitationSeasonality) <- 'PrecipitationSeasonality'


# elevation
elevation <- rast('C:/Users/user/Desktop/worldclim/wc2.1_2.5m_elev/wc2.1_2.5m_elev.tif') %>%
  resample(AnnualMeanTemperature, method="bilinear")
names(elevation) <- 'elevation'


# human population density
HPD <- rast('C:/Users/user/Desktop/Global_2020_PopulationDensity30sec_GPWv4.tiff') %>%
  resample(AnnualMeanTemperature, method="bilinear")
names(HPD) <- 'HPD'


# human footprint
HFP <- rast('C:/Users/user/Desktop/hfp2013_merisINT.tif') %>%
  terra::project('epsg:4326')%>%
  resample(AnnualMeanTemperature, method="bilinear")
names(HFP) <- 'HFP'


# merge predictors in stack
rst_predictors <- c(elevation, HFP, HPD, AnnualMeanTemperature, MinTemperatureColdestMonth,  TemperatureAnnualRange, AnnualPrecipitation, PrecipitationDriesttMonth, PrecipitationSeasonality)

writeRaster(rst_predictors, 'C:/Users/user/Desktop/rst_predictors.tiff', overwrite=TRUE) # save



# merge all protected areas shp in a single map
setwd('C:/Users/user/Desktop/protected areas')

AF_protected <- rbind(vect('WDPA_WDOECM_Oct2023_Public_AF_shp/WDPA_WDOECM_Oct2023_Public_AF_shp_0/WDPA_WDOECM_Oct2023_Public_AF_shp-polygons.shp'),
                      vect('WDPA_WDOECM_Oct2023_Public_AF_shp/WDPA_WDOECM_Oct2023_Public_AF_shp_1/WDPA_WDOECM_Oct2023_Public_AF_shp-polygons.shp'),
                      vect('WDPA_WDOECM_Oct2023_Public_AF_shp/WDPA_WDOECM_Oct2023_Public_AF_shp_2/WDPA_WDOECM_Oct2023_Public_AF_shp-polygons.shp'))

AS_protected <- rbind(vect('WDPA_WDOECM_Oct2023_Public_AS_shp/WDPA_WDOECM_Oct2023_Public_AS_shp_0/WDPA_WDOECM_Oct2023_Public_AS_shp-polygons.shp'),
                      vect('WDPA_WDOECM_Oct2023_Public_AS_shp/WDPA_WDOECM_Oct2023_Public_AS_shp_1/WDPA_WDOECM_Oct2023_Public_AS_shp-polygons.shp'),
                      vect('WDPA_WDOECM_Oct2023_Public_AS_shp/WDPA_WDOECM_Oct2023_Public_AS_shp_2/WDPA_WDOECM_Oct2023_Public_AS_shp-polygons.shp'))

NA_protected <- rbind(vect('WDPA_WDOECM_Oct2023_Public_NA_shp/WDPA_WDOECM_Oct2023_Public_NA_shp_0/WDPA_WDOECM_Oct2023_Public_NA_shp-polygons.shp'),
                      vect('WDPA_WDOECM_Oct2023_Public_NA_shp/WDPA_WDOECM_Oct2023_Public_NA_shp_1/WDPA_WDOECM_Oct2023_Public_NA_shp-polygons.shp'),
                      vect('WDPA_WDOECM_Oct2023_Public_NA_shp/WDPA_WDOECM_Oct2023_Public_NA_shp_2/WDPA_WDOECM_Oct2023_Public_NA_shp-polygons.shp'))

SA_protected <- rbind(vect('WDPA_WDOECM_Oct2023_Public_SA_shp/WDPA_WDOECM_Oct2023_Public_SA_shp_0/WDPA_WDOECM_Oct2023_Public_SA_shp-polygons.shp'),
                      vect('WDPA_WDOECM_Oct2023_Public_SA_shp/WDPA_WDOECM_Oct2023_Public_SA_shp_1/WDPA_WDOECM_Oct2023_Public_SA_shp-polygons.shp'),
                      vect('WDPA_WDOECM_Oct2023_Public_SA_shp/WDPA_WDOECM_Oct2023_Public_SA_shp_2/WDPA_WDOECM_Oct2023_Public_SA_shp-polygons.shp'))

WA_protected <- rbind(vect('WDPA_WDOECM_Oct2023_Public_WA_shp/WDPA_WDOECM_Oct2023_Public_WA_shp_0/WDPA_WDOECM_Oct2023_Public_WA_shp-polygons.shp'),
                      vect('WDPA_WDOECM_Oct2023_Public_WA_shp/WDPA_WDOECM_Oct2023_Public_WA_shp_1/WDPA_WDOECM_Oct2023_Public_WA_shp-polygons.shp'),
                      vect('WDPA_WDOECM_Oct2023_Public_WA_shp/WDPA_WDOECM_Oct2023_Public_WA_shp_2/WDPA_WDOECM_Oct2023_Public_WA_shp-polygons.shp'))

protected_world <- rbind(AF_protected['PA_DEF'], AS_protected['PA_DEF'], NA_protected['PA_DEF'], SA_protected['PA_DEF'], WA_protected['PA_DEF'])
rm(AF_protected, AS_protected, NA_protected, SA_protected, WA_protected)

# if it is near a protected area lets consider it protected
protected_rast <- rasterize(x=protected_world, y=AnnualMeanTemperature, fun='mean') 
names(protected_rast) <- 'protected_areas'

# merge to predictors
rst_predictors <- c(rst_predictors, protected_rast)
writeRaster(rst_predictors, 'C:/Users/user/Desktop/rst_predictors.tiff', overwrite=TRUE) # save



setwd('C:/Users/user/OneDrive/PUBLICACIONES/SCLERIA/Scleria_EDGE')


# check
# library(rnaturalearth)
# esp_rnat <- ne_countries("large", country = "Madagascar", returnclass = "sf") %>% vect()
# protected_rast %>% crop(esp_rnat) %>% plot(); lines(esp_rnat)


