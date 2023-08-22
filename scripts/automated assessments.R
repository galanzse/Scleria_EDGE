

# AUTOMATED ASSESSMENTS BASED ON ML METHODS (WALKER ET AL. 2023)
# WE WILL IMPUTE SPECIES WITH OCCURRENCE DATA AVAILABLE (N=67)

library(tidyverse); library(terra); library(readxl)
library(randomForest); library(caret)
library(ConR)


load('results/data_final.RData') # data



# ConR: remove S. lithosperma to avoid errors ####
occ_iucn <- data_final[['occurrences']] %>% dplyr::select(y,x,scientific_name) %>%
  subset(!(scientific_name%in%c('Scleria lithosperma (L.) Sw.','Scleria polycarpa Boeckeler')))

IUCN.eval(DATA=occ_iucn, Cell_size_AOO=2, Cell_size_locations=10, Resol_sub_pop=5, SubPop=T,
          method_locations='fixed_grid', method.range='convex.hull', exclude.area=F,
          write_file_option='excel', write_results=T)
 
# impute S. lithosperma and S. polycarpa manually
IUCN_results <- read_excel('results/IUCN_results.xlsx') # results



# RandomForest: traits, infrageneric data and coordinates ####
predictors <- data_final[['taxa']] %>%
  subset(scientific_name %in% data_final[['occurrences']]$scientific_name) %>%
  merge(data_final[['assessments']], by=c('scientific_name','section','subgenus'), all.x=T) # add iucn categories

table(is.na(predictors)) # check


# # predictors
# elevation <- rast('C:/Users/javie/Desktop/world_rasters/wc2.1_2.5m_elev.tif') # elevation
# 
# # climatic variables
# AnnualMeanTemperature <- rast('C:/Users/javie/Desktop/world_rasters/wc2.1_2.5m/wc2.1_2.5m_bio_1.tif')
# MinTemperatureColdestMonth <- rast('C:/Users/javie/Desktop/world_rasters/wc2.1_2.5m/wc2.1_2.5m_bio_6.tif')
# TemperatureAnnualRange <- rast('C:/Users/javie/Desktop/world_rasters/wc2.1_2.5m/wc2.1_2.5m_bio_7.tif')
# AnnualPrecipitation <- rast('C:/Users/javie/Desktop/world_rasters/wc2.1_2.5m/wc2.1_2.5m_bio_12.tif')
# PrecipitationWettestMonth <- rast('C:/Users/javie/Desktop/world_rasters/wc2.1_2.5m/wc2.1_2.5m_bio_13.tif')
# PrecipitationSeasonality <- rast('C:/Users/javie/Desktop/world_rasters/wc2.1_2.5m/wc2.1_2.5m_bio_15.tif')
# 
# HPD <- rast('C:/Users/javie/Desktop/world_rasters/Global_2020_PopulationDensity30sec_GPWv4.tiff')# human population density
# names(HPD) <- 'HPD'
# HPD <- terra::resample(HPD, elevation, method='average')
# 
# # human footprint
# HFP <- rast('C:/Users/javie/Desktop/world_rasters/hfp2013_merisINT.tif') %>% terra::project('epsg:4326')
# names(HFP) <- 'HFP'
# HFP <- terra::resample(HFP, elevation, method='average')
# 
# rst_predictors <- c(elevation, HFP, HPD, AnnualMeanTemperature, MinTemperatureColdestMonth,  TemperatureAnnualRange, AnnualPrecipitation, PrecipitationWettestMonth, PrecipitationSeasonality)
# writeRaster(rst_predictors, 'results/rst_predictors.tiff')

rst_predictors <- rast('C:/Users/javie/OneDrive/PUBLICACIONES/SCLERIA/Scleria_EDGE/results/rst_predictors.tiff')
names(rst_predictors)

# retrieve from ConR
temp <- IUCN_results[,c("taxa","EOO","AOO","Nbe_loc")]
colnames(temp) <- c('scientific_name',"EOO","AOO","n_locations")
predictors <- merge(predictors, temp, all.x=T)


predictors$y_mean <- NA
predictors$elevation <- NA
predictors$HPD <-  NA
predictors$HFP <-  NA
predictors$AnnualMeanTemperature <- NA
predictors$MinTemperatureColdestMonth <- NA
predictors$TemperatureAnnualRange <- NA
predictors$AnnualPrecipitation <- NA
predictors$PrecipitationWettestMonth <- NA
predictors$PrecipitationSeasonality <- NA

for (i in 1:nrow(predictors)) {
  
  temp_occ <- data_final[['occurrences']] %>% subset(scientific_name==predictors$scientific_name[i]) %>%
    dplyr::select(y,x,scientific_name)
  temp_pts <- temp_occ %>% vect(geom=c('x','y'), 'epsg:4326')
  temp_bff <- terra::buffer(temp_pts, width=5000) # take the mean value from a 5 km radius buffer around each record’s coordinates

  # latitude of range centroid
  predictors$y_mean[i] <- mean(temp_occ$y)
  
  # maximum elevation
  predictors$elevation[i] <- terra::extract(rst_predictors$wc2.1_2.5m_elev, temp_bff, na.rm=TRUE) %>%
    group_by(ID) %>%
    summarise(elev=max(wc2.1_2.5m_elev)) %>% dplyr::select(elev) %>% colMeans(na.rm=T)
  
  # minimum human population density (HPD)
  predictors$HPD[i] <- terra::extract(rst_predictors$HPD, temp_bff) %>% group_by(ID) %>%
    summarise(HPD=min(HPD)) %>% dplyr::select(HPD) %>% colMeans(na.rm=T)

  # human footprint
  predictors$HFP[i] <- terra::extract(rst_predictors$HFP, temp_bff) %>% group_by(ID) %>%
    summarise(HFP=mean(HFP, na.rm=T)) %>% dplyr::select(HFP) %>% colMeans(na.rm=T)
  
  # climatic variables
  predictors$AnnualMeanTemperature[i] <- terra::extract(rst_predictors$wc2.1_2.5m_bio_1, temp_bff, na.rm=TRUE) %>%
    group_by(ID) %>%
    summarise(clim=mean(wc2.1_2.5m_bio_1, na.rm=T)) %>% dplyr::select(clim) %>% colMeans(na.rm=T)
  
  predictors$MinTemperatureColdestMonth[i] <- terra::extract(rst_predictors$wc2.1_2.5m_bio_6, temp_bff, na.rm=TRUE) %>%
    group_by(ID) %>%
    summarise(clim=mean(wc2.1_2.5m_bio_6, na.rm=T)) %>% dplyr::select(clim) %>% colMeans(na.rm=T)
  
  predictors$TemperatureAnnualRange[i] <- terra::extract(rst_predictors$wc2.1_2.5m_bio_7, temp_bff, na.rm=TRUE) %>%
    group_by(ID) %>%
    summarise(clim=mean(wc2.1_2.5m_bio_7, na.rm=T)) %>% dplyr::select(clim) %>% colMeans(na.rm=T)
  
  predictors$AnnualPrecipitation[i] <- terra::extract(rst_predictors$wc2.1_2.5m_bio_12, temp_bff, na.rm=TRUE) %>%
    group_by(ID) %>%
    summarise(clim=mean(wc2.1_2.5m_bio_12, na.rm=T)) %>% dplyr::select(clim) %>% colMeans(na.rm=T)
  
  predictors$PrecipitationWettestMonth[i] <- terra::extract(rst_predictors$wc2.1_2.5m_bio_13, temp_bff, na.rm=TRUE) %>%
    group_by(ID) %>%
    summarise(clim=mean(wc2.1_2.5m_bio_13, na.rm=T)) %>% dplyr::select(clim) %>% colMeans(na.rm=T)
  
  predictors$PrecipitationSeasonality[i] <- terra::extract(rst_predictors$wc2.1_2.5m_bio_15, temp_bff, na.rm=TRUE) %>%
    group_by(ID) %>%
    summarise(clim=mean(wc2.1_2.5m_bio_15, na.rm=T)) %>% dplyr::select(clim) %>% colMeans(na.rm=T)
  
  # 
  print(predictors$scientific_name[i])
}

# write.table(predictors, 'results/predictors.txt')
predictors <- read.csv('results/predictors.txt', sep=' ')
predictors$section <- NULL
predictors$subgenus <- NULL

colSums(is.na(predictors))

# impute max
predictors[predictors$scientific_name=='Scleria lithosperma (L.) Sw.',c('EOO','AOO','n_locations')] <- apply(predictors[,c('EOO','AOO','n_locations')],2,max,na.rm=TRUE)
predictors[predictors$scientific_name=='Scleria polycarpa Boeckeler',c('EOO','AOO','n_locations')] <- apply(predictors[,c('EOO','AOO','n_locations')],2,max,na.rm=TRUE)

# impute mins
predictors$EOO[is.na(predictors$EOO)] <- min(predictors$EOO, na.rm=T)
predictors$AOO[is.na(predictors$AOO)] <- min(predictors$AOO, na.rm=T)


# add new variable 'threatened', see Walker et al. 2023
predictors$threatened <- NA
predictors$threatened[predictors$IUCN_category%in%c('LC','NT')] <- 'nonthreatened'
predictors$threatened[predictors$IUCN_category%in%c('CR','EN','VU')] <- 'threatened'

# exclude extinct species (S. chevaleri)
predictors <- predictors %>% subset(!(IUCN_category=='EX'))

# model
v_pred <- colnames(predictors)[-which(colnames(predictors)%in%c('scientific_name','life_form','life_form_simp','subgenus','section','IUCN_category'))]

# subset to be predicted
mod_out <- predictors %>% subset(is.na(threatened))
mod_out$threatened <- NULL

# subset to train and evaluate model
mod_in <- predictors %>% dplyr::select(v_pred) %>% subset(!is.na(threatened))

# training
v_ex <- sample(1:nrow(mod_in), nrow(mod_in)/4) # 25% of observations for external validation
data_ex <- mod_in[v_ex,]
data_tr <- mod_in[-v_ex,]

table(data_ex$threatened); table(data_tr$threatened)

# 5-fold cross validation 
control <- trainControl(method='repeatedcv', number=5, repeats=10, search='grid')
rf_random <- train(threatened~., data=data_tr[,v_pred], method='rf', metric='Accuracy', trControl=control)
rf_random
plot(rf_random)
varImp(rf_random, scale=F)
plot(varImp(rf_random, scale=F))

# confusion matrix
data_ex$rf_out <- predict(rf_random, newdata=data_ex[,v_pred], type='raw')
confusionMatrix(data=data_ex$rf_out, reference=as.factor(data_ex$threatened))

# predict
mod_out$rf_threatened <- predict(rf_random, mod_out[,v_pred[-which(v_pred=='threatened')]])
table(mod_out$rf_threatened)

# save
write.table(mod_out[,c('scientific_name','rf_threatened')], 'results/rf_results.txt')
save(rf_random, file="results/rf_random.RData")



# Model comparison ####


# RandomForest
confusionMatrix(data=data_ex$rf_out, reference=as.factor(data_ex$threatened))


# ConR
comp_trends <- merge(assessments[,c('scientific_name','IUCN_category')],
                     IUCN_results[,c('scientific_name','Category_CriteriaB')],
                     all.x=T) %>% subset(!(IUCN_category%in%c('DD','NE','EX')) & !(is.na(Category_CriteriaB)))

comp_trends$IUCN_category <- as.character(comp_trends$IUCN_category)
comp_trends$Category_CriteriaB <- as.character(comp_trends$Category_CriteriaB)
comp_trends$IUCN_category[comp_trends$IUCN_category%in%c('VU','EN','CR')] <- 'threatened'
comp_trends$IUCN_category[comp_trends$IUCN_category%in%c('LC','NT')] <- 'nonthreatened'
comp_trends$Category_CriteriaB[comp_trends$Category_CriteriaB%in%c('VU','EN','CR')] <- 'threatened'
comp_trends$Category_CriteriaB[comp_trends$Category_CriteriaB%in%c('LC or NT')] <- 'nonthreatened'

table(comp_trends$IUCN_category)
table(comp_trends$IUCN_category, comp_trends$Category_CriteriaB)


# ConR vs RF
comp_trends2 <- merge(IUCN_results[,c('scientific_name','Category_CriteriaB')],
                     mod_out[,c('scientific_name', 'rf_threatened')], all.y=T)

comp_trends2$Category_CriteriaB <- as.character(comp_trends2$Category_CriteriaB)
comp_trends2$Category_CriteriaB[comp_trends2$Category_CriteriaB%in%c('VU','EN','CR')] <- 'threatened'
comp_trends2$Category_CriteriaB[comp_trends2$Category_CriteriaB%in%c('LC or NT')] <- 'nonthreatened'

table(comp_trends2$Category_CriteriaB)
table(comp_trends2$Category_CriteriaB, comp_trends2$rf_threatened)


# RandomForest is a much better model because it sets a greater EOO, AOO and n_locations threshold to consider species as threatened

