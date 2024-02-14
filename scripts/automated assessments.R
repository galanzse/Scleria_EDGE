

# AUTOMATED ASSESSMENTS BASED ON ML METHODS (WALKER ET AL. 2023)
# WE WILL IMPUTE SPECIES WITH OCCURRENCE DATA AVAILABLE (N=67)

library(tidyverse); library(terra); library(readxl)
library(randomForest); library(caret)
# library(ConR)
# devtools::install_version("rCAT",version="0.1.6")
library(rCAT)

load('results/data_final.RData') # data
assessments <- data_final[['assessments']]



# rCAT ####
IUCN_results <- ConBatch(taxa=data_final[['occurrences']]$scientific_name,
                         lat=data_final[['occurrences']]$y,
                         long=data_final[['occurrences']]$x,
                         cellsize=2000,
                         project2gether=FALSE)
colnames(IUCN_results)[1] <- 'scientific_name'

# write.table(IUCN_results, 'results/IUCN_results.txt')
IUCN_results <- read.csv("results/IUCN_results.txt", sep="")


# RandomForest: traits, infrageneric data and coordinates ####
predictors <- data_final[['taxa']] %>%
  subset(scientific_name %in% data_final[['occurrences']]$scientific_name) %>%
  merge(data_final[['assessments']], by='scientific_name', all.x=T) %>% # add IUCN categories
  dplyr::select(scientific_name, subgenus, section, IUCN_category)

table(is.na(predictors)) # check


# import predictors
rst_predictors <- rast('C:/Users/user/Desktop/rst_predictors.tiff')
names(rst_predictors)[10] <- 'protected_areas'


# retrieve variables from ConR
predictors <- merge(predictors, IUCN_results[,c("scientific_name","EOOkm2","AOO2km")], all.x=T)


predictors$y_mean <- NA
predictors$elevation <- NA
predictors$HPD <-  NA
predictors$HFP <-  NA
predictors$AnnualMeanTemperature <- NA
predictors$MinTemperatureColdestMonth <- NA
predictors$TemperatureAnnualRange <- NA
predictors$AnnualPrecipitation <- NA
predictors$PrecipitationDriesttMonth <- NA
predictors$PrecipitationSeasonality <- NA
predictors$prop_protectedarea <- NA

temp <- predictors$scientific_name[is.na(predictors$HPD)]

for (i in 1:nrow(predictors)) {
  
  temp_occ <- data_final[['occurrences']] %>% subset(scientific_name==predictors$scientific_name[i]) %>%
    dplyr::select(y,x,scientific_name)
  temp_pts <- temp_occ %>% vect(geom=c('x','y'), 'epsg:4326')
  
  # take the mean value from a 5 km radius buffer around each record’s coordinates
  temp_bff <- terra::buffer(temp_pts, width=5000) 

  # latitude of range centroid
  predictors$y_mean[i] <- mean(temp_occ$y)
  
  # maximum elevation
  predictors$elevation[i] <- terra::extract(rst_predictors$elevation, temp_pts, ID=F, na.rm=TRUE) %>% colMeans(na.rm=T)
  
  # minimum human population density (HPD)
  predictors$HPD[i] <- terra::extract(rst_predictors$HPD, temp_bff) %>% group_by(ID) %>%
    summarise(HPD=min(HPD, na.rm=T)) %>% dplyr::select(HPD) %>% colMeans(na.rm=T)

  # human footprint
  predictors$HFP[i] <- terra::extract(rst_predictors$HFP, temp_bff) %>% group_by(ID) %>%
    summarise(HFP=mean(HFP, na.rm=T)) %>% dplyr::select(HFP) %>% colMeans(na.rm=T)
  
  # climatic variables
  predictors$AnnualMeanTemperature[i] <- terra::extract(rst_predictors$AnnualMeanTemperature, temp_pts, ID=F, na.rm=TRUE) %>% colMeans(na.rm=T)
  
  predictors$MinTemperatureColdestMonth[i] <- terra::extract(rst_predictors$MinTemperatureColdestMonth, temp_pts, ID=F, na.rm=TRUE) %>% colMeans(na.rm=T)
  
  predictors$TemperatureAnnualRange[i] <- terra::extract(rst_predictors$TemperatureAnnualRange, temp_pts, ID=F, na.rm=TRUE) %>% colMeans(na.rm=T)
  
  predictors$AnnualPrecipitation[i] <- terra::extract(rst_predictors$AnnualPrecipitation, temp_pts, ID=F, na.rm=TRUE) %>% colMeans(na.rm=T)
  
  predictors$PrecipitationDriesttMonth[i] <- terra::extract(rst_predictors$PrecipitationDriesttMonth, temp_pts, ID=F, na.rm=TRUE) %>% colMeans(na.rm=T)
  
  predictors$PrecipitationSeasonality[i] <- terra::extract(rst_predictors$PrecipitationSeasonality, temp_pts, ID=F, na.rm=TRUE) %>% colMeans(na.rm=T)
  
  pa_temp <- terra::extract(rst_predictors$protected_areas, temp_pts)
  predictors$prop_protectedarea[i] <- sum(as.numeric(pa_temp$protected_areas), na.rm=T)/nrow(pa_temp)  
  
  # 
  print(predictors$scientific_name[i])
  
}

# write.table(predictors, 'results/predictors.txt')
predictors <- read.csv('results/predictors.txt', sep=' ')
predictors$section <- NULL
predictors$subgenus <- NULL

colSums(is.na(predictors))


# add new variable 'threatened', see Walker et al. 2023
predictors$threatened <- NA
predictors$threatened[predictors$IUCN_category%in%c('LC','NT')] <- 'nonthreatened'
predictors$threatened[predictors$IUCN_category%in%c('CR','EN','VU')] <- 'threatened'

# exclude extinct species (S. chevaleri)
predictors <- predictors %>% subset(!(IUCN_category=='EX'))

# model
v_pred <- colnames(predictors)[-which(colnames(predictors)%in%c('scientific_name','IUCN_category','threatened'))]

# subset to be predicted
mod_out <- predictors %>% subset(is.na(threatened))
mod_out$threatened <- 'NE'
nrow(mod_out)

# subset to train and evaluate model
mod_in <- predictors %>% dplyr::select(all_of(v_pred), threatened) %>% subset(!is.na(threatened))

# training
v_ex <- sample(1:nrow(mod_in), nrow(mod_in)/5) # 20% of observations for external validation
data_ex <- mod_in[v_ex,]
data_tr <- mod_in[-v_ex,]

table(data_ex$threatened)
table(data_tr$threatened)

# 5-fold cross validation 
control <- trainControl(method='repeatedcv', number=5, repeats=10, search='grid')
rf_random <- train(threatened~ ., data=data_tr, method='rf', metric='Accuracy', trControl=control)

rf_random
plot(rf_random)
varImp(rf_random, scale=F)
plot(varImp(rf_random, scale=F))

# confusion matrix
data_ex$rf_out <- predict(rf_random, newdata=data_ex[,v_pred], type='raw')
confusionMatrix(data=data_ex$rf_out, reference=as.factor(data_ex$threatened))

# predict
mod_out$rf_threatened <- predict(rf_random, mod_out[,v_pred])
table(mod_out$rf_threatened)

# save
# write.table(mod_out[,c('scientific_name','rf_threatened')], 'results/rf_results.txt')
# save(rf_random, file="results/rf_random.RData")

rf_results <- read.csv("results/rf_results.txt", sep="")



# Model comparison ####


# RandomForest
table(data_ex$rf_out)
table(as.factor(data_ex$threatened), data_ex$rf_out)


# rCAT
comp_trends <- merge(assessments[,c('scientific_name','IUCN_category')],
                     IUCN_results[,c('scientific_name','EOOcat','AOOcat')], all.x=T) %>%
  merge(rf_results, all.x=T) %>%
  subset(IUCN_category!='EX' & !(is.na(EOOcat))) %>% #
  subset(!(scientific_name %in% mod_out$scientific_name))

comp_trends$IUCN_category <- as.character(comp_trends$IUCN_category)
comp_trends$IUCN_category[comp_trends$IUCN_category%in%c('VU','EN','CR')] <- 'threatened'
comp_trends$IUCN_category[comp_trends$IUCN_category%in%c('LC','NT')] <- 'nonthreatened'

comp_trends$EOOcat <- as.character(comp_trends$EOOcat)
comp_trends$EOOcat[comp_trends$EOOcat%in%c('VU','EN','CR')] <- 'threatened'
comp_trends$EOOcat[comp_trends$EOOcat%in%c('LC','NT')] <- 'nonthreatened'

# comp_trends$AOOcat <- as.character(comp_trends$AOOcat)
# comp_trends$AOOcat[comp_trends$AOOcat%in%c('VU','EN','CR')] <- 'threatened'
# comp_trends$AOOcat[comp_trends$AOOcat%in%c('LC','NT')] <- 'nonthreatened'


table(comp_trends$IUCN_category)
table(comp_trends$IUCN_category, comp_trends$EOOcat)
# table(comp_trends$IUCN_category, comp_trends$AOOcat)

comp_trends 


# ConR vs RF
comp_trends2 <- merge(rf_results[,c('scientific_name', 'rf_threatened')],
                      IUCN_results[,c('scientific_name','EOOkm2','AOO2km','EOOcat')],
                      by='scientific_name', all.y=T)

comp_trends2 <- merge(comp_trends2,
                      data_final[['assessments']][,c('scientific_name','IUCN_category')],
                      by='scientific_name', all.x=T)

comp_trends2$EOOcat <- as.character(comp_trends2$EOOcat)
comp_trends2$EOOcat[comp_trends2$EOOcat%in%c('VU','EN','CR')] <- 'threatened'
comp_trends2$EOOcat[comp_trends2$EOOcat%in%c('LC','NT')] <- 'nonthreatened'

# comp_trends2$AOOcat <- as.character(comp_trends2$AOOcat)
# comp_trends2$AOOcat[comp_trends2$AOOcat%in%c('VU','EN','CR')] <- 'threatened'
# comp_trends2$AOOcat[comp_trends2$AOOcat%in%c('LC','NT')] <- 'nonthreatened'

# table(comp_trends2$rf_threatened, comp_trends2$subgenus) # results per subgenus

comp_trends2$EOOcat <- as.character(comp_trends2$EOOcat)
comp_trends2$EOOcat[comp_trends2$EOOcat%in%c('VU','EN','CR')] <- 'threatened'
comp_trends2$EOOcat[comp_trends2$EOOcat%in%c('LC','NT')] <- 'nonthreatened'

# comp_trends2$AOOcat <- as.character(comp_trends2$AOOcat)
# comp_trends2$AOOcat[comp_trends2$AOOcat%in%c('VU','EN','CR')] <- 'threatened'
# comp_trends2$AOOcat[comp_trends2$AOOcat%in%c('LC','NT')] <- 'nonthreatened'

table(comp_trends2$EOOcat)
table(comp_trends2$EOOcat, comp_trends2$rf_threatened)

write.table(comp_trends2, 'results/aa_results.txt')

# RandomForest is a better model because it sets a greater EOO, AOO and n_locations threshold to consider species as threatened


# thresholds
comp_trends3 <- comp_trends2[!is.na(comp_trends2$rf_threatened),]
str(comp_trends3)
comp_trends3$EOOkm2 <- as.numeric(comp_trends3$EOOkm2)
comp_trends3$AOO2km <- as.numeric(comp_trends3$AOO2km)

mean(comp_trends3$EOOkm2[comp_trends3$rf_threatened=='threatened'])
mean(comp_trends3$AOO2km[comp_trends3$rf_threatened=='threatened'])
mean(comp_trends3$EOOkm2[comp_trends3$rf_threatened=='nonthreatened'])
mean(comp_trends3$AOO2km[comp_trends3$rf_threatened=='nonthreatened'])

mean(comp_trends3$EOOkm2[comp_trends3$EOOcat=='threatened'])
mean(comp_trends3$AOO2km[comp_trends3$EOOcat=='threatened'])
mean(comp_trends3$EOOkm2[comp_trends3$EOOcat=='nonthreatened'])
mean(comp_trends3$AOO2km[comp_trends3$EOOcat=='nonthreatened'])
