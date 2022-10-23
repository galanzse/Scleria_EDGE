
library(tidyverse)
library(readxl)
source('scripts/red_list.R')


# GBIF
# citation: Scleria P.J.Bergius in GBIF Secretariat (2021). GBIF Backbone Taxonomy. Checklist dataset https://doi.org/10.15468/dl.su567h accessed via GBIF.org on 2022-10-14.
occurrences_gbif <- read.delim("data/occurrences_GBIF_141022.csv")
colnames(occurrences_gbif) <- tolower(colnames(occurrences_gbif))
head(occurrences_gbif)
str(occurrences_gbif)

# remove observations with (reported) uncertainty > 10km (enough for our scope)
occurrences_gbif$coordinateuncertaintyinmeters[is.na(occurrences_gbif$coordinateuncertaintyinmeters)] <- 1000
occurrences_gbif <- occurrences_gbif %>% filter(species!='' &
                                                  decimallatitude!='' & decimallongitude!='' &
                                                  decimallatitude!=0 & decimallongitude!=0 &
                                                  coordinateuncertaintyinmeters<10000)
occurrences_gbif <- occurrences_gbif %>% dplyr::select(gbifid, species, scientificname, decimallatitude, decimallongitude)
occurrences_gbif$gbifid <- as.character(occurrences_gbif$gbifid)
occurrences_gbif$dataset <- 'new_gbif'

# GBIF 2022 tiene 25428 observaciones mas que en 2018
# hay 16040 observaciones en GBIF 2018 que no estan en 2021 (ID diferente)
# (read.delim("C:/Users/javie/OneDrive/TESIS Y PUBLICACIONES/SCLERIA/Biogeography paper/DATA/Ocurrences/versiones antiguas/0016677-181003121212138.csv"))


# Hannah
occurrences_hannah <- read_excel("data/occurrences_hannah.xlsx")
colnames(occurrences_hannah) <- tolower(colnames(occurrences_hannah))
occurrences_hannah <- occurrences_hannah %>% dplyr::select(gbifid, species, scientificname, decimallatitude, decimallongitude)
occurrences_hannah$dataset <- 'hannah'


# previous database (some occurrences were mannually added by Isabel)
occurrences_old <- read.csv("data/occurrences_jofbiogeo.txt", sep="")
occurrences_old <- occurrences_old %>% dplyr::select(gbifid, species, scientificname, decimallatitude, decimallongitude)
occurrences_old$gbifid <- as.character(occurrences_old$gbifid)
occurrences_old$dataset <- 'old'

# retain observations not on the GBIF file
occurrences_old <- occurrences_old[!(occurrences_old$gbifid %in% occurrences_gbif$gbifid),]



# rbind
occ_data <- rbind(occurrences_gbif, occurrences_old, occurrences_hannah)
colSums(is.na(occ_data))/nrow(occ_data) * 100 # % NAs

# we will not need intraspecific epithets
occ_data$species[occ_data$species=='Scleria pergracilis brachystachys'] <- 'Scleria pergracilis'
occ_data$species[occ_data$species=='Scleria distans chondrocarpa'] <-  'Scleria distans'

# all assessed species have occurrence data
table(gsub("_"," ",scleria_iucn$species) %in% occ_data$species)


# some species are present in GBIF but not assessed
setdiff(occ_data$species, gsub("_"," ",scleria_iucn$species))

occ_data <- occ_data %>% subset(!(species %in% c('Scleria Maggieville', 'Scleria Laura','Scleria Oenpelli','Scleria McMinns-Lagoon','Scleria Jabiru','Scleria Wilton-River'))) # some entries are wrong

occ_data$species[occ_data$species=='Scleria sobolifera'] <- 'Scleria sobolifer' # check spelling and synonyms in WCSP
occ_data$species[occ_data$species=='Scleria sieberi'] <- 'Scleria gaertneri'
occ_data$species[occ_data$species=='Scleria melaleuca'] <- 'Scleria gaertneri'
occ_data$species[occ_data$species=='Scleria minima'] <- 'Scleria pusilla'
occ_data$species[occ_data$species=='Scleria pterota'] <- 'Scleria gaertneri'
occ_data$species[occ_data$species=='Scleria goosenii'] <- 'Scleria goossensii'
occ_data$species[occ_data$species=='Scleria abortiva'] <- 'Scleria trialata'
occ_data$species[occ_data$species=='Scleria parallella'] <- 'Scleria parallela'
occ_data$species[occ_data$species=='Scleria valdemuricata'] <- 'Scleria tenella'
occ_data$species[occ_data$species=='Scleria griegifolia'] <- 'Scleria greigiifolia'
occ_data$species[occ_data$species=='Scleria poaeoides'] <- 'Scleria pooides'
occ_data$species[occ_data$species=='Scleria stricta'] <- 'Scleria neesii'
occ_data$species[occ_data$species=='Scleria vichadensis'] <- 'Scleria mitis'
occ_data$species[occ_data$species=='Scleria pygmaea'] <- 'Diplacrum pygmaeum'


# filter only assessed species
occ_data <- occ_data %>% filter(species %in% gsub("_"," ",scleria_iucn$species))


rm(occurrences_gbif, occurrences_hannah, occurrences_old)

