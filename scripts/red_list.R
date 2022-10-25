
library(tidyverse)
library(readxl)

# import data
scleria_iucn <- read_excel("data/Scleria_Species_List_Status.xlsx")
colnames(scleria_iucn) <- c('subgenus', 'section', 'species',
                            'in_tree', 'GBIF', 'category',
                            'reference')
scleria_iucn$in_tree[is.na(scleria_iucn$in_tree)] <- 'N'


# retain assessed species (published + submitted + about to submit)
v_categories <- c("LC","EN","DD","VU","NT","CR","EX")
scleria_iucn <- scleria_iucn %>% filter(category %in% v_categories)

# create variable with status
scleria_iucn$status <- 'published'
scleria_iucn$status[scleria_iucn$reference=='IN PREP'] <- 'in_prep'
scleria_iucn$status[scleria_iucn$reference=='SUBMITTED TO IUCN'] <- 'submitted'
scleria_iucn$status[scleria_iucn$reference=='TO BE SUBMITTED TO IUCN IN SEPT 2022'] <- 'to_be_submitted'
scleria_iucn$status[scleria_iucn$reference=='TO BE SUBMITTED TO IUCN IN SEPTEMBER 2022'] <- 'to_be_submitted'

# remove variables I wont need
scleria_iucn$reference <- NULL
scleria_iucn$GBIF <- NULL

# detect no-break spaces
scleria_iucn$species
scleria_iucn$species <- gsub("\u00A0", " ", scleria_iucn$species, fixed=TRUE)
tip2 <- str_split(scleria_iucn$species, pattern= " ", simplify = TRUE)[,1:2] %>% as.data.frame()
tip2 <- paste(tip2$V1, tip2$V2, sep=' ')
scleria_iucn$species <- tip2

# check names
scleria_iucn$section[scleria_iucn$species=='Scleria polyrrhiza'] <- 'Hypoporum'
scleria_iucn$section[scleria_iucn$species=='Scleria variegata'] <- 'Virgatae' # following Bauters et al 2019 (PlosOne)

# exploratory
par(mfrow=c(1,1),mar=c(4,4,2,2))
scleria_iucn$category <- as.factor(scleria_iucn$category)
scleria_iucn$category <- factor(scleria_iucn$category, levels = c("DD","LC","NT","VU","EN","CR","EX"))
bb <- barplot(table(scleria_iucn$category), ylab='No. of species', ylim=c(0,170), cex.names=0.8, cex.axis=0.8)
a <- as.numeric(table(scleria_iucn$category))
text(bb,a+10,labels=a,cex=0.8)

barplot(table(scleria_iucn$status))
table(scleria_iucn$in_tree)

# eliminate S. chevaleri for posterior analyses
scleria_iucn <- scleria_iucn %>% filter(category != 'EX')


rm(tip2, v_categories, a, bb)
