

# UNIFY TRAITS, PHYLOGENY, RED_LIST AND OCCURRENCE DATASETS FOR POSTERIOR ANALYSES

library(tidyverse); library(readxl)

library(ape); library(phytools); library(pez)
# https://cran.r-project.org/web/packages/ape/vignettes/DrawingPhylogenies.pdf



# import traits
LHS_means_final <- read.csv("C:/Users/user/OneDrive/PUBLICACIONES/SCLERIA/Scleria_EDGE/results/LHS_means_final.txt", sep="")

# import occurrences
occ_scleria_filt <- read.csv("C:/Users/user/OneDrive/PUBLICACIONES/SCLERIA/Scleria_EDGE/results/occ_scleria_filt.txt", sep="") %>%
  dplyr::select(species, x, y)

# We will work with species with known locations
table(occ_scleria_filt$species %in% LHS_means_final$scientific_name)
LHS_means_final <- LHS_means_final %>% subset(scientific_name %in% occ_scleria_filt$species)



# import IUCN data
# lets use Isabel's data as it includes assessments completed and waiting to be published
scleria_iucn <- read_excel("data/Scleria_Species_List_Status.xlsx") %>% dplyr::select(-Reference)
table(scleria_iucn$scientific_name %in% LHS_means_final$scientific_name)

# add taxonomic info and missing species
scleria_iucn <- merge(LHS_means_final[,c('scientific_name','subgenus','section')], scleria_iucn, by='scientific_name', all.x=T)
scleria_iucn <- unique(scleria_iucn) # S. bulbifera is assessed under two synonyms
scleria_iucn$IUCN_category[is.na(scleria_iucn$IUCN_category)] <- 'NE'

# order categories
scleria_iucn$IUCN_category <- as.factor(scleria_iucn$IUCN_category)
scleria_iucn$IUCN_category <- factor(scleria_iucn$IUCN_category, levels=c('NE', 'DD', 'LC', 'NT', 'VU', 'EN', 'CR', 'EX'))

# exploratory
par(mfrow=c(1,1),mar=c(4,4,2,2))
bb <- barplot(table(scleria_iucn$IUCN_category), ylab='No. of species', ylim=c(0,170), cex.names=0.8, cex.axis=0.8)
a <- as.numeric(table(scleria_iucn$IUCN_category))
text(bb,a+10,labels=a,cex=0.8)



# import phylogeny
scleria_tree <- read.tree('data/scleria_tree.txt'); is.ultrametric(scleria_tree)
phylogeny_tips <- read_excel("data/phylogeny_tips.xlsx")


# use scientific names with authority
scleria_tree$tip.label <- phylogeny_tips$scientific_name[order(match(phylogeny_tips$tip, scleria_tree$tip.label))]


# merge intraspecific taxa
tip2 <- which(scleria_tree$tip.label %in% c('Scleria distans Poir.','Scleria pergracilis (Nees) Kunth'))[c(1,3)]
scleria_tree <- drop.tip(scleria_tree, tip2)
# three species are not correctly located
tip2 <- which(scleria_tree$tip.label %in% c('Scleria rutenbergiana Boeckeler','Scleria williamsii Gross','Scleria gaertneri Raddi'))
scleria_tree <- drop.tip(scleria_tree, tip2)


par(mar=c(1,1,1,1))
plot(scleria_tree)



# save elements in a list
data_final <- list()
data_final[['traits']] <- LHS_means_final
data_final[['phylogeny']] <- scleria_tree
data_final[['assessments']] <- scleria_iucn
data_final[['occurrences']] <- occ_scleria_filt


save(data_final, file="results/data_final.RData")


