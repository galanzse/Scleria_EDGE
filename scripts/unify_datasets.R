

# UNIFY TRAITS, PHYLOGENY, RED_LIST AND OCCURRENCE DATASETS FOR POSTERIOR ANALYSES

library(tidyverse); library(readxl)
library(ape); library(phytools); library(pez)
# https://cran.r-project.org/web/packages/ape/vignettes/DrawingPhylogenies.pdf


# import valid species names
scleria_taxa <- read_excel("data/Kew_MNHNP_vouchers.xlsx", sheet = "traits_literature") %>%
  dplyr::select(subgenus, section, scientific_name)


# import traits
LHS_means_final <- read.csv("results/LHS_means_final.txt", sep="")


# import occurrences
occ_scleria_filt <- read.csv("results/occ_scleria_filt.txt", sep="")
colnames(occ_scleria_filt)[colnames(occ_scleria_filt)=='species'] <- 'scientific_name'


# import IUCN data: lets use Isabel's data as it includes assessments completed and waiting to be published
scleria_iucn <- read_excel("data/Scleria_Species_List_Status.xlsx") %>% dplyr::select(-Reference)
table(scleria_iucn$scientific_name %in% scleria_taxa$scientific_name)

# add taxonomic info and missing species
scleria_iucn <- merge(scleria_taxa[,c('scientific_name','subgenus','section')],
                      scleria_iucn, by='scientific_name', all.x=T)
scleria_iucn <- unique(scleria_iucn) # S. bulbifera is assessed under two synonyms
scleria_iucn$IUCN_category[is.na(scleria_iucn$IUCN_category)] <- 'NE'

# order categories
scleria_iucn$IUCN_category <- as.factor(scleria_iucn$IUCN_category)
scleria_iucn$IUCN_category <- factor(scleria_iucn$IUCN_category, levels=c('NE', 'DD', 'LC', 'NT', 'VU', 'EN', 'CR', 'EX'))

# exploratory
par(mfrow=c(1,1),mar=c(4,4,2,2))
bb <- barplot(table(scleria_iucn$IUCN_category), ylab='No. of species', ylim=c(0,170), cex.names=0.8, cex.axis=0.8,
              col=c("white","grey","forestgreen","darkolivegreen1","yellow","orange2","red","black"))
a <- as.numeric(table(scleria_iucn$IUCN_category))
text(bb,a+10,labels=a,cex=0.8)



# import phylogeny
scleria_tree <- read.tree('data/scleria_tree.txt'); is.ultrametric(scleria_tree)
phylogeny_tips <- read_excel("data/phylogeny_tips.xlsx")


# use scientific names with authority
scleria_tree$tip.label <- phylogeny_tips$scientific_name[order(match(phylogeny_tips$tip, scleria_tree$tip.label))]


# merge intraspecific taxa
dup_tip <- c('Scleria distans Poir.', 'Scleria pergracilis (Nees) Kunth', 'Scleria bulbifera Hochst. ex A.Rich.')
tip2 <- which(scleria_tree$tip.label %in% dup_tip)[c(1,3,5)]
scleria_tree <- drop.tip(scleria_tree, tip2)
# three species are not correctly located
tip2 <- which(scleria_tree$tip.label %in% c('Scleria rutenbergiana Boeckeler','Scleria williamsii Gross','Scleria gaertneri Raddi'))
scleria_tree <- drop.tip(scleria_tree, tip2)


# par(mar=c(1,1,1,1))
# plot(scleria_tree)



# summary data availability per species
scleria_taxa$trait <- NA
scleria_taxa$phylogeny <- NA
scleria_taxa$occurrence <- NA
scleria_taxa$assessment <- NA

for (i in 1:nrow(scleria_taxa)) {
  scleria_taxa$trait[i] <- scleria_taxa$scientific_name[i] %in% LHS_means_final$scientific_name
  scleria_taxa$phylogeny[i] <- scleria_taxa$scientific_name[i] %in% scleria_tree$tip.label
  scleria_taxa$occurrence[i] <- scleria_taxa$scientific_name[i] %in% occ_scleria_filt$scientific_name
  scleria_taxa$assessment[i] <- scleria_taxa$scientific_name[i] %in%  scleria_iucn$scientific_name[scleria_iucn$IUCN_category%in%c('LC','NT','VU','EN','CR','EX')]
}

# write.csv(scleria_taxa, 'results/scleria_taxa.csv')



# save elements in a list
data_final <- list()
data_final[['taxa']] <- scleria_taxa
data_final[['traits']] <- LHS_means_final[,c("scientific_name","life_form","life_form_simp","height","blade_area","nutlet_volume")]
data_final[['phylogeny']] <- scleria_tree
data_final[['assessments']] <- scleria_iucn[,c("scientific_name", "IUCN_category" )]
data_final[['occurrences']] <- occ_scleria_filt


save(data_final, file="results/data_final.RData")


