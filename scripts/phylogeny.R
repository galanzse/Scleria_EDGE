
source('scripts/red_list.R')
# library(tidyverse)
library(ape)
library(phytools)
library(pez)
# https://cran.r-project.org/web/packages/ape/vignettes/DrawingPhylogenies.pdf


# import data
scleria_tree <- read.tree('data/scleria_tree.txt')
is.ultrametric(scleria_tree)


# remove codes from tip.labels
scleria_tree$tip.label
tip2 <- str_split(scleria_tree$tip.label, pattern= "_", simplify = TRUE)[,1:2] %>% as.data.frame()
tip2 <- paste(tip2$V1, tip2$V2, sep='_')
scleria_tree$tip.label <- tip2
scleria_tree$tip.label <- gsub("_", " ", scleria_tree$tip.label) # change _ by spaces

# merge intraspecific taxa
tip2 <- which(scleria_tree$tip.label %in% c('Scleria distans','Scleria pergracilis'))[c(1,3)]
scleria_tree <- drop.tip(scleria_tree, tip2)
# three species are not correctly located
tip2 <- which(scleria_tree$tip.label %in% c('Scleria rutenbergiana','Scleria williamsii','Scleria gaertneri'))
scleria_tree <- drop.tip(scleria_tree, tip2)

# check taxonomy
scleria_tree$tip.label[scleria_tree$tip.label=='Scleria abortiva'] <- 'Scleria trialata'
scleria_tree$tip.label[scleria_tree$tip.label=='Scleria sobolifer'] <- 'Scleria sobolifera'

# 2 species in the tree do not have IUCN assessments
setdiff(scleria_tree$tip.label, scleria_iucn$species)
# I manually add these species to be considered in the analyses
not_assessed <- scleria_iucn[1:2,]
not_assessed$subgenus <- c('Hypoporum','Browniae')
not_assessed$section <- c('Hypoporum','Browniae')
not_assessed$species <- c("Scleria zambesica","Scleria depauperata")
not_assessed$in_tree <- 'Y'
not_assessed$category <- 'NE'
not_assessed$status <- 'not_assessed'
scleria_iucn <- rbind(scleria_iucn, not_assessed)


# 50 species need to be imputed in their sections
table(scleria_iucn$species %in% scleria_tree$tip.label)
scleria_iucn$species[which(!(scleria_iucn$species %in% scleria_tree$tip.label))]

# for this, I will substitute the genus by the section
scleria_iucn$sectxspp <- paste(scleria_iucn$section, str_split(scleria_iucn$species, pattern=" ", simplify = TRUE)[,2], sep=' ')

for (s in 1:length(scleria_tree$tip.label)) { 
  scleria_tree$tip.label[s] <- scleria_iucn$sectxspp[scleria_iucn$species==scleria_tree$tip.label[s]]
}

par(mar=c(1,1,1,1))
plot(scleria_tree, cex=0.5)

# vector of species to impute
v_sppximp <- scleria_iucn$sectxspp[which(!(scleria_iucn$sectxspp %in% scleria_tree$tip.label))]

# I impute the species nr times and save the results to assess the impact of imputation on the results later on

nr <- 100

imputed_trees <- list()
for (s in 1:nr) {
  imputed_trees[[s]] <- congeneric.impute(scleria_tree, species=v_sppximp, split=" ")
  print(s)
}

write.table(scleria_iucn, 'results/scleria_iucn.txt')
save(imputed_trees, file="results/imputed_trees.RData")

rm(tip2, not_assessed)


# plot an imputed tree
scleria_imputed <- imputed_trees[[1]]
scleria_imputed$tip.label <- scleria_iucn$species[order(match(scleria_iucn$sectxspp, scleria_imputed$tip.label))] # change labels

mycat <- scleria_iucn$category[order(match(scleria_iucn$species,scleria_imputed$tip.label))]
mycat <- factor(mycat)
mycol <- c("azure4","black","forestgreen","yellow","orange2","red")[mycat]

mycat2 <- scleria_iucn$in_tree[order(match(scleria_iucn$species,scleria_imputed$tip.label))]
mycat2 <- factor(mycat2)
mycex <- c(0.4,0.7)[mycat2]

par(mar=c(1,1,1,0))
plot(scleria_imputed, tip.color=mycol, cex=mycex)


