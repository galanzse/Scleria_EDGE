

# Ecologically Distinct and Globally Endangered (EcoDGE) Sclerias

library(tidyverse)
library(devtools)
library(readxl)
library(cluster)
library(ape) # library(pez)


# data
load("results/data_final.RData")
names(data_final)



# data: add RandomForest predictions to IUCN assessments ####

assessments <- data_final[['assessments']] %>% subset(!(is.na(section))) %>%
  subset(scientific_name!='Scleria chevalieri J.Raynal') # remove extinct species

rf_results <- read.csv("results/rf_results.txt", sep="") # import predictions

assessments <- merge(assessments, rf_results, all.x=T) # merge

assessments$IUCN_category <- as.character(assessments$IUCN_category) # change to character to impute
assessments$IUCN_category[assessments$IUCN_category %in% c('NE','DD')] <- assessments$rf_threatened[assessments$IUCN_category %in% c('NE','DD')]

# maintain second variable for lists
assessments$rf_threatened[assessments$IUCN_category%in%c('VU','EN','CR')] <- 'threatened'
assessments$rf_threatened[assessments$IUCN_category%in%c('LC','NT')] <- 'nonthreatened'
colnames(assessments)[which(colnames(assessments)=='rf_threatened')] <- 'thr'

assessments <- na.omit(assessments) # remove species we could not assess

# we are working with 239 species
traits <- data_final[['traits']] %>% subset(scientific_name %in% assessments$scientific_name)

# impute temporarily (I am not sure whether I measured traits on scabriuscula or schiedenana)
assessments$scientific_name[which(!(assessments$scientific_name %in% traits$scientific_name))]

# merge
traits <- merge(assessments, traits, by=c('scientific_name', 'subgenus', 'section'))

# impute by section x life form means
for (i in 1:nrow(traits)) {
  for (t in c('life_form_simp', 'height', 'blade_area', 'nutlet_volume')) {
    if (is.na(traits[i,t])) {
      temp <- traits %>% subset(section==traits$section[i] & life_form_simp==traits$life_form_simp[i]) %>%
        dplyr::select(t) %>% colMeans(na.rm=T)
      traits[i,t] <- temp
    }
  }
}

# S. khasiana is imputed from section data only
traits[traits$scientific_name=='Scleria khasiana Boeckeler',c('height', 'blade_area', 'nutlet_volume')] <- traits %>% subset(section=='elatae') %>% dplyr::select(height, blade_area, nutlet_volume) %>% colMeans(na.rm=T)


# explore functional differences between threaten categories #
table(traits$thr, traits$life_form_simp)

par(mfrow=c(1,3))
boxplot(log(height) ~ thr, data=traits)
boxplot(log(blade_area) ~ thr, data=traits)
boxplot(log(nutlet_volume) ~ thr, data=traits)



# create 15 dendrograms with LHS traits + life form: see Griffith et al 2022 (Functional Ecology) ####

v_traits <- c('life_form_simp', 'height', 'blade_area', 'nutlet_volume')
s_traits <- traits %>% dplyr::select(all_of(v_traits))
rownames(s_traits) <- traits$scientific_name
s_traits$life_form_simp <- as.factor(s_traits$life_form_simp)
str(s_traits)

# scale continuous variables
s_traits[,c('height', 'blade_area', 'nutlet_volume')] <- scale(s_traits[,c('height', 'blade_area', 'nutlet_volume')], center=T)

# distance matrix + clustering UPGMA
dend_combn <- read_excel("results/final_lists.xlsx", sheet = "dendrogram_combn")

l_dendr <- list()

for (r in 1:nrow(dend_combn)) {
  
  # select a particular combination of traits
  tr_cbmn <- dend_combn[r,] %>% t() %>% c()
  tr_cbmn <- tr_cbmn[!is.na(tr_cbmn)]
  
  # calculate dendrogram
  daisy.mat <- daisy(s_traits[,tr_cbmn], metric="gower") # %>% as.dist()
  # type=list('factor'='life_form_simp','numeric'=c('height','blade_area','nutlet_volume'))
  
  dend1 <- hclust(daisy.mat, method="average") # UPGMA
  l_dendr[[r]] <- as.phylo(dend1) # as.phylo
  
}



# Are phenotypes captured by the dendrogram?
dend1.phy <- l_EcoDGE[[11]] # dendrogram of 4 traits
mycat <- traits$subgenus[order(match(traits$scientific_name, dend1.phy$tip.label))]
mycat <- factor(mycat)
mycol <- c("brown","orange","forestgreen","blue")[mycat]

# par(mar=c(1,1,1,0))
# plot(as.phylo(dend1.phy), tip.color=mycol, cex=0.7)


# Are endangered species functionally similar?
mycat <- traits$thr[order(match(traits$thr,dend1.phy$tip.label))]
mycat <- factor(mycat)
mycol <- c("forestgreen","red")[mycat]

# par(mar=c(1,1,1,0))
# plot(as.phylo(dend1.phy), tip.color=mycol, cex=0.7)



# run EcoDGE2 ####


class(l_dendr[[1]])

# load functions from github: https://github.com/rgumbs/EDGE2/
# source_url(url='https://raw.githubusercontent.com/rgumbs/EDGE2/main/EDGE.2.calc')
# source_url(url='https://raw.githubusercontent.com/rgumbs/EDGE2/main/GE.2.calc')
source('scripts/EDGE.2.calc.R') # corrected function
prob_pext <- read.csv("results/prob_pext.txt", sep="")


# species x pext dataframe to fill in each interaction
spp_pext <- matrix(nrow=nrow(assessments), ncol=nr)
rownames(spp_pext) <- assessments$scientific_name
for (s in 1:nrow(spp_pext)) { # assign probability to each species x nr times
  temp <- assessments %>% subset(scientific_name==rownames(spp_pext)[s])
  spp_pext[s,] <- prob_pext %>% subset(RL.cat==temp$IUCN_category) %>% dplyr::select(pext) %>% deframe() %>% sample(500, replace=T)
  print(paste('---', round(s/nrow(spp_pext),2)*100, '% ---'))
}


# run EDGE nr times
l_EcoDGE <- list()
nr=500
for (r in 1:nr) {
  
  GE2_temp <- data.frame(rownames(spp_pext), spp_pext[,r])
  colnames(GE2_temp) <- c('species','pext'); rownames(GE2_temp) <- NULL
  
  dend1.phy <- l_dendr[[sample(1:length(l_dendr), 1)]] # randomly select a dendrogram
  
  l_EcoDGE[[r]] <- EDGE2_mod(dend1.phy, GE2_temp, type='ecodge')
  print(r)
}

# multiply FUD scores by 100 to balance the weighting of the two component values of the EcoDGE metric
# boxplot(l_EcoDGE[[1]][[1]]$pext, l_EcoDGE[[1]][[1]]$ED)
for (i in 1:length(l_EcoDGE)) {
  l_EcoDGE[[i]][[1]]$ED <- l_EcoDGE[[i]][[1]]$ED*100
  l_EcoDGE[[i]][[1]]$EDGE <- l_EcoDGE[[i]][[1]]$pext * l_EcoDGE[[i]][[1]]$ED
}

save(l_EcoDGE, file="results/l_EcoDGE.RData")



# explore results #### 


load("results/l_EcoDGE.RData")


# extract EDGE2 values in a df to explore distributions
EcoDGE2_values <- matrix(nrow=nrow(l_EcoDGE[[1]][[1]]), ncol=nr)
rownames(EcoDGE2_values) <- l_EDGE[[1]][[1]]$Species
FUD_values <- EcoDGE2_values
for (r in 1:nr) {
  edg1 <- l_EcoDGE[[r]][[1]] %>% dplyr::select(Species,ED,EDGE)
  edg1 <- edg1[order(match(edg1$Species,rownames(EcoDGE2_values))),]
  EcoDGE2_values[,r] <- edg1$EDGE
  FUD_values[,r] <- edg1$ED
}


# extract reference values per group
EcoDGE2_median <- matrix(nrow=length(unique(assessments$section)), ncol=5)
rownames(EcoDGE2_median) <- unique(assessments$section)
ED2_median <- EcoDGE2_median
for (s in 1:nrow(EcoDGE2_median)) {
  sp1 <- assessments$scientific_name[assessments$section==rownames(EcoDGE2_median)[s]]
  sp1 <- sp1[!is.na(sp1)]
  sp2 <- EcoDGE2_values[rownames(EcoDGE2_values)%in%sp1,] %>% as.vector() %>% unlist()
  EcoDGE2_median[s,] <- boxplot.stats(as.vector(sp2))$stats
  sp2 <- FUD_values[rownames(FUD_values)%in%sp1,] %>% as.vector() %>% unlist()
  ED2_median[s,] <- boxplot.stats(as.vector(sp2))$stats
}


par(mfrow=c(1,1),mar=c(7,4,1,1))
boxplot(t(na.omit(EcoDGE2_median)), col='white', las=2, cex.axis=0.8, ylab='EcoDGE2 score')
abline(h=median(EcoDGE2_values, na.rm=T), col='red')


# EDGE2 lists
source('scripts/EDGE.2.lists.R')

scleria_list <- EDGE2_lists(EcoDGE2_values, FUD_values, assessments[,c('scientific_name','thr')])

scleria_list <- scleria_list %>%
  merge(assessments[,c('scientific_name','section','subgenus')], all.x=T)

write.csv(scleria_list, 'results/EcoDGE2_list.csv')


