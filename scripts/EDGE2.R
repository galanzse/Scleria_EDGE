
library(tidyverse)
library(devtools)
library(ape); library(pez)


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

# S. depauperata could not be assessed because there was not occurrence data
phylogeny <- data_final[['phylogeny']] %>% drop.tip('Scleria depauperata Boeckeler')
# we are working with 239 species, 106 species need to be imputed in their sections
table(assessments$scientific_name %in% phylogeny$tip.label)


# treat sections as genera for imputation
assessments$sectxspp <- paste(assessments$section, str_split(assessments$scientific_name, pattern=" ", simplify = TRUE)[,2], sep=' ') 



# create 500 imputed trees using randtip and save in list ####

library(remotes)
remotes::install_github("iramosgutierrez/randtip")

phylogeny_sect <- phylogeny # change names in tree

intree <- assessments[assessments$scientific_name %in% phylogeny_sect$tip.label,]
phylogeny_sect$tip.label <- intree$sectxspp[order(match(intree$scientific_name, phylogeny_sect$tip.label))]

back.tree <- phylogeny_sect # backbone phylogeny
class(back.tree)
is.ultrametric(back.tree)

sp.list <- assessments$sectxspp # species list

# build info df, check and define ranks
my.info.noranks <- randtip::build_info(species=sp.list, tree=back.tree, mode='backbone', find.ranks=FALSE)
my.check <- randtip::check_info(my.info.noranks, tree=back.tree)
my.input.noranks <- randtip::info2input(my.info.noranks, back.tree)

imputed_trees <- list()
nr=500
assessments$sectxspp <- gsub(' ', '_', assessments$sectxspp) # to change to scientific names after imputation
for (i in 1:nr) {
  imp_tree <- randtip::rand_tip(input=my.input.noranks, tree=back.tree,
                                rand.type='random', # default
                                respect.mono=T, # all our clades are monophyletic so not relevant
                                prob=F, # branch selection probability is equiprobable
                                use.stem=T, # the stem branch can be considered as candidate for binding
                                prune=F, # we need the entire tree to run EDGE2
                                verbose=F)
  
  # (!!!) Nnode(imp_tree) == length(imp_tree$tip.label)

  imp_tree$tip.label <- assessments$scientific_name[order(match(assessments$sectxspp, imp_tree$tip.label))] # retrieve original names
  
  imputed_trees[[i]] <- imp_tree # save in list
  ck_n <- deframe(table(imp_tree$tip.label%in%assessments$scientific_name))
  print(paste('x',i,'   ', 'Ntip=', ck_n, sep='')) # progress
}

save(imputed_trees, file="results/imputed_trees.RData")


# plot an imputed tree
scleria_imputed <- imputed_trees[[16]]
assessments$threatened <- 'threathened'
assessments$threatened[assessments$IUCN_category%in%c('LC','NT','nonthreatened')] <- 'nonthreatened'
mycat <- assessments$threatened[order(match(assessments$scientific_name,scleria_imputed$tip.label))]
mycat <- factor(mycat)
mycol <- c("forestgreen","red")[mycat]

# par(mar=c(1,1,1,0))
# plot(scleria_imputed, tip.color=mycol, cex=0.7)



# probability of extinction ####

# load functions from github: https://github.com/rgumbs/EDGE2/
# source_url(url='https://raw.githubusercontent.com/rgumbs/EDGE2/main/GE.2.calc')

# prob_pext <- GE.2.calc(pext.vals)
# str(prob_pext)

# # add threatened/nonthreatened categories
# temp_threat <- prob_pext %>% subset(RL.cat %in% c("CR", "EN", "VU"))
# temp_threat$RL.cat <- 'threatened'
# temp_nonthreat <- prob_pext %>% subset(!(RL.cat %in% c("CR", "EN", "VU")))
# temp_nonthreat$RL.cat <- 'nonthreatened'
# prob_pext <- rbind(prob_pext, temp_threat, temp_nonthreat)
# write.table(prob_pext, 'results/prob_pext.txt')



# run EDGE2 ####


# load EDGE2 components
prob_pext <- read.csv("results/prob_pext.txt", sep="")
load("results/imputed_trees.RData")


# species x pext dataframe to fill in each interaction
spp_pext <- matrix(nrow=nrow(assessments), ncol=nr)
rownames(spp_pext) <- assessments$scientific_name
for (s in 1:nrow(spp_pext)) { # assign probability to each species x nr times
  temp <- assessments %>% subset(scientific_name==rownames(spp_pext)[s])
  spp_pext[s,] <- prob_pext %>% subset(RL.cat==temp$IUCN_category) %>% dplyr::select(pext) %>% deframe() %>% sample(500, replace=T)
  print(paste('---', round(s/nrow(spp_pext),2)*100, '% ---'))
}


# load functions from github: https://github.com/rgumbs/EDGE2/
# source_url(url='https://raw.githubusercontent.com/rgumbs/EDGE2/main/EDGE.2.calc')
source('scripts/EDGE.2.calc.R') # corrected function

# run EDGE nr times
l_EDGE <- list()
for (r in 1:nr) {
  GE2_temp <- data.frame(rownames(spp_pext), spp_pext[,r])
  colnames(GE2_temp) <- c('species','pext'); rownames(GE2_temp) <- NULL
  l_EDGE[[r]] <- EDGE2_mod(imputed_trees[[r]], GE2_temp, type='edge')
  print(paste('--- ', r/nr*100, '% ---', sep=''))
}
save(l_EDGE, file="results/l_EDGE.RData")



# explore results #### 


load("results/l_EDGE.RData")


# summary PD and ePDloss
v_PD <- vector()
v_ePDloss <- vector()
for (r in 1:length(l_EDGE)) {
  v_PD[r] <- l_EDGE[[r]][[3]][1,'PD']
  v_ePDloss[r] <- l_EDGE[[r]][[3]][1,'ePDloss']
}
mean(v_PD); sd(v_PD) # PD in MY for Scleria:
mean(v_ePDloss); sd(v_ePDloss) # Expected PD loss in MY for Scleria


# extract EDGE2 values in a df to explore distributions
EDGE2_values <- matrix(nrow=nrow(l_EDGE[[1]][[1]]), ncol=nr)
rownames(EDGE2_values) <- l_EDGE[[1]][[1]]$Species
ED2_values <- EDGE2_values
for (r in 1:nr) {
  edg1 <- l_EDGE[[r]][[1]] %>% dplyr::select(Species,ED,EDGE)
  edg1 <- edg1[order(match(edg1$Species,rownames(EDGE2_values))),]
  EDGE2_values[,r] <- edg1$EDGE
  ED2_values[,r] <- edg1$ED
}


# extract reference values per group
EDGE2_median <- matrix(nrow=length(unique(assessments$section)), ncol=5)
rownames(EDGE2_median) <- unique(assessments$section)
ED2_median <- EDGE2_median
for (s in 1:nrow(EDGE2_median)) {
  sp1 <- assessments$scientific_name[assessments$section==rownames(EDGE2_median)[s]]
  sp1 <- sp1[!is.na(sp1)]
  sp2 <- EDGE2_values[rownames(EDGE2_values)%in%sp1,] %>% as.vector() %>% unlist()
  EDGE2_median[s,] <- boxplot.stats(as.vector(sp2))$stats
  sp2 <- ED2_values[rownames(ED2_values)%in%sp1,] %>% as.vector() %>% unlist()
  ED2_median[s,] <- boxplot.stats(as.vector(sp2))$stats
}


par(mfrow=c(2,1),mar=c(7,4,1,1))
boxplot(t(na.omit(EDGE2_median)), col='white', las=2, cex.axis=0.8, ylab='EDGE2 score')
abline(h=median(EDGE2_values, na.rm=T), col='red')


# EDGE2 lists
source('scripts/EDGE.2.lists.R')

scleria_list <- EDGE2_lists(EDGE2_values, ED2_values, assessments[,c('scientific_name','thr')])

scleria_list <- scleria_list %>%
  merge(assessments[,c('scientific_name','section','subgenus')], all.x=T)

# write.csv(scleria_list, 'results/EDGE2_list.csv')

table(EDGE2_list$list)
