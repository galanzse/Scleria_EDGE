
library(tidyverse)
library(devtools)
source('scripts/functions_github.R')


# load functions from github: https://github.com/rgumbs/EDGE2/
source_url(url='https://raw.githubusercontent.com/rgumbs/EDGE2/main/EDGE.2.calc')
source_url(url='https://raw.githubusercontent.com/rgumbs/EDGE2/main/GE.2.calc')


# data
load("results/imputed_trees.RData")
scleria_iucn <- read.csv("C:/Users/user/OneDrive/TESIS Y PUBLICACIONES/SCLERIA/Scleria_EDGE/results/scleria_iucn.txt", sep="")


# probability of extinction dataframe to withdraw propabilities from
prob_pext <- GE.2.calc(pext)
str(prob_pext)


# species x pext dataframe to fill in each interaction
spp_pext <- scleria_iucn[,c('sectxspp','category')]
colnames(spp_pext) <- c('species','GE2')

# assign probability to each species x nr times
nr <- 100

l_spp_pext <- list()
for (r in 1:nr) {
  for (s in 1:nrow(spp_pext)) {
    pspp <- spp_pext$species[s]
    pcat <- scleria_iucn$category[scleria_iucn$sectxspp==pspp]
    if (pcat %in% c("LC","EN","VU","NT","CR")) {
      spp_pext$GE2[spp_pext$species==pspp] <- sample(prob_pext$pext[prob_pext$RL.cat==pcat], 1)
    } else {
      spp_pext$GE2[spp_pext$species==pspp] <- sample(prob_pext$pext, 1)
    }
  }
  spp_pext$GE2 <- as.numeric(spp_pext$GE2)
  l_spp_pext[[r]] <- spp_pext
  print(r)
}
save(l_spp_pext, file="results/l_spp_pext.RData")


# finally, I run EDGE nr times
l_EDGE <- list()
for (s in 1:nr) {
  l_EDGE[[s]] <- EDGE.2.calc(imputed_trees[[s]], l_spp_pext[[s]])
  print(s)
}
save(l_EDGE, file="results/l_EDGE.RData")


# 
v_PD <- vector()
v_ePDloss <- vector()
for (s in 1:length(l_EDGE)) {
  v_PD[s] <- l_EDGE[[s]][3][[1]][1,'PD']
  v_ePDloss[s] <- l_EDGE[[s]][3][[1]][1,'ePDloss']
}
mean(v_PD); sd(v_PD)
mean(v_ePDloss); sd(v_ePDloss)


# extract EDGE2 values in a df to explore distributions
EDGE2_values <- matrix(nrow=nrow(scleria_iucn), ncol=nr) %>% as.data.frame()
rownames(EDGE2_values) <- scleria_iucn$sectxspp
ED2_values <- EDGE2_values
for (s in 1:nr) {
  edg1 <- l_EDGE[[s]][[1]] %>% as.data.frame() %>% dplyr::select(Species,ED,EDGE); colnames(edg1)[1] <- c('sectxspp')
  edg1 <- merge(scleria_iucn[,c('sectxspp','species')], edg1, by='sectxspp', all.x=T)
  edg1 <- edg1[order(match(edg1$sectxspp,rownames(EDGE2_values))),]
  EDGE2_values[,s] <- edg1$EDGE
  ED2_values[,s] <- edg1$ED
}


# extract reference values per group
EDGE2_median <- matrix(nrow=length(unique(scleria_iucn$section)), ncol=5) %>% as.data.frame()
rownames(EDGE2_median) <- unique(scleria_iucn$section)
ED2_median <- EDGE2_median
for (s in 1:nrow(EDGE2_median)) {
  sp1 <- scleria_iucn$sectxspp[scleria_iucn$section==rownames(EDGE2_median)[s]]
  sp2 <- EDGE2_values[sp1,] %>% as.vector() %>% unlist()
  EDGE2_median[s,] <- boxplot.stats(as.vector(sp2))$stats
  sp2 <- ED2_values[sp1,] %>% as.vector() %>% unlist()
  ED2_median[s,] <- boxplot.stats(as.vector(sp2))$stats
}


par(mfrow=c(1,1),mar=c(7,4,1,1))
boxplot(t(EDGE2_median), col='white', las=2, cex.axis=0.8, ylab='EDGE2 score')
abline(h=0.116972, col='red')


# EDGE2 lists
EDGE2_list <- EDGE2_values[,1:3]
colnames(EDGE2_list) <- c('EDGE2','perc','list')
EDGE2_list$sectxspp <- rownames(EDGE2_list); rownames(EDGE2_list) <- NULL
EDGE2_list <- merge(scleria_iucn, EDGE2_list, by='sectxspp')
EDGE2_list$EDGE2 <- NA
EDGE2_list$perc <- NA
EDGE2_list$list <- NA

m1 <- as.vector(unlist(as.vector(EDGE2_values))) %>% median(na.rm=T)

for (s in 1:nrow(EDGE2_list)) {
  d1 <- EDGE2_values[EDGE2_list$sectxspp[s],] %>% t() %>% as.vector()
  # m1 <- EDGE2_median[EDGE2_list$section[s],'V3']
  
  EDGE2_list$EDGE2[s] <- median(d1, na.rm=T)
  
  t1 <- table(d1>m1)

  if (is.na(t1['TRUE'])) {
    EDGE2_list$perc[s] <- 0
  } else {
    EDGE2_list$perc[s] <- t1['TRUE']/nr
  }
  
  if (EDGE2_list$perc[s] > 0.80 & EDGE2_list$category[s] %in% c("EN","CR","VU")) {  EDGE2_list$list[s] <- 'borderline' }
  if (EDGE2_list$perc[s] > 0.95 & EDGE2_list$category[s] %in% c("EN","CR","VU")) {  EDGE2_list$list[s] <- 'main' }
  if (EDGE2_list$perc[s] > 0.80 & EDGE2_list$category[s] %in% c("NE","DD")) {  EDGE2_list$list[s] <- 'research' }
  if (EDGE2_list$perc[s] > 0.95 & EDGE2_list$category[s] %in% c("NT","LC")) {  EDGE2_list$list[s] <- 'watch' }
  
}

write.csv(EDGE2_list, 'results/EDGE2_list.csv')


