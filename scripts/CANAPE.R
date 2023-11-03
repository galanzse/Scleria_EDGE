

# Categorical Analysis of Neo- And Paleo-Endemism

library(canaper)
library(tidyverse)
library(terra)
library(picante)


# data import
load('results/data_final.RData')
load('results/imputed_trees.RData')



# import occurrences and remove species without infrageneric information
scleria_occ <- data_final[['occurrences']] %>%
  dplyr::select(scientific_name, x, y) %>%
  subset(scientific_name %in% imputed_trees[[1]]$tip.label)

scleria_occ$scientific_name[!(scleria_occ$scientific_name %in% imputed_trees[[1]]$tip.label)] %>% unique()

# points file
p_scleria_occ <- vect(scleria_occ, geom=c('x','y'), crs='epsg:4326') # points file



# import ecoregions and create presence matrix
dir_ecor <- c('C:/Users/user/Desktop/ecoregions/wwf_terr_ecos.shp')
scl_ecor_map <- terra::vect(dir_ecor) %>% terra::crop(ext(-180,180,-50,70))

scl_ecor_mat <- scl_ecor_map %>% terra::extract(p_scleria_occ) %>% dplyr::select(ECO_NAME)
scl_ecor_mat$scientific_name <- scleria_occ$scientific_name
scl_ecor_mat <- scl_ecor_mat %>% unique()
scl_ecor_mat$presence <- 1

scl_ecor_mat <- scl_ecor_mat %>% pivot_wider(names_from=scientific_name, values_from='presence') %>%
  as.data.frame()
scl_ecor_mat[is.na(scl_ecor_mat)] <- 0
rownames(scl_ecor_mat) <- scl_ecor_mat$ECO_NAME; scl_ecor_mat$ECO_NAME <- NULL
head(scl_ecor_mat)



# run analyses for every imputed tree
?vegan::commsim()
n_null1 <- "swap" # keep species richness and frequencies (same as in Mishler et al. 2014)


cpr_iter_sim(scl_ecor_mat, null_model = "swap")

canape_results <- matrix(nrow=nrow(scl_ecor_mat), ncol=length(imputed_trees))
rownames(canape_results) <- rownames(scl_ecor_mat)
for (i in 1:length(imputed_trees)) {
  canape1 <- cpr_rand_test(comm=scl_ecor_mat, null_model=n_null1, n_reps=100, n_iterations=100, phy=imputed_trees[[i]], metrics=c("pe","rpe"))
  canape1 <- cpr_classify_endem(canape1)
  canape_results[,i] <- canape1$endem_type
  print(paste('--- ', round(i/length(imputed_trees)*100,2), '% ---', sep=''))
}

write.table(canape_results, row.names=T, 'results/canape_results.txt')


unique(canape_results[])

perc_canape_results <- matrix(nrow=nrow(canape_results), ncol=5) %>% as.data.frame()
rownames(perc_canape_results) <- rownames(canape_results)
colnames(perc_canape_results) <- c('paleo', 'neo', 'not_significant', 'mixed', 'super')
for (i in 1:nrow(perc_canape_results)) {
  temp <- canape_results[i,] %>% t() %>% table()
  perc_canape_results$paleo[i] <- temp['paleo']/length(imputed_trees)*100
  perc_canape_results$neo[i] <- temp['neo']/length(imputed_trees)*100
  perc_canape_results$not_significant[i] <- temp['not significant']/length(imputed_trees)*100
  perc_canape_results$mixed[i] <- temp['mixed']/length(imputed_trees)*100
  perc_canape_results$super[i] <- temp['super']/length(imputed_trees)*100
}

perc_canape_results[is.na(perc_canape_results)] <- 0
write.table(perc_canape_results, row.names=T, 'results/perc_canape_results.txt', sep=',')


