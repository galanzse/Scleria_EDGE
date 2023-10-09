

# Categorical Analysis of Neo- And Paleo-Endemism

library(canaper)
library(tidyverse)
library(terra)


# data import
load('results/data_final.RData')
load('results/imputed_trees.RData')



# import occurrences and remove species without infrageneric information
scleria_occ <- data_final[['occurrences']] %>%
  dplyr::select(scientific_name, x, y) %>%
  merge(data_final[['assessments']][,c('scientific_name','section','subgenus')]) %>%
  subset(scientific_name %in% imputed_trees[[1]]$tip.label)

scleria_occ$scientific_name[!(scleria_occ$scientific_name %in% imputed_trees[[1]]$tip.label)] %>% unique()

# points file
p_scleria_occ <- vect(scleria_occ, geom=c('x','y'), 'epsg:4326') %>% terra::project('+proj=eqearth') # points file



# import ecoregions and create presence matrix
scl_ecor_map <- terra::vect('C:/Users/user/Desktop/ecoregions/wwf_terr_ecos.shp') %>%
  terra::crop(ext(-180,180,-50,70)) %>% terra::project('+proj=eqearth')

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
n_null1 <- "r2" # keep species richness and probabilities

canape_results <- matrix(nrow=nrow(scl_ecor_mat), ncol=length(imputed_trees))
rownames(canape_results) <- rownames(scl_ecor_mat)
for (i in 1:length(imputed_trees)) {
  canape1 <- cpr_rand_test(comm=scl_ecor_mat, null_model="r0", n_reps=100, phy=imputed_trees[[i]], metrics=c("pe","rpe"))
  canape1 <- cpr_classify_endem(canape1)
  canape_results[,i] <- canape1$endem_type
  print(paste('--- ', round(i/length(imputed_trees),2)*100, '% ---', sep=''))
}

write.table(canape_results, row.names=T, 'results/canape_results.txt')


