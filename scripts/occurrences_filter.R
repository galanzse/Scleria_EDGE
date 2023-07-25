

# CURATE OCCURRENCES SPATIALLY & CLIMATICALLY
# THE DATABASE HAS BEEN ALREADY MANUALLY CURATED AND IT IS TAXONOMICALLY CORRECT

library(tidyverse); library(readxl)
library(terra); library(maptools); library(dismo)


occ_scleria <- read_excel("data/occ_scleria_taxcorrected.xlsx") %>% dplyr::select(-source, -id) # import occurrences
colnames(occ_scleria) <- c('species','x','y') # correct variable names
pts_scleria <- vect(occ_scleria, geom=c("x","y")) # create points objects


scl_countries <- read_excel("data/Kew_MNHNP_vouchers.xlsx", sheet = "botanical_countries_wcspf") # botanical countries x species
bot_countries <- read_excel("data/occurrences/tblLevel3.xlsx") # botanical countries x political country


data(wrld_simpl) # import world map
wrld_simpl <- vect(wrld_simpl)
class(wrld_simpl)


plot(wrld_simpl); points(pts_scleria, col='red')


# correct common misidentifications of widespread species
# plot(wrld_simpl, main='Scleria secans'); points(pts_scleria[pts_scleria$species=='Scleria secans (L.) Urb.',])
occ_scleria$species[occ_scleria$species=='Scleria secans (L.) Urb.' & occ_scleria$x > 0] <- 'Scleria boivinii Steud.'
# plot(wrld_simpl, main='Scleria secans'); points(pts_scleria[pts_scleria$species=='Scleria hirtella Sw.',])
occ_scleria$species[occ_scleria$species=='Scleria hirtella Sw.' & occ_scleria$x > -30] <- 'Scleria distans Poir.'
pts_scleria <- vect(occ_scleria, geom=c("x","y")) # redo points dataframe


regions <- terra::extract(x=wrld_simpl, y=pts_scleria) %>% dplyr::select(FIPS, ISO2, ISO3) # extract regions
occ_scleria <- cbind(occ_scleria, regions)

occ_scleria <- occ_scleria %>% subset(!is.na(ISO2)) # remove points outside polygons
pts_scleria <- vect(occ_scleria, geom=c("x","y")) # redo points objects
crs(pts_scleria) <- '+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs'

unique(occ_scleria$species)[!(unique(occ_scleria$species) %in% scl_countries$scientific_name)]


MAT <- rast('C:/Users/user/Desktop/wc2.1_30s_bio/wc2.1_30s_bio_1.tif') # mean annual temperature
AP <- rast('C:/Users/user/Desktop/wc2.1_30s_bio/wc2.1_30s_bio_12.tif') # annual precipitation

GRD <- MAT %>% terra::crop(pts_scleria) %>% terra::aggregate(3) %>%
  terra::project('+proj=eqearth') # use equal area projection for spatial filter
# plot(GRD)
# cellSize(GRD, unit='km') # 2.15km2


occ_scleria$MAT <- terra::extract(MAT, pts_scleria, method='simple')[,2]
occ_scleria$AP <- terra::extract(AP, pts_scleria, method='simple')[,2]
occ_scleria$cell <- terra::extract(GRD, terra::project(pts_scleria, '+proj=eqearth'), method='simple')[,2]
# plot(GRD, main='MAT'); points(terra::project(pts_scleria, '+proj=eqearth'))


occ_scleria2 <- list() # list to store results

for (s in unique(occ_scleria$species)) {
  
  # filter observations of species s
  spp_temp <- occ_scleria %>% subset(species==s)
  # plot(wrld_simpl, main=s); points(vect(spp_temp, geom=c('x','y')), col='red')
  
  # extract botanical country information
  count_temp <- scl_countries %>% subset(scientific_name==s) %>% dplyr::select(countries) %>%
    deframe() %>% strsplit(" ") %>% deframe()
  
  # cross botanical country x iso code data
  count_temp <- bot_countries %>% subset(L3_code %in% count_temp) %>%
    dplyr::select(L3_ISOcode) %>% deframe() %>% unique()
  
  # retain observations based on iso2 code and save
  count_temp <- spp_temp[spp_temp$ISO2 %in% count_temp,]
  
  # sampling bias: select one observation per 2.15km2 cell
  ss_df <- count_temp %>% group_by(cell) %>% slice_sample(n=1)
  
  # climatic filter: only for species with more than X observations
  if (nrow(ss_df) > 5) {
    ss_df <- ss_df %>% filter(MAT > boxplot.stats(ss_df$MAT, coef=2)$stats[1] & MAT < boxplot.stats(ss_df$MAT, coef=2)$stats[5] &
                                AP > boxplot.stats(ss_df$AP, coef=2)$stats[1] & AP < boxplot.stats(ss_df$AP, coef=2)$stats[5] )
  }
  
  # save
  occ_scleria2[[s]] <- ss_df
  
  print(s)
}

occ_scleria_filt <- do.call(rbind, occ_scleria2) # paste list
pts_scleria_filt <- vect(occ_scleria_filt, geom=c("x","y")) # redo points dataframe

length(unique(occ_scleria_filt$species))
unique(occ_scleria$species)[!(unique(occ_scleria$species) %in% unique(occ_scleria_filt$species))] # check excluded species


write.table(occ_scleria_filt, 'results/occ_scleria_filt.txt') # save results
rm(AP, MAT, GRD, bot_countries, scl_countriescount_temp, occ_scleria, occ_scleria2, regions, spp_temp, ss_df, s, pts_scleria)

