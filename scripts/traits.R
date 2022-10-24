
library(tidyverse)
library(GGally)
library(readxl)
library(missMDA)


# import trait data
traits_isabel <- read_excel("data/traits_Isabel.xlsx")
traits_javi <- read_excel("data/traits.xlsx")

# check my species are in Isabel's dataset
table(traits_isabel$WCVP_name %in% traits_javi$species)
traits_isabel$WCVP_name[!(traits_isabel$WCVP_name %in% traits_javi$species)]

# merge both datasets
colnames(traits_isabel)[colnames(traits_isabel)=='WCVP_name'] <- 'species'
traits <- merge(traits_javi, traits_isabel, by='species', all.x=T)
rm(traits_isabel, traits_javi)

# check NAs per column
colSums(is.na(traits))[order(colSums(is.na(traits)))]

# remove traits with many NAs or very conserved across the phylonegy (see Bauters et al. 2016)
traits <- traits[,!(colnames(traits) %in% c('pollination - animal','pollination','bat rough generic score yes/no',
                                            'source','photosynthesis_C13','Family','climate_description',
                                            'photosynthesis_C13','phytosynthesis_type','geographic_area',
                                            'halophyte','plant_name_id','open_closed_habitat','epiphyte_lithophyte',
                                            'wet_habitat','annual','spikelet','hypogynium','leaf',
                                            'lifeform_description'))]

str(traits)
traits$inflorescence <- as.factor(traits$inflorescence)
traits$life_form <- as.factor(traits$life_form)

# log trans
traits$log.height <- log(traits$height) # height and leaf width have right-skewed distributions
traits$log.leaf_width <- log(traits$leaf_width)

# correlation
# ggpairs(traits[,c('log.height','leaf_length','log.leaf_width','nutlet_length','nutlet_width','nutlet_mass')])


# impute nutlet traits using a local regression and length as predictor
lm_nut <- loess(nutlet_width ~ nutlet_length, data=traits, span=1)
ggplot(aes(y=nutlet_width, x=nutlet_length), data=traits) +
  geom_point() + geom_smooth(method='loess', span=1)
width_pred <- predict(lm_nut, newdata=traits)
for (i in 1:nrow(traits)) {
  if (is.na(traits$nutlet_width[i])) { traits$nutlet_width[i] <- width_pred[i] }
}

# impute other cuantitative traits using PCA
traits[,c('log.height','leaf_length','log.leaf_width','nutlet_length','nutlet_width')] <- imputePCA(traits[,c('log.height','leaf_length','log.leaf_width','nutlet_length','nutlet_width')], ncp=2, scale=TRUE)$completeObs

# impute qualitative traits using mode
sp_im <- traits %>% filter(is.na(inflorescence) | is.na(life_form)) %>% dplyr::select(species) %>% deframe()
for (sp in sp_im) {
  df_sp <-traits %>% filter(species==sp)
  if (is.na(traits$inflorescence[traits$species==sp])) {
    traits$inflorescence[traits$species==sp] <- traits$inflorescence[traits$section==df_sp$section] %>% table() %>% which.max() %>% names()
  }
  if (is.na(traits$life_form[traits$species==sp])) {
    traits$life_form[traits$species==sp] <- traits$life_form[traits$section==df_sp$section] %>% table() %>% which.max() %>% names()
  }
}


# correlation nutlet size ~ mass
ggplot(aes(y=nutlet_mass, x=nutlet_length), data=traits) + 
  geom_point() +theme_classic() +
  ylim(0,80) + xlim(0,5) +
  geom_smooth(method='loess', span=4, se=F) +
  labs(x='Nutlet length (mm)', y='1000 seed weight (g)')

