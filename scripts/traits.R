
library(tidyverse)
library(GGally)
library(readxl)
library(missMDA)
library(ape)
library(cluster)
source('scripts/phylogeny.R') # need traits fro all species with molecular sequences (3 yet NE)


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

# filter species in scleria_iucn
table(scleria_iucn$species %in% traits$species)
traits <- traits %>% filter(species %in% scleria_iucn$species)

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



# create dendrogram with all considered traits: see Griffith et al 2022 (Functional Ecology) ####

# copy data
s_traits <- traits
rownames(s_traits) <- s_traits$species

# traits to use
v_traits <- c("inflorescence","life_form","leaf_length","nutlet_length","nutlet_width","log.height","log.leaf_width")
s_traits <- s_traits[,v_traits]
str(s_traits)

# scale continuous variables
s_traits[,3:7] <- scale(s_traits[,3:7])

# distance matrix + clustering UPGMA
daisy.mat <- daisy(s_traits[,v_traits], metric="gower") %>% as.dist()
dend1 <- hclust(daisy.mat, method="average") # UPGMA
dend1.phy <- as.phylo(dend1) # as.phylo

mycat <- scleria_iucn$category[order(match(scleria_iucn$species,dend1.phy$tip.label))]
mycat <- factor(mycat)
mycol <- c("azure4","black","forestgreen","yellow","orange2","red")[mycat]

par(mar=c(1,1,1,0))
plot(as.phylo(dend1.phy), tip.color=mycol, cex=0.7)

nr <- 100
imputed_dendrograms <- list()
for (r in 1:nr) { # randomly drop 1-3 traits per run
  ntr <- sample(4:6,1)
  dend.n <- s_traits %>% dplyr::select(sample(v_traits, ntr)) %>% daisy(metric="gower") %>% as.dist()
  dend.n <- hclust(dend.n, method="average") # UPGMA
  imputed_dendrograms[[r]] <- as.phylo(dend.n)
  print(r)
}

save(imputed_dendrograms, file="results/imputed_dendrograms.RData")

rm(s_traits, dend.n, ntr)
