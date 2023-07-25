

# IMPORT FUNCTIONAL TRAITS
# EXPLORE RELATIONSHIP BETWEEN LITERAUTRE AND HERBARIUM TRAITS
# AND COMPUTE LEAF AREA AND NUTLET VOLUME

library(tidyverse); library(readxl); library(GGally)
library(cluster); library(ape)


herbarium_traits <- read_excel("data/Kew_MNHNP_vouchers.xlsx", sheet='trait_data') %>% # import herbarium data
  dplyr::select(-catalogNumber, -stem_width, -vernacular, -uses, -region)
str(herbarium_traits); length(unique(herbarium_traits$scientific_name))


literature_traits <- read_excel("data/Kew_MNHNP_vouchers.xlsx", sheet='traits_literature') %>% # literature traits
  dplyr::select(-species, -nutlet_weight, -spikelet, -hypogynium, -inflorescence, -leaf, -notes)


herbarium_means <- herbarium_traits %>% group_by(scientific_name) %>% # compare literature and herbarium traits
  summarise(h.height=mean(height, na.rm=T),
            h.blade_length=mean(blade_length, na.rm=T),
            h.blade_width=mean(blade_width, na.rm=T),
            h.nutlet_length=mean(nutlet_length, na.rm=T),
            h.nutlet_width=mean(nutlet_width, na.rm=T))


comp_traits <- merge(herbarium_means, literature_traits, by='scientific_name')

# par(mfrow=c(2,3))
# plot(log(comp_traits$h.height), log(comp_traits$height),
#      main='log(height (cm))', xlab='herbarium', ylab='literature'); abline(0,1)
# plot(comp_traits$h.blade_length, comp_traits$blade_length,
#      main='blade_length (cm)', xlab='herbarium', ylab='literature'); abline(0,1)
# plot(comp_traits$h.blade_width, comp_traits$blade_width,
#      main='blade_width (cm)', xlab='herbarium', ylab='literature'); abline(0,1)
# plot(comp_traits$h.nutlet_length, comp_traits$nutlet_length,
#      main='nutlet_length (mm)', xlab='herbarium', ylab='literature'); abline(0,1)
# plot(comp_traits$h.nutlet_width, comp_traits$nutlet_width,
#      main='nutlet_width (mm)', xlab='herbarium', ylab='literature'); abline(0,1)



# lets merge species from the literature into the herbarium dataset
herbarium_traits$scientific_name[which(!(herbarium_traits$scientific_name %in% literature_traits$scientific_name))]
literature_traits$scientific_name[which(!(literature_traits$scientific_name %in% herbarium_traits$scientific_name))]

herbarium_traits <- literature_traits %>%
  dplyr::select(colnames(herbarium_traits)) %>%
  bind_rows(herbarium_traits)


# nutlet_volume	= 4/3π length*width^2	(ellipsoid)			
herbarium_traits$nutlet_volume <- 4/3 * pi * herbarium_traits$nutlet_length * herbarium_traits$nutlet_width^2	

# blade_area	= π length*width (ellipse)			
herbarium_traits$blade_area <- pi * herbarium_traits$blade_length * herbarium_traits$blade_width


# correlation nutlet volume and weight
nutlet_weight <- read_excel("data/Kew_MNHNP_vouchers.xlsx", sheet='nutlet_weight')

nutlet_means <- herbarium_traits %>% group_by(scientific_name) %>%
  summarise(nutlet_volume=mean(nutlet_volume, na.rm=T)) %>%
  dplyr::select(scientific_name, nutlet_volume)

nutlet_means <- merge(nutlet_weight[c('scientific_name','nutlet_weight')], nutlet_means)

# plot(nutlet_means$nutlet_volume, nutlet_means$nutlet_weight,
#      main='nutlet', xlab='volume (mm3)', ylab='weight (mg)')
# abline(lm(nutlet_means$nutlet_weight ~ nutlet_means$nutlet_volume))



# LHS Scleria
LHS_means <- herbarium_traits %>% dplyr::select(scientific_name, height, blade_area, nutlet_volume) %>%
  group_by(scientific_name) %>%
  summarise(height=mean(height, na.rm=T),
            blade_area=mean(blade_area, na.rm=T),
            nutlet_volume=mean(nutlet_volume, na.rm=T))
  
LHS_means <- merge(literature_traits[,c('subgenus','section','scientific_name','life_form')],
                  LHS_means, by='scientific_name')

LHS_raw <- merge(literature_traits[,c('subgenus','section','scientific_name','life_form')],
                 herbarium_traits[,c('scientific_name', 'height', 'blade_area', 'nutlet_volume')], by='scientific_name')

write.table(LHS_means, 'results/LHS_means_final.txt') # save results
write.table(LHS_raw, 'results/LHS_raw_final.txt') # save results
rm(comp_traits, herbarium_means, herbarium_traits, literature_traits, nutlet_weight, nutlet_means)


long_LHS_data <- LHS_data %>% pivot_longer(5:7, names_to='trait', values_to='value')

ggplot(aes(x=subgenus, y=value), data=long_LHS_data) +
  geom_boxplot() +
  facet_wrap(.~trait, scales="free") +
  theme_bw()
  
# library(rgl)
# mycolors <- c('brown4', 'green', 'yellow', 'blue')
# LHS_data$color <- mycolors[ as.numeric(as.factor(LHS_data$subgenus)) ]
# 
# plot3d( 
#   x=log(LHS_data$height), y=log(LHS_data$blade_area), z=log(LHS_data$nutlet_volume), 
#   xlab="height", ylab="blade_area", zlab="nutlet_volume",
#   col=LHS_data$color)
  
