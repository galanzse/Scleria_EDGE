

library(tidyverse)
library(ggplot2)


EDGE2_list <- read.csv("results/EDGE2_list.csv", row.names=1)

EcoDGE_list <- read.csv("results/EcoDGE2_list.csv", row.names=1)
colnames(EcoDGE_list)[colnames(EcoDGE_list)%in%c("EDGE2","ED2")] <- c("EcoDGE","FUD")


EDGE_cor <- merge(EDGE2_list[,c('scientific_name','EDGE2','ED2')],
                  EcoDGE_list[,c('scientific_name',"EcoDGE","FUD")])


ggplot(aes(x=log(EDGE2), y=log(EcoDGE)), data=EDGE_cor) +
  geom_point() +
  geom_smooth(method='lm')


ggplot(aes(x=ED2, y=FUD), data=EDGE_cor) +
  geom_point() +
  geom_smooth(method='lm') +
  theme_bw()

anova(lm(EDGE_cor$FUD ~ EDGE_cor$ED2))
