

# DISTRIBUTION OF EDGE2, ED2, ECODGE AND FUD SCORES PER INFRAGENERIC TAXA


library(ggpubr)



df_EDGE2_values <- as.data.frame(EDGE2_values)
df_EDGE2_values$scientific_name <- rownames(df_EDGE2_values)
df_EDGE2_values$name <- NULL
df_EDGE2_values <- df_EDGE2_values %>% pivot_longer(1:500)
df_EDGE2_values <- merge(df_EDGE2_values, assessments[,c('scientific_name','section')])
df_EDGE2_values$section <- str_to_title(df_EDGE2_values$section)

m1 <- boxplot.stats(df_EDGE2_values$value)$stats[3]

p1 <- ggplot(aes(x=section, y=value), data=df_EDGE2_values) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(limits = c(0, 3)) +
  theme_bw() +
  xlab('') + ylab('EDGE2') +
  geom_hline(yintercept=m1, col='red') +
  theme(axis.text.x = element_text(angle = -90, vjust = 1, hjust=0))


df_EcoDGE2_values <- as.data.frame(EcoDGE2_values)
df_EcoDGE2_values$scientific_name <- rownames(df_EcoDGE2_values)
df_EcoDGE2_values$name <- NULL
df_EcoDGE2_values <- df_EcoDGE2_values %>% pivot_longer(1:500)
df_EcoDGE2_values <- merge(df_EcoDGE2_values, assessments[,c('scientific_name','section')])
df_EcoDGE2_values$section <- str_to_title(df_EcoDGE2_values$section)

m1 <- boxplot.stats(df_EcoDGE2_values$value)$stats[3]

p2 <- ggplot(aes(x=section, y=value), data=df_EcoDGE2_values) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(limits = c(0, 0.2)) +
  theme_bw() +
  xlab('') + ylab('EcoDGE2') +
  geom_hline(yintercept=m1, col='red') +
  theme(axis.text.x = element_text(angle = -90, vjust = 1, hjust=0))



df_EDGE2_values <- as.data.frame(ED2_values)
df_EDGE2_values$scientific_name <- rownames(df_EDGE2_values)
df_EDGE2_values$name <- NULL
df_EDGE2_values <- df_EDGE2_values %>% pivot_longer(1:500)
df_EDGE2_values <- merge(df_EDGE2_values, assessments[,c('scientific_name','section')])
df_EDGE2_values$section <- str_to_title(df_EDGE2_values$section)

m1 <- boxplot.stats(df_EDGE2_values$value)$stats[3]

p3 <- ggplot(aes(x=section, y=value), data=df_EDGE2_values) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(limits = c(0, 25)) +
  theme_bw() +
  xlab('') + ylab('ED2') +
  geom_hline(yintercept=m1, col='red') +
  theme(axis.text.x = element_text(angle = -90, vjust = 1, hjust=0))


df_EcoDGE2_values <- as.data.frame(FUD_values)
df_EcoDGE2_values$scientific_name <- rownames(df_EcoDGE2_values)
df_EcoDGE2_values$name <- NULL
df_EcoDGE2_values <- df_EcoDGE2_values %>% pivot_longer(1:500)
df_EcoDGE2_values <- merge(df_EcoDGE2_values, assessments[,c('scientific_name','section')])
df_EcoDGE2_values$section <- str_to_title(df_EcoDGE2_values$section)

m1 <- boxplot.stats(df_EcoDGE2_values$value)$stats[3]

p4 <- ggplot(aes(x=section, y=value), data=df_EcoDGE2_values) +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(limits = c(0, 4)) +
  theme_bw() +
  xlab('') + ylab('FUD') +
  geom_hline(yintercept=m1, col='red') +
  theme(axis.text.x = element_text(angle = -90, vjust = 1, hjust=0))



ggarrange(p1, p2, p3, p4,
          labels=c('a)', 'b)', 'c)', 'd)'),
          label.x=-0.01, ncol=2, nrow=2)


