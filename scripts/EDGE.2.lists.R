

# classify species of a given clade in EDGE2 lists following standard criteria
# see 10.1371/journal.pbio.3001991 for more info


# provide two objects:
  # 1/ species x EDGE2 scores matrix, where each column includes the scores for a given run 
  # 2/ a dataframe with two columns: the first comprising species names, the second comprising their IUCN status (threatened/nonthreatened)
# function returns a dataframe with 5 columns: species, EDGE2 median and sd, percentile and list

# EDGE2_scores <- EDGE2_values
# IUCN_status <- assessments[,c('scientific_name','thr')]

EDGE2_lists <- function(EDGE2_scores, IUCN_status){

  # check format
  if (!(is.matrix(EDGE2_scores))) { warning("EDGE2_scores not a matrix")  }
  if (!(is.data.frame(IUCN_status))) { warning("IUCN_status not a data.frame")  }
  
  # new variables
  df_lists <- IUCN_status
  df_lists$EDGE2 <- NA
  df_lists$EDGE2_sd <- NA
  df_lists$perc <- NA
  df_lists$list <- NA
  
  # clade median
  m1 <- median(EDGE2_scores, na.rm=T)
  
  # classify species
  for (s in 1:nrow(df_lists)) {
    
    df_lists$EDGE2[s] <- median(EDGE2_scores[df_lists$scientific_name[s],], na.rm=T) # median
    df_lists$EDGE2_sd[s] <- sd(EDGE2_scores[df_lists$scientific_name[s],], na.rm=T) # sd
    
    t1 <- table(EDGE2_scores[df_lists$scientific_name[s],] > m1) # percentile
    
    if (is.na(t1['TRUE'])) { # check is species is always under the median
      df_lists$perc[s] <- 0
      } else {
        df_lists$perc[s] <- t1['TRUE']/nr
        }
  
  # classify according to lists
  if (df_lists$perc[s] > 0.80 & df_lists$thr[s]=='threatened') {  df_lists$list[s] <- 'borderline' }
  if (df_lists$perc[s] > 0.95 & df_lists$thr[s]=='threatened') {  df_lists$list[s] <- 'main' }
  if (df_lists$perc[s] > 0.95 & df_lists$thr[s]=='nonthreatened') {  df_lists$list[s] <- 'watch' }


  print(paste('--- ', round(s/nrow(df_lists)*100, 2), '% ---', sep='')) # status

  }
  
  return(df_lists)
}

