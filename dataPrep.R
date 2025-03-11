library(tidyverse)

load('data_2020.Rdata')

data.simp <- data.final %>%
  select(-c(source, timestamp, year, occ, occ.reg, town, pc.district.number, centroid_long, centroid_lat))

write.csv(data.simp, 'data_simp.csv', row.names=F)
