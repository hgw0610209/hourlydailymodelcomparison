library(rstan)
library(dplyr)
library(tidyverse)
library(lubridate)
library(rstanarm) # fit linear model in Bayesian

library(robustbase)
library(data.table)
library(here)


load("data_all_guangzhou.RData")

data_all_guangzhou$date = strptime(
  as.character(data_all_guangzhou$date), format = '%Y-%m-%d %H:%M:%S',
  tz='UTC'
)
data_all_guangzhou$day = trunc(data_all_guangzhou$date, units='days')

numberofStation <- length(unique(data_all_guangzhou$station))

### single pollutant missing, take 50 days out #######
set.seed(128)
valid_data=data_all_guangzhou %>% 
  group_by(pollutant,day) %>% 
  summarise(non_na_count = sum(!is.na(value)), meanValue=mean(value,na.rm = T)) %>% 
  filter(pollutant=="pm25", non_na_count==24*numberofStation) %>% 
  sample_n(50, replace=FALSE)
test_data=data_all_guangzhou %>% 
  anti_join(valid_data, by=c("pollutant","day"))
save(valid_data, file = "results/pm25missing/valid_data.RData")
save(test_data, file = "results/pm25missing/test_data.RData")

valid_data=data_all_guangzhou %>% 
  group_by(pollutant,day) %>% 
  summarise(non_na_count = sum(!is.na(value)), meanValue=mean(value,na.rm = T)) %>% 
  filter(pollutant=="pm10", non_na_count==24*numberofStation) %>% 
  sample_n(50, replace=FALSE)
test_data=data_all_guangzhou %>% 
  anti_join(valid_data, by=c("pollutant","day"))
save(valid_data, file = "results/pm10missing/valid_data.RData")
save(test_data, file = "results/pm10missing/test_data.RData")

valid_data=data_all_guangzhou %>% 
  group_by(pollutant,day) %>% 
  summarise(non_na_count = sum(!is.na(value)), meanValue=mean(value,na.rm = T)) %>% 
  filter(pollutant=="no2", non_na_count==24*numberofStation) %>% 
  sample_n(50, replace=FALSE)
test_data=data_all_guangzhou %>% 
  anti_join(valid_data, by=c("pollutant","day"))
save(valid_data, file = "results/no2missing/valid_data.RData")
save(test_data, file = "results/no2missing/test_data.RData")

valid_data=data_all_guangzhou %>% 
  group_by(pollutant,day) %>% 
  summarise(non_na_count = sum(!is.na(value)), meanValue=mean(value,na.rm = T)) %>% 
  filter(pollutant=="o3", non_na_count==24*numberofStation) %>% 
  sample_n(50, replace=FALSE)
test_data=data_all_guangzhou %>% 
  anti_join(valid_data, by=c("pollutant","day"))
save(valid_data, file = "results/o3missing/valid_data.RData")
save(test_data, file = "results/o3missing/test_data.RData")

#### two-pollutant missing ######

set.seed(128)
valid_data=data_all_guangzhou %>% 
  group_by(pollutant,day) %>% 
  summarise(non_na_count = sum(!is.na(value)), meanValue=mean(value,na.rm = T)) %>% 
  filter(pollutant %in% c("pm25","pm10"), non_na_count==24*numberofStation) %>% 
  pivot_wider(names_from = pollutant, values_from = meanValue)%>% 
  na.omit() %>% 
  sample_n(50, replace=FALSE) %>% 
  pivot_longer(cols = contains(c("pm25","pm10")),names_to = "pollutant", values_to = "value")

test_data=data_all_guangzhou %>% 
  anti_join(valid_data, by=c("pollutant","day"))

save(valid_data, file = "results/pm25pm10missing/valid_data.RData")
save(test_data, file = "results/pm25pm10missing/test_data.RData")



valid_data=data_all_guangzhou %>% 
  group_by(pollutant,day) %>% 
  summarise(non_na_count = sum(!is.na(value)), meanValue=mean(value,na.rm = T)) %>% 
  filter(pollutant %in% c("pm25","no2"), non_na_count==24*numberofStation) %>% 
  pivot_wider(names_from = pollutant, values_from = meanValue)%>% 
  na.omit() %>% 
  sample_n(50, replace=FALSE) %>% 
  pivot_longer(cols = contains(c("pm25","no2")),names_to = "pollutant", values_to = "value")

test_data=data_all_guangzhou %>% 
  anti_join(valid_data, by=c("pollutant","day"))

save(valid_data, file = "results/pm25no2missing/valid_data.RData")
save(test_data, file = "results/pm25no2missing/test_data.RData")




valid_data=data_all_guangzhou %>% 
  group_by(pollutant,day) %>% 
  summarise(non_na_count = sum(!is.na(value)), meanValue=mean(value,na.rm = T)) %>% 
  filter(pollutant %in% c("pm25","o3"), non_na_count==24*numberofStation) %>% 
  pivot_wider(names_from = pollutant, values_from = meanValue)%>% 
  na.omit() %>% 
  sample_n(50, replace=FALSE) %>% 
  pivot_longer(cols = contains(c("pm25","o3")),names_to = "pollutant", values_to = "value")

test_data=data_all_guangzhou %>% 
  anti_join(valid_data, by=c("pollutant","day"))

save(valid_data, file = "results/pm25o3missing/valid_data.RData")
save(test_data, file = "results/pm25o3missing/test_data.RData")




valid_data=data_all_guangzhou %>% 
  group_by(pollutant,day) %>% 
  summarise(non_na_count = sum(!is.na(value)), meanValue=mean(value,na.rm = T)) %>% 
  filter(pollutant %in% c("pm10","no2"), non_na_count==24*numberofStation) %>% 
  pivot_wider(names_from = pollutant, values_from = meanValue)%>% 
  na.omit() %>% 
  sample_n(50, replace=FALSE) %>% 
  pivot_longer(cols = contains(c("pm10","no2")),names_to = "pollutant", values_to = "value")

test_data=data_all_guangzhou %>% 
  anti_join(valid_data, by=c("pollutant","day"))

save(valid_data, file = "results/pm10no2missing/valid_data.RData")
save(test_data, file = "results/pm10no2missing/test_data.RData")



valid_data=data_all_guangzhou %>% 
  group_by(pollutant,day) %>% 
  summarise(non_na_count = sum(!is.na(value)), meanValue=mean(value,na.rm = T)) %>% 
  filter(pollutant %in% c("pm10","o3"), non_na_count==24*numberofStation) %>% 
  pivot_wider(names_from = pollutant, values_from = meanValue)%>% 
  na.omit() %>% 
  sample_n(50, replace=FALSE) %>% 
  pivot_longer(cols = contains(c("pm10","o3")),names_to = "pollutant", values_to = "value")

test_data=data_all_guangzhou %>% 
  anti_join(valid_data, by=c("pollutant","day"))

save(valid_data, file = "results/pm10o3missing/valid_data.RData")
save(test_data, file = "results/pm10o3missing/test_data.RData")




valid_data=data_all_guangzhou %>% 
  group_by(pollutant,day) %>% 
  summarise(non_na_count = sum(!is.na(value)), meanValue=mean(value,na.rm = T)) %>% 
  filter(pollutant %in% c("no2","o3"), non_na_count==24*numberofStation) %>% 
  pivot_wider(names_from = pollutant, values_from = meanValue)%>% 
  na.omit() %>% 
  sample_n(50, replace=FALSE) %>% 
  pivot_longer(cols = contains(c("no2","o3")),names_to = "pollutant", values_to = "value")

test_data=data_all_guangzhou %>% 
  anti_join(valid_data, by=c("pollutant","day"))

save(valid_data, file = "results/no2o3missing/valid_data.RData")
save(test_data, file = "results/no2o3missing/test_data.RData")


## we want to compare to the station level missing?? ##

set.seed(128)
valid_data=data_all_guangzhou %>% 
  group_by(pollutant,station, day) %>% 
  summarise(non_na_count = sum(!is.na(value)), meanValue=mean(value,na.rm = T)) %>% 
  filter(pollutant=="pm25", non_na_count==24) %>% 
  ungroup() %>% 
  sample_n(50, replace=FALSE)
test_data=data_all_guangzhou %>% 
  anti_join(valid_data, by=c("pollutant","station","day"))
save(valid_data, file = "results/stationpm25missing/valid_data.RData")
save(test_data, file = "results/stationpm25missing/test_data.RData")


valid_data=data_all_guangzhou %>% 
  group_by(pollutant,station, day) %>% 
  summarise(non_na_count = sum(!is.na(value)), meanValue=mean(value,na.rm = T)) %>% 
  filter(pollutant=="pm10", non_na_count==24) %>% 
  ungroup() %>% 
  sample_n(50, replace=FALSE)
test_data=data_all_guangzhou %>% 
  anti_join(valid_data, by=c("pollutant","station","day"))
save(valid_data, file = "results/stationpm10missing/valid_data.RData")
save(test_data, file = "results/stationpm10missing/test_data.RData")


valid_data=data_all_guangzhou %>% 
  group_by(pollutant,station, day) %>% 
  summarise(non_na_count = sum(!is.na(value)), meanValue=mean(value,na.rm = T)) %>% 
  filter(pollutant=="no2", non_na_count==24) %>% 
  ungroup() %>% 
  sample_n(50, replace=FALSE)
test_data=data_all_guangzhou %>% 
  anti_join(valid_data, by=c("pollutant","station","day"))
save(valid_data, file = "results/stationno2missing/valid_data.RData")
save(test_data, file = "results/stationno2missing/test_data.RData")


valid_data=data_all_guangzhou %>% 
  group_by(pollutant,station, day) %>% 
  summarise(non_na_count = sum(!is.na(value)), meanValue=mean(value,na.rm = T)) %>% 
  filter(pollutant=="o3", non_na_count==24) %>% 
  ungroup() %>% 
  sample_n(50, replace=FALSE)
test_data=data_all_guangzhou %>% 
  anti_join(valid_data, by=c("pollutant","station","day"))
save(valid_data, file = "results/stationo3missing/valid_data.RData")
save(test_data, file = "results/stationo3missing/test_data.RData")


