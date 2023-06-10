
library(sp)
library(rgdal)
library(dplyr)
library(here)
library(tidyverse)
library(lubridate)

# read in data

for(missingPollutant in c("no2o3missing"))
  # for(missingPollutant in c("no2missing","o3missing","stationpm25missing","stationpm10missing","stationno2missing","stationo3missing"))
{
  
  ######################################### load data #####################################
  # Please define your output directory, output_dir
  output_dir <- paste0("results/",missingPollutant)
  
  load(file.path(output_dir,"test_data.RData"))
  
  data_all_guangzhou <- test_data
  
  ######################################### get model input format #####################################
  
  start_year <- 2017
  end_year <- 2022
  
  data_all_guangzhou$date = strptime(
    as.character(data_all_guangzhou$date), format = '%Y-%m-%d %H:%M:%S',
    tz='UTC'
  )
  
  # the firstDay is used for day difference calculate
  firstDay = as.Date(paste0(start_year-1,"-12-31"))
  endDay = as.Date(paste0(end_year,"-12-31"))
  
  
  data_all_guangzhou$date <- trunc(data_all_guangzhou$date,units='day')
  
  # try modelling daily 
  data_2 <- data_all_guangzhou %>% 
    group_by(pollutant, date) %>% 
    summarise(pollutionValue=mean(value,na.rm=TRUE))
  
  data_2$date <- as.Date(data_2$date)
  
  # fill in full days
  
  full_date <- rbind(data.frame(pollutant="pm25", date = seq.Date(as.Date("2017-01-01"), as.Date("2022-12-31"), by="day")),
                     data.frame(pollutant="pm10", date = seq.Date(as.Date("2017-01-01"), as.Date("2022-12-31"), by="day")),
                     data.frame(pollutant="no2", date = seq.Date(as.Date("2017-01-01"), as.Date("2022-12-31"), by="day")),
                     data.frame(pollutant="o3", date = seq.Date(as.Date("2017-01-01"), as.Date("2022-12-31"), by="day")))
  
  
  data_all <- full_date %>% dplyr::left_join(data_2, by=(c("pollutant","date")))
  
  # actually, stan can not handle NA, so this full fill could be useful while doing prediction
  
  data_all <- na.omit(data_all)
  
  
  
  ###################################
  # get the input data
  ###################################
  
  library(rstan)
  library(robustbase)
  library(data.table)
  library(here)
  
  
  start_year=2017
  end_year=2022
  # the firstDay is used for day difference calculate
  firstDay = as.Date(paste0(start_year-1,"-12-31"))
  endDay = as.Date(paste0(end_year,"-12-31"))
  
  # get the covariates ready for stan model ---------------------------------
  
  data_all$day = trunc(data_all$date, units='days')
  
  data_all$cos12 = cos(2*pi*as.numeric(data_all$day)/(365.25))
  data_all$cos6 = cos(2*2*pi*as.numeric(data_all$day)/(365.25))
  data_all$sin12 = sin(2*pi*as.numeric(data_all$day)/(365.25))
  data_all$sin6 = sin(2*2*pi*as.numeric(data_all$day)/(365.25))
  
  
  # data_all$dowChinese = factor(weekdays(data_all$day), levels = 
  # weekdays(seq(ISOdate(2000,1,3), len=7, by='days')))
  data_all$dow <- lubridate::wday(data_all$date,week_start = getOption("lubridate.week.start", 1))
  
  data_all$dayInt = as.numeric(difftime(data_all$day, firstDay, units='days'))
  
  
  data_all$pollutant <- as.numeric(factor(data_all$pollutant, levels = c("pm25","pm10","no2","o3")))
  
  
  # covariates for varying covariance
  data_temp <- seq(as.Date(paste0(start_year,"-01-01")), as.Date(paste0(end_year,"-12-31")), by="day")
  
  x_cov = data.frame(cos12=cos(2*pi*as.numeric(data_temp)/(365.25))
                     ,cos6=cos(2*2*pi*as.numeric(data_temp)/(365.25))
                     ,sin12=sin(2*pi*as.numeric(data_temp)/(365.25))
                     ,sin6=sin(2*2*pi*as.numeric(data_temp)/(365.25))
  )
  
  
  
  current_data <- list(N=nrow(data_all), 
                       K=4, 
                       x=data_all[,c("cos12", "cos6", "sin12", "sin6")], 
                       y=data_all$pollutionValue,
                       Nd=nrow(x_cov), 
                       day=data_all$dayInt,
                       Nw=7, 
                       week=data_all$dow,
                       N_p=length(unique(data_all$pollutant)),
                       x_cov=x_cov,
                       N_cov=ncol(x_cov),
                       pollutant=data_all$pollutant,
                       trend_reference_day=round(nrow(x_cov)/2)
                       
  )
  
  ### fit model
  
  
  
  N_p=current_data$N_p
  
  options(mc.cores=parallel::detectCores())
  
  initf <- function() list(
    alpha=rep(2,N_p)+rnorm(N_p, mean=0, sd=0.01)
    ,miu=rnorm(N_p, mean=0, sd=0.01)
    ,arcoef_phi=rep(0.3,N_p)+rnorm(N_p, mean=0, sd=0.01)
    # ,betas=coef_store[,2:7]+rnorm(N_p*ncol(coef_store[,2:7]), mean=0, sd=0.01) # as we do not use cos sin 3
    ,betas=matrix(0, nrow=current_data$N_p, ncol = current_data$K)+rnorm(current_data$N_p*current_data$K, mean=0, sd=0.01)
    ,factor_w_pre=array(0,c(N_p,6))+rnorm(N_p*6, mean=0, sd=0.01)
    ,factor_d=matrix(0, nrow=current_data$Nd, ncol = N_p)+rnorm(current_data$Nd*N_p, mean=0, sd=0.01)
  )
  
  
  fit.code <- stanc("Daily_pollution_model.stan") # convert to C++ code
  fit.model <- stan_model(stanc_ret = fit.code) # compile C++ code
  
  rstan_options(auto_write = TRUE)
  
  
  # Sample from Stan model
  fit1 <- rstan::sampling(fit.model, data=current_data,chains = 1, iter = 1000,
                          control = list(max_treedepth = 10, adapt_delta=0.8 ,stepsize_jitter=0.6)
                          ,seed = 100
                          , pars =c("betas", "miu","alpha","factor_w","factor_d", "arcoef_phi","tau","miu_cov","betas_cov","L_save","factor_trend")
                          , verbose=T
                          ,init = initf)
  
  
  save(fit1, file = file.path(output_dir, "Daily_newchain1.RData"))
  
  fit2 <- rstan::sampling(fit.model, data=current_data,chains = 1, iter = 1000,
                          control = list(max_treedepth = 10, adapt_delta=0.8 ,stepsize_jitter=0.6)
                          ,seed = 200
                          , pars =c("betas", "miu","alpha","factor_w","factor_d", "arcoef_phi","tau","miu_cov","betas_cov","L_save","factor_trend")
                          , verbose=T
                          ,init = initf)
  
  
  save(fit2, file = file.path(output_dir, "Daily_newchain2.RData"))
}



