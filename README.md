# hourlydailymodelcomparison

Supplementary material for a paper comparing hourly and daily model strategies.

The main files are:
1. data_all_guangzhou.RData; The tidy pollution data of Guangzhou city for the study.
2. obtain_data_with_specific_missing_scenarios.R; The code to obtain sub-data with specific missing scenarios.
3. Daily_pollution_model.stan; The stan code for the new proposed daily pollution model. Note that the hourly pollution model refers to the original paper as cited in the paper.
4. fit_daily_model.R; The R code used to fit the proposed daily pollution model.
5. missingDesign.R; The R code used to take out validation data.

