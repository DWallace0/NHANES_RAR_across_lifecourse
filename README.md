# NHANES_RAR_across_lifecourse
Code for the paper: Rest-activity rhythms across the lifespan: cross-sectional findings from the U.S.-representative National Health and Nutrition Examination Survey (NHANES) 

In addition to demographic and physical exam information, this analysis used data from the physical activity monitor (PAMS) component of the 2011-2012 and 2013-2104 NHANES exams. Minute-epoch and header data for this data was used:

2011-2012 NHANES data: https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?BeginYear=2011
Minute-epoch 2011 PAMS actigraphy data: https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/PAXMIN_G.XPT
Minute-epoch 2011 PAMS actigraphy data explainer: https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/PAXMIN_G.htm
2011 actigraphy header file: https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/PAXHD_G.XPT

2013-2014 NHANES data: https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?BeginYear=2013
Minute-epoch 2013 PAMS actigraphy data: https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/PAXMIN_H.XPT
Minute-epoch 2011 PAMS actigraphy data explainer: https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/PAXMIN_H.htm
2013 actigraphy header file: https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/PAXHD_H.XPT

The provided analysis code is divided as follows for the different data processing and analysis steps:

Data download, pre-processing, imputation with “accelmissing”, and creation of RAR metrics using “nparACT” – code to do this has the filename “0_NparACT_accelmissing.R”; the resulting data has the following names:

- “accelmissing_out2” files = contains data as nested lists from raw timeseries data, formatted for accelmissing package

- “acceldata” files = activity data (PA) combined with day of week info (label), wearing / offwrist information (flag), and demographic info (age, sex, race/ethnicity) (demo) together in a dataframe

- “actdata2” files = the results after filtering the “acceldata” file through accelmissing to create a missing flag matrix and provide some information on amount missingness; data filtered for quality / number of valid days. This is the file used for imputation. 

- “sensitivity_IDs” file = subset of IDs using more stringent criteria where valid day defined as at least 20 hours of wearing and 6+ valid days 

- “accelimp” file = the imputed datasets after running accelmissing; these data will be fed into nparact to derive the RAR metrics
