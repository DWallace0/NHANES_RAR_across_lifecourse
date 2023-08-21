#########################################################
# 2023 - Danielle A. Wallace
#Rest-activity rhythms across the lifespan: 
#cross-sectional findings from the U.S.-representative 
#National Health and Nutrition Examination Survey (NHANES) 
#########################################################
#code
#########################################################
#STEP 0: Data download and pre-processing code 
#########################################################
#*******IMPORTANT*******
#The minute-epoch data on the NHANES website is provided as a single file
#in long format; that is, all IDs together are bound by row and in one large file
#for data processing, this file was downloaded to an external USB (E drive in the code below)
#and then split by SEQN (participant ID);
#a separate folder with a .csv file for each individual's data was created
#each .csv file in the folder was then looped over

#example of reading in the 2011 minute-epoch data that has been saved to a USB (E: drive)
library(haven)
act <- read_xpt("E:/PAXMIN_G.XPT") 
#example of splitting by ID and writing to individual files
for(i in unique(act$SEQN)) {
  ID <- subset(act, SEQN == i)
  write.csv(ID, file = paste0("E:/individual_files_2011/",i,".csv"))
}
#The following code assumes that you have already:
#1. downloaded the NHANES minute-epoch file and saved it somewhere (like a USB)
#2. have split the minute-epoch file by SEQN and saved each as an individual csv file in a folder
#########################################################
#set your working directory
setwd("/your path here")
#load packages
library(accelmissing)
library(fs)
library(gmodels)
library(tidyverse)
library(readxl)
library(ggplot2)
library(reshape2)
library(data.table)
library(tidyr)
library(dplyr)
library(stringr)
library(purrr)
library(haven)
library(lubridate)
library(stats)
library(zoo)
library(xfun)
library(nparACT)
#########################################################
#LOOP CODE 
####################################################
#for 2011 data
###########################################################
#load the header data; this has starttime and stoptime time stamp data
# we can use this to fill in time so that each row has a corresponding time stamp
header11 <- read_xpt("https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/PAXHD_G.XPT") 
#2011 data has n=7,821 rows, but n=904 people missing PAM data 
#provide the path info for where the individual csv files that you already created are stored:
#NOTE: replace location in "path=" to your path 
filenames <- list.files(path="E:/individual_files_2011", pattern=".*csv") 
length(filenames) #[1] n=6917 individual files
act_file <- paste0("E:/individual_files_2011/",filenames)
###########################################################
#create list objects to store output
out2<-list() 
out3<-list() 
out4<-list() 
#loop over csv files:
for(chr in act_file) {
  nhanes_act <- read.csv(chr)
  try({ 
    #if wanting to match header timestamp to individual long format data:
    nhanes_act$start_time <- header11$PAXFTIME[match(nhanes_act$SEQN, header11$SEQN)]
    nhanes_act$end_time <- header11$PAXETLDY[match(nhanes_act$SEQN, header11$SEQN)]
    start <- nhanes_act[1,18]
    nhanes_act$time <- NA
    nhanes_act$date <- NA
    #set first timestamp
    nhanes_act[1, c("time")] <- start
    nhanes_act[1, c("date")] <- "01-01-2001" #make up a date since info not available
    #create a date+time variable by combining and format to UNIX timestamp
    nhanes_act$datetime <- paste(nhanes_act$date, nhanes_act$time, sep=" ")
    nhanes_act$datetime[nhanes_act$datetime == "NA NA"] <- NA
    nhanes_act$datetime2<- parse_date_time(nhanes_act$datetime, '%m/%d/%Y %H:%M:%S')
    #create the start and end date values
    start_date <- as.Date(nhanes_act[1,23])
    end_date <- start_date+(max(nhanes_act$PAXDAYM)-1)
    end_date2 <- paste(end_date, nhanes_act[1,19], sep=" ")
    
    start2 <- as.POSIXlt(nhanes_act[1,23],format='%Y-%m-%d %H:%M:%S', tz="GMT")
    end2 <- as.POSIXlt(end_date2,format='%Y-%m-%d %H:%M:%S', tz="GMT")
    #for those with differing number of rows (NHANES documentation), create row time based on # rows in dataset
    rows <- as.numeric(nrow(nhanes_act))
    end3 <- as.POSIXlt(format(strptime(start2, format = "%Y-%m-%d %H:%M:%S") + (rows-1)*60, "%Y-%m-%d %H:%M:%S"), tz="GMT")
    #fill in 
    nhanes_act$time <- try(seq(start2, end3, by = "mins")) #may be missing rows of time (NHANES documentation)
    #for these, change end time to be based on row numbers - 1 min 
    if("try-error" %in% class(nhanes_act$time)) alternativeFunction(seq(start2, end3, by = "mins"))
   #Incorporate some data quality flags:
    #evaluate days available and data quality; select last day where data quality is good for # valid days
    nhanes_act$num_days <- max(nhanes_act$PAXDAYM)
    if(nhanes_act$num_days<8) next # if less than 8 days (<=6 full days), skip to next iteration
    ###################################
     #first, subset to only select from day 2-day 8
    nhanes_t<-subset(nhanes_act, (PAXDAYM>1& PAXDAYM<9))
    #get a measure of rows in case missing
    nhanes_t$nrows <- nrow(nhanes_t)
    #Now, use the data quality flags that CDC developed; have letters to indicate different issues
    #if flag has a letter (not NA or empty), make 1 (to indicate flagged quality)
    to_match <- paste(LETTERS, collapse = "|") #this creates letter values
    #if there are any letters in the variable, make=1
    nhanes_t$flag <- ifelse(grepl(to_match, nhanes_t$PAXFLGSM), 1, 0) 
    #create a missingness variable from the estimated non-wear variable, PAXPREDM
    nhanes_t$missing <- ifelse(nhanes_t$PAXPREDM==3, 1, 0)
    #flag combining the two prior flags
    nhanes_t$dual_flag <- ifelse(nhanes_t$flag==1 | nhanes_t$missing==1, 1, 0)
    #Recode activity data to missing if CDC quality flag is triggered or 
    #if the algorithm detects non-wear or no activity measure
    #if PAXMTSM value is -0.01, it means the activity could not be computed; recode to missing; else keep the same
    nhanes_t$activity_na <- ifelse(nhanes_t$dual_flag==1 | nhanes_t$PAXMTSM<0, NA, nhanes_t$PAXMTSM)
    nhanes_t$activity_na <- ifelse(nhanes_t$activity_na>=0, as.integer(nhanes_t$activity_na), NA)
    #needs to be integer:
    nhanes_t$offwrist <- as.integer(ifelse((is.na(nhanes_t$activity_na)), 1, 0))
    nhanes_t$day <- ((as.numeric(nhanes_t$PAXDAYM))-1)
    nhanes_t$dayofweek <- as.numeric(nhanes_t$PAXDAYWM)
    #convert row by day
    nhanes_t$hms <- substring(nhanes_t$time, 12)
    #############################################
    #output variables of interest
    #format label - condense so that each day has its own row
    label <- nhanes_t[c("SEQN", "day", "dayofweek")]
    label$rowvars <- paste0(nhanes_t$SEQN,".",nhanes_t$day)
    label2 <-  label %>% group_by(rowvars) %>% slice(1)
    label2$rowvars<- NULL
    #create an offwrist dataset that matches dimensions of light dataset
    nhanes_t3 <- nhanes_t[c("SEQN", "offwrist", "day", "hms")]
    offwrist <- reshape(nhanes_t3, idvar ="hms",v.names="offwrist", timevar = "day", direction = "wide")
    rownames(offwrist) <- NULL
    offwrist$SEQN <- NULL
    offwrist$hms <- NULL
    offwrist2 <- t(offwrist)
    rownames(offwrist2) <- NULL #so that rownames will match the other datasets
    
    #will also impute activity, so create activity dataset
    nhanes_t4 <- nhanes_t[c("SEQN", "activity_na", "day", "hms")]
    activity <- reshape(nhanes_t4, idvar ="hms",v.names="activity_na", timevar = "day", direction = "wide")
    rownames(activity) <- NULL
    activity$SEQN<- NULL
    activity$hms <- NULL
    activity2 <- t(activity)
    rownames(activity2) <- NULL #so that rownames will match the other datasets
    
  })#end brackets for "try" command
  #############
  out2[[chr]] <- label2
  out3[[chr]] <- offwrist2
  out4[[chr]] <- activity2
}
#Warning messages:
#1: In if (nhanes_act$num_days < 7) next :
#  the condition has length > 1 and only the first element will be used
####################
#save output
save(out2, out3, out4, file="accelmissing_2011_out2.RData") #6676
#########################################################
#format data
####################################################
load("accelmissing_2011_out2.RData")
label <- do.call(rbind, lapply(out2, data.frame))
offwrist <- do.call(rbind, lapply(out3, data.frame))
activity <- do.call(rbind, lapply(out4, data.frame))
#read in NHANES 2011 demographic data for imputation
#########################
demo11 <- read_xpt(url("https://wwwn.cdc.gov/Nchs/Nhanes/2011-2012/DEMO_G.XPT")) #demographics
demo11$sex <- ifelse(demo11$RIAGENDR==1, 1, 0) #male=1, female=0
demo11$race_eth <- as.factor(ifelse(demo11$RIDRETH3==1, "Mexican American",
                                    ifelse(demo11$RIDRETH3==2, "Other Hispanic",
                                           ifelse(demo11$RIDRETH3==3, "NH White",
                                                  ifelse(demo11$RIDRETH3==4, "NH Black",
                                                         ifelse(demo11$RIDRETH3==6, "NH Asian",
                                                                ifelse(demo11$RIDRETH3==7, "Other/Multiracial", NA)))))))
names(demo11)
#can't handle missing values
vars <- c("SEQN", "RIDAGEYR", "sex", "race_eth")
demo <- demo11[vars]
idv <- unique(label$SEQN) 
demo1 <- as.data.frame(demo[demo$SEQN %in% as.character(idv),])#there are 6676 from demo 
###########################
label$day <- as.integer(label$day)
label$dayofweek <- as.integer(label$dayofweek)
label$day <- NULL #getting rid of for now; only supposed to be N x 2 matrix with ID and dayofweek
#create list 
#save activity data
actdata <- list(activity, label, offwrist, demo1)
names(actdata) <- c("PA", "label", "flag", "demo") #flag=1 if missing
save(actdata, file="acceldata_11.RData") #2011 data
is.data.frame(actdata$PA)
#########################################################
#evaluate missingness and create missing flag matrix
######################################################
#accelmissing cannot handle NAs
#all NA values changed to 0 and create a flag matrix using a window of 60 minutes 
load("acceldata_11.RData") #2011 data
#replace NA values with 0
actdata$PA[is.na(actdata$PA)] <- 0
actdata$demo #demo data has age, sex, race/ethnicity
## create missing flag matrix by 60 min criterion
flag60 = create.flag(PA=actdata$PA, window=60, mark.missing=1)
mr = missing.rate(label=actdata$label, flag60, mark.missing=1, time.range=c("00:00", "23:59"))
mr$total  
#2011 = 11.1% for midnight-midnight when missing window defined as 60 min 0 counts
## missing proportion by days
mean(mr$table < 0.1) 
#PLOT wearing proportion over time 
wear.time.plot(PA=actdata$PA, label=actdata$label, flag=flag60, mark.missing=1)
# data filtering for valid days
valid.days.out = valid.days(PA=actdata$PA, label=actdata$label, flag=flag60, mark.missing=1, wear.hr=16, time.range = c("00:00", "23:59"))
ls(valid.days.out) # list with three matrix objects
# data filtering for valid subjects
x1 = list(PA=actdata$PA, label=actdata$label, flag=flag60, demo=actdata$demo) # original
x2 = valid.days.out # output of valid.days()
valid.sub.out = valid.subjects(data1=x1, data2=x2, valid.days=4)
#2011 = 6022 out of 6676 originally
length(unique(valid.sub.out$label[,1])) 
ls(valid.sub.out)
## missing rate with the filtered data
missing.rate(valid.sub.out$label, valid.sub.out$flag, mark.missing = 1, time.range = c("00:00", "23:59"))$total
# 2011 = 5.5%
# demographic data for the filtered data
idv= unique(valid.sub.out$label[,1])
matchid <- idv[idv %in% actdata$demo$SEQN]
demo1 <- as.data.frame(actdata$demo[actdata$demo$SEQN %in% as.character(idv),])
# save the data before imputation
actdata2 = list(PA=valid.sub.out$PA, label=valid.sub.out$label, flag=valid.sub.out$flag,
                demo=demo1)
save(actdata2, file="actdata2_11.RData") #2011 data
#################################################
#now let's also get a list of IDs using the stricter criteria that we will use for the sensitivity analysis
#changing valid day = 20 hours and weartime = 6 days or more
valid.days.out.x = valid.days(PA=actdata$PA, label=actdata$label, flag=flag60, mark.missing=1, wear.hr=20, time.range = c("00:00", "23:59"))
ls(valid.days.out.x) # list with three matrix objects
x1.x = list(PA=actdata$PA, label=actdata$label, flag=flag60, demo=actdata$demo) # original
x2.x = valid.days.out.x # output of valid.days()
valid.sub.out.x = valid.subjects(data1=x1.x, data2=x2.x, valid.days=6)
#2011 = 5123 out of 6676;  2013 = 5512 persons out of 7580 originally
length(unique(valid.sub.out.x$label[,1])) 
ls(valid.sub.out.x)
## missing rate with the filtered data
missing.rate(valid.sub.out.x$label, valid.sub.out.x$flag, mark.missing = 1, time.range = c("00:00", "23:59"))$total
# 2011 = 3% missing 
# demographic data for the filtered data
idv_2011= unique(valid.sub.out.x$label[,1])

sensitivity_IDs <- c(idv_2011, idv_2013)
#save(sensitivity_IDs, file="sensitivity_IDs.RData") #2011 data
#########################################################
#perform the activity imputation using accelmissing
######################################################
load("actdata2_11.RData") 
# prepare the imputation
library(mice); library(pscl)
data = actdata2
demo = as.data.frame(data$demo)
names(demo) <- c("id", "age", "sex", "race_eth")
# imputation: this takes a long time / computing power; best to run in cloud
accelimp = accel.impute(PA=data$PA, label=data$label, flag=data$flag,
                        demo=demo, time.range=c("06:00","23:00"), method="zipln.pmm",
                        K = 3, D = 5, mark.missing = 1, thresh = 600, graph.diagnostic = TRUE,
                        seed = 1234, m = 5, maxit = 6)
save(accelimp, file="accelimp_2011.RData")
###########################################################################
#using the imputed data (5 sets), derive the RAR measures with nparact
##########################################################################
library(nparACT)
library(dplyr)
library(tidyr)
load("accelimp_2011.RData")
#separate out imputed datasets to run nparact
imp1 <- as.data.frame(accelimp$imp1)
imp2 <- as.data.frame(accelimp$imp2)
imp3 <- as.data.frame(accelimp$imp3)
imp4 <- as.data.frame(accelimp$imp4)
imp5 <- as.data.frame(accelimp$imp5)
#now, to run Nparact, need to create loop to split datasets into single actigraphy files
#add a timestamp variable
#and then run nparact function
imp1$ID <- rownames(imp1)
imp1$ID = substr(imp1$ID,26,nchar(imp1$ID)-6)
head(imp1$ID)
imp1$day <- substr((rownames(imp1)),36, nchar(rownames(imp1)))
head(imp1$day)

ids <- unique(imp1$ID)

for(chr in ids) {
  try({ #try command to skip i and go to next in loop if error occurs
    #first subset by ID and create new dataframe
    imp1_sub <- imp1[which(imp1$ID==chr),]
    imp1_sub$ID <- NULL
    imp1_sub$day <- NULL
    imp1_long <- imp1_sub %>% pivot_longer(cols = colnames(imp1_sub), names_to = 'min', values_to = 'activity')
    #create time stamp
    imp1_long$Time <- format(seq(as.POSIXct("2001-01-01 00:00:00", tz="GMT"), 
                                 length.out=10080, by='1 min'), '%Y-%m-%d %H:%M:%S')
    imp1_long$Activity <- imp1_long$activity
    imp1_long$activity <- NULL
    imp1_long$min<- NULL
    ##########
    npar_a <- data.frame(imp1_long)
    # nparact function:
    npar_a_output11 <- nparACT_base("npar_a", 1/60, cutoff = 1, plot = F, fulldays = F)
    # Write RAR metrics data (from imputed activity data) to a new excel file:
    write.table(npar_a_output11, "npar_a_Results_imp1_11.csv", append = TRUE, sep = ",", row.names = paste(chr), col.names = FALSE)
  })
}

#*IMPORTANT*
#*repeat this code for additional sets of files (alter naming as necessary) - imp2, imp3, imp4, and imp5 - and save results
#*repeat steps for 2013 data
###########################################
