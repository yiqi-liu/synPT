library(dplyr)
library(tidyr)
library(foreign)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
### -------------------------------------------------------------
# The following code for data cleaning is adapted from
# https://www.openicpsr.org/openicpsr/project/146381/version/V1/view?path=/openicpsr/146381/fcr:versions/V1/synthdid-sdid-paper
# Download folder 'synthdid-sdid-paper' from the above URL and place it in the same directory as the current file. The code here is based on 'synthdid-sdid-paper/inst/extdata/placebo/cps.R'.
### -------------------------------------------------------------
source('synthdid-sdid-paper/inst/extdata/placebo/functions.R') # function used to prep the CPS aggregate data

url_place <- "https://data.nber.org/morg/annual/" # where to download CPS data

# year indices, 1979-2018, total 40 years
seq_1 <- 79:99  # 1979-1999
seq_2 <- 0:9    # 2000-2009
seq_3 <- 10:18  # 2010-2018

part_1 <- paste(paste('morg',seq_1,sep = ''),'dta',sep = '.')
part_2 <- paste(paste('morg0',seq_2,sep = ''),'dta',sep = '.')
part_3 <- paste(paste('morg',seq_3,sep = ''),'dta',sep = '.')

# Downloading the files 

fi_1 <- paste(url_place,part_1,sep = '')
# dat_1 <- lapply(fi_1,read.dta) # file too large and R throws error
files  <- file.path(tempdir(), part_1) # chatGPT solution: download first before reading
mapply(curl::curl_download, fi_1, files)
dat_1 <- lapply(files, read.dta)

fi_2 <- paste(url_place,part_2,sep = '')
# dat_2 <- lapply(fi_2,read.dta)
files  <- file.path(tempdir(), part_2)
mapply(curl::curl_download, fi_2, files)
dat_2 <- lapply(files, read.dta)

fi_3 <- paste(url_place,part_3,sep = '')
# dat_3 <- lapply(fi_3,read.dta)
files  <- file.path(tempdir(), part_3)
mapply(curl::curl_download, fi_3, files)
dat_3 <- lapply(files, read.dta)

data_assign <- as.matrix(read.table('synthdid-sdid-paper/inst/extdata/placebo/state_laws.tsv')) ==1 # auxiliary outcome data used to simulate treatment assignment

# Subsetting the data
## dictionary here: https://data.nber.org/morg/docs/cpsbjan03.ddf
index_full_1 <- c('hhid','hhnum','intmonth','age','sex','minsamp',
                  'state','year',"lineno",'earnwke','uhourse',
                  "ftpt94","ftpt89","ftpt79")
index_full_2 <- c('hhid','hhnum','intmonth','age','sex','minsamp',
                  'stfips','year',"lineno",'earnwke','uhourse',
                  "ftpt94","ftpt89","ftpt79")
index_part <- c("ftpt94","ftpt89","ftpt79")

sort_index_1 <- c("hhid", "hhnum", "lineno", "year", "minsamp", "intmonth", "state", "age")
sort_index_2 <- c("hhid", "hhnum", "lineno", "year", "minsamp", "intmonth", 'stfips', "age")

# subset to "females from 25 to 50, 4th month of interview and not from DC"
dat_1_sub <- data_const(dat_1,'state',index_full_1,sort_index_1,index_part)
dat_2_sub <- data_const(dat_2,'stfips',index_full_2,sort_index_2,index_part)
dat_3_sub <- data_const(dat_3,'stfips',index_full_2,sort_index_2,index_part)

# Average log wage
lwage_1 <- earn_const(dat_1_sub,'state')
lwage_2 <- earn_const(dat_2_sub,'stfips')
lwage_3 <- earn_const(dat_3_sub,'stfips')

data_mat_wage <- cbind(lwage_1, lwage_2, lwage_3)
colnames(data_mat_wage) <- 1979:2018

### -------------------------
# function that generates repeated cross-sectional data from the data used to construct the CPS aggregate data
gen_rc <- function(data_list, state_id){
  
  rc <- do.call( # for each year of cross sectional data in data_list, then rbind across all years
    rbind, lapply(data_list, function(dta){
      # the individual-level outcome: log wage
      ind_earn <- dta[, 'earnwke'] > 0 & !is.na(dta[, 'earnwke'])
      outcome <- round(log(dta[ind_earn, 'earnwke']), 4)
      # extract state-year corresponding to positive earnings
      state <- dta[ind_earn, state_id]
      state_num <- as.numeric(factor(as.character(state),
                                     levels=sort(unique(as.character(state)))))
      year <- dta[ind_earn, 'year']
      # get the cross-section of this year
      cross_sec <- data.frame(state_num, state, year, outcome)
      return(cross_sec)})
  )
  
  return(rc)
}
### -------------------------

# construct repeated cross-sections
rc_1 <- gen_rc(dat_1_sub, 'state')
rc_2 <- gen_rc(dat_2_sub, 'stfips')
rc_3 <- gen_rc(dat_3_sub, 'stfips')
rcs_CPS <- bind_rows(rc_1, rc_2, rc_3)
write.csv(rcs_CPS, "cps_rep_cross_sec.csv", row.names=FALSE)

## SANITY CHECK DATA CLEANING IS CORRECT
library(synthdid)
data(CPS) # load the CPS data provided by the synthdid package
CPS$state_num <- as.numeric(CPS$state)
# get (kt)-cell means and sd
cell_stats <- rcs_CPS %>%
  group_by(state_num, year) %>%
  summarise(
    n_kt=n(),
    mu=mean(outcome, na.rm = TRUE),
    sd=sd(outcome, na.rm = TRUE),
    .groups = "drop"
  )
cell_stats <- left_join(cell_stats, CPS[, c("state_num", "year", "log_wage")], by=c("state_num", "year"))

### --------- comment ---------
# there's some minor mismatch between the CPS aggregate data constructed from 'synthdid-sdid-paper/inst/extdata/placebo/cps.R'---which can also be found in their replication package 'synthdid-sdid-paper/data/CPS.csv'---and directly calling 'data(CPS)' from the synthdid package. The README.md file from 'synthdid-sdid-paper/data/README.md' says the "source data, as well as code that processes it into the csv files included here, is included in [inst/extdata](inst/extdata)," so there should be an exact match. I don't know what causes the mismatch, but I'll use the cell means I get from the original repeated cross-section to do the simulation.
### --------- comment ---------

# CPS data from the replication package 'synthdid-sdid-paper/data/CPS.csv'
CPS_rep_pkg <- readr::read_delim("synthdid-sdid-paper/data/CPS.csv", delim = ";")
CPS_rep_pkg$state_num <- as.numeric(factor(CPS_rep_pkg$state, levels = sort(unique(CPS_rep_pkg$state))))
compare_repCPSvsRCS <- left_join(cell_stats[, 1:5], CPS_rep_pkg[, c("state_num", "year", "log_wage")], by=c("state_num", "year"))

# data(CPS) doesn't match the means from repeated cross sections
quantile(round(cell_stats$mu-cell_stats$log_wage, 4))

# CPS.csv from replication package does match the means from repeated cross sections
quantile(round(compare_repCPSvsRCS$mu-compare_repCPSvsRCS$log_wage, 4))