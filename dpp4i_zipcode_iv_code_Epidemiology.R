# Title: Comparing area-level patient density and physician prescribing preference instruments
# for the effect of antidiabetics on adverse cardiovascular events among Medicare beneficiaries
# Expanded Physician Preference Instrument using Area-level Data
# Author: Jack Cordes, MS

# Libraries ####
library(tidyverse)
library(tableone)
library(pROC)
library(sf)
library(spdep)
library(rgdal)
library(cowplot)
library(survival)
library(lme4)
library(lmtest)
library(sandwich)
library(sbw)
library(quadprog)
library(OneSampleMR)
library(ivreg)
library(sp)
library(st)
library(maptools)
library(rgeos)
library(geojsonio)
library(rmapshaper)
library(scales)
library(gtable)
library(grid)
library(car)
library(ggfortify)
library(survminer)

# Zip Code Data ####

# List of pertinent ACS variable codes 2012
acs_codes_2012_1 <- c("GISJOIN", "YEAR", "ZCTA5A", "QSEE001", "QUUE046", "QSQE002", "QSQE003",
                      "QSQE004", "QSQE005", "QSQE006", "QSQE007", "QSQE008", "QSYE003", "QSYE004",
                      "QSYE005", "QSYE006", "QSYE012", "QSYE013", "QSYE014", "QSYE015", "QSYE016",
                      "QUUE047", "QUUE049", "QUUE050", "QUUE052", "QUUE054", "QUUE055", "QUUE057",
                      "QUUE059", "QUUE060", "QUUE062", "QUUE064", "QUUE065", "QUUE067", "QU1E001",
                      "QVME005", "QVUE001", "QVUE002", "QWUE001", "QXHE017", "QXHE020", "QXHE035", "QXHE038")

acs_codes_2012_2 <- c("ZCTA5A", "Q2ZE004", "Q2ZE009", "Q2ZE015", "Q2ZE020", "Q7ME015", "Q7ME016",
                      "Q7ME017", "Q7ME108", "Q7ME109", "Q7ME110", "Q7ME031", "Q7ME032", "Q7ME033",
                      "Q7ME047", "Q7ME048", "Q7ME049", "Q7ME062", "Q7ME063", "Q7ME064", "Q7ME124",
                      "Q7ME125", "Q7ME126", "Q7ME140", "Q7ME141", "Q7ME142", "Q7ME155", "Q7ME156",
                      "Q7ME157", "Q7ME092", "Q7ME093", "Q7ME094", "Q7ME185", "Q7ME186", "Q7ME187",
                      "Q7ME077", "Q7ME078", "Q7ME079", "Q7ME170", "Q7ME171", "Q7ME172", "Q8ZE036",
                      "Q8ZE037", "Q8ZE038", "Q8ZE077", "Q8ZE078", "Q8ZE079", "Q8ZE039", "Q8ZE040",
                      "Q8ZE080", "Q8ZE081", "Q8ZE041", "Q8ZE042", "Q8ZE082", "Q8ZE083", "Q9QE015",
                      "Q9QE016", "Q9QE029", "Q9QE030", "RHLE017", "RA1E016", "RA1E019", "RA1E035", 
                      "RA1E038", "RBDE013", "RBDE016", "RBDE029", "RBDE032", "RBEE013", "RBEE016",
                      "RBEE029", "RBEE032", "RCBE001", "RHLE018", "RHLE019", "RHLE020", "RHLE021")

# 5-year ACS 2008-2012
zipcode_2012_1 <- read.csv("/Volumes/PACS$/Cordes, Jack/Zip Code Data/nhgis0041_csv/nhgis0041_ds191_20125_zcta.csv")
zipcode_2012_1 <- read.csv("/Volumes/PACS$/Cordes, Jack/Zip Code Data/nhgis0041_csv/nhgis0041_ds191_20125_zcta.csv")[2:nrow(zipcode_2012_1), acs_codes_2012_1]
zipcode_2012_1[, 4:length(zipcode_2012_1)] <- sapply(zipcode_2012_1[, 4:length(zipcode_2012_1)], as.numeric)

zipcode_2012_2 <- read.csv("/Volumes/PACS$/Cordes, Jack/Zip Code Data/nhgis0041_csv/nhgis0041_ds192_20125_zcta.csv")
zipcode_2012_2 <- read.csv("/Volumes/PACS$/Cordes, Jack/Zip Code Data/nhgis0041_csv/nhgis0041_ds192_20125_zcta.csv")[2:nrow(zipcode_2012_2), acs_codes_2012_2]
zipcode_2012_2[, 2:length(zipcode_2012_2)] <- sapply(zipcode_2012_2[, 2:length(zipcode_2012_2)], as.numeric)

zipcode_2012 <- left_join(zipcode_2012_1, zipcode_2012_2, by = "ZCTA5A")

zipcode_2012$PROP_WHITE <- zipcode_2012$QSQE002/zipcode_2012$QSEE001
zipcode_2012$PROP_BLACK <- zipcode_2012$QSQE003/zipcode_2012$QSEE001
zipcode_2012$PROP_AIAN <- zipcode_2012$QSQE004/zipcode_2012$QSEE001
zipcode_2012$PROP_ASIAN <- zipcode_2012$QSQE005/zipcode_2012$QSEE001
zipcode_2012$PROP_NHP <- zipcode_2012$QSQE006/zipcode_2012$QSEE001
zipcode_2012$PROP_OTHER_RACE <- zipcode_2012$QSQE007/zipcode_2012$QSEE001
zipcode_2012$PROP_2MORE_RACE <- zipcode_2012$QSQE008/zipcode_2012$QSEE001

zipcode_2012$PROP_NONHISP_WHITE <- zipcode_2012$QSYE003/zipcode_2012$QSEE001
zipcode_2012$PROP_NONHISP_BLACK <- zipcode_2012$QSYE004/zipcode_2012$QSEE001
zipcode_2012$PROP_NONHISP_AIAN <- zipcode_2012$QSYE005/zipcode_2012$QSEE001
zipcode_2012$PROP_NONHISP_ASIAN <- zipcode_2012$QSYE006/zipcode_2012$QSEE001
zipcode_2012$PROP_HISP <- zipcode_2012$QSYE012/zipcode_2012$QSEE001
zipcode_2012$PROP_HISP_WHITE <- zipcode_2012$QSYE013/zipcode_2012$QSEE001
zipcode_2012$PROP_HISP_BLACK <- zipcode_2012$QSYE014/zipcode_2012$QSEE001
zipcode_2012$PROP_HISP_AIAN <- zipcode_2012$QSYE015/zipcode_2012$QSEE001
zipcode_2012$PROP_HISP_ASIAN <- zipcode_2012$QSYE016/zipcode_2012$QSEE001

zipcode_2012$PROP_ENGLISH_ONLY <- zipcode_2012$QUUE047/zipcode_2012$QUUE046
zipcode_2012$PROP_ENGLISH_WELL <- (zipcode_2012$QUUE049 + zipcode_2012$QUUE050 + 
                                     zipcode_2012$QUUE054 + zipcode_2012$QUUE055 + 
                                     zipcode_2012$QUUE059 + zipcode_2012$QUUE060 + 
                                     zipcode_2012$QUUE064 + zipcode_2012$QUUE065)/zipcode_2012$QUUE046
zipcode_2012$PROP_NO_ENGLISH <- (zipcode_2012$QUUE052 + zipcode_2012$QUUE057 +
                                   zipcode_2012$QUUE062 + zipcode_2012$QUUE067)/zipcode_2012$QUUE046

zipcode_2012$PROP_PUBLIC_ASST <- zipcode_2012$QVUE002/zipcode_2012$QVUE001

zipcode_2012$PROP_VET <- (zipcode_2012$QXHE017 + zipcode_2012$QXHE020 + zipcode_2012$QXHE035 + 
                            zipcode_2012$QXHE038)/zipcode_2012$QUUE046

zipcode_2012$PROP_NATIVE <- (zipcode_2012$Q2ZE004 + zipcode_2012$Q2ZE009 + zipcode_2012$Q2ZE015 +
                               zipcode_2012$Q2ZE020)/zipcode_2012$QSEE001

zipcode_2012$PROP_NEVER_MARRIED <- (zipcode_2012$Q7ME015 + zipcode_2012$Q7ME016 + zipcode_2012$Q7ME017 +
                                      zipcode_2012$Q7ME108 + zipcode_2012$Q7ME109 + zipcode_2012$Q7ME110)/zipcode_2012$QUUE046

zipcode_2012$PROP_MARRIED <- (zipcode_2012$Q7ME031 + zipcode_2012$Q7ME032 + zipcode_2012$Q7ME033 + 
                                zipcode_2012$Q7ME047 + zipcode_2012$Q7ME048 + zipcode_2012$Q7ME049 + 
                                zipcode_2012$Q7ME062 + zipcode_2012$Q7ME063 + zipcode_2012$Q7ME064 + 
                                zipcode_2012$Q7ME124 + zipcode_2012$Q7ME125 + zipcode_2012$Q7ME126 + 
                                zipcode_2012$Q7ME140 + zipcode_2012$Q7ME141 + zipcode_2012$Q7ME142 + 
                                zipcode_2012$Q7ME155 + zipcode_2012$Q7ME156 + zipcode_2012$Q7ME157)/zipcode_2012$QUUE046

zipcode_2012$PROP_DIVORCED <- (zipcode_2012$Q7ME092 + zipcode_2012$Q7ME093 + zipcode_2012$Q7ME094 + 
                                 zipcode_2012$Q7ME185 + zipcode_2012$Q7ME186 + zipcode_2012$Q7ME187)/zipcode_2012$QUUE046

zipcode_2012$PROP_WIDOWED <- (zipcode_2012$Q7ME077 + zipcode_2012$Q7ME078 + zipcode_2012$Q7ME079 + 
                                zipcode_2012$Q7ME170 + zipcode_2012$Q7ME171 + zipcode_2012$Q7ME172)/zipcode_2012$QUUE046

zipcode_2012$PROP_HS <- (zipcode_2012$Q8ZE036 + zipcode_2012$Q8ZE037 + zipcode_2012$Q8ZE038 + 
                           zipcode_2012$Q8ZE077 + zipcode_2012$Q8ZE078 + zipcode_2012$Q8ZE079)/zipcode_2012$QUUE046
zipcode_2012$PROP_ASSOC <- (zipcode_2012$Q8ZE039 + zipcode_2012$Q8ZE040 + zipcode_2012$Q8ZE080 +
                              zipcode_2012$Q8ZE081)/zipcode_2012$QUUE046
zipcode_2012$PROP_BACHELORS <- (zipcode_2012$Q8ZE041 + zipcode_2012$Q8ZE042 + zipcode_2012$Q8ZE082 +
                                  zipcode_2012$Q8ZE083)/zipcode_2012$QUUE046

zipcode_2012$PROP_POVERTY <- (zipcode_2012$Q9QE015 + zipcode_2012$Q9QE016 + zipcode_2012$Q9QE029 + 
                                zipcode_2012$Q9QE030)/zipcode_2012$QUUE046

zipcode_2012$PROP_DISABILITY <- (zipcode_2012$RA1E016 + zipcode_2012$RA1E019 + zipcode_2012$RA1E035 + 
                                   zipcode_2012$RA1E038)/zipcode_2012$RHLE017
zipcode_2012$PROP_COGNITIVE <- (zipcode_2012$RBDE013 + zipcode_2012$RBDE016 + zipcode_2012$RBDE029 + 
                                  zipcode_2012$RBDE032)/zipcode_2012$RHLE017
zipcode_2012$PROP_AMBULATORY <- (zipcode_2012$RBEE013 + zipcode_2012$RBEE016 + zipcode_2012$RBEE029 + 
                                   zipcode_2012$RBEE032)/zipcode_2012$RHLE017

# zipcode_2012$PROP_PRIVATE_INS <- (zipcode_2012$RHLE018)/zipcode_2012$RHLE017
# zipcode_2012$PROP_PUBLIC_INS <- (zipcode_2012$RHLE019)/zipcode_2012$RHLE017
# zipcode_2012$PROP_PRIVATE_PUBLIC_INS <- (zipcode_2012$RHLE020)/zipcode_2012$RHLE017
# zipcode_2012$PROP_NO_INS <- (zipcode_2012$RHLE021)/zipcode_2012$RHLE017
zipcode_2012$PROP_MEDICARE <- (zipcode_2012$RHLE019 + zipcode_2012$RHLE020)/zipcode_2012$RHLE017

zipcode_2012$MED_INC <- zipcode_2012$QU1E001
zipcode_2012$MED_INC_65 <- zipcode_2012$QVME005
zipcode_2012$GINI <- zipcode_2012$RCBE001

# List of pertinent ACS variable codes 2013
acs_codes_2013_1 <- c("GISJOIN", "YEAR", "ZCTA5A", "UEEE001", "UG6E046", "UEQE002", "UEQE003",
                      "UEQE004", "UEQE005", "UEQE006", "UEQE007", "UEQE008", "UEYE003", "UEYE004",
                      "UEYE005", "UEYE006", "UEYE012", "UEYE013", "UEYE014", "UEYE015", "UEYE016",
                      "UG6E047", "UG6E049", "UG6E050", "UG6E052", "UG6E054", "UG6E055", "UG6E057",
                      "UG6E059", "UG6E062", "UG6E060", "UG6E067", "UG6E064", "UG6E065", "UHDE001",
                      "UHYE005", "UH6E001", "UH6E002", "UJAE001", "UJXM017", "UJXM020", "UJXM035",
                      "UJXM038", "UILE001", "UM2E055", "UM2E060", "UM2E061", "UM2E062", "UM2E051")

acs_codes_2013_2 <- c("ZCTA5A", "UPIE004","UPIE009","UPIE015","UPIE020", "UT5E015","UT5E016",
                      "UT5E017","UT5E108","UT5E109","UT5E110","UT5E031","UT5E032","UT5E033",
                      "UT5E047","UT5E048","UT5E049","UT5E062","UT5E063","UT5E064","UT5E124",
                      "UT5E125","UT5E126","UT5E140","UT5E141","UT5E142","UT5E155","UT5E156",
                      "UT5E157","UT5E092","UT5E093","UT5E094","UT5E185","UT5E186","UT5E187",
                      "UT5E077","UT5E078","UT5E079","UT5E170","UT5E171","UT5E172","UVIE036",
                      "UVIE037","UVIE038","UVIE077","UVIE078","UVIE079","UVIE039","UVIE040",
                      "UVIE080","UVIE081","UVIE041","UVIE042","UVIE082","UVIE083","UV9E015",
                      "UV9E016","UV9E029","UV9E030","UXLE016","UXLE019","UXLE035","UXLE038",
                      "UXYE013","UXYE016","UXYE029","UXYE032", "UXXE013","UXXE016","UXXE029",
                      "UXXE032")

# 5-year ACS 2009-2013
zipcode_2013_1 <- read.csv("/Volumes/PACS$/Cordes, Jack/Zip Code Data/nhgis0042_csv/nhgis0042_ds201_20135_zcta.csv")
zipcode_2013_1 <- read.csv("/Volumes/PACS$/Cordes, Jack/Zip Code Data/nhgis0042_csv/nhgis0042_ds201_20135_zcta.csv")[2:nrow(zipcode_2013_1), acs_codes_2013_1]
zipcode_2013_1[, 4:length(zipcode_2013_1)] <- sapply(zipcode_2013_1[, 4:length(zipcode_2013_1)], as.numeric)

zipcode_2013_2 <- read.csv("/Volumes/PACS$/Cordes, Jack/Zip Code Data/nhgis0042_csv/nhgis0042_ds202_20135_zcta.csv")
zipcode_2013_2 <- read.csv("/Volumes/PACS$/Cordes, Jack/Zip Code Data/nhgis0042_csv/nhgis0042_ds202_20135_zcta.csv")[2:nrow(zipcode_2013_2), acs_codes_2013_2]
zipcode_2013_2[, 2:length(zipcode_2013_2)] <- sapply(zipcode_2013_2[, 2:length(zipcode_2013_2)], as.numeric)

zipcode_2013 <- left_join(zipcode_2013_1, zipcode_2013_2, by = "ZCTA5A")

zipcode_2013$PROP_WHITE <- zipcode_2013$UEQE002/zipcode_2013$UEEE001
zipcode_2013$PROP_BLACK <- zipcode_2013$UEQE003/zipcode_2013$UEEE001
zipcode_2013$PROP_AIAN <- zipcode_2013$UEQE004/zipcode_2013$UEEE001
zipcode_2013$PROP_ASIAN <- zipcode_2013$UEQE005/zipcode_2013$UEEE001
zipcode_2013$PROP_NHP <- zipcode_2013$UEQE006/zipcode_2013$UEEE001
zipcode_2013$PROP_OTHER_RACE <- zipcode_2013$UEQE007/zipcode_2013$UEEE001
zipcode_2013$PROP_2MORE_RACE <- zipcode_2013$UEQE008/zipcode_2013$UEEE001

zipcode_2013$PROP_NONHISP_WHITE <- zipcode_2013$UEYE003/zipcode_2013$UEEE001
zipcode_2013$PROP_NONHISP_BLACK <- zipcode_2013$UEYE004/zipcode_2013$UEEE001
zipcode_2013$PROP_NONHISP_AIAN <- zipcode_2013$UEYE005/zipcode_2013$UEEE001
zipcode_2013$PROP_NONHISP_ASIAN <- zipcode_2013$UEYE006/zipcode_2013$UEEE001
zipcode_2013$PROP_HISP <- zipcode_2013$UEYE012/zipcode_2013$UEEE001
zipcode_2013$PROP_HISP_WHITE <- zipcode_2013$UEYE013/zipcode_2013$UEEE001
zipcode_2013$PROP_HISP_BLACK <- zipcode_2013$UEYE014/zipcode_2013$UEEE001
zipcode_2013$PROP_HISP_AIAN <- zipcode_2013$UEYE015/zipcode_2013$UEEE001
zipcode_2013$PROP_HISP_ASIAN <- zipcode_2013$UEYE016/zipcode_2013$UEEE001

zipcode_2013$PROP_ENGLISH_ONLY <- zipcode_2013$UG6E047/zipcode_2013$UG6E046
zipcode_2013$PROP_ENGLISH_WELL <- (zipcode_2013$UG6E049 + zipcode_2013$UG6E050 + 
                                     zipcode_2013$UG6E054 + zipcode_2013$UG6E055 + 
                                     zipcode_2013$UG6E059 + zipcode_2013$UG6E060 + 
                                     zipcode_2013$UG6E064 + zipcode_2013$UG6E065)/zipcode_2013$UG6E046
zipcode_2013$PROP_NO_ENGLISH <- (zipcode_2013$UG6E052 + zipcode_2013$UG6E057 +
                                   zipcode_2013$UG6E062 + zipcode_2013$UG6E067)/zipcode_2013$UG6E046

zipcode_2013$PROP_PUBLIC_ASST <- zipcode_2013$UH6E002/zipcode_2013$UH6E001

zipcode_2013$PROP_VET <- (zipcode_2013$UJXM017 + zipcode_2013$UJXM020 + zipcode_2013$UJXM035 + 
                            zipcode_2013$UJXM038)/zipcode_2013$UG6E046

zipcode_2013$PROP_NATIVE <- (zipcode_2013$UPIE004 + zipcode_2013$UPIE009 + zipcode_2013$UPIE015 +
                               zipcode_2013$UPIE020)/zipcode_2013$UEEE001

zipcode_2013$PROP_NEVER_MARRIED <- (zipcode_2013$UT5E015 + zipcode_2013$UT5E016 + zipcode_2013$UT5E017 +
                                      zipcode_2013$UT5E108 + zipcode_2013$UT5E109 + zipcode_2013$UT5E110)/zipcode_2013$UG6E046

zipcode_2013$PROP_MARRIED <- (zipcode_2013$UT5E031 + zipcode_2013$UT5E032 + zipcode_2013$UT5E033 + 
                                zipcode_2013$UT5E047 + zipcode_2013$UT5E048 + zipcode_2013$UT5E049 + 
                                zipcode_2013$UT5E062 + zipcode_2013$UT5E063 + zipcode_2013$UT5E064 + 
                                zipcode_2013$UT5E124 + zipcode_2013$UT5E125 + zipcode_2013$UT5E126 + 
                                zipcode_2013$UT5E140 + zipcode_2013$UT5E141 + zipcode_2013$UT5E142 + 
                                zipcode_2013$UT5E155 + zipcode_2013$UT5E156 + zipcode_2013$UT5E157)/zipcode_2013$UG6E046

zipcode_2013$PROP_DIVORCED <- (zipcode_2013$UT5E092 + zipcode_2013$UT5E093 + zipcode_2013$UT5E094 + 
                                 zipcode_2013$UT5E185 + zipcode_2013$UT5E186 + zipcode_2013$UT5E187)/zipcode_2013$UG6E046

zipcode_2013$PROP_WIDOWED <- (zipcode_2013$UT5E077 + zipcode_2013$UT5E078 + zipcode_2013$UT5E079 + 
                                zipcode_2013$UT5E170 + zipcode_2013$UT5E171 + zipcode_2013$UT5E172)/zipcode_2013$UG6E046

zipcode_2013$PROP_HS <- (zipcode_2013$UVIE036 + zipcode_2013$UVIE037 + zipcode_2013$UVIE038 + 
                           zipcode_2013$UVIE077 + zipcode_2013$UVIE078 + zipcode_2013$UVIE079)/zipcode_2013$UG6E046
zipcode_2013$PROP_ASSOC <- (zipcode_2013$UVIE039 + zipcode_2013$UVIE040 + zipcode_2013$UVIE080 +
                              zipcode_2013$UVIE081)/zipcode_2013$UG6E046
zipcode_2013$PROP_BACHELORS <- (zipcode_2013$UVIE041 + zipcode_2013$UVIE042 + zipcode_2013$UVIE082 +
                                  zipcode_2013$UVIE083)/zipcode_2013$UG6E046

zipcode_2013$PROP_POVERTY <- (zipcode_2013$UV9E015 + zipcode_2013$UV9E016 + zipcode_2013$UV9E029 + 
                                zipcode_2013$UV9E030)/zipcode_2013$UG6E046

zipcode_2013$PROP_DISABILITY <- (zipcode_2013$UXLE016 + zipcode_2013$UXLE019 + zipcode_2013$UXLE035 + 
                                   zipcode_2013$UXLE038)/zipcode_2013$UM2E051
zipcode_2013$PROP_COGNITIVE <- (zipcode_2013$UXYE013 + zipcode_2013$UXYE016 + zipcode_2013$UXYE029 + 
                                  zipcode_2013$UXYE032)/zipcode_2013$UM2E051
zipcode_2013$PROP_AMBULATORY <- (zipcode_2013$UXXE013 + zipcode_2013$UXXE016 + zipcode_2013$UXXE029 + 
                                   zipcode_2013$UXXE032)/zipcode_2013$UM2E051

zipcode_2013$PROP_MEDICARE <- (zipcode_2013$UM2E055 + zipcode_2013$UM2E060 + zipcode_2013$UM2E061 + 
                                 zipcode_2013$UM2E062)/zipcode_2013$UM2E051

zipcode_2013$MED_INC <- zipcode_2013$UHDE001
zipcode_2013$MED_INC_65 <- zipcode_2013$UHYE005
zipcode_2013$GINI <- zipcode_2013$UILE001

# List of pertinent ACS variable codes 2014
acs_codes_2014_1 <- c("GISJOIN","YEAR","ZCTA5A","ABAQE001","ABDIE046","ABA2E002","ABA2E003",
                      "ABA2E004","ABA2E005","ABA2E006","ABA2E007","ABA2E008","ABBAE003",
                      "ABBAE004","ABBAE005","ABBAE006","ABBAE012","ABBAE013","ABBAE014",
                      "ABBAE015","ABBAE016","ABDIE047","ABDIE049","ABDIE050", "ABDIE052",
                      "ABDIE054","ABDIE055","ABDIE057","ABDIE059","ABDIE060","ABDIE062",
                      "ABDIE064","ABDIE065","ABDIE067","ABDPE001","ABEAE005","ABEIE001",
                      "ABEIE002","ABFIE001","ABF5E017","ABF5E020","ABF5E035","ABF5E038",
                      "ABI8E055","ABI8E060","ABI8E061","ABI8E062","ABI8E051")

acs_codes_2014_2 <- c("ZCTA5A","ABLLE004","ABLLE009","ABLLE015","ABLLE020","ABQAE015","ABQAE016",
                      "ABQAE017","ABQAE108","ABQAE109","ABQAE110","ABQAE031","ABQAE032",
                      "ABQAE033","ABQAE047","ABQAE048","ABQAE049","ABQAE062","ABQAE063",
                      "ABQAE064","ABQAE124","ABQAE125","ABQAE126","ABQAE140","ABQAE141",
                      "ABQAE142","ABQAE155","ABQAE156","ABQAE157","ABQAE092","ABQAE093",
                      "ABQAE094","ABQAE185","ABQAE186","ABQAE187","ABQAE077","ABQAE078",
                      "ABQAE079","ABQAE170","ABQAE171","ABQAE172","ABRNE036","ABRNE037",
                      "ABRNE038","ABRNE077","ABRNE078","ABRNE079","ABRNE039","ABRNE040",
                      "ABRNE080","ABRNE081","ABRNE041","ABRNE042","ABRNE082","ABRNE083",
                      "ABSEE015","ABSEE016","ABSEE029","ABSEE030","ABTQE016","ABTQE019",
                      "ABTQE035","ABTQE038","ABT2E013","ABT2E016","ABT2E029","ABT2E032",
                      "ABT3E013","ABT3E016","ABT3E029","ABT3E032","ABU0E001")

# 5-year ACS 2010-2014
zipcode_2014_1 <- read.csv("/Volumes/PACS$/Cordes, Jack/Zip Code Data/nhgis0043_csv/nhgis0043_ds206_20145_zcta.csv")
zipcode_2014_1 <- read.csv("/Volumes/PACS$/Cordes, Jack/Zip Code Data/nhgis0043_csv/nhgis0043_ds206_20145_zcta.csv")[2:nrow(zipcode_2014_1), acs_codes_2014_1]
zipcode_2014_1[, 4:length(zipcode_2014_1)] <- sapply(zipcode_2014_1[, 4:length(zipcode_2014_1)], as.numeric)

zipcode_2014_2 <- read.csv("/Volumes/PACS$/Cordes, Jack/Zip Code Data/nhgis0043_csv/nhgis0043_ds207_20145_zcta.csv")
zipcode_2014_2 <- read.csv("/Volumes/PACS$/Cordes, Jack/Zip Code Data/nhgis0043_csv/nhgis0043_ds207_20145_zcta.csv")[2:nrow(zipcode_2014_2), acs_codes_2014_2]
zipcode_2014_2[, 2:length(zipcode_2014_2)] <- sapply(zipcode_2014_2[, 2:length(zipcode_2014_2)], as.numeric)

zipcode_2014 <- left_join(zipcode_2014_1, zipcode_2014_2, by = "ZCTA5A")

zipcode_2014$PROP_WHITE <- zipcode_2014$ABA2E002/zipcode_2014$ABAQE001
zipcode_2014$PROP_BLACK <- zipcode_2014$ABA2E003/zipcode_2014$ABAQE001
zipcode_2014$PROP_AIAN <- zipcode_2014$ABA2E004/zipcode_2014$ABAQE001
zipcode_2014$PROP_ASIAN <- zipcode_2014$ABA2E005/zipcode_2014$ABAQE001
zipcode_2014$PROP_NHP <- zipcode_2014$ABA2E006/zipcode_2014$ABAQE001
zipcode_2014$PROP_OTHER_RACE <- zipcode_2014$ABA2E007/zipcode_2014$ABAQE001
zipcode_2014$PROP_2MORE_RACE <- zipcode_2014$ABA2E008/zipcode_2014$ABAQE001

zipcode_2014$PROP_NONHISP_WHITE <- zipcode_2014$ABBAE003/zipcode_2014$ABAQE001
zipcode_2014$PROP_NONHISP_BLACK <- zipcode_2014$ABBAE004/zipcode_2014$ABAQE001
zipcode_2014$PROP_NONHISP_AIAN <- zipcode_2014$ABBAE005/zipcode_2014$ABAQE001
zipcode_2014$PROP_NONHISP_ASIAN <- zipcode_2014$ABBAE006/zipcode_2014$ABAQE001
zipcode_2014$PROP_HISP <- zipcode_2014$ABBAE012/zipcode_2014$ABAQE001
zipcode_2014$PROP_HISP_WHITE <- zipcode_2014$ABBAE013/zipcode_2014$ABAQE001
zipcode_2014$PROP_HISP_BLACK <- zipcode_2014$ABBAE014/zipcode_2014$ABAQE001
zipcode_2014$PROP_HISP_AIAN <- zipcode_2014$ABBAE015/zipcode_2014$ABAQE001
zipcode_2014$PROP_HISP_ASIAN <- zipcode_2014$ABBAE016/zipcode_2014$ABAQE001

zipcode_2014$PROP_ENGLISH_ONLY <- zipcode_2014$ABDIE047/zipcode_2014$ABDIE046
zipcode_2014$PROP_ENGLISH_WELL <- (zipcode_2014$ABDIE049 + zipcode_2014$ABDIE050 + 
                                     zipcode_2014$ABDIE054 + zipcode_2014$ABDIE055 + 
                                     zipcode_2014$ABDIE059 + zipcode_2014$ABDIE060 + 
                                     zipcode_2014$ABDIE064 + zipcode_2014$ABDIE065)/zipcode_2014$ABDIE046
zipcode_2014$PROP_NO_ENGLISH <- (zipcode_2014$ABDIE052 + zipcode_2014$ABDIE057 +
                                   zipcode_2014$ABDIE062 + zipcode_2014$ABDIE067)/zipcode_2014$ABDIE046

zipcode_2014$PROP_PUBLIC_ASST <- zipcode_2014$ABEIE002/zipcode_2014$ABEIE001

zipcode_2014$PROP_VET <- (zipcode_2014$ABF5E017 + zipcode_2014$ABF5E020 + zipcode_2014$ABF5E035 + 
                            zipcode_2014$ABF5E038)/zipcode_2014$ABDIE046

zipcode_2014$PROP_NATIVE <- (zipcode_2014$ABLLE004 + zipcode_2014$ABLLE009 + zipcode_2014$ABLLE015 +
                               zipcode_2014$ABLLE020)/zipcode_2014$ABAQE001

zipcode_2014$PROP_NEVER_MARRIED <- (zipcode_2014$ABQAE015 + zipcode_2014$ABQAE016 + zipcode_2014$ABQAE017 +
                                      zipcode_2014$ABQAE108 + zipcode_2014$ABQAE109 + zipcode_2014$ABQAE110)/zipcode_2014$ABDIE046

zipcode_2014$PROP_MARRIED <- (zipcode_2014$ABQAE031 + zipcode_2014$ABQAE032 + zipcode_2014$ABQAE033 + 
                                zipcode_2014$ABQAE047 + zipcode_2014$ABQAE048 + zipcode_2014$ABQAE049 + 
                                zipcode_2014$ABQAE062 + zipcode_2014$ABQAE063 + zipcode_2014$ABQAE064 + 
                                zipcode_2014$ABQAE124 + zipcode_2014$ABQAE125 + zipcode_2014$ABQAE126 + 
                                zipcode_2014$ABQAE140 + zipcode_2014$ABQAE141 + zipcode_2014$ABQAE142 + 
                                zipcode_2014$ABQAE155 + zipcode_2014$ABQAE156 + zipcode_2014$ABQAE157)/zipcode_2014$ABDIE046

zipcode_2014$PROP_DIVORCED <- (zipcode_2014$ABQAE092 + zipcode_2014$ABQAE093 + zipcode_2014$ABQAE094 + 
                                 zipcode_2014$ABQAE185 + zipcode_2014$ABQAE186 + zipcode_2014$ABQAE187)/zipcode_2014$ABDIE046

zipcode_2014$PROP_WIDOWED <- (zipcode_2014$ABQAE077 + zipcode_2014$ABQAE078 + zipcode_2014$ABQAE079 + 
                                zipcode_2014$ABQAE170 + zipcode_2014$ABQAE171 + zipcode_2014$ABQAE172)/zipcode_2014$ABDIE046

zipcode_2014$PROP_HS <- (zipcode_2014$ABRNE036 + zipcode_2014$ABRNE037 + zipcode_2014$ABRNE038 + 
                           zipcode_2014$ABRNE077 + zipcode_2014$ABRNE078 + zipcode_2014$ABRNE079)/zipcode_2014$ABDIE046
zipcode_2014$PROP_ASSOC <- (zipcode_2014$ABRNE039 + zipcode_2014$ABRNE040 + zipcode_2014$ABRNE080 +
                              zipcode_2014$ABRNE081)/zipcode_2014$ABDIE046
zipcode_2014$PROP_BACHELORS <- (zipcode_2014$ABRNE041 + zipcode_2014$ABRNE042 + zipcode_2014$ABRNE082 +
                                  zipcode_2014$ABRNE083)/zipcode_2014$ABDIE046

zipcode_2014$PROP_POVERTY <- (zipcode_2014$ABSEE015 + zipcode_2014$ABSEE016 + zipcode_2014$ABSEE029 + 
                                zipcode_2014$ABSEE030)/zipcode_2014$ABDIE046

zipcode_2014$PROP_DISABILITY <- (zipcode_2014$ABTQE016 + zipcode_2014$ABTQE019 + zipcode_2014$ABTQE035 + 
                                   zipcode_2014$ABTQE038)/zipcode_2014$ABI8E051
zipcode_2014$PROP_COGNITIVE <- (zipcode_2014$ABT2E013 + zipcode_2014$ABT2E016 + zipcode_2014$ABT2E029 + 
                                  zipcode_2014$ABT2E032)/zipcode_2014$ABI8E051
zipcode_2014$PROP_AMBULATORY <- (zipcode_2014$ABT3E013 + zipcode_2014$ABT3E016 + zipcode_2014$ABT3E029 + 
                                   zipcode_2014$ABT3E032)/zipcode_2014$ABI8E051

zipcode_2014$PROP_MEDICARE <- (zipcode_2014$ABI8E055 + zipcode_2014$ABI8E060 + zipcode_2014$ABI8E061 + 
                                 zipcode_2014$ABI8E062)/zipcode_2014$ABI8E051

zipcode_2014$MED_INC <- zipcode_2014$ABDPE001
zipcode_2014$MED_INC_65 <- zipcode_2014$ABEAE005
zipcode_2014$GINI <- zipcode_2014$ABU0E001

# List of pertinent ACS variable codes 2015
acs_codes_2015_1 <- c("GISJOIN","YEAR","ZCTA5A","ADKLE001","ADNDE046","ADKXE002","ADKXE003",
                      "ADKXE004","ADKXE005","ADKXE006","ADKXE007","ADKXE008","ADK5E003",
                      "ADK5E004","ADK5E005","ADK5E006","ADK5E012","ADK5E013","ADK5E014",
                      "ADK5E015","ADK5E016","ADNDE047","ADNDE049","ADNDE050","ADNDE052",
                      "ADNDE054","ADNDE055","ADNDE057","ADNDE059","ADNDE060","ADNDE062",
                      "ADNDE064","ADNDE065","ADNDE067","ADNKE001","ADNWE005","ADN4E001",
                      "ADN4E002","ADOLE001","ADO8E017","ADO8E020","ADO8E035","ADO8E038",
                      "ADSBE055","ADSBE060","ADSBE061","ADSBE062","ADSBE051")

acs_codes_2015_2 <- c("ZCTA5A","ADURE004","ADURE009","ADURE015","ADURE020","ADZEE015","ADZEE016",
                      "ADZEE017","ADZEE108","ADZEE109","ADZEE110","ADZEE031","ADZEE032",
                      "ADZEE033","ADZEE047","ADZEE048","ADZEE049","ADZEE062","ADZEE063",
                      "ADZEE064","ADZEE124","ADZEE125","ADZEE126","ADZEE140","ADZEE141",
                      "ADZEE142","ADZEE155","ADZEE156","ADZEE157","ADZEE092","ADZEE093",
                      "ADZEE094","ADZEE185","ADZEE186","ADZEE187","ADZEE077","ADZEE078",
                      "ADZEE079","ADZEE170","ADZEE171","ADZEE172","AD0PE036","AD0PE037",
                      "AD0PE038","AD0PE077","AD0PE078","AD0PE079","AD0PE039","AD0PE040",
                      "AD0PE080","AD0PE081","AD0PE041","AD0PE042","AD0PE082","AD0PE083",
                      "AD1GE015","AD1GE016","AD1GE029","AD1GE030","AD2SE016","AD2SE019",
                      "AD2SE035","AD2SE038","AD24E013","AD24E016","AD24E029","AD24E032",
                      "AD25E013","AD25E016","AD25E029","AD25E032","AD4BE001")

# 5-year ACS 2011-2015
zipcode_2015_1 <- read.csv("/Volumes/PACS$/Cordes, Jack/Zip Code Data/nhgis0044_csv/nhgis0044_ds215_20155_zcta.csv")
zipcode_2015_1 <- read.csv("/Volumes/PACS$/Cordes, Jack/Zip Code Data/nhgis0044_csv/nhgis0044_ds215_20155_zcta.csv")[2:nrow(zipcode_2015_1), acs_codes_2015_1]
zipcode_2015_1[, 4:length(zipcode_2015_1)] <- sapply(zipcode_2015_1[, 4:length(zipcode_2015_1)], as.numeric)

zipcode_2015_2 <- read.csv("/Volumes/PACS$/Cordes, Jack/Zip Code Data/nhgis0044_csv/nhgis0044_ds216_20155_zcta.csv")
zipcode_2015_2 <- read.csv("/Volumes/PACS$/Cordes, Jack/Zip Code Data/nhgis0044_csv/nhgis0044_ds216_20155_zcta.csv")[2:nrow(zipcode_2015_2), acs_codes_2015_2]
zipcode_2015_2[, 2:length(zipcode_2015_2)] <- sapply(zipcode_2015_2[, 2:length(zipcode_2015_2)], as.numeric)

zipcode_2015 <- left_join(zipcode_2015_1, zipcode_2015_2, by = "ZCTA5A")

zipcode_2015$PROP_WHITE <- zipcode_2015$ADKXE002/zipcode_2015$ADKLE001
zipcode_2015$PROP_BLACK <- zipcode_2015$ADKXE003/zipcode_2015$ADKLE001
zipcode_2015$PROP_AIAN <- zipcode_2015$ADKXE004/zipcode_2015$ADKLE001
zipcode_2015$PROP_ASIAN <- zipcode_2015$ADKXE005/zipcode_2015$ADKLE001
zipcode_2015$PROP_NHP <- zipcode_2015$ADKXE006/zipcode_2015$ADKLE001
zipcode_2015$PROP_OTHER_RACE <- zipcode_2015$ADKXE007/zipcode_2015$ADKLE001
zipcode_2015$PROP_2MORE_RACE <- zipcode_2015$ADKXE008/zipcode_2015$ADKLE001

zipcode_2015$PROP_NONHISP_WHITE <- zipcode_2015$ADK5E003/zipcode_2015$ADKLE001
zipcode_2015$PROP_NONHISP_BLACK <- zipcode_2015$ADK5E004/zipcode_2015$ADKLE001
zipcode_2015$PROP_NONHISP_AIAN <- zipcode_2015$ADK5E005/zipcode_2015$ADKLE001
zipcode_2015$PROP_NONHISP_ASIAN <- zipcode_2015$ADK5E006/zipcode_2015$ADKLE001
zipcode_2015$PROP_HISP <- zipcode_2015$ADK5E012/zipcode_2015$ADKLE001
zipcode_2015$PROP_HISP_WHITE <- zipcode_2015$ADK5E013/zipcode_2015$ADKLE001
zipcode_2015$PROP_HISP_BLACK <- zipcode_2015$ADK5E014/zipcode_2015$ADKLE001
zipcode_2015$PROP_HISP_AIAN <- zipcode_2015$ADK5E015/zipcode_2015$ADKLE001
zipcode_2015$PROP_HISP_ASIAN <- zipcode_2015$ADK5E016/zipcode_2015$ADKLE001

zipcode_2015$PROP_ENGLISH_ONLY <- zipcode_2015$ADNDE047/zipcode_2015$ADNDE046
zipcode_2015$PROP_ENGLISH_WELL <- (zipcode_2015$ADNDE049 + zipcode_2015$ADNDE050 + 
                                     zipcode_2015$ADNDE054 + zipcode_2015$ADNDE055 + 
                                     zipcode_2015$ADNDE059 + zipcode_2015$ADNDE060 + 
                                     zipcode_2015$ADNDE064 + zipcode_2015$ADNDE065)/zipcode_2015$ADNDE046
zipcode_2015$PROP_NO_ENGLISH <- (zipcode_2015$ADNDE052 + zipcode_2015$ADNDE057 + 
                                   zipcode_2015$ADNDE062 + zipcode_2015$ADNDE067)/zipcode_2015$ADNDE046

zipcode_2015$PROP_PUBLIC_ASST <- zipcode_2015$ADN4E002/zipcode_2015$ADN4E001

zipcode_2015$PROP_VET <- (zipcode_2015$ADO8E017 + zipcode_2015$ADO8E020 + zipcode_2015$ADO8E035 + 
                            zipcode_2015$ADO8E038)/zipcode_2015$ADNDE046

zipcode_2015$PROP_NATIVE <- (zipcode_2015$ADURE004 + zipcode_2015$ADURE009 + zipcode_2015$ADURE015 +
                               zipcode_2015$ADURE020)/zipcode_2015$ADKLE001

zipcode_2015$PROP_NEVER_MARRIED <- (zipcode_2015$ADZEE015 + zipcode_2015$ADZEE016 + zipcode_2015$ADZEE017 +
                                      zipcode_2015$ADZEE108 + zipcode_2015$ADZEE109 + zipcode_2015$ADZEE110)/zipcode_2015$ADNDE046

zipcode_2015$PROP_MARRIED <- (zipcode_2015$ADZEE031 + zipcode_2015$ADZEE032 + zipcode_2015$ADZEE033 + 
                                zipcode_2015$ADZEE047 + zipcode_2015$ADZEE048 + zipcode_2015$ADZEE049 + 
                                zipcode_2015$ADZEE062 + zipcode_2015$ADZEE063 + zipcode_2015$ADZEE064 + 
                                zipcode_2015$ADZEE124 + zipcode_2015$ADZEE125 + zipcode_2015$ADZEE126 + 
                                zipcode_2015$ADZEE140 + zipcode_2015$ADZEE141 + zipcode_2015$ADZEE142 + 
                                zipcode_2015$ADZEE155 + zipcode_2015$ADZEE156 + zipcode_2015$ADZEE157)/zipcode_2015$ADNDE046

zipcode_2015$PROP_DIVORCED <- (zipcode_2015$ADZEE092 + zipcode_2015$ADZEE093 + zipcode_2015$ADZEE094 + 
                                 zipcode_2015$ADZEE185 + zipcode_2015$ADZEE186 + zipcode_2015$ADZEE187)/zipcode_2015$ADNDE046

zipcode_2015$PROP_WIDOWED <- (zipcode_2015$ADZEE077 + zipcode_2015$ADZEE078 + zipcode_2015$ADZEE079 + 
                                zipcode_2015$ADZEE170 + zipcode_2015$ADZEE171 + zipcode_2015$ADZEE172)/zipcode_2015$ADNDE046

zipcode_2015$PROP_HS <- (zipcode_2015$AD0PE036 + zipcode_2015$AD0PE037 + zipcode_2015$AD0PE038 + 
                           zipcode_2015$AD0PE077 + zipcode_2015$AD0PE078 + zipcode_2015$AD0PE079)/zipcode_2015$ADNDE046
zipcode_2015$PROP_ASSOC <- (zipcode_2015$AD0PE039 + zipcode_2015$AD0PE040 + zipcode_2015$AD0PE080 +
                              zipcode_2015$AD0PE081)/zipcode_2015$ADNDE046
zipcode_2015$PROP_BACHELORS <- (zipcode_2015$AD0PE041 + zipcode_2015$AD0PE042 + zipcode_2015$AD0PE082 +
                                  zipcode_2015$AD0PE083)/zipcode_2015$ADNDE046

zipcode_2015$PROP_POVERTY <- (zipcode_2015$AD1GE015 + zipcode_2015$AD1GE016 + zipcode_2015$AD1GE029 + 
                                zipcode_2015$AD1GE030)/zipcode_2015$ADNDE046

zipcode_2015$PROP_DISABILITY <- (zipcode_2015$AD2SE016 + zipcode_2015$AD2SE019 + zipcode_2015$AD2SE035 + 
                                   zipcode_2015$AD2SE038)/zipcode_2015$ADSBE051
zipcode_2015$PROP_COGNITIVE <- (zipcode_2015$AD24E013 + zipcode_2015$AD24E016 + zipcode_2015$AD24E029 + 
                                  zipcode_2015$AD24E032)/zipcode_2015$ADSBE051
zipcode_2015$PROP_AMBULATORY <- (zipcode_2015$AD25E013 + zipcode_2015$AD25E016 + zipcode_2015$AD25E029 + 
                                   zipcode_2015$AD25E032)/zipcode_2015$ADSBE051

zipcode_2015$PROP_MEDICARE <- (zipcode_2015$ADSBE055 + zipcode_2015$ADSBE060 + zipcode_2015$ADSBE061 + 
                                 zipcode_2015$ADSBE062)/zipcode_2015$ADSBE051

zipcode_2015$MED_INC <- zipcode_2015$ADNKE001
zipcode_2015$MED_INC_65 <- zipcode_2015$ADNWE005
zipcode_2015$GINI <- zipcode_2015$AD4BE001

# List of pertinent ACS variable codes 2016
acs_codes_2016_1 <- c("GISJOIN","YEAR","ZCTA5A","AF2AE001","AF42E046","AF2ME002","AF2ME003",
                      "AF2ME004","AF2ME005","AF2ME006","AF2ME007","AF2ME008","AF2UE003",
                      "AF2UE004","AF2UE005","AF2UE006","AF2UE012","AF2UE013","AF2UE014",
                      "AF2UE015","AF2UE016","AF42E047","AF42E049","AF42E050","AF42E052",
                      "AF42E054","AF42E055","AF42E057","AF42E059","AF42E060","AF42E062",
                      "AF42E064","AF42E065","AF42E067","AF49E001","AF5LE005","AF5TE001",
                      "AF5TE002","AF6AE001","AF6XE017","AF6XE020","AF6XE035","AF6XE038",
                      "AF90E055","AF90E060","AF90E061","AF90E062","AF90E051")

acs_codes_2016_2 <- c("ZCTA5A","AGCGE004","AGCGE009","AGCGE015","AGCGE020","AGG3E015",
                      "AGG3E016","AGG3E017","AGG3E108","AGG3E109","AGG3E110","AGG3E031",
                      "AGG3E032","AGG3E033","AGG3E047","AGG3E048","AGG3E049","AGG3E062",
                      "AGG3E063","AGG3E064","AGG3E124","AGG3E125","AGG3E126","AGG3E140",
                      "AGG3E141","AGG3E142","AGG3E155","AGG3E156","AGG3E157","AGG3E092",
                      "AGG3E093","AGG3E094","AGG3E185","AGG3E186","AGG3E187","AGG3E077",
                      "AGG3E078","AGG3E079","AGG3E170","AGG3E171","AGG3E172","AGIEE036",
                      "AGIEE037","AGIEE038","AGIEE077","AGIEE078","AGIEE079","AGIEE039",
                      "AGIEE040","AGIEE080","AGIEE081","AGIEE041","AGIEE042","AGIEE082",
                      "AGIEE083","AGI6E015","AGI6E016","AGI6E029","AGI6E030","AGKIE016",
                      "AGKIE019","AGKIE035","AGKIE038","AGKUE013","AGKUE016","AGKUE029",
                      "AGKUE032","AGKVE013","AGKVE016","AGKVE029","AGKVE032","AGL1E001")

# 5-year ACS 2012-2016
zipcode_2016_1 <- read.csv("/Volumes/PACS$/Cordes, Jack/Zip Code Data/nhgis0045_csv/nhgis0045_ds225_20165_zcta.csv")
zipcode_2016_1 <- read.csv("/Volumes/PACS$/Cordes, Jack/Zip Code Data/nhgis0045_csv/nhgis0045_ds225_20165_zcta.csv")[2:nrow(zipcode_2016_1), acs_codes_2016_1]
zipcode_2016_1[, 4:length(zipcode_2016_1)] <- sapply(zipcode_2016_1[, 4:length(zipcode_2016_1)], as.numeric)

zipcode_2016_2 <- read.csv("/Volumes/PACS$/Cordes, Jack/Zip Code Data/nhgis0045_csv/nhgis0045_ds226_20165_zcta.csv")
zipcode_2016_2 <- read.csv("/Volumes/PACS$/Cordes, Jack/Zip Code Data/nhgis0045_csv/nhgis0045_ds226_20165_zcta.csv")[2:nrow(zipcode_2016_2), acs_codes_2016_2]
zipcode_2016_2[, 2:length(zipcode_2016_2)] <- sapply(zipcode_2016_2[, 2:length(zipcode_2016_2)], as.numeric)

zipcode_2016 <- left_join(zipcode_2016_1, zipcode_2016_2, by = "ZCTA5A")

zipcode_2016$PROP_WHITE <- zipcode_2016$AF2ME002/zipcode_2016$AF2AE001
zipcode_2016$PROP_BLACK <- zipcode_2016$AF2ME003/zipcode_2016$AF2AE001
zipcode_2016$PROP_AIAN <- zipcode_2016$AF2ME004/zipcode_2016$AF2AE001
zipcode_2016$PROP_ASIAN <- zipcode_2016$AF2ME005/zipcode_2016$AF2AE001
zipcode_2016$PROP_NHP <- zipcode_2016$AF2ME006/zipcode_2016$AF2AE001
zipcode_2016$PROP_OTHER_RACE <- zipcode_2016$AF2ME007/zipcode_2016$AF2AE001
zipcode_2016$PROP_2MORE_RACE <- zipcode_2016$AF2ME008/zipcode_2016$AF2AE001

zipcode_2016$PROP_NONHISP_WHITE <- zipcode_2016$AF2UE003/zipcode_2016$AF2AE001
zipcode_2016$PROP_NONHISP_BLACK <- zipcode_2016$AF2UE004/zipcode_2016$AF2AE001
zipcode_2016$PROP_NONHISP_AIAN <- zipcode_2016$AF2UE005/zipcode_2016$AF2AE001
zipcode_2016$PROP_NONHISP_ASIAN <- zipcode_2016$AF2UE006/zipcode_2016$AF2AE001
zipcode_2016$PROP_HISP <- zipcode_2016$AF2UE012/zipcode_2016$AF2AE001
zipcode_2016$PROP_HISP_WHITE <- zipcode_2016$AF2UE013/zipcode_2016$AF2AE001
zipcode_2016$PROP_HISP_BLACK <- zipcode_2016$AF2UE014/zipcode_2016$AF2AE001
zipcode_2016$PROP_HISP_AIAN <- zipcode_2016$AF2UE015/zipcode_2016$AF2AE001
zipcode_2016$PROP_HISP_ASIAN <- zipcode_2016$AF2UE016/zipcode_2016$AF2AE001

zipcode_2016$PROP_ENGLISH_ONLY <- zipcode_2016$AF42E047/zipcode_2016$AF42E046
zipcode_2016$PROP_ENGLISH_WELL <- (zipcode_2016$AF42E049 + zipcode_2016$AF42E050 + 
                                     zipcode_2016$AF42E054 + zipcode_2016$AF42E055 + 
                                     zipcode_2016$AF42E059 + zipcode_2016$AF42E060 + 
                                     zipcode_2016$AF42E064 + zipcode_2016$AF42E065)/zipcode_2016$AF42E046
zipcode_2016$PROP_NO_ENGLISH <- (zipcode_2016$AF42E052 + zipcode_2016$AF42E057 +
                                   zipcode_2016$AF42E062 + zipcode_2016$AF42E067)/zipcode_2016$AF42E046

zipcode_2016$PROP_PUBLIC_ASST <- zipcode_2016$AF5TE002/zipcode_2016$AF5TE001

zipcode_2016$PROP_VET <- (zipcode_2016$AF6XE017 + zipcode_2016$AF6XE020 + zipcode_2016$AF6XE035 + 
                            zipcode_2016$AF6XE038)/zipcode_2016$AF42E046

zipcode_2016$PROP_NATIVE <- (zipcode_2016$AGCGE004 + zipcode_2016$AGCGE009 + zipcode_2016$AGCGE015 +
                               zipcode_2016$AGCGE020)/zipcode_2016$AF2AE001

zipcode_2016$PROP_NEVER_MARRIED <- (zipcode_2016$AGG3E015 + zipcode_2016$AGG3E016 + zipcode_2016$AGG3E017 +
                                      zipcode_2016$AGG3E108 + zipcode_2016$AGG3E109 + zipcode_2016$AGG3E110)/zipcode_2016$AF42E046

zipcode_2016$PROP_MARRIED <- (zipcode_2016$AGG3E031 + zipcode_2016$AGG3E032 + zipcode_2016$AGG3E033 + 
                                zipcode_2016$AGG3E047 + zipcode_2016$AGG3E048 + zipcode_2016$AGG3E049 + 
                                zipcode_2016$AGG3E062 + zipcode_2016$AGG3E063 + zipcode_2016$AGG3E064 + 
                                zipcode_2016$AGG3E124 + zipcode_2016$AGG3E125 + zipcode_2016$AGG3E126 + 
                                zipcode_2016$AGG3E140 + zipcode_2016$AGG3E141 + zipcode_2016$AGG3E142 + 
                                zipcode_2016$AGG3E155 + zipcode_2016$AGG3E156 + zipcode_2016$AGG3E157)/zipcode_2016$AF42E046

zipcode_2016$PROP_DIVORCED <- (zipcode_2016$AGG3E092 + zipcode_2016$AGG3E093 + zipcode_2016$AGG3E094 + 
                                 zipcode_2016$AGG3E185 + zipcode_2016$AGG3E186 + zipcode_2016$AGG3E187)/zipcode_2016$AF42E046

zipcode_2016$PROP_WIDOWED <- (zipcode_2016$AGG3E077 + zipcode_2016$AGG3E078 + zipcode_2016$AGG3E079 + 
                                zipcode_2016$AGG3E170 + zipcode_2016$AGG3E171 + zipcode_2016$AGG3E172)/zipcode_2016$AF42E046

zipcode_2016$PROP_HS <- (zipcode_2016$AGIEE036 + zipcode_2016$AGIEE037 + zipcode_2016$AGIEE038 + 
                           zipcode_2016$AGIEE077 + zipcode_2016$AGIEE078 + zipcode_2016$AGIEE079)/zipcode_2016$AF42E046
zipcode_2016$PROP_ASSOC <- (zipcode_2016$AGIEE039 + zipcode_2016$AGIEE040 + zipcode_2016$AGIEE080 +
                              zipcode_2016$AGIEE081)/zipcode_2016$AF42E046
zipcode_2016$PROP_BACHELORS <- (zipcode_2016$AGIEE041 + zipcode_2016$AGIEE042 + zipcode_2016$AGIEE082 +
                                  zipcode_2016$AGIEE083)/zipcode_2016$AF42E046

zipcode_2016$PROP_POVERTY <- (zipcode_2016$AGI6E015 + zipcode_2016$AGI6E016 + zipcode_2016$AGI6E029 + 
                                zipcode_2016$AGI6E030)/zipcode_2016$AF42E046

zipcode_2016$PROP_DISABILITY <- (zipcode_2016$AGKIE016 + zipcode_2016$AGKIE019 + zipcode_2016$AGKIE035 + 
                                   zipcode_2016$AGKIE038)/zipcode_2016$AF90E051
zipcode_2016$PROP_COGNITIVE <- (zipcode_2016$AGKUE013 + zipcode_2016$AGKUE016 + zipcode_2016$AGKUE029 + 
                                  zipcode_2016$AGKUE032)/zipcode_2016$AF90E051
zipcode_2016$PROP_AMBULATORY <- (zipcode_2016$AGKVE013 + zipcode_2016$AGKVE016 + zipcode_2016$AGKVE029 + 
                                   zipcode_2016$AGKVE032)/zipcode_2016$AF90E051

zipcode_2016$PROP_MEDICARE <- (zipcode_2016$AF90E055 + zipcode_2016$AF90E060 + zipcode_2016$AF90E061 + 
                                 zipcode_2016$AF90E062)/zipcode_2016$AF90E051

zipcode_2016$MED_INC <- zipcode_2016$AF49E001
zipcode_2016$MED_INC_65 <- zipcode_2016$AF5LE005
zipcode_2016$GINI <- zipcode_2016$AGL1E001

# List of pertinent ACS variable codes 2017
acs_codes_2017_1 <- c("GISJOIN","YEAR","ZCTA5A","AHYQE001","AH1IE046","AHY2E002","AHY2E003",
                      "AHY2E004","AHY2E005","AHY2E006","AHY2E007","AHY2E008","AHZAE003",
                      "AHZAE004","AHZAE005","AHZAE006","AHZAE012","AHZAE013","AHZAE014",
                      "AHZAE015","AHZAE016","AH1IE047","AH1IE049","AH1IE050","AH1IE052",
                      "AH1IE054","AH1IE055","AH1IE057","AH1IE059","AH1IE060","AH1IE062",
                      "AH1IE064","AH1IE065","AH1IE067","AH1PE001","AH11E005","AH19E001",
                      "AH19E002","AH2RE001","AH3FE017","AH3FE020","AH3FE035","AH3FE038",
                      "AH6IE055","AH6IE060","AH6IE061","AH6IE062","AH6IE051")

acs_codes_2017_2 <- c("ZCTA5A","AH8YE004","AH8YE009","AH8YE015","AH8YE020","AIDLE015",
                      "AIDLE016","AIDLE017","AIDLE108","AIDLE109","AIDLE110","AIDLE031",
                      "AIDLE032","AIDLE033","AIDLE047","AIDLE048","AIDLE049","AIDLE062",
                      "AIDLE063","AIDLE064","AIDLE124","AIDLE125","AIDLE126","AIDLE140",
                      "AIDLE141","AIDLE142","AIDLE155","AIDLE156","AIDLE157","AIDLE092",
                      "AIDLE093","AIDLE094","AIDLE185","AIDLE186","AIDLE187","AIDLE077",
                      "AIDLE078","AIDLE079","AIDLE170","AIDLE171","AIDLE172","AIEWE036",
                      "AIEWE037","AIEWE038","AIEWE077","AIEWE078","AIEWE079","AIEWE039",
                      "AIEWE040","AIEWE080","AIEWE081","AIEWE041","AIEWE042","AIEWE082",
                      "AIEWE083","AIFOE015","AIFOE016","AIFOE029","AIFOE030","AIG0E016",
                      "AIG0E019","AIG0E035","AIG0E038","AIHCE013","AIHCE016","AIHCE029",
                      "AIHCE032","AIHDE013","AIHDE016","AIHDE029","AIHDE032","AIIJE001")

# 5-year ACS 2013-2017
zipcode_2017_1 <- read.csv("/Volumes/PACS$/Cordes, Jack/Zip Code Data/nhgis0046_csv/nhgis0046_ds233_20175_zcta.csv")
zipcode_2017_1 <- read.csv("/Volumes/PACS$/Cordes, Jack/Zip Code Data/nhgis0046_csv/nhgis0046_ds233_20175_zcta.csv")[2:nrow(zipcode_2017_1), acs_codes_2017_1]
zipcode_2017_1[, 4:length(zipcode_2017_1)] <- sapply(zipcode_2017_1[, 4:length(zipcode_2017_1)], as.numeric)

zipcode_2017_2 <- read.csv("/Volumes/PACS$/Cordes, Jack/Zip Code Data/nhgis0046_csv/nhgis0046_ds234_20175_zcta.csv")
zipcode_2017_2 <- read.csv("/Volumes/PACS$/Cordes, Jack/Zip Code Data/nhgis0046_csv/nhgis0046_ds234_20175_zcta.csv")[2:nrow(zipcode_2017_2), acs_codes_2017_2]
zipcode_2017_2[, 2:length(zipcode_2017_2)] <- sapply(zipcode_2017_2[, 2:length(zipcode_2017_2)], as.numeric)

zipcode_2017 <- left_join(zipcode_2017_1, zipcode_2017_2, by = "ZCTA5A")

zipcode_2017$PROP_WHITE <- zipcode_2017$AHY2E002/zipcode_2017$AHYQE001
zipcode_2017$PROP_BLACK <- zipcode_2017$AHY2E003/zipcode_2017$AHYQE001
zipcode_2017$PROP_AIAN <- zipcode_2017$AHY2E004/zipcode_2017$AHYQE001
zipcode_2017$PROP_ASIAN <- zipcode_2017$AHY2E005/zipcode_2017$AHYQE001
zipcode_2017$PROP_NHP <- zipcode_2017$AHY2E006/zipcode_2017$AHYQE001
zipcode_2017$PROP_OTHER_RACE <- zipcode_2017$AHY2E007/zipcode_2017$AHYQE001
zipcode_2017$PROP_2MORE_RACE <- zipcode_2017$AHY2E008/zipcode_2017$AHYQE001

zipcode_2017$PROP_NONHISP_WHITE <- zipcode_2017$AHZAE003/zipcode_2017$AHYQE001
zipcode_2017$PROP_NONHISP_BLACK <- zipcode_2017$AHZAE004/zipcode_2017$AHYQE001
zipcode_2017$PROP_NONHISP_AIAN <- zipcode_2017$AHZAE005/zipcode_2017$AHYQE001
zipcode_2017$PROP_NONHISP_ASIAN <- zipcode_2017$AHZAE006/zipcode_2017$AHYQE001
zipcode_2017$PROP_HISP <- zipcode_2017$AHZAE012/zipcode_2017$AHYQE001
zipcode_2017$PROP_HISP_WHITE <- zipcode_2017$AHZAE013/zipcode_2017$AHYQE001
zipcode_2017$PROP_HISP_BLACK <- zipcode_2017$AHZAE014/zipcode_2017$AHYQE001
zipcode_2017$PROP_HISP_AIAN <- zipcode_2017$AHZAE015/zipcode_2017$AHYQE001
zipcode_2017$PROP_HISP_ASIAN <- zipcode_2017$AHZAE016/zipcode_2017$AHYQE001

zipcode_2017$PROP_ENGLISH_ONLY <- zipcode_2017$AH1IE047/zipcode_2017$AH1IE046
zipcode_2017$PROP_ENGLISH_WELL <- (zipcode_2017$AH1IE049 + zipcode_2017$AH1IE050 + 
                                     zipcode_2017$AH1IE054 + zipcode_2017$AH1IE055 + 
                                     zipcode_2017$AH1IE059 + zipcode_2017$AH1IE060 + 
                                     zipcode_2017$AH1IE064 + zipcode_2017$AH1IE065)/zipcode_2017$AH1IE046
zipcode_2017$PROP_NO_ENGLISH <- (zipcode_2017$AH1IE052 + zipcode_2017$AH1IE057 +
                                   zipcode_2017$AH1IE062 + zipcode_2017$AH1IE067)/zipcode_2017$AH1IE046

zipcode_2017$PROP_PUBLIC_ASST <- zipcode_2017$AH19E002/zipcode_2017$AH19E001

zipcode_2017$PROP_VET <- (zipcode_2017$AH3FE017 + zipcode_2017$AH3FE020 + zipcode_2017$AH3FE035 + 
                            zipcode_2017$AH3FE038)/zipcode_2017$AH1IE046

zipcode_2017$PROP_NATIVE <- (zipcode_2017$AH8YE004 + zipcode_2017$AH8YE009 + zipcode_2017$AH8YE015 +
                               zipcode_2017$AH8YE020)/zipcode_2017$AHYQE001

zipcode_2017$PROP_NEVER_MARRIED <- (zipcode_2017$AIDLE015 + zipcode_2017$AIDLE016 + zipcode_2017$AIDLE017 +
                                      zipcode_2017$AIDLE108 + zipcode_2017$AIDLE109 + zipcode_2017$AIDLE110)/zipcode_2017$AH1IE046

zipcode_2017$PROP_MARRIED <- (zipcode_2017$AIDLE031 + zipcode_2017$AIDLE032 + zipcode_2017$AIDLE033 + 
                                zipcode_2017$AIDLE047 + zipcode_2017$AIDLE048 + zipcode_2017$AIDLE049 + 
                                zipcode_2017$AIDLE062 + zipcode_2017$AIDLE063 + zipcode_2017$AIDLE064 + 
                                zipcode_2017$AIDLE124 + zipcode_2017$AIDLE125 + zipcode_2017$AIDLE126 + 
                                zipcode_2017$AIDLE140 + zipcode_2017$AIDLE141 + zipcode_2017$AIDLE142 + 
                                zipcode_2017$AIDLE155 + zipcode_2017$AIDLE156 + zipcode_2017$AIDLE157)/zipcode_2017$AH1IE046

zipcode_2017$PROP_DIVORCED <- (zipcode_2017$AIDLE092 + zipcode_2017$AIDLE093 + zipcode_2017$AIDLE094 + 
                                 zipcode_2017$AIDLE185 + zipcode_2017$AIDLE186 + zipcode_2017$AIDLE187)/zipcode_2017$AH1IE046

zipcode_2017$PROP_WIDOWED <- (zipcode_2017$AIDLE077 + zipcode_2017$AIDLE078 + zipcode_2017$AIDLE079 + 
                                zipcode_2017$AIDLE170 + zipcode_2017$AIDLE171 + zipcode_2017$AIDLE172)/zipcode_2017$AH1IE046

zipcode_2017$PROP_HS <- (zipcode_2017$AIEWE036 + zipcode_2017$AIEWE037 + zipcode_2017$AIEWE038 + 
                           zipcode_2017$AIEWE077 + zipcode_2017$AIEWE078 + zipcode_2017$AIEWE079)/zipcode_2017$AH1IE046
zipcode_2017$PROP_ASSOC <- (zipcode_2017$AIEWE039 + zipcode_2017$AIEWE040 + zipcode_2017$AIEWE080 +
                              zipcode_2017$AIEWE081)/zipcode_2017$AH1IE046
zipcode_2017$PROP_BACHELORS <- (zipcode_2017$AIEWE041 + zipcode_2017$AIEWE042 + zipcode_2017$AIEWE082 +
                                  zipcode_2017$AIEWE083)/zipcode_2017$AH1IE046

zipcode_2017$PROP_POVERTY <- (zipcode_2017$AIFOE015 + zipcode_2017$AIFOE016 + zipcode_2017$AIFOE029 + 
                                zipcode_2017$AIFOE030)/zipcode_2017$AH1IE046

zipcode_2017$PROP_DISABILITY <- (zipcode_2017$AIG0E016 + zipcode_2017$AIG0E019 + zipcode_2017$AIG0E035 + 
                                   zipcode_2017$AIG0E038)/zipcode_2017$AH6IE051
zipcode_2017$PROP_COGNITIVE <- (zipcode_2017$AIHCE013 + zipcode_2017$AIHCE016 + zipcode_2017$AIHCE029 + 
                                  zipcode_2017$AIHCE032)/zipcode_2017$AH6IE051
zipcode_2017$PROP_AMBULATORY <- (zipcode_2017$AIHDE013 + zipcode_2017$AIHDE016 + zipcode_2017$AIHDE029 + 
                                   zipcode_2017$AIHDE032)/zipcode_2017$AH6IE051

zipcode_2017$PROP_MEDICARE <- (zipcode_2017$AH6IE055 + zipcode_2017$AH6IE060 + zipcode_2017$AH6IE061 + 
                                 zipcode_2017$AH6IE062)/zipcode_2017$AH6IE051

zipcode_2017$MED_INC <- zipcode_2017$AH1PE001
zipcode_2017$MED_INC_65 <- zipcode_2017$AH11E005
zipcode_2017$GINI <- zipcode_2017$AIIJE001

# 5-year ACS All Years 2012-2017
zipcode_2012_2017 <- rbind(zipcode_2012[,c(1:3,(length(zipcode_2012)-36):length(zipcode_2012))], 
                           zipcode_2013[,c(1:3,(length(zipcode_2013)-36):length(zipcode_2013))], 
                           zipcode_2014[,c(1:3,(length(zipcode_2014)-36):length(zipcode_2014))], 
                           zipcode_2015[,c(1:3,(length(zipcode_2015)-36):length(zipcode_2015))], 
                           zipcode_2016[,c(1:3,(length(zipcode_2016)-36):length(zipcode_2016))],
                           zipcode_2017[,c(1:3,(length(zipcode_2017)-36):length(zipcode_2017))])

zipcode_2012_2017$YEAR <- as.numeric(substr(zipcode_2012_2017$YEAR, 6, 10))
zipcode_2012_2017$ZCTA_USE <- zipcode_2012_2017$ZCTA5A

# Load zip code to ZCTA crosswalk files (source: UDS Mapper) ####
zip_zcta_crosswalk_2012 <- read.csv("/Volumes/PACS$/Cordes, Jack/Zip Code Data/Zip Code to ZCTA Crosswalks/ZIPCodetoZCTACrosswalk2012.csv",
                                    colClasses = c(ZIP = "character", ZCTA_USE = "character"))
zip_zcta_crosswalk_2013 <- read.csv("/Volumes/PACS$/Cordes, Jack/Zip Code Data/Zip Code to ZCTA Crosswalks/ZIPCodetoZCTACrosswalk2013.csv",
                                    colClasses = c(ZIP = "character", ZCTA_USE = "character"))
zip_zcta_crosswalk_2014 <- read.csv("/Volumes/PACS$/Cordes, Jack/Zip Code Data/Zip Code to ZCTA Crosswalks/ZIPCodetoZCTACrosswalk2014.csv",
                                    colClasses = c(ZIP = "character", ZCTA_USE = "character"))
zip_zcta_crosswalk_2015 <- read.csv("/Volumes/PACS$/Cordes, Jack/Zip Code Data/Zip Code to ZCTA Crosswalks/ZIPCodetoZCTACrosswalk2015.csv",
                                    colClasses = c(ZIP = "character", ZCTA = "character"))
zip_zcta_crosswalk_2015$ZCTA_USE <- zip_zcta_crosswalk_2015$ZCTA
zip_zcta_crosswalk_2016 <- read.csv("/Volumes/PACS$/Cordes, Jack/Zip Code Data/Zip Code to ZCTA Crosswalks/ZIPCodetoZCTACrosswalk2016.csv",
                                    colClasses = c(ZIP = "character", ZCTA_USE = "character"))
zip_zcta_crosswalk_2017 <- read.csv("/Volumes/PACS$/Cordes, Jack/Zip Code Data/Zip Code to ZCTA Crosswalks/ZIPCodetoZCTACrosswalk2017.csv",
                                    colClasses = c(ZIP_CODE = "character", ZCTA = "character"))
zip_zcta_crosswalk_2017$ZIP <- zip_zcta_crosswalk_2017$ZIP_CODE
zip_zcta_crosswalk_2017$ZCTA_USE <- zip_zcta_crosswalk_2017$ZCTA

zip_zcta_crosswalk_2012 <- zip_zcta_crosswalk_2012 %>% select(c(ZIP, ZCTA_USE))
zip_zcta_crosswalk_2013 <- zip_zcta_crosswalk_2012 %>% select(c(ZIP, ZCTA_USE))
zip_zcta_crosswalk_2014 <- zip_zcta_crosswalk_2012 %>% select(c(ZIP, ZCTA_USE))
zip_zcta_crosswalk_2015 <- zip_zcta_crosswalk_2012 %>% select(c(ZIP, ZCTA_USE))
zip_zcta_crosswalk_2016 <- zip_zcta_crosswalk_2012 %>% select(c(ZIP, ZCTA_USE))
zip_zcta_crosswalk_2017 <- zip_zcta_crosswalk_2012 %>% select(c(ZIP, ZCTA_USE))

zip_zcta_crosswalk_2012$YEAR <- 2012
zip_zcta_crosswalk_2013$YEAR <- 2013
zip_zcta_crosswalk_2014$YEAR <- 2014
zip_zcta_crosswalk_2015$YEAR <- 2015
zip_zcta_crosswalk_2016$YEAR <- 2016
zip_zcta_crosswalk_2017$YEAR <- 2017

zip_zcta_crosswalk <- rbind(zip_zcta_crosswalk_2012, zip_zcta_crosswalk_2013,
                            zip_zcta_crosswalk_2014, zip_zcta_crosswalk_2015,
                            zip_zcta_crosswalk_2016, zip_zcta_crosswalk_2017)

# Load shapefiles ####
zcta_shp <- rgdal::readOGR("/Volumes/PACS$/Cordes, Jack/Zip Code Data/nhgis0046_shape/nhgis0046_shapefile_tl2017_us_zcta_2017", layer = "US_zcta_2017",
                           verbose = F)

#zcta_sf <- as(zcta_shp, "sf")

#zcta_sf <- st_sfc(zcta_sf, crs = 2163)

#zcta_sf <- zcta_sf %>% st_set_crs(2831) %>% st_transform(crs=4617)

state_shp <- rgdal::readOGR("/Volumes/PACS$/Cordes, Jack/Zip Code Data/nhgis0047_shape/nhgis0047_shapefile_tl2017_us_state_2017", layer = "US_state_2017",
                            verbose = F)

# TECOS Primary ####
tecos <- read.csv("/Volumes/PACS$/Cordes, Jack/Data_without_Prescriber_ID/JCA01_TECOS_primary_all_followup.csv", colClasses = c(ZIP_CD = "character"))

tecos$OUTCOME <- ifelse(tecos$OUTCOME == "false", 0, 1)

tecos_exposure <- read.csv("/Volumes/PACS$/Cordes, Jack/TECOS/outcome_Primary_Composite_Outcome__Angina__stroke__MI__all-cause_mortality_from_vital_/primary/all_followup.csv")[,1:2]

tecos_exposure$EXPOSURE <- ifelse(tecos_exposure$EXPOSURE == "R", 0, 1)

tecos <- left_join(tecos, tecos_exposure, by = "PID")

tecos$YEAR <- as.numeric(substr(tecos$ENTRY_CLASS, 1,4))
tecos$ZIP <- tecos$ZIP_CD

# Add prescriber ID
tecos_rx_id <- read.csv("/Volumes/PACS$/Cordes, Jack/Files with Zip Codes/jca02_TECOS_primary_all_followup.csv")[,c(1,217:219)]

tecos <- left_join(tecos, tecos_rx_id, by = "PID")

# Convert cohort zip codes to ZCTA
tecos <- merge(tecos, zip_zcta_crosswalk, by=c("ZIP","YEAR"))

# Join ACS zip code data to cohort
tecos <- merge(tecos, zipcode_2012_2017, by=c("ZCTA_USE","YEAR"))

# Check for missingness in ACS data and remove individuals with missing zipcode information
colSums(is.na(tecos))

# Remove individuals younger than 65 years old
# Total of 6935 individuals comprising 1.69% of the cohort
# 99 individuals had both missing zipcode data and were <65 years old
tecos <- tecos %>% filter(RUN1_ENTRY_COVARIATE_147 >= 65)

# Remove individuals with missing zipcode ACS data
# Total of 13050 individuals comprising 3.18% of the cohort
tecos <- tecos[complete.cases(tecos), ]

# Remove individuals with zipcodes from Puerto Rico
# Total of 24 individuals comprising 0.0062% of the cohort
tecos <- tecos %>% filter(as.numeric(ZCTA_USE) >= 1000)

# Rename Aetion variables
tecos <- tecos %>% rename(AGE_CAT	=	RUN1_ENTRY_COVARIATE_1,
                          SEX	=	RUN1_ENTRY_COVARIATE_2,
                          RACE	=	RUN1_ENTRY_COVARIATE_3,
                          REGION	=	RUN1_ENTRY_COVARIATE_4,
                          OBESITY	=	RUN1_ENTRY_COVARIATE_5,
                          OVERWEIGHT	=	RUN1_ENTRY_COVARIATE_6,
                          SMOKE	=	RUN1_ENTRY_COVARIATE_7,
                          ALCOHOL	=	RUN1_ENTRY_COVARIATE_8,
                          DRUG	=	RUN1_ENTRY_COVARIATE_9,
                          DM_RETINOPATHY	=	RUN1_ENTRY_COVARIATE_10,
                          DM_OPHTHALMIC	=	RUN1_ENTRY_COVARIATE_11,
                          RETINAL_DETACH	=	RUN1_ENTRY_COVARIATE_12,
                          RETINAL_LASER	=	RUN1_ENTRY_COVARIATE_13,
                          DM_NEUROPATHY2	=	RUN1_ENTRY_COVARIATE_14,
                          DM_NEUROPATHY3	=	RUN1_ENTRY_COVARIATE_15,
                          HYPOGLYCEMIA	=	RUN1_ENTRY_COVARIATE_16,
                          HYPERGLYCEMIA	=	RUN1_ENTRY_COVARIATE_17,
                          ELECTROLYTE_ACID_BASE	=	RUN1_ENTRY_COVARIATE_18,
                          DM_KETOACIDOSIS	=	RUN1_ENTRY_COVARIATE_19,
                          HONK	=	RUN1_ENTRY_COVARIATE_20,
                          DM_PERIPHERAL	=	RUN1_ENTRY_COVARIATE_21,
                          DM_FOOT	=	RUN1_ENTRY_COVARIATE_22,
                          GANGRENE	=	RUN1_ENTRY_COVARIATE_23,
                          LOWER_AMPUTATION	=	RUN1_ENTRY_COVARIATE_24,
                          OSTEOMYELITIS	=	RUN1_ENTRY_COVARIATE_25,
                          SKIN_INFECTION	=	RUN1_ENTRY_COVARIATE_26,
                          ED	=	RUN1_ENTRY_COVARIATE_27,
                          DM_COMPLICATION	=	RUN1_ENTRY_COVARIATE_28,
                          DM_NO_COMPLICATION	=	RUN1_ENTRY_COVARIATE_29,
                          HYPERTENSION	=	RUN1_ENTRY_COVARIATE_30,
                          HYPERLIPIDEMIA	=	RUN1_ENTRY_COVARIATE_31,
                          IHD	=	RUN1_ENTRY_COVARIATE_32,
                          ACUTE_MI	=	RUN1_ENTRY_COVARIATE_33,
                          UNSTABLE_ANGINA	=	RUN1_ENTRY_COVARIATE_34,
                          OLD_MI	=	RUN1_ENTRY_COVARIATE_35,
                          STABLE_ANGINA	=	RUN1_ENTRY_COVARIATE_36,
                          CORONARY_ATHEROSCLEROSIS	=	RUN1_ENTRY_COVARIATE_37,
                          OTHER_ATHEROSCLEROSIS	=	RUN1_ENTRY_COVARIATE_38,
                          CARDIAC_PROCEDURE	=	RUN1_ENTRY_COVARIATE_39,
                          CABG_PTCA	=	RUN1_ENTRY_COVARIATE_40,
                          STROKE	=	RUN1_ENTRY_COVARIATE_41,
                          ISCHEMIC_STROKE	=	RUN1_ENTRY_COVARIATE_42,
                          HEMORRHAGIC_STROKE	=	RUN1_ENTRY_COVARIATE_43,
                          TIA	=	RUN1_ENTRY_COVARIATE_44,
                          OTHER_CEREBROVASCULAR	=	RUN1_ENTRY_COVARIATE_45,
                          LATE_CEREBROVASCULAR	=	RUN1_ENTRY_COVARIATE_46,
                          CEREBROVASCULAR_PROCEDURE	=	RUN1_ENTRY_COVARIATE_47,
                          CHF	=	RUN1_ENTRY_COVARIATE_48,
                          PVD	=	RUN1_ENTRY_COVARIATE_49,
                          ATRIAL_FIBRILLATION	=	RUN1_ENTRY_COVARIATE_50,
                          CARDIAC_DYSRHYTHMIA	=	RUN1_ENTRY_COVARIATE_51,
                          CARDIAC_CONDUCTION	=	RUN1_ENTRY_COVARIATE_52,
                          OTHER_CVD	=	RUN1_ENTRY_COVARIATE_53,
                          EDEMA	=	RUN1_ENTRY_COVARIATE_54,
                          COPD	=	RUN1_ENTRY_COVARIATE_55,
                          ASTHMA	=	RUN1_ENTRY_COVARIATE_56,
                          OBS_SLEEP_APNEA	=	RUN1_ENTRY_COVARIATE_57,
                          PNEUMONIA	=	RUN1_ENTRY_COVARIATE_58,
                          RENAL_DYSFUNCTION	=	RUN1_ENTRY_COVARIATE_59,
                          ACUTE_RENAL_DISEASE	=	RUN1_ENTRY_COVARIATE_60,
                          CHRONIC_RENAL_INSUFFICIENCY	=	RUN1_ENTRY_COVARIATE_61,
                          CKD	=	RUN1_ENTRY_COVARIATE_62,
                          CKD_3_4	=	RUN1_ENTRY_COVARIATE_63,
                          HYPERTENSIVE_NEPHROPATHY	=	RUN1_ENTRY_COVARIATE_64,
                          OTHER_RENAL_INSUFFICIENCY	=	RUN1_ENTRY_COVARIATE_65,
                          LIVER_DISEASE	=	RUN1_ENTRY_COVARIATE_66,
                          OSTEOARTHRITIS	=	RUN1_ENTRY_COVARIATE_67,
                          OTHER_ARTHRITIS	=	RUN1_ENTRY_COVARIATE_68,
                          DORSOPATHIES	=	RUN1_ENTRY_COVARIATE_69,
                          FRACTURE	=	RUN1_ENTRY_COVARIATE_70,
                          FALL	=	RUN1_ENTRY_COVARIATE_71,
                          OSTEOPOROSIS	=	RUN1_ENTRY_COVARIATE_72,
                          HYPERTHYROIDISM	=	RUN1_ENTRY_COVARIATE_73,
                          HYPOTHYROIDISM	=	RUN1_ENTRY_COVARIATE_74,
                          OTHER_THYROID	=	RUN1_ENTRY_COVARIATE_75,
                          DEPRESSION	=	RUN1_ENTRY_COVARIATE_76,
                          ANXIETY	=	RUN1_ENTRY_COVARIATE_77,
                          SLEEP_DISORDER	=	RUN1_ENTRY_COVARIATE_78,
                          DEMENTIA	=	RUN1_ENTRY_COVARIATE_79,
                          DELIRIUM	=	RUN1_ENTRY_COVARIATE_80,
                          PSYCHOSIS	=	RUN1_ENTRY_COVARIATE_81,
                          FRAILTY_QUALITATIVE	=	RUN1_ENTRY_COVARIATE_82,
                          FRAILTY_EMPIRICAL_V3	=	RUN1_ENTRY_COVARIATE_83,
                          NON_FRAILTY	=	RUN1_ENTRY_COVARIATE_84,
                          N_ANTIDIABETICS	=	RUN1_ENTRY_COVARIATE_85,
                          NAIVE_NEW_USER	=	RUN1_ENTRY_COVARIATE_86,
                          ACE_INHIBITORS	=	RUN1_ENTRY_COVARIATE_87,
                          ARB	=	RUN1_ENTRY_COVARIATE_88,
                          LOOP_DIURETICS	=	RUN1_ENTRY_COVARIATE_89,
                          OTHER_DIURETICS	=	RUN1_ENTRY_COVARIATE_90,
                          NITRATES	=	RUN1_ENTRY_COVARIATE_91,
                          OTHER_ANTIHYPERTENSIVE	=	RUN1_ENTRY_COVARIATE_92,
                          DIGOXIN	=	RUN1_ENTRY_COVARIATE_93,
                          ANTI_ARRHYTHMICS	=	RUN1_ENTRY_COVARIATE_94,
                          COPD_ASTHMA_DRUG	=	RUN1_ENTRY_COVARIATE_95,
                          STATIN	=	RUN1_ENTRY_COVARIATE_96,
                          OTHER_LIPID_DRUG	=	RUN1_ENTRY_COVARIATE_97,
                          ANTIPLATELET	=	RUN1_ENTRY_COVARIATE_98,
                          ORAL_ANTICOAGULANT	=	RUN1_ENTRY_COVARIATE_99,
                          HEPARIN	=	RUN1_ENTRY_COVARIATE_100,
                          NSAID	=	RUN1_ENTRY_COVARIATE_101,
                          ORAL_CORTICOSTEROID	=	RUN1_ENTRY_COVARIATE_102,
                          BISPHOSPHONATE	=	RUN1_ENTRY_COVARIATE_103,
                          OPIOID	=	RUN1_ENTRY_COVARIATE_104,
                          ANTIDEPRESSANT	=	RUN1_ENTRY_COVARIATE_105,
                          ANTIPSYCHOTIC	=	RUN1_ENTRY_COVARIATE_106,
                          ANTICONVULSANT	=	RUN1_ENTRY_COVARIATE_107,
                          LITHIUM	=	RUN1_ENTRY_COVARIATE_108,
                          BENZODIAZAPINE	=	RUN1_ENTRY_COVARIATE_109,
                          ANXIOLYTIC_HYPNOTIC	=	RUN1_ENTRY_COVARIATE_110,
                          DEMENTIA_DRUG	=	RUN1_ENTRY_COVARIATE_111,
                          PARKINSON_DRUG	=	RUN1_ENTRY_COVARIATE_112,
                          N_DIAGNOSES	=	RUN1_ENTRY_COVARIATE_113,
                          N_DRUG_RX	=	RUN1_ENTRY_COVARIATE_114,
                          HOSPITALIZATION	=	RUN1_ENTRY_COVARIATE_115,
                          ENDOCRINOLOGIST	=	RUN1_ENTRY_COVARIATE_116,
                          INTERNAL	=	RUN1_ENTRY_COVARIATE_117,
                          CARDIOLOGIST	=	RUN1_ENTRY_COVARIATE_118,
                          ELECTROCARDIOGRAM	=	RUN1_ENTRY_COVARIATE_119,
                          GLUCOSE_TESTS	=	RUN1_ENTRY_COVARIATE_120,
                          HOSPITALIZATION_30	=	RUN1_ENTRY_COVARIATE_121,
                          HOSPITALIZATION_31_180	=	RUN1_ENTRY_COVARIATE_122,
                          N_HOSPITALIZATION	=	RUN1_ENTRY_COVARIATE_123,
                          N_HOSPITAL_DAYS	=	RUN1_ENTRY_COVARIATE_124,
                          N_ED	=	RUN1_ENTRY_COVARIATE_125,
                          N_OFFICE	=	RUN1_ENTRY_COVARIATE_126,
                          N_ENDOCRINOLOGIST	=	RUN1_ENTRY_COVARIATE_127,
                          N_INTERNAL	=	RUN1_ENTRY_COVARIATE_128,
                          N_CARDIOLOGIST	=	RUN1_ENTRY_COVARIATE_129,
                          N_ELECTROCARDIOGRAM	=	RUN1_ENTRY_COVARIATE_130,
                          N_HBA1C	=	RUN1_ENTRY_COVARIATE_131,
                          N_GLUCOSE_TEST	=	RUN1_ENTRY_COVARIATE_132,
                          N_LIPID_TEST	=	RUN1_ENTRY_COVARIATE_133,
                          N_CREATININE_TEST	=	RUN1_ENTRY_COVARIATE_134,
                          N_BUN_TEST	=	RUN1_ENTRY_COVARIATE_135,
                          N_MICROALBUMINURIA_TEST	=	RUN1_ENTRY_COVARIATE_136,
                          DM_DRUG_AGI	=	RUN1_ENTRY_COVARIATE_137,
                          DM_DRUG_GLITAZONE	=	RUN1_ENTRY_COVARIATE_138,
                          DM_DRUG_GLP1	=	RUN1_ENTRY_COVARIATE_139,
                          DM_DRUG_INSULIN	=	RUN1_ENTRY_COVARIATE_140,
                          DM_DRUG_MEGLITINIDE	=	RUN1_ENTRY_COVARIATE_141,
                          DM_DRUG_METFORMIN	=	RUN1_ENTRY_COVARIATE_142,
                          PRAMLINTIDE	=	RUN1_ENTRY_COVARIATE_143,
                          GEN1_SU	=	RUN1_ENTRY_COVARIATE_144,
                          MONOTHERAPY_INITIATION	=	RUN1_ENTRY_COVARIATE_145,
                          METFORMIN_DUAL	=	RUN1_ENTRY_COVARIATE_146,
                          AGE	=	RUN1_ENTRY_COVARIATE_147,
                          CEREBROVASCULAR_HEM_STROKE	=	RUN1_ENTRY_COVARIATE_148,
                          BLADDER_STONE	=	RUN1_ENTRY_COVARIATE_149,
                          KIDNEY_STONE	=	RUN1_ENTRY_COVARIATE_150,
                          UTI	=	RUN1_ENTRY_COVARIATE_151,
                          DIPSTICK_URINALYSIS	=	RUN1_ENTRY_COVARIATE_152,
                          NON_DIPSTICK_URINALYSIS	=	RUN1_ENTRY_COVARIATE_153,
                          URINE_FUNCTION_TEST	=	RUN1_ENTRY_COVARIATE_154,
                          CYTOLOGY	=	RUN1_ENTRY_COVARIATE_155,
                          CYSTOSCOPY	=	RUN1_ENTRY_COVARIATE_156,
                          FRAILTY_EMPIRICAL	=	RUN1_ENTRY_COVARIATE_157,
                          CREATININE_TEST	=	RUN1_ENTRY_COVARIATE_158,
                          BUN_TEST	=	RUN1_ENTRY_COVARIATE_159,
                          CRI_NO_CKD	=	RUN1_ENTRY_COVARIATE_160,
                          CKD_1_2	=	RUN1_ENTRY_COVARIATE_161,
                          CKD_3_6	=	RUN1_ENTRY_COVARIATE_162,
                          CONCOMITANT_SGLT2I	=	RUN1_ENTRY_COVARIATE_163,
                          CONCOMITANT_AGI	=	RUN1_ENTRY_COVARIATE_164,
                          CONCOMITANT_GLITAZONE	=	RUN1_ENTRY_COVARIATE_165,
                          CONCOMITANT_GLP1	=	RUN1_ENTRY_COVARIATE_166,
                          CONCOMITANT_INSULIN	=	RUN1_ENTRY_COVARIATE_167,
                          CONCOMITANT_MEGLITINIDE	=	RUN1_ENTRY_COVARIATE_168,
                          CONCOMITANT_METFORMIN	=	RUN1_ENTRY_COVARIATE_169,
                          FRAILTY_QUALITATIVE_V1	=	RUN1_ENTRY_COVARIATE_170,
                          PAST_SGLT2I	=	RUN1_ENTRY_COVARIATE_171,
                          PAST_AGI	=	RUN1_ENTRY_COVARIATE_172,
                          PAST_GLITAZONE	=	RUN1_ENTRY_COVARIATE_173,
                          PAST_GLP1	=	RUN1_ENTRY_COVARIATE_174,
                          PAST_INSULIN	=	RUN1_ENTRY_COVARIATE_175,
                          PAST_MEGLITINIDE	=	RUN1_ENTRY_COVARIATE_176,
                          PAST_METFORMIN	=	RUN1_ENTRY_COVARIATE_177,
                          CALENDAR_TIME_DAY	=	RUN1_ENTRY_COVARIATE_178,
                          BLADDER_KIDNEY_STONE	=	RUN1_ENTRY_COVARIATE_179,
                          PERIPHERAL_GANGRENE_OSTEOMYELITIS	=	RUN1_ENTRY_COVARIATE_180,
                          AGE_DECILE	=	RUN1_ENTRY_COVARIATE_181,
                          ALCOHOL_DRUG	=	RUN1_ENTRY_COVARIATE_182,
                          DM_EYE	=	RUN1_ENTRY_COVARIATE_183,
                          COMPOSITE_CVD	=	RUN1_ENTRY_COVARIATE_184,
                          COMPOSITE_CARDIAC_PROCEDURE	=	RUN1_ENTRY_COVARIATE_185,
                          THYROID	=	RUN1_ENTRY_COVARIATE_186,
                          DELIRIUM_PSYCHOSIS	=	RUN1_ENTRY_COVARIATE_187,
                          MEGLITINIDE	=	RUN1_ENTRY_COVARIATE_188,
                          AGI	=	RUN1_ENTRY_COVARIATE_189,
                          GLAUCOMA_CATARACTS	=	RUN1_ENTRY_COVARIATE_190,
                          IMAGING	=	RUN1_ENTRY_COVARIATE_191,
                          CELLULITIS_ABCESS_TOE	=	RUN1_ENTRY_COVARIATE_192,
                          FOOT_ULCER	=	RUN1_ENTRY_COVARIATE_193,
                          ENTRESTO	=	RUN1_ENTRY_COVARIATE_194,
                          ENDOCRINOLOGIST_30	=	RUN1_ENTRY_COVARIATE_195,
                          ENDOCRINOLOGIST_31_180	=	RUN1_ENTRY_COVARIATE_196,
                          CARDIOLOGIST_30	=	RUN1_ENTRY_COVARIATE_197,
                          CARDIOLOGIST_31_180	=	RUN1_ENTRY_COVARIATE_198,
                          INTERNAL_30	=	RUN1_ENTRY_COVARIATE_199,
                          INTERNAL_31_180	=	RUN1_ENTRY_COVARIATE_200,
                          DIALYSIS	=	RUN1_ENTRY_COVARIATE_201,
                          CCI_180	=	RUN1_ENTRY_COVARIATE_202,
                          CKD_3_6_DIALYSIS	=	RUN1_ENTRY_COVARIATE_203,
                          BASELINE_CV	=	RUN1_ENTRY_COVARIATE_204,
                          THIAZIDE	=	RUN1_ENTRY_COVARIATE_205,
                          BETA_BLOCKER	=	RUN1_ENTRY_COVARIATE_206,
                          CA_CHANNEL_BLOCKER	=	RUN1_ENTRY_COVARIATE_207)

# Convert all true/false logicals to 1/0 numeric
tecos[,c(9:85, 88, 90:116, 119:126, 141:150, 152:160, 162:173, 175:181, 183:184, 186:205, 207:211)] <- lapply(tecos[,c(9:85, 88, 90:116, 119:126, 141:150, 152:160, 162:173, 175:181, 183:184, 186:205, 207:211)], function(x) as.numeric(as.logical(x)))

# Convert all true/false logicals to 1/0 numeric if no zipcode data included
#tecos[,c(6:82, 85, 87:113, 116:123, 138:147, 149:157, 159:170, 172:178, 180:181, 183:202, 204:208)] <- lapply(tecos[,c(6:82, 85, 87:113, 116:123, 138:147, 149:157, 159:170, 172:178, 180:181, 183:202, 204:208)], function(x) as.numeric(as.logical(x)))

# Create normalized versions of MED_INC and MED_INC_65 to [0, 1] for use in analysis
tecos <- tecos %>% mutate(MED_INC_STD = (MED_INC - min(MED_INC))/(max(MED_INC) - min(MED_INC)))
tecos <- tecos %>% mutate(MED_INC_65_STD =
                            (MED_INC_65 - min(MED_INC_65))/(max(MED_INC_65) - min(MED_INC_65)))

# Add state variable based on ZCTA
tecos <- tecos %>% mutate(
  STATE_FIPS = case_when(
    as.numeric(ZCTA_USE) >= 35000 & as.numeric(ZCTA_USE) < 37000 ~ "01",
    as.numeric(ZCTA_USE) >= 99500 ~ "02",
    as.numeric(ZCTA_USE) >= 85000 & as.numeric(ZCTA_USE) < 87000 ~ "04",
    as.numeric(ZCTA_USE) >= 71600 & as.numeric(ZCTA_USE) < 73000 ~ "05",
    as.numeric(ZCTA_USE) >= 90000 & as.numeric(ZCTA_USE) < 96200 ~ "06",
    as.numeric(ZCTA_USE) >= 80000 & as.numeric(ZCTA_USE) < 82000 ~ "08",
    as.numeric(ZCTA_USE) >= 6000 & as.numeric(ZCTA_USE) < 7000 ~ "09",
    as.numeric(ZCTA_USE) >= 19700 & as.numeric(ZCTA_USE) < 20000 ~ "10",
    as.numeric(ZCTA_USE) >= 20000 & as.numeric(ZCTA_USE) < 20600 ~ "11",
    as.numeric(ZCTA_USE) >= 32000 & as.numeric(ZCTA_USE) < 35000 ~ "12",
    as.numeric(ZCTA_USE) >= 30000 & as.numeric(ZCTA_USE) < 32000 ~ "13",
    as.numeric(ZCTA_USE) >= 96700 & as.numeric(ZCTA_USE) < 96900 ~ "15",
    as.numeric(ZCTA_USE) >= 83200 & as.numeric(ZCTA_USE) < 84000 ~ "16",
    as.numeric(ZCTA_USE) >= 60000 & as.numeric(ZCTA_USE) < 63000 ~ "17",
    as.numeric(ZCTA_USE) >= 46000 & as.numeric(ZCTA_USE) < 48000 ~ "18",
    as.numeric(ZCTA_USE) >= 50000 & as.numeric(ZCTA_USE) < 53000 ~ "19",
    as.numeric(ZCTA_USE) >= 66000 & as.numeric(ZCTA_USE) < 68000 ~ "20",
    as.numeric(ZCTA_USE) >= 40000 & as.numeric(ZCTA_USE) < 43000 ~ "21",
    as.numeric(ZCTA_USE) >= 70000 & as.numeric(ZCTA_USE) < 71600 ~ "22",
    as.numeric(ZCTA_USE) >= 3900 & as.numeric(ZCTA_USE) < 5000 ~ "23",
    as.numeric(ZCTA_USE) >= 20600 & as.numeric(ZCTA_USE) < 22000 ~ "24",
    as.numeric(ZCTA_USE) >= 1000 & as.numeric(ZCTA_USE) < 2800 ~ "25",
    as.numeric(ZCTA_USE) >= 48000 & as.numeric(ZCTA_USE) < 50000 ~ "26",
    as.numeric(ZCTA_USE) >= 55000 & as.numeric(ZCTA_USE) < 56800 ~ "27",
    as.numeric(ZCTA_USE) >= 38600 & as.numeric(ZCTA_USE) < 40000 ~ "28",
    as.numeric(ZCTA_USE) >= 63000 & as.numeric(ZCTA_USE) < 66000 ~ "29",
    as.numeric(ZCTA_USE) >= 59000 & as.numeric(ZCTA_USE) < 60000 ~ "30",
    as.numeric(ZCTA_USE) >= 68000 & as.numeric(ZCTA_USE) < 70000 ~ "31",
    as.numeric(ZCTA_USE) >= 88900 & as.numeric(ZCTA_USE) < 90000 ~ "32",
    as.numeric(ZCTA_USE) >= 3000 & as.numeric(ZCTA_USE) < 3900 ~ "33",
    as.numeric(ZCTA_USE) >= 7000 & as.numeric(ZCTA_USE) < 9000 ~ "34",
    as.numeric(ZCTA_USE) >= 87000 & as.numeric(ZCTA_USE) < 88500 ~ "35",
    as.numeric(ZCTA_USE) >= 10000 & as.numeric(ZCTA_USE) < 15000 ~ "36",
    as.numeric(ZCTA_USE) >= 27000 & as.numeric(ZCTA_USE) < 29000 ~ "37",
    as.numeric(ZCTA_USE) >= 58000 & as.numeric(ZCTA_USE) < 59000 ~ "38",
    as.numeric(ZCTA_USE) >= 43000 & as.numeric(ZCTA_USE) < 46000 ~ "39",
    as.numeric(ZCTA_USE) >= 73000 & as.numeric(ZCTA_USE) < 75000 ~ "40",
    as.numeric(ZCTA_USE) >= 97000 & as.numeric(ZCTA_USE) < 98000 ~ "41",
    as.numeric(ZCTA_USE) >= 15000 & as.numeric(ZCTA_USE) < 19700 ~ "42",
    as.numeric(ZCTA_USE) >= 2800 & as.numeric(ZCTA_USE) < 3000 ~ "44",
    as.numeric(ZCTA_USE) >= 29000 & as.numeric(ZCTA_USE) < 30000 ~ "45",
    as.numeric(ZCTA_USE) >= 57000 & as.numeric(ZCTA_USE) < 58000 ~ "46",
    as.numeric(ZCTA_USE) >= 37000 & as.numeric(ZCTA_USE) < 38600 ~ "47",
    as.numeric(ZCTA_USE) >= 75000 & as.numeric(ZCTA_USE) < 80000 ~ "48",
    as.numeric(ZCTA_USE) >= 84000 & as.numeric(ZCTA_USE) < 85000 ~ "49",
    as.numeric(ZCTA_USE) >= 5000 & as.numeric(ZCTA_USE) < 6000 ~ "50",
    as.numeric(ZCTA_USE) >= 22000 & as.numeric(ZCTA_USE) < 24700 ~ "51",
    as.numeric(ZCTA_USE) >= 98000 & as.numeric(ZCTA_USE) < 99500 ~ "53",
    as.numeric(ZCTA_USE) >= 24700 & as.numeric(ZCTA_USE) < 27000 ~ "54",
    as.numeric(ZCTA_USE) >= 53000 & as.numeric(ZCTA_USE) < 55000 ~ "55",
    TRUE ~ "56"
  )
)

tecos <- tecos %>% mutate(
  STATE = case_when(
    as.numeric(ZCTA_USE) >= 35000 & as.numeric(ZCTA_USE) < 37000 ~ "AL",
    as.numeric(ZCTA_USE) >= 99500 ~ "AK",
    as.numeric(ZCTA_USE) >= 85000 & as.numeric(ZCTA_USE) < 87000 ~ "AZ",
    as.numeric(ZCTA_USE) >= 71600 & as.numeric(ZCTA_USE) < 73000 ~ "AR",
    as.numeric(ZCTA_USE) >= 90000 & as.numeric(ZCTA_USE) < 96200 ~ "CA",
    as.numeric(ZCTA_USE) >= 80000 & as.numeric(ZCTA_USE) < 82000 ~ "CO",
    as.numeric(ZCTA_USE) >= 6000 & as.numeric(ZCTA_USE) < 7000 ~ "CT",
    as.numeric(ZCTA_USE) >= 19700 & as.numeric(ZCTA_USE) < 20000 ~ "DE",
    (as.numeric(ZCTA_USE) >= 20000 & as.numeric(ZCTA_USE) < 20100) | (as.numeric(ZCTA_USE) >= 20200 & as.numeric(ZCTA_USE) < 20600) ~ "DC",
    as.numeric(ZCTA_USE) >= 32000 & as.numeric(ZCTA_USE) < 35000 ~ "FL",
    (as.numeric(ZCTA_USE) >= 30000 & as.numeric(ZCTA_USE) < 32000) | (as.numeric(ZCTA_USE) >= 39800 & as.numeric(ZCTA_USE) < 40000) ~ "GA",
    as.numeric(ZCTA_USE) >= 96700 & as.numeric(ZCTA_USE) < 96900 ~ "HI",
    as.numeric(ZCTA_USE) >= 83200 & as.numeric(ZCTA_USE) < 84000 ~ "ID",
    as.numeric(ZCTA_USE) >= 60000 & as.numeric(ZCTA_USE) < 63000 ~ "IL",
    as.numeric(ZCTA_USE) >= 46000 & as.numeric(ZCTA_USE) < 48000 ~ "IN",
    as.numeric(ZCTA_USE) >= 50000 & as.numeric(ZCTA_USE) < 53000 ~ "IA",
    as.numeric(ZCTA_USE) >= 66000 & as.numeric(ZCTA_USE) < 68000 ~ "KS",
    as.numeric(ZCTA_USE) >= 40000 & as.numeric(ZCTA_USE) < 43000 ~ "KY",
    as.numeric(ZCTA_USE) >= 70000 & as.numeric(ZCTA_USE) < 71600 ~ "LA",
    as.numeric(ZCTA_USE) >= 3900 & as.numeric(ZCTA_USE) < 5000 ~ "ME",
    as.numeric(ZCTA_USE) >= 20600 & as.numeric(ZCTA_USE) < 22000 ~ "MD",
    as.numeric(ZCTA_USE) >= 1000 & as.numeric(ZCTA_USE) < 2800 ~ "MA",
    as.numeric(ZCTA_USE) >= 48000 & as.numeric(ZCTA_USE) < 50000 ~ "MI",
    as.numeric(ZCTA_USE) >= 55000 & as.numeric(ZCTA_USE) < 56800 ~ "MN",
    (as.numeric(ZCTA_USE) >= 38600 & as.numeric(ZCTA_USE) < 39800) ~ "MS",
    as.numeric(ZCTA_USE) >= 63000 & as.numeric(ZCTA_USE) < 66000 ~ "MO",
    as.numeric(ZCTA_USE) >= 59000 & as.numeric(ZCTA_USE) < 60000 ~ "MT",
    as.numeric(ZCTA_USE) >= 68000 & as.numeric(ZCTA_USE) < 70000 ~ "NE",
    as.numeric(ZCTA_USE) >= 88900 & as.numeric(ZCTA_USE) < 90000 ~ "NV",
    as.numeric(ZCTA_USE) >= 3000 & as.numeric(ZCTA_USE) < 3900 ~ "NH",
    as.numeric(ZCTA_USE) >= 7000 & as.numeric(ZCTA_USE) < 9000 ~ "NJ",
    as.numeric(ZCTA_USE) >= 87000 & as.numeric(ZCTA_USE) < 88500 ~ "NM",
    as.numeric(ZCTA_USE) >= 10000 & as.numeric(ZCTA_USE) < 15000 ~ "NY",
    as.numeric(ZCTA_USE) >= 27000 & as.numeric(ZCTA_USE) < 29000 ~ "NC",
    as.numeric(ZCTA_USE) >= 58000 & as.numeric(ZCTA_USE) < 59000 ~ "ND",
    as.numeric(ZCTA_USE) >= 43000 & as.numeric(ZCTA_USE) < 46000 ~ "OH",
    as.numeric(ZCTA_USE) >= 73000 & as.numeric(ZCTA_USE) < 75000 ~ "OK",
    as.numeric(ZCTA_USE) >= 97000 & as.numeric(ZCTA_USE) < 98000 ~ "OR",
    as.numeric(ZCTA_USE) >= 15000 & as.numeric(ZCTA_USE) < 19700 ~ "PA",
    as.numeric(ZCTA_USE) >= 2800 & as.numeric(ZCTA_USE) < 3000 ~ "RI",
    as.numeric(ZCTA_USE) >= 29000 & as.numeric(ZCTA_USE) < 30000 ~ "SC",
    as.numeric(ZCTA_USE) >= 57000 & as.numeric(ZCTA_USE) < 58000 ~ "SD",
    as.numeric(ZCTA_USE) >= 37000 & as.numeric(ZCTA_USE) < 38600 ~ "TN",
    as.numeric(ZCTA_USE) >= 75000 & as.numeric(ZCTA_USE) < 80000 ~ "TX",
    as.numeric(ZCTA_USE) >= 84000 & as.numeric(ZCTA_USE) < 85000 ~ "UT",
    as.numeric(ZCTA_USE) >= 5000 & as.numeric(ZCTA_USE) < 6000 ~ "VT",
    (as.numeric(ZCTA_USE) >= 22000 & as.numeric(ZCTA_USE) < 24700) | (as.numeric(ZCTA_USE) >= 20100 & as.numeric(ZCTA_USE) < 20200) ~ "VA",
    as.numeric(ZCTA_USE) >= 98000 & as.numeric(ZCTA_USE) < 99500 ~ "WA",
    as.numeric(ZCTA_USE) >= 24700 & as.numeric(ZCTA_USE) < 27000 ~ "WV",
    as.numeric(ZCTA_USE) >= 53000 & as.numeric(ZCTA_USE) < 55000 ~ "WI",
    TRUE ~ "WY"
  )
)

# Add region variable based on ZCTA
tecos <- tecos %>% mutate(
  REGION_ZCTA = case_when(
    (STATE == "ME" |
       STATE == "NH" |
       STATE == "VT" |
       STATE == "MA" |
       STATE == "RI" |
       STATE == "CT" |
       STATE == "NY" |
       STATE == "NJ" |
       STATE == "PA") ~ "0",
    (STATE == "DE" |
       STATE == "MD" |
       STATE == "DC" |
       STATE == "WV" |
       STATE == "VA" |
       STATE == "NC" |
       STATE == "SC" |
       STATE == "GA" |
       STATE == "FL" |
       STATE == "AL" |
       STATE == "TN" |
       STATE == "KY" |
       STATE == "MS" |
       STATE == "AR" |
       STATE == "LA" |
       STATE == "OK" |
       STATE == "TX") ~ "1",
    (STATE == "OH" |
       STATE == "MI" |
       STATE == "IN" |
       STATE == "IL" |
       STATE == "WI" |
       STATE == "MN" |
       STATE == "IA" |
       STATE == "MO" |
       STATE == "KS" |
       STATE == "NE" |
       STATE == "SD" |
       STATE == "ND") ~ "2",
    TRUE ~ "3"
  )
)

# Descriptive Statistics
# Key non-zipcode variables
summary(tecos$AGE)
a_label <- gTree("A", children = gList(textGrob("A", x = 0, y = 0.8,
                                                just = c("left", "top"),
                                                gp = gpar(fontsize = 32, col =  "black"))))
tecos_age <- tecos %>%
  ggplot() + 
  geom_histogram(aes(AGE), binwidth = 1, color = "black", fill = "white") + 
  geom_vline(xintercept = 75.96, color = "red") +
  scale_y_continuous(labels = label_comma()) +
  coord_cartesian(xlim = c(65, 110)) +
  labs(title = "Age Distribution\nSitagliptin", x = "Age (years)", y = "Count") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))
tecos_age_grob <- gtable_add_grob(ggplot_gtable(ggplot_build(tecos_age)), a_label, t = 1, l = 4, b = 6)
ggsave("tecos_age_dist.pdf", tecos_age_grob)

table(tecos$RACE)
b_label <- gTree("B", children = gList(textGrob("B", x = 0, y = 0.8,
                                                just = c("left", "top"),
                                                gp = gpar(fontsize = 32, col =  "black"))))
tecos_race <- tecos %>%
  ggplot() + 
  geom_bar(aes(RACE), color = "black", fill = "white") +
  scale_x_continuous(breaks = 0:5,
                     labels = c("White","Black","Asian","Hispanic","Native\nAmerican","Other")) +
  scale_y_continuous(labels = label_comma()) +
  labs(title = "Race Distribution\nSitagliptin", x = "Race", y = "Count") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(hjust = 0.5))
tecos_race_grob <- gtable_add_grob(ggplot_gtable(ggplot_build(tecos_race)), b_label, t = 1, l = 4, b = 6)
ggsave("tecos_race_dist.pdf", tecos_race_grob)

table(tecos$REGION_ZCTA)
c_label <- gTree("C", children = gList(textGrob("C", x = 0, y = 0.8,
                                                just = c("left", "top"),
                                                gp = gpar(fontsize = 32, col =  "black"))))
tecos_region <- tecos %>%
  ggplot() + 
  geom_bar(aes(as.numeric(REGION_ZCTA)), color = "black", fill = "white") +
  scale_x_continuous(breaks = 0:3,
                     labels = c("Northeast",
                                "South",
                                "Midwest",
                                "West")) +
  scale_y_continuous(labels = label_comma()) +
  labs(title = "Region Distribution\nSitagliptin", x = "Region", y = "Count") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(hjust = 0.5))
tecos_region_grob <- gtable_add_grob(ggplot_gtable(ggplot_build(tecos_region)), c_label, t = 1, l = 4, b = 6)
ggsave("tecos_region_dist.pdf", tecos_region_grob)

# Histogram of individuals per zip code
tecos_zip_frequency <- tecos %>% group_by(ZCTA_USE) %>% dplyr::summarize(count=n())

summary(tecos_zip_frequency$count)
table(tecos_zip_frequency$count==1)
table(tecos_zip_frequency$count==2)
table(tecos_zip_frequency$count>=5)
table(tecos_zip_frequency$count>=10)
table(tecos_zip_frequency$count>=15)

ggplot(tecos_zip_frequency) + 
  geom_histogram(aes(count), binwidth = 10, color = "black", fill = "white") + 
  geom_vline(xintercept = 15.19, color = "red") +
  labs(title = "ZCTA Population Distribution (Sitagliptin)", x = "Number of Individuals in ZCTA", y = "Count") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("tecos_zcta_frequency.pdf")

filter(tecos_zip_frequency, count>=5) %>%
  ggplot() + 
  geom_histogram(aes(count), binwidth = 10, color = "black", fill = "white") + 
  geom_vline(xintercept = 23.65, color = "red") +
  labs(title = "ZCTA Population Distribution (Sitagliptin)\n5+ in ZCTA", x ="Number of Individuals in ZCTA", y = "Count") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("tecos_zcta_frequency_5.pdf")

filter(tecos_zip_frequency, count>=10) %>%
  ggplot() + 
  geom_histogram(aes(count), binwidth = 10, color = "black", fill = "white") + 
  geom_vline(xintercept = 30.72, color = "red") +
  labs(title = "ZCTA Population Distribution (Sitagliptin)\n10+ in ZCTA", x = "Number of Individuals in ZCTA", y = "Count") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("tecos_zcta_frequency_10.pdf")

filter(tecos_zip_frequency, count>=15) %>%
  ggplot() + 
  geom_histogram(aes(count), binwidth = 10, color = "black", fill = "white") + 
  geom_vline(xintercept = 36.47, color = "red") +
  labs(title = "ZCTA Population Distribution (Sitagliptin)\n15+ in ZCTA", x = "Number of Individuals in ZCTA", y = "Count") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("tecos_zcta_frequency_15.pdf")

# Add zip code counts to individuals
tecos <- merge(tecos, tecos_zip_frequency, by=c("ZCTA_USE"))

summary(tecos$count)
table(tecos$count==1)
table(tecos$count==2)
table(tecos$count>=5)
table(tecos$count>=10)
table(tecos$count>=15)

ggplot(tecos) + 
  geom_histogram(aes(count), binwidth = 10, color = "black", fill = "white") + 
  geom_vline(xintercept = 43.9, color = "red") +
  labs(title = "Individuals by ZCTA Population Distribution\n(Sitagliptin)", x = "Number of Individuals in ZCTA", y = "Count") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("tecos_individuals_by_zcta_pop.pdf")

# Zip code exposure distribution: 5 or more in zip code
tecos_exp_dist_5 <- tecos %>%
  filter(count>=5) %>%
  group_by(ZCTA_USE) %>%
  dplyr::summarize(prop_exp = mean(EXPOSURE))

summary(tecos_exp_dist_5$prop_exp) # Mean exposure proportion = 0.24, Median = 0.22

# Re-add zip code counts of individuals
tecos_exp_dist_5 <- merge(tecos[,c(1,ncol(tecos))], tecos_exp_dist_5, by = 'ZCTA_USE') %>%
  distinct()

tecos_exp_dist_5_plot <- tecos_exp_dist_5 %>%
  ggplot() + 
  geom_histogram(aes(prop_exp), binwidth = 0.05, color = "black", fill = "white") + 
  geom_vline(xintercept = 0.24, color = "red") +
  scale_y_continuous(labels = label_comma()) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(title = "Sitagliptin Distribution\n5+ in ZCTA", x = "Proportion Sitagliptin", y = "Count of ZCTAs") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(hjust = 0.5))
tecos_exp_dist_5_grob <- gtable_add_grob(ggplot_gtable(ggplot_build(tecos_exp_dist_5_plot)), a_label, t = 1, l = 4, b = 6)
ggsave("tecos_zcta_exp_dist_5.pdf", tecos_exp_dist_5_grob)

# Zip code exposure distribution: 10 or more in zip code
tecos_exp_dist_10 <- tecos %>%
  filter(count>=10) %>%
  group_by(ZCTA_USE) %>%
  dplyr::summarize(prop_exp = mean(EXPOSURE))

summary(tecos_exp_dist_10$prop_exp) # Mean exposure proportion = 0.25, Median = 0.23

tecos_exp_dist_10 %>%
  ggplot() + 
  geom_histogram(aes(prop_exp), binwidth = 0.05, color = "black", fill = "white") + 
  geom_vline(xintercept = 0.25, color = "red") +
  coord_cartesian(xlim = c(0, 1)) +
  labs(title = "Sitagliptin Distribution\n10+ in ZCTA", x = "Proportion Sitagliptin", y = "Count") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("tecos_zcta_exp_dist_10.pdf")

# Zip code exposure distribution: 15 or more in zip code
tecos_exp_dist_15 <- tecos %>%
  filter(count>=15) %>%
  group_by(ZCTA_USE) %>%
  dplyr::summarize(prop_exp = mean(EXPOSURE))

summary(tecos_exp_dist_15$prop_exp) # Mean exposure proportion = 0.25, Median = 0.24

tecos_exp_dist_15 %>%
  ggplot() + 
  geom_histogram(aes(prop_exp), binwidth = 0.05, color = "black", fill = "white") + 
  geom_vline(xintercept = 0.25, color = "red") +
  coord_cartesian(xlim = c(0, 1)) +
  labs(title = "Sitagliptin Distribution\n15+ in ZCTA", x = "Proportion Sitagliptin", y = "Count") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("tecos_zcta_exp_dist_15.pdf")

# Calculate zipcode sitagliptin prescribing proportion for all zipcodes
tecos_exp_dist <- tecos %>%
  group_by(ZCTA_USE) %>%
  dplyr::summarize(prop_exp = mean(EXPOSURE))

# Add zipcode sitagliptin prescribing proportion as new variable for all individuals
tecos <- tecos %>% merge(tecos_exp_dist, by = "ZCTA_USE")

# Descriptive statistics of zipcode sitagliptin exposure proportion
summary(tecos$prop_exp) # Mean exposure proportion = 0.25, Median = 0.24

tecos_exp_dist_5_plot <- filter(tecos, count > 4) %>%
  ggplot() + 
  geom_histogram(aes(prop_exp), binwidth = 0.05, color = "black", fill = "white") + 
  geom_vline(xintercept = 0.26, color = "red") +
  scale_y_continuous(labels = label_comma()) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(title = "Distribution of Individuals by ZCTA Sitagliptin Dispensing Proportion\n5+ in ZCTA", x = "Proportion Sitagliptin in ZCTA", y = "Count of Individuals") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(hjust = 0.5))
tecos_exp_dist_5_grob <- gtable_add_grob(ggplot_gtable(ggplot_build(tecos_exp_dist_5_plot)), a_label, t = 1, l = 4, b = 6)
ggsave("tecos_exp_dist_5.pdf", tecos_exp_dist_5_grob)


# Table 1 of characteristics for all individuals in cohort
tecos_factor_vars <- c("AGE_CAT",
                       "SEX",
                       "RACE",
                       "REGION",
                       "OBESITY",
                       "OVERWEIGHT",
                       "SMOKE",
                       "ALCOHOL",
                       "DRUG",
                       "DM_RETINOPATHY",
                       "DM_OPHTHALMIC",
                       "RETINAL_DETACH",
                       "RETINAL_LASER",
                       "DM_NEUROPATHY2",
                       "DM_NEUROPATHY3",
                       "HYPOGLYCEMIA",
                       "HYPERGLYCEMIA",
                       "ELECTROLYTE_ACID_BASE",
                       "DM_KETOACIDOSIS",
                       "HONK",
                       "DM_PERIPHERAL",
                       "DM_FOOT",
                       "GANGRENE",
                       "LOWER_AMPUTATION",
                       "OSTEOMYELITIS",
                       "SKIN_INFECTION",
                       "ED",
                       "DM_COMPLICATION",
                       "DM_NO_COMPLICATION",
                       "HYPERTENSION",
                       "HYPERLIPIDEMIA",
                       "IHD",
                       "ACUTE_MI",
                       "UNSTABLE_ANGINA",
                       "OLD_MI",
                       "STABLE_ANGINA",
                       "CORONARY_ATHEROSCLEROSIS",
                       "OTHER_ATHEROSCLEROSIS",
                       "CARDIAC_PROCEDURE",
                       "CABG_PTCA",
                       "STROKE",
                       "ISCHEMIC_STROKE",
                       "HEMORRHAGIC_STROKE",
                       "TIA",
                       "OTHER_CEREBROVASCULAR",
                       "LATE_CEREBROVASCULAR",
                       "CEREBROVASCULAR_PROCEDURE",
                       "CHF",
                       "PVD",
                       "ATRIAL_FIBRILLATION",
                       "CARDIAC_DYSRHYTHMIA",
                       "CARDIAC_CONDUCTION",
                       "OTHER_CVD",
                       "EDEMA",
                       "COPD",
                       "ASTHMA",
                       "OBS_SLEEP_APNEA",
                       "PNEUMONIA",
                       "RENAL_DYSFUNCTION",
                       "ACUTE_RENAL_DISEASE",
                       "CHRONIC_RENAL_INSUFFICIENCY",
                       "CKD",
                       "CKD_3_4",
                       "HYPERTENSIVE_NEPHROPATHY",
                       "OTHER_RENAL_INSUFFICIENCY",
                       "LIVER_DISEASE",
                       "OSTEOARTHRITIS",
                       "OTHER_ARTHRITIS",
                       "DORSOPATHIES",
                       "FRACTURE",
                       "FALL",
                       "OSTEOPOROSIS",
                       "HYPERTHYROIDISM",
                       "HYPOTHYROIDISM",
                       "OTHER_THYROID",
                       "DEPRESSION",
                       "ANXIETY",
                       "SLEEP_DISORDER",
                       "DEMENTIA",
                       "DELIRIUM",
                       "PSYCHOSIS",
                       "FRAILTY_QUALITATIVE",
                       "FRAILTY_EMPIRICAL_V3",
                       "NON_FRAILTY",
                       "NAIVE_NEW_USER",
                       "ACE_INHIBITORS",
                       "ARB",
                       "LOOP_DIURETICS",
                       "OTHER_DIURETICS",
                       "NITRATES",
                       "OTHER_ANTIHYPERTENSIVE",
                       "DIGOXIN",
                       "ANTI_ARRHYTHMICS",
                       "COPD_ASTHMA_DRUG",
                       "STATIN",
                       "OTHER_LIPID_DRUG",
                       "ANTIPLATELET",
                       "ORAL_ANTICOAGULANT",
                       "HEPARIN",
                       "NSAID",
                       "ORAL_CORTICOSTEROID",
                       "BISPHOSPHONATE",
                       "OPIOID",
                       "ANTIDEPRESSANT",
                       "ANTIPSYCHOTIC",
                       "ANTICONVULSANT",
                       "LITHIUM",
                       "BENZODIAZAPINE",
                       "ANXIOLYTIC_HYPNOTIC",
                       "DEMENTIA_DRUG",
                       "PARKINSON_DRUG",
                       "HOSPITALIZATION",
                       "ENDOCRINOLOGIST",
                       "INTERNAL",
                       "CARDIOLOGIST",
                       "ELECTROCARDIOGRAM",
                       "GLUCOSE_TESTS",
                       "HOSPITALIZATION_30",
                       "HOSPITALIZATION_31_180",
                       "DM_DRUG_AGI",
                       "DM_DRUG_GLITAZONE",
                       "DM_DRUG_GLP1",
                       "DM_DRUG_INSULIN",
                       "DM_DRUG_MEGLITINIDE",
                       "DM_DRUG_METFORMIN",
                       "PRAMLINTIDE",
                       "GEN1_SU",
                       "MONOTHERAPY_INITIATION",
                       "METFORMIN_DUAL",
                       "CEREBROVASCULAR_HEM_STROKE",
                       "BLADDER_STONE",
                       "KIDNEY_STONE",
                       "UTI",
                       "DIPSTICK_URINALYSIS",
                       "NON_DIPSTICK_URINALYSIS",
                       "URINE_FUNCTION_TEST",
                       "CYTOLOGY",
                       "CYSTOSCOPY",
                       "CREATININE_TEST",
                       "BUN_TEST",
                       "CRI_NO_CKD",
                       "CKD_1_2",
                       "CKD_3_6",
                       "CONCOMITANT_SGLT2I",
                       "CONCOMITANT_AGI",
                       "CONCOMITANT_GLITAZONE",
                       "CONCOMITANT_GLP1",
                       "CONCOMITANT_INSULIN",
                       "CONCOMITANT_MEGLITINIDE",
                       "CONCOMITANT_METFORMIN",
                       "PAST_SGLT2I",
                       "PAST_AGI",
                       "PAST_GLITAZONE",
                       "PAST_GLP1",
                       "PAST_INSULIN",
                       "PAST_MEGLITINIDE",
                       "PAST_METFORMIN",
                       "BLADDER_KIDNEY_STONE",
                       "PERIPHERAL_GANGRENE_OSTEOMYELITIS",
                       "AGE_DECILE",
                       "ALCOHOL_DRUG",
                       "DM_EYE",
                       "COMPOSITE_CVD",
                       "COMPOSITE_CARDIAC_PROCEDURE",
                       "THYROID",
                       "DELIRIUM_PSYCHOSIS",
                       "MEGLITINIDE",
                       "AGI",
                       "GLAUCOMA_CATARACTS",
                       "IMAGING",
                       "CELLULITIS_ABCESS_TOE",
                       "FOOT_ULCER",
                       "ENTRESTO",
                       "ENDOCRINOLOGIST_30",
                       "ENDOCRINOLOGIST_31_180",
                       "CARDIOLOGIST_30",
                       "CARDIOLOGIST_31_180",
                       "INTERNAL_30",
                       "INTERNAL_31_180",
                       "DIALYSIS",
                       "CCI_180",
                       "CKD_3_6_DIALYSIS",
                       "BASELINE_CV",
                       "THIAZIDE",
                       "BETA_BLOCKER",
                       "CA_CHANNEL_BLOCKER")

tecos_vars <- c("AGE_CAT",
                "SEX",
                "RACE",
                "REGION",
                "OBESITY",
                "OVERWEIGHT",
                "SMOKE",
                "ALCOHOL",
                "DRUG",
                "DM_RETINOPATHY",
                "DM_OPHTHALMIC",
                "RETINAL_DETACH",
                "RETINAL_LASER",
                "DM_NEUROPATHY2",
                "DM_NEUROPATHY3",
                "HYPOGLYCEMIA",
                "HYPERGLYCEMIA",
                "ELECTROLYTE_ACID_BASE",
                "DM_KETOACIDOSIS",
                "HONK",
                "DM_PERIPHERAL",
                "DM_FOOT",
                "GANGRENE",
                "LOWER_AMPUTATION",
                "OSTEOMYELITIS",
                "SKIN_INFECTION",
                "ED",
                "DM_COMPLICATION",
                "DM_NO_COMPLICATION",
                "HYPERTENSION",
                "HYPERLIPIDEMIA",
                "IHD",
                "ACUTE_MI",
                "UNSTABLE_ANGINA",
                "OLD_MI",
                "STABLE_ANGINA",
                "CORONARY_ATHEROSCLEROSIS",
                "OTHER_ATHEROSCLEROSIS",
                "CARDIAC_PROCEDURE",
                "CABG_PTCA",
                "STROKE",
                "ISCHEMIC_STROKE",
                "HEMORRHAGIC_STROKE",
                "TIA",
                "OTHER_CEREBROVASCULAR",
                "LATE_CEREBROVASCULAR",
                "CEREBROVASCULAR_PROCEDURE",
                "CHF",
                "PVD",
                "ATRIAL_FIBRILLATION",
                "CARDIAC_DYSRHYTHMIA",
                "CARDIAC_CONDUCTION",
                "OTHER_CVD",
                "EDEMA",
                "COPD",
                "ASTHMA",
                "OBS_SLEEP_APNEA",
                "PNEUMONIA",
                "RENAL_DYSFUNCTION",
                "ACUTE_RENAL_DISEASE",
                "CHRONIC_RENAL_INSUFFICIENCY",
                "CKD",
                "CKD_3_4",
                "HYPERTENSIVE_NEPHROPATHY",
                "OTHER_RENAL_INSUFFICIENCY",
                "LIVER_DISEASE",
                "OSTEOARTHRITIS",
                "OTHER_ARTHRITIS",
                "DORSOPATHIES",
                "FRACTURE",
                "FALL",
                "OSTEOPOROSIS",
                "HYPERTHYROIDISM",
                "HYPOTHYROIDISM",
                "OTHER_THYROID",
                "DEPRESSION",
                "ANXIETY",
                "SLEEP_DISORDER",
                "DEMENTIA",
                "DELIRIUM",
                "PSYCHOSIS",
                "FRAILTY_QUALITATIVE",
                "FRAILTY_EMPIRICAL_V3",
                "NON_FRAILTY",
                "N_ANTIDIABETICS",
                "NAIVE_NEW_USER",
                "ACE_INHIBITORS",
                "ARB",
                "LOOP_DIURETICS",
                "OTHER_DIURETICS",
                "NITRATES",
                "OTHER_ANTIHYPERTENSIVE",
                "DIGOXIN",
                "ANTI_ARRHYTHMICS",
                "COPD_ASTHMA_DRUG",
                "STATIN",
                "OTHER_LIPID_DRUG",
                "ANTIPLATELET",
                "ORAL_ANTICOAGULANT",
                "HEPARIN",
                "NSAID",
                "ORAL_CORTICOSTEROID",
                "BISPHOSPHONATE",
                "OPIOID",
                "ANTIDEPRESSANT",
                "ANTIPSYCHOTIC",
                "ANTICONVULSANT",
                "LITHIUM",
                "BENZODIAZAPINE",
                "ANXIOLYTIC_HYPNOTIC",
                "DEMENTIA_DRUG",
                "PARKINSON_DRUG",
                "N_DIAGNOSES",
                "N_DRUG_RX",
                "HOSPITALIZATION",
                "ENDOCRINOLOGIST",
                "INTERNAL",
                "CARDIOLOGIST",
                "ELECTROCARDIOGRAM",
                "GLUCOSE_TESTS",
                "HOSPITALIZATION_30",
                "HOSPITALIZATION_31_180",
                "N_HOSPITALIZATION",
                "N_HOSPITAL_DAYS",
                "N_ED",
                "N_OFFICE",
                "N_ENDOCRINOLOGIST",
                "N_INTERNAL",
                "N_CARDIOLOGIST",
                "N_ELECTROCARDIOGRAM",
                "N_HBA1C",
                "N_GLUCOSE_TEST",
                "N_LIPID_TEST",
                "N_CREATININE_TEST",
                "N_BUN_TEST",
                "N_MICROALBUMINURIA_TEST",
                "DM_DRUG_AGI",
                "DM_DRUG_GLITAZONE",
                "DM_DRUG_GLP1",
                "DM_DRUG_INSULIN",
                "DM_DRUG_MEGLITINIDE",
                "DM_DRUG_METFORMIN",
                "PRAMLINTIDE",
                "GEN1_SU",
                "MONOTHERAPY_INITIATION",
                "METFORMIN_DUAL",
                "AGE",
                "CEREBROVASCULAR_HEM_STROKE",
                "BLADDER_STONE",
                "KIDNEY_STONE",
                "UTI",
                "DIPSTICK_URINALYSIS",
                "NON_DIPSTICK_URINALYSIS",
                "URINE_FUNCTION_TEST",
                "CYTOLOGY",
                "CYSTOSCOPY",
                "FRAILTY_EMPIRICAL",
                "CREATININE_TEST",
                "BUN_TEST",
                "CRI_NO_CKD",
                "CKD_1_2",
                "CKD_3_6",
                "CONCOMITANT_SGLT2I",
                "CONCOMITANT_AGI",
                "CONCOMITANT_GLITAZONE",
                "CONCOMITANT_GLP1",
                "CONCOMITANT_INSULIN",
                "CONCOMITANT_MEGLITINIDE",
                "CONCOMITANT_METFORMIN",
                "FRAILTY_QUALITATIVE_V1",
                "PAST_SGLT2I",
                "PAST_AGI",
                "PAST_GLITAZONE",
                "PAST_GLP1",
                "PAST_INSULIN",
                "PAST_MEGLITINIDE",
                "PAST_METFORMIN",
                "CALENDAR_TIME_DAY",
                "BLADDER_KIDNEY_STONE",
                "PERIPHERAL_GANGRENE_OSTEOMYELITIS",
                "ALCOHOL_DRUG",
                "DM_EYE",
                "COMPOSITE_CVD",
                "COMPOSITE_CARDIAC_PROCEDURE",
                "THYROID",
                "DELIRIUM_PSYCHOSIS",
                "MEGLITINIDE",
                "AGI",
                "GLAUCOMA_CATARACTS",
                "IMAGING",
                "CELLULITIS_ABCESS_TOE",
                "FOOT_ULCER",
                "ENTRESTO",
                "ENDOCRINOLOGIST_30",
                "ENDOCRINOLOGIST_31_180",
                "CARDIOLOGIST_30",
                "CARDIOLOGIST_31_180",
                "INTERNAL_30",
                "INTERNAL_31_180",
                "DIALYSIS",
                "CCI_180",
                "CKD_3_6_DIALYSIS",
                "BASELINE_CV",
                "THIAZIDE",
                "BETA_BLOCKER",
                "CA_CHANNEL_BLOCKER",
                "PROP_WHITE",
                "PROP_BLACK",
                "PROP_AIAN",
                "PROP_ASIAN",
                "PROP_NHP",
                "PROP_OTHER_RACE",
                "PROP_2MORE_RACE",
                "PROP_NONHISP_WHITE",
                "PROP_NONHISP_BLACK",
                "PROP_NONHISP_AIAN",
                "PROP_NONHISP_ASIAN",
                "PROP_HISP",
                "PROP_HISP_WHITE",
                "PROP_HISP_BLACK",
                "PROP_HISP_ASIAN",
                "PROP_ENGLISH_ONLY",
                "PROP_ENGLISH_WELL",
                "PROP_NO_ENGLISH",
                "PROP_PUBLIC_ASST",
                "PROP_VET",
                "PROP_NATIVE",
                "PROP_NEVER_MARRIED",
                "PROP_MARRIED",
                "PROP_DIVORCED",
                "PROP_WIDOWED",
                "PROP_HS",
                "PROP_ASSOC",
                "PROP_BACHELORS",
                "PROP_POVERTY",
                "PROP_DISABILITY",
                "PROP_COGNITIVE",
                "PROP_AMBULATORY",
                "MED_INC",
                "MED_INC_65",
                "GINI")

tecos_tableone_all <- CreateTableOne(vars = tecos_vars, data = tecos, 
                                     factorVars = tecos_factor_vars,
                                     strata = "EXPOSURE", test = FALSE, smd = TRUE)
tecos_tableone_all_export <- print(tecos_tableone_all, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(tecos_tableone_all_export, file = "/Volumes/PACS$/Cordes, Jack/tecos_tableone_all.csv")

# Table 1 of characteristics for individuals with at least 5 in zip code by exposure group
tecos_5 <- tecos %>%
  filter(count>=5)

tecos_tableone_5_exp <- CreateTableOne(vars = tecos_vars, data = tecos_5, 
                                       factorVars = tecos_factor_vars,
                                       strata = "EXPOSURE", test = FALSE, smd = TRUE)
tecos_tableone_5_exp_export <- print(tecos_tableone_5_exp, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(tecos_tableone_5_exp_export, file = "/Volumes/PACS$/Cordes, Jack/tecos_tableone_5_exp.csv")

# Table 1 of characteristics for individuals with at least 5 in zip code vs. fewer
tecos$AT_LEAST_5 <- ifelse(tecos$count > 4,1,0)

tecos_tableone_5_vs_less <- CreateTableOne(vars = c(tecos_vars, "EXPOSURE"), data = tecos, 
                                           factorVars = c(tecos_factor_vars, "EXPOSURE"),
                                           strata = "AT_LEAST_5", test = FALSE, smd = TRUE)
tecos_tableone_5_vs_less_export <- print(tecos_tableone_5_vs_less, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(tecos_tableone_5_vs_less_export, file = "/Volumes/PACS$/Cordes, Jack/tecos_tableone_5_vs_less.csv")

# Table 1 of characteristics for all individuals in cohort limited to PS variables
tecos_factor_vars_ps <- c("AGE_CAT",
                          "SEX",
                          "RACE",
                          "REGION",
                          "OBESITY",
                          "OVERWEIGHT",
                          "SMOKE",
                          "ALCOHOL",
                          "DRUG",
                          "DM_NEUROPATHY3",
                          "HYPOGLYCEMIA",
                          "HYPERGLYCEMIA",
                          "DM_KETOACIDOSIS",
                          "HYPERTENSION",
                          "HYPERLIPIDEMIA",
                          "UNSTABLE_ANGINA",
                          "STABLE_ANGINA",
                          "CORONARY_ATHEROSCLEROSIS",
                          "TIA",
                          "CHF",
                          "PVD",
                          "CARDIAC_DYSRHYTHMIA",
                          "EDEMA",
                          "COPD",
                          "ASTHMA",
                          "CKD",
                          "LIVER_DISEASE",
                          "DEPRESSION",
                          "ANXIETY",
                          "FRAILTY_EMPIRICAL_V3",
                          "ACE_INHIBITORS",
                          "ARB",
                          "STATIN",
                          "ANTIPLATELET",
                          "NSAID",
                          "ORAL_CORTICOSTEROID",
                          "OPIOID",
                          "BENZODIAZAPINE",
                          "HOSPITALIZATION",
                          "ENDOCRINOLOGIST",
                          "CARDIOLOGIST",
                          "DM_DRUG_AGI",
                          "DM_DRUG_GLITAZONE",
                          "DM_DRUG_GLP1",
                          "DM_DRUG_INSULIN",
                          "DM_DRUG_MEGLITINIDE",
                          "DM_DRUG_METFORMIN",
                          "PRAMLINTIDE",
                          "GEN1_SU",
                          "MONOTHERAPY_INITIATION",
                          "METFORMIN_DUAL",
                          "CEREBROVASCULAR_HEM_STROKE",
                          "URINE_FUNCTION_TEST",
                          "CRI_NO_CKD",
                          "CONCOMITANT_SGLT2I",
                          "CONCOMITANT_AGI",
                          "CONCOMITANT_GLITAZONE",
                          "CONCOMITANT_GLP1",
                          "CONCOMITANT_INSULIN",
                          "CONCOMITANT_MEGLITINIDE",
                          "CONCOMITANT_METFORMIN",
                          "PAST_SGLT2I",
                          "PAST_AGI",
                          "PAST_GLITAZONE",
                          "PAST_GLP1",
                          "PAST_INSULIN",
                          "PAST_MEGLITINIDE",
                          "PAST_METFORMIN",
                          "DM_EYE",
                          "COMPOSITE_CVD",
                          "COMPOSITE_CARDIAC_PROCEDURE",
                          "THYROID",
                          "DIALYSIS",
                          "BETA_BLOCKER",
                          "CA_CHANNEL_BLOCKER")

tecos_vars_ps <- c("AGE_CAT",
                   "SEX",
                   "RACE",
                   "REGION",
                   "OBESITY",
                   "OVERWEIGHT",
                   "SMOKE",
                   "ALCOHOL",
                   "DRUG",
                   "DM_NEUROPATHY3",
                   "HYPOGLYCEMIA",
                   "HYPERGLYCEMIA",
                   "DM_KETOACIDOSIS",
                   "HYPERTENSION",
                   "HYPERLIPIDEMIA",
                   "UNSTABLE_ANGINA",
                   "STABLE_ANGINA",
                   "CORONARY_ATHEROSCLEROSIS",
                   "TIA",
                   "CHF",
                   "PVD",
                   "CARDIAC_DYSRHYTHMIA",
                   "EDEMA",
                   "COPD",
                   "ASTHMA",
                   "CKD",
                   "LIVER_DISEASE",
                   "DEPRESSION",
                   "ANXIETY",
                   "FRAILTY_EMPIRICAL_V3",
                   "N_ANTIDIABETICS",
                   "ACE_INHIBITORS",
                   "ARB",
                   "STATIN",
                   "ANTIPLATELET",
                   "NSAID",
                   "ORAL_CORTICOSTEROID",
                   "OPIOID",
                   "BENZODIAZAPINE",
                   "N_DIAGNOSES",
                   "N_DRUG_RX",
                   "HOSPITALIZATION",
                   "ENDOCRINOLOGIST",
                   "CARDIOLOGIST",
                   "N_OFFICE",
                   "N_HBA1C",
                   "N_GLUCOSE_TEST",
                   "DM_DRUG_AGI",
                   "DM_DRUG_GLITAZONE",
                   "DM_DRUG_GLP1",
                   "DM_DRUG_INSULIN",
                   "DM_DRUG_MEGLITINIDE",
                   "DM_DRUG_METFORMIN",
                   "PRAMLINTIDE",
                   "GEN1_SU",
                   "MONOTHERAPY_INITIATION",
                   "METFORMIN_DUAL",
                   "CEREBROVASCULAR_HEM_STROKE",
                   "URINE_FUNCTION_TEST",
                   "CRI_NO_CKD",
                   "CONCOMITANT_SGLT2I",
                   "CONCOMITANT_AGI",
                   "CONCOMITANT_GLITAZONE",
                   "CONCOMITANT_GLP1",
                   "CONCOMITANT_INSULIN",
                   "CONCOMITANT_MEGLITINIDE",
                   "CONCOMITANT_METFORMIN",
                   "PAST_SGLT2I",
                   "PAST_AGI",
                   "PAST_GLITAZONE",
                   "PAST_GLP1",
                   "PAST_INSULIN",
                   "PAST_MEGLITINIDE",
                   "PAST_METFORMIN",
                   "CALENDAR_TIME_DAY",
                   "DM_EYE",
                   "COMPOSITE_CVD",
                   "COMPOSITE_CARDIAC_PROCEDURE",
                   "THYROID",
                   "DIALYSIS",
                   "CCI_180",
                   "BETA_BLOCKER",
                   "CA_CHANNEL_BLOCKER",
                   "PROP_WHITE",
                   "PROP_BLACK",
                   "PROP_AIAN",
                   "PROP_ASIAN",
                   "PROP_OTHER_RACE",
                   "PROP_HISP",
                   "PROP_NO_ENGLISH",
                   "PROP_PUBLIC_ASST",
                   "PROP_NATIVE",
                   "PROP_MARRIED",
                   "PROP_BACHELORS",
                   "PROP_POVERTY",
                   "PROP_DISABILITY",
                   "PROP_COGNITIVE",
                   "PROP_AMBULATORY",
                   "MED_INC",
                   "MED_INC_65",
                   "GINI")

# Vector of zipcode level variables
tecos_vars_zip <- c("PROP_WHITE",
                    "PROP_BLACK",
                    "PROP_AIAN",
                    "PROP_ASIAN",
                    "PROP_OTHER_RACE",
                    "PROP_HISP",
                    "PROP_NO_ENGLISH",
                    "PROP_PUBLIC_ASST",
                    "PROP_NATIVE",
                    "PROP_MARRIED",
                    "PROP_BACHELORS",
                    "PROP_POVERTY",
                    "PROP_DISABILITY",
                    "PROP_COGNITIVE",
                    "PROP_AMBULATORY",
                    "MED_INC",
                    "MED_INC_65",
                    "GINI")

tecos_tableone_all_ps <- CreateTableOne(vars = tecos_vars_ps, data = tecos, 
                                        factorVars = tecos_factor_vars_ps,
                                        strata = "EXPOSURE", test = FALSE, smd = TRUE)
tecos_tableone_all_ps_export <- print(tecos_tableone_all_ps, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(tecos_tableone_all_ps_export, file = "/Volumes/PACS$/Cordes, Jack/tecos_tableone_all_ps_14_17.csv")

# Table 1 of characteristics for individuals with at least 5 in zip code by exposure group
tecos_5 <- tecos %>%
  filter(count>=5)

tecos_tableone_5_exp_ps <- CreateTableOne(vars = tecos_vars_ps, data = tecos_5, 
                                          factorVars = tecos_factor_vars_ps,
                                          strata = "EXPOSURE", test = FALSE, smd = TRUE)
tecos_tableone_5_exp_ps_export <- print(tecos_tableone_5_exp_ps, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(tecos_tableone_5_exp_ps_export, file = "/Volumes/PACS$/Cordes, Jack/tecos_tableone_5_exp_ps.csv")


# Table 1 of characteristics for individuals with at least 5 in zip code vs. fewer, ps variables only
tecos$AT_LEAST_5 <- ifelse(tecos$count > 4,1,0)

tecos_tableone_5_vs_less_ps <- CreateTableOne(vars = c(tecos_vars_ps, "EXPOSURE"), data = tecos, 
                                              factorVars = c(tecos_factor_vars_ps, "EXPOSURE"),
                                              strata = "AT_LEAST_5", test = FALSE, smd = TRUE)
tecos_tableone_5_vs_less_ps_export <- print(tecos_tableone_5_vs_less_ps, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(tecos_tableone_5_vs_less_ps_export, file = "/Volumes/PACS$/Cordes, Jack/tecos_tableone_5_vs_less_ps.csv")

# Mahalanobis distance of PS variables by treatment group for all individuals
tecos_mah_dist_total <- tecos[,c(1:13,19:21,23,34:35,38,40:41,48,52:53,55,58:60,66,70,80:81,87,89,91:92,100,102,105:106,108,113,117:120,122,130,135:136,141:150,152,158,164,167:173,175:182,187:190,205:206,210:211,220,226:229,231,237,244:245,247,249,254:258,262:264)]

tecos_mah_dist_total_treatment <- mahalanobis(colMeans(tecos_mah_dist_total %>% filter(EXPOSURE == 1) %>% select(c(5:87,89:106))),
                                              colMeans(tecos_mah_dist_total %>% filter(EXPOSURE == 0) %>% select(c(5:87,89:106))),
                                              cov(tecos_mah_dist_total %>% select(c(5:87,89:106))),
                                              tol = 1e-30)
print(tecos_mah_dist_total_treatment)

# Mahalanobis distance of PS variables by treatment group at least 5 in ZCTA
tecos_mah_dist_total_5 <- filter(tecos[,c(1:13,19:21,23,34:35,38,40:41,48,52:53,55,58:60,66,70,80:81,87,89,91:92,100,102,105:106,108,113,117:120,122,130,135:136,141:150,152,158,164,167:173,175:182,187:190,205:206,210:211,220,226:229,231,237,244:245,247,249,254:258,262:264,268)], count > 4)

tecos_mah_dist_total_treatment_5 <- mahalanobis(colMeans(tecos_mah_dist_total_5 %>% filter(EXPOSURE == 1) %>% select(c(5:87,89:106))),
                                                colMeans(tecos_mah_dist_total_5 %>% filter(EXPOSURE == 0) %>% select(c(5:87,89:106))),
                                                cov(tecos_mah_dist_total_5 %>% select(c(5:87,89:106))),
                                                tol = 1e-30)
print(tecos_mah_dist_total_treatment_5)

# Mahalanobis distance of ZCTA variables by treatment group at least 5 in ZCTA
tecos_mah_dist_total_5_zcta <- filter(tecos[,c(1:4,220,226:229,231,237,244:245,247,249,254:258,262:264,268)], count > 4)

tecos_mah_dist_total_treatment_5_zcta <- mahalanobis(colMeans(tecos_mah_dist_total_5_zcta %>% filter(EXPOSURE == 1) %>% select(c(6:23))),
                                                     colMeans(tecos_mah_dist_total_5_zcta %>% filter(EXPOSURE == 0) %>% select(c(6:23))),
                                                     cov(tecos_mah_dist_total_5_zcta %>% select(c(6:23))),
                                                     tol = 1e-30)
print(tecos_mah_dist_total_treatment_5_zcta)


# Risk factor analysis
# Among Sitagliptin
summary(
  glm(data = filter(tecos, EXPOSURE == 1), formula = OUTCOME ~ PROP_WHITE,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 1), formula = OUTCOME ~ PROP_BLACK,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 1), formula = OUTCOME ~ PROP_AIAN,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 1), formula = OUTCOME ~ PROP_ASIAN,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 1), formula = OUTCOME ~ PROP_OTHER_RACE,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 1), formula = OUTCOME ~ PROP_HISP,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 1), formula = OUTCOME ~ PROP_NO_ENGLISH,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 1), formula = OUTCOME ~ PROP_PUBLIC_ASST,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 1), formula = OUTCOME ~ PROP_NATIVE,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 1), formula = OUTCOME ~ PROP_MARRIED,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 1), formula = OUTCOME ~ PROP_BACHELORS,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 1), formula = OUTCOME ~ PROP_POVERTY,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 1), formula = OUTCOME ~ PROP_DISABILITY,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 1), formula = OUTCOME ~ PROP_COGNITIVE,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 1), formula = OUTCOME ~ PROP_AMBULATORY,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 1), formula = OUTCOME ~ MED_INC_STD,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 1), formula = OUTCOME ~ MED_INC_65_STD,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 1), formula = OUTCOME ~ GINI,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)

# Among SU
summary(
  glm(data = filter(tecos, EXPOSURE == 0), formula = OUTCOME ~ PROP_WHITE,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 0), formula = OUTCOME ~ PROP_BLACK,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 0), formula = OUTCOME ~ PROP_AIAN,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 0), formula = OUTCOME ~ PROP_ASIAN,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 0), formula = OUTCOME ~ PROP_OTHER_RACE,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 0), formula = OUTCOME ~ PROP_HISP,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 0), formula = OUTCOME ~ PROP_NO_ENGLISH,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 0), formula = OUTCOME ~ PROP_PUBLIC_ASST,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 0), formula = OUTCOME ~ PROP_NATIVE,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 0), formula = OUTCOME ~ PROP_MARRIED,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 0), formula = OUTCOME ~ PROP_BACHELORS,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 0), formula = OUTCOME ~ PROP_POVERTY,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 0), formula = OUTCOME ~ PROP_DISABILITY,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 0), formula = OUTCOME ~ PROP_COGNITIVE,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 0), formula = OUTCOME ~ PROP_AMBULATORY,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 0), formula = OUTCOME ~ MED_INC_STD,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 0), formula = OUTCOME ~ MED_INC_65_STD,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(tecos, EXPOSURE == 0), formula = OUTCOME ~ GINI,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)

# Full cohort
summary(
  glm(data = tecos, formula = OUTCOME ~ PROP_WHITE,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = tecos, formula = OUTCOME ~ PROP_BLACK,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = tecos, formula = OUTCOME ~ PROP_AIAN,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = tecos, formula = OUTCOME ~ PROP_ASIAN,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = tecos, formula = OUTCOME ~ PROP_OTHER_RACE,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = tecos, formula = OUTCOME ~ PROP_HISP,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = tecos, formula = OUTCOME ~ PROP_NO_ENGLISH,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = tecos, formula = OUTCOME ~ PROP_PUBLIC_ASST,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = tecos, formula = OUTCOME ~ PROP_NATIVE,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = tecos, formula = OUTCOME ~ PROP_MARRIED,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = tecos, formula = OUTCOME ~ PROP_BACHELORS,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = tecos, formula = OUTCOME ~ PROP_POVERTY,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = tecos, formula = OUTCOME ~ PROP_DISABILITY,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = tecos, formula = OUTCOME ~ PROP_COGNITIVE,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = tecos, formula = OUTCOME ~ PROP_AMBULATORY,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = tecos, formula = OUTCOME ~ MED_INC_STD,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = tecos, formula = OUTCOME ~ MED_INC_65_STD,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = tecos, formula = OUTCOME ~ GINI,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)

# Unadjusted Association
# Distribution of follow-up times
summary(tecos$FUP_TIME)

tecos %>%
  ggplot() + 
  geom_histogram(aes(FUP_TIME), binwidth = 30, color = "black", fill = "white") + 
  geom_vline(xintercept = 247.4, color = "red") +
  scale_y_continuous(labels = label_comma()) +
  coord_cartesian(xlim = c(0, 2200)) +
  labs(title = "Follow-up Time Distribution\nSitagliptin", x = "Days", y = "Count") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))
ggsave(filename = "tecos_followup_dist.pdf")

# Raw outcome counts by exposure status
table(tecos$EXPOSURE, tecos$OUTCOME)
tecos %>% group_by(EXPOSURE) %>% summarise(OUTCOME_SUM = sum(OUTCOME))

# Person-years by exposure status
tecos %>% group_by(EXPOSURE) %>% summarise(FUP_TIME_SUM = sum(FUP_TIME)/365.25)

# Risk in the exposed and unexposed
tecos_exp_risk <- 5243/(5243 + 88947)
tecos_unexp_risk <- 19434/(19434 + 255053)

# Unadjusted risk ratio and odds ratio
tecos_rr <- tecos_exp_risk/tecos_unexp_risk

tecos_or <- (284588*5785)/(21714*97892)

# Unadjusted ITT cumulative incidence analysis using logistic regression
tecos_model_crude_itt <- glm(data = tecos, formula = OUTCOME ~ EXPOSURE,
                             family = binomial(link = "logit"),
                             na.action = na.exclude)
summary(tecos_model_crude_itt)
exp(tecos_model_crude_itt$coefficients)

# Unadjusted survival analysis
tecos_model_crude_survival <- coxph(data = tecos, formula = Surv(FUP_TIME, OUTCOME) ~ EXPOSURE,
                                    na.action = na.exclude)
summary(tecos_model_crude_survival)
cox.zph(tecos_model_crude_survival)
ggcoxzph(cox.zph(tecos_model_crude_survival))
ggsave(filename = "tecos_model_crude_survival.pdf")

# Kaplan-Meier curves
tecos_km <- survfit(data = tecos, formula = Surv(FUP_TIME, OUTCOME) ~ EXPOSURE)

center_title <- function() {
  theme_survminer() %+replace%
    theme(plot.title = element_text(size = 18, hjust = 0),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18))
}

tecos_km_plot <- ggsurvplot(fit = tecos_km,
                            palette = c("Red", "Blue"),
                            censor = FALSE,
                            title = "A                                      Survival Curves\n                                      Sitagliptin Cohort",
                            xlab = "Days",
                            ylab = "Survival Probability",
                            legend.title = "",
                            legend.labs = c("SU", "Sitagliptin"),
                            legend = "right",
                            conf.int = T,
                            ggtheme = center_title())
ggsave("tecos_km.pdf", tecos_km_plot$plot)

# Generate propensity scores based on 2014-2017 cohort
# Subset data to 2014-2017
tecos_phys <- filter(tecos, YEAR >= 2014)

# Raw outcome counts by exposure status
table(tecos_phys$EXPOSURE, tecos_phys$OUTCOME)
tecos_phys %>% group_by(EXPOSURE) %>% summarise(OUTCOME_SUM = sum(OUTCOME))

# Unadjusted risk ratio analysis linear regression 6-month follow-up
tecos_model_crude_linear_6mo <- glm(data = tecos_phys, formula = OUTCOME_6MO ~ EXPOSURE,
                                    na.action = na.exclude)
summary(tecos_model_crude_linear_6mo)

# Unadjusted risk ratio analysis linear regression 1-year follow-up
tecos_model_crude_linear_1yr <- glm(data = tecos_phys, formula = OUTCOME_1YR ~ EXPOSURE,
                                    na.action = na.exclude)
summary(tecos_model_crude_linear_1yr)

# Unadjusted risk ratio analysis linear regression all follow-up
tecos_model_crude_linear <- glm(data = tecos_phys, formula = OUTCOME ~ EXPOSURE,
                                na.action = na.exclude)
summary(tecos_model_crude_linear)

# Analysis adjusted for zipcode variables
tecos_ps_model_zip_14_17 <- glm(data = tecos_phys, formula = EXPOSURE ~ 
                                  as.factor(AGE_CAT) +
                                  as.factor(SEX) +
                                  as.factor(RACE) +
                                  as.factor(REGION) +
                                  as.factor(OBESITY) +
                                  as.factor(OVERWEIGHT) +
                                  as.factor(SMOKE) +
                                  as.factor(ALCOHOL) +
                                  as.factor(DRUG) +
                                  as.factor(DM_NEUROPATHY3) +
                                  as.factor(HYPOGLYCEMIA) +
                                  as.factor(HYPERGLYCEMIA) +
                                  as.factor(DM_KETOACIDOSIS) +
                                  as.factor(HYPERTENSION) +
                                  as.factor(HYPERLIPIDEMIA) +
                                  as.factor(UNSTABLE_ANGINA) +
                                  as.factor(STABLE_ANGINA) +
                                  as.factor(CORONARY_ATHEROSCLEROSIS) +
                                  as.factor(TIA) +
                                  as.factor(CHF) +
                                  as.factor(PVD) +
                                  as.factor(CARDIAC_DYSRHYTHMIA) +
                                  as.factor(EDEMA) +
                                  as.factor(COPD) +
                                  as.factor(ASTHMA) +
                                  as.factor(CKD) +
                                  as.factor(LIVER_DISEASE) +
                                  as.factor(DEPRESSION) +
                                  as.factor(ANXIETY) +
                                  as.factor(FRAILTY_EMPIRICAL_V3) +
                                  N_ANTIDIABETICS +
                                  as.factor(ACE_INHIBITORS) +
                                  as.factor(ARB) +
                                  as.factor(STATIN) +
                                  as.factor(ANTIPLATELET) +
                                  as.factor(NSAID) +
                                  as.factor(ORAL_CORTICOSTEROID) +
                                  as.factor(OPIOID) +
                                  as.factor(BENZODIAZAPINE) +
                                  N_DIAGNOSES +
                                  N_DRUG_RX +
                                  as.factor(HOSPITALIZATION) +
                                  as.factor(ENDOCRINOLOGIST) +
                                  as.factor(CARDIOLOGIST) +
                                  N_OFFICE +
                                  N_HBA1C +
                                  N_GLUCOSE_TEST +
                                  as.factor(DM_DRUG_AGI) +
                                  as.factor(DM_DRUG_GLITAZONE) +
                                  as.factor(DM_DRUG_GLP1) +
                                  as.factor(DM_DRUG_INSULIN) +
                                  as.factor(DM_DRUG_MEGLITINIDE) +
                                  as.factor(DM_DRUG_METFORMIN) +
                                  as.factor(PRAMLINTIDE) +
                                  as.factor(GEN1_SU) +
                                  as.factor(MONOTHERAPY_INITIATION) +
                                  as.factor(METFORMIN_DUAL) +
                                  as.factor(CEREBROVASCULAR_HEM_STROKE) +
                                  as.factor(URINE_FUNCTION_TEST) +
                                  as.factor(CRI_NO_CKD) +
                                  as.factor(CONCOMITANT_SGLT2I) +
                                  as.factor(CONCOMITANT_AGI) +
                                  as.factor(CONCOMITANT_GLITAZONE) +
                                  as.factor(CONCOMITANT_GLP1) +
                                  as.factor(CONCOMITANT_INSULIN) +
                                  as.factor(CONCOMITANT_MEGLITINIDE) +
                                  as.factor(CONCOMITANT_METFORMIN) +
                                  as.factor(PAST_SGLT2I) +
                                  as.factor(PAST_AGI) +
                                  as.factor(PAST_GLITAZONE) +
                                  as.factor(PAST_GLP1) +
                                  as.factor(PAST_INSULIN) +
                                  as.factor(PAST_MEGLITINIDE) +
                                  as.factor(PAST_METFORMIN) +
                                  CALENDAR_TIME_DAY +
                                  as.factor(DM_EYE) +
                                  as.factor(COMPOSITE_CVD) +
                                  as.factor(COMPOSITE_CARDIAC_PROCEDURE) +
                                  as.factor(THYROID) +
                                  as.factor(DIALYSIS) +
                                  CCI_180 +
                                  as.factor(BETA_BLOCKER) +
                                  as.factor(CA_CHANNEL_BLOCKER) +
                                  PROP_WHITE +
                                  PROP_BLACK +
                                  PROP_AIAN +
                                  PROP_ASIAN +
                                  PROP_OTHER_RACE +
                                  PROP_HISP +
                                  PROP_NO_ENGLISH +
                                  PROP_PUBLIC_ASST +
                                  PROP_NATIVE +
                                  PROP_MARRIED +
                                  PROP_BACHELORS +
                                  PROP_POVERTY +
                                  PROP_DISABILITY +
                                  PROP_COGNITIVE +
                                  PROP_AMBULATORY +
                                  MED_INC_STD +
                                  MED_INC_65_STD +
                                  GINI, family = binomial(link = "logit"),
                                na.action = na.exclude)

tecos_ps_model_zip_14_17_pred <- predict(tecos_ps_model_zip_14_17, type = "response")
roc(tecos_phys$EXPOSURE, tecos_ps_model_zip_14_17_pred)

# Calculate PS
tecos_phys$PS_ZIP_14_17 <- fitted(tecos_ps_model_zip_14_17)
summary(tecos_phys$PS_ZIP_14_17)

# PS quintiles
tecos_phys %>%
  dplyr::summarize(min = quantile(PS_ZIP_14_17, probs = 0),
                   q1 = quantile(PS_ZIP_14_17, probs = 0.2),
                   q2 = quantile(PS_ZIP_14_17, probs = 0.4),
                   q3 = quantile(PS_ZIP_14_17, probs = 0.6),
                   q4 = quantile(PS_ZIP_14_17, probs = 0.8),
                   max = quantile(PS_ZIP_14_17, probs = 1))

# PS distribution by exposure
tecos_ps_dist_zip_14_17_plot <- ggplot(tecos_phys) +
  geom_density(aes(PS_ZIP_14_17, group = factor(EXPOSURE), color = factor(EXPOSURE))) + 
  labs(title = "Sitagliptin Cohort Propensity Scores\nwith ZCTA Variables 2014-2017", x = "Propensity Score", y = "Density") +
  scale_color_manual(name = "", values = c("red", "blue"), labels = c("SU", "Sitagliptin")) +
  scale_x_continuous(limits = c(0, 1)) +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(hjust = 0.5),
        legend.text = element_text(size = 18))
ggsave("tecos_ps_dist_zip_14_17.pdf", tecos_ps_dist_zip_14_17_plot)

# Create PS quintiles variable
tecos_phys <- tecos_phys %>% mutate(
  PS_ZIP_QUINTILE_14_17 = case_when(
    PS_ZIP_14_17 < 0.181776 ~ "1",
    PS_ZIP_14_17 >= 0.181776 & PS_ZIP_14_17 < 0.2261959 ~ "2",
    PS_ZIP_14_17 >= 0.2261959 & PS_ZIP_14_17 < 0.2703549 ~ "3",
    PS_ZIP_14_17 >= 0.2703549 & PS_ZIP_14_17 < 0.3370492 ~ "4",
    TRUE ~ "5"
  )
)

# Global historical zcta preference IV analysis
# Sort on zcta and CED
tecos_phys <- tecos_phys[order(tecos_phys$ZCTA_USE, tecos_phys$ENTRYDATE),]

# Create empty global historical zcta preference column
tecos_phys <- tecos_phys %>% mutate(prop_exp_hist = NA)

# Create empty counter column
tecos_phys <- tecos_phys %>% mutate(counter = NA)

# Create empty cumulative sum column
tecos_phys <- tecos_phys %>% mutate(cum_sum = NA)

# Create global historical zcta preference variable based on average of all previous prescriptions
for (i in 2:nrow(tecos_phys)){
  if ((tecos_phys$ZCTA_USE[i] == tecos_phys$ZCTA_USE[i-1]) & (is.na(tecos_phys$prop_exp_hist[i-1]) == TRUE)){
    tecos_phys$prop_exp_hist[i] <- tecos_phys$EXPOSURE[i-1]
  }
  if ((tecos_phys$ZCTA_USE[i] == tecos_phys$ZCTA_USE[i-1]) & (is.na(tecos_phys$prop_exp_hist[i-1]) == TRUE)){
    tecos_phys$counter[i] <- 1
  }
  if ((tecos_phys$ZCTA_USE[i] == tecos_phys$ZCTA_USE[i-1]) & (is.na(tecos_phys$prop_exp_hist[i-1]) == TRUE)){
    tecos_phys$cum_sum[i] <- tecos_phys$EXPOSURE[i-1]
  }
  if ((tecos_phys$ZCTA_USE[i] == tecos_phys$ZCTA_USE[i-1]) & (is.na(tecos_phys$prop_exp_hist[i-1]) == FALSE)){
    tecos_phys$counter[i] <- tecos_phys$counter[i-1]+1
  }
  if ((tecos_phys$ZCTA_USE[i] == tecos_phys$ZCTA_USE[i-1]) & (is.na(tecos_phys$prop_exp_hist[i-1]) == FALSE)){
    tecos_phys$cum_sum[i] <- tecos_phys$cum_sum[i-1] + tecos_phys$EXPOSURE[i-1]
  }
  if ((tecos_phys$ZCTA_USE[i] == tecos_phys$ZCTA_USE[i-1]) & (is.na(tecos_phys$prop_exp_hist[i-1]) == FALSE)){
    tecos_phys$prop_exp_hist[i] <- tecos_phys$cum_sum[i]/tecos_phys$counter[i]
  }
}

# Distribution of global historical zcta proportions
summary(tecos_phys$prop_exp_hist)

# Create dataset limited to non-missing zcta historical IV
tecos_zip_hist <- filter(tecos_phys, is.na(prop_exp_hist) == FALSE)

tecos_zip_hist_dist_plot <- tecos_zip_hist %>%
  ggplot() + 
  geom_histogram(aes(prop_exp_hist), binwidth = 0.05, color = "black", fill = "white") + 
  geom_vline(xintercept = 0.258, color = "red") +
  scale_y_continuous(labels = label_comma()) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(title = "Sitagliptin Cohort\nCumulative ZCTA Proportion for Each Individual", x = "Proportion Sitagliptin", y = "Count of Individuals") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(hjust = 0.5))
tecos_zip_hist_dist_grob <- gtable_add_grob(ggplot_gtable(ggplot_build(tecos_zip_hist_dist_plot)), a_label, t = 1, l = 4, b = 6)
ggsave("tecos_zip_hist_dist.pdf", tecos_zip_hist_dist_grob)

# For IV with proportion sitagliptin cutoffs 100% vs. 0%, instrument is collapsed with exposure
# because no variation in individual-level exposure
# So becomes simple regression comparing individuals with 100% to 0% prescribers

# ZCTA proportion sitagliptin 100% vs. 0%
tecos_phys <- tecos_phys %>%
  mutate(ZIP_HIST_0_100 = case_when(
    prop_exp_hist == 1 ~ "1",
    prop_exp_hist == 0 ~ "0",
    TRUE ~ NA_character_
  )
  )

# Ordinary least squares - no adjustment for confounders
tecos_model_zip_hist_ols_0_100 <- glm(data = filter(tecos_phys, (prop_exp_hist == 0 | prop_exp_hist == 1)),
                                      formula = OUTCOME ~ EXPOSURE,
                                      na.action = na.exclude)
summary(tecos_model_zip_hist_ols_0_100)

# Ordinary least squares - adjusted for PS
tecos_model_zip_hist_ols_0_100_adj <- glm(data = filter(tecos_phys, (prop_exp_hist == 0 | prop_exp_hist == 1)),
                                          formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                          na.action = na.exclude)
summary(tecos_model_zip_hist_ols_0_100_adj)

# Logistic regression - no adjustment for confounders
tecos_model_zip_hist_logistic_0_100 <- glm(data = filter(tecos_phys, (prop_exp_hist == 0 | prop_exp_hist == 1)),
                                           formula = OUTCOME ~ EXPOSURE,
                                           family = binomial(link = "logit"),
                                           na.action = na.exclude)
summary(tecos_model_zip_hist_logistic_0_100)
exp(tecos_model_zip_hist_logistic_0_100$coefficients)

# Logistic regression - adjusted for PS
tecos_model_zip_hist_logistic_0_100_adj <- glm(data = filter(tecos_phys, (prop_exp_hist == 0 | prop_exp_hist == 1)),
                                               formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                               family = binomial(link = "logit"),
                                               na.action = na.exclude)
summary(tecos_model_zip_hist_logistic_0_100_adj)
exp(tecos_model_zip_hist_logistic_0_100_adj$coefficients)

# Table 1 at the individual level
tecos_tableone_all_iv_zip_hist_0_100 <- CreateTableOne(vars = tecos_vars_ps, data = filter(tecos_phys, is.na(ZIP_HIST_0_100) == FALSE),
                                                       factorVars = tecos_factor_vars_ps,
                                                       strata = "ZIP_HIST_0_100", test = FALSE, smd = TRUE)
tecos_tableone_all_iv_zip_hist_0_100_export <- print(tecos_tableone_all_iv_zip_hist_0_100, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(tecos_tableone_all_iv_zip_hist_0_100_export, file = "/Volumes/PACS$/Cordes, Jack/tecos_tableone_all_iv_zip_hist_0_100.csv")

# Check balance improvement by percent change in Mahalanobis distance using treatment
tecos_zip_hist_iv_0_100_mah_dist <- tecos_phys[,c(1:13,19:21,23,34:35,38,40:41,48,52:53,55,58:60,66,70,80:81,87,89,91:92,100,102,105:106,108,113,117:120,122,130,135:136,141:150,152,158,164,167:173,175:182,187:190,205:206,210:211,220,226:229,231,237,244:245,247,249,254:258,262:264,271,295)]

tecos_mah_dist_total_zip_hist_iv_0_100 <- mahalanobis(colMeans(tecos_zip_hist_iv_0_100_mah_dist %>% filter(ZIP_HIST_0_100 == 1) %>% select(c(5:87,89:106))),
                                                      colMeans(tecos_zip_hist_iv_0_100_mah_dist %>% filter(ZIP_HIST_0_100 == 0) %>% select(c(5:87,89:106))),
                                                      cov(tecos_zip_hist_iv_0_100_mah_dist %>% select(c(5:87,89:106))),
                                                      tol = 1e-30)

((tecos_mah_dist_total_zip_hist_iv_0_100-tecos_mah_dist_total_14_17_treatment)/(tecos_mah_dist_total_14_17_treatment))*100

# Distribution of outcomes by zipcode exposure and by actual exposure
tecos_iv_zip_hist_0_100_table <- filter(tecos_phys, is.na(ZIP_HIST_0_100) == F & quadrant == 0)
table(tecos_iv_zip_hist_0_100_table$OUTCOME,
      tecos_iv_zip_hist_0_100_table$EXPOSURE,
      tecos_iv_zip_hist_0_100_table$ZIP_HIST_0_100)

# ZCTA proportion sitagliptin ≥90% vs. ≤10%
tecos_phys <- tecos_phys %>%
  mutate(ZIP_HIST_10_90 = case_when(
    prop_exp_hist >= 0.9 ~ "1",
    prop_exp_hist <= 0.1 ~ "0",
    TRUE ~ NA_character_
  )
  )

# Check strength of instrument - linear regression
tecos_model_iv_zip_hist_strength_10_90 <- glm(data = tecos_phys,
                                              formula = EXPOSURE ~ ZIP_HIST_10_90 + as.factor(PS_ZIP_QUINTILE_14_17),
                                              na.action = na.exclude)
summary(tecos_model_iv_zip_hist_strength_10_90)
linearHypothesis(tecos_model_iv_zip_hist_strength_10_90, "ZIP_HIST_10_901 = 0", test = "F")

# Ordinary least squares - no adjustment for confounders
tecos_model_zip_hist_ols_10_90 <- glm(data = filter(tecos_phys, (prop_exp_hist <= 0.1 | prop_exp_hist >= 0.9)),
                                      formula = OUTCOME ~ EXPOSURE,
                                      na.action = na.exclude)
summary(tecos_model_zip_hist_ols_10_90)

# 2 stage least squares - no adjustment for confounders
tecos_model_zip_hist_tsls_10_90 <- ivreg(OUTCOME ~ EXPOSURE | ZIP_HIST_10_90,
                                         data = tecos_phys)
summary(tecos_model_zip_hist_tsls_10_90, vcov = sandwich, diagnostics = TRUE)

# Ordinary least squares - adjusted for PS
tecos_model_zip_hist_ols_10_90_adj <- glm(data = filter(tecos_phys, (prop_exp_hist <= 0.1 | prop_exp_hist >= 0.9)),
                                          formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                          na.action = na.exclude)
summary(tecos_model_zip_hist_ols_10_90_adj)

# 2 stage least squares - adjusted for PS
tecos_model_zip_hist_tsls_10_90_adj <- ivreg(OUTCOME ~ as.factor(PS_ZIP_QUINTILE_14_17) | EXPOSURE | ZIP_HIST_10_90,
                                             data = tecos_phys)
summary(tecos_model_zip_hist_tsls_10_90_adj, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - no adjustment for confounders
tecos_model_zip_hist_logistic_10_90 <- glm(data = filter(tecos_phys, (prop_exp_hist <= 0.1 | prop_exp_hist >= 0.9)),
                                           formula = OUTCOME ~ EXPOSURE,
                                           family = binomial(link = "logit"),
                                           na.action = na.exclude)
summary(tecos_model_zip_hist_logistic_10_90)
exp(tecos_model_zip_hist_logistic_10_90$coefficients)

# 2 stage residual inclusion - no adjustment for confounders
tecos_model_zip_hist_tsri_10_90 <- tsri(OUTCOME ~ EXPOSURE | ZIP_HIST_10_90,
                                        data = tecos_phys,
                                        link = "logit")
summary(tecos_model_zip_hist_tsri_10_90, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - adjusted for PS
tecos_model_zip_hist_logistic_10_90_adj <- glm(data = filter(tecos_phys, (prop_exp_hist <= 0.1 | prop_exp_hist >= 0.9)),
                                               formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                               family = binomial(link = "logit"),
                                               na.action = na.exclude)
summary(tecos_model_zip_hist_logistic_10_90_adj)
exp(tecos_model_zip_hist_logistic_10_90_adj$coefficients)

# 2 stage residual inclusion - adjusted for PS
tecos_model_zip_hist_tsri_10_90_adj <- tsri(OUTCOME ~ PS_ZIP_QUINTILE_14_17 | EXPOSURE | ZIP_HIST_10_90,
                                            data = tecos_phys,
                                            link = "logit")
summary(tecos_model_zip_hist_tsri_10_90_adj, vcov = sandwich, diagnostics = TRUE)

# Table 1 at the individual level
tecos_tableone_all_iv_zip_hist_10_90 <- CreateTableOne(vars = tecos_vars_ps, data = filter(tecos_phys, is.na(ZIP_HIST_10_90) == FALSE),
                                                       factorVars = tecos_factor_vars_ps,
                                                       strata = "ZIP_HIST_10_90", test = FALSE, smd = TRUE)
tecos_tableone_all_iv_zip_hist_10_90_export <- print(tecos_tableone_all_iv_zip_hist_10_90, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(tecos_tableone_all_iv_zip_hist_10_90_export, file = "/Volumes/PACS$/Cordes, Jack/tecos_tableone_all_iv_zip_hist_10_90.csv")

# Check balance improvement by percent change in Mahalanobis distance using treatment
tecos_zip_hist_iv_10_90_mah_dist <- tecos_phys[,c(1:13,19:21,23,34:35,38,40:41,48,52:53,55,58:60,66,70,80:81,87,89,91:92,100,102,105:106,108,113,117:120,122,130,135:136,141:150,152,158,164,167:173,175:182,187:190,205:206,210:211,220,226:229,231,237,244:245,247,249,254:258,262:264,271,296)]

tecos_mah_dist_total_zip_hist_iv_10_90 <- mahalanobis(colMeans(tecos_zip_hist_iv_10_90_mah_dist %>% filter(ZIP_HIST_10_90 == 1) %>% select(c(5:87,89:106))),
                                                      colMeans(tecos_zip_hist_iv_10_90_mah_dist %>% filter(ZIP_HIST_10_90 == 0) %>% select(c(5:87,89:106))),
                                                      cov(tecos_zip_hist_iv_10_90_mah_dist %>% select(c(5:87,89:106))),
                                                      tol = 1e-30)

((tecos_mah_dist_total_zip_hist_iv_10_90-tecos_mah_dist_total_14_17_treatment)/(tecos_mah_dist_total_14_17_treatment))*100

# Distribution of outcomes by ZCTA exposure and by actual exposure
tecos_iv_zip_hist_10_90_table <- filter(tecos_phys, is.na(ZIP_HIST_10_90) == F & quadrant == 0)
table(tecos_iv_zip_hist_10_90_table$OUTCOME,
      tecos_iv_zip_hist_10_90_table$EXPOSURE,
      tecos_iv_zip_hist_10_90_table$ZIP_HIST_10_90)

# ZCTA proportion sitagliptin ≥80% vs. ≤20%
tecos_phys <- tecos_phys %>%
  mutate(ZIP_HIST_20_80 = case_when(
    prop_exp_hist >= 0.8 ~ "1",
    prop_exp_hist <= 0.2 ~ "0",
    TRUE ~ NA_character_
  )
  )

# Check strength of instrument - linear regression
tecos_model_iv_zip_hist_strength_20_80 <- glm(data = tecos_phys,
                                              formula = EXPOSURE ~ ZIP_HIST_20_80 + as.factor(PS_ZIP_QUINTILE_14_17),
                                              na.action = na.exclude)
summary(tecos_model_iv_zip_hist_strength_20_80)
linearHypothesis(tecos_model_iv_zip_hist_strength_20_80, "ZIP_HIST_20_801 = 0", test = "F")

# Ordinary least squares - no adjustment for confounders
tecos_model_zip_hist_ols_20_80 <- glm(data = filter(tecos_phys, (prop_exp_hist <= 0.2 | prop_exp_hist >= 0.8)),
                                      formula = OUTCOME ~ EXPOSURE,
                                      na.action = na.exclude)
summary(tecos_model_zip_hist_ols_20_80)

# 2 stage least squares - no adjustment for confounders
tecos_model_zip_hist_tsls_20_80 <- ivreg(OUTCOME ~ EXPOSURE | ZIP_HIST_20_80,
                                         data = tecos_phys)
summary(tecos_model_zip_hist_tsls_20_80, vcov = sandwich, diagnostics = TRUE)

# Ordinary least squares - adjusted for PS
tecos_model_zip_hist_ols_20_80_adj <- glm(data = filter(tecos_phys, (prop_exp_hist <= 0.2 | prop_exp_hist >= 0.8)),
                                          formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                          na.action = na.exclude)
summary(tecos_model_zip_hist_ols_20_80_adj)

# 2 stage least squares - adjusted for PS
tecos_model_zip_hist_tsls_20_80_adj <- ivreg(OUTCOME ~ as.factor(PS_ZIP_QUINTILE_14_17) | EXPOSURE | ZIP_HIST_20_80,
                                             data = tecos_phys)
summary(tecos_model_zip_hist_tsls_20_80_adj, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - no adjustment for confounders
tecos_model_zip_hist_logistic_20_80 <- glm(data = filter(tecos_phys, (prop_exp_hist <= 0.2 | prop_exp_hist >= 0.8)),
                                           formula = OUTCOME ~ EXPOSURE,
                                           family = binomial(link = "logit"),
                                           na.action = na.exclude)
summary(tecos_model_zip_hist_logistic_20_80)
exp(tecos_model_zip_hist_logistic_20_80$coefficients)

# 2 stage residual inclusion - no adjustment for confounders
tecos_model_zip_hist_tsri_20_80 <- tsri(OUTCOME ~ EXPOSURE | ZIP_HIST_20_80,
                                        data = tecos_phys,
                                        link = "logit")
summary(tecos_model_zip_hist_tsri_20_80, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - adjusted for PS
tecos_model_zip_hist_logistic_20_80_adj <- glm(data = filter(tecos_phys, (prop_exp_hist <= 0.2 | prop_exp_hist >= 0.8)),
                                               formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                               family = binomial(link = "logit"),
                                               na.action = na.exclude)
summary(tecos_model_zip_hist_logistic_20_80_adj)
exp(tecos_model_zip_hist_logistic_20_80_adj$coefficients)

# 2 stage residual inclusion - adjusted for PS
tecos_model_zip_hist_tsri_20_80_adj <- tsri(OUTCOME ~ PS_ZIP_QUINTILE_14_17 | EXPOSURE | ZIP_HIST_20_80,
                                            data = tecos_phys,
                                            link = "logit")
summary(tecos_model_zip_hist_tsri_20_80_adj, vcov = sandwich, diagnostics = TRUE)

# Table 1 at the individual level
tecos_tableone_all_iv_zip_hist_20_80 <- CreateTableOne(vars = tecos_vars_ps, data = filter(tecos_phys, is.na(ZIP_HIST_20_80) == FALSE),
                                                       factorVars = tecos_factor_vars_ps,
                                                       strata = "ZIP_HIST_20_80", test = FALSE, smd = TRUE)
tecos_tableone_all_iv_zip_hist_20_80_export <- print(tecos_tableone_all_iv_zip_hist_20_80, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(tecos_tableone_all_iv_zip_hist_20_80_export, file = "/Volumes/PACS$/Cordes, Jack/tecos_tableone_all_iv_zip_hist_20_80.csv")

# Check balance improvement by percent change in Mahalanobis distance using treatment
tecos_zip_hist_iv_20_80_mah_dist <- tecos_phys[,c(1:13,19:21,23,34:35,38,40:41,48,52:53,55,58:60,66,70,80:81,87,89,91:92,100,102,105:106,108,113,117:120,122,130,135:136,141:150,152,158,164,167:173,175:182,187:190,205:206,210:211,220,226:229,231,237,244:245,247,249,254:258,262:264,271,297)]

tecos_mah_dist_total_zip_hist_iv_20_80 <- mahalanobis(colMeans(tecos_zip_hist_iv_20_80_mah_dist %>% filter(ZIP_HIST_20_80 == 1) %>% select(c(5:87,89:106))),
                                                      colMeans(tecos_zip_hist_iv_20_80_mah_dist %>% filter(ZIP_HIST_20_80 == 0) %>% select(c(5:87,89:106))),
                                                      cov(tecos_zip_hist_iv_20_80_mah_dist %>% select(c(5:87,89:106))),
                                                      tol = 1e-30)

((tecos_mah_dist_total_zip_hist_iv_20_80-tecos_mah_dist_total_14_17_treatment)/(tecos_mah_dist_total_14_17_treatment))*100

# Distribution of outcomes by ZCTA exposure and by actual exposure
tecos_iv_zip_hist_20_80_table <- filter(tecos_phys, is.na(ZIP_HIST_20_80) == F & quadrant == 0)
table(tecos_iv_zip_hist_20_80_table$OUTCOME,
      tecos_iv_zip_hist_20_80_table$EXPOSURE,
      tecos_iv_zip_hist_20_80_table$ZIP_HIST_20_80)

# ZCTA proportion sitagliptin ≥70% vs. ≤30%
tecos_phys <- tecos_phys %>%
  mutate(ZIP_HIST_30_70 = case_when(
    prop_exp_hist >= 0.7 ~ "1",
    prop_exp_hist <= 0.3 ~ "0",
    TRUE ~ NA_character_
  )
  )

# Check strength of instrument - linear regression
tecos_model_iv_zip_hist_strength_30_70 <- glm(data = tecos_phys,
                                              formula = EXPOSURE ~ ZIP_HIST_30_70 + as.factor(PS_ZIP_QUINTILE_14_17),
                                              na.action = na.exclude)
summary(tecos_model_iv_zip_hist_strength_30_70)
linearHypothesis(tecos_model_iv_zip_hist_strength_30_70, "ZIP_HIST_30_701 = 0", test = "F")

# Ordinary least squares - no adjustment for confounders
tecos_model_zip_hist_ols_30_70 <- glm(data = filter(tecos_phys, (prop_exp_hist <= 0.3 | prop_exp_hist >= 0.7)),
                                      formula = OUTCOME ~ EXPOSURE,
                                      na.action = na.exclude)
summary(tecos_model_zip_hist_ols_30_70)

# 2 stage least squares - no adjustment for confounders
tecos_model_zip_hist_tsls_30_70 <- ivreg(OUTCOME ~ EXPOSURE | ZIP_HIST_30_70,
                                         data = tecos_phys)
summary(tecos_model_zip_hist_tsls_30_70, vcov = sandwich, diagnostics = TRUE)

# Ordinary least squares - adjusted for PS
tecos_model_zip_hist_ols_30_70_adj <- glm(data = filter(tecos_phys, (prop_exp_hist <= 0.3 | prop_exp_hist >= 0.7)),
                                          formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                          na.action = na.exclude)
summary(tecos_model_zip_hist_ols_30_70_adj)

# 2 stage least squares - adjusted for PS
tecos_model_zip_hist_tsls_30_70_adj <- ivreg(OUTCOME ~ as.factor(PS_ZIP_QUINTILE_14_17) | EXPOSURE | ZIP_HIST_30_70,
                                             data = tecos_phys)
summary(tecos_model_zip_hist_tsls_30_70_adj, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - no adjustment for confounders
tecos_model_zip_hist_logistic_30_70 <- glm(data = filter(tecos_phys, (prop_exp_hist <= 0.3 | prop_exp_hist >= 0.7)),
                                           formula = OUTCOME ~ EXPOSURE,
                                           family = binomial(link = "logit"),
                                           na.action = na.exclude)
summary(tecos_model_zip_hist_logistic_30_70)
exp(tecos_model_zip_hist_logistic_30_70$coefficients)

# 2 stage residual inclusion - no adjustment for confounders
tecos_model_zip_hist_tsri_30_70 <- tsri(OUTCOME ~ EXPOSURE | ZIP_HIST_30_70,
                                        data = tecos_phys,
                                        link = "logit")
summary(tecos_model_zip_hist_tsri_30_70, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - adjusted for PS
tecos_model_zip_hist_logistic_30_70_adj <- glm(data = filter(tecos_phys, (prop_exp_hist <= 0.3 | prop_exp_hist >= 0.7)),
                                               formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                               family = binomial(link = "logit"),
                                               na.action = na.exclude)
summary(tecos_model_zip_hist_logistic_30_70_adj)
exp(tecos_model_zip_hist_logistic_30_70_adj$coefficients)

# 2 stage residual inclusion - adjusted for PS
tecos_model_zip_hist_tsri_30_70_adj <- tsri(OUTCOME ~ PS_ZIP_QUINTILE_14_17 | EXPOSURE | ZIP_HIST_30_70,
                                            data = tecos_phys,
                                            link = "logit")
summary(tecos_model_zip_hist_tsri_30_70_adj, vcov = sandwich, diagnostics = TRUE)

# Table 1 at the individual level
tecos_tableone_all_iv_zip_hist_30_70 <- CreateTableOne(vars = tecos_vars_ps, data = filter(tecos_phys, is.na(ZIP_HIST_30_70) == FALSE),
                                                       factorVars = tecos_factor_vars_ps,
                                                       strata = "ZIP_HIST_30_70", test = FALSE, smd = TRUE)
tecos_tableone_all_iv_zip_hist_30_70_export <- print(tecos_tableone_all_iv_zip_hist_30_70, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(tecos_tableone_all_iv_zip_hist_30_70_export, file = "/Volumes/PACS$/Cordes, Jack/tecos_tableone_all_iv_zip_hist_30_70.csv")

# Check balance improvement by percent change in Mahalanobis distance using treatment
tecos_zip_hist_iv_30_70_mah_dist <- tecos_phys[,c(1:13,19:21,23,34:35,38,40:41,48,52:53,55,58:60,66,70,80:81,87,89,91:92,100,102,105:106,108,113,117:120,122,130,135:136,141:150,152,158,164,167:173,175:182,187:190,205:206,210:211,220,226:229,231,237,244:245,247,249,254:258,262:264,271,298)]

tecos_mah_dist_total_zip_hist_iv_30_70 <- mahalanobis(colMeans(tecos_zip_hist_iv_30_70_mah_dist %>% filter(ZIP_HIST_30_70 == 1) %>% select(c(5:87,89:106))),
                                                      colMeans(tecos_zip_hist_iv_30_70_mah_dist %>% filter(ZIP_HIST_30_70 == 0) %>% select(c(5:87,89:106))),
                                                      cov(tecos_zip_hist_iv_30_70_mah_dist %>% select(c(5:87,89:106))),
                                                      tol = 1e-30)

((tecos_mah_dist_total_zip_hist_iv_30_70-tecos_mah_dist_total_14_17_treatment)/(tecos_mah_dist_total_14_17_treatment))*100

# Distribution of outcomes by ZCTA exposure and by actual exposure
tecos_iv_zip_hist_30_70_table <- filter(tecos_phys, is.na(ZIP_HIST_30_70) == F & quadrant == 0)
table(tecos_iv_zip_hist_30_70_table$OUTCOME,
      tecos_iv_zip_hist_30_70_table$EXPOSURE,
      tecos_iv_zip_hist_30_70_table$ZIP_HIST_30_70)

# ZCTA proportion sitagliptin ≥60% vs. ≤40%
tecos_phys <- tecos_phys %>%
  mutate(ZIP_HIST_40_60 = case_when(
    prop_exp_hist >= 0.6 ~ "1",
    prop_exp_hist <= 0.4 ~ "0",
    TRUE ~ NA_character_
  )
  )

# Check strength of instrument - linear regression
tecos_model_iv_zip_hist_strength_40_60 <- glm(data = tecos_phys,
                                              formula = EXPOSURE ~ ZIP_HIST_40_60 + as.factor(PS_ZIP_QUINTILE_14_17),
                                              na.action = na.exclude)
summary(tecos_model_iv_zip_hist_strength_40_60)
linearHypothesis(tecos_model_iv_zip_hist_strength_40_60, "ZIP_HIST_40_601 = 0", test = "F")

# Ordinary least squares - no adjustment for confounders
tecos_model_zip_hist_ols_40_60 <- glm(data = filter(tecos_phys, (prop_exp_hist <= 0.4 | prop_exp_hist >= 0.6)),
                                      formula = OUTCOME ~ EXPOSURE,
                                      na.action = na.exclude)
summary(tecos_model_zip_hist_ols_40_60)

# 2 stage least squares - no adjustment for confounders
tecos_model_zip_hist_tsls_40_60 <- ivreg(OUTCOME ~ EXPOSURE | ZIP_HIST_40_60,
                                         data = tecos_phys)
summary(tecos_model_zip_hist_tsls_40_60, vcov = sandwich, diagnostics = TRUE)

# Ordinary least squares - adjusted for PS
tecos_model_zip_hist_ols_40_60_adj <- glm(data = filter(tecos_phys, (prop_exp_hist <= 0.4 | prop_exp_hist >= 0.6)),
                                          formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                          na.action = na.exclude)
summary(tecos_model_zip_hist_ols_40_60_adj)

# 2 stage least squares - adjusted for PS
tecos_model_zip_hist_tsls_40_60_adj <- ivreg(OUTCOME ~ as.factor(PS_ZIP_QUINTILE_14_17) | EXPOSURE | ZIP_HIST_40_60,
                                             data = tecos_phys)
summary(tecos_model_zip_hist_tsls_40_60_adj, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - no adjustment for confounders
tecos_model_zip_hist_logistic_40_60 <- glm(data = filter(tecos_phys, (prop_exp_hist <= 0.4 | prop_exp_hist >= 0.6)),
                                           formula = OUTCOME ~ EXPOSURE,
                                           family = binomial(link = "logit"),
                                           na.action = na.exclude)
summary(tecos_model_zip_hist_logistic_40_60)
exp(tecos_model_zip_hist_logistic_40_60$coefficients)

# 2 stage residual inclusion - no adjustment for confounders
tecos_model_zip_hist_tsri_40_60 <- tsri(OUTCOME ~ EXPOSURE | ZIP_HIST_40_60,
                                        data = tecos_phys,
                                        link = "logit")
summary(tecos_model_zip_hist_tsri_40_60, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - adjusted for PS
tecos_model_zip_hist_logistic_40_60_adj <- glm(data = filter(tecos_phys, (prop_exp_hist <= 0.4 | prop_exp_hist >= 0.6)),
                                               formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                               family = binomial(link = "logit"),
                                               na.action = na.exclude)
summary(tecos_model_zip_hist_logistic_40_60_adj)
exp(tecos_model_zip_hist_logistic_40_60_adj$coefficients)

# 2 stage residual inclusion - adjusted for PS
tecos_model_zip_hist_tsri_40_60_adj <- tsri(OUTCOME ~ PS_ZIP_QUINTILE_14_17 | EXPOSURE | ZIP_HIST_40_60,
                                            data = tecos_phys,
                                            link = "logit")
summary(tecos_model_zip_hist_tsri_40_60_adj, vcov = sandwich, diagnostics = TRUE)

# Table 1 at the individual level
tecos_tableone_all_iv_zip_hist_40_60 <- CreateTableOne(vars = tecos_vars_ps, data = filter(tecos_phys, is.na(ZIP_HIST_40_60) == FALSE),
                                                       factorVars = tecos_factor_vars_ps,
                                                       strata = "ZIP_HIST_40_60", test = FALSE, smd = TRUE)
tecos_tableone_all_iv_zip_hist_40_60_export <- print(tecos_tableone_all_iv_zip_hist_40_60, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(tecos_tableone_all_iv_zip_hist_40_60_export, file = "/Volumes/PACS$/Cordes, Jack/tecos_tableone_all_iv_zip_hist_40_60.csv")

# Check balance improvement by percent change in Mahalanobis distance using treatment
tecos_zip_hist_iv_40_60_mah_dist <- tecos_phys[,c(1:13,19:21,23,34:35,38,40:41,48,52:53,55,58:60,66,70,80:81,87,89,91:92,100,102,105:106,108,113,117:120,122,130,135:136,141:150,152,158,164,167:173,175:182,187:190,205:206,210:211,220,226:229,231,237,244:245,247,249,254:258,262:264,271,299)]

tecos_mah_dist_total_zip_hist_iv_40_60 <- mahalanobis(colMeans(tecos_zip_hist_iv_40_60_mah_dist %>% filter(ZIP_HIST_40_60 == 1) %>% select(c(5:87,89:106))),
                                                      colMeans(tecos_zip_hist_iv_40_60_mah_dist %>% filter(ZIP_HIST_40_60 == 0) %>% select(c(5:87,89:106))),
                                                      cov(tecos_zip_hist_iv_40_60_mah_dist %>% select(c(5:87,89:106))),
                                                      tol = 1e-30)

((tecos_mah_dist_total_zip_hist_iv_40_60-tecos_mah_dist_total_14_17_treatment)/(tecos_mah_dist_total_14_17_treatment))*100

# Distribution of outcomes by ZCTA exposure and by actual exposure
tecos_iv_zip_hist_40_60_table <- filter(tecos_phys, is.na(ZIP_HIST_40_60) == F & quadrant == 0)
table(tecos_iv_zip_hist_40_60_table$OUTCOME,
      tecos_iv_zip_hist_40_60_table$EXPOSURE,
      tecos_iv_zip_hist_40_60_table$ZIP_HIST_40_60)

# ZCTA proportion sitagliptin ≥50% vs. <50%
tecos_phys <- tecos_phys %>%
  mutate(ZIP_HIST_50_50 = case_when(
    prop_exp_hist >= 0.5 ~ "1",
    prop_exp_hist < 0.5 ~ "0",
    TRUE ~ NA_character_
  )
  )

# Check strength of instrument - linear regression
tecos_model_iv_zip_hist_strength_50_50 <- glm(data = tecos_phys,
                                              formula = EXPOSURE ~ ZIP_HIST_50_50 + as.factor(PS_ZIP_QUINTILE_14_17),
                                              na.action = na.exclude)
summary(tecos_model_iv_zip_hist_strength_50_50)
linearHypothesis(tecos_model_iv_zip_hist_strength_50_50, "ZIP_HIST_50_501 = 0", test = "F")

# Ordinary least squares - no adjustment for confounders
tecos_model_zip_hist_ols_50_50 <- glm(data = filter(tecos_phys, (prop_exp_hist < 0.5 | prop_exp_hist >= 0.5)),
                                      formula = OUTCOME ~ EXPOSURE,
                                      na.action = na.exclude)
summary(tecos_model_zip_hist_ols_50_50)

# 2 stage least squares - no adjustment for confounders
tecos_model_zip_hist_tsls_50_50 <- ivreg(OUTCOME ~ EXPOSURE | ZIP_HIST_50_50,
                                         data = tecos_phys)
summary(tecos_model_zip_hist_tsls_50_50, vcov = sandwich, diagnostics = TRUE)

# Ordinary least squares - adjusted for PS
tecos_model_zip_hist_ols_50_50_adj <- glm(data = filter(tecos_phys, (prop_exp_hist < 0.5 | prop_exp_hist >= 0.5)),
                                          formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                          na.action = na.exclude)
summary(tecos_model_zip_hist_ols_50_50_adj)

# 2 stage least squares - adjusted for PS
tecos_model_zip_hist_tsls_50_50_adj <- ivreg(OUTCOME ~ as.factor(PS_ZIP_QUINTILE_14_17) | EXPOSURE | ZIP_HIST_50_50,
                                             data = tecos_phys)
summary(tecos_model_zip_hist_tsls_50_50_adj, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - no adjustment for confounders
tecos_model_zip_hist_logistic_50_50 <- glm(data = filter(tecos_phys, (prop_exp_hist < 0.5 | prop_exp_hist >= 0.5)),
                                           formula = OUTCOME ~ EXPOSURE,
                                           family = binomial(link = "logit"),
                                           na.action = na.exclude)
summary(tecos_model_zip_hist_logistic_50_50)
exp(tecos_model_zip_hist_logistic_50_50$coefficients)

# 2 stage residual inclusion - no adjustment for confounders
tecos_model_zip_hist_tsri_50_50 <- tsri(OUTCOME ~ EXPOSURE | ZIP_HIST_50_50,
                                        data = tecos_phys,
                                        link = "logit")
summary(tecos_model_zip_hist_tsri_50_50, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - adjusted for PS
tecos_model_zip_hist_logistic_50_50_adj <- glm(data = filter(tecos_phys, (prop_exp_hist < 0.5 | prop_exp_hist >= 0.5)),
                                               formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                               family = binomial(link = "logit"),
                                               na.action = na.exclude)
summary(tecos_model_zip_hist_logistic_50_50_adj)
exp(tecos_model_zip_hist_logistic_50_50_adj$coefficients)

# 2 stage residual inclusion - adjusted for PS
tecos_model_zip_hist_tsri_50_50_adj <- tsri(OUTCOME ~ PS_ZIP_QUINTILE_14_17 | EXPOSURE | ZIP_HIST_50_50,
                                            data = tecos_phys,
                                            link = "logit")
summary(tecos_model_zip_hist_tsri_50_50_adj, vcov = sandwich, diagnostics = TRUE)

# Table 1 at the individual level
tecos_tableone_all_iv_zip_hist_50_50 <- CreateTableOne(vars = tecos_vars_ps, data = filter(tecos_phys, is.na(ZIP_HIST_50_50) == FALSE),
                                                       factorVars = tecos_factor_vars_ps,
                                                       strata = "ZIP_HIST_50_50", test = FALSE, smd = TRUE)
tecos_tableone_all_iv_zip_hist_50_50_export <- print(tecos_tableone_all_iv_zip_hist_50_50, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(tecos_tableone_all_iv_zip_hist_50_50_export, file = "/Volumes/PACS$/Cordes, Jack/tecos_tableone_all_iv_zip_hist_50_50.csv")

# Check balance improvement by percent change in Mahalanobis distance using treatment
tecos_zip_hist_iv_50_50_mah_dist <- tecos_phys[,c(1:13,19:21,23,34:35,38,40:41,48,52:53,55,58:60,66,70,80:81,87,89,91:92,100,102,105:106,108,113,117:120,122,130,135:136,141:150,152,158,164,167:173,175:182,187:190,205:206,210:211,220,226:229,231,237,244:245,247,249,254:258,262:264,271,300)]

tecos_mah_dist_total_zip_hist_iv_50_50 <- mahalanobis(colMeans(tecos_zip_hist_iv_50_50_mah_dist %>% filter(ZIP_HIST_50_50 == 1) %>% select(c(5:87,89:106))),
                                                      colMeans(tecos_zip_hist_iv_50_50_mah_dist %>% filter(ZIP_HIST_50_50 == 0) %>% select(c(5:87,89:106))),
                                                      cov(tecos_zip_hist_iv_50_50_mah_dist %>% select(c(5:87,89:106))),
                                                      tol = 1e-30)

((tecos_mah_dist_total_zip_hist_iv_50_50-tecos_mah_dist_total_14_17_treatment)/(tecos_mah_dist_total_14_17_treatment))*100

# Distribution of outcomes by ZCTA exposure and by actual exposure
tecos_iv_zip_hist_50_50_table <- filter(tecos_phys, is.na(ZIP_HIST_50_50) == F & quadrant == 0)
table(tecos_iv_zip_hist_50_50_table$OUTCOME,
      tecos_iv_zip_hist_50_50_table$EXPOSURE,
      tecos_iv_zip_hist_50_50_table$ZIP_HIST_50_50)

# IV using physician preference
# Distribution of physicians across ZCTAs
tecos_zcta_phys_frequency <- tecos_phys %>%
  group_by(ZCTA_USE) %>%
  dplyr::summarize(count_zcta_phys = n_distinct(prscrbr_id))

summary(tecos_zcta_phys_frequency$count_zcta_phys)

# Add prescriber ZCTA counts to individuals
tecos_phys <- left_join(tecos_phys, tecos_zcta_phys_frequency, by = "ZCTA_USE")

# Distribution of prescribers across ZCTAs with at least 5 individuals
tecos_phys_5 <- filter(tecos_phys, count > 4)

# Filter dataset to only ZCTAs
tecos_phys_5_zcta <- distinct(tecos_phys_5, ZCTA_USE, .keep_all = TRUE)
summary(tecos_phys_5_zcta$count_zcta_phys)

tecos_zcta_phys_dist_plot <- tecos_phys_5_zcta %>%
  ggplot() + 
  geom_histogram(aes(count_zcta_phys), binwidth = 5, color = "black", fill = "white") + 
  geom_vline(xintercept = 12.11, color = "red") +
  scale_y_continuous(labels = label_comma()) +
  coord_cartesian(xlim = c(0, 175)) +
  labs(title = "Sitagliptin Cohort Prescriber Distribution across ZCTAs\n5+ Individuals in ZCTA", x = "Count of Prescribers", y = "Count of ZCTAs") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(hjust = 0.5))
tecos_zcta_phys_dist_grob <- gtable_add_grob(ggplot_gtable(ggplot_build(tecos_zcta_phys_dist_plot)), a_label, t = 1, l = 4, b = 6)
ggsave("tecos_zcta_phys_dist.pdf", tecos_zcta_phys_dist_grob)

# Spatial distribution of prescribers by ZCTA with at least 5 cohort individuals in ZCTA
# Join frequency data to ZCTA polygons
zcta_shp_tecos <- merge(zcta_shp_tecos, tecos_phys_5_zcta, by = "ZCTA_USE")

# Create Lower 48 ZCTA and state shapefiles
zcta_shp_tecos_48 <- st_transform(zcta_shp_tecos %>% filter(str_starts(ZCTA_USE, "0") |
                                                              str_starts(ZCTA_USE, "1") |
                                                              str_starts(ZCTA_USE, "2") |
                                                              str_starts(ZCTA_USE, "3") |
                                                              str_starts(ZCTA_USE, "4") |
                                                              str_starts(ZCTA_USE, "5") |
                                                              str_starts(ZCTA_USE, "6") |
                                                              str_starts(ZCTA_USE, "7") |
                                                              str_starts(ZCTA_USE, "8") |
                                                              str_starts(ZCTA_USE, "90") |
                                                              str_starts(ZCTA_USE, "91") |
                                                              str_starts(ZCTA_USE, "92") |
                                                              str_starts(ZCTA_USE, "93") |
                                                              str_starts(ZCTA_USE, "94") |
                                                              str_starts(ZCTA_USE, "95") |
                                                              str_starts(ZCTA_USE, "961") |
                                                              str_starts(ZCTA_USE, "97") |
                                                              str_starts(ZCTA_USE, "980") |
                                                              str_starts(ZCTA_USE, "981") |
                                                              str_starts(ZCTA_USE, "982") |
                                                              str_starts(ZCTA_USE, "983") |
                                                              str_starts(ZCTA_USE, "984") |
                                                              str_starts(ZCTA_USE, "985") |
                                                              str_starts(ZCTA_USE, "986") |
                                                              str_starts(ZCTA_USE, "987") |
                                                              str_starts(ZCTA_USE, "988") |
                                                              str_starts(ZCTA_USE, "989") |
                                                              str_starts(ZCTA_USE, "990") |
                                                              str_starts(ZCTA_USE, "991") |
                                                              str_starts(ZCTA_USE, "992") |
                                                              str_starts(ZCTA_USE, "993") |
                                                              str_starts(ZCTA_USE, "994")), 5070)
state_shp_48 <- st_transform(state_shp %>% filter(str_starts(STATEFP, "01") |
                                                    str_starts(STATEFP, "04") |
                                                    str_starts(STATEFP, "05") |
                                                    str_starts(STATEFP, "06") |
                                                    str_starts(STATEFP, "08") |
                                                    str_starts(STATEFP, "09") |
                                                    str_starts(STATEFP, "10") |
                                                    str_starts(STATEFP, "11") |
                                                    str_starts(STATEFP, "12") |
                                                    str_starts(STATEFP, "13") |
                                                    str_starts(STATEFP, "16") |
                                                    str_starts(STATEFP, "17") |
                                                    str_starts(STATEFP, "18") |
                                                    str_starts(STATEFP, "19") |
                                                    str_starts(STATEFP, "20") |
                                                    str_starts(STATEFP, "21") |
                                                    str_starts(STATEFP, "22") |
                                                    str_starts(STATEFP, "23") |
                                                    str_starts(STATEFP, "24") |
                                                    str_starts(STATEFP, "25") |
                                                    str_starts(STATEFP, "26") |
                                                    str_starts(STATEFP, "27") |
                                                    str_starts(STATEFP, "28") |
                                                    str_starts(STATEFP, "29") |
                                                    str_starts(STATEFP, "30") |
                                                    str_starts(STATEFP, "31") |
                                                    str_starts(STATEFP, "32") |
                                                    str_starts(STATEFP, "33") |
                                                    str_starts(STATEFP, "34") |
                                                    str_starts(STATEFP, "35") |
                                                    str_starts(STATEFP, "36") |
                                                    str_starts(STATEFP, "37") |
                                                    str_starts(STATEFP, "38") |
                                                    str_starts(STATEFP, "39") |
                                                    str_starts(STATEFP, "40") |
                                                    str_starts(STATEFP, "41") |
                                                    str_starts(STATEFP, "42") |
                                                    str_starts(STATEFP, "44") |
                                                    str_starts(STATEFP, "45") |
                                                    str_starts(STATEFP, "46") |
                                                    str_starts(STATEFP, "47") |
                                                    str_starts(STATEFP, "48") |
                                                    str_starts(STATEFP, "49") |
                                                    str_starts(STATEFP, "50") |
                                                    str_starts(STATEFP, "51") |
                                                    str_starts(STATEFP, "53") |
                                                    str_starts(STATEFP, "54") |
                                                    str_starts(STATEFP, "55") |
                                                    str_starts(STATEFP, "56")), 5070)
state_shp_48 <- st_buffer(state_shp_48, dist = 0)

# Create Alaska ZCTA and state shapefiles
zcta_shp_tecos_AK <- st_transform(zcta_shp_tecos %>% filter(str_starts(ZCTA_USE, "995") | 
                                                              str_starts(ZCTA_USE, "996") |
                                                              str_starts(ZCTA_USE, "997") |
                                                              str_starts(ZCTA_USE, "998") |
                                                              str_starts(ZCTA_USE, "999")), 3338)
state_shp_AK <- st_transform(state_shp %>% filter(str_starts(STATEFP, "02")), 3338)
state_shp_AK <- st_buffer(state_shp_AK, dist = 0)

# Create Hawaii ZCTA and state shapefiles
zcta_shp_tecos_HI <- st_transform(zcta_shp_tecos %>% filter(str_starts(ZCTA_USE, "967") | 
                                                              str_starts(ZCTA_USE, "968")), 32603)
state_shp_HI <- st_transform(state_shp %>% filter(str_starts(STATEFP, "15")), 32603)
state_shp_HI <- st_buffer(state_shp_HI, dist = 0)

# Map of number of prescribers by ZCTA
#options(device = "X11")

zcta_map_tecos_48_phys <- ggplot() +
  geom_sf(data = zcta_shp_tecos_48, aes(fill = cut(count_zcta_phys, c(0,1,5,10,25,175))), color = NA, show.legend = TRUE) +
  scale_fill_brewer(palette = "Spectral", name = "Sitagliptin Cohort\nCount of Prescribers\nper ZCTA", direction = -1) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  geom_sf(data = state_shp_48, fill = NA, color = "black", size = 0.01)
#zcta_map_tecos_48_phys

zcta_map_tecos_AK_phys <- ggplot() +
  geom_sf(data = zcta_shp_tecos_AK, aes(fill = cut(count_zcta_phys, c(0,1,5,10,25,175))), color = NA, show.legend = FALSE) +
  scale_fill_brewer(palette = "Spectral", direction = -1) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_blank()) + 
  geom_sf(data = st_crop(state_shp_AK, xmin = -1000000,
                         xmax = 1491822,
                         ymin = 414200,
                         ymax = 2378425), fill = NA, color = "black", size = 0.01)
#zcta_map_tecos_AK_phys

zcta_map_tecos_HI_phys <- ggplot() +
  geom_sf(data = zcta_shp_tecos_HI, aes(fill = cut(count_zcta_phys, c(0,1,5,10,25,175))), color = NA, show.legend = FALSE) +
  scale_fill_brewer(palette = "Spectral", direction = -1) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_blank()) +
  geom_sf(data = st_crop(state_shp_HI, xmin = 1000000,
                         xmax = 1573922,
                         ymin = 2117066,
                         ymax = 3214762), fill = NA, color = "black", size = 0.01)
#zcta_map_tecos_HI_phys

zcta_map_tecos_phys <- ggdraw() +
  draw_plot(zcta_map_tecos_48_phys) +
  draw_plot(zcta_map_tecos_AK_phys, x = 0, y = 0, width = 0.3, height = 0.28) +
  draw_plot(zcta_map_tecos_HI_phys, x = 0.4, y = 0, width = 0.3, height = 0.18) +
  draw_grob(a_label, x = 0.03, y = 0.15)
#zcta_map_tecos_phys
ggsave("/Volumes/PACS$/Cordes, Jack/tecos_zcta_map_phys.pdf", zcta_map_tecos_phys)

# Instantaneous physician preference IV analysis
# Convert PS_ZIP_QUINTILE_14_17 to factor for tsri
tecos_phys$PS_ZIP_QUINTILE_14_17 <- as.factor(tecos_phys$PS_ZIP_QUINTILE_14_17)

# Sort on prescriber ID and CED
tecos_phys <- tecos_phys[order(tecos_phys$prscrbr_id, tecos_phys$ENTRYDATE),]

# Create empty instantaneous physician preference column
tecos_phys <- tecos_phys %>% mutate(PHYS_PREF_INST = NA)

# Create instantaneous physician preference variable based on most recent prescription in cohort
for (i in 2:nrow(tecos_phys)){
  if (tecos_phys$prscrbr_id[i] == tecos_phys$prscrbr_id[i-1]){
    tecos_phys$PHYS_PREF_INST[i] <- tecos_phys$EXPOSURE[i-1]
  }
}

# Distribution of global historical physician preference proportions
summary(tecos_phys$PHYS_PREF_INST)

# Check strength of instrument - linear regression
tecos_model_iv_strength_phys_pref_inst <- glm(data = tecos_phys,
                                              formula = EXPOSURE ~ PHYS_PREF_INST + as.factor(PS_ZIP_QUINTILE_14_17),
                                              na.action = na.exclude)
summary(tecos_model_iv_strength_phys_pref_inst)
linearHypothesis(tecos_model_iv_strength_phys_pref_inst, "PHYS_PREF_INST = 0", test = "F")

# Ordinary least squares - no adjustment for confounders
tecos_model_ols_phys_pref_inst <- glm(data = filter(tecos_phys, is.na(PHYS_PREF_INST) == FALSE),
                                      formula = OUTCOME ~ EXPOSURE,
                                      na.action = na.exclude)
summary(tecos_model_ols_phys_pref_inst)

# 2 stage least squares - no adjustment for confounders
tecos_model_tsls_phys_pref_inst <- ivreg(OUTCOME ~ EXPOSURE | PHYS_PREF_INST,
                                         data = filter(tecos_phys, is.na(PHYS_PREF_INST) == FALSE))
summary(tecos_model_tsls_phys_pref_inst, vcov = sandwich, diagnostics = TRUE)

# Ordinary least squares - adjusted for PS
tecos_model_ols_phys_pref_inst_adj <- glm(data = filter(tecos_phys, is.na(PHYS_PREF_INST) == FALSE),
                                          formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                          na.action = na.exclude)
summary(tecos_model_ols_phys_pref_inst_adj)

# 2 stage least squares - adjusted for PS
tecos_model_tsls_phys_pref_inst_adj <- ivreg(OUTCOME ~ as.factor(PS_ZIP_QUINTILE_14_17) | EXPOSURE | PHYS_PREF_INST,
                                             data = filter(tecos_phys, is.na(PHYS_PREF_INST) == FALSE))
summary(tecos_model_tsls_phys_pref_inst_adj, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - no adjustment for confounders
tecos_model_logistic_phys_pref_inst <- glm(data = filter(tecos_phys, is.na(PHYS_PREF_INST) == FALSE),
                                           formula = OUTCOME ~ EXPOSURE,
                                           family = binomial(link = "logit"),
                                           na.action = na.exclude)
summary(tecos_model_logistic_phys_pref_inst)
exp(tecos_model_logistic_phys_pref_inst$coefficients)

# 2 stage residual inclusion - no adjustment for confounders
tecos_model_tsri_phys_pref_inst <- tsri(OUTCOME ~ EXPOSURE | PHYS_PREF_INST,
                                        data = filter(tecos_phys, is.na(PHYS_PREF_INST) == FALSE),
                                        link = "logit")
summary(tecos_model_tsri_phys_pref_inst, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - adjusted for PS
tecos_model_logistic_phys_pref_inst_adj <- glm(data = filter(tecos_phys, is.na(PHYS_PREF_INST) == FALSE),
                                               formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                               family = binomial(link = "logit"),
                                               na.action = na.exclude)
summary(tecos_model_logistic_phys_pref_inst_adj)
exp(tecos_model_logistic_phys_pref_inst_adj$coefficients)

# 2 stage residual inclusion - adjusted for PS
tecos_model_tsri_phys_pref_inst_adj <- tsri(OUTCOME ~ PS_ZIP_QUINTILE_14_17 | EXPOSURE | PHYS_PREF_INST,
                                            data = filter(tecos_phys, is.na(PHYS_PREF_INST) == FALSE),
                                            link = "logit")
summary(tecos_model_tsri_phys_pref_inst_adj, vcov = sandwich, diagnostics = TRUE)

# Table 1 at the individual level
tecos_tableone_all_iv_phys_pref_inst <- CreateTableOne(vars = tecos_vars_ps, data = filter(tecos_phys, is.na(PHYS_PREF_INST) == FALSE),
                                                       factorVars = tecos_factor_vars_ps,
                                                       strata = "PHYS_PREF_INST", test = FALSE, smd = TRUE)
tecos_tableone_all_iv_phys_pref_inst_export <- print(tecos_tableone_all_iv_phys_pref_inst, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(tecos_tableone_all_iv_phys_pref_inst_export, file = "/Volumes/PACS$/Cordes, Jack/tecos_tableone_all_iv_phys_pref_inst.csv")

# Check balance improvement by percent change in Mahalanobis distance using treatment
tecos_phys_iv_inst_mah_dist <- tecos_phys[,c(1:13,19:21,23,34:35,38,40:41,48,52:53,55,58:60,66,70,80:81,87,89,91:92,100,102,105:106,108,113,117:120,122,130,135:136,141:150,152,158,164,167:173,175:182,187:190,205:206,210:211,220,226:229,231,237,244:245,247,249,254:258,262:264,302)]

tecos_mah_dist_total_phys_iv_inst <- mahalanobis(colMeans(tecos_phys_iv_inst_mah_dist %>% filter(PHYS_PREF_INST == 1) %>% select(c(5:87,89:106))),
                                                 colMeans(tecos_phys_iv_inst_mah_dist %>% filter(PHYS_PREF_INST == 0) %>% select(c(5:87,89:106))),
                                                 cov(tecos_phys_iv_inst_mah_dist %>% select(c(5:87,89:106))),
                                                 tol = 1e-30)

((tecos_mah_dist_total_phys_iv_inst-tecos_mah_dist_total_14_17_treatment)/(tecos_mah_dist_total_14_17_treatment))*100

# Distribution of outcomes by physician preference and by actual exposure
tecos_iv_phys_pref_inst_table <- filter(tecos_phys, is.na(PHYS_PREF_INST) == FALSE)
table(tecos_iv_phys_pref_inst_table$OUTCOME,
      tecos_iv_phys_pref_inst_table$EXPOSURE,
      tecos_iv_phys_pref_inst_table$PHYS_PREF_INST)

# Global historical physician preference IV analysis
# Sort on prescriber ID and CED
tecos_phys <- tecos_phys[order(tecos_phys$prscrbr_id, tecos_phys$ENTRYDATE),]

# Create empty global historical physician preference column
tecos_phys <- tecos_phys %>% mutate(prop_exp_phys = NA)

# Create empty counter column
tecos_phys <- tecos_phys %>% mutate(counter = NA)

# Create empty cumulative sum column
tecos_phys <- tecos_phys %>% mutate(cum_sum = NA)

# Create global historical physician preference variable based on average of all previous prescriptions
for (i in 2:nrow(tecos_phys)){
  if ((tecos_phys$prscrbr_id[i] == tecos_phys$prscrbr_id[i-1]) & (is.na(tecos_phys$prop_exp_phys[i-1]) == TRUE)){
    tecos_phys$prop_exp_phys[i] <- tecos_phys$EXPOSURE[i-1]
  }
  if ((tecos_phys$prscrbr_id[i] == tecos_phys$prscrbr_id[i-1]) & (is.na(tecos_phys$prop_exp_phys[i-1]) == TRUE)){
    tecos_phys$counter[i] <- 1
  }
  if ((tecos_phys$prscrbr_id[i] == tecos_phys$prscrbr_id[i-1]) & (is.na(tecos_phys$prop_exp_phys[i-1]) == TRUE)){
    tecos_phys$cum_sum[i] <- tecos_phys$EXPOSURE[i-1]
  }
  if ((tecos_phys$prscrbr_id[i] == tecos_phys$prscrbr_id[i-1]) & (is.na(tecos_phys$prop_exp_phys[i-1]) == FALSE)){
    tecos_phys$counter[i] <- tecos_phys$counter[i-1]+1
  }
  if ((tecos_phys$prscrbr_id[i] == tecos_phys$prscrbr_id[i-1]) & (is.na(tecos_phys$prop_exp_phys[i-1]) == FALSE)){
    tecos_phys$cum_sum[i] <- tecos_phys$cum_sum[i-1] + tecos_phys$EXPOSURE[i-1]
  }
  if ((tecos_phys$prscrbr_id[i] == tecos_phys$prscrbr_id[i-1]) & (is.na(tecos_phys$prop_exp_phys[i-1]) == FALSE)){
    tecos_phys$prop_exp_phys[i] <- tecos_phys$cum_sum[i]/tecos_phys$counter[i]
  }
}

# Distribution of global historical physician preference proportions
summary(tecos_phys$prop_exp_phys)

# Create dataset limited to non-missing physician preference IV
tecos_phys_global <- filter(tecos_phys, is.na(prop_exp_phys) == FALSE)

tecos_phys_global_pref_dist_plot <- tecos_phys_global %>%
  ggplot() + 
  geom_histogram(aes(prop_exp_phys), binwidth = 0.05, color = "black", fill = "white") + 
  geom_vline(xintercept = 0.27, color = "red") +
  scale_y_continuous(labels = label_comma()) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(title = "Sitagliptin Cohort\nCumulative Prescriber Preference for Each Individual ", x = "Proportion Sitagliptin", y = "Count of Individuals") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(hjust = 0.5))
tecos_phys_global_pref_dist_grob <- gtable_add_grob(ggplot_gtable(ggplot_build(tecos_phys_global_pref_dist_plot)), a_label, t = 1, l = 4, b = 6)
ggsave("tecos_phys_global_pref_dist.pdf", tecos_phys_global_pref_dist_grob)

# For IV with proportion sitagliptin cutoffs 100% vs. 0%, instrument is collapsed with exposure
# because no variation in individual-level exposure
# So becomes simple regression comparing individuals with 100% to 0% prescribers

# Prescriber proportion sitagliptin 100% vs. 0%
tecos_phys <- tecos_phys %>%
  mutate(PHYS_PREF_GLOBAL_0_100 = case_when(
    prop_exp_phys == 1 ~ "1",
    prop_exp_phys == 0 ~ "0",
    TRUE ~ NA_character_
  )
  )

# Ordinary least squares - no adjustment for confounders
tecos_model_phys_pref_global_ols_0_100 <- glm(data = filter(tecos_phys, (prop_exp_phys == 0 | prop_exp_phys == 1)),
                                              formula = OUTCOME ~ EXPOSURE,
                                              na.action = na.exclude)
summary(tecos_model_phys_pref_global_ols_0_100)

# Ordinary least squares - adjusted for PS
tecos_model_phys_pref_global_ols_0_100_adj <- glm(data = filter(tecos_phys, (prop_exp_phys == 0 | prop_exp_phys == 1)),
                                                  formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                                  na.action = na.exclude)
summary(tecos_model_phys_pref_global_ols_0_100_adj)

# Logistic regression - no adjustment for confounders
tecos_model_phys_pref_global_logistic_0_100 <- glm(data = filter(tecos_phys, (prop_exp_phys == 0 | prop_exp_phys == 1)),
                                                   formula = OUTCOME ~ EXPOSURE,
                                                   family = binomial(link = "logit"),
                                                   na.action = na.exclude)
summary(tecos_model_phys_pref_global_logistic_0_100)
exp(tecos_model_phys_pref_global_logistic_0_100$coefficients)

# Logistic regression - adjusted for PS
tecos_model_phys_pref_global_logistic_0_100_adj <- glm(data = filter(tecos_phys, (prop_exp_phys == 0 | prop_exp_phys == 1)),
                                                       formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                                       family = binomial(link = "logit"),
                                                       na.action = na.exclude)
summary(tecos_model_phys_pref_global_logistic_0_100_adj)
exp(tecos_model_phys_pref_global_logistic_0_100_adj$coefficients)

# Table 1 at the individual level
tecos_tableone_all_iv_phys_pref_global_0_100 <- CreateTableOne(vars = tecos_vars_ps, data = filter(tecos_phys, is.na(PHYS_PREF_GLOBAL_0_100) == FALSE),
                                                               factorVars = tecos_factor_vars_ps,
                                                               strata = "PHYS_PREF_GLOBAL_0_100", test = FALSE, smd = TRUE)
tecos_tableone_all_iv_phys_pref_global_0_100_export <- print(tecos_tableone_all_iv_phys_pref_global_0_100, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(tecos_tableone_all_iv_phys_pref_global_0_100_export, file = "/Volumes/PACS$/Cordes, Jack/tecos_tableone_all_iv_phys_pref_global_0_100.csv")

# Check balance improvement by percent change in Mahalanobis distance using treatment
tecos_phys_iv_0_100_mah_dist <- tecos_phys[,c(1:13,19:21,23,34:35,38,40:41,48,52:53,55,58:60,66,70,80:81,87,89,91:92,100,102,105:106,108,113,117:120,122,130,135:136,141:150,152,158,164,167:173,175:182,187:190,205:206,210:211,220,226:229,231,237,244:245,247,249,254:258,262:264,304)]

tecos_mah_dist_total_phys_iv_0_100 <- mahalanobis(colMeans(tecos_phys_iv_0_100_mah_dist %>% filter(PHYS_PREF_GLOBAL_0_100 == 1) %>% select(c(5:87,89:106))),
                                                  colMeans(tecos_phys_iv_0_100_mah_dist %>% filter(PHYS_PREF_GLOBAL_0_100 == 0) %>% select(c(5:87,89:106))),
                                                  cov(tecos_phys_iv_0_100_mah_dist %>% select(c(5:87,89:106))),
                                                  tol = 1e-30)

((tecos_mah_dist_total_phys_iv_0_100-tecos_mah_dist_total_14_17_treatment)/(tecos_mah_dist_total_14_17_treatment))*100

# Distribution of outcomes by zipcode exposure and by actual exposure
tecos_iv_phys_pref_global_0_100_table <- filter(tecos_phys, is.na(PHYS_PREF_GLOBAL_0_100) == FALSE)
table(tecos_iv_phys_pref_global_0_100_table$OUTCOME,
      tecos_iv_phys_pref_global_0_100_table$EXPOSURE,
      tecos_iv_phys_pref_global_0_100_table$PHYS_PREF_GLOBAL_0_100)

# Prescriber proportion sitagliptin ≥90% vs. ≤10%
tecos_phys <- tecos_phys %>%
  mutate(PHYS_PREF_GLOBAL_10_90 = case_when(
    prop_exp_phys >= 0.9 ~ "1",
    prop_exp_phys <= 0.1 ~ "0",
    TRUE ~ NA_character_
  )
  )

# Check strength of instrument - linear regression
tecos_model_iv_phys_pref_global_strength_10_90 <- glm(data = tecos_phys,
                                                      formula = EXPOSURE ~ PHYS_PREF_GLOBAL_10_90 + as.factor(PS_ZIP_QUINTILE_14_17),
                                                      na.action = na.exclude)
summary(tecos_model_iv_phys_pref_global_strength_10_90)
linearHypothesis(tecos_model_iv_phys_pref_global_strength_10_90, "PHYS_PREF_GLOBAL_10_901 = 0", test = "F")

# Ordinary least squares - no adjustment for confounders
tecos_model_phys_pref_global_ols_10_90 <- glm(data = filter(tecos_phys, (prop_exp_phys <= 0.1 | prop_exp_phys >= 0.9)),
                                              formula = OUTCOME ~ EXPOSURE,
                                              na.action = na.exclude)
summary(tecos_model_phys_pref_global_ols_10_90)

# 2 stage least squares - no adjustment for confounders
tecos_model_phys_pref_global_tsls_10_90 <- ivreg(OUTCOME ~ EXPOSURE | PHYS_PREF_GLOBAL_10_90,
                                                 data = tecos_phys)
summary(tecos_model_phys_pref_global_tsls_10_90, vcov = sandwich, diagnostics = TRUE)

# Ordinary least squares - adjusted for PS
tecos_model_phys_pref_global_ols_10_90_adj <- glm(data = filter(tecos_phys, (prop_exp_phys <= 0.1 | prop_exp_phys >= 0.9)),
                                                  formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                                  na.action = na.exclude)
summary(tecos_model_phys_pref_global_ols_10_90_adj)

# 2 stage least squares - adjusted for PS
tecos_model_phys_pref_global_tsls_10_90_adj <- ivreg(OUTCOME ~ as.factor(PS_ZIP_QUINTILE_14_17) | EXPOSURE | PHYS_PREF_GLOBAL_10_90,
                                                     data = tecos_phys)
summary(tecos_model_phys_pref_global_tsls_10_90_adj, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - no adjustment for confounders
tecos_model_phys_pref_global_logistic_10_90 <- glm(data = filter(tecos_phys, (prop_exp_phys <= 0.1 | prop_exp_phys >= 0.9)),
                                                   formula = OUTCOME ~ EXPOSURE,
                                                   family = binomial(link = "logit"),
                                                   na.action = na.exclude)
summary(tecos_model_phys_pref_global_logistic_10_90)
exp(tecos_model_phys_pref_global_logistic_10_90$coefficients)

# 2 stage residual inclusion - no adjustment for confounders
tecos_model_phys_pref_global_tsri_10_90 <- tsri(OUTCOME ~ EXPOSURE | PHYS_PREF_GLOBAL_10_90,
                                                data = tecos_phys,
                                                link = "logit")
summary(tecos_model_phys_pref_global_tsri_10_90, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - adjusted for PS
tecos_model_phys_pref_global_logistic_10_90_adj <- glm(data = filter(tecos_phys, (prop_exp_phys <= 0.1 | prop_exp_phys >= 0.9)),
                                                       formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                                       family = binomial(link = "logit"),
                                                       na.action = na.exclude)
summary(tecos_model_phys_pref_global_logistic_10_90_adj)
exp(tecos_model_phys_pref_global_logistic_10_90_adj$coefficients)

# 2 stage residual inclusion - adjusted for PS
tecos_model_phys_pref_global_tsri_10_90_adj <- tsri(OUTCOME ~ PS_ZIP_QUINTILE_14_17 | EXPOSURE | PHYS_PREF_GLOBAL_10_90,
                                                    data = tecos_phys,
                                                    link = "logit")
summary(tecos_model_phys_pref_global_tsri_10_90_adj, vcov = sandwich, diagnostics = TRUE)

# Table 1 at the individual level
tecos_tableone_all_iv_phys_pref_global_10_90 <- CreateTableOne(vars = tecos_vars_ps, data = filter(tecos_phys, is.na(PHYS_PREF_GLOBAL_10_90) == FALSE),
                                                               factorVars = tecos_factor_vars_ps,
                                                               strata = "PHYS_PREF_GLOBAL_10_90", test = FALSE, smd = TRUE)
tecos_tableone_all_iv_phys_pref_global_10_90_export <- print(tecos_tableone_all_iv_phys_pref_global_10_90, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(tecos_tableone_all_iv_phys_pref_global_10_90_export, file = "/Volumes/PACS$/Cordes, Jack/tecos_tableone_all_iv_phys_pref_global_10_90.csv")

# Check balance improvement by percent change in Mahalanobis distance using treatment
tecos_phys_iv_10_90_mah_dist <- tecos_phys[,c(1:13,19:21,23,34:35,38,40:41,48,52:53,55,58:60,66,70,80:81,87,89,91:92,100,102,105:106,108,113,117:120,122,130,135:136,141:150,152,158,164,167:173,175:182,187:190,205:206,210:211,220,226:229,231,237,244:245,247,249,254:258,262:264,305)]

tecos_mah_dist_total_phys_iv_10_90 <- mahalanobis(colMeans(tecos_phys_iv_10_90_mah_dist %>% filter(PHYS_PREF_GLOBAL_10_90 == 1) %>% select(c(5:87,89:106))),
                                                  colMeans(tecos_phys_iv_10_90_mah_dist %>% filter(PHYS_PREF_GLOBAL_10_90 == 0) %>% select(c(5:87,89:106))),
                                                  cov(tecos_phys_iv_10_90_mah_dist %>% select(c(5:87,89:106))),
                                                  tol = 1e-30)

((tecos_mah_dist_total_phys_iv_10_90-tecos_mah_dist_total_14_17_treatment)/(tecos_mah_dist_total_14_17_treatment))*100

# Distribution of outcomes by prescriber exposure and by actual exposure
tecos_iv_phys_pref_global_10_90_table <- filter(tecos_phys, is.na(PHYS_PREF_GLOBAL_10_90) == FALSE)
table(tecos_iv_phys_pref_global_10_90_table$OUTCOME,
      tecos_iv_phys_pref_global_10_90_table$EXPOSURE,
      tecos_iv_phys_pref_global_10_90_table$PHYS_PREF_GLOBAL_10_90)

# Prescriber proportion sitagliptin ≥80% vs. ≤20%
tecos_phys <- tecos_phys %>%
  mutate(PHYS_PREF_GLOBAL_20_80 = case_when(
    prop_exp_phys >= 0.8 ~ "1",
    prop_exp_phys <= 0.2 ~ "0",
    TRUE ~ NA_character_
  )
  )

# Check strength of instrument - linear regression
tecos_model_iv_phys_pref_global_strength_20_80 <- glm(data = tecos_phys,
                                                      formula = EXPOSURE ~ PHYS_PREF_GLOBAL_20_80 + as.factor(PS_ZIP_QUINTILE_14_17),
                                                      na.action = na.exclude)
summary(tecos_model_iv_phys_pref_global_strength_20_80)
linearHypothesis(tecos_model_iv_phys_pref_global_strength_20_80, "PHYS_PREF_GLOBAL_20_801 = 0", test = "F")

# Ordinary least squares - no adjustment for confounders
tecos_model_phys_pref_global_ols_20_80 <- glm(data = filter(tecos_phys, (prop_exp_phys <= 0.2 | prop_exp_phys >= 0.8)),
                                              formula = OUTCOME ~ EXPOSURE,
                                              na.action = na.exclude)
summary(tecos_model_phys_pref_global_ols_20_80)

# 2 stage least squares - no adjustment for confounders
tecos_model_phys_pref_global_tsls_20_80 <- ivreg(OUTCOME ~ EXPOSURE | PHYS_PREF_GLOBAL_20_80,
                                                 data = tecos_phys)
summary(tecos_model_phys_pref_global_tsls_20_80, vcov = sandwich, diagnostics = TRUE)

# Ordinary least squares - adjusted for PS
tecos_model_phys_pref_global_ols_20_80_adj <- glm(data = filter(tecos_phys, (prop_exp_phys <= 0.2 | prop_exp_phys >= 0.8)),
                                                  formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                                  na.action = na.exclude)
summary(tecos_model_phys_pref_global_ols_20_80_adj)

# 2 stage least squares - adjusted for PS
tecos_model_phys_pref_global_tsls_20_80_adj <- ivreg(OUTCOME ~ as.factor(PS_ZIP_QUINTILE_14_17) | EXPOSURE | PHYS_PREF_GLOBAL_20_80,
                                                     data = tecos_phys)
summary(tecos_model_phys_pref_global_tsls_20_80_adj, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - no adjustment for confounders
tecos_model_phys_pref_global_logistic_20_80 <- glm(data = filter(tecos_phys, (prop_exp_phys <= 0.2 | prop_exp_phys >= 0.8)),
                                                   formula = OUTCOME ~ EXPOSURE,
                                                   family = binomial(link = "logit"),
                                                   na.action = na.exclude)
summary(tecos_model_phys_pref_global_logistic_20_80)
exp(tecos_model_phys_pref_global_logistic_20_80$coefficients)

# 2 stage residual inclusion - no adjustment for confounders
tecos_model_phys_pref_global_tsri_20_80 <- tsri(OUTCOME ~ EXPOSURE | PHYS_PREF_GLOBAL_20_80,
                                                data = tecos_phys,
                                                link = "logit")
summary(tecos_model_phys_pref_global_tsri_20_80, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - adjusted for PS
tecos_model_phys_pref_global_logistic_20_80_adj <- glm(data = filter(tecos_phys, (prop_exp_phys <= 0.2 | prop_exp_phys >= 0.8)),
                                                       formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                                       family = binomial(link = "logit"),
                                                       na.action = na.exclude)
summary(tecos_model_phys_pref_global_logistic_20_80_adj)
exp(tecos_model_phys_pref_global_logistic_20_80_adj$coefficients)

# 2 stage residual inclusion - adjusted for PS
tecos_model_phys_pref_global_tsri_20_80_adj <- tsri(OUTCOME ~ PS_ZIP_QUINTILE_14_17 | EXPOSURE | PHYS_PREF_GLOBAL_20_80,
                                                    data = tecos_phys,
                                                    link = "logit")
summary(tecos_model_phys_pref_global_tsri_20_80_adj, vcov = sandwich, diagnostics = TRUE)

# Table 1 at the individual level
tecos_tableone_all_iv_phys_pref_global_20_80 <- CreateTableOne(vars = tecos_vars_ps, data = filter(tecos_phys, is.na(PHYS_PREF_GLOBAL_20_80) == FALSE),
                                                               factorVars = tecos_factor_vars_ps,
                                                               strata = "PHYS_PREF_GLOBAL_20_80", test = FALSE, smd = TRUE)
tecos_tableone_all_iv_phys_pref_global_20_80_export <- print(tecos_tableone_all_iv_phys_pref_global_20_80, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(tecos_tableone_all_iv_phys_pref_global_20_80_export, file = "/Volumes/PACS$/Cordes, Jack/tecos_tableone_all_iv_phys_pref_global_20_80.csv")

# Check balance improvement by percent change in Mahalanobis distance using treatment
tecos_phys_iv_20_80_mah_dist <- tecos_phys[,c(1:13,19:21,23,34:35,38,40:41,48,52:53,55,58:60,66,70,80:81,87,89,91:92,100,102,105:106,108,113,117:120,122,130,135:136,141:150,152,158,164,167:173,175:182,187:190,205:206,210:211,220,226:229,231,237,244:245,247,249,254:258,262:264,306)]

tecos_mah_dist_total_phys_iv_20_80 <- mahalanobis(colMeans(tecos_phys_iv_20_80_mah_dist %>% filter(PHYS_PREF_GLOBAL_20_80 == 1) %>% select(c(5:87,89:106))),
                                                  colMeans(tecos_phys_iv_20_80_mah_dist %>% filter(PHYS_PREF_GLOBAL_20_80 == 0) %>% select(c(5:87,89:106))),
                                                  cov(tecos_phys_iv_20_80_mah_dist %>% select(c(5:87,89:106))),
                                                  tol = 1e-30)

((tecos_mah_dist_total_phys_iv_20_80-tecos_mah_dist_total_14_17_treatment)/(tecos_mah_dist_total_14_17_treatment))*100

# Distribution of outcomes by prescriber exposure and by actual exposure
tecos_iv_phys_pref_global_20_80_table <- filter(tecos_phys, is.na(PHYS_PREF_GLOBAL_20_80) == FALSE)
table(tecos_iv_phys_pref_global_20_80_table$OUTCOME,
      tecos_iv_phys_pref_global_20_80_table$EXPOSURE,
      tecos_iv_phys_pref_global_20_80_table$PHYS_PREF_GLOBAL_20_80)

# Prescriber proportion sitagliptin ≥70% vs. ≤30%
tecos_phys <- tecos_phys %>%
  mutate(PHYS_PREF_GLOBAL_30_70 = case_when(
    prop_exp_phys >= 0.7 ~ "1",
    prop_exp_phys <= 0.3 ~ "0",
    TRUE ~ NA_character_
  )
  )

# Check strength of instrument - linear regression
tecos_model_iv_phys_pref_global_strength_30_70 <- glm(data = tecos_phys,
                                                      formula = EXPOSURE ~ PHYS_PREF_GLOBAL_30_70 + as.factor(PS_ZIP_QUINTILE_14_17),
                                                      na.action = na.exclude)
summary(tecos_model_iv_phys_pref_global_strength_30_70)
linearHypothesis(tecos_model_iv_phys_pref_global_strength_30_70, "PHYS_PREF_GLOBAL_30_701 = 0", test = "F")

# Ordinary least squares - no adjustment for confounders
tecos_model_phys_pref_global_ols_30_70 <- glm(data = filter(tecos_phys, (prop_exp_phys <= 0.3 | prop_exp_phys >= 0.7)),
                                              formula = OUTCOME ~ EXPOSURE,
                                              na.action = na.exclude)
summary(tecos_model_phys_pref_global_ols_30_70)

# 2 stage least squares - no adjustment for confounders
tecos_model_phys_pref_global_tsls_30_70 <- ivreg(OUTCOME ~ EXPOSURE | PHYS_PREF_GLOBAL_30_70,
                                                 data = tecos_phys)
summary(tecos_model_phys_pref_global_tsls_30_70, vcov = sandwich, diagnostics = TRUE)

# Ordinary least squares - adjusted for PS
tecos_model_phys_pref_global_ols_30_70_adj <- glm(data = filter(tecos_phys, (prop_exp_phys <= 0.3 | prop_exp_phys >= 0.7)),
                                                  formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                                  na.action = na.exclude)
summary(tecos_model_phys_pref_global_ols_30_70_adj)

# 2 stage least squares - adjusted for PS
tecos_model_phys_pref_global_tsls_30_70_adj <- ivreg(OUTCOME ~ as.factor(PS_ZIP_QUINTILE_14_17) | EXPOSURE | PHYS_PREF_GLOBAL_30_70,
                                                     data = tecos_phys)
summary(tecos_model_phys_pref_global_tsls_30_70_adj, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - no adjustment for confounders
tecos_model_phys_pref_global_logistic_30_70 <- glm(data = filter(tecos_phys, (prop_exp_phys <= 0.3 | prop_exp_phys >= 0.7)),
                                                   formula = OUTCOME ~ EXPOSURE,
                                                   family = binomial(link = "logit"),
                                                   na.action = na.exclude)
summary(tecos_model_phys_pref_global_logistic_30_70)
exp(tecos_model_phys_pref_global_logistic_30_70$coefficients)

# 2 stage residual inclusion - no adjustment for confounders
tecos_model_phys_pref_global_tsri_30_70 <- tsri(OUTCOME ~ EXPOSURE | PHYS_PREF_GLOBAL_30_70,
                                                data = tecos_phys,
                                                link = "logit")
summary(tecos_model_phys_pref_global_tsri_30_70, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - adjusted for PS
tecos_model_phys_pref_global_logistic_30_70_adj <- glm(data = filter(tecos_phys, (prop_exp_phys <= 0.3 | prop_exp_phys >= 0.7)),
                                                       formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                                       family = binomial(link = "logit"),
                                                       na.action = na.exclude)
summary(tecos_model_phys_pref_global_logistic_30_70_adj)
exp(tecos_model_phys_pref_global_logistic_30_70_adj$coefficients)

# 2 stage residual inclusion - prescribers adjusted for PS
tecos_model_phys_pref_global_tsri_30_70_adj <- tsri(OUTCOME ~ PS_ZIP_QUINTILE_14_17 | EXPOSURE | PHYS_PREF_GLOBAL_30_70,
                                                    data = tecos_phys,
                                                    link = "logit")
summary(tecos_model_phys_pref_global_tsri_30_70_adj, vcov = sandwich, diagnostics = TRUE)

# Table 1 at the individual level
tecos_tableone_all_iv_phys_pref_global_30_70 <- CreateTableOne(vars = tecos_vars_ps, data = filter(tecos_phys, is.na(PHYS_PREF_GLOBAL_30_70) == FALSE),
                                                               factorVars = tecos_factor_vars_ps,
                                                               strata = "PHYS_PREF_GLOBAL_30_70", test = FALSE, smd = TRUE)
tecos_tableone_all_iv_phys_pref_global_30_70_export <- print(tecos_tableone_all_iv_phys_pref_global_30_70, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(tecos_tableone_all_iv_phys_pref_global_30_70_export, file = "/Volumes/PACS$/Cordes, Jack/tecos_tableone_all_iv_phys_pref_global_30_70.csv")

# Check balance improvement by percent change in Mahalanobis distance using treatment
tecos_phys_iv_30_70_mah_dist <- tecos_phys[,c(1:13,19:21,23,34:35,38,40:41,48,52:53,55,58:60,66,70,80:81,87,89,91:92,100,102,105:106,108,113,117:120,122,130,135:136,141:150,152,158,164,167:173,175:182,187:190,205:206,210:211,220,226:229,231,237,244:245,247,249,254:258,262:264,307)]

tecos_mah_dist_total_phys_iv_30_70 <- mahalanobis(colMeans(tecos_phys_iv_30_70_mah_dist %>% filter(PHYS_PREF_GLOBAL_30_70 == 1) %>% select(c(5:87,89:106))),
                                                  colMeans(tecos_phys_iv_30_70_mah_dist %>% filter(PHYS_PREF_GLOBAL_30_70 == 0) %>% select(c(5:87,89:106))),
                                                  cov(tecos_phys_iv_30_70_mah_dist %>% select(c(5:87,89:106))),
                                                  tol = 1e-30)

((tecos_mah_dist_total_phys_iv_30_70-tecos_mah_dist_total_14_17_treatment)/(tecos_mah_dist_total_14_17_treatment))*100

# Distribution of outcomes by prescriber exposure and by actual exposure
tecos_iv_phys_pref_global_30_70_table <- filter(tecos_phys, is.na(PHYS_PREF_GLOBAL_30_70) == FALSE)
table(tecos_iv_phys_pref_global_30_70_table$OUTCOME,
      tecos_iv_phys_pref_global_30_70_table$EXPOSURE,
      tecos_iv_phys_pref_global_30_70_table$PHYS_PREF_GLOBAL_30_70)

# Prescriber proportion sitagliptin ≥60% vs. ≤40%
tecos_phys <- tecos_phys %>%
  mutate(PHYS_PREF_GLOBAL_40_60 = case_when(
    prop_exp_phys >= 0.6 ~ "1",
    prop_exp_phys <= 0.4 ~ "0",
    TRUE ~ NA_character_
  )
  )

# Check strength of instrument - linear regression
tecos_model_iv_phys_pref_global_strength_40_60 <- glm(data = tecos_phys,
                                                      formula = EXPOSURE ~ PHYS_PREF_GLOBAL_40_60 + as.factor(PS_ZIP_QUINTILE_14_17),
                                                      na.action = na.exclude)
summary(tecos_model_iv_phys_pref_global_strength_40_60)
linearHypothesis(tecos_model_iv_phys_pref_global_strength_40_60, "PHYS_PREF_GLOBAL_40_601 = 0", test = "F")

# Ordinary least squares - no adjustment for confounders
tecos_model_phys_pref_global_ols_40_60 <- glm(data = filter(tecos_phys, (prop_exp_phys <= 0.4 | prop_exp_phys >= 0.6)),
                                              formula = OUTCOME ~ EXPOSURE,
                                              na.action = na.exclude)
summary(tecos_model_phys_pref_global_ols_40_60)

# 2 stage least squares - no adjustment for confounders
tecos_model_phys_pref_global_tsls_40_60 <- ivreg(OUTCOME ~ EXPOSURE | PHYS_PREF_GLOBAL_40_60,
                                                 data = tecos_phys)
summary(tecos_model_phys_pref_global_tsls_40_60, vcov = sandwich, diagnostics = TRUE)

# Ordinary least squares - adjusted for PS
tecos_model_phys_pref_global_ols_40_60_adj <- glm(data = filter(tecos_phys, (prop_exp_phys <= 0.4 | prop_exp_phys >= 0.6)),
                                                  formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                                  na.action = na.exclude)
summary(tecos_model_phys_pref_global_ols_40_60_adj)

# 2 stage least squares - adjusted for PS
tecos_model_phys_pref_global_tsls_40_60_adj <- ivreg(OUTCOME ~ as.factor(PS_ZIP_QUINTILE_14_17) | EXPOSURE | PHYS_PREF_GLOBAL_40_60,
                                                     data = tecos_phys)
summary(tecos_model_phys_pref_global_tsls_40_60_adj, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - no adjustment for confounders
tecos_model_phys_pref_global_logistic_40_60 <- glm(data = filter(tecos_phys, (prop_exp_phys <= 0.4 | prop_exp_phys >= 0.6)),
                                                   formula = OUTCOME ~ EXPOSURE,
                                                   family = binomial(link = "logit"),
                                                   na.action = na.exclude)
summary(tecos_model_phys_pref_global_logistic_40_60)
exp(tecos_model_phys_pref_global_logistic_40_60$coefficients)

# 2 stage residual inclusion - no adjustment for confounders
tecos_model_phys_pref_global_tsri_40_60 <- tsri(OUTCOME ~ EXPOSURE | PHYS_PREF_GLOBAL_40_60,
                                                data = tecos_phys,
                                                link = "logit")
summary(tecos_model_phys_pref_global_tsri_40_60, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - adjusted for PS
tecos_model_phys_pref_global_logistic_40_60_adj <- glm(data = filter(tecos_phys, (prop_exp_phys <= 0.4 | prop_exp_phys >= 0.6)),
                                                       formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                                       family = binomial(link = "logit"),
                                                       na.action = na.exclude)
summary(tecos_model_phys_pref_global_logistic_40_60_adj)
exp(tecos_model_phys_pref_global_logistic_40_60_adj$coefficients)

# 2 stage residual inclusion - adjusted for PS
tecos_model_phys_pref_global_tsri_40_60_adj <- tsri(OUTCOME ~ PS_ZIP_QUINTILE_14_17 | EXPOSURE | PHYS_PREF_GLOBAL_40_60,
                                                    data = tecos_phys,
                                                    link = "logit")
summary(tecos_model_phys_pref_global_tsri_40_60_adj, vcov = sandwich, diagnostics = TRUE)

# Table 1 at the individual level
tecos_tableone_all_iv_phys_pref_global_40_60 <- CreateTableOne(vars = tecos_vars_ps, data = filter(tecos_phys, is.na(PHYS_PREF_GLOBAL_40_60) == FALSE),
                                                               factorVars = tecos_factor_vars_ps,
                                                               strata = "PHYS_PREF_GLOBAL_40_60", test = FALSE, smd = TRUE)
tecos_tableone_all_iv_phys_pref_global_40_60_export <- print(tecos_tableone_all_iv_phys_pref_global_40_60, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(tecos_tableone_all_iv_phys_pref_global_40_60_export, file = "/Volumes/PACS$/Cordes, Jack/tecos_tableone_all_iv_phys_pref_global_40_60.csv")

# Check balance improvement by percent change in Mahalanobis distance using treatment
tecos_phys_iv_40_60_mah_dist <- tecos_phys[,c(1:13,19:21,23,34:35,38,40:41,48,52:53,55,58:60,66,70,80:81,87,89,91:92,100,102,105:106,108,113,117:120,122,130,135:136,141:150,152,158,164,167:173,175:182,187:190,205:206,210:211,220,226:229,231,237,244:245,247,249,254:258,262:264,308)]

tecos_mah_dist_total_phys_iv_40_60 <- mahalanobis(colMeans(tecos_phys_iv_40_60_mah_dist %>% filter(PHYS_PREF_GLOBAL_40_60 == 1) %>% select(c(5:87,89:106))),
                                                  colMeans(tecos_phys_iv_40_60_mah_dist %>% filter(PHYS_PREF_GLOBAL_40_60 == 0) %>% select(c(5:87,89:106))),
                                                  cov(tecos_phys_iv_40_60_mah_dist %>% select(c(5:87,89:106))),
                                                  tol = 1e-30)

((tecos_mah_dist_total_phys_iv_40_60-tecos_mah_dist_total_14_17_treatment)/(tecos_mah_dist_total_14_17_treatment))*100

# Distribution of outcomes by prescriber exposure and by actual exposure
tecos_iv_phys_pref_global_40_60_table <- filter(tecos_phys, is.na(PHYS_PREF_GLOBAL_40_60) == FALSE)
table(tecos_iv_phys_pref_global_40_60_table$OUTCOME,
      tecos_iv_phys_pref_global_40_60_table$EXPOSURE,
      tecos_iv_phys_pref_global_40_60_table$PHYS_PREF_GLOBAL_40_60)

# Prescriber proportion sitagliptin ≥50% vs. <50%
tecos_phys <- tecos_phys %>%
  mutate(PHYS_PREF_GLOBAL_50_50 = case_when(
    prop_exp_phys >= 0.5 ~ "1",
    prop_exp_phys < 0.5 ~ "0",
    TRUE ~ NA_character_
  )
  )

# Check strength of instrument - linear regression
tecos_model_iv_phys_pref_global_strength_50_50 <- glm(data = tecos_phys,
                                                      formula = EXPOSURE ~ PHYS_PREF_GLOBAL_50_50 + as.factor(PS_ZIP_QUINTILE_14_17),
                                                      na.action = na.exclude)
summary(tecos_model_iv_phys_pref_global_strength_50_50)
linearHypothesis(tecos_model_iv_phys_pref_global_strength_50_50, "PHYS_PREF_GLOBAL_50_501 = 0", test = "F")

# Ordinary least squares - no adjustment for confounders
tecos_model_phys_pref_global_ols_50_50 <- glm(data = filter(tecos_phys, (prop_exp_phys < 0.5 | prop_exp_phys >= 0.5)),
                                              formula = OUTCOME ~ EXPOSURE,
                                              na.action = na.exclude)
summary(tecos_model_phys_pref_global_ols_50_50)

# 2 stage least squares - no adjustment for confounders
tecos_model_phys_pref_global_tsls_50_50 <- ivreg(OUTCOME ~ EXPOSURE | PHYS_PREF_GLOBAL_50_50,
                                                 data = tecos_phys)
summary(tecos_model_phys_pref_global_tsls_50_50, vcov = sandwich, diagnostics = TRUE)

# Ordinary least squares - adjusted for PS
tecos_model_phys_pref_global_ols_50_50_adj <- glm(data = filter(tecos_phys, (prop_exp_phys < 0.5 | prop_exp_phys >= 0.5)),
                                                  formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                                  na.action = na.exclude)
summary(tecos_model_phys_pref_global_ols_50_50_adj)

# 2 stage least squares - adjusted for PS
tecos_model_phys_pref_global_tsls_50_50_adj <- ivreg(OUTCOME ~ as.factor(PS_ZIP_QUINTILE_14_17) | EXPOSURE | PHYS_PREF_GLOBAL_50_50,
                                                     data = tecos_phys)
summary(tecos_model_phys_pref_global_tsls_50_50_adj, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - no adjustment for confounders
tecos_model_phys_pref_global_logistic_50_50 <- glm(data = filter(tecos_phys, (prop_exp_phys < 0.5 | prop_exp_phys >= 0.5)),
                                                   formula = OUTCOME ~ EXPOSURE,
                                                   family = binomial(link = "logit"),
                                                   na.action = na.exclude)
summary(tecos_model_phys_pref_global_logistic_50_50)
exp(tecos_model_phys_pref_global_logistic_50_50$coefficients)

# 2 stage residual inclusion - no adjustment for confounders
tecos_model_phys_pref_global_tsri_50_50 <- tsri(OUTCOME ~ EXPOSURE | PHYS_PREF_GLOBAL_50_50,
                                                data = tecos_phys,
                                                link = "logit")
summary(tecos_model_phys_pref_global_tsri_50_50, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - adjusted for PS
tecos_model_phys_pref_global_logistic_50_50_adj <- glm(data = filter(tecos_phys, (prop_exp_phys < 0.5 | prop_exp_phys >= 0.5)),
                                                       formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                                       family = binomial(link = "logit"),
                                                       na.action = na.exclude)
summary(tecos_model_phys_pref_global_logistic_50_50_adj)
exp(tecos_model_phys_pref_global_logistic_50_50_adj$coefficients)

# 2 stage residual inclusion - adjusted for PS
tecos_model_phys_pref_global_tsri_50_50_adj <- tsri(OUTCOME ~ PS_ZIP_QUINTILE_14_17 | EXPOSURE | PHYS_PREF_GLOBAL_50_50,
                                                    data = tecos_phys,
                                                    link = "logit")
summary(tecos_model_phys_pref_global_tsri_50_50_adj, vcov = sandwich, diagnostics = TRUE)

# Table 1 at the individual level
tecos_tableone_all_iv_phys_pref_global_50_50 <- CreateTableOne(vars = tecos_vars_ps, data = filter(tecos_phys, is.na(PHYS_PREF_GLOBAL_50_50) == FALSE),
                                                               factorVars = tecos_factor_vars_ps,
                                                               strata = "PHYS_PREF_GLOBAL_50_50", test = FALSE, smd = TRUE)
tecos_tableone_all_iv_phys_pref_global_50_50_export <- print(tecos_tableone_all_iv_phys_pref_global_50_50, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(tecos_tableone_all_iv_phys_pref_global_50_50_export, file = "/Volumes/PACS$/Cordes, Jack/tecos_tableone_all_iv_phys_pref_global_50_50.csv")

# Check balance improvement by percent change in Mahalanobis distance using treatment
tecos_phys_iv_50_50_mah_dist <- tecos_phys[,c(1:13,19:21,23,34:35,38,40:41,48,52:53,55,58:60,66,70,80:81,87,89,91:92,100,102,105:106,108,113,117:120,122,130,135:136,141:150,152,158,164,167:173,175:182,187:190,205:206,210:211,220,226:229,231,237,244:245,247,249,254:258,262:264,309)]

tecos_mah_dist_total_phys_iv_50_50 <- mahalanobis(colMeans(tecos_phys_iv_50_50_mah_dist %>% filter(PHYS_PREF_GLOBAL_50_50 == 1) %>% select(c(5:87,89:106))),
                                                  colMeans(tecos_phys_iv_50_50_mah_dist %>% filter(PHYS_PREF_GLOBAL_50_50 == 0) %>% select(c(5:87,89:106))),
                                                  cov(tecos_phys_iv_50_50_mah_dist %>% select(c(5:87,89:106))),
                                                  tol = 1e-30)

((tecos_mah_dist_total_phys_iv_50_50-tecos_mah_dist_total_14_17_treatment)/(tecos_mah_dist_total_14_17_treatment))*100

# Distribution of outcomes by prescriber exposure and by actual exposure
tecos_iv_phys_pref_global_50_50_table <- filter(tecos_phys, is.na(PHYS_PREF_GLOBAL_50_50) == FALSE)
table(tecos_iv_phys_pref_global_50_50_table$OUTCOME,
      tecos_iv_phys_pref_global_50_50_table$EXPOSURE,
      tecos_iv_phys_pref_global_50_50_table$PHYS_PREF_GLOBAL_50_50)

# SAVOR-TIMI Primary ####
savor <- read.csv("/Volumes/PACS$/Cordes, Jack/Data_without_Prescriber_ID/JCA01_TIMI_primary_all_followup.csv", colClasses = c(ZIP_CD = "character"))

savor$OUTCOME <- ifelse(savor$OUTCOME == "false", 0, 1)

savor_exposure <- read.csv("/Volumes/PACS$/Cordes, Jack/SAVOR-TIMI/outcome_Primary_Composite_Outcome__stroke__MI__all-cause_mortality_from_vital_/primary/all_followup.csv")[,1:2]

savor_exposure$EXPOSURE <- ifelse(savor_exposure$EXPOSURE == "R", 0, 1)

savor <- left_join(savor, savor_exposure, by = "PID")

savor$YEAR <- as.numeric(substr(savor$ENTRY_CLASS, 1,4))
savor$ZIP <- savor$ZIP_CD

# Add prescriber ID
savor_rx_id <- read.csv("/Volumes/PACS$/Cordes, Jack/Files with Zip Codes/jca02_TIMI_primary_all_followup.csv")[,c(1,217:219)]

savor <- left_join(savor, savor_rx_id, by = "PID")

# Convert cohort zip codes to ZCTA
savor <- merge(savor, zip_zcta_crosswalk, by=c("ZIP","YEAR"))

# Join ACS zip code data to cohort
savor <- merge(savor, zipcode_2012_2017, by=c("ZCTA_USE","YEAR"))

# Check for missingness in ACS data and remove individuals with missing zipcode information
colSums(is.na(savor))

# Remove individuals younger than 65 years old
# Total of 15130 individuals comprising 1.73% of the cohort
# 237 individuals had both missing zipcode data and were <65 years old
savor <- savor %>% filter(RUN1_ENTRY_COVARIATE_147 >= 65)

# Remove individuals with missing zipcode ACS data
# Total of 26812 individuals comprising 3.07% of the cohort
savor <- savor[complete.cases(savor), ]

# Remove individuals with zipcodes from Puerto Rico
# Total of 39 individuals comprising 0.0047% of the cohort
savor <- savor %>% filter(as.numeric(ZCTA_USE) >= 1000)

# Rename Aetion variables
savor <- savor %>% rename(AGE_CAT	=	RUN1_ENTRY_COVARIATE_1,
                          SEX	=	RUN1_ENTRY_COVARIATE_2,
                          RACE	=	RUN1_ENTRY_COVARIATE_3,
                          REGION	=	RUN1_ENTRY_COVARIATE_4,
                          OBESITY	=	RUN1_ENTRY_COVARIATE_5,
                          OVERWEIGHT	=	RUN1_ENTRY_COVARIATE_6,
                          SMOKE	=	RUN1_ENTRY_COVARIATE_7,
                          ALCOHOL	=	RUN1_ENTRY_COVARIATE_8,
                          DRUG	=	RUN1_ENTRY_COVARIATE_9,
                          DM_RETINOPATHY	=	RUN1_ENTRY_COVARIATE_10,
                          DM_OPHTHALMIC	=	RUN1_ENTRY_COVARIATE_11,
                          RETINAL_DETACH	=	RUN1_ENTRY_COVARIATE_12,
                          RETINAL_LASER	=	RUN1_ENTRY_COVARIATE_13,
                          DM_NEUROPATHY2	=	RUN1_ENTRY_COVARIATE_14,
                          DM_NEUROPATHY3	=	RUN1_ENTRY_COVARIATE_15,
                          HYPOGLYCEMIA	=	RUN1_ENTRY_COVARIATE_16,
                          HYPERGLYCEMIA	=	RUN1_ENTRY_COVARIATE_17,
                          ELECTROLYTE_ACID_BASE	=	RUN1_ENTRY_COVARIATE_18,
                          DM_KETOACIDOSIS	=	RUN1_ENTRY_COVARIATE_19,
                          HONK	=	RUN1_ENTRY_COVARIATE_20,
                          DM_PERIPHERAL	=	RUN1_ENTRY_COVARIATE_21,
                          DM_FOOT	=	RUN1_ENTRY_COVARIATE_22,
                          GANGRENE	=	RUN1_ENTRY_COVARIATE_23,
                          LOWER_AMPUTATION	=	RUN1_ENTRY_COVARIATE_24,
                          OSTEOMYELITIS	=	RUN1_ENTRY_COVARIATE_25,
                          SKIN_INFECTION	=	RUN1_ENTRY_COVARIATE_26,
                          ED	=	RUN1_ENTRY_COVARIATE_27,
                          DM_COMPLICATION	=	RUN1_ENTRY_COVARIATE_28,
                          DM_NO_COMPLICATION	=	RUN1_ENTRY_COVARIATE_29,
                          HYPERTENSION	=	RUN1_ENTRY_COVARIATE_30,
                          HYPERLIPIDEMIA	=	RUN1_ENTRY_COVARIATE_31,
                          IHD	=	RUN1_ENTRY_COVARIATE_32,
                          ACUTE_MI	=	RUN1_ENTRY_COVARIATE_33,
                          UNSTABLE_ANGINA	=	RUN1_ENTRY_COVARIATE_34,
                          OLD_MI	=	RUN1_ENTRY_COVARIATE_35,
                          STABLE_ANGINA	=	RUN1_ENTRY_COVARIATE_36,
                          CORONARY_ATHEROSCLEROSIS	=	RUN1_ENTRY_COVARIATE_37,
                          OTHER_ATHEROSCLEROSIS	=	RUN1_ENTRY_COVARIATE_38,
                          CARDIAC_PROCEDURE	=	RUN1_ENTRY_COVARIATE_39,
                          CABG_PTCA	=	RUN1_ENTRY_COVARIATE_40,
                          STROKE	=	RUN1_ENTRY_COVARIATE_41,
                          ISCHEMIC_STROKE	=	RUN1_ENTRY_COVARIATE_42,
                          HEMORRHAGIC_STROKE	=	RUN1_ENTRY_COVARIATE_43,
                          TIA	=	RUN1_ENTRY_COVARIATE_44,
                          OTHER_CEREBROVASCULAR	=	RUN1_ENTRY_COVARIATE_45,
                          LATE_CEREBROVASCULAR	=	RUN1_ENTRY_COVARIATE_46,
                          CEREBROVASCULAR_PROCEDURE	=	RUN1_ENTRY_COVARIATE_47,
                          CHF	=	RUN1_ENTRY_COVARIATE_48,
                          PVD	=	RUN1_ENTRY_COVARIATE_49,
                          ATRIAL_FIBRILLATION	=	RUN1_ENTRY_COVARIATE_50,
                          CARDIAC_DYSRHYTHMIA	=	RUN1_ENTRY_COVARIATE_51,
                          CARDIAC_CONDUCTION	=	RUN1_ENTRY_COVARIATE_52,
                          OTHER_CVD	=	RUN1_ENTRY_COVARIATE_53,
                          EDEMA	=	RUN1_ENTRY_COVARIATE_54,
                          COPD	=	RUN1_ENTRY_COVARIATE_55,
                          ASTHMA	=	RUN1_ENTRY_COVARIATE_56,
                          OBS_SLEEP_APNEA	=	RUN1_ENTRY_COVARIATE_57,
                          PNEUMONIA	=	RUN1_ENTRY_COVARIATE_58,
                          RENAL_DYSFUNCTION	=	RUN1_ENTRY_COVARIATE_59,
                          ACUTE_RENAL_DISEASE	=	RUN1_ENTRY_COVARIATE_60,
                          CHRONIC_RENAL_INSUFFICIENCY	=	RUN1_ENTRY_COVARIATE_61,
                          CKD	=	RUN1_ENTRY_COVARIATE_62,
                          CKD_3_4	=	RUN1_ENTRY_COVARIATE_63,
                          HYPERTENSIVE_NEPHROPATHY	=	RUN1_ENTRY_COVARIATE_64,
                          OTHER_RENAL_INSUFFICIENCY	=	RUN1_ENTRY_COVARIATE_65,
                          LIVER_DISEASE	=	RUN1_ENTRY_COVARIATE_66,
                          OSTEOARTHRITIS	=	RUN1_ENTRY_COVARIATE_67,
                          OTHER_ARTHRITIS	=	RUN1_ENTRY_COVARIATE_68,
                          DORSOPATHIES	=	RUN1_ENTRY_COVARIATE_69,
                          FRACTURE	=	RUN1_ENTRY_COVARIATE_70,
                          FALL	=	RUN1_ENTRY_COVARIATE_71,
                          OSTEOPOROSIS	=	RUN1_ENTRY_COVARIATE_72,
                          HYPERTHYROIDISM	=	RUN1_ENTRY_COVARIATE_73,
                          HYPOTHYROIDISM	=	RUN1_ENTRY_COVARIATE_74,
                          OTHER_THYROID	=	RUN1_ENTRY_COVARIATE_75,
                          DEPRESSION	=	RUN1_ENTRY_COVARIATE_76,
                          ANXIETY	=	RUN1_ENTRY_COVARIATE_77,
                          SLEEP_DISORDER	=	RUN1_ENTRY_COVARIATE_78,
                          DEMENTIA	=	RUN1_ENTRY_COVARIATE_79,
                          DELIRIUM	=	RUN1_ENTRY_COVARIATE_80,
                          PSYCHOSIS	=	RUN1_ENTRY_COVARIATE_81,
                          FRAILTY_QUALITATIVE	=	RUN1_ENTRY_COVARIATE_82,
                          FRAILTY_EMPIRICAL_V3	=	RUN1_ENTRY_COVARIATE_83,
                          NON_FRAILTY	=	RUN1_ENTRY_COVARIATE_84,
                          N_ANTIDIABETICS	=	RUN1_ENTRY_COVARIATE_85,
                          NAIVE_NEW_USER	=	RUN1_ENTRY_COVARIATE_86,
                          ACE_INHIBITORS	=	RUN1_ENTRY_COVARIATE_87,
                          ARB	=	RUN1_ENTRY_COVARIATE_88,
                          LOOP_DIURETICS	=	RUN1_ENTRY_COVARIATE_89,
                          OTHER_DIURETICS	=	RUN1_ENTRY_COVARIATE_90,
                          NITRATES	=	RUN1_ENTRY_COVARIATE_91,
                          OTHER_ANTIHYPERTENSIVE	=	RUN1_ENTRY_COVARIATE_92,
                          DIGOXIN	=	RUN1_ENTRY_COVARIATE_93,
                          ANTI_ARRHYTHMICS	=	RUN1_ENTRY_COVARIATE_94,
                          COPD_ASTHMA_DRUG	=	RUN1_ENTRY_COVARIATE_95,
                          STATIN	=	RUN1_ENTRY_COVARIATE_96,
                          OTHER_LIPID_DRUG	=	RUN1_ENTRY_COVARIATE_97,
                          ANTIPLATELET	=	RUN1_ENTRY_COVARIATE_98,
                          ORAL_ANTICOAGULANT	=	RUN1_ENTRY_COVARIATE_99,
                          HEPARIN	=	RUN1_ENTRY_COVARIATE_100,
                          NSAID	=	RUN1_ENTRY_COVARIATE_101,
                          ORAL_CORTICOSTEROID	=	RUN1_ENTRY_COVARIATE_102,
                          BISPHOSPHONATE	=	RUN1_ENTRY_COVARIATE_103,
                          OPIOID	=	RUN1_ENTRY_COVARIATE_104,
                          ANTIDEPRESSANT	=	RUN1_ENTRY_COVARIATE_105,
                          ANTIPSYCHOTIC	=	RUN1_ENTRY_COVARIATE_106,
                          ANTICONVULSANT	=	RUN1_ENTRY_COVARIATE_107,
                          LITHIUM	=	RUN1_ENTRY_COVARIATE_108,
                          BENZODIAZAPINE	=	RUN1_ENTRY_COVARIATE_109,
                          ANXIOLYTIC_HYPNOTIC	=	RUN1_ENTRY_COVARIATE_110,
                          DEMENTIA_DRUG	=	RUN1_ENTRY_COVARIATE_111,
                          PARKINSON_DRUG	=	RUN1_ENTRY_COVARIATE_112,
                          N_DIAGNOSES	=	RUN1_ENTRY_COVARIATE_113,
                          N_DRUG_RX	=	RUN1_ENTRY_COVARIATE_114,
                          HOSPITALIZATION	=	RUN1_ENTRY_COVARIATE_115,
                          ENDOCRINOLOGIST	=	RUN1_ENTRY_COVARIATE_116,
                          INTERNAL	=	RUN1_ENTRY_COVARIATE_117,
                          CARDIOLOGIST	=	RUN1_ENTRY_COVARIATE_118,
                          ELECTROCARDIOGRAM	=	RUN1_ENTRY_COVARIATE_119,
                          GLUCOSE_TESTS	=	RUN1_ENTRY_COVARIATE_120,
                          HOSPITALIZATION_30	=	RUN1_ENTRY_COVARIATE_121,
                          HOSPITALIZATION_31_180	=	RUN1_ENTRY_COVARIATE_122,
                          N_HOSPITALIZATION	=	RUN1_ENTRY_COVARIATE_123,
                          N_HOSPITAL_DAYS	=	RUN1_ENTRY_COVARIATE_124,
                          N_ED	=	RUN1_ENTRY_COVARIATE_125,
                          N_OFFICE	=	RUN1_ENTRY_COVARIATE_126,
                          N_ENDOCRINOLOGIST	=	RUN1_ENTRY_COVARIATE_127,
                          N_INTERNAL	=	RUN1_ENTRY_COVARIATE_128,
                          N_CARDIOLOGIST	=	RUN1_ENTRY_COVARIATE_129,
                          N_ELECTROCARDIOGRAM	=	RUN1_ENTRY_COVARIATE_130,
                          N_HBA1C	=	RUN1_ENTRY_COVARIATE_131,
                          N_GLUCOSE_TEST	=	RUN1_ENTRY_COVARIATE_132,
                          N_LIPID_TEST	=	RUN1_ENTRY_COVARIATE_133,
                          N_CREATININE_TEST	=	RUN1_ENTRY_COVARIATE_134,
                          N_BUN_TEST	=	RUN1_ENTRY_COVARIATE_135,
                          N_MICROALBUMINURIA_TEST	=	RUN1_ENTRY_COVARIATE_136,
                          DM_DRUG_AGI	=	RUN1_ENTRY_COVARIATE_137,
                          DM_DRUG_GLITAZONE	=	RUN1_ENTRY_COVARIATE_138,
                          DM_DRUG_GLP1	=	RUN1_ENTRY_COVARIATE_139,
                          DM_DRUG_INSULIN	=	RUN1_ENTRY_COVARIATE_140,
                          DM_DRUG_MEGLITINIDE	=	RUN1_ENTRY_COVARIATE_141,
                          DM_DRUG_METFORMIN	=	RUN1_ENTRY_COVARIATE_142,
                          PRAMLINTIDE	=	RUN1_ENTRY_COVARIATE_143,
                          GEN1_SU	=	RUN1_ENTRY_COVARIATE_144,
                          MONOTHERAPY_INITIATION	=	RUN1_ENTRY_COVARIATE_145,
                          METFORMIN_DUAL	=	RUN1_ENTRY_COVARIATE_146,
                          AGE	=	RUN1_ENTRY_COVARIATE_147,
                          CEREBROVASCULAR_HEM_STROKE	=	RUN1_ENTRY_COVARIATE_148,
                          BLADDER_STONE	=	RUN1_ENTRY_COVARIATE_149,
                          KIDNEY_STONE	=	RUN1_ENTRY_COVARIATE_150,
                          UTI	=	RUN1_ENTRY_COVARIATE_151,
                          DIPSTICK_URINALYSIS	=	RUN1_ENTRY_COVARIATE_152,
                          NON_DIPSTICK_URINALYSIS	=	RUN1_ENTRY_COVARIATE_153,
                          URINE_FUNCTION_TEST	=	RUN1_ENTRY_COVARIATE_154,
                          CYTOLOGY	=	RUN1_ENTRY_COVARIATE_155,
                          CYSTOSCOPY	=	RUN1_ENTRY_COVARIATE_156,
                          FRAILTY_EMPIRICAL	=	RUN1_ENTRY_COVARIATE_157,
                          CREATININE_TEST	=	RUN1_ENTRY_COVARIATE_158,
                          BUN_TEST	=	RUN1_ENTRY_COVARIATE_159,
                          CRI_NO_CKD	=	RUN1_ENTRY_COVARIATE_160,
                          CKD_1_2	=	RUN1_ENTRY_COVARIATE_161,
                          CKD_3_6	=	RUN1_ENTRY_COVARIATE_162,
                          CONCOMITANT_SGLT2I	=	RUN1_ENTRY_COVARIATE_163,
                          CONCOMITANT_AGI	=	RUN1_ENTRY_COVARIATE_164,
                          CONCOMITANT_GLITAZONE	=	RUN1_ENTRY_COVARIATE_165,
                          CONCOMITANT_GLP1	=	RUN1_ENTRY_COVARIATE_166,
                          CONCOMITANT_INSULIN	=	RUN1_ENTRY_COVARIATE_167,
                          CONCOMITANT_MEGLITINIDE	=	RUN1_ENTRY_COVARIATE_168,
                          CONCOMITANT_METFORMIN	=	RUN1_ENTRY_COVARIATE_169,
                          FRAILTY_QUALITATIVE_V1	=	RUN1_ENTRY_COVARIATE_170,
                          PAST_SGLT2I	=	RUN1_ENTRY_COVARIATE_171,
                          PAST_AGI	=	RUN1_ENTRY_COVARIATE_172,
                          PAST_GLITAZONE	=	RUN1_ENTRY_COVARIATE_173,
                          PAST_GLP1	=	RUN1_ENTRY_COVARIATE_174,
                          PAST_INSULIN	=	RUN1_ENTRY_COVARIATE_175,
                          PAST_MEGLITINIDE	=	RUN1_ENTRY_COVARIATE_176,
                          PAST_METFORMIN	=	RUN1_ENTRY_COVARIATE_177,
                          CALENDAR_TIME_DAY	=	RUN1_ENTRY_COVARIATE_178,
                          BLADDER_KIDNEY_STONE	=	RUN1_ENTRY_COVARIATE_179,
                          PERIPHERAL_GANGRENE_OSTEOMYELITIS	=	RUN1_ENTRY_COVARIATE_180,
                          AGE_DECILE	=	RUN1_ENTRY_COVARIATE_181,
                          ALCOHOL_DRUG	=	RUN1_ENTRY_COVARIATE_182,
                          DM_EYE	=	RUN1_ENTRY_COVARIATE_183,
                          COMPOSITE_CVD	=	RUN1_ENTRY_COVARIATE_184,
                          COMPOSITE_CARDIAC_PROCEDURE	=	RUN1_ENTRY_COVARIATE_185,
                          THYROID	=	RUN1_ENTRY_COVARIATE_186,
                          DELIRIUM_PSYCHOSIS	=	RUN1_ENTRY_COVARIATE_187,
                          MEGLITINIDE	=	RUN1_ENTRY_COVARIATE_188,
                          AGI	=	RUN1_ENTRY_COVARIATE_189,
                          GLAUCOMA_CATARACTS	=	RUN1_ENTRY_COVARIATE_190,
                          IMAGING	=	RUN1_ENTRY_COVARIATE_191,
                          CELLULITIS_ABCESS_TOE	=	RUN1_ENTRY_COVARIATE_192,
                          FOOT_ULCER	=	RUN1_ENTRY_COVARIATE_193,
                          ENTRESTO	=	RUN1_ENTRY_COVARIATE_194,
                          ENDOCRINOLOGIST_30	=	RUN1_ENTRY_COVARIATE_195,
                          ENDOCRINOLOGIST_31_180	=	RUN1_ENTRY_COVARIATE_196,
                          CARDIOLOGIST_30	=	RUN1_ENTRY_COVARIATE_197,
                          CARDIOLOGIST_31_180	=	RUN1_ENTRY_COVARIATE_198,
                          INTERNAL_30	=	RUN1_ENTRY_COVARIATE_199,
                          INTERNAL_31_180	=	RUN1_ENTRY_COVARIATE_200,
                          DIALYSIS	=	RUN1_ENTRY_COVARIATE_201,
                          CCI_180	=	RUN1_ENTRY_COVARIATE_202,
                          CKD_3_6_DIALYSIS	=	RUN1_ENTRY_COVARIATE_203,
                          BASELINE_CV	=	RUN1_ENTRY_COVARIATE_204,
                          THIAZIDE	=	RUN1_ENTRY_COVARIATE_205,
                          BETA_BLOCKER	=	RUN1_ENTRY_COVARIATE_206,
                          CA_CHANNEL_BLOCKER	=	RUN1_ENTRY_COVARIATE_207)

# Convert all true/false logicals to 1/0 numeric
savor[,c(9:85, 88, 90:116, 119:126, 141:150, 152:160, 162:173, 175:181, 183:184, 186:205, 207:211)] <- lapply(savor[,c(9:85, 88, 90:116, 119:126, 141:150, 152:160, 162:173, 175:181, 183:184, 186:205, 207:211)], function(x) as.numeric(as.logical(x)))

# Create normalized versions of MED_INC and MED_INC_65 to [0, 1] for use in analysis
savor <- savor %>% mutate(MED_INC_STD = (MED_INC - min(MED_INC))/(max(MED_INC) - min(MED_INC)))
savor <- savor %>% mutate(MED_INC_65_STD =
                            (MED_INC_65 - min(MED_INC_65))/(max(MED_INC_65) - min(MED_INC_65)))

# Add state variable based on ZCTA
savor <- savor %>% mutate(
  STATE_FIPS = case_when(
    as.numeric(ZCTA_USE) >= 35000 & as.numeric(ZCTA_USE) < 37000 ~ "01",
    as.numeric(ZCTA_USE) >= 99500 ~ "02",
    as.numeric(ZCTA_USE) >= 85000 & as.numeric(ZCTA_USE) < 87000 ~ "04",
    as.numeric(ZCTA_USE) >= 71600 & as.numeric(ZCTA_USE) < 73000 ~ "05",
    as.numeric(ZCTA_USE) >= 90000 & as.numeric(ZCTA_USE) < 96200 ~ "06",
    as.numeric(ZCTA_USE) >= 80000 & as.numeric(ZCTA_USE) < 82000 ~ "08",
    as.numeric(ZCTA_USE) >= 6000 & as.numeric(ZCTA_USE) < 7000 ~ "09",
    as.numeric(ZCTA_USE) >= 19700 & as.numeric(ZCTA_USE) < 20000 ~ "10",
    as.numeric(ZCTA_USE) >= 20000 & as.numeric(ZCTA_USE) < 20600 ~ "11",
    as.numeric(ZCTA_USE) >= 32000 & as.numeric(ZCTA_USE) < 35000 ~ "12",
    as.numeric(ZCTA_USE) >= 30000 & as.numeric(ZCTA_USE) < 32000 ~ "13",
    as.numeric(ZCTA_USE) >= 96700 & as.numeric(ZCTA_USE) < 96900 ~ "15",
    as.numeric(ZCTA_USE) >= 83200 & as.numeric(ZCTA_USE) < 84000 ~ "16",
    as.numeric(ZCTA_USE) >= 60000 & as.numeric(ZCTA_USE) < 63000 ~ "17",
    as.numeric(ZCTA_USE) >= 46000 & as.numeric(ZCTA_USE) < 48000 ~ "18",
    as.numeric(ZCTA_USE) >= 50000 & as.numeric(ZCTA_USE) < 53000 ~ "19",
    as.numeric(ZCTA_USE) >= 66000 & as.numeric(ZCTA_USE) < 68000 ~ "20",
    as.numeric(ZCTA_USE) >= 40000 & as.numeric(ZCTA_USE) < 43000 ~ "21",
    as.numeric(ZCTA_USE) >= 70000 & as.numeric(ZCTA_USE) < 71600 ~ "22",
    as.numeric(ZCTA_USE) >= 3900 & as.numeric(ZCTA_USE) < 5000 ~ "23",
    as.numeric(ZCTA_USE) >= 20600 & as.numeric(ZCTA_USE) < 22000 ~ "24",
    as.numeric(ZCTA_USE) >= 1000 & as.numeric(ZCTA_USE) < 2800 ~ "25",
    as.numeric(ZCTA_USE) >= 48000 & as.numeric(ZCTA_USE) < 50000 ~ "26",
    as.numeric(ZCTA_USE) >= 55000 & as.numeric(ZCTA_USE) < 56800 ~ "27",
    as.numeric(ZCTA_USE) >= 38600 & as.numeric(ZCTA_USE) < 40000 ~ "28",
    as.numeric(ZCTA_USE) >= 63000 & as.numeric(ZCTA_USE) < 66000 ~ "29",
    as.numeric(ZCTA_USE) >= 59000 & as.numeric(ZCTA_USE) < 60000 ~ "30",
    as.numeric(ZCTA_USE) >= 68000 & as.numeric(ZCTA_USE) < 70000 ~ "31",
    as.numeric(ZCTA_USE) >= 88900 & as.numeric(ZCTA_USE) < 90000 ~ "32",
    as.numeric(ZCTA_USE) >= 3000 & as.numeric(ZCTA_USE) < 3900 ~ "33",
    as.numeric(ZCTA_USE) >= 7000 & as.numeric(ZCTA_USE) < 9000 ~ "34",
    as.numeric(ZCTA_USE) >= 87000 & as.numeric(ZCTA_USE) < 88500 ~ "35",
    as.numeric(ZCTA_USE) >= 10000 & as.numeric(ZCTA_USE) < 15000 ~ "36",
    as.numeric(ZCTA_USE) >= 27000 & as.numeric(ZCTA_USE) < 29000 ~ "37",
    as.numeric(ZCTA_USE) >= 58000 & as.numeric(ZCTA_USE) < 59000 ~ "38",
    as.numeric(ZCTA_USE) >= 43000 & as.numeric(ZCTA_USE) < 46000 ~ "39",
    as.numeric(ZCTA_USE) >= 73000 & as.numeric(ZCTA_USE) < 75000 ~ "40",
    as.numeric(ZCTA_USE) >= 97000 & as.numeric(ZCTA_USE) < 98000 ~ "41",
    as.numeric(ZCTA_USE) >= 15000 & as.numeric(ZCTA_USE) < 19700 ~ "42",
    as.numeric(ZCTA_USE) >= 2800 & as.numeric(ZCTA_USE) < 3000 ~ "44",
    as.numeric(ZCTA_USE) >= 29000 & as.numeric(ZCTA_USE) < 30000 ~ "45",
    as.numeric(ZCTA_USE) >= 57000 & as.numeric(ZCTA_USE) < 58000 ~ "46",
    as.numeric(ZCTA_USE) >= 37000 & as.numeric(ZCTA_USE) < 38600 ~ "47",
    as.numeric(ZCTA_USE) >= 75000 & as.numeric(ZCTA_USE) < 80000 ~ "48",
    as.numeric(ZCTA_USE) >= 84000 & as.numeric(ZCTA_USE) < 85000 ~ "49",
    as.numeric(ZCTA_USE) >= 5000 & as.numeric(ZCTA_USE) < 6000 ~ "50",
    as.numeric(ZCTA_USE) >= 22000 & as.numeric(ZCTA_USE) < 24700 ~ "51",
    as.numeric(ZCTA_USE) >= 98000 & as.numeric(ZCTA_USE) < 99500 ~ "53",
    as.numeric(ZCTA_USE) >= 24700 & as.numeric(ZCTA_USE) < 27000 ~ "54",
    as.numeric(ZCTA_USE) >= 53000 & as.numeric(ZCTA_USE) < 55000 ~ "55",
    TRUE ~ "56"
  )
)

savor <- savor %>% mutate(
  STATE = case_when(
    as.numeric(ZCTA_USE) >= 35000 & as.numeric(ZCTA_USE) < 37000 ~ "AL",
    as.numeric(ZCTA_USE) >= 99500 ~ "AK",
    as.numeric(ZCTA_USE) >= 85000 & as.numeric(ZCTA_USE) < 87000 ~ "AZ",
    as.numeric(ZCTA_USE) >= 71600 & as.numeric(ZCTA_USE) < 73000 ~ "AR",
    as.numeric(ZCTA_USE) >= 90000 & as.numeric(ZCTA_USE) < 96200 ~ "CA",
    as.numeric(ZCTA_USE) >= 80000 & as.numeric(ZCTA_USE) < 82000 ~ "CO",
    as.numeric(ZCTA_USE) >= 6000 & as.numeric(ZCTA_USE) < 7000 ~ "CT",
    as.numeric(ZCTA_USE) >= 19700 & as.numeric(ZCTA_USE) < 20000 ~ "DE",
    (as.numeric(ZCTA_USE) >= 20000 & as.numeric(ZCTA_USE) < 20100) | (as.numeric(ZCTA_USE) >= 20200 & as.numeric(ZCTA_USE) < 20600) ~ "DC",
    as.numeric(ZCTA_USE) >= 32000 & as.numeric(ZCTA_USE) < 35000 ~ "FL",
    (as.numeric(ZCTA_USE) >= 30000 & as.numeric(ZCTA_USE) < 32000) | (as.numeric(ZCTA_USE) >= 39800 & as.numeric(ZCTA_USE) < 40000) ~ "GA",
    as.numeric(ZCTA_USE) >= 96700 & as.numeric(ZCTA_USE) < 96900 ~ "HI",
    as.numeric(ZCTA_USE) >= 83200 & as.numeric(ZCTA_USE) < 84000 ~ "ID",
    as.numeric(ZCTA_USE) >= 60000 & as.numeric(ZCTA_USE) < 63000 ~ "IL",
    as.numeric(ZCTA_USE) >= 46000 & as.numeric(ZCTA_USE) < 48000 ~ "IN",
    as.numeric(ZCTA_USE) >= 50000 & as.numeric(ZCTA_USE) < 53000 ~ "IA",
    as.numeric(ZCTA_USE) >= 66000 & as.numeric(ZCTA_USE) < 68000 ~ "KS",
    as.numeric(ZCTA_USE) >= 40000 & as.numeric(ZCTA_USE) < 43000 ~ "KY",
    as.numeric(ZCTA_USE) >= 70000 & as.numeric(ZCTA_USE) < 71600 ~ "LA",
    as.numeric(ZCTA_USE) >= 3900 & as.numeric(ZCTA_USE) < 5000 ~ "ME",
    as.numeric(ZCTA_USE) >= 20600 & as.numeric(ZCTA_USE) < 22000 ~ "MD",
    as.numeric(ZCTA_USE) >= 1000 & as.numeric(ZCTA_USE) < 2800 ~ "MA",
    as.numeric(ZCTA_USE) >= 48000 & as.numeric(ZCTA_USE) < 50000 ~ "MI",
    as.numeric(ZCTA_USE) >= 55000 & as.numeric(ZCTA_USE) < 56800 ~ "MN",
    (as.numeric(ZCTA_USE) >= 38600 & as.numeric(ZCTA_USE) < 39800) ~ "MS",
    as.numeric(ZCTA_USE) >= 63000 & as.numeric(ZCTA_USE) < 66000 ~ "MO",
    as.numeric(ZCTA_USE) >= 59000 & as.numeric(ZCTA_USE) < 60000 ~ "MT",
    as.numeric(ZCTA_USE) >= 68000 & as.numeric(ZCTA_USE) < 70000 ~ "NE",
    as.numeric(ZCTA_USE) >= 88900 & as.numeric(ZCTA_USE) < 90000 ~ "NV",
    as.numeric(ZCTA_USE) >= 3000 & as.numeric(ZCTA_USE) < 3900 ~ "NH",
    as.numeric(ZCTA_USE) >= 7000 & as.numeric(ZCTA_USE) < 9000 ~ "NJ",
    as.numeric(ZCTA_USE) >= 87000 & as.numeric(ZCTA_USE) < 88500 ~ "NM",
    as.numeric(ZCTA_USE) >= 10000 & as.numeric(ZCTA_USE) < 15000 ~ "NY",
    as.numeric(ZCTA_USE) >= 27000 & as.numeric(ZCTA_USE) < 29000 ~ "NC",
    as.numeric(ZCTA_USE) >= 58000 & as.numeric(ZCTA_USE) < 59000 ~ "ND",
    as.numeric(ZCTA_USE) >= 43000 & as.numeric(ZCTA_USE) < 46000 ~ "OH",
    as.numeric(ZCTA_USE) >= 73000 & as.numeric(ZCTA_USE) < 75000 ~ "OK",
    as.numeric(ZCTA_USE) >= 97000 & as.numeric(ZCTA_USE) < 98000 ~ "OR",
    as.numeric(ZCTA_USE) >= 15000 & as.numeric(ZCTA_USE) < 19700 ~ "PA",
    as.numeric(ZCTA_USE) >= 2800 & as.numeric(ZCTA_USE) < 3000 ~ "RI",
    as.numeric(ZCTA_USE) >= 29000 & as.numeric(ZCTA_USE) < 30000 ~ "SC",
    as.numeric(ZCTA_USE) >= 57000 & as.numeric(ZCTA_USE) < 58000 ~ "SD",
    as.numeric(ZCTA_USE) >= 37000 & as.numeric(ZCTA_USE) < 38600 ~ "TN",
    as.numeric(ZCTA_USE) >= 75000 & as.numeric(ZCTA_USE) < 80000 ~ "TX",
    as.numeric(ZCTA_USE) >= 84000 & as.numeric(ZCTA_USE) < 85000 ~ "UT",
    as.numeric(ZCTA_USE) >= 5000 & as.numeric(ZCTA_USE) < 6000 ~ "VT",
    (as.numeric(ZCTA_USE) >= 22000 & as.numeric(ZCTA_USE) < 24700) | (as.numeric(ZCTA_USE) >= 20100 & as.numeric(ZCTA_USE) < 20200) ~ "VA",
    as.numeric(ZCTA_USE) >= 98000 & as.numeric(ZCTA_USE) < 99500 ~ "WA",
    as.numeric(ZCTA_USE) >= 24700 & as.numeric(ZCTA_USE) < 27000 ~ "WV",
    as.numeric(ZCTA_USE) >= 53000 & as.numeric(ZCTA_USE) < 55000 ~ "WI",
    TRUE ~ "WY"
  )
)

# Add region variable based on ZCTA
savor <- savor %>% mutate(
  REGION_ZCTA = case_when(
    (STATE == "ME" |
       STATE == "NH" |
       STATE == "VT" |
       STATE == "MA" |
       STATE == "RI" |
       STATE == "CT" |
       STATE == "NY" |
       STATE == "NJ" |
       STATE == "PA") ~ "0",
    (STATE == "DE" |
       STATE == "MD" |
       STATE == "DC" |
       STATE == "WV" |
       STATE == "VA" |
       STATE == "NC" |
       STATE == "SC" |
       STATE == "GA" |
       STATE == "FL" |
       STATE == "AL" |
       STATE == "TN" |
       STATE == "KY" |
       STATE == "MS" |
       STATE == "AR" |
       STATE == "LA" |
       STATE == "OK" |
       STATE == "TX") ~ "1",
    (STATE == "OH" |
       STATE == "MI" |
       STATE == "IN" |
       STATE == "IL" |
       STATE == "WI" |
       STATE == "MN" |
       STATE == "IA" |
       STATE == "MO" |
       STATE == "KS" |
       STATE == "NE" |
       STATE == "SD" |
       STATE == "ND") ~ "2",
    TRUE ~ "3"
  )
)

# Descriptive Statistics
# Key non-zipcode variables
summary(savor$AGE)
d_label <- gTree("D", children = gList(textGrob("D", x = 0, y = 0.8,
                                                just = c("left", "top"),
                                                gp = gpar(fontsize = 32, col =  "black"))))
savor_age <- savor %>%
  ggplot() + 
  geom_histogram(aes(AGE), binwidth = 1, color = "black", fill = "white") + 
  geom_vline(xintercept = 74.64, color = "red") +
  scale_y_continuous(labels = label_comma()) +
  coord_cartesian(xlim = c(65, 115)) +
  labs(title = "Age Distribution\nSaxagliptin", x = "Age (years)", y = "Count") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))
savor_age_grob <- gtable_add_grob(ggplot_gtable(ggplot_build(savor_age)), d_label, t = 1, l = 4, b = 6)
ggsave("savor_age_dist.pdf", savor_age_grob)

table(savor$RACE)
e_label <- gTree("E", children = gList(textGrob("E", x = 0, y = 0.8,
                                                just = c("left", "top"),
                                                gp = gpar(fontsize = 32, col =  "black"))))
savor_race <- savor %>%
  ggplot() + 
  geom_bar(aes(RACE), color = "black", fill = "white") +
  scale_x_continuous(breaks = 0:5,
                     labels = c("White","Black","Asian","Hispanic","Native\nAmerican","Other")) +
  scale_y_continuous(labels = label_comma()) +
  labs(title = "Race Distribution\nSaxagliptin", x = "Race", y = "Count") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(hjust = 0.5))
savor_race_grob <- gtable_add_grob(ggplot_gtable(ggplot_build(savor_race)), e_label, t = 1, l = 4, b = 6)
ggsave("savor_race_dist.pdf", savor_race_grob)

table(savor$REGION_ZCTA)
f_label <- gTree("F", children = gList(textGrob("F", x = 0, y = 0.8,
                                                just = c("left", "top"),
                                                gp = gpar(fontsize = 32, col =  "black"))))
savor_region <- savor %>%
  ggplot() + 
  geom_bar(aes(as.numeric(REGION_ZCTA)), color = "black", fill = "white") +
  scale_x_continuous(breaks = 0:3,
                     labels = c("Northeast",
                                "South",
                                "Midwest",
                                "West")) +
  scale_y_continuous(labels = label_comma()) +
  labs(title = "Region Distribution\nSaxagliptin", x = "Region", y = "Count") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(hjust = 0.5))
savor_region_grob <- gtable_add_grob(ggplot_gtable(ggplot_build(savor_region)), f_label, t = 1, l = 4, b = 6)
ggsave("savor_region_dist.pdf", savor_region_grob)

# Histogram of individuals per zip code
savor_zip_frequency <- savor %>% group_by(ZCTA_USE) %>% dplyr::summarize(count=n())

summary(savor_zip_frequency$count) 
table(savor_zip_frequency$count==1)
table(savor_zip_frequency$count==2)
table(savor_zip_frequency$count>=5)
table(savor_zip_frequency$count>=10)
table(savor_zip_frequency$count>=15)

ggplot(savor_zip_frequency) + 
  geom_histogram(aes(count), binwidth = 10, color = "black", fill = "white") + 
  geom_vline(xintercept = 29.35, color = "red") +
  labs(title = "ZCTA Population Distribution (Saxagliptin)", x = "Number of Individuals in ZCTA", y = "Count") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("savor_zcta_frequency.pdf")

filter(savor_zip_frequency, count>=5) %>%
  ggplot() + 
  geom_histogram(aes(count), binwidth = 10, color = "black", fill = "white") + 
  geom_vline(xintercept = 23.65, color = "red") +
  labs(title = "ZCTA Population Distribution (Saxagliptin)\n5+ in ZCTA", x = "Number of Individuals in ZCTA", y = "Count") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("savor_zcta_frequency_5.pdf")

filter(savor_zip_frequency, count>=10) %>%
  ggplot() + 
  geom_histogram(aes(count), binwidth = 10, color = "black", fill = "white") + 
  geom_vline(xintercept = 30.72, color = "red") +
  labs(title = "ZCTA Population Distribution (Saxagliptin)\n10+ in ZCTA", x = "Number of Individuals in ZCTA", y = "Count") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("savor_zcta_frequency_10.pdf")

filter(savor_zip_frequency, count>=15) %>%
  ggplot() + 
  geom_histogram(aes(count), binwidth = 10, color = "black", fill = "white") + 
  geom_vline(xintercept = 36.47, color = "red") +
  labs(title = "ZCTA Population Distribution (Saxagliptin)\n15+ in ZCTA", x = "Number of Individuals in ZCTA", y = "Count") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("savor_zcta_frequency_15.pdf")

# Add zip code counts to individuals
savor <- merge(savor, savor_zip_frequency, by=c("ZCTA_USE"))

summary(savor$count)
table(savor$count==1)
table(savor$count==2)
table(savor$count>=5)
table(savor$count>=10)
table(savor$count>=15)

ggplot(savor) + 
  geom_histogram(aes(count), binwidth = 10, color = "black", fill = "white") + 
  geom_vline(xintercept = 82.63, color = "red") +
  labs(title = "Individuals by ZCTA Population Distribution\n(Saxagliptin)", x = "Number of Individuals in ZCTA", y = "Count") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("savor_individuals_by_zcta_pop.pdf")

# Zip code exposure distribution: 5 or more in zip code
savor_exp_dist_5 <- savor %>%
  filter(count>=5) %>%
  group_by(ZCTA_USE) %>%
  dplyr::summarize(prop_exp = mean(EXPOSURE))

summary(savor_exp_dist_5$prop_exp) # Mean exposure proportion = 0.044, Median = 0.025

# Re-add zip code counts of individuals
savor_exp_dist_5 <- merge(savor[,c(1,ncol(savor))], savor_exp_dist_5, by = 'ZCTA_USE') %>%
  distinct()

b_label <- gTree("B", children = gList(textGrob("B", x = 0, y = 0.8,
                                                just = c("left", "top"),
                                                gp = gpar(fontsize = 32, col =  "black"))))

savor_exp_dist_5_plot <- savor_exp_dist_5 %>%
  ggplot() + 
  geom_histogram(aes(prop_exp), binwidth = 0.05, color = "black", fill = "white") + 
  geom_vline(xintercept = 0.044, color = "red") +
  scale_y_continuous(labels = label_comma()) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(title = "Saxagliptin Distribution\n5+ in ZCTA", x = "Proportion Saxagliptin", y = "Count") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))
savor_exp_dist_5_grob <- gtable_add_grob(ggplot_gtable(ggplot_build(savor_exp_dist_5_plot)), b_label, t = 1, l = 4, b = 6)
ggsave("savor_zcta_exp_dist_5.pdf", savor_exp_dist_5_grob)

# Zip code exposure distribution: 10 or more in zip code
savor_exp_dist_10 <- savor %>%
  filter(count>=10) %>%
  group_by(ZCTA_USE) %>%
  dplyr::summarize(prop_exp = mean(EXPOSURE))

summary(savor_exp_dist_10$prop_exp) # Mean exposure proportion = 0.045, Median = 0.035

savor_exp_dist_10 %>%
  ggplot() + 
  geom_histogram(aes(prop_exp), binwidth = 0.05, color = "black", fill = "white") + 
  geom_vline(xintercept = 0.045, color = "red") +
  coord_cartesian(xlim = c(0, 1)) +
  labs(title = "Saxagliptin Distribution\n10+ in ZCTA", x = "Proportion Saxagliptin", y = "Count") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("savor_zcta_exp_dist_10.pdf")

# Zip code exposure distribution: 15 or more in zip code
savor_exp_dist_15 <- savor %>%
  filter(count>=15) %>%
  group_by(ZCTA_USE) %>%
  dplyr::summarize(prop_exp = mean(EXPOSURE))

summary(savor_exp_dist_15$prop_exp) # Mean exposure proportion = 0.046, Median = 0.038

savor_exp_dist_15 %>%
  ggplot() + 
  geom_histogram(aes(prop_exp), binwidth = 0.05, color = "black", fill = "white") + 
  geom_vline(xintercept = 0.046, color = "red") +
  coord_cartesian(xlim = c(0, 1)) +
  labs(title = "Saxagliptin Distribution\n15+ in ZCTA", x = "Proportion Saxagliptin", y = "Count") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("savor_zcta_exp_dist_15.pdf")

# Calculate zipcode saxagliptin prescribing proportion for all zipcodes
savor_exp_dist <- savor %>%
  group_by(ZCTA_USE) %>%
  dplyr::summarize(prop_exp = mean(EXPOSURE))

# Add zipcode saxagliptin prescribing proportion as new variable for all individuals
savor <- savor %>% merge(savor_exp_dist, by = "ZCTA_USE")

# Descriptive statistics of zipcode saxagliptin exposure proportion
summary(savor$prop_exp) # Mean exposure proportion = 0.047, Median = 0.038

savor_exp_dist_5_plot <- filter(savor, count > 4) %>%
  ggplot() + 
  geom_histogram(aes(prop_exp), binwidth = 0.05, color = "black", fill = "white") + 
  geom_vline(xintercept = 0.047, color = "red") +
  scale_y_continuous(labels = label_comma()) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(title = "Distribution of Individuals by ZCTA Saxagliptin Dispensing Proportion\n5+ in ZCTA", x = "Proportion Saxagliptin in ZCTA", y = "Count of Individuals") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(hjust = 0.5))
savor_exp_dist_5_grob <- gtable_add_grob(ggplot_gtable(ggplot_build(savor_exp_dist_5_plot)), b_label, t = 1, l = 4, b = 6)
ggsave("savor_exp_dist_5.pdf", savor_exp_dist_5_grob)

# Table 1 of characteristics for all individuals in cohort
savor_factor_vars <- c("AGE_CAT",
                       "SEX",
                       "RACE",
                       "REGION",
                       "OBESITY",
                       "OVERWEIGHT",
                       "SMOKE",
                       "ALCOHOL",
                       "DRUG",
                       "DM_RETINOPATHY",
                       "DM_OPHTHALMIC",
                       "RETINAL_DETACH",
                       "RETINAL_LASER",
                       "DM_NEUROPATHY2",
                       "DM_NEUROPATHY3",
                       "HYPOGLYCEMIA",
                       "HYPERGLYCEMIA",
                       "ELECTROLYTE_ACID_BASE",
                       "DM_KETOACIDOSIS",
                       "HONK",
                       "DM_PERIPHERAL",
                       "DM_FOOT",
                       "GANGRENE",
                       "LOWER_AMPUTATION",
                       "OSTEOMYELITIS",
                       "SKIN_INFECTION",
                       "ED",
                       "DM_COMPLICATION",
                       "DM_NO_COMPLICATION",
                       "HYPERTENSION",
                       "HYPERLIPIDEMIA",
                       "IHD",
                       "ACUTE_MI",
                       "UNSTABLE_ANGINA",
                       "OLD_MI",
                       "STABLE_ANGINA",
                       "CORONARY_ATHEROSCLEROSIS",
                       "OTHER_ATHEROSCLEROSIS",
                       "CARDIAC_PROCEDURE",
                       "CABG_PTCA",
                       "STROKE",
                       "ISCHEMIC_STROKE",
                       "HEMORRHAGIC_STROKE",
                       "TIA",
                       "OTHER_CEREBROVASCULAR",
                       "LATE_CEREBROVASCULAR",
                       "CEREBROVASCULAR_PROCEDURE",
                       "CHF",
                       "PVD",
                       "ATRIAL_FIBRILLATION",
                       "CARDIAC_DYSRHYTHMIA",
                       "CARDIAC_CONDUCTION",
                       "OTHER_CVD",
                       "EDEMA",
                       "COPD",
                       "ASTHMA",
                       "OBS_SLEEP_APNEA",
                       "PNEUMONIA",
                       "RENAL_DYSFUNCTION",
                       "ACUTE_RENAL_DISEASE",
                       "CHRONIC_RENAL_INSUFFICIENCY",
                       "CKD",
                       "CKD_3_4",
                       "HYPERTENSIVE_NEPHROPATHY",
                       "OTHER_RENAL_INSUFFICIENCY",
                       "LIVER_DISEASE",
                       "OSTEOARTHRITIS",
                       "OTHER_ARTHRITIS",
                       "DORSOPATHIES",
                       "FRACTURE",
                       "FALL",
                       "OSTEOPOROSIS",
                       "HYPERTHYROIDISM",
                       "HYPOTHYROIDISM",
                       "OTHER_THYROID",
                       "DEPRESSION",
                       "ANXIETY",
                       "SLEEP_DISORDER",
                       "DEMENTIA",
                       "DELIRIUM",
                       "PSYCHOSIS",
                       "FRAILTY_QUALITATIVE",
                       "FRAILTY_EMPIRICAL_V3",
                       "NON_FRAILTY",
                       "NAIVE_NEW_USER",
                       "ACE_INHIBITORS",
                       "ARB",
                       "LOOP_DIURETICS",
                       "OTHER_DIURETICS",
                       "NITRATES",
                       "OTHER_ANTIHYPERTENSIVE",
                       "DIGOXIN",
                       "ANTI_ARRHYTHMICS",
                       "COPD_ASTHMA_DRUG",
                       "STATIN",
                       "OTHER_LIPID_DRUG",
                       "ANTIPLATELET",
                       "ORAL_ANTICOAGULANT",
                       "HEPARIN",
                       "NSAID",
                       "ORAL_CORTICOSTEROID",
                       "BISPHOSPHONATE",
                       "OPIOID",
                       "ANTIDEPRESSANT",
                       "ANTIPSYCHOTIC",
                       "ANTICONVULSANT",
                       "LITHIUM",
                       "BENZODIAZAPINE",
                       "ANXIOLYTIC_HYPNOTIC",
                       "DEMENTIA_DRUG",
                       "PARKINSON_DRUG",
                       "HOSPITALIZATION",
                       "ENDOCRINOLOGIST",
                       "INTERNAL",
                       "CARDIOLOGIST",
                       "ELECTROCARDIOGRAM",
                       "GLUCOSE_TESTS",
                       "HOSPITALIZATION_30",
                       "HOSPITALIZATION_31_180",
                       "DM_DRUG_AGI",
                       "DM_DRUG_GLITAZONE",
                       "DM_DRUG_GLP1",
                       "DM_DRUG_INSULIN",
                       "DM_DRUG_MEGLITINIDE",
                       "DM_DRUG_METFORMIN",
                       "PRAMLINTIDE",
                       "GEN1_SU",
                       "MONOTHERAPY_INITIATION",
                       "METFORMIN_DUAL",
                       "CEREBROVASCULAR_HEM_STROKE",
                       "BLADDER_STONE",
                       "KIDNEY_STONE",
                       "UTI",
                       "DIPSTICK_URINALYSIS",
                       "NON_DIPSTICK_URINALYSIS",
                       "URINE_FUNCTION_TEST",
                       "CYTOLOGY",
                       "CYSTOSCOPY",
                       "CREATININE_TEST",
                       "BUN_TEST",
                       "CRI_NO_CKD",
                       "CKD_1_2",
                       "CKD_3_6",
                       "CONCOMITANT_SGLT2I",
                       "CONCOMITANT_AGI",
                       "CONCOMITANT_GLITAZONE",
                       "CONCOMITANT_GLP1",
                       "CONCOMITANT_INSULIN",
                       "CONCOMITANT_MEGLITINIDE",
                       "CONCOMITANT_METFORMIN",
                       "PAST_SGLT2I",
                       "PAST_AGI",
                       "PAST_GLITAZONE",
                       "PAST_GLP1",
                       "PAST_INSULIN",
                       "PAST_MEGLITINIDE",
                       "PAST_METFORMIN",
                       "BLADDER_KIDNEY_STONE",
                       "PERIPHERAL_GANGRENE_OSTEOMYELITIS",
                       "AGE_DECILE",
                       "ALCOHOL_DRUG",
                       "DM_EYE",
                       "COMPOSITE_CVD",
                       "COMPOSITE_CARDIAC_PROCEDURE",
                       "THYROID",
                       "DELIRIUM_PSYCHOSIS",
                       "MEGLITINIDE",
                       "AGI",
                       "GLAUCOMA_CATARACTS",
                       "IMAGING",
                       "CELLULITIS_ABCESS_TOE",
                       "FOOT_ULCER",
                       "ENTRESTO",
                       "ENDOCRINOLOGIST_30",
                       "ENDOCRINOLOGIST_31_180",
                       "CARDIOLOGIST_30",
                       "CARDIOLOGIST_31_180",
                       "INTERNAL_30",
                       "INTERNAL_31_180",
                       "DIALYSIS",
                       "CCI_180",
                       "CKD_3_6_DIALYSIS",
                       "BASELINE_CV",
                       "THIAZIDE",
                       "BETA_BLOCKER",
                       "CA_CHANNEL_BLOCKER")

savor_vars <- c("AGE_CAT",
                "SEX",
                "RACE",
                "REGION",
                "OBESITY",
                "OVERWEIGHT",
                "SMOKE",
                "ALCOHOL",
                "DRUG",
                "DM_RETINOPATHY",
                "DM_OPHTHALMIC",
                "RETINAL_DETACH",
                "RETINAL_LASER",
                "DM_NEUROPATHY2",
                "DM_NEUROPATHY3",
                "HYPOGLYCEMIA",
                "HYPERGLYCEMIA",
                "ELECTROLYTE_ACID_BASE",
                "DM_KETOACIDOSIS",
                "HONK",
                "DM_PERIPHERAL",
                "DM_FOOT",
                "GANGRENE",
                "LOWER_AMPUTATION",
                "OSTEOMYELITIS",
                "SKIN_INFECTION",
                "ED",
                "DM_COMPLICATION",
                "DM_NO_COMPLICATION",
                "HYPERTENSION",
                "HYPERLIPIDEMIA",
                "IHD",
                "ACUTE_MI",
                "UNSTABLE_ANGINA",
                "OLD_MI",
                "STABLE_ANGINA",
                "CORONARY_ATHEROSCLEROSIS",
                "OTHER_ATHEROSCLEROSIS",
                "CARDIAC_PROCEDURE",
                "CABG_PTCA",
                "STROKE",
                "ISCHEMIC_STROKE",
                "HEMORRHAGIC_STROKE",
                "TIA",
                "OTHER_CEREBROVASCULAR",
                "LATE_CEREBROVASCULAR",
                "CEREBROVASCULAR_PROCEDURE",
                "CHF",
                "PVD",
                "ATRIAL_FIBRILLATION",
                "CARDIAC_DYSRHYTHMIA",
                "CARDIAC_CONDUCTION",
                "OTHER_CVD",
                "EDEMA",
                "COPD",
                "ASTHMA",
                "OBS_SLEEP_APNEA",
                "PNEUMONIA",
                "RENAL_DYSFUNCTION",
                "ACUTE_RENAL_DISEASE",
                "CHRONIC_RENAL_INSUFFICIENCY",
                "CKD",
                "CKD_3_4",
                "HYPERTENSIVE_NEPHROPATHY",
                "OTHER_RENAL_INSUFFICIENCY",
                "LIVER_DISEASE",
                "OSTEOARTHRITIS",
                "OTHER_ARTHRITIS",
                "DORSOPATHIES",
                "FRACTURE",
                "FALL",
                "OSTEOPOROSIS",
                "HYPERTHYROIDISM",
                "HYPOTHYROIDISM",
                "OTHER_THYROID",
                "DEPRESSION",
                "ANXIETY",
                "SLEEP_DISORDER",
                "DEMENTIA",
                "DELIRIUM",
                "PSYCHOSIS",
                "FRAILTY_QUALITATIVE",
                "FRAILTY_EMPIRICAL_V3",
                "NON_FRAILTY",
                "N_ANTIDIABETICS",
                "NAIVE_NEW_USER",
                "ACE_INHIBITORS",
                "ARB",
                "LOOP_DIURETICS",
                "OTHER_DIURETICS",
                "NITRATES",
                "OTHER_ANTIHYPERTENSIVE",
                "DIGOXIN",
                "ANTI_ARRHYTHMICS",
                "COPD_ASTHMA_DRUG",
                "STATIN",
                "OTHER_LIPID_DRUG",
                "ANTIPLATELET",
                "ORAL_ANTICOAGULANT",
                "HEPARIN",
                "NSAID",
                "ORAL_CORTICOSTEROID",
                "BISPHOSPHONATE",
                "OPIOID",
                "ANTIDEPRESSANT",
                "ANTIPSYCHOTIC",
                "ANTICONVULSANT",
                "LITHIUM",
                "BENZODIAZAPINE",
                "ANXIOLYTIC_HYPNOTIC",
                "DEMENTIA_DRUG",
                "PARKINSON_DRUG",
                "N_DIAGNOSES",
                "N_DRUG_RX",
                "HOSPITALIZATION",
                "ENDOCRINOLOGIST",
                "INTERNAL",
                "CARDIOLOGIST",
                "ELECTROCARDIOGRAM",
                "GLUCOSE_TESTS",
                "HOSPITALIZATION_30",
                "HOSPITALIZATION_31_180",
                "N_HOSPITALIZATION",
                "N_HOSPITAL_DAYS",
                "N_ED",
                "N_OFFICE",
                "N_ENDOCRINOLOGIST",
                "N_INTERNAL",
                "N_CARDIOLOGIST",
                "N_ELECTROCARDIOGRAM",
                "N_HBA1C",
                "N_GLUCOSE_TEST",
                "N_LIPID_TEST",
                "N_CREATININE_TEST",
                "N_BUN_TEST",
                "N_MICROALBUMINURIA_TEST",
                "DM_DRUG_AGI",
                "DM_DRUG_GLITAZONE",
                "DM_DRUG_GLP1",
                "DM_DRUG_INSULIN",
                "DM_DRUG_MEGLITINIDE",
                "DM_DRUG_METFORMIN",
                "PRAMLINTIDE",
                "GEN1_SU",
                "MONOTHERAPY_INITIATION",
                "METFORMIN_DUAL",
                "AGE",
                "CEREBROVASCULAR_HEM_STROKE",
                "BLADDER_STONE",
                "KIDNEY_STONE",
                "UTI",
                "DIPSTICK_URINALYSIS",
                "NON_DIPSTICK_URINALYSIS",
                "URINE_FUNCTION_TEST",
                "CYTOLOGY",
                "CYSTOSCOPY",
                "FRAILTY_EMPIRICAL",
                "CREATININE_TEST",
                "BUN_TEST",
                "CRI_NO_CKD",
                "CKD_1_2",
                "CKD_3_6",
                "CONCOMITANT_SGLT2I",
                "CONCOMITANT_AGI",
                "CONCOMITANT_GLITAZONE",
                "CONCOMITANT_GLP1",
                "CONCOMITANT_INSULIN",
                "CONCOMITANT_MEGLITINIDE",
                "CONCOMITANT_METFORMIN",
                "FRAILTY_QUALITATIVE_V1",
                "PAST_SGLT2I",
                "PAST_AGI",
                "PAST_GLITAZONE",
                "PAST_GLP1",
                "PAST_INSULIN",
                "PAST_MEGLITINIDE",
                "PAST_METFORMIN",
                "CALENDAR_TIME_DAY",
                "BLADDER_KIDNEY_STONE",
                "PERIPHERAL_GANGRENE_OSTEOMYELITIS",
                "ALCOHOL_DRUG",
                "DM_EYE",
                "COMPOSITE_CVD",
                "COMPOSITE_CARDIAC_PROCEDURE",
                "THYROID",
                "DELIRIUM_PSYCHOSIS",
                "MEGLITINIDE",
                "AGI",
                "GLAUCOMA_CATARACTS",
                "IMAGING",
                "CELLULITIS_ABCESS_TOE",
                "FOOT_ULCER",
                "ENTRESTO",
                "ENDOCRINOLOGIST_30",
                "ENDOCRINOLOGIST_31_180",
                "CARDIOLOGIST_30",
                "CARDIOLOGIST_31_180",
                "INTERNAL_30",
                "INTERNAL_31_180",
                "DIALYSIS",
                "CCI_180",
                "CKD_3_6_DIALYSIS",
                "BASELINE_CV",
                "THIAZIDE",
                "BETA_BLOCKER",
                "CA_CHANNEL_BLOCKER",
                "PROP_WHITE",
                "PROP_BLACK",
                "PROP_AIAN",
                "PROP_ASIAN",
                "PROP_NHP",
                "PROP_OTHER_RACE",
                "PROP_2MORE_RACE",
                "PROP_NONHISP_WHITE",
                "PROP_NONHISP_BLACK",
                "PROP_NONHISP_AIAN",
                "PROP_NONHISP_ASIAN",
                "PROP_HISP",
                "PROP_HISP_WHITE",
                "PROP_HISP_BLACK",
                "PROP_HISP_ASIAN",
                "PROP_ENGLISH_ONLY",
                "PROP_ENGLISH_WELL",
                "PROP_NO_ENGLISH",
                "PROP_PUBLIC_ASST",
                "PROP_VET",
                "PROP_NATIVE",
                "PROP_NEVER_MARRIED",
                "PROP_MARRIED",
                "PROP_DIVORCED",
                "PROP_WIDOWED",
                "PROP_HS",
                "PROP_ASSOC",
                "PROP_BACHELORS",
                "PROP_POVERTY",
                "PROP_DISABILITY",
                "PROP_COGNITIVE",
                "PROP_AMBULATORY",
                "MED_INC",
                "MED_INC_65",
                "GINI")

savor_tableone_all <- CreateTableOne(vars = savor_vars, data = savor, 
                                     factorVars = savor_factor_vars,
                                     strata = "EXPOSURE", test = FALSE, smd = TRUE)
savor_tableone_all_export <- print(savor_tableone_all, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(savor_tableone_all_export, file = "/Volumes/PACS$/Cordes, Jack/savor_tableone_all.csv")

# Table 1 of characteristics for individuals with at least 5 in zip code by exposure group
savor_5 <- savor %>%
  filter(count>=5)

savor_tableone_5_exp <- CreateTableOne(vars = savor_vars, data = savor_5, 
                                       factorVars = savor_factor_vars,
                                       strata = "EXPOSURE", test = FALSE, smd = TRUE)
savor_tableone_5_exp_export <- print(savor_tableone_5_exp, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(savor_tableone_5_exp_export, file = "/Volumes/PACS$/Cordes, Jack/savor_tableone_5_exp.csv")


# Table 1 of characteristics for individuals with at least 5 in zip code vs. fewer
savor$AT_LEAST_5 <- ifelse(savor$count>=5,1,0)

savor_tableone_5_vs_less <- CreateTableOne(vars = c(savor_vars, "EXPOSURE"), data = savor, 
                                           factorVars = c(savor_factor_vars, "EXPOSURE"),
                                           strata = "AT_LEAST_5", test = FALSE, smd = TRUE)
savor_tableone_5_vs_less_export <- print(savor_tableone_5_vs_less, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(savor_tableone_5_vs_less_export, file = "/Volumes/PACS$/Cordes, Jack/savor_tableone_5_vs_less.csv")

# Table 1 of characteristics for all individuals in cohort limited to PS variables
savor_factor_vars_ps <- c("AGE_CAT",
                          "SEX",
                          "RACE",
                          "REGION",
                          "OBESITY",
                          "OVERWEIGHT",
                          "SMOKE",
                          "ALCOHOL",
                          "DRUG",
                          "DM_NEUROPATHY3",
                          "HYPOGLYCEMIA",
                          "HYPERGLYCEMIA",
                          "DM_KETOACIDOSIS",
                          "HYPERTENSION",
                          "HYPERLIPIDEMIA",
                          "UNSTABLE_ANGINA",
                          "STABLE_ANGINA",
                          "CORONARY_ATHEROSCLEROSIS",
                          "TIA",
                          "CHF",
                          "PVD",
                          "CARDIAC_DYSRHYTHMIA",
                          "EDEMA",
                          "COPD",
                          "ASTHMA",
                          "CKD",
                          "LIVER_DISEASE",
                          "DEPRESSION",
                          "ANXIETY",
                          "FRAILTY_EMPIRICAL_V3",
                          "ACE_INHIBITORS",
                          "ARB",
                          "STATIN",
                          "ANTIPLATELET",
                          "NSAID",
                          "ORAL_CORTICOSTEROID",
                          "OPIOID",
                          "BENZODIAZAPINE",
                          "HOSPITALIZATION",
                          "ENDOCRINOLOGIST",
                          "CARDIOLOGIST",
                          "DM_DRUG_AGI",
                          "DM_DRUG_GLITAZONE",
                          "DM_DRUG_GLP1",
                          "DM_DRUG_INSULIN",
                          "DM_DRUG_MEGLITINIDE",
                          "DM_DRUG_METFORMIN",
                          "PRAMLINTIDE",
                          "GEN1_SU",
                          "MONOTHERAPY_INITIATION",
                          "METFORMIN_DUAL",
                          "CEREBROVASCULAR_HEM_STROKE",
                          "URINE_FUNCTION_TEST",
                          "CRI_NO_CKD",
                          "CONCOMITANT_SGLT2I",
                          "CONCOMITANT_AGI",
                          "CONCOMITANT_GLITAZONE",
                          "CONCOMITANT_GLP1",
                          "CONCOMITANT_INSULIN",
                          "CONCOMITANT_MEGLITINIDE",
                          "CONCOMITANT_METFORMIN",
                          "PAST_SGLT2I",
                          "PAST_AGI",
                          "PAST_GLITAZONE",
                          "PAST_GLP1",
                          "PAST_INSULIN",
                          "PAST_MEGLITINIDE",
                          "PAST_METFORMIN",
                          "DM_EYE",
                          "COMPOSITE_CVD",
                          "COMPOSITE_CARDIAC_PROCEDURE",
                          "THYROID",
                          "DIALYSIS",
                          "BETA_BLOCKER",
                          "CA_CHANNEL_BLOCKER")

savor_vars_ps <- c("AGE_CAT",
                   "SEX",
                   "RACE",
                   "REGION",
                   "OBESITY",
                   "OVERWEIGHT",
                   "SMOKE",
                   "ALCOHOL",
                   "DRUG",
                   "DM_NEUROPATHY3",
                   "HYPOGLYCEMIA",
                   "HYPERGLYCEMIA",
                   "DM_KETOACIDOSIS",
                   "HYPERTENSION",
                   "HYPERLIPIDEMIA",
                   "UNSTABLE_ANGINA",
                   "STABLE_ANGINA",
                   "CORONARY_ATHEROSCLEROSIS",
                   "TIA",
                   "CHF",
                   "PVD",
                   "CARDIAC_DYSRHYTHMIA",
                   "EDEMA",
                   "COPD",
                   "ASTHMA",
                   "CKD",
                   "LIVER_DISEASE",
                   "DEPRESSION",
                   "ANXIETY",
                   "FRAILTY_EMPIRICAL_V3",
                   "N_ANTIDIABETICS",
                   "ACE_INHIBITORS",
                   "ARB",
                   "STATIN",
                   "ANTIPLATELET",
                   "NSAID",
                   "ORAL_CORTICOSTEROID",
                   "OPIOID",
                   "BENZODIAZAPINE",
                   "N_DIAGNOSES",
                   "N_DRUG_RX",
                   "HOSPITALIZATION",
                   "ENDOCRINOLOGIST",
                   "CARDIOLOGIST",
                   "N_OFFICE",
                   "N_HBA1C",
                   "N_GLUCOSE_TEST",
                   "DM_DRUG_AGI",
                   "DM_DRUG_GLITAZONE",
                   "DM_DRUG_GLP1",
                   "DM_DRUG_INSULIN",
                   "DM_DRUG_MEGLITINIDE",
                   "DM_DRUG_METFORMIN",
                   "PRAMLINTIDE",
                   "GEN1_SU",
                   "MONOTHERAPY_INITIATION",
                   "METFORMIN_DUAL",
                   "CEREBROVASCULAR_HEM_STROKE",
                   "URINE_FUNCTION_TEST",
                   "CRI_NO_CKD",
                   "CONCOMITANT_SGLT2I",
                   "CONCOMITANT_AGI",
                   "CONCOMITANT_GLITAZONE",
                   "CONCOMITANT_GLP1",
                   "CONCOMITANT_INSULIN",
                   "CONCOMITANT_MEGLITINIDE",
                   "CONCOMITANT_METFORMIN",
                   "PAST_SGLT2I",
                   "PAST_AGI",
                   "PAST_GLITAZONE",
                   "PAST_GLP1",
                   "PAST_INSULIN",
                   "PAST_MEGLITINIDE",
                   "PAST_METFORMIN",
                   "CALENDAR_TIME_DAY",
                   "DM_EYE",
                   "COMPOSITE_CVD",
                   "COMPOSITE_CARDIAC_PROCEDURE",
                   "THYROID",
                   "DIALYSIS",
                   "CCI_180",
                   "BETA_BLOCKER",
                   "CA_CHANNEL_BLOCKER",
                   "PROP_WHITE",
                   "PROP_BLACK",
                   "PROP_AIAN",
                   "PROP_ASIAN",
                   "PROP_OTHER_RACE",
                   "PROP_HISP",
                   "PROP_NO_ENGLISH",
                   "PROP_PUBLIC_ASST",
                   "PROP_NATIVE",
                   "PROP_MARRIED",
                   "PROP_BACHELORS",
                   "PROP_POVERTY",
                   "PROP_DISABILITY",
                   "PROP_COGNITIVE",
                   "PROP_AMBULATORY",
                   "MED_INC",
                   "MED_INC_65",
                   "GINI")

savor_vars_zip <- c("PROP_WHITE",
                    "PROP_BLACK",
                    "PROP_AIAN",
                    "PROP_ASIAN",
                    "PROP_OTHER_RACE",
                    "PROP_HISP",
                    "PROP_NO_ENGLISH",
                    "PROP_PUBLIC_ASST",
                    "PROP_NATIVE",
                    "PROP_MARRIED",
                    "PROP_BACHELORS",
                    "PROP_POVERTY",
                    "PROP_DISABILITY",
                    "PROP_COGNITIVE",
                    "PROP_AMBULATORY",
                    "MED_INC",
                    "MED_INC_65",
                    "GINI")

savor_tableone_all_ps <- CreateTableOne(vars = savor_vars_ps, data = savor, 
                                        factorVars = savor_factor_vars_ps,
                                        strata = "EXPOSURE", test = FALSE, smd = TRUE)
savor_tableone_all_ps_export <- print(savor_tableone_all_ps, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(savor_tableone_all_ps_export, file = "/Volumes/PACS$/Cordes, Jack/savor_tableone_all_ps_14_17.csv")

# Table 1 of characteristics for individuals with at least 5 in zip code by exposure group
savor_5 <- savor %>%
  filter(count>=5)

savor_tableone_5_exp_ps <- CreateTableOne(vars = savor_vars_ps, data = savor_5, 
                                          factorVars = savor_factor_vars_ps,
                                          strata = "EXPOSURE", test = FALSE, smd = TRUE)
savor_tableone_5_exp_ps_export <- print(savor_tableone_5_exp_ps, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(savor_tableone_5_exp_ps_export, file = "/Volumes/PACS$/Cordes, Jack/savor_tableone_5_exp_ps.csv")


# Table 1 of characteristics for individuals with at least 5 in zip code vs. fewer, ps variables only
savor$AT_LEAST_5 <- ifelse(savor$count > 4,1,0)

savor_tableone_5_vs_less_ps <- CreateTableOne(vars = c(savor_vars_ps, "EXPOSURE"), data = savor, 
                                              factorVars = c(savor_factor_vars_ps, "EXPOSURE"),
                                              strata = "AT_LEAST_5", test = FALSE, smd = TRUE)
savor_tableone_5_vs_less_ps_export <- print(savor_tableone_5_vs_less_ps, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(savor_tableone_5_vs_less_ps_export, file = "/Volumes/PACS$/Cordes, Jack/savor_tableone_5_vs_less_ps.csv")

# Mahalanobis distance of PS variables by treatment group for all individuals
savor_mah_dist_total <- savor[,c(1:13,19:21,23,34:35,38,40:41,48,52:53,55,58:60,66,70,80:81,87,89,91:92,100,102,105:106,108,113,117:120,122,130,135:136,141:142,144:150,152,158,164,167:169,171:173,175:177,179:182,187:190,206,210:211,220,226:229,231,237,244:245,247,249,254:258,262:264)]

savor_mah_dist_total_treatment <- mahalanobis(colMeans(savor_mah_dist_total %>% filter(EXPOSURE == 1) %>% select(c(5:83,85:102))),
                                              colMeans(savor_mah_dist_total %>% filter(EXPOSURE == 0) %>% select(c(5:83,85:102))),
                                              cov(savor_mah_dist_total %>% select(c(5:83,85:102))),
                                              tol = 1e-30)
print(savor_mah_dist_total_treatment)

# Mahalanobis distance of PS variables by treatment group at least 5 in ZCTA
savor_mah_dist_total_5 <- filter(savor[,c(1:13,19:21,23,34:35,38,40:41,48,52:53,55,58:60,66,70,80:81,87,89,91:92,100,102,105:106,108,113,117:120,122,130,135:136,141:142,144:150,152,158,164,167:169,171:173,175:177,179:182,187:190,206,210:211,220,226:229,231,237,244:245,247,249,254:258,262:264,268)], count > 4)

savor_mah_dist_total_treatment_5 <- mahalanobis(colMeans(savor_mah_dist_total_5 %>% filter(EXPOSURE == 1) %>% select(c(5:83,85:102))),
                                                colMeans(savor_mah_dist_total_5 %>% filter(EXPOSURE == 0) %>% select(c(5:83,85:102))),
                                                cov(savor_mah_dist_total_5 %>% select(c(5:83,85:102))),
                                                tol = 1e-30)
print(savor_mah_dist_total_treatment_5)

# Mahalanobis distance of ZCTA variables by treatment group at least 5 in ZCTA
savor_mah_dist_total_5_zcta <- filter(savor[,c(1:4,220,226:229,231,237,244:245,247,249,254:258,262:264,268)], count > 4)

savor_mah_dist_total_treatment_5_zcta <- mahalanobis(colMeans(savor_mah_dist_total_5_zcta %>% filter(EXPOSURE == 1) %>% select(c(6:23))),
                                                     colMeans(savor_mah_dist_total_5_zcta %>% filter(EXPOSURE == 0) %>% select(c(6:23))),
                                                     cov(savor_mah_dist_total_5_zcta %>% select(c(6:23))),
                                                     tol = 1e-30)
print(savor_mah_dist_total_treatment_5_zcta)

# Add LISA quadrant designations to allow limiting analyses to non-clustered zipcodes
savor <- left_join(savor, st_drop_geometry(zcta_shp_savor[,c(1,37)]), by = "ZCTA_USE")

# Mahalanobis distance of PS variables only by treatment group with at least 5 in ZCTA limited to nonclustered
savor_mah_dist_total_5_noclusters <- filter(savor[,c(1:13,19:21,23,34:35,38,40:41,48,52:53,55,58:60,66,70,80:81,87,89,91:92,100,102,105:106,108,113,117:120,122,130,135:136,141:142,144:150,152,158,164,167:169,171:173,175:177,179:182,187:190,206,210:211,220,226:229,231,237,244:245,247,249,254:258,262:264,268,271)], count > 4 & quadrant == 0)
savor_mah_dist_total_treatment_5_noclusters <- mahalanobis(colMeans(savor_mah_dist_total_5_noclusters %>% filter(EXPOSURE == 1) %>% select(c(6:23))),
                                                           colMeans(savor_mah_dist_total_5_noclusters %>% filter(EXPOSURE == 0) %>% select(c(6:23))),
                                                           cov(savor_mah_dist_total_5_noclusters %>% select(c(6:23))),
                                                           tol = 1e-30)
print(savor_mah_dist_total_treatment_5_noclusters)

# Mahalanobis distance of PS variables by treatment group for 2014-2017
savor_mah_dist_total_14_17 <- filter(savor[,c(1:13,19:21,23,34:35,38,40:41,48,52:53,55,58:60,66,70,80:81,87,89,91:92,100,102,105:106,108,113,117:120,122,130,135:136,141:142,144:150,152,158,164,167:169,171:173,175:177,179:182,187:190,206,210:211,220,226:229,231,237,244:245,247,249,254:258,262:264)], YEAR >= 2014)

savor_mah_dist_total_14_17_treatment <- mahalanobis(colMeans(savor_mah_dist_total_14_17 %>% filter(EXPOSURE == 1) %>% select(c(5:83,85:102))),
                                                    colMeans(savor_mah_dist_total_14_17 %>% filter(EXPOSURE == 0) %>% select(c(5:83,85:102))),
                                                    cov(savor_mah_dist_total_14_17 %>% select(c(5:83,85:102))),
                                                    tol = 1e-30)
print(savor_mah_dist_total_14_17_treatment)

# Risk factor analysis
# Among Saxagliptin
summary(
  glm(data = filter(savor, EXPOSURE == 1), formula = OUTCOME ~ PROP_WHITE,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 1), formula = OUTCOME ~ PROP_BLACK,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 1), formula = OUTCOME ~ PROP_AIAN,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 1), formula = OUTCOME ~ PROP_ASIAN,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 1), formula = OUTCOME ~ PROP_OTHER_RACE,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 1), formula = OUTCOME ~ PROP_HISP,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 1), formula = OUTCOME ~ PROP_NO_ENGLISH,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 1), formula = OUTCOME ~ PROP_PUBLIC_ASST,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 1), formula = OUTCOME ~ PROP_NATIVE,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 1), formula = OUTCOME ~ PROP_MARRIED,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 1), formula = OUTCOME ~ PROP_BACHELORS,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 1), formula = OUTCOME ~ PROP_POVERTY,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 1), formula = OUTCOME ~ PROP_DISABILITY,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 1), formula = OUTCOME ~ PROP_COGNITIVE,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 1), formula = OUTCOME ~ PROP_AMBULATORY,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 1), formula = OUTCOME ~ MED_INC_STD,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 1), formula = OUTCOME ~ MED_INC_65_STD,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 1), formula = OUTCOME ~ GINI,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)

# Among SU
summary(
  glm(data = filter(savor, EXPOSURE == 0), formula = OUTCOME ~ PROP_WHITE,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 0), formula = OUTCOME ~ PROP_BLACK,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 0), formula = OUTCOME ~ PROP_AIAN,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 0), formula = OUTCOME ~ PROP_ASIAN,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 0), formula = OUTCOME ~ PROP_OTHER_RACE,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 0), formula = OUTCOME ~ PROP_HISP,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 0), formula = OUTCOME ~ PROP_NO_ENGLISH,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 0), formula = OUTCOME ~ PROP_PUBLIC_ASST,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 0), formula = OUTCOME ~ PROP_NATIVE,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 0), formula = OUTCOME ~ PROP_MARRIED,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 0), formula = OUTCOME ~ PROP_BACHELORS,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 0), formula = OUTCOME ~ PROP_POVERTY,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 0), formula = OUTCOME ~ PROP_DISABILITY,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 0), formula = OUTCOME ~ PROP_COGNITIVE,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 0), formula = OUTCOME ~ PROP_AMBULATORY,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 0), formula = OUTCOME ~ MED_INC_STD,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 0), formula = OUTCOME ~ MED_INC_65_STD,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = filter(savor, EXPOSURE == 0), formula = OUTCOME ~ GINI,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)

# Full cohort
summary(
  glm(data = savor, formula = OUTCOME ~ PROP_WHITE,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = savor, formula = OUTCOME ~ PROP_BLACK,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = savor, formula = OUTCOME ~ PROP_AIAN,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = savor, formula = OUTCOME ~ PROP_ASIAN,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = savor, formula = OUTCOME ~ PROP_OTHER_RACE,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = savor, formula = OUTCOME ~ PROP_HISP,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = savor, formula = OUTCOME ~ PROP_NO_ENGLISH,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = savor, formula = OUTCOME ~ PROP_PUBLIC_ASST,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = savor, formula = OUTCOME ~ PROP_NATIVE,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = savor, formula = OUTCOME ~ PROP_MARRIED,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = savor, formula = OUTCOME ~ PROP_BACHELORS,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = savor, formula = OUTCOME ~ PROP_POVERTY,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = savor, formula = OUTCOME ~ PROP_DISABILITY,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = savor, formula = OUTCOME ~ PROP_COGNITIVE,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = savor, formula = OUTCOME ~ PROP_AMBULATORY,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = savor, formula = OUTCOME ~ MED_INC_STD,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = savor, formula = OUTCOME ~ MED_INC_65_STD,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)
summary(
  glm(data = savor, formula = OUTCOME ~ GINI,
      family = binomial(link = "logit"),
      na.action = na.exclude)
)

# Unadjusted association
# Distribution of follow-up times
summary(savor$FUP_TIME)

savor %>%
  ggplot() + 
  geom_histogram(aes(FUP_TIME), binwidth = 30, color = "black", fill = "white") + 
  geom_vline(xintercept = 299.5, color = "red") +
  scale_y_continuous(labels = label_comma()) +
  coord_cartesian(xlim = c(0, 2200)) +
  labs(title = "Follow-up Time Distribution\nSaxagliptin", x = "Days", y = "Count") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18))
ggsave(filename = "savor_followup_dist.pdf")

# Raw outcome counts by exposure status
table(savor$EXPOSURE, savor$OUTCOME)
savor %>% group_by(EXPOSURE) %>% summarise(OUTCOME_SUM = sum(OUTCOME))

# Person-years by exposure status
savor %>% group_by(EXPOSURE) %>% summarise(FUP_TIME_SUM = sum(FUP_TIME)/365.25)

# Risk in the exposed and unexposed
savor_exp_risk <- 996/(996 + 37715)
savor_unexp_risk <- 33119/(33119 + 744827)

# Unadjusted risk ratio and odds ratio
savor_rr <- savor_exp_risk/savor_unexp_risk

savor_or <- (284588*5785)/(21714*97892)

# Unadjusted ITT cumulative incidence analysis using logistic regression
savor_model_crude_itt <- glm(data = savor, formula = OUTCOME ~ EXPOSURE,
                             family = binomial(link = "logit"),
                             na.action = na.exclude)
summary(savor_model_crude_itt)
exp(savor_model_crude_itt$coefficients)

# Unadjusted survival analysis
savor_model_crude_survival <- coxph(data = savor, formula = Surv(FUP_TIME, OUTCOME) ~ EXPOSURE,
                                    na.action = na.exclude)
summary(savor_model_crude_survival)
cox.zph(savor_model_crude_survival)
ggcoxzph(cox.zph(savor_model_crude_survival))
ggsave(filename = "savor_model_crude_survival.pdf")

# Kaplan-Meier curves
savor_km <- survfit(data = savor, formula = Surv(FUP_TIME, OUTCOME) ~ EXPOSURE)

center_title <- function() {
  theme_survminer() %+replace%
    theme(plot.title = element_text(size = 18, hjust = 0),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18))
}

savor_km_plot <- ggsurvplot(fit = savor_km,
                            palette = c("Red", "Blue"),
                            censor = FALSE,
                            title = "B                                      Survival Curves\n                                      Saxagliptin Cohort",
                            xlab = "Days",
                            ylab = "Survival Probability",
                            legend.title = "",
                            legend.labs = c("SU", "Saxagliptin"),
                            legend = "right",
                            conf.int = T,
                            ggtheme = center_title())
ggsave("savor_km.pdf", savor_km_plot$plot)

# Generate propensity scores based on 2014-2017 cohort
# Subset data to 2014-2017
savor_phys <- filter(savor, YEAR >= 2014)

# Raw outcome counts by exposure status
table(savor_phys$EXPOSURE, savor_phys$OUTCOME)
savor_phys %>% group_by(EXPOSURE) %>% summarise(OUTCOME_SUM = sum(OUTCOME))

# Unadjusted risk ratio analysis linear regression all follow-up
savor_model_crude_linear <- glm(data = savor_phys, formula = OUTCOME ~ EXPOSURE,
                                na.action = na.exclude)
summary(savor_model_crude_linear)

# Analysis adjusted for zipcode variables
savor_ps_model_zip_14_17 <- glm(data = savor_phys, formula = EXPOSURE ~ 
                                  as.factor(AGE_CAT) +
                                  as.factor(SEX) +
                                  as.factor(RACE) +
                                  as.factor(REGION) +
                                  as.factor(OBESITY) +
                                  as.factor(OVERWEIGHT) +
                                  as.factor(SMOKE) +
                                  as.factor(ALCOHOL) +
                                  as.factor(DRUG) +
                                  as.factor(DM_NEUROPATHY3) +
                                  as.factor(HYPOGLYCEMIA) +
                                  as.factor(HYPERGLYCEMIA) +
                                  as.factor(DM_KETOACIDOSIS) +
                                  as.factor(HYPERTENSION) +
                                  as.factor(HYPERLIPIDEMIA) +
                                  as.factor(UNSTABLE_ANGINA) +
                                  as.factor(STABLE_ANGINA) +
                                  as.factor(CORONARY_ATHEROSCLEROSIS) +
                                  as.factor(TIA) +
                                  as.factor(CHF) +
                                  as.factor(PVD) +
                                  as.factor(CARDIAC_DYSRHYTHMIA) +
                                  as.factor(EDEMA) +
                                  as.factor(COPD) +
                                  as.factor(ASTHMA) +
                                  as.factor(CKD) +
                                  as.factor(LIVER_DISEASE) +
                                  as.factor(DEPRESSION) +
                                  as.factor(ANXIETY) +
                                  as.factor(FRAILTY_EMPIRICAL_V3) +
                                  N_ANTIDIABETICS +
                                  as.factor(ACE_INHIBITORS) +
                                  as.factor(ARB) +
                                  as.factor(STATIN) +
                                  as.factor(ANTIPLATELET) +
                                  as.factor(NSAID) +
                                  as.factor(ORAL_CORTICOSTEROID) +
                                  as.factor(OPIOID) +
                                  as.factor(BENZODIAZAPINE) +
                                  N_DIAGNOSES +
                                  N_DRUG_RX +
                                  as.factor(HOSPITALIZATION) +
                                  as.factor(ENDOCRINOLOGIST) +
                                  as.factor(CARDIOLOGIST) +
                                  N_OFFICE +
                                  N_HBA1C +
                                  N_GLUCOSE_TEST +
                                  as.factor(DM_DRUG_AGI) +
                                  as.factor(DM_DRUG_GLITAZONE) +
                                  #                                  as.factor(DM_DRUG_GLP1) +
                                  as.factor(DM_DRUG_INSULIN) +
                                  as.factor(DM_DRUG_MEGLITINIDE) +
                                  as.factor(DM_DRUG_METFORMIN) +
                                  as.factor(PRAMLINTIDE) +
                                  as.factor(GEN1_SU) +
                                  as.factor(MONOTHERAPY_INITIATION) +
                                  as.factor(METFORMIN_DUAL) +
                                  as.factor(CEREBROVASCULAR_HEM_STROKE) +
                                  as.factor(URINE_FUNCTION_TEST) +
                                  as.factor(CRI_NO_CKD) +
                                  as.factor(CONCOMITANT_SGLT2I) +
                                  as.factor(CONCOMITANT_AGI) +
                                  as.factor(CONCOMITANT_GLITAZONE) +
                                  #                                  as.factor(CONCOMITANT_GLP1) +
                                  as.factor(CONCOMITANT_INSULIN) +
                                  as.factor(CONCOMITANT_MEGLITINIDE) +
                                  as.factor(CONCOMITANT_METFORMIN) +
                                  as.factor(PAST_SGLT2I) +
                                  as.factor(PAST_AGI) +
                                  as.factor(PAST_GLITAZONE) +
                                  #                                  as.factor(PAST_GLP1) +
                                  as.factor(PAST_INSULIN) +
                                  as.factor(PAST_MEGLITINIDE) +
                                  as.factor(PAST_METFORMIN) +
                                  CALENDAR_TIME_DAY +
                                  as.factor(DM_EYE) +
                                  as.factor(COMPOSITE_CVD) +
                                  as.factor(COMPOSITE_CARDIAC_PROCEDURE) +
                                  as.factor(THYROID) +
                                  #                                  as.factor(DIALYSIS) +
                                  CCI_180 +
                                  as.factor(BETA_BLOCKER) +
                                  as.factor(CA_CHANNEL_BLOCKER) +
                                  PROP_WHITE +
                                  PROP_BLACK +
                                  PROP_AIAN +
                                  PROP_ASIAN +
                                  PROP_OTHER_RACE +
                                  PROP_HISP +
                                  PROP_NO_ENGLISH +
                                  PROP_PUBLIC_ASST +
                                  PROP_NATIVE +
                                  PROP_MARRIED +
                                  PROP_BACHELORS +
                                  PROP_POVERTY +
                                  PROP_DISABILITY +
                                  PROP_COGNITIVE +
                                  PROP_AMBULATORY +
                                  MED_INC_STD +
                                  MED_INC_65_STD +
                                  GINI, family = binomial(link = "logit"),
                                na.action = na.exclude)

savor_ps_model_zip_14_17_pred <- predict(savor_ps_model_zip_14_17, type = "response")
roc(savor_phys$EXPOSURE, savor_ps_model_zip_14_17_pred)

# Calculate PS
savor_phys$PS_ZIP_14_17 <- fitted(savor_ps_model_zip_14_17)
summary(savor_phys$PS_ZIP_14_17)

# PS quintiles
savor_phys %>%
  dplyr::summarize(min = quantile(PS_ZIP_14_17, probs = 0),
                   q1 = quantile(PS_ZIP_14_17, probs = 0.2),
                   q2 = quantile(PS_ZIP_14_17, probs = 0.4),
                   q3 = quantile(PS_ZIP_14_17, probs = 0.6),
                   q4 = quantile(PS_ZIP_14_17, probs = 0.8),
                   max = quantile(PS_ZIP_14_17, probs = 1))

# PS distribution by exposure
savor_ps_dist_zip_14_17_plot <- ggplot(savor_phys) +
  geom_density(aes(PS_ZIP_14_17, group = factor(EXPOSURE), color = factor(EXPOSURE))) + 
  labs(title = "Saxagliptin Cohort Propensity Scores\nwith ZCTA Variables 2014-2017", x = "Propensity Score", y = "Density") +
  scale_color_manual(name = "", values = c("red", "blue"), labels = c("SU", "Saxagliptin")) +
  scale_x_continuous(limits = c(0, 1)) +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(hjust = 0.5),
        legend.text = element_text(size = 18))
ggsave("savor_ps_dist_zip_14_17.pdf", savor_ps_dist_zip_14_17_plot)

# Create PS quintiles variable
savor_phys <- savor_phys %>% mutate(
  PS_ZIP_QUINTILE_14_17 = case_when(
    PS_ZIP_14_17 < 0.01856576 ~ "1",
    PS_ZIP_14_17 >= 0.01856576 & PS_ZIP_14_17 < 0.02601198 ~ "2",
    PS_ZIP_14_17 >= 0.02601198 & PS_ZIP_14_17 < 0.03494436 ~ "3",
    PS_ZIP_14_17 >= 0.03494436 & PS_ZIP_14_17 < 0.04919348 ~ "4",
    TRUE ~ "5"
  )
)

# Global historical zcta preference IV analysis
# Sort on zcta and CED
savor_phys <- savor_phys[order(savor_phys$ZCTA_USE, savor_phys$ENTRYDATE),]

# Create empty global historical zcta preference column
savor_phys <- savor_phys %>% mutate(prop_exp_hist = NA)

# Create empty counter column
savor_phys <- savor_phys %>% mutate(counter = NA)

# Create empty cumulative sum column
savor_phys <- savor_phys %>% mutate(cum_sum = NA)

# Create global historical zcta preference variable based on average of all previous prescriptions
for (i in 2:nrow(savor_phys)){
  if ((savor_phys$ZCTA_USE[i] == savor_phys$ZCTA_USE[i-1]) & (is.na(savor_phys$prop_exp_hist[i-1]) == TRUE)){
    savor_phys$prop_exp_hist[i] <- savor_phys$EXPOSURE[i-1]
  }
  if ((savor_phys$ZCTA_USE[i] == savor_phys$ZCTA_USE[i-1]) & (is.na(savor_phys$prop_exp_hist[i-1]) == TRUE)){
    savor_phys$counter[i] <- 1
  }
  if ((savor_phys$ZCTA_USE[i] == savor_phys$ZCTA_USE[i-1]) & (is.na(savor_phys$prop_exp_hist[i-1]) == TRUE)){
    savor_phys$cum_sum[i] <- savor_phys$EXPOSURE[i-1]
  }
  if ((savor_phys$ZCTA_USE[i] == savor_phys$ZCTA_USE[i-1]) & (is.na(savor_phys$prop_exp_hist[i-1]) == FALSE)){
    savor_phys$counter[i] <- savor_phys$counter[i-1]+1
  }
  if ((savor_phys$ZCTA_USE[i] == savor_phys$ZCTA_USE[i-1]) & (is.na(savor_phys$prop_exp_hist[i-1]) == FALSE)){
    savor_phys$cum_sum[i] <- savor_phys$cum_sum[i-1] + savor_phys$EXPOSURE[i-1]
  }
  if ((savor_phys$ZCTA_USE[i] == savor_phys$ZCTA_USE[i-1]) & (is.na(savor_phys$prop_exp_hist[i-1]) == FALSE)){
    savor_phys$prop_exp_hist[i] <- savor_phys$cum_sum[i]/savor_phys$counter[i]
  }
}

# Distribution of global historical zcta proportions
summary(savor_phys$prop_exp_hist)

# Create dataset limited to non-missing zcta historical IV
savor_zip_hist <- filter(savor_phys, is.na(prop_exp_hist) == FALSE)

savor_zip_hist_dist_plot <- savor_zip_hist %>%
  ggplot() + 
  geom_histogram(aes(prop_exp_hist), binwidth = 0.05, color = "black", fill = "white") + 
  geom_vline(xintercept = 0.046, color = "red") +
  scale_y_continuous(labels = label_comma()) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(title = "Saxagliptin Cohort\nCumulative ZCTA Proportion for Each Individual", x = "Proportion Saxagliptin", y = "Count of Individuals") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(hjust = 0.5))
savor_zip_hist_dist_grob <- gtable_add_grob(ggplot_gtable(ggplot_build(savor_zip_hist_dist_plot)), b_label, t = 1, l = 4, b = 6)
ggsave("savor_zip_hist_dist.pdf", savor_zip_hist_dist_grob)

# For IV with proportion saxagliptin cutoffs 100% vs. 0%, instrument is collapsed with exposure
# because no variation in individual-level exposure
# So becomes simple regression comparing individuals with 100% to 0% prescribers

# ZCTA proportion saxagliptin 100% vs. 0%
savor_phys <- savor_phys %>%
  mutate(ZIP_HIST_0_100 = case_when(
    prop_exp_hist == 1 ~ "1",
    prop_exp_hist == 0 ~ "0",
    TRUE ~ NA_character_
  )
  )

# Ordinary least squares - no adjustment for confounders
savor_model_zip_hist_ols_0_100 <- glm(data = filter(savor_phys, (prop_exp_hist == 0 | prop_exp_hist == 1)),
                                      formula = OUTCOME ~ EXPOSURE,
                                      na.action = na.exclude)
summary(savor_model_zip_hist_ols_0_100)

# Ordinary least squares - adjusted for PS
savor_model_zip_hist_ols_0_100_adj <- glm(data = filter(savor_phys, (prop_exp_hist == 0 | prop_exp_hist == 1)),
                                          formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                          na.action = na.exclude)
summary(savor_model_zip_hist_ols_0_100_adj)

# Logistic regression - no adjustment for confounders
savor_model_zip_hist_logistic_0_100 <- glm(data = filter(savor_phys, (prop_exp_hist == 0 | prop_exp_hist == 1)),
                                           formula = OUTCOME ~ EXPOSURE,
                                           family = binomial(link = "logit"),
                                           na.action = na.exclude)
summary(savor_model_zip_hist_logistic_0_100)
exp(savor_model_zip_hist_logistic_0_100$coefficients)

# Logistic regression - adjusted for PS
savor_model_zip_hist_logistic_0_100_adj <- glm(data = filter(savor_phys, (prop_exp_hist == 0 | prop_exp_hist == 1)),
                                               formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                               family = binomial(link = "logit"),
                                               na.action = na.exclude)
summary(savor_model_zip_hist_logistic_0_100_adj)
exp(savor_model_zip_hist_logistic_0_100_adj$coefficients)

# Table 1 at the individual level
savor_tableone_all_iv_zip_hist_0_100 <- CreateTableOne(vars = savor_vars_ps, data = filter(savor_phys, is.na(ZIP_HIST_0_100) == FALSE),
                                                       factorVars = savor_factor_vars_ps,
                                                       strata = "ZIP_HIST_0_100", test = FALSE, smd = TRUE)
savor_tableone_all_iv_zip_hist_0_100_export <- print(savor_tableone_all_iv_zip_hist_0_100, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(savor_tableone_all_iv_zip_hist_0_100_export, file = "/Volumes/PACS$/Cordes, Jack/savor_tableone_all_iv_zip_hist_0_100.csv")

# Check balance improvement by percent change in Mahalanobis distance using treatment
savor_zip_hist_iv_0_100_mah_dist <- savor_phys[,c(1:13,19:21,23,34:35,38,40:41,48,52:53,55,58:60,66,70,80:81,87,89,91:92,100,102,105:106,108,113,117:120,122,130,135:136,141:142,144:150,152,158,164,167:169,171:173,175:177,179:182,187:190,206,210:211,220,226:229,231,237,244:245,247,249,254:258,262:264,271,291)]

savor_mah_dist_total_zip_hist_iv_0_100 <- mahalanobis(colMeans(savor_zip_hist_iv_0_100_mah_dist %>% filter(ZIP_HIST_0_100 == 1) %>% select(c(5:83,85:102))),
                                                      colMeans(savor_zip_hist_iv_0_100_mah_dist %>% filter(ZIP_HIST_0_100 == 0) %>% select(c(5:83,85:102))),
                                                      cov(savor_zip_hist_iv_0_100_mah_dist %>% select(c(5:83,85:102))),
                                                      tol = 1e-30)

((savor_mah_dist_total_zip_hist_iv_0_100-savor_mah_dist_total_14_17_treatment)/(savor_mah_dist_total_14_17_treatment))*100

# Distribution of outcomes by zipcode exposure and by actual exposure
savor_iv_zip_hist_0_100_table <- filter(savor_phys, is.na(ZIP_HIST_0_100) == F)
table(savor_iv_zip_hist_0_100_table$OUTCOME,
      savor_iv_zip_hist_0_100_table$EXPOSURE,
      savor_iv_zip_hist_0_100_table$ZIP_HIST_0_100)

# ZCTA proportion saxagliptin ≥90% vs. ≤10%
savor_phys <- savor_phys %>%
  mutate(ZIP_HIST_10_90 = case_when(
    prop_exp_hist >= 0.9 ~ "1",
    prop_exp_hist <= 0.1 ~ "0",
    TRUE ~ NA_character_
  )
  )

# Check strength of instrument - linear regression
savor_model_iv_zip_hist_strength_10_90 <- glm(data = savor_phys,
                                              formula = EXPOSURE ~ ZIP_HIST_10_90 + as.factor(PS_ZIP_QUINTILE_14_17),
                                              na.action = na.exclude)
summary(savor_model_iv_zip_hist_strength_10_90)
linearHypothesis(savor_model_iv_zip_hist_strength_10_90, "ZIP_HIST_10_901 = 0", test = "F")

# Ordinary least squares - no adjustment for confounders
savor_model_zip_hist_ols_10_90 <- glm(data = filter(savor_phys, (prop_exp_hist <= 0.1 | prop_exp_hist >= 0.9)),
                                      formula = OUTCOME ~ EXPOSURE,
                                      na.action = na.exclude)
summary(savor_model_zip_hist_ols_10_90)

# 2 stage least squares - no adjustment for confounders
savor_model_zip_hist_tsls_10_90 <- ivreg(OUTCOME ~ EXPOSURE | ZIP_HIST_10_90,
                                         data = savor_phys)
summary(savor_model_zip_hist_tsls_10_90, vcov = sandwich, diagnostics = TRUE)

# Ordinary least squares - adjusted for PS
savor_model_zip_hist_ols_10_90_adj <- glm(data = filter(savor_phys, (prop_exp_hist <= 0.1 | prop_exp_hist >= 0.9)),
                                          formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                          na.action = na.exclude)
summary(savor_model_zip_hist_ols_10_90_adj)

# 2 stage least squares - adjusted for PS
savor_model_zip_hist_tsls_10_90_adj <- ivreg(OUTCOME ~ as.factor(PS_ZIP_QUINTILE_14_17) | EXPOSURE | ZIP_HIST_10_90,
                                             data = savor_phys)
summary(savor_model_zip_hist_tsls_10_90_adj, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - no adjustment for confounders
savor_model_zip_hist_logistic_10_90 <- glm(data = filter(savor_phys, (prop_exp_hist <= 0.1 | prop_exp_hist >= 0.9)),
                                           formula = OUTCOME ~ EXPOSURE,
                                           family = binomial(link = "logit"),
                                           na.action = na.exclude)
summary(savor_model_zip_hist_logistic_10_90)
exp(savor_model_zip_hist_logistic_10_90$coefficients)

# 2 stage residual inclusion - no adjustment for confounders
savor_model_zip_hist_tsri_10_90 <- tsri(OUTCOME ~ EXPOSURE | ZIP_HIST_10_90,
                                        data = savor_phys,
                                        link = "logit")
summary(savor_model_zip_hist_tsri_10_90, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - adjusted for PS
savor_model_zip_hist_logistic_10_90_adj <- glm(data = filter(savor_phys, (prop_exp_hist <= 0.1 | prop_exp_hist >= 0.9)),
                                               formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                               family = binomial(link = "logit"),
                                               na.action = na.exclude)
summary(savor_model_zip_hist_logistic_10_90_adj)
exp(savor_model_zip_hist_logistic_10_90_adj$coefficients)

# 2 stage residual inclusion - adjusted for PS
savor_model_zip_hist_tsri_10_90_adj <- tsri(OUTCOME ~ PS_ZIP_QUINTILE_14_17 | EXPOSURE | ZIP_HIST_10_90,
                                            data = savor_phys,
                                            link = "logit")
summary(savor_model_zip_hist_tsri_10_90_adj, vcov = sandwich, diagnostics = TRUE)

# Table 1 at the individual level
savor_tableone_all_iv_zip_hist_10_90 <- CreateTableOne(vars = savor_vars_ps, data = filter(savor_phys, is.na(ZIP_HIST_10_90) == FALSE),
                                                       factorVars = savor_factor_vars_ps,
                                                       strata = "ZIP_HIST_10_90", test = FALSE, smd = TRUE)
savor_tableone_all_iv_zip_hist_10_90_export <- print(savor_tableone_all_iv_zip_hist_10_90, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(savor_tableone_all_iv_zip_hist_10_90_export, file = "/Volumes/PACS$/Cordes, Jack/savor_tableone_all_iv_zip_hist_10_90.csv")

# Check balance improvement by percent change in Mahalanobis distance using treatment
savor_zip_hist_iv_10_90_mah_dist <- savor_phys[,c(1:13,19:21,23,34:35,38,40:41,48,52:53,55,58:60,66,70,80:81,87,89,91:92,100,102,105:106,108,113,117:120,122,130,135:136,141:142,144:150,152,158,164,167:169,171:173,175:177,179:182,187:190,206,210:211,220,226:229,231,237,244:245,247,249,254:258,262:264,271,292)]

savor_mah_dist_total_zip_hist_iv_10_90 <- mahalanobis(colMeans(savor_zip_hist_iv_10_90_mah_dist %>% filter(ZIP_HIST_10_90 == 1) %>% select(c(5:83,85:102))),
                                                      colMeans(savor_zip_hist_iv_10_90_mah_dist %>% filter(ZIP_HIST_10_90 == 0) %>% select(c(5:83,85:102))),
                                                      cov(savor_zip_hist_iv_10_90_mah_dist %>% select(c(5:83,85:102))),
                                                      tol = 1e-30)

((savor_mah_dist_total_zip_hist_iv_10_90-savor_mah_dist_total_14_17_treatment)/(savor_mah_dist_total_14_17_treatment))*100

# Distribution of outcomes by ZCTA exposure and by actual exposure
savor_iv_zip_hist_10_90_table <- filter(savor_phys, is.na(ZIP_HIST_10_90) == F)
table(savor_iv_zip_hist_10_90_table$OUTCOME,
      savor_iv_zip_hist_10_90_table$EXPOSURE,
      savor_iv_zip_hist_10_90_table$ZIP_HIST_10_90)

# ZCTA proportion saxagliptin ≥80% vs. ≤20%
savor_phys <- savor_phys %>%
  mutate(ZIP_HIST_20_80 = case_when(
    prop_exp_hist >= 0.8 ~ "1",
    prop_exp_hist <= 0.2 ~ "0",
    TRUE ~ NA_character_
  )
  )

# Check strength of instrument - linear regression
savor_model_iv_zip_hist_strength_20_80 <- glm(data = savor_phys,
                                              formula = EXPOSURE ~ ZIP_HIST_20_80 + as.factor(PS_ZIP_QUINTILE_14_17),
                                              na.action = na.exclude)
summary(savor_model_iv_zip_hist_strength_20_80)
linearHypothesis(savor_model_iv_zip_hist_strength_20_80, "ZIP_HIST_20_801 = 0", test = "F")

# Ordinary least squares - no adjustment for confounders
savor_model_zip_hist_ols_20_80 <- glm(data = filter(savor_phys, (prop_exp_hist <= 0.2 | prop_exp_hist >= 0.8)),
                                      formula = OUTCOME ~ EXPOSURE,
                                      na.action = na.exclude)
summary(savor_model_zip_hist_ols_20_80)

# 2 stage least squares - no adjustment for confounders
savor_model_zip_hist_tsls_20_80 <- ivreg(OUTCOME ~ EXPOSURE | ZIP_HIST_20_80,
                                         data = savor_phys)
summary(savor_model_zip_hist_tsls_20_80, vcov = sandwich, diagnostics = TRUE)

# Ordinary least squares - adjusted for PS
savor_model_zip_hist_ols_20_80_adj <- glm(data = filter(savor_phys, (prop_exp_hist <= 0.2 | prop_exp_hist >= 0.8)),
                                          formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                          na.action = na.exclude)
summary(savor_model_zip_hist_ols_20_80_adj)

# 2 stage least squares - adjusted for PS
savor_model_zip_hist_tsls_20_80_adj <- ivreg(OUTCOME ~ as.factor(PS_ZIP_QUINTILE_14_17) | EXPOSURE | ZIP_HIST_20_80,
                                             data = savor_phys)
summary(savor_model_zip_hist_tsls_20_80_adj, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - no adjustment for confounders
savor_model_zip_hist_logistic_20_80 <- glm(data = filter(savor_phys, (prop_exp_hist <= 0.2 | prop_exp_hist >= 0.8)),
                                           formula = OUTCOME ~ EXPOSURE,
                                           family = binomial(link = "logit"),
                                           na.action = na.exclude)
summary(savor_model_zip_hist_logistic_20_80)
exp(savor_model_zip_hist_logistic_20_80$coefficients)

# 2 stage residual inclusion - no adjustment for confounders
savor_model_zip_hist_tsri_20_80 <- tsri(OUTCOME ~ EXPOSURE | ZIP_HIST_20_80,
                                        data = savor_phys,
                                        link = "logit")
summary(savor_model_zip_hist_tsri_20_80, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - adjusted for PS
savor_model_zip_hist_logistic_20_80_adj <- glm(data = filter(savor_phys, (prop_exp_hist <= 0.2 | prop_exp_hist >= 0.8)),
                                               formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                               family = binomial(link = "logit"),
                                               na.action = na.exclude)
summary(savor_model_zip_hist_logistic_20_80_adj)
exp(savor_model_zip_hist_logistic_20_80_adj$coefficients)

# 2 stage residual inclusion - adjusted for PS
savor_model_zip_hist_tsri_20_80_adj <- tsri(OUTCOME ~ PS_ZIP_QUINTILE_14_17 | EXPOSURE | ZIP_HIST_20_80,
                                            data = savor_phys,
                                            link = "logit")
summary(savor_model_zip_hist_tsri_20_80_adj, vcov = sandwich, diagnostics = TRUE)

# Table 1 at the individual level
savor_tableone_all_iv_zip_hist_20_80 <- CreateTableOne(vars = savor_vars_ps, data = filter(savor_phys, is.na(ZIP_HIST_20_80) == FALSE),
                                                       factorVars = savor_factor_vars_ps,
                                                       strata = "ZIP_HIST_20_80", test = FALSE, smd = TRUE)
savor_tableone_all_iv_zip_hist_20_80_export <- print(savor_tableone_all_iv_zip_hist_20_80, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(savor_tableone_all_iv_zip_hist_20_80_export, file = "/Volumes/PACS$/Cordes, Jack/savor_tableone_all_iv_zip_hist_20_80.csv")

# Check balance improvement by percent change in Mahalanobis distance using treatment
savor_zip_hist_iv_20_80_mah_dist <- savor_phys[,c(1:13,19:21,23,34:35,38,40:41,48,52:53,55,58:60,66,70,80:81,87,89,91:92,100,102,105:106,108,113,117:120,122,130,135:136,141:142,144:150,152,158,164,167:169,171:173,175:177,179:182,187:190,206,210:211,220,226:229,231,237,244:245,247,249,254:258,262:264,271,293)]

savor_mah_dist_total_zip_hist_iv_20_80 <- mahalanobis(colMeans(savor_zip_hist_iv_20_80_mah_dist %>% filter(ZIP_HIST_20_80 == 1) %>% select(c(5:83,85:102))),
                                                      colMeans(savor_zip_hist_iv_20_80_mah_dist %>% filter(ZIP_HIST_20_80 == 0) %>% select(c(5:83,85:102))),
                                                      cov(savor_zip_hist_iv_20_80_mah_dist %>% select(c(5:83,85:102))),
                                                      tol = 1e-30)

((savor_mah_dist_total_zip_hist_iv_20_80-savor_mah_dist_total_14_17_treatment)/(savor_mah_dist_total_14_17_treatment))*100

# Distribution of outcomes by ZCTA exposure and by actual exposure
savor_iv_zip_hist_20_80_table <- filter(savor_phys, is.na(ZIP_HIST_20_80) == F)
table(savor_iv_zip_hist_20_80_table$OUTCOME,
      savor_iv_zip_hist_20_80_table$EXPOSURE,
      savor_iv_zip_hist_20_80_table$ZIP_HIST_20_80)

# ZCTA proportion saxagliptin ≥70% vs. ≤30%
savor_phys <- savor_phys %>%
  mutate(ZIP_HIST_30_70 = case_when(
    prop_exp_hist >= 0.7 ~ "1",
    prop_exp_hist <= 0.3 ~ "0",
    TRUE ~ NA_character_
  )
  )

# Check strength of instrument - linear regression
savor_model_iv_zip_hist_strength_30_70 <- glm(data = savor_phys,
                                              formula = EXPOSURE ~ ZIP_HIST_30_70 + as.factor(PS_ZIP_QUINTILE_14_17),
                                              na.action = na.exclude)
summary(savor_model_iv_zip_hist_strength_30_70)
linearHypothesis(savor_model_iv_zip_hist_strength_30_70, "ZIP_HIST_30_701 = 0", test = "F")

# Ordinary least squares - no adjustment for confounders
savor_model_zip_hist_ols_30_70 <- glm(data = filter(savor_phys, (prop_exp_hist <= 0.3 | prop_exp_hist >= 0.7)),
                                      formula = OUTCOME ~ EXPOSURE,
                                      na.action = na.exclude)
summary(savor_model_zip_hist_ols_30_70)

# 2 stage least squares - no adjustment for confounders
savor_model_zip_hist_tsls_30_70 <- ivreg(OUTCOME ~ EXPOSURE | ZIP_HIST_30_70,
                                         data = savor_phys)
summary(savor_model_zip_hist_tsls_30_70, vcov = sandwich, diagnostics = TRUE)

# Ordinary least squares - adjusted for PS
savor_model_zip_hist_ols_30_70_adj <- glm(data = filter(savor_phys, (prop_exp_hist <= 0.3 | prop_exp_hist >= 0.7)),
                                          formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                          na.action = na.exclude)
summary(savor_model_zip_hist_ols_30_70_adj)

# 2 stage least squares - adjusted for PS
savor_model_zip_hist_tsls_30_70_adj <- ivreg(OUTCOME ~ as.factor(PS_ZIP_QUINTILE_14_17) | EXPOSURE | ZIP_HIST_30_70,
                                             data = savor_phys)
summary(savor_model_zip_hist_tsls_30_70_adj, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - no adjustment for confounders
savor_model_zip_hist_logistic_30_70 <- glm(data = filter(savor_phys, (prop_exp_hist <= 0.3 | prop_exp_hist >= 0.7)),
                                           formula = OUTCOME ~ EXPOSURE,
                                           family = binomial(link = "logit"),
                                           na.action = na.exclude)
summary(savor_model_zip_hist_logistic_30_70)
exp(savor_model_zip_hist_logistic_30_70$coefficients)

# 2 stage residual inclusion - no adjustment for confounders
savor_model_zip_hist_tsri_30_70 <- tsri(OUTCOME ~ EXPOSURE | ZIP_HIST_30_70,
                                        data = savor_phys,
                                        link = "logit")
summary(savor_model_zip_hist_tsri_30_70, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - adjusted for PS
savor_model_zip_hist_logistic_30_70_adj <- glm(data = filter(savor_phys, (prop_exp_hist <= 0.3 | prop_exp_hist >= 0.7)),
                                               formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                               family = binomial(link = "logit"),
                                               na.action = na.exclude)
summary(savor_model_zip_hist_logistic_30_70_adj)
exp(savor_model_zip_hist_logistic_30_70_adj$coefficients)

# 2 stage residual inclusion - adjusted for PS
savor_model_zip_hist_tsri_30_70_adj <- tsri(OUTCOME ~ PS_ZIP_QUINTILE_14_17 | EXPOSURE | ZIP_HIST_30_70,
                                            data = savor_phys,
                                            link = "logit")
summary(savor_model_zip_hist_tsri_30_70_adj, vcov = sandwich, diagnostics = TRUE)

# Table 1 at the individual level
savor_tableone_all_iv_zip_hist_30_70 <- CreateTableOne(vars = savor_vars_ps, data = filter(savor_phys, is.na(ZIP_HIST_30_70) == FALSE),
                                                       factorVars = savor_factor_vars_ps,
                                                       strata = "ZIP_HIST_30_70", test = FALSE, smd = TRUE)
savor_tableone_all_iv_zip_hist_30_70_export <- print(savor_tableone_all_iv_zip_hist_30_70, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(savor_tableone_all_iv_zip_hist_30_70_export, file = "/Volumes/PACS$/Cordes, Jack/savor_tableone_all_iv_zip_hist_30_70.csv")

# Check balance improvement by percent change in Mahalanobis distance using treatment
savor_zip_hist_iv_30_70_mah_dist <- savor_phys[,c(1:13,19:21,23,34:35,38,40:41,48,52:53,55,58:60,66,70,80:81,87,89,91:92,100,102,105:106,108,113,117:120,122,130,135:136,141:142,144:150,152,158,164,167:169,171:173,175:177,179:182,187:190,206,210:211,220,226:229,231,237,244:245,247,249,254:258,262:264,271,294)]

savor_mah_dist_total_zip_hist_iv_30_70 <- mahalanobis(colMeans(savor_zip_hist_iv_30_70_mah_dist %>% filter(ZIP_HIST_30_70 == 1) %>% select(c(5:83,85:102))),
                                                      colMeans(savor_zip_hist_iv_30_70_mah_dist %>% filter(ZIP_HIST_30_70 == 0) %>% select(c(5:83,85:102))),
                                                      cov(savor_zip_hist_iv_30_70_mah_dist %>% select(c(5:83,85:102))),
                                                      tol = 1e-30)

((savor_mah_dist_total_zip_hist_iv_30_70-savor_mah_dist_total_14_17_treatment)/(savor_mah_dist_total_14_17_treatment))*100

# Distribution of outcomes by ZCTA exposure and by actual exposure
savor_iv_zip_hist_30_70_table <- filter(savor_phys, is.na(ZIP_HIST_30_70) == F)
table(savor_iv_zip_hist_30_70_table$OUTCOME,
      savor_iv_zip_hist_30_70_table$EXPOSURE,
      savor_iv_zip_hist_30_70_table$ZIP_HIST_30_70)

# ZCTA proportion saxagliptin ≥60% vs. ≤40%
savor_phys <- savor_phys %>%
  mutate(ZIP_HIST_40_60 = case_when(
    prop_exp_hist >= 0.6 ~ "1",
    prop_exp_hist <= 0.4 ~ "0",
    TRUE ~ NA_character_
  )
  )

# Check strength of instrument - linear regression
savor_model_iv_zip_hist_strength_40_60 <- glm(data = savor_phys,
                                              formula = EXPOSURE ~ ZIP_HIST_40_60 + as.factor(PS_ZIP_QUINTILE_14_17),
                                              na.action = na.exclude)
summary(savor_model_iv_zip_hist_strength_40_60)
linearHypothesis(savor_model_iv_zip_hist_strength_40_60, "ZIP_HIST_40_601 = 0", test = "F")

# Ordinary least squares - no adjustment for confounders
savor_model_zip_hist_ols_40_60 <- glm(data = filter(savor_phys, (prop_exp_hist <= 0.4 | prop_exp_hist >= 0.6)),
                                      formula = OUTCOME ~ EXPOSURE,
                                      na.action = na.exclude)
summary(savor_model_zip_hist_ols_40_60)

# 2 stage least squares - no adjustment for confounders
savor_model_zip_hist_tsls_40_60 <- ivreg(OUTCOME ~ EXPOSURE | ZIP_HIST_40_60,
                                         data = savor_phys)
summary(savor_model_zip_hist_tsls_40_60, vcov = sandwich, diagnostics = TRUE)

# Ordinary least squares - adjusted for PS
savor_model_zip_hist_ols_40_60_adj <- glm(data = filter(savor_phys, (prop_exp_hist <= 0.4 | prop_exp_hist >= 0.6)),
                                          formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                          na.action = na.exclude)
summary(savor_model_zip_hist_ols_40_60_adj)

# 2 stage least squares - adjusted for PS
savor_model_zip_hist_tsls_40_60_adj <- ivreg(OUTCOME ~ as.factor(PS_ZIP_QUINTILE_14_17) | EXPOSURE | ZIP_HIST_40_60,
                                             data = savor_phys)
summary(savor_model_zip_hist_tsls_40_60_adj, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - no adjustment for confounders
savor_model_zip_hist_logistic_40_60 <- glm(data = filter(savor_phys, (prop_exp_hist <= 0.4 | prop_exp_hist >= 0.6)),
                                           formula = OUTCOME ~ EXPOSURE,
                                           family = binomial(link = "logit"),
                                           na.action = na.exclude)
summary(savor_model_zip_hist_logistic_40_60)
exp(savor_model_zip_hist_logistic_40_60$coefficients)

# 2 stage residual inclusion - no adjustment for confounders
savor_model_zip_hist_tsri_40_60 <- tsri(OUTCOME ~ EXPOSURE | ZIP_HIST_40_60,
                                        data = savor_phys,
                                        link = "logit")
summary(savor_model_zip_hist_tsri_40_60, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - adjusted for PS
savor_model_zip_hist_logistic_40_60_adj <- glm(data = filter(savor_phys, (prop_exp_hist <= 0.4 | prop_exp_hist >= 0.6)),
                                               formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                               family = binomial(link = "logit"),
                                               na.action = na.exclude)
summary(savor_model_zip_hist_logistic_40_60_adj)
exp(savor_model_zip_hist_logistic_40_60_adj$coefficients)

# 2 stage residual inclusion - adjusted for PS
savor_model_zip_hist_tsri_40_60_adj <- tsri(OUTCOME ~ PS_ZIP_QUINTILE_14_17 | EXPOSURE | ZIP_HIST_40_60,
                                            data = savor_phys,
                                            link = "logit")
summary(savor_model_zip_hist_tsri_40_60_adj, vcov = sandwich, diagnostics = TRUE)

# Table 1 at the individual level
savor_tableone_all_iv_zip_hist_40_60 <- CreateTableOne(vars = savor_vars_ps, data = filter(savor_phys, is.na(ZIP_HIST_40_60) == FALSE),
                                                       factorVars = savor_factor_vars_ps,
                                                       strata = "ZIP_HIST_40_60", test = FALSE, smd = TRUE)
savor_tableone_all_iv_zip_hist_40_60_export <- print(savor_tableone_all_iv_zip_hist_40_60, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(savor_tableone_all_iv_zip_hist_40_60_export, file = "/Volumes/PACS$/Cordes, Jack/savor_tableone_all_iv_zip_hist_40_60.csv")

# Check balance improvement by percent change in Mahalanobis distance using treatment
savor_zip_hist_iv_40_60_mah_dist <- savor_phys[,c(1:13,19:21,23,34:35,38,40:41,48,52:53,55,58:60,66,70,80:81,87,89,91:92,100,102,105:106,108,113,117:120,122,130,135:136,141:142,144:150,152,158,164,167:169,171:173,175:177,179:182,187:190,206,210:211,220,226:229,231,237,244:245,247,249,254:258,262:264,271,295)]

savor_mah_dist_total_zip_hist_iv_40_60 <- mahalanobis(colMeans(savor_zip_hist_iv_40_60_mah_dist %>% filter(ZIP_HIST_40_60 == 1) %>% select(c(5:83,85:102))),
                                                      colMeans(savor_zip_hist_iv_40_60_mah_dist %>% filter(ZIP_HIST_40_60 == 0) %>% select(c(5:83,85:102))),
                                                      cov(savor_zip_hist_iv_40_60_mah_dist %>% select(c(5:83,85:102))),
                                                      tol = 1e-30)

((savor_mah_dist_total_zip_hist_iv_40_60-savor_mah_dist_total_14_17_treatment)/(savor_mah_dist_total_14_17_treatment))*100

# Distribution of outcomes by ZCTA exposure and by actual exposure
savor_iv_zip_hist_40_60_table <- filter(savor_phys, is.na(ZIP_HIST_40_60) == F)
table(savor_iv_zip_hist_40_60_table$OUTCOME,
      savor_iv_zip_hist_40_60_table$EXPOSURE,
      savor_iv_zip_hist_40_60_table$ZIP_HIST_40_60)

# ZCTA proportion saxagliptin ≥50% vs. <50%
savor_phys <- savor_phys %>%
  mutate(ZIP_HIST_50_50 = case_when(
    prop_exp_hist >= 0.5 ~ "1",
    prop_exp_hist < 0.5 ~ "0",
    TRUE ~ NA_character_
  )
  )

# Check strength of instrument - linear regression
savor_model_iv_zip_hist_strength_50_50 <- glm(data = savor_phys,
                                              formula = EXPOSURE ~ ZIP_HIST_50_50 + as.factor(PS_ZIP_QUINTILE_14_17),
                                              na.action = na.exclude)
summary(savor_model_iv_zip_hist_strength_50_50)
linearHypothesis(savor_model_iv_zip_hist_strength_50_50, "ZIP_HIST_50_501 = 0", test = "F")

# Ordinary least squares - no adjustment for confounders
savor_model_zip_hist_ols_50_50 <- glm(data = filter(savor_phys, (prop_exp_hist < 0.5 | prop_exp_hist >= 0.5)),
                                      formula = OUTCOME ~ EXPOSURE,
                                      na.action = na.exclude)
summary(savor_model_zip_hist_ols_50_50)

# 2 stage least squares - no adjustment for confounders
savor_model_zip_hist_tsls_50_50 <- ivreg(OUTCOME ~ EXPOSURE | ZIP_HIST_50_50,
                                         data = savor_phys)
summary(savor_model_zip_hist_tsls_50_50, vcov = sandwich, diagnostics = TRUE)

# Ordinary least squares - adjusted for PS
savor_model_zip_hist_ols_50_50_adj <- glm(data = filter(savor_phys, (prop_exp_hist < 0.5 | prop_exp_hist >= 0.5)),
                                          formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                          na.action = na.exclude)
summary(savor_model_zip_hist_ols_50_50_adj)

# 2 stage least squares - adjusted for PS
savor_model_zip_hist_tsls_50_50_adj <- ivreg(OUTCOME ~ as.factor(PS_ZIP_QUINTILE_14_17) | EXPOSURE | ZIP_HIST_50_50,
                                             data = savor_phys)
summary(savor_model_zip_hist_tsls_50_50_adj, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - no adjustment for confounders
savor_model_zip_hist_logistic_50_50 <- glm(data = filter(savor_phys, (prop_exp_hist < 0.5 | prop_exp_hist >= 0.5)),
                                           formula = OUTCOME ~ EXPOSURE,
                                           family = binomial(link = "logit"),
                                           na.action = na.exclude)
summary(savor_model_zip_hist_logistic_50_50)
exp(savor_model_zip_hist_logistic_50_50$coefficients)

# 2 stage residual inclusion - no adjustment for confounders
savor_model_zip_hist_tsri_50_50 <- tsri(OUTCOME ~ EXPOSURE | ZIP_HIST_50_50,
                                        data = savor_phys,
                                        link = "logit")
summary(savor_model_zip_hist_tsri_50_50, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - adjusted for PS
savor_model_zip_hist_logistic_50_50_adj <- glm(data = filter(savor_phys, (prop_exp_hist < 0.5 | prop_exp_hist >= 0.5)),
                                               formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                               family = binomial(link = "logit"),
                                               na.action = na.exclude)
summary(savor_model_zip_hist_logistic_50_50_adj)
exp(savor_model_zip_hist_logistic_50_50_adj$coefficients)

# 2 stage residual inclusion - adjusted for PS
savor_model_zip_hist_tsri_50_50_adj <- tsri(OUTCOME ~ PS_ZIP_QUINTILE_14_17 | EXPOSURE | ZIP_HIST_50_50,
                                            data = savor_phys,
                                            link = "logit")
summary(savor_model_zip_hist_tsri_50_50_adj, vcov = sandwich, diagnostics = TRUE)

# Table 1 at the individual level
savor_tableone_all_iv_zip_hist_50_50 <- CreateTableOne(vars = savor_vars_ps, data = filter(savor_phys, is.na(ZIP_HIST_50_50) == FALSE),
                                                       factorVars = savor_factor_vars_ps,
                                                       strata = "ZIP_HIST_50_50", test = FALSE, smd = TRUE)
savor_tableone_all_iv_zip_hist_50_50_export <- print(savor_tableone_all_iv_zip_hist_50_50, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(savor_tableone_all_iv_zip_hist_50_50_export, file = "/Volumes/PACS$/Cordes, Jack/savor_tableone_all_iv_zip_hist_50_50.csv")

# Check balance improvement by percent change in Mahalanobis distance using treatment
savor_zip_hist_iv_50_50_mah_dist <- savor_phys[,c(1:13,19:21,23,34:35,38,40:41,48,52:53,55,58:60,66,70,80:81,87,89,91:92,100,102,105:106,108,113,117:120,122,130,135:136,141:142,144:150,152,158,164,167:169,171:173,175:177,179:182,187:190,206,210:211,220,226:229,231,237,244:245,247,249,254:258,262:264,271,296)]

savor_mah_dist_total_zip_hist_iv_50_50 <- mahalanobis(colMeans(savor_zip_hist_iv_50_50_mah_dist %>% filter(ZIP_HIST_50_50 == 1) %>% select(c(5:83,85:102))),
                                                      colMeans(savor_zip_hist_iv_50_50_mah_dist %>% filter(ZIP_HIST_50_50 == 0) %>% select(c(5:83,85:102))),
                                                      cov(savor_zip_hist_iv_50_50_mah_dist %>% select(c(5:83,85:102))),
                                                      tol = 1e-30)

((savor_mah_dist_total_zip_hist_iv_50_50-savor_mah_dist_total_14_17_treatment)/(savor_mah_dist_total_14_17_treatment))*100

# Distribution of outcomes by ZCTA exposure and by actual exposure
savor_iv_zip_hist_50_50_table <- filter(savor_phys, is.na(ZIP_HIST_50_50) == F)
table(savor_iv_zip_hist_50_50_table$OUTCOME,
      savor_iv_zip_hist_50_50_table$EXPOSURE,
      savor_iv_zip_hist_50_50_table$ZIP_HIST_50_50)

# IV using physician preference
# Distribution of physicians across ZCTAs
savor_zcta_phys_frequency <- savor_phys %>%
  group_by(ZCTA_USE) %>%
  dplyr::summarize(count_zcta_phys = n_distinct(prscrbr_id))

summary(savor_zcta_phys_frequency$count_zcta_phys)

# Add prescriber ZCTA counts to individuals
savor_phys <- left_join(savor_phys, savor_zcta_phys_frequency, by = "ZCTA_USE")

# Distribution of prescribers across ZCTAs with at least 5 individuals
savor_phys_5 <- filter(savor_phys, count > 4)

# Filter dataset to only ZCTAs
savor_phys_5_zcta <- distinct(savor_phys_5, ZCTA_USE, .keep_all = TRUE)
summary(savor_phys_5_zcta$count_zcta_phys)

savor_zcta_phys_dist_plot <- savor_phys_5_zcta %>%
  ggplot() + 
  geom_histogram(aes(count_zcta_phys), binwidth = 5, color = "black", fill = "white") + 
  geom_vline(xintercept = 17.66, color = "red") +
  scale_y_continuous(labels = label_comma()) +
  coord_cartesian(xlim = c(0, 180)) +
  labs(title = "Saxagliptin Cohort Prescriber Distribution across ZCTAs\n5+ Individuals in ZCTA", x = "Count of Prescribers", y = "Count of ZCTAs") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(hjust = 0.5))
savor_zcta_phys_dist_grob <- gtable_add_grob(ggplot_gtable(ggplot_build(savor_zcta_phys_dist_plot)), b_label, t = 1, l = 4, b = 6)
ggsave("savor_zcta_phys_dist.pdf", savor_zcta_phys_dist_grob)

# Spatial distribution of prescribers by ZCTA with at least 5 cohort individuals in ZCTA
# Join frequency data to ZCTA polygons
zcta_shp_savor <- merge(zcta_shp_savor, savor_phys_5_zcta, by = "ZCTA_USE")

# Create Lower 48 ZCTA and state shapefiles
zcta_shp_savor_48 <- st_transform(zcta_shp_savor %>% filter(str_starts(ZCTA_USE, "0") |
                                                              str_starts(ZCTA_USE, "1") |
                                                              str_starts(ZCTA_USE, "2") |
                                                              str_starts(ZCTA_USE, "3") |
                                                              str_starts(ZCTA_USE, "4") |
                                                              str_starts(ZCTA_USE, "5") |
                                                              str_starts(ZCTA_USE, "6") |
                                                              str_starts(ZCTA_USE, "7") |
                                                              str_starts(ZCTA_USE, "8") |
                                                              str_starts(ZCTA_USE, "90") |
                                                              str_starts(ZCTA_USE, "91") |
                                                              str_starts(ZCTA_USE, "92") |
                                                              str_starts(ZCTA_USE, "93") |
                                                              str_starts(ZCTA_USE, "94") |
                                                              str_starts(ZCTA_USE, "95") |
                                                              str_starts(ZCTA_USE, "961") |
                                                              str_starts(ZCTA_USE, "97") |
                                                              str_starts(ZCTA_USE, "980") |
                                                              str_starts(ZCTA_USE, "981") |
                                                              str_starts(ZCTA_USE, "982") |
                                                              str_starts(ZCTA_USE, "983") |
                                                              str_starts(ZCTA_USE, "984") |
                                                              str_starts(ZCTA_USE, "985") |
                                                              str_starts(ZCTA_USE, "986") |
                                                              str_starts(ZCTA_USE, "987") |
                                                              str_starts(ZCTA_USE, "988") |
                                                              str_starts(ZCTA_USE, "989") |
                                                              str_starts(ZCTA_USE, "990") |
                                                              str_starts(ZCTA_USE, "991") |
                                                              str_starts(ZCTA_USE, "992") |
                                                              str_starts(ZCTA_USE, "993") |
                                                              str_starts(ZCTA_USE, "994")), 5070)
state_shp_48 <- st_transform(state_shp %>% filter(str_starts(STATEFP, "01") |
                                                    str_starts(STATEFP, "04") |
                                                    str_starts(STATEFP, "05") |
                                                    str_starts(STATEFP, "06") |
                                                    str_starts(STATEFP, "08") |
                                                    str_starts(STATEFP, "09") |
                                                    str_starts(STATEFP, "10") |
                                                    str_starts(STATEFP, "11") |
                                                    str_starts(STATEFP, "12") |
                                                    str_starts(STATEFP, "13") |
                                                    str_starts(STATEFP, "16") |
                                                    str_starts(STATEFP, "17") |
                                                    str_starts(STATEFP, "18") |
                                                    str_starts(STATEFP, "19") |
                                                    str_starts(STATEFP, "20") |
                                                    str_starts(STATEFP, "21") |
                                                    str_starts(STATEFP, "22") |
                                                    str_starts(STATEFP, "23") |
                                                    str_starts(STATEFP, "24") |
                                                    str_starts(STATEFP, "25") |
                                                    str_starts(STATEFP, "26") |
                                                    str_starts(STATEFP, "27") |
                                                    str_starts(STATEFP, "28") |
                                                    str_starts(STATEFP, "29") |
                                                    str_starts(STATEFP, "30") |
                                                    str_starts(STATEFP, "31") |
                                                    str_starts(STATEFP, "32") |
                                                    str_starts(STATEFP, "33") |
                                                    str_starts(STATEFP, "34") |
                                                    str_starts(STATEFP, "35") |
                                                    str_starts(STATEFP, "36") |
                                                    str_starts(STATEFP, "37") |
                                                    str_starts(STATEFP, "38") |
                                                    str_starts(STATEFP, "39") |
                                                    str_starts(STATEFP, "40") |
                                                    str_starts(STATEFP, "41") |
                                                    str_starts(STATEFP, "42") |
                                                    str_starts(STATEFP, "44") |
                                                    str_starts(STATEFP, "45") |
                                                    str_starts(STATEFP, "46") |
                                                    str_starts(STATEFP, "47") |
                                                    str_starts(STATEFP, "48") |
                                                    str_starts(STATEFP, "49") |
                                                    str_starts(STATEFP, "50") |
                                                    str_starts(STATEFP, "51") |
                                                    str_starts(STATEFP, "53") |
                                                    str_starts(STATEFP, "54") |
                                                    str_starts(STATEFP, "55") |
                                                    str_starts(STATEFP, "56")), 5070)
state_shp_48 <- st_buffer(state_shp_48, dist = 0)

# Create Alaska ZCTA and state shapefiles
zcta_shp_savor_AK <- st_transform(zcta_shp_savor %>% filter(str_starts(ZCTA_USE, "995") | 
                                                              str_starts(ZCTA_USE, "996") |
                                                              str_starts(ZCTA_USE, "997") |
                                                              str_starts(ZCTA_USE, "998") |
                                                              str_starts(ZCTA_USE, "999")), 3338)
state_shp_AK <- st_transform(state_shp %>% filter(str_starts(STATEFP, "02")), 3338)
state_shp_AK <- st_buffer(state_shp_AK, dist = 0)

# Create Hawaii ZCTA and state shapefiles
zcta_shp_savor_HI <- st_transform(zcta_shp_savor %>% filter(str_starts(ZCTA_USE, "967") | 
                                                              str_starts(ZCTA_USE, "968")), 32603)
state_shp_HI <- st_transform(state_shp %>% filter(str_starts(STATEFP, "15")), 32603)
state_shp_HI <- st_buffer(state_shp_HI, dist = 0)

# Map of number of prescribers by ZCTA
#options(device = "X11")

zcta_map_savor_48_phys <- ggplot() +
  geom_sf(data = zcta_shp_savor_48, aes(fill = cut(count_zcta_phys, c(0,1,5,10,25,180))), color = NA, show.legend = TRUE) +
  scale_fill_brewer(palette = "Spectral", name = "Saxagliptin Cohort\nCount of Prescribers\nper ZCTA", direction = -1) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + 
  geom_sf(data = state_shp_48, fill = NA, color = "black", size = 0.01)
#zcta_map_savor_48_phys

zcta_map_savor_AK_phys <- ggplot() +
  geom_sf(data = zcta_shp_savor_AK, aes(fill = cut(count_zcta_phys, c(0,1,5,10,25,180))), color = NA, show.legend = FALSE) +
  scale_fill_brewer(palette = "Spectral", direction = -1) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_blank()) + 
  geom_sf(data = st_crop(state_shp_AK, xmin = -1000000,
                         xmax = 1491822,
                         ymin = 414200,
                         ymax = 2378425), fill = NA, color = "black", size = 0.01)
#zcta_map_savor_AK_phys

zcta_map_savor_HI_phys <- ggplot() +
  geom_sf(data = zcta_shp_savor_HI, aes(fill = cut(count_zcta_phys, c(0,1,5,10,25,180))), color = NA, show.legend = FALSE) +
  scale_fill_brewer(palette = "Spectral", direction = -1) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.background = element_blank()) +
  geom_sf(data = st_crop(state_shp_HI, xmin = 1000000,
                         xmax = 1573922,
                         ymin = 2117066,
                         ymax = 3214762), fill = NA, color = "black", size = 0.01)
#zcta_map_savor_HI_phys

zcta_map_savor_phys <- ggdraw() +
  draw_plot(zcta_map_savor_48_phys) +
  draw_plot(zcta_map_savor_AK_phys, x = 0, y = 0, width = 0.3, height = 0.28) +
  draw_plot(zcta_map_savor_HI_phys, x = 0.4, y = 0, width = 0.3, height = 0.18) +
  draw_grob(b_label, x = 0.03, y = 0.15)
#zcta_map_savor_phys
ggsave("/Volumes/PACS$/Cordes, Jack/savor_zcta_map_phys.pdf", zcta_map_savor_phys)

# Instantaneous physician preference IV analysis
# Convert PS_ZIP_QUINTILE_14_17 to factor for tsri
savor_phys$PS_ZIP_QUINTILE_14_17 <- as.factor(savor_phys$PS_ZIP_QUINTILE_14_17)

# Sort on prescriber ID and CED
savor_phys <- savor_phys[order(savor_phys$prscrbr_id, savor_phys$ENTRYDATE),]

# Create empty instantaneous physician preference column
savor_phys <- savor_phys %>% mutate(PHYS_PREF_INST = NA)

# Create instantaneous physician preference variable based on most recent prescription in cohort
for (i in 2:nrow(savor_phys)){
  if (savor_phys$prscrbr_id[i] == savor_phys$prscrbr_id[i-1]){
    savor_phys$PHYS_PREF_INST[i] <- savor_phys$EXPOSURE[i-1]
  }
}

# Distribution of global historical physician preference proportions
summary(savor_phys$PHYS_PREF_INST)

# Check strength of instrument - linear regression
savor_model_iv_strength_phys_pref_inst <- glm(data = savor_phys,
                                              formula = EXPOSURE ~ PHYS_PREF_INST + as.factor(PS_ZIP_QUINTILE_14_17),
                                              na.action = na.exclude)
summary(savor_model_iv_strength_phys_pref_inst)
linearHypothesis(savor_model_iv_strength_phys_pref_inst, "PHYS_PREF_INST = 0", test = "F")

# Ordinary least squares - no adjustment for confounders
savor_model_ols_phys_pref_inst <- glm(data = filter(savor_phys, is.na(PHYS_PREF_INST) == FALSE),
                                      formula = OUTCOME ~ EXPOSURE,
                                      na.action = na.exclude)
summary(savor_model_ols_phys_pref_inst)

# 2 stage least squares - no adjustment for confounders
savor_model_tsls_phys_pref_inst <- ivreg(OUTCOME ~ EXPOSURE | PHYS_PREF_INST,
                                         data = filter(savor_phys, is.na(PHYS_PREF_INST) == FALSE))
summary(savor_model_tsls_phys_pref_inst, vcov = sandwich, diagnostics = TRUE)

# Ordinary least squares - adjusted for PS
savor_model_ols_phys_pref_inst_adj <- glm(data = filter(savor_phys, is.na(PHYS_PREF_INST) == FALSE),
                                          formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                          na.action = na.exclude)
summary(savor_model_ols_phys_pref_inst_adj)

# 2 stage least squares - adjusted for PS
savor_model_tsls_phys_pref_inst_adj <- ivreg(OUTCOME ~ as.factor(PS_ZIP_QUINTILE_14_17) | EXPOSURE | PHYS_PREF_INST,
                                             data = filter(savor_phys, is.na(PHYS_PREF_INST) == FALSE))
summary(savor_model_tsls_phys_pref_inst_adj, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - no adjustment for confounders
savor_model_logistic_phys_pref_inst <- glm(data = filter(savor_phys, is.na(PHYS_PREF_INST) == FALSE),
                                           formula = OUTCOME ~ EXPOSURE,
                                           family = binomial(link = "logit"),
                                           na.action = na.exclude)
summary(savor_model_logistic_phys_pref_inst)
exp(savor_model_logistic_phys_pref_inst$coefficients)

# 2 stage residual inclusion - no adjustment for confounders
savor_model_tsri_phys_pref_inst <- tsri(OUTCOME ~ EXPOSURE | PHYS_PREF_INST,
                                        data = filter(savor_phys, is.na(PHYS_PREF_INST) == FALSE),
                                        link = "logit")
summary(savor_model_tsri_phys_pref_inst, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - adjusted for PS
savor_model_logistic_phys_pref_inst_adj <- glm(data = filter(savor_phys, is.na(PHYS_PREF_INST) == FALSE),
                                               formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                               family = binomial(link = "logit"),
                                               na.action = na.exclude)
summary(savor_model_logistic_phys_pref_inst_adj)
exp(savor_model_logistic_phys_pref_inst_adj$coefficients)

# 2 stage residual inclusion - adjusted for PS
savor_model_tsri_phys_pref_inst_adj <- tsri(OUTCOME ~ PS_ZIP_QUINTILE_14_17 | EXPOSURE | PHYS_PREF_INST,
                                            data = filter(savor_phys, is.na(PHYS_PREF_INST) == FALSE),
                                            link = "logit")
summary(savor_model_tsri_phys_pref_inst_adj, vcov = sandwich, diagnostics = TRUE)

# Table 1 at the individual level
savor_tableone_all_iv_phys_pref_inst <- CreateTableOne(vars = savor_vars_ps, data = filter(savor_phys, is.na(PHYS_PREF_INST) == FALSE),
                                                       factorVars = savor_factor_vars_ps,
                                                       strata = "PHYS_PREF_INST", test = FALSE, smd = TRUE)
savor_tableone_all_iv_phys_pref_inst_export <- print(savor_tableone_all_iv_phys_pref_inst, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(savor_tableone_all_iv_phys_pref_inst_export, file = "/Volumes/PACS$/Cordes, Jack/savor_tableone_all_iv_phys_pref_inst.csv")

# Check balance improvement by percent change in Mahalanobis distance using treatment
savor_phys_iv_inst_mah_dist <- savor_phys[,c(1:13,19:21,23,34:35,38,40:41,48,52:53,55,58:60,66,70,80:81,87,89,91:92,100,102,105:106,108,113,117:120,122,130,135:136,141:142,144:150,152,158,164,167:169,171:173,175:177,179:182,187:190,206,210:211,220,226:229,231,237,244:245,247,249,254:258,262:264,298)]

savor_mah_dist_total_phys_iv_inst <- mahalanobis(colMeans(savor_phys_iv_inst_mah_dist %>% filter(PHYS_PREF_INST == 1) %>% select(c(5:83,85:102))),
                                                 colMeans(savor_phys_iv_inst_mah_dist %>% filter(PHYS_PREF_INST == 0) %>% select(c(5:83,85:102))),
                                                 cov(savor_phys_iv_inst_mah_dist %>% select(c(5:83,85:102))),
                                                 tol = 1e-30)

((savor_mah_dist_total_phys_iv_inst-savor_mah_dist_total_14_17_treatment)/(savor_mah_dist_total_14_17_treatment))*100

# Distribution of outcomes by physician preference and by actual exposure
savor_iv_phys_pref_inst_table <- filter(savor_phys, is.na(PHYS_PREF_INST) == FALSE)
table(savor_iv_phys_pref_inst_table$OUTCOME,
      savor_iv_phys_pref_inst_table$EXPOSURE,
      savor_iv_phys_pref_inst_table$PHYS_PREF_INST)

# Global historical physician preference IV analysis
# Sort on prescriber ID and CED
savor_phys <- savor_phys[order(savor_phys$prscrbr_id, savor_phys$ENTRYDATE),]

# Create empty global historical physician preference column
savor_phys <- savor_phys %>% mutate(prop_exp_phys = NA)

# Create empty counter column
savor_phys <- savor_phys %>% mutate(counter = NA)

# Create empty cumulative sum column
savor_phys <- savor_phys %>% mutate(cum_sum = NA)

# Create global historical physician preference variable based on average of all previous prescriptions
for (i in 2:nrow(savor_phys)){
  if ((savor_phys$prscrbr_id[i] == savor_phys$prscrbr_id[i-1]) & (is.na(savor_phys$prop_exp_phys[i-1]) == TRUE)){
    savor_phys$prop_exp_phys[i] <- savor_phys$EXPOSURE[i-1]
  }
  if ((savor_phys$prscrbr_id[i] == savor_phys$prscrbr_id[i-1]) & (is.na(savor_phys$prop_exp_phys[i-1]) == TRUE)){
    savor_phys$counter[i] <- 1
  }
  if ((savor_phys$prscrbr_id[i] == savor_phys$prscrbr_id[i-1]) & (is.na(savor_phys$prop_exp_phys[i-1]) == TRUE)){
    savor_phys$cum_sum[i] <- savor_phys$EXPOSURE[i-1]
  }
  if ((savor_phys$prscrbr_id[i] == savor_phys$prscrbr_id[i-1]) & (is.na(savor_phys$prop_exp_phys[i-1]) == FALSE)){
    savor_phys$counter[i] <- savor_phys$counter[i-1]+1
  }
  if ((savor_phys$prscrbr_id[i] == savor_phys$prscrbr_id[i-1]) & (is.na(savor_phys$prop_exp_phys[i-1]) == FALSE)){
    savor_phys$cum_sum[i] <- savor_phys$cum_sum[i-1] + savor_phys$EXPOSURE[i-1]
  }
  if ((savor_phys$prscrbr_id[i] == savor_phys$prscrbr_id[i-1]) & (is.na(savor_phys$prop_exp_phys[i-1]) == FALSE)){
    savor_phys$prop_exp_phys[i] <- savor_phys$cum_sum[i]/savor_phys$counter[i]
  }
}

# Distribution of global historical physician preference proportions
summary(savor_phys$prop_exp_phys)

# Create dataset limited to non-missing physician preference IV
savor_phys_global <- filter(savor_phys, is.na(prop_exp_phys) == FALSE)

savor_phys_global_pref_dist_plot <- savor_phys_global %>%
  ggplot() + 
  geom_histogram(aes(prop_exp_phys), binwidth = 0.05, color = "black", fill = "white") + 
  geom_vline(xintercept = 0.04, color = "red") +
  scale_y_continuous(labels = label_comma()) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(title = "Saxagliptin Cohort\nCumulative Prescriber Preference for Each Individual ", x = "Proportion Saxagliptin", y = "Count of Individuals") +
  theme(plot.title = element_text(size = 18, hjust = 0.5),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(hjust = 0.5))
savor_phys_global_pref_dist_grob <- gtable_add_grob(ggplot_gtable(ggplot_build(savor_phys_global_pref_dist_plot)), b_label, t = 1, l = 4, b = 6)
ggsave("savor_phys_global_pref_dist.pdf", savor_phys_global_pref_dist_grob)

# For IV with proportion saxagliptin cutoffs 100% vs. 0%, instrument is collapsed with exposure
# because no variation in individual-level exposure
# So becomes simple regression comparing individuals with 100% to 0% prescribers

# Prescriber proportion saxagliptin 100% vs. 0%
savor_phys <- savor_phys %>%
  mutate(PHYS_PREF_GLOBAL_0_100 = case_when(
    prop_exp_phys == 1 ~ "1",
    prop_exp_phys == 0 ~ "0",
    TRUE ~ NA_character_
  )
  )

# Ordinary least squares - no adjustment for confounders
savor_model_phys_pref_global_ols_0_100 <- glm(data = filter(savor_phys, (prop_exp_phys == 0 | prop_exp_phys == 1)),
                                              formula = OUTCOME ~ EXPOSURE,
                                              na.action = na.exclude)
summary(savor_model_phys_pref_global_ols_0_100)

# Ordinary least squares - adjusted for PS
savor_model_phys_pref_global_ols_0_100_adj <- glm(data = filter(savor_phys, (prop_exp_phys == 0 | prop_exp_phys == 1)),
                                                  formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                                  na.action = na.exclude)
summary(savor_model_phys_pref_global_ols_0_100_adj)

# Logistic regression - no adjustment for confounders
savor_model_phys_pref_global_logistic_0_100 <- glm(data = filter(savor_phys, (prop_exp_phys == 0 | prop_exp_phys == 1)),
                                                   formula = OUTCOME ~ EXPOSURE,
                                                   family = binomial(link = "logit"),
                                                   na.action = na.exclude)
summary(savor_model_phys_pref_global_logistic_0_100)
exp(savor_model_phys_pref_global_logistic_0_100$coefficients)

# Logistic regression - adjusted for PS
savor_model_phys_pref_global_logistic_0_100_adj <- glm(data = filter(savor_phys, (prop_exp_phys == 0 | prop_exp_phys == 1)),
                                                       formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                                       family = binomial(link = "logit"),
                                                       na.action = na.exclude)
summary(savor_model_phys_pref_global_logistic_0_100_adj)
exp(savor_model_phys_pref_global_logistic_0_100_adj$coefficients)

# Table 1 at the individual level
savor_tableone_all_iv_phys_pref_global_0_100 <- CreateTableOne(vars = savor_vars_ps, data = filter(savor_phys, is.na(PHYS_PREF_GLOBAL_0_100) == FALSE),
                                                               factorVars = savor_factor_vars_ps,
                                                               strata = "PHYS_PREF_GLOBAL_0_100", test = FALSE, smd = TRUE)
savor_tableone_all_iv_phys_pref_global_0_100_export <- print(savor_tableone_all_iv_phys_pref_global_0_100, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(savor_tableone_all_iv_phys_pref_global_0_100_export, file = "/Volumes/PACS$/Cordes, Jack/savor_tableone_all_iv_phys_pref_global_0_100.csv")

# Check balance improvement by percent change in Mahalanobis distance using treatment
savor_phys_iv_0_100_mah_dist <- savor_phys[,c(1:13,19:21,23,34:35,38,40:41,48,52:53,55,58:60,66,70,80:81,87,89,91:92,100,102,105:106,108,113,117:120,122,130,135:136,141:142,144:150,152,158,164,167:169,171:173,175:177,179:182,187:190,206,210:211,220,226:229,231,237,244:245,247,249,254:258,262:264,300)]

savor_mah_dist_total_phys_iv_0_100 <- mahalanobis(colMeans(savor_phys_iv_0_100_mah_dist %>% filter(PHYS_PREF_GLOBAL_0_100 == 1) %>% select(c(5:83,85:102))),
                                                  colMeans(savor_phys_iv_0_100_mah_dist %>% filter(PHYS_PREF_GLOBAL_0_100 == 0) %>% select(c(5:83,85:102))),
                                                  cov(savor_phys_iv_0_100_mah_dist %>% select(c(5:83,85:102))),
                                                  tol = 1e-30)

((savor_mah_dist_total_phys_iv_0_100-savor_mah_dist_total_14_17_treatment)/(savor_mah_dist_total_14_17_treatment))*100

# Distribution of outcomes by zipcode exposure and by actual exposure
savor_iv_phys_pref_global_0_100_table <- filter(savor_phys, is.na(PHYS_PREF_GLOBAL_0_100) == F)
table(savor_iv_phys_pref_global_0_100_table$OUTCOME,
      savor_iv_phys_pref_global_0_100_table$EXPOSURE,
      savor_iv_phys_pref_global_0_100_table$ZCTA_EXP_0_100)

# Prescriber proportion saxagliptin ≥90% vs. ≤10%
savor_phys <- savor_phys %>%
  mutate(PHYS_PREF_GLOBAL_10_90 = case_when(
    prop_exp_phys >= 0.9 ~ "1",
    prop_exp_phys <= 0.1 ~ "0",
    TRUE ~ NA_character_
  )
  )

# Check strength of instrument - linear regression
savor_model_iv_phys_pref_global_strength_10_90 <- glm(data = savor_phys,
                                                      formula = EXPOSURE ~ PHYS_PREF_GLOBAL_10_90 + as.factor(PS_ZIP_QUINTILE_14_17),
                                                      na.action = na.exclude)
summary(savor_model_iv_phys_pref_global_strength_10_90)
linearHypothesis(savor_model_iv_phys_pref_global_strength_10_90, "PHYS_PREF_GLOBAL_10_901 = 0", test = "F")

# Ordinary least squares - no adjustment for confounders
savor_model_phys_pref_global_ols_10_90 <- glm(data = filter(savor_phys, (prop_exp_phys <= 0.1 | prop_exp_phys >= 0.9)),
                                              formula = OUTCOME ~ EXPOSURE,
                                              na.action = na.exclude)
summary(savor_model_phys_pref_global_ols_10_90)

# 2 stage least squares - no adjustment for confounders
savor_model_phys_pref_global_tsls_10_90 <- ivreg(OUTCOME ~ EXPOSURE | PHYS_PREF_GLOBAL_10_90,
                                                 data = savor_phys)
summary(savor_model_phys_pref_global_tsls_10_90, vcov = sandwich, diagnostics = TRUE)

# Ordinary least squares - adjusted for PS
savor_model_phys_pref_global_ols_10_90_adj <- glm(data = filter(savor_phys, (prop_exp_phys <= 0.1 | prop_exp_phys >= 0.9)),
                                                  formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                                  na.action = na.exclude)
summary(savor_model_phys_pref_global_ols_10_90_adj)

# 2 stage least squares - adjusted for PS
savor_model_phys_pref_global_tsls_10_90_adj <- ivreg(OUTCOME ~ as.factor(PS_ZIP_QUINTILE_14_17) | EXPOSURE | PHYS_PREF_GLOBAL_10_90,
                                                     data = savor_phys)
summary(savor_model_phys_pref_global_tsls_10_90_adj, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - no adjustment for confounders
savor_model_phys_pref_global_logistic_10_90 <- glm(data = filter(savor_phys, (prop_exp_phys <= 0.1 | prop_exp_phys >= 0.9)),
                                                   formula = OUTCOME ~ EXPOSURE,
                                                   family = binomial(link = "logit"),
                                                   na.action = na.exclude)
summary(savor_model_phys_pref_global_logistic_10_90)
exp(savor_model_phys_pref_global_logistic_10_90$coefficients)

# 2 stage residual inclusion - no adjustment for confounders
savor_model_phys_pref_global_tsri_10_90 <- tsri(OUTCOME ~ EXPOSURE | PHYS_PREF_GLOBAL_10_90,
                                                data = savor_phys,
                                                link = "logit")
summary(savor_model_phys_pref_global_tsri_10_90, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - adjusted for PS
savor_model_phys_pref_global_logistic_10_90_adj <- glm(data = filter(savor_phys, (prop_exp_phys <= 0.1 | prop_exp_phys >= 0.9)),
                                                       formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                                       family = binomial(link = "logit"),
                                                       na.action = na.exclude)
summary(savor_model_phys_pref_global_logistic_10_90_adj)
exp(savor_model_phys_pref_global_logistic_10_90_adj$coefficients)

# 2 stage residual inclusion - adjusted for PS
savor_model_phys_pref_global_tsri_10_90_adj <- tsri(OUTCOME ~ PS_ZIP_QUINTILE_14_17 | EXPOSURE | PHYS_PREF_GLOBAL_10_90,
                                                    data = savor_phys,
                                                    link = "logit")
summary(savor_model_phys_pref_global_tsri_10_90_adj, vcov = sandwich, diagnostics = TRUE)

# Table 1 at the individual level
savor_tableone_all_iv_phys_pref_global_10_90 <- CreateTableOne(vars = savor_vars_ps, data = filter(savor_phys, is.na(PHYS_PREF_GLOBAL_10_90) == FALSE),
                                                               factorVars = savor_factor_vars_ps,
                                                               strata = "PHYS_PREF_GLOBAL_10_90", test = FALSE, smd = TRUE)
savor_tableone_all_iv_phys_pref_global_10_90_export <- print(savor_tableone_all_iv_phys_pref_global_10_90, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(savor_tableone_all_iv_phys_pref_global_10_90_export, file = "/Volumes/PACS$/Cordes, Jack/savor_tableone_all_iv_phys_pref_global_10_90.csv")

# Check balance improvement by percent change in Mahalanobis distance using treatment
savor_phys_iv_10_90_mah_dist <- savor_phys[,c(1:13,19:21,23,34:35,38,40:41,48,52:53,55,58:60,66,70,80:81,87,89,91:92,100,102,105:106,108,113,117:120,122,130,135:136,141:142,144:150,152,158,164,167:169,171:173,175:177,179:182,187:190,206,210:211,220,226:229,231,237,244:245,247,249,254:258,262:264,301)]

savor_mah_dist_total_phys_iv_10_90 <- mahalanobis(colMeans(savor_phys_iv_10_90_mah_dist %>% filter(PHYS_PREF_GLOBAL_10_90 == 1) %>% select(c(5:83,85:102))),
                                                  colMeans(savor_phys_iv_10_90_mah_dist %>% filter(PHYS_PREF_GLOBAL_10_90 == 0) %>% select(c(5:83,85:102))),
                                                  cov(savor_phys_iv_10_90_mah_dist %>% select(c(5:83,85:102))),
                                                  tol = 1e-30)

((savor_mah_dist_total_phys_iv_10_90-savor_mah_dist_total_14_17_treatment)/(savor_mah_dist_total_14_17_treatment))*100

# Distribution of outcomes by prescriber exposure and by actual exposure
savor_iv_phys_pref_global_10_90_table <- filter(savor_phys, is.na(PHYS_PREF_GLOBAL_10_90) == F)
table(savor_iv_phys_pref_global_10_90_table$OUTCOME,
      savor_iv_phys_pref_global_10_90_table$EXPOSURE,
      savor_iv_phys_pref_global_10_90_table$PHYS_PREF_GLOBAL_10_90)

# Prescriber proportion saxagliptin ≥80% vs. ≤20%
savor_phys <- savor_phys %>%
  mutate(PHYS_PREF_GLOBAL_20_80 = case_when(
    prop_exp_phys >= 0.8 ~ "1",
    prop_exp_phys <= 0.2 ~ "0",
    TRUE ~ NA_character_
  )
  )

# Check strength of instrument - linear regression
savor_model_iv_phys_pref_global_strength_20_80 <- glm(data = savor_phys,
                                                      formula = EXPOSURE ~ PHYS_PREF_GLOBAL_20_80 + as.factor(PS_ZIP_QUINTILE_14_17),
                                                      na.action = na.exclude)
summary(savor_model_iv_phys_pref_global_strength_20_80)
linearHypothesis(savor_model_iv_phys_pref_global_strength_20_80, "PHYS_PREF_GLOBAL_20_801 = 0", test = "F")

# Ordinary least squares - no adjustment for confounders
savor_model_phys_pref_global_ols_20_80 <- glm(data = filter(savor_phys, (prop_exp_phys <= 0.2 | prop_exp_phys >= 0.8)),
                                              formula = OUTCOME ~ EXPOSURE,
                                              na.action = na.exclude)
summary(savor_model_phys_pref_global_ols_20_80)

# 2 stage least squares - no adjustment for confounders
savor_model_phys_pref_global_tsls_20_80 <- ivreg(OUTCOME ~ EXPOSURE | PHYS_PREF_GLOBAL_20_80,
                                                 data = savor_phys)
summary(savor_model_phys_pref_global_tsls_20_80, vcov = sandwich, diagnostics = TRUE)

# Ordinary least squares - adjusted for PS
savor_model_phys_pref_global_ols_20_80_adj <- glm(data = filter(savor_phys, (prop_exp_phys <= 0.2 | prop_exp_phys >= 0.8)),
                                                  formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                                  na.action = na.exclude)
summary(savor_model_phys_pref_global_ols_20_80_adj)

# 2 stage least squares - adjusted for PS
savor_model_phys_pref_global_tsls_20_80_adj <- ivreg(OUTCOME ~ as.factor(PS_ZIP_QUINTILE_14_17) | EXPOSURE | PHYS_PREF_GLOBAL_20_80,
                                                     data = savor_phys)
summary(savor_model_phys_pref_global_tsls_20_80_adj, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - no adjustment for confounders
savor_model_phys_pref_global_logistic_20_80 <- glm(data = filter(savor_phys, (prop_exp_phys <= 0.2 | prop_exp_phys >= 0.8)),
                                                   formula = OUTCOME ~ EXPOSURE,
                                                   family = binomial(link = "logit"),
                                                   na.action = na.exclude)
summary(savor_model_phys_pref_global_logistic_20_80)
exp(savor_model_phys_pref_global_logistic_20_80$coefficients)

# 2 stage residual inclusion - no adjustment for confounders
savor_model_phys_pref_global_tsri_20_80 <- tsri(OUTCOME ~ EXPOSURE | PHYS_PREF_GLOBAL_20_80,
                                                data = savor_phys,
                                                link = "logit")
summary(savor_model_phys_pref_global_tsri_20_80, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - adjusted for PS
savor_model_phys_pref_global_logistic_20_80_adj <- glm(data = filter(savor_phys, (prop_exp_phys <= 0.2 | prop_exp_phys >= 0.8)),
                                                       formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                                       family = binomial(link = "logit"),
                                                       na.action = na.exclude)
summary(savor_model_phys_pref_global_logistic_20_80_adj)
exp(savor_model_phys_pref_global_logistic_20_80_adj$coefficients)

# 2 stage residual inclusion - adjusted for PS
savor_model_phys_pref_global_tsri_20_80_adj <- tsri(OUTCOME ~ PS_ZIP_QUINTILE_14_17 | EXPOSURE | PHYS_PREF_GLOBAL_20_80,
                                                    data = savor_phys,
                                                    link = "logit")
summary(savor_model_phys_pref_global_tsri_20_80_adj, vcov = sandwich, diagnostics = TRUE)

# Table 1 at the individual level
savor_tableone_all_iv_phys_pref_global_20_80 <- CreateTableOne(vars = savor_vars_ps, data = filter(savor_phys, is.na(PHYS_PREF_GLOBAL_20_80) == FALSE),
                                                               factorVars = savor_factor_vars_ps,
                                                               strata = "PHYS_PREF_GLOBAL_20_80", test = FALSE, smd = TRUE)
savor_tableone_all_iv_phys_pref_global_20_80_export <- print(savor_tableone_all_iv_phys_pref_global_20_80, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(savor_tableone_all_iv_phys_pref_global_20_80_export, file = "/Volumes/PACS$/Cordes, Jack/savor_tableone_all_iv_phys_pref_global_20_80.csv")

# Check balance improvement by percent change in Mahalanobis distance using treatment
savor_phys_iv_20_80_mah_dist <- savor_phys[,c(1:13,19:21,23,34:35,38,40:41,48,52:53,55,58:60,66,70,80:81,87,89,91:92,100,102,105:106,108,113,117:120,122,130,135:136,141:142,144:150,152,158,164,167:169,171:173,175:177,179:182,187:190,206,210:211,220,226:229,231,237,244:245,247,249,254:258,262:264,302)]

savor_mah_dist_total_phys_iv_20_80 <- mahalanobis(colMeans(savor_phys_iv_20_80_mah_dist %>% filter(PHYS_PREF_GLOBAL_20_80 == 1) %>% select(c(5:83,85:102))),
                                                  colMeans(savor_phys_iv_20_80_mah_dist %>% filter(PHYS_PREF_GLOBAL_20_80 == 0) %>% select(c(5:83,85:102))),
                                                  cov(savor_phys_iv_20_80_mah_dist %>% select(c(5:83,85:102))),
                                                  tol = 1e-30)

((savor_mah_dist_total_phys_iv_20_80-savor_mah_dist_total_14_17_treatment)/(savor_mah_dist_total_14_17_treatment))*100

# Distribution of outcomes by prescriber exposure and by actual exposure
savor_iv_phys_pref_global_20_80_table <- filter(savor_phys, is.na(PHYS_PREF_GLOBAL_20_80) == F)
table(savor_iv_phys_pref_global_20_80_table$OUTCOME,
      savor_iv_phys_pref_global_20_80_table$EXPOSURE,
      savor_iv_phys_pref_global_20_80_table$PHYS_PREF_GLOBAL_20_80)

# Prescriber proportion saxagliptin ≥70% vs. ≤30%
savor_phys <- savor_phys %>%
  mutate(PHYS_PREF_GLOBAL_30_70 = case_when(
    prop_exp_phys >= 0.7 ~ "1",
    prop_exp_phys <= 0.3 ~ "0",
    TRUE ~ NA_character_
  )
  )

# Check strength of instrument - linear regression
savor_model_iv_phys_pref_global_strength_30_70 <- glm(data = savor_phys,
                                                      formula = EXPOSURE ~ PHYS_PREF_GLOBAL_30_70 + as.factor(PS_ZIP_QUINTILE_14_17),
                                                      na.action = na.exclude)
summary(savor_model_iv_phys_pref_global_strength_30_70)
linearHypothesis(savor_model_iv_phys_pref_global_strength_30_70, "PHYS_PREF_GLOBAL_30_701 = 0", test = "F")

# Ordinary least squares - no adjustment for confounders
savor_model_phys_pref_global_ols_30_70 <- glm(data = filter(savor_phys, (prop_exp_phys <= 0.3 | prop_exp_phys >= 0.7)),
                                              formula = OUTCOME ~ EXPOSURE,
                                              na.action = na.exclude)
summary(savor_model_phys_pref_global_ols_30_70)

# 2 stage least squares - no adjustment for confounders
savor_model_phys_pref_global_tsls_30_70 <- ivreg(OUTCOME ~ EXPOSURE | PHYS_PREF_GLOBAL_30_70,
                                                 data = savor_phys)
summary(savor_model_phys_pref_global_tsls_30_70, vcov = sandwich, diagnostics = TRUE)

# Ordinary least squares - adjusted for PS
savor_model_phys_pref_global_ols_30_70_adj <- glm(data = filter(savor_phys, (prop_exp_phys <= 0.3 | prop_exp_phys >= 0.7)),
                                                  formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                                  na.action = na.exclude)
summary(savor_model_phys_pref_global_ols_30_70_adj)

# 2 stage least squares - adjusted for PS
savor_model_phys_pref_global_tsls_30_70_adj <- ivreg(OUTCOME ~ as.factor(PS_ZIP_QUINTILE_14_17) | EXPOSURE | PHYS_PREF_GLOBAL_30_70,
                                                     data = savor_phys)
summary(savor_model_phys_pref_global_tsls_30_70_adj, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - no adjustment for confounders
savor_model_phys_pref_global_logistic_30_70 <- glm(data = filter(savor_phys, (prop_exp_phys <= 0.3 | prop_exp_phys >= 0.7)),
                                                   formula = OUTCOME ~ EXPOSURE,
                                                   family = binomial(link = "logit"),
                                                   na.action = na.exclude)
summary(savor_model_phys_pref_global_logistic_30_70)
exp(savor_model_phys_pref_global_logistic_30_70$coefficients)

# 2 stage residual inclusion - no adjustment for confounders
savor_model_phys_pref_global_tsri_30_70 <- tsri(OUTCOME ~ EXPOSURE | PHYS_PREF_GLOBAL_30_70,
                                                data = savor_phys,
                                                link = "logit")
summary(savor_model_phys_pref_global_tsri_30_70, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - adjusted for PS
savor_model_phys_pref_global_logistic_30_70_adj <- glm(data = filter(savor_phys, (prop_exp_phys <= 0.3 | prop_exp_phys >= 0.7)),
                                                       formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                                       family = binomial(link = "logit"),
                                                       na.action = na.exclude)
summary(savor_model_phys_pref_global_logistic_30_70_adj)
exp(savor_model_phys_pref_global_logistic_30_70_adj$coefficients)

# 2 stage residual inclusion - prescribers adjusted for PS
savor_model_phys_pref_global_tsri_30_70_adj <- tsri(OUTCOME ~ PS_ZIP_QUINTILE_14_17 | EXPOSURE | PHYS_PREF_GLOBAL_30_70,
                                                    data = savor_phys,
                                                    link = "logit")
summary(savor_model_phys_pref_global_tsri_30_70_adj, vcov = sandwich, diagnostics = TRUE)

# Table 1 at the individual level
savor_tableone_all_iv_phys_pref_global_30_70 <- CreateTableOne(vars = savor_vars_ps, data = filter(savor_phys, is.na(PHYS_PREF_GLOBAL_30_70) == FALSE),
                                                               factorVars = savor_factor_vars_ps,
                                                               strata = "PHYS_PREF_GLOBAL_30_70", test = FALSE, smd = TRUE)
savor_tableone_all_iv_phys_pref_global_30_70_export <- print(savor_tableone_all_iv_phys_pref_global_30_70, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(savor_tableone_all_iv_phys_pref_global_30_70_export, file = "/Volumes/PACS$/Cordes, Jack/savor_tableone_all_iv_phys_pref_global_30_70.csv")

# Check balance improvement by percent change in Mahalanobis distance using treatment
savor_phys_iv_30_70_mah_dist <- savor_phys[,c(1:13,19:21,23,34:35,38,40:41,48,52:53,55,58:60,66,70,80:81,87,89,91:92,100,102,105:106,108,113,117:120,122,130,135:136,141:142,144:150,152,158,164,167:169,171:173,175:177,179:182,187:190,206,210:211,220,226:229,231,237,244:245,247,249,254:258,262:264,303)]

savor_mah_dist_total_phys_iv_30_70 <- mahalanobis(colMeans(savor_phys_iv_30_70_mah_dist %>% filter(PHYS_PREF_GLOBAL_30_70 == 1) %>% select(c(5:83,85:102))),
                                                  colMeans(savor_phys_iv_30_70_mah_dist %>% filter(PHYS_PREF_GLOBAL_30_70 == 0) %>% select(c(5:83,85:102))),
                                                  cov(savor_phys_iv_30_70_mah_dist %>% select(c(5:83,85:102))),
                                                  tol = 1e-30)

((savor_mah_dist_total_phys_iv_30_70-savor_mah_dist_total_14_17_treatment)/(savor_mah_dist_total_14_17_treatment))*100

# Distribution of outcomes by prescriber exposure and by actual exposure
savor_iv_phys_pref_global_30_70_table <- filter(savor_phys, is.na(PHYS_PREF_GLOBAL_30_70) == F)
table(savor_iv_phys_pref_global_30_70_table$OUTCOME,
      savor_iv_phys_pref_global_30_70_table$EXPOSURE,
      savor_iv_phys_pref_global_30_70_table$PHYS_PREF_GLOBAL_30_70)

# Prescriber proportion saxagliptin ≥60% vs. ≤40%
savor_phys <- savor_phys %>%
  mutate(PHYS_PREF_GLOBAL_40_60 = case_when(
    prop_exp_phys >= 0.6 ~ "1",
    prop_exp_phys <= 0.4 ~ "0",
    TRUE ~ NA_character_
  )
  )

# Check strength of instrument - linear regression
savor_model_iv_phys_pref_global_strength_40_60 <- glm(data = savor_phys,
                                                      formula = EXPOSURE ~ PHYS_PREF_GLOBAL_40_60 + as.factor(PS_ZIP_QUINTILE_14_17),
                                                      na.action = na.exclude)
summary(savor_model_iv_phys_pref_global_strength_40_60)
linearHypothesis(savor_model_iv_phys_pref_global_strength_40_60, "PHYS_PREF_GLOBAL_40_601 = 0", test = "F")

# Ordinary least squares - no adjustment for confounders
savor_model_phys_pref_global_ols_40_60 <- glm(data = filter(savor_phys, (prop_exp_phys <= 0.4 | prop_exp_phys >= 0.6)),
                                              formula = OUTCOME ~ EXPOSURE,
                                              na.action = na.exclude)
summary(savor_model_phys_pref_global_ols_40_60)

# 2 stage least squares - no adjustment for confounders
savor_model_phys_pref_global_tsls_40_60 <- ivreg(OUTCOME ~ EXPOSURE | PHYS_PREF_GLOBAL_40_60,
                                                 data = savor_phys)
summary(savor_model_phys_pref_global_tsls_40_60, vcov = sandwich, diagnostics = TRUE)

# Ordinary least squares - adjusted for PS
savor_model_phys_pref_global_ols_40_60_adj <- glm(data = filter(savor_phys, (prop_exp_phys <= 0.4 | prop_exp_phys >= 0.6)),
                                                  formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                                  na.action = na.exclude)
summary(savor_model_phys_pref_global_ols_40_60_adj)

# 2 stage least squares - adjusted for PS
savor_model_phys_pref_global_tsls_40_60_adj <- ivreg(OUTCOME ~ as.factor(PS_ZIP_QUINTILE_14_17) | EXPOSURE | PHYS_PREF_GLOBAL_40_60,
                                                     data = savor_phys)
summary(savor_model_phys_pref_global_tsls_40_60_adj, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - no adjustment for confounders
savor_model_phys_pref_global_logistic_40_60 <- glm(data = filter(savor_phys, (prop_exp_phys <= 0.4 | prop_exp_phys >= 0.6)),
                                                   formula = OUTCOME ~ EXPOSURE,
                                                   family = binomial(link = "logit"),
                                                   na.action = na.exclude)
summary(savor_model_phys_pref_global_logistic_40_60)
exp(savor_model_phys_pref_global_logistic_40_60$coefficients)

# 2 stage residual inclusion - no adjustment for confounders
savor_model_phys_pref_global_tsri_40_60 <- tsri(OUTCOME ~ EXPOSURE | PHYS_PREF_GLOBAL_40_60,
                                                data = savor_phys,
                                                link = "logit")
summary(savor_model_phys_pref_global_tsri_40_60, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - adjusted for PS
savor_model_phys_pref_global_logistic_40_60_adj <- glm(data = filter(savor_phys, (prop_exp_phys <= 0.4 | prop_exp_phys >= 0.6)),
                                                       formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                                       family = binomial(link = "logit"),
                                                       na.action = na.exclude)
summary(savor_model_phys_pref_global_logistic_40_60_adj)
exp(savor_model_phys_pref_global_logistic_40_60_adj$coefficients)

# 2 stage residual inclusion - adjusted for PS
savor_model_phys_pref_global_tsri_40_60_adj <- tsri(OUTCOME ~ PS_ZIP_QUINTILE_14_17 | EXPOSURE | PHYS_PREF_GLOBAL_40_60,
                                                    data = savor_phys,
                                                    link = "logit")
summary(savor_model_phys_pref_global_tsri_40_60_adj, vcov = sandwich, diagnostics = TRUE)

# Table 1 at the individual level
savor_tableone_all_iv_phys_pref_global_40_60 <- CreateTableOne(vars = savor_vars_ps, data = filter(savor_phys, is.na(PHYS_PREF_GLOBAL_40_60) == FALSE),
                                                               factorVars = savor_factor_vars_ps,
                                                               strata = "PHYS_PREF_GLOBAL_40_60", test = FALSE, smd = TRUE)
savor_tableone_all_iv_phys_pref_global_40_60_export <- print(savor_tableone_all_iv_phys_pref_global_40_60, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(savor_tableone_all_iv_phys_pref_global_40_60_export, file = "/Volumes/PACS$/Cordes, Jack/savor_tableone_all_iv_phys_pref_global_40_60.csv")

# Check balance improvement by percent change in Mahalanobis distance using treatment
savor_phys_iv_40_60_mah_dist <- savor_phys[,c(1:13,19:21,23,34:35,38,40:41,48,52:53,55,58:60,66,70,80:81,87,89,91:92,100,102,105:106,108,113,117:120,122,130,135:136,141:142,144:150,152,158,164,167:169,171:173,175:177,179:182,187:190,206,210:211,220,226:229,231,237,244:245,247,249,254:258,262:264,304)]

savor_mah_dist_total_phys_iv_40_60 <- mahalanobis(colMeans(savor_phys_iv_40_60_mah_dist %>% filter(PHYS_PREF_GLOBAL_40_60 == 1) %>% select(c(5:83,85:102))),
                                                  colMeans(savor_phys_iv_40_60_mah_dist %>% filter(PHYS_PREF_GLOBAL_40_60 == 0) %>% select(c(5:83,85:102))),
                                                  cov(savor_phys_iv_40_60_mah_dist %>% select(c(5:83,85:102))),
                                                  tol = 1e-30)

((savor_mah_dist_total_phys_iv_40_60-savor_mah_dist_total_14_17_treatment)/(savor_mah_dist_total_14_17_treatment))*100

# Distribution of outcomes by prescriber exposure and by actual exposure
savor_iv_phys_pref_global_40_60_table <- filter(savor_phys, is.na(PHYS_PREF_GLOBAL_40_60) == F)
table(savor_iv_phys_pref_global_40_60_table$OUTCOME,
      savor_iv_phys_pref_global_40_60_table$EXPOSURE,
      savor_iv_phys_pref_global_40_60_table$PHYS_PREF_GLOBAL_40_60)

# Prescriber proportion saxagliptin ≥50% vs. <50%
savor_phys <- savor_phys %>%
  mutate(PHYS_PREF_GLOBAL_50_50 = case_when(
    prop_exp_phys >= 0.5 ~ "1",
    prop_exp_phys < 0.5 ~ "0",
    TRUE ~ NA_character_
  )
  )

# Check strength of instrument - linear regression
savor_model_iv_phys_pref_global_strength_50_50 <- glm(data = savor_phys,
                                                      formula = EXPOSURE ~ PHYS_PREF_GLOBAL_50_50 + as.factor(PS_ZIP_QUINTILE_14_17),
                                                      na.action = na.exclude)
summary(savor_model_iv_phys_pref_global_strength_50_50)
linearHypothesis(savor_model_iv_phys_pref_global_strength_50_50, "PHYS_PREF_GLOBAL_50_501 = 0", test = "F")

# Ordinary least squares - no adjustment for confounders
savor_model_phys_pref_global_ols_50_50 <- glm(data = filter(savor_phys, (prop_exp_phys < 0.5 | prop_exp_phys >= 0.5)),
                                              formula = OUTCOME ~ EXPOSURE,
                                              na.action = na.exclude)
summary(savor_model_phys_pref_global_ols_50_50)

# 2 stage least squares - no adjustment for confounders
savor_model_phys_pref_global_tsls_50_50 <- ivreg(OUTCOME ~ EXPOSURE | PHYS_PREF_GLOBAL_50_50,
                                                 data = savor_phys)
summary(savor_model_phys_pref_global_tsls_50_50, vcov = sandwich, diagnostics = TRUE)

# Ordinary least squares - adjusted for PS
savor_model_phys_pref_global_ols_50_50_adj <- glm(data = filter(savor_phys, (prop_exp_phys < 0.5 | prop_exp_phys >= 0.5)),
                                                  formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                                  na.action = na.exclude)
summary(savor_model_phys_pref_global_ols_50_50_adj)

# 2 stage least squares - adjusted for PS
savor_model_phys_pref_global_tsls_50_50_adj <- ivreg(OUTCOME ~ as.factor(PS_ZIP_QUINTILE_14_17) | EXPOSURE | PHYS_PREF_GLOBAL_50_50,
                                                     data = savor_phys)
summary(savor_model_phys_pref_global_tsls_50_50_adj, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - no adjustment for confounders
savor_model_phys_pref_global_logistic_50_50 <- glm(data = filter(savor_phys, (prop_exp_phys < 0.5 | prop_exp_phys >= 0.5)),
                                                   formula = OUTCOME ~ EXPOSURE,
                                                   family = binomial(link = "logit"),
                                                   na.action = na.exclude)
summary(savor_model_phys_pref_global_logistic_50_50)
exp(savor_model_phys_pref_global_logistic_50_50$coefficients)

# 2 stage residual inclusion - no adjustment for confounders
savor_model_phys_pref_global_tsri_50_50 <- tsri(OUTCOME ~ EXPOSURE | PHYS_PREF_GLOBAL_50_50,
                                                data = savor_phys,
                                                link = "logit")
summary(savor_model_phys_pref_global_tsri_50_50, vcov = sandwich, diagnostics = TRUE)

# Logistic regression - adjusted for PS
savor_model_phys_pref_global_logistic_50_50_adj <- glm(data = filter(savor_phys, (prop_exp_phys < 0.5 | prop_exp_phys >= 0.5)),
                                                       formula = OUTCOME ~ EXPOSURE + as.factor(PS_ZIP_QUINTILE_14_17),
                                                       family = binomial(link = "logit"),
                                                       na.action = na.exclude)
summary(savor_model_phys_pref_global_logistic_50_50_adj)
exp(savor_model_phys_pref_global_logistic_50_50_adj$coefficients)

# 2 stage residual inclusion - adjusted for PS
savor_model_phys_pref_global_tsri_50_50_adj <- tsri(OUTCOME ~ PS_ZIP_QUINTILE_14_17 | EXPOSURE | PHYS_PREF_GLOBAL_50_50,
                                                    data = savor_phys,
                                                    link = "logit")
summary(savor_model_phys_pref_global_tsri_50_50_adj, vcov = sandwich, diagnostics = TRUE)

# Table 1 at the individual level
savor_tableone_all_iv_phys_pref_global_50_50 <- CreateTableOne(vars = savor_vars_ps, data = filter(savor_phys, is.na(PHYS_PREF_GLOBAL_50_50) == FALSE),
                                                               factorVars = savor_factor_vars_ps,
                                                               strata = "PHYS_PREF_GLOBAL_50_50", test = FALSE, smd = TRUE)
savor_tableone_all_iv_phys_pref_global_50_50_export <- print(savor_tableone_all_iv_phys_pref_global_50_50, smd = TRUE, quote = TRUE, noSpaces = TRUE)
write.csv(savor_tableone_all_iv_phys_pref_global_50_50_export, file = "/Volumes/PACS$/Cordes, Jack/savor_tableone_all_iv_phys_pref_global_50_50.csv")

# Check balance improvement by percent change in Mahalanobis distance using treatment
savor_phys_iv_50_50_mah_dist <- savor_phys[,c(1:13,19:21,23,34:35,38,40:41,48,52:53,55,58:60,66,70,80:81,87,89,91:92,100,102,105:106,108,113,117:120,122,130,135:136,141:142,144:150,152,158,164,167:169,171:173,175:177,179:182,187:190,206,210:211,220,226:229,231,237,244:245,247,249,254:258,262:264,305)]

savor_mah_dist_total_phys_iv_50_50 <- mahalanobis(colMeans(savor_phys_iv_50_50_mah_dist %>% filter(PHYS_PREF_GLOBAL_50_50 == 1) %>% select(c(5:83,85:102))),
                                                  colMeans(savor_phys_iv_50_50_mah_dist %>% filter(PHYS_PREF_GLOBAL_50_50 == 0) %>% select(c(5:83,85:102))),
                                                  cov(savor_phys_iv_50_50_mah_dist %>% select(c(5:83,85:102))),
                                                  tol = 1e-30)

((savor_mah_dist_total_phys_iv_50_50-savor_mah_dist_total_14_17_treatment)/(savor_mah_dist_total_14_17_treatment))*100

# Distribution of outcomes by prescriber exposure and by actual exposure
savor_iv_phys_pref_global_50_50_table <- filter(savor_phys, is.na(PHYS_PREF_GLOBAL_50_50) == F)
table(savor_iv_phys_pref_global_50_50_table$OUTCOME,
      savor_iv_phys_pref_global_50_50_table$EXPOSURE,
      savor_iv_phys_pref_global_50_50_table$PHYS_PREF_GLOBAL_50_50)
