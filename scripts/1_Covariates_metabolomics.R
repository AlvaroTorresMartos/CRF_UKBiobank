
## Load packages ------
if (!require("pacman")) install.packages("pacman")

pacman::p_load(dplyr, readr, descr, stringr, tidyr)

## Import data -----
t1 = Sys.time()

demo = read_csv("/mnt/project/data/raw/demographics.csv")

hes = read_csv("/mnt/project/data/raw/health_outcomes.csv")

bio = read_csv("/mnt/project/data/raw/biomarkers_281125.csv")

diet = read_csv("/mnt/project/data/raw/diet.csv")

cardio = read_csv("/mnt/project/data/raw/cardiorespiratory.csv")

t2 = Sys.time()
t1-t2
#Time difference of 2.781869 mins


t1 = Sys.time()
# Rename Participant ID to eid in diet
diet_rename <- diet %>% 
  rename(eid = `Participant ID`)

# Helper: merge evitando columnas duplicadas (excepto eid)
merge_no_dups <- function(x, y, by = "eid") {
  dup_cols <- setdiff(intersect(names(x), names(y)), by)
  if (length(dup_cols) > 0) {
    y <- y[, !(names(y) %in% dup_cols)]
  }
  merge(x, y, by = by)
}

## Merge datasets sin generar _x / _y ------
data_merged <- demo %>%
  merge_no_dups(hes, by = "eid") %>%
  merge_no_dups(bio, by = "eid") %>%
  merge_no_dups(diet_rename, by = "eid") %>%
  merge_no_dups(cardio, by = "eid")

# Remove original datasets
rm(bio, demo, hes, diet, diet_rename, cardio)

# Standardize column names
data_merged <- data_merged %>%
  rename_with(~ .x %>%
                str_replace_all("-", "_") %>%
                str_replace_all("\\.", "_"))

# Save
saveRDS(data_merged, 
        file = "./upload/data_merged_17022026.rds")

t2 = Sys.time()
t1-t2
# Time difference of -18.58748 mins



#PART2: Covariates-----
#FROM THIS POINT MAKE A COPY OF DATA_MERGED.RDS AND OPERATE WITH THE CODES PER BLOCK
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, readr, stringr, tidyr, descr)

# Load merged data
data_merged_copy <- readRDS("/mnt/project/data/processed/data_merged_17022026.rds")

## 21022: Age at recruitment----
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=21022

data_merged_copy <- data_merged_copy %>%
  dplyr::mutate(age_group = case_when(
    `21022_0_0` >= 39 & `21022_0_0` < 45 ~ 1, 
    `21022_0_0` >= 45 & `21022_0_0` < 50 ~ 2, 
    `21022_0_0` >= 50 & `21022_0_0` < 55 ~ 3, 
    `21022_0_0` >= 55 & `21022_0_0` < 60 ~ 4, 
    `21022_0_0` >= 60 & `21022_0_0` < 65 ~ 5, 
    `21022_0_0` >= 65 & `21022_0_0` < 70 ~ 6, 
    `21022_0_0` >= 70                ~ 7,
    is.na(`21022_0_0`) ~ 9
  ))

data_merged_copy$age_group %>% freq()
#      Frequency   Percent
#1         51704   10.3010
#2         65959   13.1410
#3         76215   15.1843
#4         90705   18.0711
#5        121379   24.1823
#6         93550   18.6379
#7          2422   0.4825
#NA's          2               
#Total    501936   100.0000

## 21001: Body mass index (BMI) -----
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=21001

data_merged_copy <- data_merged_copy %>%
  dplyr::mutate(bmi_group = case_when(
    `21001_0_0` < 18.5 ~ 0,
    `21001_0_0` >= 18.5 & `21001_0_0` < 25 ~ 1,
    `21001_0_0` >= 25 & `21001_0_0` < 30 ~ 2,
    `21001_0_0` >= 30 & `21001_0_0` < 35 ~ 3,
    `21001_0_0` >= 35 & `21001_0_0` < 40 ~ 4,
    `21001_0_0` >= 40 ~ 5, 
    is.na(`21001_0_0`) ~ 9
  ))

data_merged_copy$bmi_group %>% freq()
#      Frequency  Percent
#0          2624   0.5228
#1        162184  32.3117
#2        211888  42.2141
#3         87475  17.4275
#4         24966   4.9739
#5          9696   1.9317
#9          3103   0.6182
#Total    501936 100.0000      

## 21000: Ethnic background -----
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?tk=MozhW45lvgR3ND07YxPpY0fVTVY7dwBZ483974&id=21000
data_merged_copy = data_merged_copy %>%
  dplyr::rename(ethnicity0 = `21000_0_0`, ethnicity1 = `21000_1_0`,
                ethnicity2 = `21000_2_0`, ethnicity3 = `21000_3_0`) %>% 
  dplyr::mutate(across(starts_with("ethnicity"), ~ case_when(
    . %in% c("1", "1001", "1002", "1003") ~ 1,
    . %in% c("4", "4001", "4002", "4003") ~ 2, 
    . %in% c("3", "3001", "3002", "3003", "3004", "5") ~ 3,
    . %in% c("2", "2001",  "2002",  "2003", "2004", "6") ~ 4,
    . %in% c("-3", "-1") ~ 9,  
    TRUE ~ .))) %>% 
  dplyr::mutate(ethnicity = coalesce(ethnicity0, ethnicity1, 
                                     ethnicity2, ethnicity3), 
                ethnicity = if_else(is.na(ethnicity), 9, ethnicity)
  )

data_merged_copy$ethnicity %>% freq()
#Frequency  Percent
#1        472188  94.0733 #white
#2          8048   1.6034 #black
#3         11440   2.2792 #asian
#4          7501   1.4944 #mixed
#9          2759   0.5497 #unknown
#Total    501936 100.0000


## 22189: Townsend index of deprivation -----
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?tk=MozhW45lvgR3ND07YxPpY0fVTVY7dwBZ483974&id=22189
data_merged_copy = data_merged_copy %>%
  dplyr::mutate(townsend = `22189_0_0`, 
                # Create quintiles for Townsend (5 quantiles)
                townsend_q5 = ntile(townsend, 5),
                # Recode missing values to 9 (replace NA with 9)
                townsend_q5 = if_else(is.na(townsend_q5), 9, townsend_q5))

data_merged_copy$townsend_q5 %>% freq()
#      Frequency  Percent
#1        100263  19.9753
#2        100263  19.9753
#3        100263  19.9753
#4        100263  19.9753
#5        100263  19.9753
#9           621   0.1237
#Total    501936 100.0000



## 6138: Qualifications ------
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?tk=MozhW45lvgR3ND07YxPpY0fVTVY7dwBZ483974&id=6138
data_merged_copy = data_merged_copy %>%
  dplyr::mutate(`6138_0_0` = coalesce(`6138_0_0`, `6138_1_0`, `6138_2_0`, `6138_3_0`)) %>% 
  dplyr::rename(edu_group0 = `6138_0_0`) %>% 
  dplyr::mutate(across(starts_with("edu_group0"), 
                       ~ str_split(.x, "\\|"))) %>%
  unnest_wider(where(is.list), names_sep = "_") %>% 
  dplyr::mutate(education_group =  case_when(
    edu_group0_1 %in% c(1, 6) | edu_group0_2 %in% c(1, 6) | 
      edu_group0_3 %in% c(1, 6) | edu_group0_4 %in% c(1, 6) | 
      edu_group0_5 %in% c(1, 6) | edu_group0_6 %in% c(1, 6) ~ 1, #Group 1
    edu_group0_1 %in% c(2, 3, 4) | edu_group0_2 %in% c(2, 3, 4) | 
      edu_group0_3 %in% c(2, 3, 4) | edu_group0_4 %in% c(2, 3, 4) | 
      edu_group0_5 %in% c(2, 3, 4) | edu_group0_6 %in% c(2, 3, 4) ~ 2, #Group 2
    edu_group0_1 %in% c(5) | edu_group0_2 %in% c(5) | 
      edu_group0_3 %in% c(5) | edu_group0_4 %in% c(5) | 
      edu_group0_5 %in% c(5) | edu_group0_6 %in% c(5) ~ 3, #Group 3
    edu_group0_1 %in% c(-7) | edu_group0_2 %in% c(-7) | 
      edu_group0_3 %in% c(-7) | edu_group0_4 %in% c(-7) | 
      edu_group0_5 %in% c(-7) | edu_group0_6 %in% c(-7) ~ 4, #Group 4
    edu_group0_1 %in% c(-3) | edu_group0_2 %in% c(-3) | 
      edu_group0_3 %in% c(-3) | edu_group0_4 %in% c(-3) | 
      edu_group0_5 %in% c(-3) | edu_group0_6 %in% c(-3) ~ 9, #No answer or no data
    is.na(edu_group0_1) |  is.na(edu_group0_2) |  is.na(edu_group0_3) | 
      is.na(edu_group0_4) | is.na(edu_group0_5) |  is.na(edu_group0_6) ~ 9 #Group 9
  ))

data_merged_copy$education_group %>% freq()
#      Frequency Percent
#1        233996  46.619
#2        144215  28.732
#3         29444   5.866
#4         85301  16.994
#9          8980   1.789
#Total    501936 100.000

## 54: Geographical and location ----
#ID 54 BUT, encoded using ID 10 https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=10
# 1 "East Midlands" 
# 2 "London" 
# 3 "North East" 
# 4 "North West" 
# 5 "Scotland" 
# 6 "South East" 
# 7 "South West"
# 8 "Wales" 
# 9 "West Midlands" 
# 10 "York

data_merged_copy = data_merged_copy %>% 
  dplyr::mutate(region = case_when(
    `54_0_0` %in% c(11013) ~ 1, 
    `54_0_0` %in% c(11012, 11018, 11020) ~ 2, 
    `54_0_0` %in% c(11009, 11017) ~ 3, 
    `54_0_0` %in% c(10003, 11001, 11008, 11016) ~ 4, 
    `54_0_0` %in% c(11004, 11005) ~ 5, 
    `54_0_0` %in% c(11002, 11007) ~ 6, 
    `54_0_0` %in% c(11011) ~ 7, 
    `54_0_0` %in% c(11003, 11022, 11023) ~ 8, 
    `54_0_0` %in% c(11006, 11021) ~ 9, 
    `54_0_0` %in% c(11010, 11014) ~ 10))

data_merged_copy$region %>% freq()
#    Frequency Percent
#1         33838   6.741
#2         68728  13.693
#3         58239  11.603
#4         78772  15.694
#5         35815   7.135
#6         43400   8.647
#7         42954   8.558
#8         20802   4.144
#9         44880   8.941
#10        74508  14.844
#Total    501936 100.000

### Country (derived from the previous variable) -----
data_merged_copy = data_merged_copy %>% 
  dplyr::mutate(country = case_when(
    region == 5 ~ 1, # Scotland
    region == 8 ~ 2, # Wales
    TRUE ~ 3)) # England

data_merged_copy$country %>% freq()
#      Frequency Percent
#1         35815   7.135
#2         20802   4.144
#3        445319  88.720
#Total    501936 100.000

## 20116: Smoking status ----- 
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20116

data_merged_copy = data_merged_copy %>%
  mutate(smoking_status = case_when(
    `20116_0_0` == -3 ~ 9, # No answer
    `20116_0_0` == 0 ~ 0,  # Never
    `20116_0_0` == 1 ~ 1,  # Previous
    `20116_0_0` == 2 ~ 2,  # Current
    TRUE ~ 9))

data_merged_copy$smoking_status %>% freq()
#      Frequency  Percent
#0        273220  54.4332
#1        172839  34.4345
#2         52929  10.5450
#9          2948   0.5873
#Total    501936 100.0000

## 22038: MET minutes per week for moderate activity ----
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=22038
## 22039: MET minutes per week for vigorous activity -----
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=22039
## 22037: MET minutes per week for walking (included in 22040)----
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=22037
## 22040: Summed MET minutes per week for all activity----
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=22040

data_merged_copy <- data_merged_copy %>%
  mutate(
    total_MET_hrwk = `22040_0_0` / 60,
    minsperweek_mod = `22038_0_0`,
    minsperweek_vig = `22039_0_0`
  ) %>%
  mutate(
    #Recode total MET-hours/week into meaningful categories
    pa_group = case_when(
      is.na(total_MET_hrwk) ~ 9,
      total_MET_hrwk < 10 ~ 1,     # Low
      total_MET_hrwk >= 10 & total_MET_hrwk < 50 ~ 2,  # Moderate
      total_MET_hrwk >= 50 ~ 3     # High
    )
  ) %>%
  # Categorize PA based on recommendations
  mutate(
    pa_rec = case_when(
      minsperweek_vig >= 75 ~ 1,  
      minsperweek_mod >= 150 & minsperweek_vig < 75 ~ 1,  
      is.na(minsperweek_vig) & is.na(minsperweek_mod) ~ 9,
      TRUE ~ 0          
    )
  )

data_merged_copy$pa_group %>% freq()
#      Frequency Percent
#1         70902   14.13
#2        194943   38.84
#3        118991   23.71
#9        117100   23.33
#Total    501936  100.00
data_merged_copy$pa_rec %>%  freq()
#      Frequency Percent
#0         73825   14.71
#1        311011   61.96
#9        117100   23.33
#Total    501936  100.00

## 1558: Alcohol intake frequency  -----
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?tk=MozhW45lvgR3ND07YxPpY0fVTVY7dwBZ483974&id=1558

data_merged_copy <- data_merged_copy %>%
  mutate(alc_group = case_when(
    `1558_0_0` == 6 ~ 0,  # None
    `1558_0_0` %in% c(4, 5) ~ 1,  # Occasional
    `1558_0_0` %in% c(2, 3) ~ 2,  # Moderate
    `1558_0_0` == 1 ~ 3,  # High intake
    `1558_0_0` %in% c(-3,NA) ~ 9,  # Prefer not to answer or NA
  ))

data_merged_copy$alc_group %>% freq()
#      Frequency  Percent
#0         40584   8.0855
#1        113732  22.6587
#2        244444  48.7002
#3        101676  20.2568
#9          1500   0.2988
#Total    501936 100.0000



## Dietary exposures from FFQ
# Fruits
# Fruit/veg variables have specific coding: -10 (less than one),
# -3 (prefer not to answer), -1 (do not know)
# 1 portion of cooked (1289)/raw veg (1299) = 3 tablespoons/day
# 1 portion of fresh fruit (1309) = 1 piece/day
# Portion sizes based on: 
# https://www.nhs.uk/live-well/eat-well/5-a-day-portion-sizes/#5-a-day-fruit-portions


## 1289: Cooked vegetable intake ----
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=1289

## 1299: Salad / raw vegetable intake ----
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=1299

## 1309: Fresh fruit intake ----
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=1309

# Step 1 - Recoding variables
data_merged_copy = data_merged_copy %>%
  dplyr::mutate(cooked_veg  = case_when(
    `Cooked vegetable intake | Instance 0` == -10 ~ 0.5,  # recode -10 to 0.5
    `Cooked vegetable intake | Instance 0`  %in% c(-3, -1) ~ NA_real_,  # recode -3 and -1 to NA
    TRUE ~  `Cooked vegetable intake | Instance 0`), 
    raw_veg  = case_when(
      `Salad / raw vegetable intake | Instance 0` == -10 ~ 0.5,  # recode -10 to 0.5
      `Salad / raw vegetable intake | Instance 0`  %in% c(-3, -1) ~ NA_real_,  # recode -3 and -1 to NA
      TRUE ~  `Salad / raw vegetable intake | Instance 0`), 
    fruit  = case_when(
      `Fresh fruit intake | Instance 0` == -10 ~ 0.5,  # recode -10 to 0.5
      `Fresh fruit intake | Instance 0`  %in% c(-3, -1) ~ NA_real_,  # recode -3 and -1 to NA
      TRUE ~  `Fresh fruit intake | Instance 0`)
  )

# Step 2 - Generate fruit-veg portions per day
data_merged_copy = data_merged_copy %>%
  dplyr::mutate(fruit_veg = rowSums(select(., cooked_veg, raw_veg, fruit)), 
                fruit_veg = round(fruit_veg, digits = 0),
                fruit_veg_3 = ntile(fruit_veg, 3), 
                fruit_veg_3 = case_when(
                  is.na(fruit) & is.na(raw_veg) &
                    is.na(cooked_veg) ~ 9,
                  fruit_veg == 0 ~ 0, 
                  is.na(fruit_veg_3) ~ 9, 
                  TRUE ~ fruit_veg_3))


# Step 3 - Indicate whether meeting the fruit-veg recommendations
data_merged_copy = data_merged_copy %>%
  dplyr::mutate(fruit_veg_rec = case_when(
    fruit_veg >= 5 ~ 1,  # meeting recommendation
    fruit_veg < 5 ~ 0,   # not meeting recommendation
    TRUE ~ 9))  # missing values recoded as 9

data_merged_copy$fruit_veg_3 %>% freq()
#      Frequency  Percent
#0          3223   0.6421
#1        160179  31.9122
#2        163401  32.5542
#3        163401  32.5542
#9         11732   2.3373
#Total    501936 100.0000

data_merged_copy$fruit_veg_rec %>% freq()
#Frequency Percent
#0        111703  22.254 # not meeting recommendation
#1        378501  75.408 # meeting recommendation
#9        11732   2.337 # missing values
#Total    501936 100.000

## 2724: Had Menopause  ------
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=2724

# Recode variables for menopausal status
data_merged_copy = data_merged_copy %>%
  dplyr::mutate(
    menopause = case_when(`2724_0_0` == 1 ~ 1,
                          `2724_0_0` == 0 ~ 0,
                          `2724_0_0`%in% c(2, 3, -3) ~ 0,
                          is.na(`2724_0_0`) ~ 9)) 

data_merged_copy$menopause %>% freq()
#      Frequency Percent
#0        107321   21.38
#1        165240   32.92
#9        229375   45.70
#Total    501936  100.00

## 31: Sex -----
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?tk=gcADdik0CdtQfaqQOFUeh3WFjm8biyRO471674&id=31

data_merged_copy = data_merged_copy %>%
  dplyr::mutate(
    sex = case_when(`31_0_0` == 1 ~ 1, #male
                    `31_0_0` == 0 ~ 0, #female
                    is.na(`31_0_0`) ~ 9)) 

data_merged_copy$Sex %>% freq()
#Frequency Percent
#      Frequency Percent
#0        273036    54.4
#1        228900    45.6
#Total    501936   100.0

##1160: Sleep duration-----
#https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=1160
#Based on https://www.sciencedirect.com/science/article/pii/S2352721815000157?via%3Dihub recommendations

data_merged_copy <- data_merged_copy %>%
  mutate(
    sleep = case_when(
      is.na(`1160_0_0`) | `1160_0_0` %in% c(-1, -3) ~ NA,
      `1160_0_0` > 24 ~ NA,
      TRUE ~ round(`1160_0_0`, 1))
  )

data_merged_copy <- data_merged_copy %>%
  mutate(
    sleep_group = case_when(
      is.na(sleep) ~ 0, #Unknown/Invalid,
      sleep < 7 ~ 1, #Insufficient sleep,
      sleep >= 7 & sleep <= 9 ~ 2, #Recommended range",
      sleep > 10 ~ 3 #Excessive sleep
    )
  )

data_merged_copy$sleep_group %>% freq()
#      Frequency  Percent Valid Percent
#0          4210   0.8388        0.8508
#1        123111  24.5272       24.8787
#2        365378  72.7937       73.8370
#3          2145   0.4273        0.4335
#NA's       7092   1.4129              
#Total    501936 100.0000      100.0000    

##3140: Pregnancy-----
#https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=3140

data_merged_copy <- data_merged_copy %>%
  mutate(
    pregnancy = case_when(
      `3140_0_0` %in% c(0, 1, 2) ~ `3140_0_0`,
      TRUE ~ NA_real_
    )
  )

data_merged_copy$pregnancy %>% freq()
#Frequency   Percent Valid
#0        271824  54.15511 #no
#1           150   0.02988 #yes
#2           220   0.04383 #unsure
#NA's     229742  45.77117 #nodata           
#Total    501936 100.00000

##2453: Cancer----
#https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=2453
data_merged_copy <- data_merged_copy %>%
  mutate(
    cancer = case_when(
      `2453_0_0` == 1 ~ 1,
      `2453_0_0` < 0 | is.na("2453_0_0") ~ 9,
      TRUE ~ 0
    ),
  )

data_merged_copy$cancer %>% freq()
#      Frequency  Percent
#0        461514  91.9468 #no
#1         38579   7.6860 #yes
#9          1843   0.3672 #NAs
#Total    501936 100.0000





## 20107: Illnesses of father (CVD,COPD, diabetes and Cancer)----
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20107

## 20110: Illnesses of mother (CVD,COPD, diabetes)----
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20110

## 20111: Illnesses of siblings (CVD,COPD, diabetes)-----
# https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=20111

#Preprocessing
data_merged_copy = data_merged_copy %>% 
  dplyr::mutate(illness_father = `20107_0_0`, illness_mother = `20110_0_0`, 
                illness_siblings = `20111_0_0`) 

data_merged_copy = data_merged_copy %>% 
  dplyr::mutate(across(starts_with("illness_") , 
                       ~ str_split(.x, "\\|"))) %>% 
  unnest_wider(where(is.list), names_sep = "_") 

#Cardiovascular diseases

data_merged_copy = data_merged_copy %>% 
  dplyr::mutate(across(starts_with("illness_") , 
                       ~ as.numeric(.))) %>% 
  dplyr::mutate(across(starts_with("illness_") , 
                       ~ case_when(
                         . %in% c(1) ~ 1, #heartdisease
                         . %in% c(2) ~ 1, #stroke
                         . %in% c(8) ~ 1, #highbloodpressure
                         is.na(.) ~ NA_real_, 
                         TRUE ~ 0), .names = "cvd_{col}")) 

data_merged_copy = data_merged_copy %>% 
  dplyr::mutate(
    cvd_fam = case_when(
      cvd_illness_father_1 == 1 | cvd_illness_father_2 == 1 |  cvd_illness_father_3 == 1 | 
        cvd_illness_father_4 == 1 | cvd_illness_father_5 == 1 | cvd_illness_father_6 == 1 | 
        cvd_illness_father_7 == 1 | cvd_illness_father_8 == 1 | cvd_illness_father_9 == 1 | 
        cvd_illness_father_10 == 1 | 
        cvd_illness_mother_1 == 1 | cvd_illness_mother_2 == 1 |  cvd_illness_mother_3 == 1 | 
        cvd_illness_mother_4 == 1 | cvd_illness_mother_5 == 1 | cvd_illness_mother_6 == 1 | 
        cvd_illness_mother_7 == 1 | cvd_illness_mother_8 == 1 | cvd_illness_mother_9 == 1 | 
        cvd_illness_mother_10 == 1 | cvd_illness_mother_11 == 1 | 
        cvd_illness_siblings_1 == 1 | cvd_illness_siblings_2 == 1 |  cvd_illness_siblings_3 == 1 | 
        cvd_illness_siblings_4 == 1 | cvd_illness_siblings_5 == 1 | cvd_illness_siblings_6 == 1 | 
        cvd_illness_siblings_7 == 1 | cvd_illness_siblings_8 == 1 | cvd_illness_siblings_9 == 1 | 
        cvd_illness_siblings_10 == 1  ~ 1, 
      cvd_illness_father_1 == 0 | cvd_illness_father_2 == 0 |  cvd_illness_father_3 == 0 | 
        cvd_illness_father_4 == 0 | cvd_illness_father_5 == 0 | cvd_illness_father_6 == 0 | 
        cvd_illness_father_7 == 0 | cvd_illness_father_8 == 0 | cvd_illness_father_9 == 0 | 
        cvd_illness_father_10 == 0 | 
        cvd_illness_mother_1 == 0 | cvd_illness_mother_2 == 0 |  cvd_illness_mother_3 == 0 | 
        cvd_illness_mother_4 == 0 | cvd_illness_mother_5 == 0 | cvd_illness_mother_6 == 0 | 
        cvd_illness_mother_7 == 0 | cvd_illness_mother_8 == 0 | cvd_illness_mother_9 == 0 | 
        cvd_illness_mother_10 == 0 | cvd_illness_mother_11 == 0 | 
        cvd_illness_siblings_1 == 0 | cvd_illness_siblings_2 == 0 |  cvd_illness_siblings_3 == 0 | 
        cvd_illness_siblings_4 == 0 | cvd_illness_siblings_5 == 0 | cvd_illness_siblings_6 == 0 | 
        cvd_illness_siblings_7 == 0 | cvd_illness_siblings_8 == 0 | cvd_illness_siblings_9 == 0 | 
        cvd_illness_siblings_10 == 0 | cvd_illness_siblings_11 == 0 | cvd_illness_siblings_12 == 0  ~ 0, 
      TRUE ~ NA_real_)) 

data_merged_copy$cvd_fam %>% freq()  
#      Frequency Percent Valid Percent
#0        125018  24.907         25.34
#1        368354  73.387         74.66
#NA's       8564   1.706              
#Total    501936 100.000        100.00

#Diabetes
data_merged_copy = data_merged_copy %>% 
  dplyr::mutate(across(starts_with("illness_") , 
                       ~ case_when(
                         . == 9 ~ 1, #diabetes
                         is.na(.) ~ NA_real_, 
                         TRUE ~ 0),
                       .names = "diab_{col}")) %>% 
  dplyr::mutate(
    diab_fam = case_when(
      diab_illness_father_1 == 1 | diab_illness_father_2 == 1 |  diab_illness_father_3 == 1 | 
        diab_illness_father_4 == 1 | diab_illness_father_5 == 1 | diab_illness_father_6 == 1 | 
        diab_illness_father_7 == 1 | diab_illness_father_8 == 1 | diab_illness_father_9 == 1 | 
        diab_illness_father_10 == 1 | 
        diab_illness_mother_1 == 1 | diab_illness_mother_2 == 1 |  diab_illness_mother_3 == 1 | 
        diab_illness_mother_4 == 1 | diab_illness_mother_5 == 1 | diab_illness_mother_6 == 1 | 
        diab_illness_mother_7 == 1 | diab_illness_mother_8 == 1 | diab_illness_mother_9 == 1 | 
        diab_illness_mother_10 == 1 | diab_illness_mother_11 == 1 | 
        diab_illness_siblings_1 == 1 | diab_illness_siblings_2 == 1 |  diab_illness_siblings_3 == 1 | 
        diab_illness_siblings_4 == 1 | diab_illness_siblings_5 == 1 | diab_illness_siblings_6 == 1 | 
        diab_illness_siblings_7 == 1 | diab_illness_siblings_8 == 1 | diab_illness_siblings_9 == 1 | 
        diab_illness_siblings_10 == 1  ~ 1, 
      diab_illness_father_1 == 0 | diab_illness_father_2 == 0 |  diab_illness_father_3 == 0 | 
        diab_illness_father_4 == 0 | diab_illness_father_5 == 0 | diab_illness_father_6 == 0 | 
        diab_illness_father_7 == 0 | diab_illness_father_8 == 0 | diab_illness_father_9 == 0 | 
        diab_illness_father_10 == 0 | 
        diab_illness_mother_1 == 0 | diab_illness_mother_2 == 0 |  diab_illness_mother_3 == 0 | 
        diab_illness_mother_4 == 0 | diab_illness_mother_5 == 0 | diab_illness_mother_6 == 0 | 
        diab_illness_mother_7 == 0 | diab_illness_mother_8 == 0 | diab_illness_mother_9 == 0 | 
        diab_illness_mother_10 == 0 | diab_illness_mother_11 == 0 | 
        diab_illness_siblings_1 == 0 | diab_illness_siblings_2 == 0 |  diab_illness_siblings_3 == 0 | 
        diab_illness_siblings_4 == 0 | diab_illness_siblings_5 == 0 | diab_illness_siblings_6 == 0 | 
        diab_illness_siblings_7 == 0 | diab_illness_siblings_8 == 0 | diab_illness_siblings_9 == 0 | 
        diab_illness_siblings_10 == 0 | diab_illness_siblings_11 == 0 | diab_illness_siblings_12 == 0  ~ 0, 
      TRUE ~ NA_real_))


data_merged_copy$diab_fam %>% freq()
#      Frequency Percent
#0        384903  76.684
#1        108469  21.610
#NA's       8564   1.706              
#Total    501936 100.000

##Chronic bronchitis/emphysema
data_merged_copy = data_merged_copy %>% 
  dplyr::mutate(across(starts_with("illness_"), 
                       ~ case_when(
                         . == 6 ~ 1,  # COPD
                         is.na(.) ~ NA_real_,
                         TRUE ~ 0),
                       .names = "copd_{col}")) %>% 
  dplyr::mutate(
    copd_fam = case_when(
      copd_illness_father_1 == 1 | copd_illness_father_2 == 1 | copd_illness_father_3 == 1 | 
        copd_illness_father_4 == 1 | copd_illness_father_5 == 1 | copd_illness_father_6 == 1 | 
        copd_illness_father_7 == 1 | copd_illness_father_8 == 1 | copd_illness_father_9 == 1 | 
        copd_illness_father_10 == 1 | 
        copd_illness_mother_1 == 1 | copd_illness_mother_2 == 1 | copd_illness_mother_3 == 1 | 
        copd_illness_mother_4 == 1 | copd_illness_mother_5 == 1 | copd_illness_mother_6 == 1 | 
        copd_illness_mother_7 == 1 | copd_illness_mother_8 == 1 | copd_illness_mother_9 == 1 | 
        copd_illness_mother_10 == 1 | copd_illness_mother_11 == 1 | 
        copd_illness_siblings_1 == 1 | copd_illness_siblings_2 == 1 | copd_illness_siblings_3 == 1 | 
        copd_illness_siblings_4 == 1 | copd_illness_siblings_5 == 1 | copd_illness_siblings_6 == 1 | 
        copd_illness_siblings_7 == 1 | copd_illness_siblings_8 == 1 | copd_illness_siblings_9 == 1 | 
        copd_illness_siblings_10 == 1 ~ 1,
      
      copd_illness_father_1 == 0 | copd_illness_father_2 == 0 | copd_illness_father_3 == 0 | 
        copd_illness_father_4 == 0 | copd_illness_father_5 == 0 | copd_illness_father_6 == 0 | 
        copd_illness_father_7 == 0 | copd_illness_father_8 == 0 | copd_illness_father_9 == 0 | 
        copd_illness_father_10 == 0 | 
        copd_illness_mother_1 == 0 | copd_illness_mother_2 == 0 | copd_illness_mother_3 == 0 | 
        copd_illness_mother_4 == 0 | copd_illness_mother_5 == 0 | copd_illness_mother_6 == 0 | 
        copd_illness_mother_7 == 0 | copd_illness_mother_8 == 0 | copd_illness_mother_9 == 0 | 
        copd_illness_mother_10 == 0 | copd_illness_mother_11 == 0 | 
        copd_illness_siblings_1 == 0 | copd_illness_siblings_2 == 0 | copd_illness_siblings_3 == 0 | 
        copd_illness_siblings_4 == 0 | copd_illness_siblings_5 == 0 | copd_illness_siblings_6 == 0 | 
        copd_illness_siblings_7 == 0 | copd_illness_siblings_8 == 0 | copd_illness_siblings_9 == 0 | 
        copd_illness_siblings_10 == 0 | copd_illness_siblings_11 == 0 | copd_illness_siblings_12 == 0 ~ 0,
      
      TRUE ~ NA_real_
    )
  )

data_merged_copy$copd_fam %>% freq()
#      Frequency Percent
#0        415458  82.771         
#1         77914  15.523        
#NA's       8564   1.706              
#Total    501936 100.000

##Cancer
data_merged_copy = data_merged_copy %>% 
  dplyr::mutate(across(starts_with("illness_"), 
                       ~ case_when(
                         . %in% c(13, 5, 4, 3) ~ 1,  # Cáncer
                         is.na(.) ~ NA_real_,
                         TRUE ~ 0),
                       .names = "cancer_{col}")) %>% 
  dplyr::mutate(
    cancer_fam = case_when(
      cancer_illness_father_1 == 1 | cancer_illness_father_2 == 1 | cancer_illness_father_3 == 1 | 
        cancer_illness_father_4 == 1 | cancer_illness_father_5 == 1 | cancer_illness_father_6 == 1 | 
        cancer_illness_father_7 == 1 | cancer_illness_father_8 == 1 | cancer_illness_father_9 == 1 | 
        cancer_illness_father_10 == 1 | 
        cancer_illness_mother_1 == 1 | cancer_illness_mother_2 == 1 | cancer_illness_mother_3 == 1 | 
        cancer_illness_mother_4 == 1 | cancer_illness_mother_5 == 1 | cancer_illness_mother_6 == 1 | 
        cancer_illness_mother_7 == 1 | cancer_illness_mother_8 == 1 | cancer_illness_mother_9 == 1 | 
        cancer_illness_mother_10 == 1 | cancer_illness_mother_11 == 1 | 
        cancer_illness_siblings_1 == 1 | cancer_illness_siblings_2 == 1 | cancer_illness_siblings_3 == 1 | 
        cancer_illness_siblings_4 == 1 | cancer_illness_siblings_5 == 1 | cancer_illness_siblings_6 == 1 | 
        cancer_illness_siblings_7 == 1 | cancer_illness_siblings_8 == 1 | cancer_illness_siblings_9 == 1 | 
        cancer_illness_siblings_10 == 1 ~ 1,
      
      cancer_illness_father_1 == 0 | cancer_illness_father_2 == 0 | cancer_illness_father_3 == 0 | 
        cancer_illness_father_4 == 0 | cancer_illness_father_5 == 0 | cancer_illness_father_6 == 0 | 
        cancer_illness_father_7 == 0 | cancer_illness_father_8 == 0 | cancer_illness_father_9 == 0 | 
        cancer_illness_father_10 == 0 | 
        cancer_illness_mother_1 == 0 | cancer_illness_mother_2 == 0 | cancer_illness_mother_3 == 0 | 
        cancer_illness_mother_4 == 0 | cancer_illness_mother_5 == 0 | cancer_illness_mother_6 == 0 | 
        cancer_illness_mother_7 == 0 | cancer_illness_mother_8 == 0 | cancer_illness_mother_9 == 0 | 
        cancer_illness_mother_10 == 0 | cancer_illness_mother_11 == 0 | 
        cancer_illness_siblings_1 == 0 | cancer_illness_siblings_2 == 0 | cancer_illness_siblings_3 == 0 | 
        cancer_illness_siblings_4 == 0 | cancer_illness_siblings_5 == 0 | cancer_illness_siblings_6 == 0 | 
        cancer_illness_siblings_7 == 0 | cancer_illness_siblings_8 == 0 | cancer_illness_siblings_9 == 0 | 
        cancer_illness_siblings_10 == 0 | cancer_illness_siblings_11 == 0 | cancer_illness_siblings_12 == 0 ~ 0,
      
      TRUE ~ NA_real_
    )
  )

data_merged_copy$cancer_fam %>% freq()
#     Frequency Percent
#0        317838  63.322
#1        175534  34.971
#NA's       8564   1.706              



# Covariates to retain in the final dataset
covariates_clean <- c(
  "eid", "age_group", "bmi_group", "ethnicity", "townsend_q5",
  "education_group", "region", "country", "smoking_status",
  "total_MET_hrwk", "minsperweek_mod", "minsperweek_vig", "pa_group", "pa_rec",
  "alc_group", "cooked_veg", "raw_veg", "fruit", "fruit_veg",
  "fruit_veg_3", "fruit_veg_rec", "menopause", "sex", "sleep",
  "sleep_group", "pregnancy", "cancer", "cvd_fam", "diab_fam", "copd_fam", "cancer_fam"
)

# Metabolomics field IDs (prefixes)
metabo_ids <- c(
  23474, 23475, 23476, 23477, 23460, 23479, 23440, 23439, 23441, 23433,
  23432, 23431, 23484, 23526, 23561, 23533, 23498, 23568, 23540, 23505,
  23575, 23547, 23512, 23554, 23491, 23519, 23580, 23610, 23635, 23615,
  23590, 23640, 23620, 23595, 23645, 23625, 23600, 23630, 23585, 23605,
  23485, 23418, 23527, 23417, 23562, 23534, 23499, 23569, 23541, 23506,
  23576, 23548, 23513, 23416, 23555, 23492, 23520, 23581, 23611, 23636,
  23616, 23591, 23641, 23621, 23596, 23646, 23626, 23601, 23631, 23586,
  23606, 23473, 23404, 23481, 23430, 23523, 23429, 23558, 23530, 23495,
  23565, 23537, 23502, 23572, 23544, 23509, 23428, 23551, 23488, 23516,
  23478, 23443, 23450, 23457, 23486, 23422, 23528, 23421, 23563, 23535,
  23500, 23570, 23542, 23507, 23577, 23549, 23514, 23420, 23556, 23493,
  23521, 23582, 23612, 23637, 23617, 23592, 23642, 23622, 23597, 23647,
  23627, 23602, 23632, 23587, 23607, 23470, 23461, 23462, 23480, 23406,
  23463, 23465, 23405, 23471, 23466, 23449, 23456, 23447, 23454, 23444,
  23451, 23445, 23459, 23452, 23468, 23437, 23434, 23483, 23414, 23525,
  23413, 23560, 23532, 23497, 23567, 23539, 23504, 23574, 23546, 23511,
  23412, 23553, 23490, 23518, 23579, 23609, 23634, 23614, 23589, 23639,
  23619, 23594, 23644, 23624, 23599, 23629, 23584, 23604, 23446, 23458,
  23453, 23472, 23402, 23448, 23455, 23438, 23400, 23401, 23436, 23464,
  23427, 23415, 23442, 23419, 23482, 23426, 23524, 23425, 23559, 23531,
  23496, 23423, 23566, 23538, 23503, 23573, 23545, 23510, 23424, 23552,
  23489, 23517, 23411, 23407, 23487, 23410, 23529, 23409, 23564, 23536,
  23501, 23571, 23543, 23508, 23578, 23550, 23515, 23408, 23557, 23494,
  23522, 23435, 23583, 23613, 23638, 23618, 23593, 23643, 23623, 23598,
  23648, 23628, 23603, 23633, 23588, 23608, 23469, 23403, 23467, 20281,20280
)
for (id in metabo_ids) {
  pattern <- paste0("^", id, "_")
  cols <- grep(pattern, names(data_merged_copy), value = TRUE)
  for (col in cols) {
    new_name <- paste0("metabo_", col)
    data_merged_copy[[new_name]] <- data_merged_copy[[col]]
  }
}

# Add metabolomics variables to the final variable list
metabo_vars <- grep("^metabo_", names(data_merged_copy), value = TRUE)

# NCD field IDs (non-communicable diseases)
ncd_ids <- c(
  130792, 130793, 131286, 131287, 131294, 131295, 131306, 131307,
  131354, 131355, 22150, 130708, 130709, 130710, 130711, 130712,
  130713, 130714, 130715
)

for (id in ncd_ids) {
  pattern <- paste0("^", id, "_")
  cols <- grep(pattern, names(data_merged_copy), value = TRUE)
  for (col in cols) {
    new_name <- paste0("ncd_", col)
    data_merged_copy[[new_name]] <- data_merged_copy[[col]]
  }
}

ncd_vars <- grep("^ncd_", names(data_merged_copy), value = TRUE)

vo2max_vars <- c("30038_0_0", "30038_1_0")

# Final list of variables to include
final_vars <- c(covariates_clean, metabo_vars, ncd_vars, vo2max_vars)

# Subset the merged dataset to include only selected variables
data_merged_copy <- data_merged_copy %>% select(all_of(final_vars))

# Save the final dataset to file, next time prevalent discard ncd diseases for elastic-net construction
saveRDS(data_merged_copy, "./upload/data_merged_copy_17_02_2026_needed_to_merge_by_ids.rds")


data_merged_copy <- readRDS("./upload/data_merged_copy_17_02_2026_needed_to_merge_by_ids.rds")



# Identify renamed metabolomics columns from instance 0_0_0
metabo_cols <- grep("^metabo_.*_0_0$", names(data_merged_copy), value = TRUE)

# Keep subjects with at least one non-missing metabolomics value (instance 0). # Save the final dataset to file, next time prevalent discard ncd diseases for elastic-net construction
data_merged_copy <- data_merged_copy[apply(data_merged_copy[, metabo_cols, drop = FALSE], 1, function(x) any(!is.na(x))), ]



saveRDS(data_merged_copy, "./upload/data_merged_copy_17_02_2026_at_least_1_missing_metabolite_needed_to_merge_by_ids.rds")

# Upload the files to RAP UK Bibank----
## Locate all the files generated in the upload folder 

# setwd("./upload")

# system("dx ls data/processed")

# system("pwd")

system("dx upload /home/rstudio-server/upload/*")