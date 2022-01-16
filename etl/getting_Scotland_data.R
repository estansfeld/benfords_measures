# some background
# http://statistics.gov.scot/data/deaths-involving-coronavirus-covid-19
# https://medium.swirrl.com/using-r-to-analyse-linked-data-7225eefe2eb8
#
# training resources
# https://guides.statistics.gov.scot/article/22-querying-data-with-sparql
# https://guides.statistics.gov.scot/article/34-understanding-the-data-structure


library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(here)

source(here("etl", "SG_SparQL_Queries.R"))

pd <- get_SG_Population_SparQL()
pop_full <- pd$results
pop_boards <- pop_full %>% 
  dplyr::filter(sex == "All", age == "All") 
  

# get the data
qd <- get_SG_CovidDeaths_SparQL()
covid_full <- data.table(qd$results)
summary(covid_full)
# 
# cd <- get_SG_CovidDeaths_Region_SparQL()
# covid_region <- data.table(cd$results)

# add a date field to use instead of the week commencing field
# get the area code from the URI
# convert Cause into a factor
covid_region <- covid_full[
  sex == "All" & cause == "COVID-19 related" & location == "All" & age == "All"
][
  periodname != 2020 & periodname != 2021
][
  j= c("week_beginning", "areauri"):=.(
    ymd(str_extract(string = periodname, pattern = "\\d{4}[-]\\d{2}[-]\\d{2}")),
    str_extract(areauri, "\\S\\d\\d*")
  )
][
  j = .(deaths = sum(nDeaths)),
  by = .(week_beginning, areauri, region = area)
]

write.csv(covid_region, "covid_20_scotland.csv")

summary(Covid_full)
Up_to_week_beginning <- max(Covid_full$week_beginning, na.rm = T)

# add the Scotland Population Projection for 2020
Covid_full <- Covid_full %>% 
  left_join(pop_boards %>% select(area, nPeople)) %>% 
  mutate(deaths_per_thousand = nDeaths / (nPeople/1000)) 

#split the data into three tables:
#  - Scotland_Location_YearToDate
#  - Scotland_Location_by_Week
#  - Scotland_All_by_Week
#  - BoardsCause_YearToDate
#  - BoardsCause_by_Week
#  - Demographics_YearToDate
#  - Demographics_by_Week

Scotland_Location_YearToDate <- Covid_full %>% 
  dplyr::filter(periodname == 2020 & area == "Scotland" & sex == "All" & age == "All") %>% 
  select(-week_beginning, -sex, -age)

Scotland_Location_by_Week <- Covid_full %>%
  dplyr::filter(area == "Scotland" & periodname != 2020 & sex == "All" & age == "All") %>% 
  select( -age, -sex)  

Scotland_5Year_Mean_by_Week <- Covid_full %>%
  dplyr::filter(periodname != 2020 & cause == "All causes - average of corresponding week over previous 5 years") %>% 
  select(-area, -areauri, -age, -sex, -location, -deaths_per_thousand, -nPeople) %>% 
  mutate(cum_deaths = cumsum(nDeaths))

Scotland_nDeaths_Total <- Covid_full %>%
  dplyr::filter(periodname == 2020 & sex != "All") %>% 
  select(-area, -areauri, -location, -week_beginning, -deaths_per_thousand, -nPeople) %>% 
  group_by(cause) %>% 
  summarise(total_deaths = sum(nDeaths))

Scotland_All_by_Week <- Covid_full %>%
  dplyr::filter(area == "Scotland" & 
                  periodname != 2020 & 
                  sex == "All" & 
                  age == "All" &
                  location ==  "All") %>% 
  select(-area, -areauri, -age, -sex, -location, -deaths_per_thousand) %>% 
  pivot_wider(names_from = cause, values_from = c(nDeaths)) %>% 
  rename(deaths_COVID_19 = `COVID-19 related`) %>% 
  rename(All = `All causes`) %>% 
  rename(deaths_5_Year_Av = `All causes - average of corresponding week over previous 5 years`) %>% 
  mutate(deaths_weekly_COVID_19 = replace_na(deaths_COVID_19, 0)) %>% 
  mutate(deaths_weekly_Not_Covid = All - deaths_weekly_COVID_19) %>% 
  mutate(deaths_cum_COVID = cumsum(deaths_weekly_COVID_19)) %>% 
  mutate(deaths_cum_Not_COVID = cumsum(deaths_weekly_Not_Covid)) %>% 
  mutate(deaths_YTD_2020 = cumsum(All)) %>% 
  mutate(deaths_YTD_5_Year_Av = cumsum(deaths_5_Year_Av)) %>%
  mutate(Excess_Weekly_Deaths = deaths_YTD_2020 - deaths_YTD_5_Year_Av) %>% 
  pivot_longer(col = starts_with("deaths"), names_to = "cause")

Scotland_Overall <- Covid_full %>%
  dplyr::filter(area == "Scotland" & 
                  periodname == 2020 & 
                  sex == "All" & 
                  age == "All" &
                  location ==  "All") %>% 
  select(-area, -areauri, -age, -sex, -location, -deaths_per_thousand) %>%
  pivot_wider(names_from = cause, values_from = c(nDeaths)) %>%
  rename(`YTD_Deaths_5_Year_Avg` = `All causes - average of corresponding week over previous 5 years`) %>%
  select(- `COVID-19 related`) %>% 
  rename(`YTD_Deaths 2020` = `All causes`) %>% 
  mutate(Excess_Deaths = `YTD_Deaths 2020` - `YTD_Deaths_5_Year_Avg`) 

Scotland_Age_YearToDate <- Covid_full %>% 
  dplyr::filter(periodname == 2020 & area == "Scotland" & sex == "All" & location == "All") %>% 
  select(-week_beginning, -sex, -location) 

Scotland_Age_by_Week <- Covid_full %>%
  dplyr::filter(area == "Scotland" & periodname != 2020 & sex == "All" & location == "All") %>% 
  select( -location, -sex) 

BoardsCause_YearToDate<- Covid_full %>% 
  dplyr::filter(periodname == 2020 & areauri == "S08") %>% 
  select(-week_beginning, -sex, -age)

BoardsCause_by_Week<- Covid_full %>% 
  dplyr::filter(periodname != 2020 & areauri == "S08") %>% 
  select(-age, -sex, -areauri, -location)

CouncilsCause_YearToDate<- Covid_full %>% 
  dplyr::filter(periodname == 2020 & areauri == "S12") %>% 
  select(-week_beginning, -sex, -age)

CouncilsCause_by_Week<- Covid_full %>% 
  dplyr::filter(periodname != 2020 & areauri == "S12") %>% 
  select(-age, -sex, -location)

Demographics_YearToDate<- Covid_full %>% 
  dplyr::filter(periodname == 2020 & sex != "All") %>% 
  select(-week_beginning, -area, -areauri, -location) %>% 
  left_join(Scotland_nDeaths_Total) %>% 
  mutate(share_of_deaths = nDeaths / total_deaths)

Demographics_by_Week<- Covid_full %>% 
  dplyr::filter(periodname != 2020 & sex != "All" & age !="All") %>% 
  select(-area, -areauri, -location)

External_events <- tribble(
  ~ date, ~ event,
  "2020-01-30", "Global health\nemergency declared",
  "2020-03-11", "Pandemic\ndeclared",
  "2020-02-13", "China reporting\nchange",
  "2020-03-23", "UK lockdown\nstarted",
  "2020-06-18", "Scotland phase 2 started"
) %>%
  mutate(date = as.Date(date))

write.csv(Covid_full, 
          file = here("results", str_c("NRS_Scotland_All_to_wb_", Up_to_week_beginning,".csv")))
  
  # export data
write.csv(Scotland_All_by_Week, 
          file = here("results", "Scotland_All_by_Week.csv"))

write.csv(BoardsCause_by_Week, 
          file = here("results", "BoardsCause_by_Week.csv"))

write.csv(CouncilsCause_by_Week, 
          file = here("results", "CouncilsCause_by_Week.csv"))

write.csv(Demographics_by_Week, 
          file = here("results", "Demographics_by_Week.csv"))

write.csv(Scotland_All_by_Week, 
          file = here("results", str_c("Scotland_All_by_Week_to_wb_", Up_to_week_beginning,".csv")))

write.csv(BoardsCause_by_Week, 
          file = here("results", str_c("BoardsCause_by_Week_to_wb_", Up_to_week_beginning,".csv")))

write.csv(CouncilsCause_by_Week, 
          file = here("results", str_c("CouncilsCause_by_Week_to_wb_", Up_to_week_beginning,".csv")))

write.csv(Demographics_by_Week, 
          file = here("results", str_c("Demographics_by_Week_to_wb_", Up_to_week_beginning,".csv")))

