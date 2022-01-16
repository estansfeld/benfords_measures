library(data.table)
library(stringr)
library(here)

india_case_data <- fread(here("etl", "covid_19_india.csv"))
setnames(india_case_data, c("State/UnionTerritory", "Confirmed", "Deaths", "Cured"), c("state", "cum_confirmed", "cum_deaths", "cum_cured"))

# fix typos for Telangana state
india_case_data <- india_case_data[
  j = state:=fifelse(str_starts(state, "Tel"), "Telangana", state)
][
  j = state:=fifelse(str_starts(state, "Karn"), "Karnataka", state)
]
setkey(india_case_data, state, Date)

# summarise states and remove states with fewer than 100 records
india_states <- india_case_data[
  j = .(N=.N, T = sum(as.numeric(cum_confirmed), na.rm=T)), by = .(state)
][
  i = N>100, .(state)]


india_daily_cases <- india_case_data[india_states, on = .(state)][,
                                                                  j = .(Date, state, cum_confirmed, cum_deaths, cum_cured)
][
  j = c("confirmed", "deaths", "cured") := lapply(.SD, function(x) c(NA, diff(x))), 
  .SDcols = c("cum_confirmed", "cum_deaths", "cum_cured"), by = .(state)
][
  i = !is.na(confirmed)
][
  j = c("cum_confirmed", "cum_deaths", "cum_cured"):=NULL
][j = c("magnitude_confirmed", "magnitude_deaths", "magnitude_cured"):=.(
  fifelse(confirmed > 0, floor(log10(confirmed)), NA_integer_), 
  fifelse(deaths >0, floor(log10(deaths)),NA_integer_),
  fifelse(cured >0, floor(log10(cured)),NA_integer_))
]
setkey(india_daily_cases, state, Date)

whole_scotland_MIS <- fread(here("etl", "SG_covid_whole_data_set.csv"))
scotland_confirmed_cases <- whole_scotland_MIS[Variable == "Testing - New cases reported", .(date = DateCode, areauri = GeographyCode, confirmed = as.numeric(Value))][j = c("magnitude_confirmed"):=.(
  fifelse(confirmed >0, floor(log10(confirmed)),NA_integer_))
]# [!is.na(magnitude_confirmed)]

scotland_case_data <- fread(here("etl", "covid_20_scotland.csv"))
setnames(scotland_case_data, "week_beginning", "Date")
setkey(scotland_case_data, region, Date)
scotland_weekly_deaths <- scotland_case_data[region !="Scotland"][j = c("magnitude_deaths"):=.(
  fifelse(deaths >0, floor(log10(deaths)),NA_integer_))
]# [!is.na(magnitude_deaths)]
# summarise deaths by state
scotland_state <- scotland_weekly_deaths[i = str_starts(areauri, "S08"),
                                         j = .(N=.N, T = sum(as.numeric(deaths), na.rm=T)), by = region
][j=.(region)]
# summarise deaths by state
scotland_states <- scotland_weekly_deaths[
  j = .(N=.N, T = sum(as.numeric(deaths), na.rm=T)), by = region
][j = .(region)]
setkey(scotland_states, "region")
# daily data
MIS_pos_test <- whole_scotland_MIS[
  Variable == "Testing - Cumulative people tested for COVID-19 - Positive", .(date = DateCode, areauri = GeographyCode, cum_pos = as.numeric(Value))
][
  j = c("positive") := lapply(.SD, function(x) c(NA, diff(x))), 
  .SDcols = c("cum_pos"), by = .(areauri)
][
  j = c("magnitude_positive"):=.(
    fifelse(positive !=0, floor(log10(abs(positive))),NA_integer_))
][
  !is.na(positive) & positive !=0
][
  , cum_pos:=NULL
]
# get the uri names
MIS_pos_test <- scotland_case_data[
  str_starts(areauri, "S08")
][
  j=.(.N), by = .(region, areauri)
][
  MIS_pos_test, on = .(areauri)
][
  j = c("N", "areauri"):=NULL
]
# summarise MIS_pos_test by variable
postives <- MIS_pos_test[i = positive != 0,
                         j = .(N=.N), by = .(region)
]
# summarise MIS by variable
variable <- whole_scotland_MIS[
  j = .(N=.N), by = .(GeographyCode, Variable)
]
scotland_deaths <- scotland_weekly_deaths[deaths != 0]$deaths
scotland_deaths_10 <- scotland_deaths[scotland_deaths>=10]
scotland_confirmed <- scotland_confirmed_cases[confirmed!=0]$confirmed
scotland_confirmed_10 <- scotland_confirmed[scotland_confirmed>=10]