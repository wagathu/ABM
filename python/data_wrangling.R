
# Importing packages ------------------------------------------------------

pacman::p_load(
  dplyr,
  ggplot2,
  data.table,
  lubridate,
  ggalt,
  patchwork,
  stringr,
  rKenyaCensus,
  INLA,
  inlabru,
  tidyr,
  readxl
)


# Importing data ----------------------------------------------------------

k <- fread("data/kenya_measles.csv")
co <- fread("data/county_measles.csv")
bdf2 <- readRDS("data/year_pop.rds")
who_cases <- read_excel("data/who_cases.xlsx")
who_month_cases <- read_excel("data/who_monthly_cases.xlsx")
who_vaccs <- read_excel("data/who_vaccs.xlsx")
who_inci <- read_excel("data/who_inci.xlsx")
week_measles <- fread("data/weekly_measles.csv")

# Data wrangling ----------------------------------------------------------

# Cleaning WHO cases
wc <- data.frame(t(who_cases))[-c(1:2),1] %>%
  data.frame("cases" = ., year = rev(1980:2023)) |> 
  mutate(cases = as.numeric(str_remove_all(cases, "\\,")))
  
wi <- data.frame(t(who_inci))[-c(1:3),1] %>%
  data.frame("inci" = ., year = rev(1980:2023)) |> 
  mutate(inci = as.numeric(str_remove_all(inci, "\\,")))

wmc <- who_month_cases %>% 
  setNames(str_to_lower(colnames(.))) |> 
  dplyr::select(-region, -iso3, -country) |> 
  data.table() |> 
  melt(id.vars = "year") |> 
  mutate(date = ymd(paste(year, str_to_title(variable), "01", sep = "-"))) |> 
  dplyr::select(date, cases = value) |> 
  mutate(cases = as.numeric(cases))

wk_measles <- week_measles |>
  _[, date := periodname] |>
  _[, `:=`(
    county = orgunitlevel2,
    subcounty = orgunitlevel3,
    cases = `IDSR Measles Total`
  )] |> 
  _[county == "", ] |>
  _[, end_of_week := as.Date(sub(".*-\\s*(\\d{4}-\\d{2}-\\d{2})", "\\1", date))] |> 
  _[, .(date = end_of_week, cases)]
fwrite(wk_measles, "data/wk_measles.csv", row.names = F)

k2 <- k |> 
  _[, .(
    month = str_split(periodname, " ", simplify = T)[,1],
    year = str_split(periodname, " ", simplify = T)[,2],
    deaths = `IDSR Measles Deaths`,
    cases = `IDSR Measles Total`,
    mcv1 = `Proportion of under 1 year receiving vaccine against Measles and Rubella 1`,
    mcv2 = `Proportion of under two years receiving  vaccine against Measles and Rubella 2`
  )] |> 
  _[, date := ymd(paste(year, month, "01", sep = "-"))] |> 
  _[, !c("month", "year")] |> 
  setcolorder(c("date", "cases", "deaths", "mcv1", "mcv2"))
fwrite(k2, "data/clean_kenya_measles.csv", row.names = F)

# Population dataset ------------------------------------------------------

# Age and population
pop_age <- rKenyaCensus::V3_T2.2 %>%
  setNames(c(str_to_lower(colnames(.)))) |>
  data.table() |> 
  _[!age %in% c("Total", "NotStated") & !str_detect(age, "\\-"),] |> 
  _[, .(age = str_replace_all(age, "\\+", "") |> as.numeric(), value = total)]
fwrite(pop_age, "data/pop_age.csv", row.names = F)


# Population projections using INLA ---------------------------------------

mbdf2 <- merge(
  bdf2 |> mutate(month = "Jun"),
  expand.grid(
    county = bdf2 |> pull(county) |> unique(),
    Age = bdf2 |> pull(Age) |> unique(),
    gender = bdf2 |> pull(gender) |> unique(),
    year = c(bdf2 |> pull(year) |> unique(), 2026:2030),
    month = month.abb
  ),
  by = c('county', 'Age', 'gender', 'year', 'month'),
  all = T
)

mbdf3 <- mbdf2 |>
  mutate(date = paste0(year, '-', month) |> ym()) |>
  group_by(county, Age, gender) |>
  arrange(date) |>
  mutate(time = row_number()) |>
  ungroup()

mbdf4 <- mbdf3 |>
  filter(year > 2019) |>
  mutate(Age = as.factor(Age),
         gender = as.factor(gender),
         county = as.factor(county),
         countynum = as.integer(county),
         gtime = time,
         gmat = paste0(county, "-", Age, "-", gender) |> factor(),
         gmatnum = as.integer(gmat),
         logpop = log(population))

i2 <- inla(
  logpop ~ 
    f(gtime, model = "linear") + # global
    f(gmat, model = "iid") + # unstructured effects
    f(gmatnum, time, model = "iid"), # group specific slopes
  family = "gaussian",
  data = mbdf4,
  control.predictor = list(link = 1)
  # verbose = T
)
beepr::beep(2)

mbdf5 <- mbdf4 |>
  mutate(projections = i2$summary.fitted.values$mean |> exp()) |>
  select(county, Age, gender, year, month, date, population, projections)

mbdf5 |>
  filter(county == "Kajiado" & gender == "Male") |>
  ggplot(aes(x = date, y = population)) + 
  geom_point() +
  geom_line(aes(x = date, y = projections)) +
  facet_wrap(~Age, scales = "free")

ky <- mbdf5 %>%
  filter(month(date) == 12 & gender == "Total") |> 
  summarise(n_alive = sum(projections), .by = c(year)) |> 
  mutate(year = as.numeric(year)) |> 
  arrange(year)
fwrite(ky, "data/ky.csv")

cy <- mbdf5 %>%
  filter(month(date) == 12 & gender == "Total") |> 
  summarise(n_alive = sum(projections), .by = c(year, county)) |> 
  mutate(year = as.numeric(year)) |> 
  arrange(year)
fwrite(cy, "data/cy.csv")

# Plotting ----------------------------------------------------------------


wc |> 
  filter(year >= 2000) |> 
  ggplot(aes(x = year)) +
  geom_xspline(aes(y = cases)) +
  scale_y_continuous(labels = scales::label_comma()) +
  theme_light()
wi |> 
  filter(year >= 2000) |> 
  ggplot(aes(x = year)) +
  geom_xspline(aes(y = inci)) +
  scale_y_continuous(labels = scales::label_comma()) +
  theme_light()

wmc |> 
  filter(date >= as.Date("2020-01-01")) |> 
  ggplot(aes(x = date)) +
  geom_line(aes(y = cases)) +
  scale_y_continuous(labels = scales::label_comma()) +
  theme_light()

wk_measles |> 
  ggplot(aes(x = date)) +
  geom_xspline(aes(y = cases))+
  theme_light()

wi
k2 |> 
  ggplot(aes(x = date)) +
  geom_xspline(aes(y = cases)) +
  theme_minimal()

pop_age |> 
  ggplot(aes(x = age)) +
  geom_col(aes(y = total)) +
  theme_minimal() +
  scale_y_continuous(labels = scales::label_comma())
