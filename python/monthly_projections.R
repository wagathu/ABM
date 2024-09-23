

require(tidyr)

# working on one dataset
bdf2 <- readRDS("data/year_pop.rds")


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

# plotting
mbdf4 |> mutate(year = factor(year)) |>
  # na.omit() |>
  filter(county == "Baringo" & gender == "Male") |>
  ggplot(aes(x = date, y = population)) + 
  geom_point() +
  facet_wrap(~Age, scales = "free")

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

saveRDS(mbdf5, "monthly_projections.rds")
