
# Importing data -------------------------------------------------------------------------


if (require(pacman))
  p_load(data.table,
         dplyr,
         ggplot2,
         scales,
         dplyr,
         ggalt,
         lubridate,
         purrr,
         patchwork,
         plotly
         )

# Imporitng data -------------------------------------------------------------------------

df.mcv1 <- fread('data/mcv1_only_scenario.csv') |> 
  data.table()
df.mcv1_mcv2 <- fread('data/mcv1_at_95_and_mcv2_wvar.csv')


# Cleaning data --------------------------------------------------------------------------

cols.numeric <- c("prevalence", "new_infections", "susceptible", 
                 "infected", "recovered", "exposed", "cum_infections")

cols <- colnames(df.mcv1)

df1 <- copy(df.mcv1) |> 
  setDT() |> 
  _[, map_at(.SD, 4:10, as.numeric) ] |> 
  _[, time := floor(time)] |> 
  _[as.numeric(time) <= 2049, ] |> 
  _[, month := rep(c(1:12), 300)] |> 
  _[, date := ymd(paste(time, month, '01', sep = '-'))] |> 
  _[, mcv1 := as.factor(mcv1)] |> 
  _[, !c('time', 'month', 'V1')] |> 
  _[, melt(
    .SD,
    id.vars = c('date', 'mcv1'),
    variable.name = 'results'
    )]

# Plotting -------------------------------------------------------------------------------


plts <- df1 %>%
  split(.$results) |> 
  map(
    ~. |> 
      ggplot(aes(x = date)) +
      geom_xspline(aes(y = value, color = mcv1)) +
      theme_light() +
      ggtitle(.$results)
  ) |> 
  setNames(c("prevalence", "new_infections", "susceptible", 
             "infected", "recovered", "exposed", "cum_infections"))

map(plts, ggplotly)
