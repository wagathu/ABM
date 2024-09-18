
# Importing packages ---------------------------------------------------------------------

if(require(pacman))
  p_load(data.table, dplyr, readxl, reticulate, ggplot2, tidyverse)

# Importing the needed datasets ----------------------------------------------------------

cmy <- fread('data/measles_df.csv')
brt <- fread('data/demo.csv') 

# Creating datasets ----------------------------------------------------------------------

brt <- brt |> 
  _[, county := str_to_title(county)] |> 
  _[, !('V1')] |> 
  _[, county := case_when(
    county == "Murang'a" ~ 'Muranga',
    county == 'Tharaka-Nithi' ~ 'Tharaka Nithi',
    TRUE ~ county
    
  )]

scenarios <- data.table(
  mcv1 = rep(.95, 10),
  mcv2 = c(.95, .80, .70, .60, .50, .40, .30, .20, .10, 1e-10)
)

# The parameter dataframes
mm = cmy |> 
  setDT() |> 
  _[, .( county, measles_incidence, mcv1 = prop_mcv1, mcv2 = prop_mcv2)] |> 
  _[, map_if(.SD, is.numeric, mean), by = 'county'] 

 pars_df <- brt |> 
   _[, merge(.SD,
             mm,
             by = 'county'
             )] |> 
   _[, `:=`(
     initial_prev = measles_incidence/1e6,
     initial_immunity = ((mcv1/1e2)*under_fives*.85)/under_fives,
     death_rate = cdr
     
   )]
 
 write.csv(pars_df, 'pars_df.csv', row.names = F)
 