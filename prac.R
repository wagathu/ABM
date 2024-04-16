require(reticulate);require(data.table);require(ggplot2);require(ggalt) 
ss <- import('starsim')
plt <- import("matplotlib.pyplot")

# Set the matplotlib backend to inline

pars <- list(
  n_agents = 5000,
  birth_rate = 20,
  death_rate = 15,
  
  networks = list(
    type = 'randomnet',
    n_contacts = 4
  ),
  diseases = list(
    type = 'sir',
    dur_inf = 10,
    beta = 1.8
  )
)

mcv1 = ss$sir_vaccine(pars = list(efficacy = 0.7))
mcv2 = ss$sir_vaccine(pars = list(efficacy = 0.5))

myIntervention1 = ss$routine_vx(
    start_year = 2000,
    prob = 0.95,
    product = mcv1
)

myIntervention2 = ss$routine_vx(
  start_year = 2015,
  prob = 0.1,
  product = mcv2
)

myIntervention3 = c(myIntervention1, myIntervention2)
sim_base <- ss$Sim(pars = pars)
base = sim_base$run()
sim_intervention <- ss$Sim(pars = pars, interventions = myIntervention3)
interv = sim_intervention$run()

df <- data.table(
  time = base$yearvec,
  base_prevalence = base$results$sir$prevalence,
  inter_prevalence = sim_intervention$results$sir$prevalence
) |> 
  _[, melt(.SD,
           id.vars = 'time'
           )]

df |> 
  ggplot(aes(x = time)) +
  geom_xspline(aes(y = value, col = variable)) +
  geom_vline(aes(xintercept = 2000), linetype = 2) +
  theme_light() +
  theme(legend.position = 'bottom') 

# immunization status  = Efficacy * vaccination coverage
