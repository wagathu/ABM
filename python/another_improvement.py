# Modules
import numpy as np
import sciris as sc
import matplotlib.pyplot as pl
import starsim as ss
from plotnine import *
import pandas as pd
import matplotlib.pyplot as plt


# Data
kenya_popsize = pd.read_csv("data/ky.csv")
pop_age = pd.read_csv("data/pop_age.csv")
pop_age['value'] = pop_age['value'].astype(float)

# Measles class -------------------------------------------------------------------------------------------

class SEIR(ss.Infection):
    """
    Example SEIR model
    This class implements a basic SEIR model with states for susceptible,
    exposed, infected/infectious, and recovered. It also includes deaths and basic
    results.
    """
    
    def __init__(self, pars=None, **kwargs):
        super().__init__()
        self.default_pars(
            beta = 18/9 ,  # Mean transmission rate: R0 / D
            init_prev = ss.bernoulli(p=0.005),
            init_immune = ss.bernoulli(p=0.01)
            dur_exp = ss.lognorm_ex(mean=10/12, stdev=2),
            dur_inf = ss.lognorm_ex(mean=9/12, stdev=2),
            p_death = ss.bernoulli(p=0.018)

        )
        self.update_pars(pars, **kwargs)

        self.add_states(
            ss.BoolArr('exposed', label='Exposed'),
            ss.BoolArr('recovered', label='Recovered'),
            ss.FloatArr('ti_exposed', label='Time of exposure'),
            ss.FloatArr('ti_infectious', label='Time of becoming infectious'),
            ss.FloatArr('ti_recovered', label='Time of recovery'),
            ss.FloatArr('ti_dead', label='Time of death'),
        )
        return
    
    @property
    def infectious(self):
        return self.infected
      
     def set_initial_states(self, sim):
        super().set_initial_states(sim)
        imm_prob = self.pars.imm_prob
        alive_uids = ss.true(sim.people.alive & self.susceptible)
        initially_immune = imm_prob > np.random.rand()

        self.susceptible[initially_immune] = False
        self.recovered[initially_immune] = True
        

    def update_pre(self):
        """ Update states before the next time step """
        sim = self.sim
        p = self.pars
        ti = sim.ti
        dt = sim.dt
    
        # # # handle beta here: at start of infection prior to transmission
        # beta_mean = 0.15  # Mean transmission rate
        # beta_amplitude = 0.21   #Amplitude of seasonal forcing
        # beta_rate = beta_mean * (1 + beta_amplitude * np.cos(4 * np.pi * ti/12))
        # beta_prob = 1 - np.exp(-beta_rate)
        # 
        # # # Dynamically get the network keys from the simulation
        # network_keys = self.sim.networks.keys()
        # self.pars.beta = {key: [beta_prob, beta_prob] for key in network_keys}

        # conditions for all older people above 20 years to never get measles
        all_ids_above_20 = sim.people.uid[sim.people.age > 20]
        self.susceptible[all_ids_above_20] = False
        self.recovered[all_ids_above_20] = True

        # Progress exposed -> infectious
        new_infectious = (self.exposed & (self.ti_infectious <= ti) ).uids
        self.exposed[new_infectious] = False
        self.infected[new_infectious] = True

        # Progress infectious -> recovered
        recovered = (self.infected & (self.ti_recovered <= ti)).uids
        self.infected[recovered] = False
        self.recovered[recovered] = True

        # Trigger deaths
        deaths = (self.ti_dead <= ti).uids
        if len(deaths):
            sim.people.request_death(deaths)

        return

    def set_prognoses(self, uids, source_uids=None):
        """ Set prognoses """
        super().set_prognoses(uids, source_uids)
        ti = self.sim.ti
        dt = self.sim.dt
        self.susceptible[uids] = False
        self.exposed[uids] = True
        self.ti_exposed[uids] = ti
        p = self.pars
        
        # Sample durations, being careful to only sample from the
        # distributions once per timestep.
        dur_exp = p.dur_exp.rvs(uids)
        dur_inf = p.dur_inf.rvs(uids)

        # Set time of becoming infectious
        self.ti_infectious[uids] = ti + dur_exp / dt

        # Determine who dies and who recovers and when
        will_die = p.p_death.rvs(uids)
        dead_uids = uids[will_die]
        rec_uids = uids[~will_die]
        self.ti_dead[dead_uids] = ti + (dur_exp[will_die] + dur_inf[will_die]) / dt
        self.ti_recovered[rec_uids] = ti + (dur_exp[~will_die] + dur_inf[~will_die]) / dt
        return

    def update_death(self, uids):
        """ Reset exposed/infected/recovered flags for dead agents """
        self.susceptible[uids] = False
        self.exposed[uids] = False
        self.infected[uids] = False
        self.recovered[uids] = False
        return
      
    def update_results(self):
        super().update_results()
        res = self.results
        ti = self.sim.ti
        
        # Count the number of new exposures during this timestep
        new_exposures = np.count_nonzero((self.ti_exposed == ti) & self.exposed)
        
        # Update results accordingly
        res.new_infections[ti] = new_exposures
        res.cum_infections[ti] = np.sum(res['new_infections'][:ti+1])
        res.prevalence[ti] = (res.n_infected[ti] + res.n_exposed[ti]) / np.count_nonzero(self.sim.people.alive)
        return

    def plot(self, plot_kw=None):
        """ Default plot for SEIR model """
        fig = pl.figure()
        plot_kw = sc.mergedicts(dict(lw=2, alpha=0.8), plot_kw)
        for rkey in ['n_susceptible', 'n_exposed', 'n_infected', 'n_recovered']:
            pl.plot(self.sim.results.yearvec, self.results[rkey], label=self.results[rkey].label, **plot_kw)
        pl.legend(frameon=False)
        pl.xlabel('Year')
        pl.ylabel('Number of people')
        sc.boxoff()
        sc.commaticks()
        return fig
measles = SEIR()

# Interventions -------------------------------------------------------------------------------------------

my_vax1 = measles_vaccine(name='vax1', pars=dict(efficacy=0.85))
my_vax2 = measles_vaccine(name='vax2', pars=dict(efficacy=0.99))

intv1 = measles_routine_vx(
    name='routine1', 
    start_year=2020,
    product=my_vax1,
    prob=.95,
    dose = "mcv1"
    
)

intv2 = measles_routine_vx(
    name='routine2',
    start_year=2020,
    product=my_vax2,
    prob=.95,
    dose = "mcv2"
)

intv = [intv1, intv2]

# Parameters -------------------------------------------------------------------------------------------

pars = dict(
    n_agents = 30_000,     # Number of agents to simulate
    total_pop = 48e6,
    birth_rate = 27.58,    # parameters for monthly: birth rate 2022 is 27.58
    death_rate = 7.8,        # parameters for monthly: death rate 2022 is 7.8
    networks = ss.RandomNet(pars={'n_contacts': 4})
)

ppl = ss.People(
    n_agents = 30_000,     # Number of agents to simulate
    age_data = pop_age
   
    )

# Simulations -------------------------------------------------------------------------------------------

mysim = ss.Sim(
    pars = pars,
    start = 2020, 
    people = ppl, 
    diseases = measles, 
    rand_seed = 765,
    n_years = 20,
    dt = 1/12
)

mysim_Intv = ss.Sim(
    pars = pars,
    start = 2020, 
    people = ppl, 
    diseases = measles, 
    interventions = intv,
    rand_seed = 765,
    n_years = 20,
    dt = 1/12
)

# Running the model -------------------------------------------------------------------------------------------

mysim.run()
mysim_Intv.run()
mysim.plot()
plt.show()

# Consolidating the results -------------------------------------------------------------------------------------------

res = pd.DataFrame({
    'year': mysim.yearvec,
    'susceptible': mysim.results.seir.n_susceptible,
    "Exposed": mysim.results.seir.n_exposed,
    "Infected": mysim.results.seir.n_infected,
    "new_infections": mysim.results.seir.new_infections,
    "Recovered": mysim.results.seir.n_recovered,
    "No Intv Prevalence":mysim.results.seir.prevalence,
    "Intv Prevalence":mysim_Intv.results.seir.prevalence
})
res_long = pd.melt(res, id_vars=['year'], var_name='state', value_name='count')

# Filter for specific states
filtered_res = res_long.loc[res_long['state'].isin(['Intv Prevalence', 'No Intv Prevalence'])]

# New infections
res2 = pd.DataFrame({
    'year': mysim.yearvec,
    "new_infections": mysim.results.seir.new_infections,
    "new_infections_intv": mysim_Intv.results.seir.new_infections,
})

res2_long = pd.melt(res2, id_vars=['year'], var_name='state', value_name='count')

# Age at vaccination
age_vaccinated = pd.DataFrame({
    "age_mcv1": mysim_Intv.interventions.routine1.age_at_vaccination*12,
    "age_mcv2" : mysim_Intv.interventions.routine2.age_at_vaccination*12
    
    })


# Plotting -------------------------------------------------------------------------------------------

(
 ggplot(res_long, aes(x='year', y='count')) +
 geom_line(aes(color='state'), size = 1) +
 theme_light()
)


(
  ggplot(res2_long, aes(x="year", y = "count", color = "state")) +
   geom_line(size = 1) +
 theme_light()
)

(
  ggplot(filtered_res, aes(x="year", y = "count", color = "state")) +
  geom_line(size = 1) +
  theme_light() +
  labs(title = "Prevalence")
)


(
  ggplot(age_vaccinated, aes(x = "age_mcv2")) +
  geom_histogram()
)

mysim_Intv.interventions.routine2.n_doses.max()
