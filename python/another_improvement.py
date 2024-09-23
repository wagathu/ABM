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

# Measles class
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
            beta = .8, # 1 - np.exp(-(15/7)),  # Mean transmission rate: R0 / D
            
            init_prev = ss.bernoulli(p=.05),
            dur_exp = ss.lognorm_ex(mean=10/30, stdev=2),
            dur_inf = ss.lognorm_ex(mean=10/30, stdev=2),
            p_death = ss.bernoulli(p=0.018),
            p_immunity = ss.bernoulli(p = .0)

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
        return self.infected | self.exposed


    def update_pre(self):
        """ Update states before the next time step """
        sim = self.sim
        p = self.pars
        ti = sim.ti
        dt = sim.dt
        
        # handle beta here: at start of infection prior to transmission
        beta_mean = .3  # Mean transmission rate
        beta_amplitude = 1  # Amplitude of seasonal forcing
        #beta_phase = 100  # Phase shift of seasonal forcing (in months)
        #beta_period = 20  # Period of seasonal forcing (in months)
        beta_rate = beta_mean * (1 + beta_amplitude * np.cos(2 * np.pi * ti/365) * np.sin(2 * np.pi * ti/365))
        beta_prob = 1 - np.exp(-beta_rate)
        
        # Dynamically get the network keys from the simulation
        network_keys = self.sim.networks.keys()
        self.pars.beta = {key: [beta_prob, beta_prob] for key in network_keys}
        
        # add a condition here to make children less than/eq to 6 months of age immune
        p_immunity = np.where((sim.people.age < .5) | (sim.people.age > 50), 1, 0)

        
        # conditions for all older people above 20 years to never get measles
        all_ids_above_20 = sim.people.uid[sim.people.age > 20]
        self.susceptible[all_ids_above_20] = False
        self.recovered[all_ids_above_20] = True
        
        # Handle initially immune individuals
        if ti == 0:  
            all_susceptibles = self.susceptible.uids
            initially_immune = p.p_immunity.rvs(len(all_susceptibles))
            immune_uids = all_susceptibles[initially_immune]
            self.susceptible[immune_uids] = False
            self.recovered[immune_uids] = True
            self.ti_recovered[immune_uids] = ti
            

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
# Pars
pars = dict(
    n_agents = 20_000,     # Number of agents to simulate
    birth_rate = 27.58,    # parameters for monthly: birth rate 2022 is 27.58
    death_rate = 8,        # parameters for monthly: death rate 2022 is 7.8
    networks = ss.RandomNet(pars={'n_contacts': 40})
)

ppl = ss.People(
    20_000, 
    age_data = pop_age
    )

mysim = ss.Sim(
    pars = pars, 
    total_pop = 47563609, 
    start = 2020, 
    people = ppl, 
    diseases = measles, 
    rand_seed = 765,
    n_years = 5, 
    dt=1/365
)

mysim.run()
mysim.plot()
plt.show()


res = pd.DataFrame({
    'year': mysim.yearvec,
    #'susceptible': mysim.results.seir.n_susceptible,
    #"Exposed": mysim.results.seir.n_exposed,
    "Infected": mysim.results.seir.n_infected,
    #"Recovered": mysim.results.seir.n_recovered
})
res_long = pd.melt(res, id_vars=['year'], var_name='state', value_name='count')

(
 ggplot(res_long.iloc[10:, ], aes(x='year', y='count')) +
 geom_line(aes(color='state')) +
 geom_point(aes(color='state')) +
 theme_light()
)

# Age population
x = mysim.people.age
x_df = pd.DataFrame({"age":x})
(
 ggplot(x_df, aes(x = "age")) +
 geom_histogram(color = "white")
)
