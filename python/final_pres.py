
# Importing modules
import sciris as sc
import numpy as np
import starsim as ss
import matplotlib.pyplot as plt
import pandas as pd
import pylab as pl
from plotnine import *

# Creating the measles class
class Measles(ss.SIR):
    def __init__(self, pars=None, par_dists=None, *args, **kwargs):
        """ Initialize with parameters """

        pars = ss.omergeleft(pars,
            # Natural history parameters, all specified in days
            dur_exp = 10,       # (days) - source: US CDCv3
            dur_inf = 8,      # 4 days before rash and 4 days after rash
            p_death = 0.002,   # Probability of death

            # Initial conditions and beta
            init_prev = None, # 0.000495,
            imm_prob = None,  # Could replace this with a function that depends on age
            beta = None,
        )
        # initial immunity turkana: .578
        # initial immunity kisumu: .75
        # initial immunity makueni: .8
        

        par_dists = ss.omergeleft(par_dists,
            dur_exp   = ss.normal,
            dur_inf   = ss.normal,
            init_prev = ss.bernoulli,
            #imm_prob  = ss.bernoulli,
            p_death   = ss.bernoulli,
        )

        super().__init__(pars=pars, par_dists=par_dists, *args, **kwargs)

        # SIR are added automatically, here we add E
        self.add_states(
            ss.State('exposed', bool, False),
            ss.State('ti_exposed', float, np.nan),
        )
    

        return
    
    # initial immunity
    # def set_imm_prob(self, sim):
    #     imm_prob = self.imm_prob # ... .578
    #     return imm_prob

    @property
    def infectious(self):
        return self.infected | self.exposed
    
    def set_initial_states(self, sim):
        super().set_initial_states(sim)
        imm_prob = self.pars.imm_prob
        alive_uids = ss.true(sim.people.alive & self.susceptible)
        initially_immune = imm_prob > np.random.rand()

     # initially_immune = self.pars.init_imm.filter(alive_uids)
        self.susceptible[initially_immune] = False
        self.recovered[initially_immune] = True
        

    def update_pre(self, sim):
        # Progress exposed -> infected
        infected = ss.true(self.exposed & (self.ti_infected <= sim.ti)) #  & sim.people.age < 20
        self.exposed[infected] = False
        self.infected[infected] = True

        # Progress infected -> recovered
        recovered = ss.true(self.infected & (self.ti_recovered <= sim.ti))
        self.infected[recovered] = False
        self.recovered[recovered] = True

        # Trigger deaths
        deaths = ss.true(self.ti_dead <= sim.ti)
        if len(deaths):
            sim.people.request_death(deaths)
        return

    def set_prognoses(self, sim, uids, source_uids=None):
        """ Set prognoses for those who get infected """
        # Do not call set_prognosis on parent
        # super().set_prognoses(sim, uids, source_uids)
        

        self.susceptible[uids] = False
        self.exposed[uids] = True
        self.ti_exposed[uids] = sim.ti

        p = self.pars

        # Determine when exposed become infected
        self.ti_infected[uids] = sim.ti + p.dur_exp.rvs(uids) / sim.dt

        # Sample duration of infection, being careful to only sample from the
        # distribution once per timestep.
        dur_inf = p.dur_inf.rvs(uids)

        # Determine who dies and who recovers and when
        will_die = p.p_death.rvs(uids)
        dead_uids = uids[will_die]
        rec_uids = uids[~will_die]
        self.ti_dead[dead_uids] = self.ti_infected[dead_uids] + dur_inf[will_die] / sim.dt
        self.ti_recovered[rec_uids] = self.ti_infected[rec_uids] + dur_inf[~will_die] / sim.dt

        return

    def update_death(self, sim, uids):
        # Reset infected/recovered flags for dead agents
        for state in ['susceptible', 'exposed', 'infected', 'recovered']:
            self.statesdict[state][uids] = False
        return
    
# Creating measles vaccine
class measles_vaccine(ss.sir_vaccine):
    """
    Create a vaccine product that changes susceptible people to recovered (i.e., perfect immunity)
    """
    def administer(self, people, uids):
        people.measles.rel_sus[uids] *= 1-self.pars.efficacy
        return
    
# Routine measles vaccine
class routine_measles_vx(ss.routine_vx):
    pass



# The function for simulating
def county_sim(pars_df, scs, **kwargs):
    all_new_infections = {}  # Dictionary to accumulate data for all counties and scs combinations
    
    pars = {}
    
    eligibility_mcv2 = lambda sim: (sim.interventions.routine1.n_doses == 1)
    my_vax1 = measles_vaccine(name='vax1', pars=dict(efficacy=0.85))
    my_vax2 = measles_vaccine(name='vax2', pars=dict(efficacy=0.99))
    my_vax3 = measles_vaccine(name='vax3', pars=dict(efficacy=0.95))
    
    for _, county_row in pars_df.iterrows():
        county_name = county_row['county']
        all_new_infections[county_name] = {}  # Create a sub-dictionary for the county
        
        pars[county_name] = sc.objdict(
            n_agents = 5000,
            birth_rate = county_row['birth_rate'], 
            death_rate = county_row['death_rate'],
            networks = dict(
                type = 'randomnet',
                n_contacts = 4 
            )
        )
        
        for _, sc_row in scs.iterrows():
            # The interventions from the vaccine
            intv1 = routine_measles_vx(
                name='routine1', 
                start_year=2020,
                product=my_vax1,
                prob=sc_row['mcv1']
            )
            intv2 = routine_measles_vx(
                name='routine2',
                start_year=2020,
                eligibility=eligibility_mcv2,
                product=my_vax2,
                prob=sc_row['mcv2']
            )
            
            # The interventions
            intv = [intv1, intv2]
            
            # The simulations
            print(f"\033[92mRunning simulations for county: {county_name} with scs: {sc_row}\033[0m")
            
            sim_intv = ss.Sim(
                pars=pars[county_name],
                diseases=Measles(beta=.9, init_prev=county_row['initial_prev'], imm_prob=county_row['initial_immunity']), 
                interventions=intv, 
                start=2020, 
                end=2050
            )
            sim_intv.run()
            
            # Collecting data for data_new_infections
            data_new_infections = sim_intv.results.measles.new_infections
            df_new_infections = pd.DataFrame({'new_infections': data_new_infections}, index=sim_intv.yearvec)  # Add index
            all_new_infections[county_name][tuple(sc_row)] = df_new_infections  # Store data for this scs combination
            
    return all_new_infections

# res = county_sim(pars_df, scs)
res_with_immunity = county_sim(pars_df, scs)


# For kenya at large
eligibility_mcv2 = lambda sim: (sim.interventions.routine1.n_doses == 1)
my_vax1 = measles_vaccine(name='vax1', pars=dict(efficacy=0.85))
my_vax2 = measles_vaccine(name='vax2', pars=dict(efficacy=0.99))
my_vax3 = measles_vaccine(name='vax3', pars=dict(efficacy=0.95))

pars = sc.objdict(
n_agents = 5000,
birth_rate = 21.3, 
death_rate = 6.2,
networks = dict(
    type = 'randomnet',
    n_contacts = 4 
)
)

intv1 = routine_measles_vx(
    name='routine1', 
    start_year=2020,
    product=my_vax1,
    prob=.95
)
intv2 = routine_measles_vx(
    name='routine2',
    start_year=2020,
    eligibility=eligibility_mcv2,
    product=my_vax2,
    prob=.95
)

# np.arange(2021, 2030.5, 0.5)
# 3y = 2023, 2026, 2029
# 2y 2022, 2024, 2026, 2028, 2030 # .85 original
sia = ss.campaign_vx(name='SIA', years=np.arange(2021, 2030.5, 2) , prob=0.95, product = my_vax3) # The campaigns representing the SIAs


# The interventions
intv = [intv1, intv2, sia]

# The simulations
#print(f"\033[92mRunning simulations for county: {county_name} with scs: {sc_row}\033[0m")

sim_intv = ss.Sim(
    pars=pars,
    diseases=Measles(beta=.9, init_prev=0.0004945954, imm_prob=0),  # 
    interventions=intv, 
    start=2020, 
    end=2050
)
sim_intv.run()

# Collecting data for data_new_infections
data_new_infections = sim_intv.results.measles.new_infections
df_new_infections = pd.DataFrame({'new_infections': data_new_infections}, index=sim_intv.yearvec)  # Add index
all_new_infections = df_new_infections  # Store data for this scs combination

sim_intv.plot()
plt.show()
# plotting
import matplotlib.pyplot as plt
plt.figure()
plt.plot(sim_intv.yearvec, sim_intv.results.measles.prevalence, label='Baseline')
plt.plot(sim_intv.yearvec, sim_intv.results.measles.prevalence, label='Vax')
#plt.axvline(x=2015, color='k', ls='--')
plt.title('Prevalence')
plt.legend()
plt.show();


