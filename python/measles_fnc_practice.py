
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
                             dur_exp=10,       # (days) - source: US CDCv3
                             dur_inf=8,      # 4 days before rash and 4 days after rash
                             p_death=0.002,   # Probability of death

                             # Initial conditions and beta
                             init_prev=None,  # 0.000495,
                             imm_prob=None,  # Could replace this with a function that depends on age
                             beta=None,
                             )
        # initial immunity turkana: .578
        # initial immunity kisumu: .75
        # initial immunity makueni: .8

        par_dists = ss.omergeleft(par_dists,
                                  dur_exp=ss.normal,
                                  dur_inf=ss.normal,
                                  init_prev=ss.bernoulli,
                                  # imm_prob  = ss.bernoulli,
                                  p_death=ss.bernoulli,
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
        # & sim.people.age < 20
        infected = ss.true(self.exposed & (self.ti_infected <= sim.ti))
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
        self.ti_dead[dead_uids] = self.ti_infected[dead_uids] + \
            dur_inf[will_die] / sim.dt
        self.ti_recovered[rec_uids] = self.ti_infected[rec_uids] + \
            dur_inf[~will_die] / sim.dt

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
    # Dictionary to accumulate data for all counties and scs combinations
    all_new_infections = []

    pars = {}

    def eligibility_mcv2(sim): return (sim.interventions.routine1.n_doses == 1)
    my_vax1 = measles_vaccine(name='vax1', pars=dict(efficacy=0.85))
    my_vax2 = measles_vaccine(name='vax2', pars=dict(efficacy=0.99))
    my_vax3 = measles_vaccine(name='vax3', pars=dict(efficacy=0.95))

    for _, county_row in pars_df.iterrows():
        county_name = county_row['county']

        pars[county_name] = sc.objdict(
            n_agents=5000,
            birth_rate=county_row['birth_rate'],
            death_rate=county_row['death_rate'],
            networks=dict(
                type='randomnet',
                n_contacts=4
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
            print(
                f"\033[92mRunning simulations for county: {county_name} with scs: {sc_row}\033[0m")

            sim_intv = ss.Sim(
                pars=pars[county_name],
                diseases=Measles(
                    beta=.9, init_prev=county_row['initial_prev'], imm_prob=county_row['initial_immunity']),
                interventions=intv,
                start=2020,
                end=2050
            )
            sim_intv.run()

            # Collecting data for data_new_infections
            # Collecting results: number of S, I, R, and new infections
            data_s = sim_intv.results.measles.n_susceptible
            data_i = sim_intv.results.measles.n_infected
            data_r = sim_intv.results.measles.n_recovered
            data_new_infections = sim_intv.results.measles.new_infections

            # Creating a DataFrame with S, I, R, and new infections
            df_results = pd.DataFrame({
                'year': sim_intv.yearvec,
                'susceptible': data_s,
                'infected': data_i,
                'recovered': data_r,
                'new_infections': data_new_infections,
                'county': county_name,
                'mcv1': sc_row['mcv1'],
                'mcv2': sc_row['mcv2']
            })

            # Accumulating results
            all_new_infections.append(df_results)
        final_results = pd.concat(all_new_infections).reset_index(drop=True)
    return final_results


# res = county_sim(pars_df, scs)
res_with_immunity = county_sim(pars_df, scs)

pars_df = pd.read_csv("/Users/macuser/Documents/GitHub/ABM/python/pars_df.csv")
# Extend the range to just above 1, so it includes 0.95
mcv2_values = np.arange(0, 1.05, 0.1).tolist()
mcv2_values = [min(v, 0.95) for v in mcv2_values]  # Limit the values to 0.95


# Create the DataFrame
scs = pd.DataFrame({
    'mcv1': [0.95] * len(mcv2_values),  # Fill the mcv1 column with 0.95
    'mcv2': mcv2_values                # Use the generated values for mcv2
})
