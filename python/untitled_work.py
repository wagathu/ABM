
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import starsim as ss
from starsim.diseases.sir import SIR

__all__ = ['Measles']



class Measles(SIR):
    def __init__(self, pars=None, *args, **kwargs):
        """Initialize with parameters"""
        super().__init__(*args, **kwargs)

        # Assuming self.default_pars is a method that sets default parameters
        self.default_pars(
            beta=1e3,  # Placeholder value for transmission rate
            init_prev=ss.bernoulli(p=0.005),  # Initial prevalence
            imm_prob=0.5,  # Immunity probability

            # Natural history parameters, all specified in days
            dur_exp=ss.normal(loc=8/365),        # Incubation period (days)
            dur_inf=ss.normal(loc=11/365),       # Infectious period (days)
            p_death=ss.bernoulli(p=0.005),   # Probability of death
        )

        # Update parameters with provided values
        self.update_pars(pars=pars, **kwargs)

        # Add additional states for the Measles model
        self.add_states(
            ss.BoolArr('exposed', label='Exposed'),
            ss.FloatArr('ti_exposed', label='Time of exposure'),
        )

        return

    @property
    def infectious(self):
        return self.infected | self.exposed
    
    def set_initial_states(self, sim):
        super().set_initial_states(sim)
        imm_prob = self.default_pars().imm_prob
        alive_uids = ss.true(sim.people.alive & self.susceptible)
        initially_immune = imm_prob > np.random.rand()

     # initially_immune = self.pars.init_imm.filter(alive_uids)
        self.susceptible[initially_immune] = False
        self.recovered[initially_immune] = True


    def update_pre(self):
        # Progress exposed -> infected
        ti = self.sim.ti
        infected = (self.exposed & (self.ti_infected <= ti)).uids
        self.exposed[infected] = False
        self.infected[infected] = True

        # Progress infected -> recovered
        recovered = (self.infected & (self.ti_recovered <= ti)).uids
        self.infected[recovered] = False
        self.recovered[recovered] = True

        # Trigger deaths
        deaths = (self.ti_dead <= ti).uids
        if len(deaths):
            self.sim.people.request_death(deaths)
        return

    def set_prognoses(self, uids, source_uids=None):
        """ Set prognoses for those who get infected """
        super().set_prognoses(uids, source_uids)
        ti = self.sim.ti
        dt = self.sim.dt

        self.susceptible[uids] = False
        self.exposed[uids] = True
        self.ti_exposed[uids] = ti

        p = self.pars

        # Determine when exposed become infected
        self.ti_infected[uids] = ti + p.dur_exp.rvs(uids) / dt

        # Sample duration of infection, being careful to only sample from the
        # distribution once per timestep.
        dur_inf = p.dur_inf.rvs(uids)

        # Determine who dies and who recovers and when
        will_die = p.p_death.rvs(uids)
        dead_uids = uids[will_die]
        rec_uids = uids[~will_die]
        self.ti_dead[dead_uids] = self.ti_infected[dead_uids] + dur_inf[will_die] / dt
        self.ti_recovered[rec_uids] = self.ti_infected[rec_uids] + dur_inf[~will_die] / dt

        return

    def update_death(self, uids):
        # Reset infected/recovered flags for dead agents
        for state in ['susceptible', 'exposed', 'infected', 'recovered']:
            self.statesdict[state][uids] = False
        return

pars = sc.objdict(
    n_agents=5000,
    birth_rate=0,
    death_rate=0,
    networks=ss.RandomNet()
)

measles = Measles()
sim = ss.Sim(pars = pars, diseases = measles )
sim.run()
sim.results

