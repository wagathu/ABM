#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 16:25:15 2024

@author: macuser
"""
# Importing the Modules
import sciris as sc
import starsim as ss
import matplotlib.pyplot as plt

# Creating the class measles
class Measles(SIR):

    def __init__(self, pars=None, par_dists=None, *args, **kwargs):
        """ Initialize with parameters """

        pars = ss.omergeleft(pars,
            # Natural history parameters, all specified in days
            dur_exp = 8,       # (days) - source: US CDC
            dur_inf = 11,      # (days) - source: US CDC
            p_death = 0.005,   # Probability of death

            # Initial conditions and beta
            init_prev = 0.005,
            beta = None,
        )

        par_dists = ss.omergeleft(par_dists,
            dur_exp   = ss.normal,
            dur_inf   = ss.normal,
            init_prev = ss.bernoulli,
            p_death   = ss.bernoulli,
        )

        super().__init__(pars=pars, par_dists=par_dists, *args, **kwargs)

        # SIR are added automatically, here we add E
        self.add_states(
            ss.State('exposed', bool, False),
            ss.State('ti_exposed', float, np.nan),
        )

        return

    @property
    def infectious(self):
        return self.infected | self.exposed

    def update_pre(self, sim):
        # Progress exposed -> infected
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
        self.ti_dead[dead_uids] = self.ti_infected[dead_uids] + dur_inf[will_die] / sim.dt
        self.ti_recovered[rec_uids] = self.ti_infected[rec_uids] + dur_inf[~will_die] / sim.dt

        return

    def update_death(self, sim, uids):
        # Reset infected/recovered flags for dead agents
        for state in ['susceptible', 'exposed', 'infected', 'recovered']:
            self.statesdict[state][uids] = False
        return

# Creating the class for the two doses
class measles_vaccine(ss.sir_vaccine):
    """
    Create a vaccine product that changes susceptible people to recovered (i.e., perfect immunity)
    """
    def administer(self, people, uids):
        people.measles.rel_sus[uids] *= 1-self.pars.efficacy
        return
      
# Creating the class for the SIAs
class measles_sia(ss.campaign_vx):
    """
    Create a vaccine product that changes susceptible people to recovered (i.e., perfect immunity)
    """
    def administer(self, people, uids):
        people.measles.rel_sus[uids] *= 1-self.pars.efficacy
        return

# The parameters
pars = sc.objdict(
    n_agents = 5_000,
    birth_rate = 20,
    death_rate = 15,
    networks = dict(
        type = 'randomnet',
        contacts = 4
        ),
    diseases = dict(
        name = 'measles',
        beta = 0.2
        )
    )

# The vaccines
my_vax1 = measles_vaccine(pars=dict(efficacy=0.7))
my_vax2 = measles_vaccine(pars=dict(efficacy=0.8))
my_vax3 = measles_vaccine(pars=dict(efficacy=0.8))

# The interventions from the vaccine
intv1 = ss.routine_vx(name='vax1', start_year=2000, prob=0, product=my_vax1)
intv2 = ss.routine_vx(name='vax1', start_year=2013, prob=0.8, product=my_vax2)
sia = ss.campaign_vx(name='SIA', years= [2015, 2017, 2019, 2021], prob=0.8, product=my_vax3) # The campaigns representing the SIAs
intv = [intv1, intv2]

# The simulations
sim_base = ss.Sim(pars=pars)
sim_intv = ss.Sim(pars = pars, interventions = intv)
sim_base.run()
sim_intv.run()

# Plotting
plt.figure()
plt.plot(sim_base.yearvec, sim_base.results.measles.prevalence, label='Baseline')
plt.plot(sim_intv.yearvec, sim_intv.results.measles.prevalence, label='Vax')
plt.axvline(x=2013, color='k', ls='--')
plt.title('Prevalence')
plt.legend()
plt.show();





