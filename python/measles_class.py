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
class MyMeasles(ss.Measles):
    pass


# Creating the class for the two doses
class measles_vaccine(ss.sir_vaccine):
    """
    Create a vaccine product that changes susceptible people to recovered (i.e., perfect immunity)
    """
    def administer(self, people, uids):
        people.mymeasles.rel_sus[uids] *= 1-self.pars.efficacy
        return
      
# Creating the class for the SIAs
class measles_sia(ss.campaign_vx):
    """
    Create a vaccine product that changes susceptible people to recovered (i.e., perfect immunity)
    """
    def administer(self, people, uids):
        people.mymeasles.rel_sus[uids] *= 1-self.pars.efficacy
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
        name = 'mymeasles',
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
plt.plot(sim_base.yearvec, sim_base.results.mymeasles.prevalence, label='Baseline')
plt.plot(sim_intv.yearvec, sim_intv.results.mymeasles.prevalence, label='Vax')
plt.axvline(x=2013, color='k', ls='--')
plt.title('Prevalence')
plt.legend()
plt.show();





