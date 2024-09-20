#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 16:25:15 2024

@author: macuser
"""

import starsim as ss
import matplotlib.pyplot as plt


class MyMeasles(ss.Measles):
    pass


class measles_vaccine(ss.sir_vaccine):
    """
    Create a vaccine product that changes susceptible people to recovered (i.e., perfect immunity)
    """

    def administer(self, people, uids):
        people.mymeasles.rel_sus[uids] *= 1-self.pars.efficacy
        return


pars = sc.objdict(
    n_agents=1000,
    birth_rate=20,
    death_rate=15,
    networks='random',
)


measles = MyMeasles()

my_vax = measles_vaccine(pars=dict(efficacy=0.5))
intv = ss.routine_vx(name='vax1', start_year=2020, prob=0.4, product=my_vax)
sim = ss.Sim(pars=pars, diseases=measles, interventions=intv)
sim.run()
sim.plot()
plt.show()
