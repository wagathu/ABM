
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 16:25:15 2024

@author: macuser
"""
# Importing the Modules
import sciris as sc
import numpy as np
import starsim as ss
import matplotlib.pyplot as plt

# Creating the class measles
class Measles2(ss.SIR):

    def __init__(self, pars=None, par_dists=None, *args, **kwargs):
        """ Initialize with parameters """

        pars = ss.omergeleft(pars,
            # Natural history parameters, all specified in days
            dur_exp = 8,       # (days) - source: US CDCv3
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
        
        self.age_bins = [0, 1, 5]

        return

    @property
    def infectious(self):
        return self.infected | self.exposed

    def update_pre(self, sim):
        # Progress exposed -> infected
        infected = ss.true(self.exposed & (self.ti_infected <= sim.ti) & sim.people.age < 20)
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
    
    def get_age_bins(self, sim, uids):
        ages = sim.people.age[uids]
        age_bins = np.digitize(ages, self.age_bins)
        return age_bins

    def set_prognoses(self, sim, uids, source_uids=None):
        """ Set prognoses for those who get infected """
        # Do not call set_prognosis on parent
        # super().set_prognoses(sim, uids, source_uids)
        
        age_bins = self.get_age_bins(sim, uids)

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

class routine_measles_vx(ss.routine_vx):
    """
    Routine vaccination - an instance of base vaccination combined with routine delivery.
    See base classes for a description of input arguments.
    """

    def __init__(self, product=None, prob=None, eligibility=None, mean = None,
                 start_year=None, end_year=None, years=None, uid = None, age_lim1=9/12, age_lim2=12/12, **kwargs):
        super().__init__(product=product, prob=prob, eligibility=eligibility, 
                         start_year=start_year, end_year=end_year, years=years, **kwargs)
        self.age_lim1 = age_lim1
        self.age_lim2 = age_lim2
        self.uid = uid,
        self.age = self.people.age
        return

    def check_eligibility(self, sim):
        eligible = (sim.people.age >= self.age_lim1) & (sim.people.age <= self.age_lim2)
        #eligible_uid = sim.people.uid[eligible]
        return eligible
    
    def probs_vac(self, mean):
        ages = self.people.age
        prob = np.random.normal(loc=ages/mean, scale=0.2)
        prob = np.abs(prob)
        if prob > 1:
            prob = 1/prob
        return prob
        


#measles = measles(beta = .6)

# The parameters
pars = sc.objdict(
    n_agents = 5000,
    dt = 1/12,
    birth_rate = 27,
    death_rate = 8,
    networks = dict(
        type = 'randomnet',
        n_contacts = 10
        ),
    diseases = dict(
        type = 'measles',
        
        )
    )

# The vaccines
my_vax1 = measles_vaccine(name='vax1', pars=dict(efficacy=0.85))
my_vax2 = measles_vaccine(name='vax2', pars=dict(efficacy=0.95))
my_vax3 = measles_vaccine(name='vax3', pars=dict(efficacy=0.9))



# The interventions from the vaccine
intv1 =  routine_measles_vx(name='routine1', 
                            start_year=2009,
                            #prob = probs_vac(self, mean = 12),
                            product=my_vax1,
                            age_lim1=9/12,
                            age_lim2=12/12
                            )
                            
intv2 =  routine_measles_vx(name='routine2',
                            start_year=2006, 
                            #prob = probs_vac(self, mean = 21),
                            product=my_vax2,
                            age_lim1=18/12,
                            age_lim2=24/12
                            )
                            
sia = ss.campaign_vx(name='SIA', years= [2006], prob=0.9, product=my_vax3) # The campaigns representing the SIAs
intv = [intv1, intv2, sia]

# The simulations
sim_base = ss.Sim(pars=pars, start = 2005, end = 2030)
sim_intv = ss.Sim(pars = pars, interventions = intv, start = 2005, end = 2030)
sim_base.run()
sim_intv.run()

# Plotting
plt.figure()
plt.plot(sim_base.yearvec, sim_base.results.measles.prevalence, label='Baseline')
plt.plot(sim_intv.yearvec, sim_intv.results.measles.prevalence, label='Vax')
plt.axvline(x=2013, color='k', ls='--')
plt.axvline(x=2010, color='k', ls='--')

plt.title('Prevalence')
plt.legend()
plt.show();


sim_base.plot()

