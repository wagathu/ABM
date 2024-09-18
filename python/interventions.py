
import starsim as ss
pars = dict(
    n_agents = 5000,
    # birth_rate = 20,
    # death_rate = 15,
    networks = dict(
        type = 'randomnet',
        n_contacts = 4
    ),
    diseases = dict(
        type = 'sir',
        dur_inf = 10,
        beta = 0.1,
    )
)

# Create the product - a vaccine with 50% efficacy
MyVaccine1 = ss.sir_vaccine(pars=dict(efficacy=0.7), name='vax1')
MyVaccine2 = ss.sir_vaccine(pars=dict(efficacy=0.2), name='vax2') # You can add other vaccine -specific paramters


# Create the intervention
MyIntervention1 = ss.routine_vx(
    name = 'routine1',
    start_year = 2013,    # Begin vaccination in 2015
    prob = 0.6,           # 20% coverage
    product = MyVaccine1   # Use the MyVaccine product
)

MyIntervention2 = ss.routine_vx(
    name = 'routine2',
    start_year = 2014,    # Begin vaccination in 2015
    prob = 0.2,           # 20% coverage
    product=MyVaccine2   # Use the MyVaccine product
)




# Now create two sims: a baseline sim and one with the intervention
sim_base = ss.Sim(pars=pars)
sim_base = sim_base.run()
sim_intv = ss.Sim(pars=pars, interventions=[MyIntervention1, MyIntervention2])
sim_intv.run()

import matplotlib.pyplot as plt
plt.figure()
plt.plot(sim_base.yearvec, sim_base.results.sir.prevalence, label='Baseline')
plt.plot(sim_intv.yearvec, sim_intv.results.sir.prevalence, label='Vax')
plt.axvline(x=2015, color='k', ls='--')
plt.title('Prevalence')
plt.legend()
plt.show();
