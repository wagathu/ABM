import starsim as ss
import sciris as sc

class measlesIntervention(ss.Plugin):
    """
    Base class for interventions.
    
    The key method of the measlesIntervention is ``apply()``, which is called with the sim
    on each timestep.
    """

    def __init__(self, eligibility=None, dose = None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.eligibility = eligibility
        self.dose = dose
        return

    def _parse_product(self, product):
        """
        Parse the product input
        """
        if isinstance(product, ss.Product):  # No need to do anything
            self.product = product
        elif isinstance(product, str):
            self.product = self._parse_product_str(product)
        else:
            errormsg = f'Cannot understand {product} - please provide it as a Product.'
            raise ValueError(errormsg)
        return

    def _parse_product_str(self, product):
        raise NotImplementedError

    def check_eligibility(self, sim):
        """
        Return an array of indices of agents eligible for screening at time t
        """
        if self.eligibility is not None:
            is_eligible = self.eligibility(sim)
            if is_eligible is not None and len(is_eligible): # Only worry if non-None/nonzero length
                if isinstance(is_eligible, ss.BoolArr):
                    is_eligible = is_eligible.uids
                if not isinstance(is_eligible, ss.uids):
                    errormsg = f'Eligibility function must return BoolArr or UIDs, not {type(is_eligible)} {is_eligible}'
                    raise TypeError(errormsg)
        else:
            is_eligible = sim.people.auids # Everyone
        return is_eligible
    # dummy function for checking dose
    def check_dose(self, sim):
        if self.dose == "mcv1":
            dose = "mcv1"
        else:
            dose = "mcv2"
        return dose
    
class measlesBaseVaccination(measlesIntervention):
    """
    Base vaccination class for determining who will receive a vaccine.
    """
    def __init__(self, product=None, prob=None, label=None, **kwargs):
        super().__init__(**kwargs)
        self.prob = sc.promotetoarray(prob)
        self.label = label
        self._parse_product(product)
        self.vaccinated = ss.BoolArr('vaccinated')
        self.n_doses = ss.FloatArr('doses', default=0)
        self.ti_vaccinated = ss.FloatArr('ti_vaccinated')
        self.age_at_vaccination = ss.FloatArr('age_at_vaccination')  # New array to store age at vaccination
        self.coverage_dist = ss.bernoulli(p=0)  # Placeholder
        return

    def apply(self, sim):
        """
        Deliver the diagnostics by finding who's eligible, finding who accepts, and applying the product.
        """
        accept_uids = np.array([])
        if sim.ti in self.timepoints:

            ti = sc.findinds(self.timepoints, sim.ti)[0]
            prob = self.prob[ti]  # Get the proportion of people who will be tested this timestep
            #is_eligible = self.check_eligibility(sim)  # Check eligibility
            #self.coverage_dist.set(p=prob)
            #accept_uids = self.coverage_dist.filter(is_eligible)
            
            check_dose = self.check_dose(sim)
            if check_dose == "mcv1":
                eligible_accept_uids = sim.people.auids[(sim.people.age >= 9/12) & (sim.people.age <= 12/12) & (sim.interventions.routine1.n_doses == 0)]
            else:
                eligible_accept_uids = sim.people.auids[
                    (sim.people.age >= 18/12) & 
                    (sim.people.age <= 24/12) & 
                    (sim.interventions.routine1.n_doses == 1) &
                    (sim.interventions.routine2.n_doses == 0) 
                ]
                
            accept_uids = eligible_accept_uids[ss.bernoulli(p=prob, strict = False).rvs(eligible_accept_uids.shape)]

            if len(accept_uids):
                self.product.administer(sim.people, accept_uids)

                # Update people's state and dates
                self.vaccinated[accept_uids] = True
                self.ti_vaccinated[accept_uids] = sim.ti
                self.n_doses[accept_uids] += 1
                
                self.age_at_vaccination[accept_uids] = sim.people.age[accept_uids]

        return accept_uids

class measles_routine_vx(measlesBaseVaccination, ss.RoutineDelivery):
    """
    Routine vaccination - an instance of base vaccination combined with routine delivery.
    See base classes for a description of input arguments.
    """

    def __init__(self, product=None, prob=None, eligibility=None,
                 start_year=None, end_year=None, years=None, **kwargs):

        measlesBaseVaccination.__init__(self, product=product, eligibility=eligibility, **kwargs)
        ss.RoutineDelivery.__init__(self, prob=prob, start_year=start_year, end_year=end_year, years=years)
        return

    def init_pre(self, sim):
        ss.RoutineDelivery.init_pre(self, sim)  # Initialize this first, as it ensures that prob is interpolated properly
        measlesBaseVaccination.init_pre(self, sim)  # Initialize this next

class measles_vaccine(ss.Vx):
    """
    Create a vaccine product that affects the probability of infection.
    
    The vaccine can be either "leaky", in which everyone who receives the vaccine 
    receives the same amount of protection (specified by the efficacy parameter) 
    each time they are exposed to an infection. The alternative (leaky=False) is
    that the efficacy is the probability that the vaccine "takes", in which case
    that person is 100% protected (and the remaining people are 0% protected).
    
    Args:
        efficacy (float): efficacy of the vaccine (0<=efficacy<=1)
        leaky (bool): see above
    """
    def __init__(self, pars=None, *args, **kwargs):
        super().__init__()
        self.default_pars(
            efficacy = 0.9,
            leaky = True
        )
        self.update_pars(pars, **kwargs)
        return

    def administer(self, people, uids):        
        if self.pars.leaky:
            people.seir.rel_sus[uids] *= 1-self.pars.efficacy
        else:
            people.seir.rel_sus[uids] *= np.random.binomial(1, 1-self.pars.efficacy, len(uids))
        return
   
