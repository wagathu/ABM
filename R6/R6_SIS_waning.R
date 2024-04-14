if (require(pacman))
  p_load(R6, ggplot2, data.table, purrr, ggalt)

default_pars = list(
  beta = 0.3,                    # Infection rate per contact per unit time
  gamma = 0.5,                   # Recovery rate
  waning_rate = 1.0,             # Rate at which immunity wanes
  n_contacts = 10,               # Number of people each person is connected to
  distance = 0.1,                # The distance over which people form contacts
  I0 = 3,                        # Number of people initially infected
  N = 1000,                      # Total population size
  maxtime = 40,                  # How long to simulate for
  npts = 100,                    # Number of time points during the simulation
  seed = 1,                      # Random seed to use
  colors = list(S = 'darkgreen', I = 'gold', R = 'skyblue'),
  save_movie = FALSE             # Whether to save the movie (slow)
)


# Define the Person class
Person <- R6Class(
  'Person',
  public = list(
    pars = NULL,
    S = NULL,                    # Initialize S and I as public fields
    I = NULL,
    t_I = NULL,
    imm = NULL,
    x = NULL,
    y = NULL,
    
    initialize = \(pars, ...) {
      self$pars = pars
      self$S = TRUE              # Initialize S and I as public fields
      self$I = FALSE
      self$t_I = NULL
      self$imm = 0
      self$x = runif(1)
      self$y = runif(1)
      return(self)
    },
    
    infect = \(t) {
      self$S = F
      self$I = T
      self$t_I = t
      return(self)
    },
    
    recover = \(t) {
      self$I = F
      self$S = T
      self$imm = 1
      return(self)
    },
    
    check_immunity = \(t) {
      pars = self$pars
      if (self$S & !is.null(self$t_I)) {
        time_diff = t - self$t_I
        waning_factor = time_diff / pars$waning_rate
        self$imm = exp(-waning_factor)
      }
      return(self)
    },
    
    check_infection = \(other, t) {
      pars = self$pars
      if (self$S & other$I) {
        susceptibility = 1 - self$imm
        if (runif(1) < (pars$beta * pars$dt * susceptibility)) {
          self = self$infect(t)
        }
      }
      return(self)
    },
    check_recovery = \(t) {
      pars = self$pars
      if (self$I) {
        if (runif(1) < (pars$gamma * pars$dt)) {
          self = self$recover(t)
        }
      }
      return(self)
    }
    
  )
)


# Define the Sim class
Sim <- R6Class(
  'Sim',
  public = list(
    S = NULL,
    I = NULL,
    Ti = NULL,
    time = NULL,
    S_full = NULL,
    I_full = NULL,
    contacts = NULL,
    people = list(),
    pars = NULL,
    
    initialize = \(pars = default_pars, ...) {
      pars = as.list(default_pars, list(...))
      self$pars = pars
      self$pars$dt = self$pars$maxtime / self$pars$npts
      self$Ti =  1:self$pars$npts
      self$time = self$Ti * self$pars$dt
      
      self$initialize_objects()
      return(self)
    },
    
    initialize_objects = \() {
      pars = self$pars
      
      set.seed(pars$seed)
      self$people = self$people = purrr::map(1:pars$N, ~ Person$new(pars))
      self$people[1:pars$I0] = purrr::map(1:pars$I0, ~ self$people[[.x]]$infect(t = 0))
      self$make_network()
      
      self$S = rep(0, pars$npts)
      self$I = rep(0, pars$npts)
      
      self$S_full = list()
      self$I_full = list()
      return(self)
    },
    
    
    get_xy = \() {
      # Get the location of each agent
      x <- map_dbl(self$people, ~ .x$x)
      y <- map_dbl(self$people, ~ .x$y)
      return(list(x = x, y = y))
    },
    
    make_network = function() {
      pars = self$pars
      xy = self$get_xy()
      x = xy$x
      y = xy$y
      dist = matrix(0, nrow = pars$N, ncol = pars$N)
      
      for (i in 1:pars$N) {
        dist[i, ] = 1 + ((x - x[i]) ^ 2 + (y - y[i]) ^ 2) ^ 0.5 / pars$distance
        dist[i, i] = Inf
      }
      
      rnds = matrix(runif(pars$N * pars$N), nrow = pars$N, ncol = pars$N)
      ratios = dist / rnds
      order = order(ratios, na.last = NA)
      inds = order[0:(pars$N * pars$n_contacts / 2)]
      contacts = arrayInd(inds, .dim = dim(ratios))
      self$contacts = t(cbind(contacts))
      return(self)
    },
    
    check_infections = \(t) {
      for (i in 1:ncol(self$contacts)) {
        p1 <- self$contacts[1, i]
        p2 <- self$contacts[2, i]
        
        person1 = self$people[[p1]]
        person2 = self$people[[p2]]
        
        person1$check_infection(person2, t)
        person2$check_infection(person1, t)
      }
      return(self)
    },
    
    check_recoveries = \(t) {
      for (person in self$people) {
        person$check_recovery(t)
      }
      return(self)
    },
    
    check_immunities = \(t) {
      for (person in self$people) {
        person$check_immunity(t)
      }
      return(self)
    },
    
    count = \(t) {
      this_S = list()
      this_I = list()
      
      for (i in seq_along(self$people)) {
        person = self$people[[i]]
        if (person$S)
          this_S = c(this_S, i)
        if (person$I)
          this_I = c(this_I, i)
      }
      
      self$S[t] = self$S[t] + length(this_S)
      self$I[t] = self$I[t] + length(this_I)
      
      self$S_full <- c(self$S_full, list(this_S))
      self$I_full <- c(self$I_full, list(this_I))
      
      return(self)
    },
    
    run = \() {
      for (t in self$Ti) {
        self$check_infections(t) # Check which infectious occur
        self$check_recoveries(t) # Check which recoveries occur
        self$check_immunities(t) # Check which recoveries occur
        self$count(t)            # Store results
      }
      
      return(self)
    },
    
    plot = \() {
      pars = self$pars
      df <- data.frame(t = self$time,
                       S = self$S,
                       I = self$I) |>
        setDT() |>
        _[, melt(
          .SD,
          id.vars = 't',
          variable.name = 'state',
          value.name = 'number'
        )]
      
      df |>
        ggplot() +
        geom_xspline(aes(x = t, y = number, col = state)) +
        theme_light() +
        theme(
          axis.title = element_text(color = 'black'),
          axis.text = element_text(color = 'black'),
          axis.ticks = element_line(color = 'black'),
          legend.position = 'bottom'
        ) +
        scale_color_manual(
          values = c('S' = pars$colors$S,
                     'I' = pars$colors$I),
          labels = c('S' = 'Susceptible',
                     'I' = 'Infected')
        ) +
        ylab('Number of people') +
        xlab('Time') +
        labs(color = 'State:')
      
    }
  )
)



# Example usage
sim <- Sim$new()
sim$run()
sim$plot()
