

if (require(pacman))
  p_load(R6, ggplot2, data.table, purrr, ggalt)

# Define parameters
beta <- 2.5           # Infection rate
gamma <- 1.0          # Recovery rate
contact_rate <- 0.5   # Fraction of population each person connects to
I0 <- 5               # Number of people initially infected
N <- 100              # Total population size
maxtime <- 10         # How long to simulate for
npts <- 100           # Number of time points during the simulation
dt <- maxtime / npts  # Timestep length

Person <- R6Class(
  'Person',
  public = list(
    S = T,
    I = F,
    R = F,
    
    infect = \() {
      self$S = F
      self$I = T
      return(self)
    },
    
    recover = \() {
      self$I = F
      self$R = T
    },
    
    check_infection = \(other) {
      if (self$S &
          other$I & (runif(1) < beta / (N * contact_rate) * dt)) {
        self <- self$infect()
        return(self)
      }
    },
    
    check_recovery = \() {
      if (self$I & (runif(1) < (gamma * dt))) {
        self <- self$recover()
        return(self)
      }
    }
    
  )
)

Sim <- R6Class(
  'Sim',
  public = list(
    x = seq(0, npts - 1),
    S = rep(0, npts),
    I = rep(0, npts),
    R = rep(0, npts),
    time = NULL,
    people = list(),
    
    # Setting the initial conditions
    initialize = \() {
      self$people = purrr::map(1:N, ~ Person$new())
      self$people[1:I0] = map(1:I0, ~ self$people[[.x]]$infect())
      self$time = self$x * dt
    },
    
    count = \(self) {
      S = sum(unlist(map(self$people, \(x) x$S)))
      I = sum(unlist(map(self$people, \(x) x$I)))
      R = sum(unlist(map(self$people, \(x) x$R)))
      
      return(list(S = S, I = I, R = R))
    },
    
    check_infections = \() {
      contacts <- sample(1:N, size = floor(N * contact_rate))
      map(self$people, \(person1) {
        map(contacts, \(contact) {
          person2 <- self$people[[contact]]
          person1$check_infection(person2)
        })
      })
      return(self)
    },
    
    check_recoveries = function() {
      map(self$people, ~ .x$check_recovery())
      return(self)
    },
    
    run = function() {
      for (t in self$x) {
        self$check_infections()
        self$check_recoveries()
        
        # Update results
        counts <- self$count(self)
        self$S[[t + 1]] <- counts$S
        self$I[[t + 1]] <- counts$I
        self$R[[t + 1]] <- counts$R
      }
      cat("Run finished\n")
      return(self)
    },
    plot = \() {
      df <- data.frame(
        t = self$time,
        S = self$S,
        I = self$I,
        R = self$R
      ) |>
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
          values = c(
            'S' = 'blue',
            'I' = 'red',
            'R' = 'green'
          ),
          labels = c('S' = 'Susceptible',
                     'I' = 'Infectious',
                     'R' = 'Recovered')
        ) +
        ylab('Number of people') +
        xlab('Time') +
        labs(color = 'State:')
      
    }
  )
  
)
sim <- Sim$new()
sim$run()
sim$plot()


