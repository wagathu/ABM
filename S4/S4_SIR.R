

beta <- 2.5  # Infection rate
gamma <- 1.0  # Recovery rate
contact_rate <- 0.5  # Fraction of population each person connects to
I0 <- 5  # Number of people initially infected
N <- 100  # Total population size
maxtime <- 10  # How long to simulate for
npts <- 100  # Number of time points during the simulation
dt <- maxtime / npts  # Timestep length




setClass('Person',
         slots = list(
           S = 'numeric',
           I = 'numeric',
           R = 'numeric'
         )
)

setMethod('initialize',
          'Person',
          \(.Object){
            .Object <-callNextMethod(.Object)
            .Object@S <- 1
            .Object@I <- 0
            .Object@R <- 0
            return(.Object)
          })

setGeneric('Infect',
           \(.Object) {
             standardGeneric('Infect')
           })

setMethod('Infect',
          'Person',
          \(.Object) {
            .Object@S <- 0
            .Object@I <- 1
            return(.Object)
          })

setGeneric('Recover',
           \(.Object) {
             standardGeneric('Recover')
           })

setMethod('Recover',
          'Person',
          \(.Object) {
            .Object@I <- 0
            .Object@R <- 1
            return(.Object)
          })

setGeneric('check_infection',
           \(.Object, other) {
             standardGeneric('check_infection')
           })

setMethod('check_infection',
          signature = 'Person',
          \(.Object, other) {
            if (.Object@S & other@I & runif(1) < beta / (N * contact_rate) * dt) {
              .Object <- Infect(.Object)
            }
         return(.Object)
          }
)

setGeneric('check_recovery',
           \(.Object) {
             standardGeneric('check_recovery')
           })

setMethod('check_recovery',
          'Person',
          \(.Object) {
            if (.Object@I & runif(1) < (gamma * dt)) {
              .Object <- Recover(.Object)
            }
            return(.Object)
          })


# Define the Simulation class
setClass(
  'Sim',
   slots = list(
     x = 'numeric',
     S = 'numeric',
     I = 'numeric',
     R = 'numeric',
     time = 'numeric',
     people = 'list'
   )
)


setMethod('initialize',
          'Sim',
          \(.Object) {
            .Object <- callNextMethod(.Object)
              .Object@x = seq(0, npts - 1)
              .Object@S = rep(0, npts)
              .Object@I = rep(0, npts)
              .Object@R = rep(0, npts)
              .Object@time = NA_real_
              .Object@people = list()
 
          
            for (i in 1:N) {
              .Object@people[[i]] <- new('Person')
            }
            
            for (i in 1:I0) {
              .Object@people[[i]] <- Infect(.Object@people[[i]] )
            }
              .Object@time <- .Object@x * dt
            return(.Object) 
          })

setGeneric(
  'Count',
  \(.Object) {
  standardGeneric('Count')
  }
)

setMethod(
  'Count',
  'Sim',
  \(.Object){
    S = 0
    I = 0
    R = 0
    
    for (person in .Object@people) {
      S = S + person@S
      I = I + person@I
      R = R + person@R
    }
    return(list(S = S, I = I, R = R))
    
  }
)

setGeneric(
  'check_infections',
  \(.Object) {
    standardGeneric('check_infections')
  }
)

setMethod(  
  'check_infections',
  'Sim',
  \(.Object){
    for (person1 in .Object@people) {
      contacts <- sample(1:N, size = floor(N * contact_rate))
      for (contact in contacts) {
        person2 <- .Object@people[[contact]]
        check_infection(person1, person2) # Here is where the error is
      }
     
    }
  return(.Object)
  }
)

setGeneric(
  'check_recoveries',
  \(.Object){
    standardGeneric('check_recoveries')
  }
  )

setMethod(
  'check_recoveries',
  'Sim',
  \(.Object){
    for (person in .Object@people) {
      check_recovery(person)
      }
    return(.Object)
  }
  )

setGeneric(
  'Run',
  \(.Object){
    standardGeneric('Run')
  }
)

setMethod(
  'Run',
  'Sim',
  \(.Object){
    for (t in .Object@x) {
      .Object <- check_infections(.Object)
      .Object <- check_recoveries(.Object)
      
      # Update results
      counts <- Count(.Object)
      .Object@S[t + 1] <- counts$S
      .Object@I[t + 1] <- counts$I
      .Object@R[t + 1] <- counts$R
    }
    
    cat("Run finished\n")
    return(.Object)
  }
)

sim.init <- new('Sim')
res <- Run(sim.init)

jj <- new('Sim')
.Object <- jj




