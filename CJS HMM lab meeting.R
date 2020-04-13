###
# CJS data simulation from Kery and Shaub, Bayesian Populatin Analysis Using WinBUGS
###

#Phi(.)p(.)

# Define parameter values
n.occasions <- 6                   # Number of capture occasions
marked <-
  rep(50, n.occasions - 1)   # Annual number of newly marked individuals
phi <- rep(0.65, n.occasions - 1)
p <- rep(0.4, n.occasions - 1)

# Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol = n.occasions - 1, nrow = sum(marked))
P <- matrix(p, ncol = n.occasions - 1, nrow = sum(marked))

# Define function to simulate a capture-history (CH) matrix
simul.cjs <- function(PHI, P, marked) {
  n.occasions <- dim(PHI)[2] + 1
  CH <- matrix(0, ncol = n.occasions, nrow = sum(marked))
  # Define a vector with the occasion of marking
  mark.occ <- rep(1:length(marked), marked[1:length(marked)])
  # Fill the CH matrix
  for (i in 1:sum(marked)) {
    CH[i, mark.occ[i]] <- 1       # Write an 1 at the release occasion
    if (mark.occ[i] == n.occasions)
      next
    for (t in (mark.occ[i] + 1):n.occasions) {
      # Bernoulli trial: does individual survive occasion?
      sur <- rbinom(1, 1, PHI[i, t - 1])
      if (sur == 0)
        break		# If dead, move to next individual
      # Bernoulli trial: is individual recaptured?
      rp <- rbinom(1, 1, P[i, t - 1])
      if (rp == 1)
        CH[i, t] <- 1
    } #t
  } #i
  return(CH)
}

# Execute function
set.seed(2020)
CH <- simul.cjs(PHI, P, marked)

# Create vector with occasion of marking
get.first <- function(x)
  min(which(x != 0))
first <- apply(CH, 1, get.first)


##############################################################
#       Begin HMM model 
##############################################################

#function to return negative log probability of data | parameters (a.k.a negative log likelihood)
HMM_CJS <- function(params, CH, first) {
  phi <- plogis(params[1]) #take survival parameter from parameter vector and transform from logit to probabilty scale
  p <- plogis(params[2])   #take detection and transform to probability
  
  #create initial matrix
  delta <-
    matrix(c(1, 0), nrow = 1, ncol = 2) #all animals start alive when tagged
  
  #create transition probabilities matrix
  Gamma <- matrix(0, nrow = 2, ncol = 2)
  Gamma[1, 1] <- phi    #survive
  Gamma[1, 2] <- 1 - phi#die
  Gamma[2, 1] <- 0      #rise from the dead (funny that I'm writing this on Easter)
  Gamma[2, 2] <- 1      #remain dead
  
  #observation matrix
  obs <- matrix(0, nrow = 2, ncol = 2)
  obs[1, 1] <- (1 - p) # not seen given alive
  obs[1, 2] <- p       # seen given alive
  obs[2, 1] <- 1       # not seen given dead
  obs[2, 2] <- 0       # seen given dead
  
  ##
  # Likelihood
  ##
  
  #code from Pg. 39 of Zucchini et al 2016 generalized to multiple animals 
  
  #Initialize negative log likelihood of data for all animals at 0
  lscale <- 0
  
  for (i in 1:nrow(CH)) { #loop through animals
    
    alpha <- delta        #start with initial distribution (i.e. conditioning on alive at capture)
    for (j in (first[i] + 1):(n.occasions)) {     # loop through each time step
      
      alpha <- alpha %*% Gamma * obs[, (CH[i, j] + 1)] #probabilities through t-1 * transition probability matrix * observation (seen or not seen) probability given state probabilities
    }
    lscale <- lscale + log(sum(alpha)) # add log likelihood of for a given animal to the total log likelihood
  }
  return(-lscale) # return negative log likelihood
}

#########################################################################
#          Alternative model algorithem designed to prevent underflow
#########################################################################

HMM_CJS_2<- function(params, CH, first) {
  phi <- plogis(params[1]) #take survival parameter from parameter vector and transform from logit to probabilty scale
  p <- plogis(params[2])  #take detection and transform to probability
  
  #create initial matrix
  delta <-
    matrix(c(1, 0), nrow = 1, ncol = 2) #all animals start alive when tagged
  
  #create transition probabilities matrix
  Gamma <- matrix(0, nrow = 2, ncol = 2)
  Gamma[1, 1] <- phi    #survive
  Gamma[1, 2] <- 1 - phi  #die
  Gamma[2, 1] <-
    0      #rise from the dead (funny that I'm writing this on Easter)
  Gamma[2, 2] <- 1      #remain dead
  
  #observation matrix
  obs <- matrix(0, nrow = 2, ncol = 2)
  obs[1, 1] <- (1 - p)    # not seen given alive
  obs[1, 2] <- p        # seen given alive
  obs[2, 1] <- 1       # not seen given dead
  obs[2, 2] <- 0       # seen given dead
  
  ##
  # Likelihood
  ##
  
  #code from Pg. 49 of Zucchini et al 2016 generalized to multiple animals 
  
  #Initialize negative log likelihood at 0
  lscale <- 0
  
  for (i in 1:nrow(CH)) { #loop through each animal
    alpha <- delta  #start with initial distribution (i.e. conditioning on alive at capture)
    for (j in (first[i] + 1):(n.occasions)) {     # loop through each time step
      
      alpha <- alpha %*% Gamma * obs[, (CH[i, j] + 1)] #probs at t-1  (scaled) * transition probability matrix * observation (seen or not seen) probability given state probabilities
      lscale <- lscale + log(sum(alpha))  #add log probabilities to the total log likelihood
      alpha <- alpha / sum(alpha)         #rescale probabilties
    }
  }
  return(-lscale) #return negative log likelihood
}

#########################################################################
#          function to fit model and return delta method standard deviations and 95% confidence intercal
#########################################################################
HMM_opt<-function(func,         # function to optimize
                  inits,        # startihng parmeters
                  CH,           # function to optimize
                  first){       # time at first capture (from data)
  
out<-nlm(f=func,                # function to optimize
         p=inits,               # startihng parmeters
         CH=CH,                 # capture history (data)
         first=first,           # time at first capture (from data)
         hessian = TRUE)        # return hessian

params<-out$estimate            # paramaters in logit space

pars.real<-plogis(out$estimate) # parameters in probability scale
pars.real

l.var<-diag(solve(out$hessian)) # Variance on logit scale

var.params.real<-l.var * 
  ((exp(params))/((exp(params)+1)^2))^2 #delta method variance on probability scale


sd.params.real<-sqrt(var.params.real)   # Standard deviation on probabilty scale

#data table of results
return(data.frame(pars=pars.real,sd=sd.params.real,lower.95=plogis(qnorm(c(.025),params,sqrt(l.var))),upper.95=plogis(qnorm(c(.975),params,sqrt(l.var))),row.names =c(if(length(pars.real)>2) {paste(rep("Phi",5),1:5,sep=".")}else{"Phi"},"p"))) #c("Phi","p")))

}

#########################################################################
#          optimize model
#########################################################################
HMM_opt(func= HMM_CJS,          # function to optimize
        inits= c(0,0),          # startihng parmeters
        CH = CH,                # capture history (data)
        first=first)

#second algorithm
HMM_opt(func= HMM_CJS_2,          # function to optimize
        inits= c(0,0),          # startihng parmeters
        CH = CH,                # capture history (data)
        first=first)

##########################################################
##########################################################
##########################################################
#       Fixed time effects on survival
##########################################################

##simulate data

# Define parameter values
n.occasions <- 6                   # Number of capture occasions
marked <-
  rep(200, n.occasions - 1)        # Annual number of newly marked individuals 
phi <- seq(0.75,by=-.1, length=n.occasions - 1) # decreasing survival
p <- rep(0.4, n.occasions - 1)      # constant detection prob


# Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol = n.occasions - 1, nrow = sum(marked),byrow = T)
P <- matrix(p, ncol = n.occasions - 1, nrow = sum(marked))

# Execute function
set.seed(2020)
CH <- simul.cjs(PHI, P, marked)

# Create vector with occasion of marking
get.first <- function(x)
  min(which(x != 0))
first <- apply(CH, 1, get.first)

#########################################################################
#          Alternative model algorithem designed to prevent underflow
#          Update to allow fixed time effects on survival 
#########################################################################
HMM_CJS_2<- function(params, CH, first) {
  phi <- plogis(params[1:(n.occasions-1)]) #take survival parameterz from parameter vector and transform from logit to probabilty scale
  p <- plogis(params[n.occasions])  #take detection and transform to probability
  
  #create initial matrix
  delta <-
    matrix(c(1, 0), nrow = 1, ncol = 2) #all animals start alive when tagged
  
  #create transition probabilities matrix
  Gamma <- matrix(0, nrow = 2, ncol = 2)
  Gamma[2, 1] <-
    0      #rise from the dead (funny that I'm writing this on Easter)
  Gamma[2, 2] <- 1      #remain dead
  
  #observation matrix
  obs <- matrix(0, nrow = 2, ncol = 2)
  obs[1, 1] <- (1 - p)    # not seen given alive
  obs[1, 2] <- p        # seen given alive
  obs[2, 1] <- 1       # not seen given dead
  obs[2, 2] <- 0       # seen given dead
  
  ##
  # Likelihood
  ##
  
  #code from Pg. 49 of Zucchini et al 2016 generalized to multiple animals 
  
  #Initialize negative log likelihood at 0
  lscale <- 0
  
  for (i in 1:nrow(CH)) { #loop through each animal
    alpha <- delta  #start with initial distribution (i.e. conditioning on alive at capture)
    for (j in (first[i] + 1):(n.occasions)) {     # loop through each time step
      Gamma[1, 1] <- phi[j-1]    #survive
      Gamma[1, 2] <- 1 - phi[j-1]  #die
      
      alpha <- alpha %*% Gamma * obs[, (CH[i, j] + 1)] #probs at t-1  (scaled) * transition probability matrix * observation (seen or not seen) probability given state probabilities
      lscale <- lscale + log(sum(alpha))  #add log probabilities to the total log likelihood
      alpha <- alpha / sum(alpha)         #rescale probabilties
    }
  }
  return(-lscale) #return negative log likelihood
}


#second algorithm
HMM_opt(func= HMM_CJS_2,          # function to optimize
        inits= c(rep(0,n.occasions-1),0),          # startihng parmeters
        CH = CH,                # capture history (data)
        first=first)


##########################################################
##########################################################
##########################################################
#       Random time effects on survival
##########################################################

# 7.4.2. Random time effects
# Define parameter values
n.occasions <- 20                  # Number of capture occasions
marked <- rep(30, n.occasions-1)   # Annual number of newly marked individuals
mean.phi <- 0.65
var.phi <- 1                       # Temporal variance of survival
p <- rep(0.4, n.occasions-1)

# Determine annual survival probabilities
logit.phi <- rnorm(n.occasions-1, qlogis(mean.phi), var.phi^0.5)
phi <- plogis(logit.phi)

# Define matrices with survival and recapture probabilities
PHI <- matrix(phi, ncol = n.occasions-1, nrow = sum(marked), byrow = TRUE)
P <- matrix(p, ncol = n.occasions-1, nrow = sum(marked))

# Simulate capture-histories
CH <- simul.cjs(PHI, P, marked)

# Create vector with occasion of marking
get.first <- function(x) min(which(x!=0))
first <- apply(CH, 1, get.first)

##########################################################
library(TMB)

#compile and load TMB model
TMB::compile("CJS_HMM_lab_meeting.cpp")
dyn.load(dynlib("CJS_HMM_lab_meeting"))

# create data list
data_CJS<-list(CH=CH,first=first)

#initial parameters list
params<-list(logit_phi=0,logit_p=0,log_sigma=0,epsilon=rep(0,n.occasions-1))

#construct model
CJS_HMM<-MakeADFun(data_CJS,params,random = c("epsilon"),
                DLL="CJS_HMM_lab_meeting",silent=T,map=list())

# fit model (TMBhelper is available from Kasper Kristensen gitHub)
out<-TMBhelper::fit_tmb(CJS_HMM,CJS_HMM$fn, CJS_HMM$gr,CJS_HMM$par,getsd=TRUE, newtonsteps=1  )

df_out<-data.frame(MLE=out$SD$value,
SD=out$SD$sd)[c(2:1,22),]
rownames(df_out)<-c("phi","p","sigma")
df_out

