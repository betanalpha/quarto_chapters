############################################################
# Setup
############################################################

par(family="serif", las=1, bty="l",
    cex.axis=1, cex.lab=1, cex.main=1,
    xaxs="i", yaxs="i", mar = c(5, 5, 3, 1))

library(rstan)
rstan_options(auto_write = TRUE)            # Cache compiled Stan programs
options(mc.cores = parallel::detectCores()) # Parallelize chains
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

util <- new.env()
source('mcmc_analysis_tools_rstan.R', local=util)
source('mcmc_visualization_tools.R', local=util)

############################################################
# Simulate Data
############################################################

# True Conversion Behavior
theta <- function(x, psi1, psi2) {
  psi1 * (-expm1(-psi2 * x))
}

psi1 <- 0.47
psi2 <- 5e-4

par(mfrow=c(1, 1), mar=c(5, 5, 3, 1))

xs <- seq(0, 5000, 10)
ys <- sapply(xs, function(x) theta(x, psi1, psi2))
plot(xs, ys, type="l", lwd=2, col=util$c_dark,
     xlim=c(0, 5000), xlab="USD",
     ylim=c(0, 1), ylab="Conversion Probability")

lambda <- 0.3
theta_VIP <- 0.95
  
ys <- sapply(xs, function(x) lambda * theta_VIP + 
                             (1 -lambda) * theta(x, psi1, psi2))
lines(xs, ys, type="l", lwd=2, col=util$c_light)


# Simulate Data
simu <- stan(file="stan_programs/simu_data.stan",
             algorithm="Fixed_param",
             warmup=0, iter=1, chains=1, 
             seed=4838282, refresh=0)

samples <- util$extract_expectands(simu)

N <- 1500
x <- sapply(1:N, function(n)
                 samples[[paste0('x[', n,  ']')]][1, 1])
y <- sapply(1:N, function(n)
                 samples[[paste0('y[', n,  ']')]][1, 1])

stan_rdump(c('N', 'x', 'y'), 'logs.data.R')

N_aux <- 200
y_aux <- sapply(1:N_aux, function(n)
                         samples[[paste0('y_aux[', n,  ']')]][1, 1])

stan_rdump(c('N_aux', 'y_aux'), 'aux_logs.data.R')
