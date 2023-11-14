############################################################
# Initial setup
############################################################

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

util <- new.env()
source('stan_utility_rstan.R', local=util)

library(colormap)

c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

c_light_teal <- c("#6B8E8E")
c_mid_teal <- c("#487575")
c_dark_teal <- c("#1D4F4F")

par(family="serif", las=1, bty="l", cex.axis=1, cex.lab=1, cex.main=1,
    xaxs="i", yaxs="i", mar = c(5, 5, 3, 1))

############################################################
# Visualize Prior Model
############################################################

disc_colors <- c("#FFFFFF", "#DCBCBC", "#C79999", "#B97C7C",
                 "#A25050", "#8F2727", "#7C0000")
cont_colors <- colormap(colormap=disc_colors, nshades=100)

# Implementation of log(x^alpha) with safe handling of alpha = 0
lmult <- function(x, a) {
  if (a == 0) {
    (0)
  } else {
    (a * log(x))
  }
}

# Dirichlet probability density functions
dirichlet_pdf <- function(p1, p2, p3, a1, a2, a3) {
  exp( lmult(p1, a1 - 1) + lmult(p2, a2 - 1) + lmult(p3, a3 - 1) 
       + lgamma(a1 + a2 + a3) - lgamma(a1) - lgamma(a2) - lgamma(a3) )
}

# Plot Dirichlet probability density function across a simplex
plot_dirichlet <- function(a1, a2, a3) {
  N <- 200
  C <- 1 / sqrt(3)
  xs <- seq(-C - 0.2, C + 0.2, (2 * C + 0.4) / (N - 1))
  ys <- seq(0 - 0.1, 1 + 0.1, 1.2 / (N - 1))
  zs <- matrix(0, nrow=N, ncol=N)

  for (n in 1:N) {
    for (m in 1:N) {
      p3 <- ys[m]
      p2 <- (1 + xs[n] / C - ys[m]) / 2
      p1 <- (1 - xs[n] / C - ys[m]) / 2
    
      if (p1 >= 0 & p2 >= 0 & p3 >= 0 & p1 + p2 + p3 <= 1) {
        zs[n, m] <- dirichlet_pdf(p1, p2, p3, a1, a2, a3)
      } else {
        zs[n, m] <- NA
      }
    }
  }

  par(mar = c(0, 0, 0, 0)) 
  image(xs, ys, zs, col=rev(cont_colors), axes=FALSE, ann=FALSE)
  
  lines( c(-C, 0), c(0, 1), lwd=3)
  lines( c(+C, 0), c(0, 1), lwd=3)
  lines( c(-C, +C), c(0, 0), lwd=3)

  text_delta <- 0.05
  text( 0, 1 + text_delta, "(0, 0, 1)", cex=2)
  text(-C - text_delta, -text_delta, "(1, 0, 0)", cex=2)
  text(+C + text_delta, -text_delta, "(0, 1, 0)", cex=2)
  
  tick_delta <- 0.025
  lines( c(0, 0), c(0, tick_delta), lwd=3)
  text(0, 0 - text_delta, "(1/2, 1/2, 0)", cex=2)
  
  lines( c(+C * 0.5, +C * 0.5 - tick_delta * 0.5 * sqrt(3)), 
         c(0.5, 0.5 - tick_delta * 0.5), lwd=3)
  text(C * 0.5 + text_delta * 0.5 * sqrt(3) + 2.5 * text_delta, 
       0.5 + text_delta * 0.5, "(0, 1/2, 1/2)", cex=2)

  lines( c(-C * 0.5, -C * 0.5 + tick_delta * 0.5 * sqrt(3)), 
         c(0.5, 0.5 - tick_delta * 0.5), lwd=3)
  text(-C * 0.5 - text_delta * 0.5 * sqrt(3) - 2.5 * text_delta, 
       0.5 + text_delta * 0.5, "(1/2, 0, 1/2)", cex=2)
  
  points(0, 1/3, col="white", pch=16, cex=1.5)
  points(0, 1/3, col="black", pch=16, cex=1)
  text(0, 1/3 - 1.5 * text_delta, "(1/3, 1/3, 1/3)", cex=2)
}

# Uniform probability density function
plot_dirichlet(1, 1, 1) 

# Probability density function concentrates around equal probability 
# configuration (1/3, 1/3, 1/3)
plot_dirichlet(5, 5, 5)

# Probability density function peaks around (1, 0, 0) but tails are 
# multimodal, concentrating on configurations with either p2 >> p3 or 
# p3 >> p2.
plot_dirichlet(1, 0.8, 0.8)

# Probability density function peaks around (1, 0, 0), but not tails 
# are symmetric for p2 and p3.
plot_dirichlet(2, 1, 1)
plot_dirichlet(3, 1, 1)
plot_dirichlet(4, 1, 1)
plot_dirichlet(5, 1, 1)

############################################################
# Simulate Data
############################################################

simu <- stan(file="stan_programs/simu.stan",
             iter=1, warmup=0, chains=1,
             seed=4838282, algorithm="Fixed_param")

# Separate into observations with unknown correct answers
# and a smaller calibration data set where the the correct
# answers are known.

data <- list('K' = 3, 'N_reviewers' = 5, 'N_assessments' = 100,
             'y' = extract(simu)$y[1,1:100,])

true_zs <- extract(simu)$z[1,1:100]

calib_data <- list('N_calib_assessments' = 10,
                   'y_calib' = extract(simu)$y[1,101:110,],
                   'z_calib' = extract(simu)$z[1,101:110])

############################################################
# Display Data
############################################################

# We can visualize the data for each assessment by 
# histogramming the reviewer responses.

plot_line_hist <- function(s, bin_min, bin_max, delta, xlab, ylim=NA, main="") {
  bins <- seq(bin_min, bin_max, delta)
  B <- length(bins) - 1
  idx <- rep(1:B, each=2)
  x <- sapply(1:length(idx),
              function(b) if(b %% 2 == 1) bins[idx[b]] else bins[idx[b] + 1])
  x <- c(bin_min - 10, x, bin_max + 10)
  
  counts <- hist(s, breaks=bins, plot=FALSE)$counts
  y <- counts[idx]
  y <- c(0, y, 0)
  
  ymax <- max(y) + 1
  
  if (any(is.na(ylim))) {
    ylim <- c(0, ymax)
  }
  
  plot(x, y, type="l", main=main, col="black", lwd=2,
       xlab=xlab, xlim=c(bin_min, bin_max),
       ylab="Counts", ylim=ylim)
}

# We'll just look at the first 25 assessments.
par(mfrow=c(5, 5), mar = c(5, 4, 2, 1))

for (n in 1:25) {
  plot_line_hist(data$y[n,], -0.5, data$K + 1.5, 1, "Answer", 
                 c(0, 5.5), paste("Assesment", n))
  abline(v=true_zs[n], lwd=3, col=c_light)
}

# There is usually some consensus amongst the reviewers 
# but in some cases the agreement is much weaker leaving 
# us to question what the true answer might be.

############################################################
# Model 1
############################################################

# We'll start with flat priors for all of the simplices involved.

# Quantify posterior distribution
fit <- stan(file="stan_programs/fit1.stan", 
            data=data, seed=8438338,
            warmup=1000, iter=2024, refresh=0)

# Diagnostics reveal serious R-hat issues for the reviewer 
# fidelities.  Because there are not accompanying empirical
# effective sample size warnings this suggests multimodality
# in the posterior distribution.

diagnostics <- util$extract_hmc_diagnostics(fit) 
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectands(fit)
names <- grep('pred', names(samples), value=TRUE, invert=TRUE)
base_samples <- samples[names]
util$check_all_expectand_diagnostics(base_samples)

# Before exploring the multimodality, however, let's examine
# the posterior retrodicive performance.  We'll use the histogram 
# of reviewer responses for each assessment as our summary
# statistic.

hist_retro <- function(obs, samples, pred_names, bin_min, bin_max, delta, xlab="", ylim=NA, title="") {
  if (is.na(bin_min)) bin_min <- min(pred)
  if (is.na(bin_max)) bin_max <- max(pred)
  breaks <- seq(bin_min, bin_max, delta)
  B <- length(breaks) - 1
  
  idx <- rep(1:B, each=2)
  xs <- sapply(1:length(idx), 
               function(b) if(b %% 2 == 0) breaks[idx[b] + 1] else breaks[idx[b]] )
  
  obs_counts <- hist(obs, breaks=breaks, plot=FALSE)$counts
  pad_obs_counts <- do.call(cbind, lapply(idx, function(n) obs_counts[n]))
  
  pred <- sapply(pred_names, function(name) c(t(samples[[name]]), recursive=TRUE))
  N <- dim(pred)[1]
  pred_counts <- sapply(1:N, 
                        function(n) hist(pred[n,], breaks=breaks, plot=FALSE)$counts)
  
  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  cred <- sapply(1:B, function(b) quantile(pred_counts[b,], probs=probs))
  pad_cred <- do.call(cbind, lapply(idx, function(n) cred[1:9, n]))
  
  if (any(is.na(ylim))) {
    ylim <- c(0, max(c(obs_counts, cred[9,])))
  }
  
  plot(1, type="n", main=title,
       xlim=c(bin_min, bin_max), xlab=xlab,
       ylim=ylim, ylab="Counts")
  
  polygon(c(xs, rev(xs)), c(pad_cred[1,], rev(pad_cred[9,])),
          col = c_light, border = NA)
  polygon(c(xs, rev(xs)), c(pad_cred[2,], rev(pad_cred[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(xs, rev(xs)), c(pad_cred[3,], rev(pad_cred[7,])),
          col = c_mid, border = NA)
  polygon(c(xs, rev(xs)), c(pad_cred[4,], rev(pad_cred[6,])),
          col = c_mid_highlight, border = NA)
  lines(xs, pad_cred[5,], col=c_dark, lwd=2)
  
  lines(xs, pad_obs_counts, col="white", lty=1, lw=2.5)
  lines(xs, pad_obs_counts, col="black", lty=1, lw=2)
}

# Again we'll consider only the first 25 assessments to simplify
# the presentation.

# Overall there doesn't seem to be any retrodictive tension, but 
# only because the posterior retrodictive distribution is so
# diffuse.

par(mfrow=c(5, 5), mar = c(5, 4, 2, 1))

for (n in 1:25) {
  pred_names <- grep(paste0('y_pred\\[', n, ','), names(samples), value=TRUE)
  hist_retro(data$y[n,], samples, pred_names, -0.5, data$K + 1.5, 1, 
             "Answer", c(0, 5), paste("Assesment", n))
}

# Let's see if we can find the multimodality hinted at by the
# diagnostic warnings.  We'll focus on the first reviewer.

# Very strong multimodalities!

par(mfrow=c(3, 3), mar = c(5, 4, 2, 1))

r <- 1
for (k_answer in 1:3) {
  for (k_true in 1:3) {
    name <- paste0('fidelity[', r, ',', k_true, ',', k_answer, ']')
    util$plot_expectand_pushforward(samples[[name]], 25, name, flim=c(0, 1))
  }
}

# We have to be careful, however, interpreting the number and relative sizes of 
# each of these peaks.  If there are many peaks then the four Markov chains we 
# ran may not have encountered them.  Moreover if each Markov chain has been
# captured by a particular mode then we won't have any meaningful information 
# about the relative posterior probabilities of those modes.

r <- 1
k_true <- 1
name1 <- paste0('fidelity[', r, ',', k_true, ',', 1, ']')
name2 <- paste0('fidelity[', r, ',', k_true, ',', 2, ']')
name3 <- paste0('fidelity[', r, ',', k_true, ',', 3, ']')

util$plot_pairs_by_chain(samples[[name1]], name1, samples[[name3]], name3)

# Our visualizations can be made much more efficient if we take advantage of
# the natural geometry of the 3-dimensoional simplex, allowing us to plot
# all three probabilities at once.

plot_simplex_samples_by_chain <- function(samples, names, title="") {
  vals <- lapply(names, function(name) samples[[name]])
  C <- dim(vals[[1]])[1]
  N <- dim(vals[[1]])[2]
  
  D <- 1 / sqrt(3)
  xlim <- c(-D - 0.2, +D+ 0.2)
  ylim <- c(-0.1, 1.1)
  
  par(mar = c(0, 0, 5, 0)) 
  plot(1, type="n", main=title, 
       xlim=xlim, ylim=ylim,
       axes=FALSE)
  
  lines( c(-D, 0),  c(0, 1), lwd=3)
  lines( c(+D, 0),  c(0, 1), lwd=3)
  lines( c(-D, +D), c(0, 0), lwd=3)
  
  text_delta <- 0.05
  text( 0, 1 + text_delta, "(0, 0, 1)", cex=1)
  text(-D - text_delta, -text_delta, "(1, 0, 0)", cex=1)
  text(+D + text_delta, -text_delta, "(0, 1, 0)", cex=1)
  
  nom_colors <- c("#DCBCBC", "#C79999", "#B97C7C",
                  "#A25050", "#8F2727", "#7C0000")
  chain_colors <- colormap(colormap=nom_colors, nshades=C)
  
  for (c in 1:C) {
    points(D * (vals[[2]][c,] - vals[[1]][c,]), vals[[3]][c,], 
           pch=16, cex=0.75, col=chain_colors[c])
  }
  
  points(0, 1/3, col="white", pch=16, cex=1.5)
  points(0, 1/3, col="black", pch=16, cex=1)
  text(0, 1/3 - 1.5 * text_delta, "(1/3, 1/3, 1/3)", cex=1)
}

r <- 1
par(mfrow=c(3, 1))

for (k_true in 1:3) {
  names <- sapply(1:data$K, function(k) paste0('fidelity[', r, ',', k_true, ',', k, ']'))
  plot_simplex_samples_by_chain(samples, names, paste0("k_true = ", k_true))
}

# Now the underlying problem is coming into focus.  The data are consistent with
# all of the reviewers being accurate, providing the correct answer more often 
# than not.  They are also consistent, however, with the reviewers all being 
# biased in the same way.  In other words there is a fundamental non-identifiability
# in the rows of the fidelity matrix for each reviewer.

# Concensus doesn't necessarily mean correctness!

############################################################
# Model 2
############################################################

# There are a few ways to break this degeneracy.

# For example we could incorporate calibration trials, where we give reviewers 
# assessments with knownanswers.  Alternatively if the true answer has other
# observational consequences then we could integrate that second data generating 
# process into a single joint model.

# Here we will use a more informative prior model that breaks the degeneracy 
# by assuming that reviewers should be reasonably accurate.

par(mfrow=c(3, 1))

for (k_true in 1:3) {
  plot_dirichlet(4, 1, 1)
  plot_dirichlet(1, 4, 1)
  plot_dirichlet(1, 1, 4)
}

# Initialize from the prior model to prevent the Markov chains from 
# being captured by extreme modes.

dirichlet_rng <- function(a1, a2, a3) {
  gammas <- c(rgamma(1, a1, 1), rgamma(1, a2, 1), rgamma(1, a3, 1))
  gammas / sum(gammas)
}

set.seed(48383499)

inits <- list()

for (c in 1:4) {
  chain_inits <- list()
  
  for (n in 1:data$N_assessments) {
    name <- paste0('p[', n, ']')
    chain_inits[[name]] <- dirichlet_rng(1, 1, 1)
  }
  
  for (r in 1:data$N_reviewers) {
    for (k in 1:data$K) {
      name <- paste0('fidelity[', r, ',', k, ']')
      if (k == 1)
        chain_inits[[name]] <- dirichlet_rng(4, 1, 1)
      else if (k == 2)
        chain_inits[[name]] <- dirichlet_rng(1, 4, 1)
      else if (k == 3)
        chain_inits[[name]] <- dirichlet_rng(1, 1, 4)
    }
  }
  
  inits[[c]] <- chain_inits
}


plot_simplex_inits <- function(inits, name, title="") {
  C <- length(inits)
  
  D <- 1 / sqrt(3)
  xlim <- c(-D - 0.2, +D+ 0.2)
  ylim <- c(-0.1, 1.1)
  
  par(mar = c(0, 0, 5, 0)) 
  plot(1, type="n", main=title, 
       xlim=xlim, ylim=ylim,
       axes=FALSE)
  
  lines( c(-D, 0),  c(0, 1), lwd=3)
  lines( c(+D, 0),  c(0, 1), lwd=3)
  lines( c(-D, +D), c(0, 0), lwd=3)
  
  text_delta <- 0.05
  text( 0, 1 + text_delta, "(0, 0, 1)", cex=1)
  text(-D - text_delta, -text_delta, "(1, 0, 0)", cex=1)
  text(+D + text_delta, -text_delta, "(0, 1, 0)", cex=1)
  
  nom_colors <- c("#DCBCBC", "#C79999", "#B97C7C",
                  "#A25050", "#8F2727", "#7C0000")
  chain_colors <- colormap(colormap=nom_colors, nshades=C)
  
  for (c in 1:C) {
    val <- inits[[c]][[name]]
    points(D * (val[2] - val[1]), val[3], 
           pch=16, cex=0.75, col=chain_colors[c])
  }
  
  points(0, 1/3, col="white", pch=16, cex=1.5)
  points(0, 1/3, col="black", pch=16, cex=1)
  text(0, 1/3 - 1.5 * text_delta, "(1/3, 1/3, 1/3)", cex=1)
}

r <- 1
par(mfrow=c(3, 1))

for (k_true in 1:3) {
  name <- paste0('fidelity[', r, ',', k_true, ']')
  plot_simplex_inits(inits, name, paste0("k_true = ", k_true))
}

# Quantify posterior distribution
fit <- stan(file="stan_programs/fit2.stan", 
            data=data, seed=8438338,
            warmup=1000, iter=2024, refresh=0, 
            init=inits)

# Check diagnostics
diagnostics <- util$extract_hmc_diagnostics(fit) 
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectands(fit)
names <- grep('pred', names(samples), value=TRUE, invert=TRUE)
base_samples <- samples[names]
util$check_all_expectand_diagnostics(base_samples)

# Retrodictive checks
par(mfrow=c(5, 5), mar = c(5, 4, 2, 1))

for (n in 1:25) {
  pred_names <- grep(paste0('y_pred\\[', n, ','), names(samples), value=TRUE)
  hist_retro(data$y[n,], samples, pred_names, -0.5, data$K + 1.5, 1,
             "Answer", c(0, 5), paste("Assesment", n))
}

# Posterior inferences
par(mfrow=c(1, 1))

for (k_true in c(2)) {
  names <- sapply(1:data$K, function(k) paste0('fidelity[', r, ',', k_true, ',', k, ']'))
  plot_simplex_samples_by_chain(samples, names, paste0("k_true = ", k_true))
}

r <- 1
par(mfrow=c(3, 1))

for (k_true in 1:3) {
  names <- sapply(1:data$K, function(k) paste0('fidelity[', r, ',', k_true, ',', k, ']'))
  plot_simplex_samples_by_chain(samples, names, paste0("k_true = ", k_true))
}

par(mfrow=c(3, 3))
for (n in 1:9) {
  names <- sapply(1:data$K, function(k) paste0('p[', n, ',', k, ']'))
  plot_simplex_samples_by_chain(samples, names, paste0("Assessment = ", n))
}


plot_simplex_samples <- function(samples, names, title="") {
  vals <- sapply(names, function(name) c(t(samples[[name]]), recursive=TRUE))
  N <- dim(vals)[1]
  
  C <- 1 / sqrt(3)
  xlim <- c(-C - 0.2, +C + 0.2)
  ylim <- c(-0.1, 1.1)
  
  par(mar = c(0, 0, 5, 0)) 
  plot(1, type="n", main=title, 
       xlim=xlim, ylim=ylim,
       axes=FALSE)
  
  lines( c(-C, 0), c(0, 1), lwd=3)
  lines( c(+C, 0), c(0, 1), lwd=3)
  lines( c(-C, +C), c(0, 0), lwd=3)
  
  text_delta <- 0.05
  text( 0, 1 + text_delta, "(0, 0, 1)", cex=1)
  text(-C - text_delta, -text_delta, "(1, 0, 0)", cex=1)
  text(+C + text_delta, -text_delta, "(0, 1, 0)", cex=1)
  
  points(C * (vals[,2] - vals[,1]), vals[,3], pch=16, cex=0.75, col=c_dark)
  
  points(0, 1/3, col="white", pch=16, cex=1.5)
  points(0, 1/3, col="black", pch=16, cex=1)
  text(0, 1/3 - 1.5 * text_delta, "(1/3, 1/3, 1/3)", cex=1)
}

r <- 2
par(mfrow=c(1, 3))
    
for (k_true in 1:3) {
  names <- sapply(1:data$K, function(k) paste0('fidelity[', r, ',', k_true, ',', k, ']'))
  plot_simplex_samples(samples, names, paste0("k_true = ", k_true))
}



plot_disc_marginal_quantiles <- function(samples, names, x_name="", title="") {
  N <- length(names)
  idx <- rep(1:N, each=2)
  xs <- sapply(1:length(idx), function(k) if(k %% 2 == 0) idx[k] + 0.5 else idx[k] - 0.5)
  
  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  cred <- sapply(1:N, function(n) quantile(c(t(samples[[names[n]]]), recursive=TRUE), probs=probs))
  pad_cred <- do.call(cbind, lapply(idx, function(n) cred[1:9,n]))
  
  ylims <- c(min(cred[1,]), max(cred[9,]))
  
  plot(1, type="n", main=title,
       xlim=c(0.5, N + 0.5), xlab=x_name,
       ylim=ylims, ylab="Marginal Posterior Quantiles")
  
  polygon(c(xs, rev(xs)), c(pad_cred[1,], rev(pad_cred[9,])),
          col = c_light, border = NA)
  polygon(c(xs, rev(xs)), c(pad_cred[2,], rev(pad_cred[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(xs, rev(xs)), c(pad_cred[3,], rev(pad_cred[7,])),
          col = c_mid, border = NA)
  polygon(c(xs, rev(xs)), c(pad_cred[4,], rev(pad_cred[6,])),
          col = c_mid_highlight, border = NA)
  for (n in 1:N) {
    lines(xs[(2 * n - 1):(2 * n)], pad_cred[5,(2 * n - 1):(2 * n)], col=c_dark, lwd=2)
  }
}

par(mfrow=c(5, 5), mar = c(5, 4, 2, 1))
for (n in 1:25) {
  names <- sapply(1:data$K, function(k) paste0('p[', n, ',', k, ']'))
  plot_disc_marginal_quantiles(samples, names, 
                               x_name="p_correct_answer", 
                               title=paste("Assesment", n))
  abline(v=true_zs[n], lwd=3, col="white")
  abline(v=true_zs[n], lwd=2, col="black")
}


plot_disc_realizations <- function(samples, names, N, xlab="", ylab="", ylim=NA, title="") {
  vals <- sapply(names, function(name) c(t(samples[[name]]), recursive=TRUE))
  
  I <- dim(vals)[1]
  J <- min(N, I)
  
  plot_idx <- sapply(1:J, function(j) (I %/% J) * (j - 1) + 1)
  
  nom_colors <- c("#DCBCBC", "#C79999", "#B97C7C",
                  "#A25050", "#8F2727", "#7C0000")
  line_colors <- colormap(colormap=nom_colors, nshades=J)
  
  K <- dim(vals)[2]
  
  idxs <- rep(1:K, each=2)
  xs <- sapply(1:length(idxs),
               function(b) if(b %% 2 == 1) idxs[b] - 0.5 else idxs[b] + 0.5)
  
  xlim <- c(0.5, K + 0.5)
  if (any(is.na(ylim))) {
    ylim <- c(0, max(vals))
  }
  
  plot(1, type="n", main=title,
       xlab=xlab, xlim=xlim,
       ylab=ylab, ylim=ylim)
  for (j in 1:J) {
    ys <- sapply(idxs, function(idx) vals[plot_idx[j], idx])
    lines(xs, ys, col=line_colors[j], lwd=3)
  }
}

par(mfrow=c(5, 5), mar = c(5, 4, 2, 1))
for (n in 1:25) {
  names <- sapply(1:data$K, function(k) paste0('p[', n, ',', k, ']'))
  plot_disc_realizations(samples, names, 5,
                         "p_answer", "p_correct_answer", c(0, 1),
                         title=paste("Assesment", n))
}

############################################################
# Model 3
############################################################

# The stronger prior model alone didn't seem to be sufficient
# to avoid the "off-diagonal" modes.  If we could elicit more
# precise domain expertise then we could push the prior model 
# even further, but that's not always possible.

# Here we'll introduce a small calibration data set.

data <- c(data, calib_data)

# Quantify posterior distribution
fit <- stan(file="stan_programs/fit3.stan", 
            data=data, seed=8438338,
            warmup=1000, iter=2024, refresh=0, 
            init=inits)

# Those frustrating Rhat warnings are gone.  Have we
# actually resolved the multimodality?

diagnostics <- util$extract_hmc_diagnostics(fit) 
util$check_all_hmc_diagnostics(diagnostics)

samples <- util$extract_expectands(fit)
names <- grep('pred', names(samples), value=TRUE, invert=TRUE)
base_samples <- samples[names]
util$check_all_expectand_diagnostics(base_samples)

# Retrodictive checks
par(mfrow=c(5, 5), mar = c(5, 4, 2, 1))

for (n in 1:25) {
  pred_names <- grep(paste0('y_pred\\[', n, ','), names(samples), value=TRUE)
  hist_retro(data$y[n,], samples, pred_names, -0.5, data$K + 1.5, 1,
             "Answer", c(0, 5), paste("Assesment", n))
}

# Posterior inferences
par(mfrow=c(1, 1))

# The right corners but still strong uncertainties.
r <- 1
par(mfrow=c(3, 1))

for (k_true in 1:3) {
  names <- sapply(1:data$K, function(k) paste0('fidelity[', r, ',', k_true, ',', k, ']'))
  plot_simplex_samples_by_chain(samples, names, paste0("k_true = ", k_true))
}

# Which results in very weakly informative posterior inferences
# for the correct answer simplices.
par(mfrow=c(1, 1))
for (n in c(1)) {
  names <- sapply(1:data$K, function(k) paste0('p[', n, ',', k, ']'))
  plot_simplex_samples_by_chain(samples, names, paste0("Assessment = ", n))
}


plot_disc_marginal_quantiles <- function(samples, names, x_name="", title="") {
  N <- length(names)
  idx <- rep(1:N, each=2)
  xs <- sapply(1:length(idx), function(k) if(k %% 2 == 0) idx[k] + 0.5 else idx[k] - 0.5)
  
  probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  cred <- sapply(1:N, function(n) quantile(c(t(samples[[names[n]]]), recursive=TRUE), probs=probs))
  pad_cred <- do.call(cbind, lapply(idx, function(n) cred[1:9,n]))
  
  ylims <- c(min(cred[1,]), max(cred[9,]))
  
  plot(1, type="n", main=title,
       xlim=c(0.5, N + 0.5), xlab=x_name,
       ylim=ylims, ylab="Marginal Posterior Quantiles")
  
  polygon(c(xs, rev(xs)), c(pad_cred[1,], rev(pad_cred[9,])),
          col = c_light, border = NA)
  polygon(c(xs, rev(xs)), c(pad_cred[2,], rev(pad_cred[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(xs, rev(xs)), c(pad_cred[3,], rev(pad_cred[7,])),
          col = c_mid, border = NA)
  polygon(c(xs, rev(xs)), c(pad_cred[4,], rev(pad_cred[6,])),
          col = c_mid_highlight, border = NA)
  for (n in 1:N) {
    lines(xs[(2 * n - 1):(2 * n)], pad_cred[5,(2 * n - 1):(2 * n)], col=c_dark, lwd=2)
  }
}

# That said the simplices do concentrate towards the correct 
# answers more than the incorrect answers.

par(mfrow=c(5, 5), mar = c(5, 4, 2, 1))
for (n in 1:25) {
  names <- sapply(1:data$K, function(k) paste0('p[', n, ',', k, ']'))
  plot_disc_marginal_quantiles(samples, names, 
                               x_name="p_correct_answer", 
                               title=paste("Assesment", n))
  abline(v=true_zs[n], lwd=3, col="white")
  abline(v=true_zs[n], lwd=2, col="black")
}


# Unused for the moment.

plot_disc_realizations <- function(samples, names, N, xlab="", ylab="", ylim=NA, title="") {
  vals <- sapply(names, function(name) c(t(samples[[name]]), recursive=TRUE))
  
  I <- dim(vals)[1]
  J <- min(N, I)
  
  plot_idx <- sapply(1:J, function(j) (I %/% J) * (j - 1) + 1)
  
  nom_colors <- c("#DCBCBC", "#C79999", "#B97C7C",
                  "#A25050", "#8F2727", "#7C0000")
  line_colors <- colormap(colormap=nom_colors, nshades=J)
  
  K <- dim(vals)[2]
  
  idxs <- rep(1:K, each=2)
  xs <- sapply(1:length(idxs),
               function(b) if(b %% 2 == 1) idxs[b] - 0.5 else idxs[b] + 0.5)
  
  xlim <- c(0.5, K + 0.5)
  if (any(is.na(ylim))) {
    ylim <- c(0, max(vals))
  }
  
  plot(1, type="n", main=title,
       xlab=xlab, xlim=xlim,
       ylab=ylab, ylim=ylim)
  for (j in 1:J) {
    ys <- sapply(idxs, function(idx) vals[plot_idx[j], idx])
    lines(xs, ys, col=line_colors[j], lwd=3)
  }
}