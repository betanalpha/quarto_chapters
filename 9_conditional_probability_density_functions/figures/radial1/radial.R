######################################################################
# Configure Environment
######################################################################

library(MSCquartets)
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


disc_colors <- c("#FFFFFF", "#DCBCBC", "#C79999", "#B97C7C",
                 "#A25050", "#8F2727", "#7C0000")
cont_colors <- colormap(colormap=disc_colors, nshades=100)

######################################################################
# Initial Probability Density Function
######################################################################

multinormal_pdf <- function(x, y, s, rho) {
  Q <- ( x**2 - 2 * rho * x * y + y**2 ) / ( s**2 * (1 - rho**2) )
  exp(-0.5 * Q) / (6.283185307179586 * sx * sy * sqrt(1 - rho**2) )
}

s <- 1.5
rho <- 0.75

N <- 200
xs <- seq(0, 5, 5 / (N - 1))
ys <- seq(0, 5, 5 / (N - 1))
zs <- matrix(0, nrow=N, ncol=N)

for (n in 1:N) {
  for (m in 1:N) {
    zs[n, m] <- multinormal_pdf(xs[n], ys[m], s, rho)
  }
}

image(xs, ys, zs, col=rev(cont_colors), xlab="x", ylab="y")

# Estimate normalization under truncation to the upper-right 
# quadrant or R^2.
delta <- 1e-3
xs <- seq(0, 5, delta)
ys <- seq(0, 5, delta)

norm <- 0
for (x in xs) {
  for (y in ys) {
    norm <- norm + multinormal_pdf(x, y, s, rho)
  }
}
norm <- delta**2 * norm

######################################################################
# Pushforward Probability Density Function 
######################################################################

# Generate samples from initial probability distribution function
# and pushforward along a radial function.
Sigma <- matrix(0, nrow=2, ncol=2)
Sigma[1, 1] <- s**2
Sigma[1, 2] <- rho * s**2
Sigma[2, 1] <- rho * s**2
Sigma[2, 2] <- s**2

L <- t(chol(Sigma))

S <- 100000
sxs <- rep(NA, S)
sys <- rep(NA, S)
srs <- rep(NA, S)

idx <- 0
while (idx < S) {
  z <- rnorm(2, 0, 1)
  p <- L %*% z
  
  if (p[1] >=0 & p[2] >= 0) {
    idx <- idx + 1
    sxs[idx] <- p[1]
    sys[idx] <- p[2]
    srs[idx] <- sqrt(p[1]**2 + p[2]**2)
  }
}

image(xs, ys, zs, col=rev(cont_colors), xlab="x", ylab="y")
points(sxs, sys, pch=16, cex=1, col="#0000FF03")

# Compare analytic pushforward probability density function to 
# histogram of pushforward samples
radial_pushforward_pdf <- function(r, s, rho) {
  C <- 1 / (4 * s**2 * sqrt(1 - rho**2)) * (1 / norm)
  alpha <- 0.5 * (1 / (1 - rho**2)) * (r / s)**2
  C * r * exp(-alpha) * ( 2 * besselI(alpha * rho, 0) + M0(alpha * rho) )
}

hist(srs, breaks=seq(0, 10, 0.2), 
     col=c_mid, border=c_dark_highlight, prob=TRUE, xlab="r")

rs <- seq(0, 10, 0.01)
pds <- sapply(rs, function(r) radial_pushforward_pdf(r, s, rho))
lines(rs, pds, col="blue")

# Save density values for TikZ visualization
rs <- seq(0, 10, 0.05)
pds <- sapply(rs, function(r) radial_pushforward_pdf(r, s, rho))
cat(sprintf("%.3f/%.3f, ", rs, pds), "\n")

######################################################################
# Conditional Probability Density Functions
######################################################################

conditional_pdf <- function(theta, r, s, rho) {
  Q <- ( 1 - rho * sin(2 * theta) ) * r**2 / ( s**2 * (1 - rho**2) )
  input_pd <- r * exp(-0.5 * Q) / (6.283185307179586 * sx * sy * sqrt(1 - rho**2) ) * (1 / norm)
  output_pd <- radial_pushforward_pdf(r, s, rho)
  input_pd / output_pd
}

par(mfrow=c(5, 5))

delta <- 7 / 25
rs <- seq(delta, 7, delta)
thetas <- seq(0, pi / 2, 0.1)

for (r in rs) {
  cpds <- sapply(thetas, function(theta) conditional_pdf(theta, r, s, rho) )
  
  plot(thetas, cpds, col=c_dark, type="l", 
       xlab="theta", ylim=c(0, 3), main=paste0("r = ", r))
}

# Save density values for TikZ visualization
rs <- seq(0.5, 6.5, 0.5)
for (r in rs) {
  thetas <- seq(0, pi / 2, (pi / 2) / 200)
  xs <- sapply(thetas, function(theta) r * cos(theta) )
  ys <- sapply(thetas, function(theta) r * sin(theta) )
  cpds <- sapply(thetas, function(theta) conditional_pdf(theta, r, s, rho) )

  xs <- c(r, xs, 0.0)
  ys <- c(0.0, ys, r)
  cpds <- c(0, cpds, 0)

  cat(sprintf("(%.3f, %.3f, %.3f)\n", xs, ys, cpds), 
      file=paste0('conditional/data/', as.integer(10 * r), ".dat"))
}
