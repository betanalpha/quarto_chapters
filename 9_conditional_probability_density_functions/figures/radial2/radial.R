######################################################################
# Configure Environment
######################################################################

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
# Marginal Probability Density Function 
######################################################################

alpha <- 0.8 * 5
beta <- 5

par(mfrow=c(1, 1), mar = c(5, 5, 2, 1))

rs <- seq(0, 3, 0.025)
ps <- dgamma(rs, alpha, beta)

plot(rs, ps)

# Save density values for TikZ visualization
cat(sprintf('%.3f/%.3f,', rs, ps))

######################################################################
# Conditional Probability Density Function
######################################################################

kappa <- 5

par(mfrow=c(4, 1), mar = c(5, 5, 2, 1))

thetas <- seq(0, 2 * pi, 0.01)

for (r in c(0.5, 1, 1.5, 2)) {
  mu <- pi * r
  ps <- exp( kappa * cos(thetas - mu) ) / (2 * pi * besselI(kappa, 0))
  
  plot(thetas, ps)
}

rs <- seq(0.2, 1.8, 0.4)
for (r in rs) {
  theta0 <- pi * r
  theta1 <- pi * (r - 1)
  theta2 <- pi * (r + 1)
  
  thetas <- seq(theta1, theta2, 2 * pi / 300)
  xs <- sapply(thetas, function(theta) r * cos(theta) )
  ys <- sapply(thetas, function(theta) r * sin(theta) )
  cpds <- sapply(thetas, function(theta) exp( kappa * cos(theta - pi * r) ) / (2 * pi * besselI(kappa, 0)))
  
  cat(sprintf("(%.3f, %.3f, %.3f)\n", xs, ys, cpds), 
      file=paste0('', as.integer(100 * r), ".dat"))
}

######################################################################
# Joint Probability Density Function
######################################################################

C <- alpha * log(beta) - lgamma(alpha)
D <- 1 / (2 * pi * besselI(kappa, 0))

cat(sprintf('\\a = %.8f;\n', alpha))
cat(sprintf('\\b = %.8f;\n', beta))
cat(sprintf('\\C = %.8f;\n', C))
cat(sprintf('\\k = %.8f;\n', kappa))
cat(sprintf('\\D = %.8f;\n', D))

N <- 500
xs <- seq(-2, 2, 4 / N)
ys <- seq(-2, 2, 4 / N)
zs <- matrix(0, nrow=N, ncol=N)

for (n1 in 1:N) {
  for (n2 in 1:N) {
    x <- xs[n1]
    y <- ys[n2]
    
    r <- sqrt(x**2 + y**2)
    theta <- atan2(y, x)
    
    zs[n1, n2] <- D * dgamma(r, alpha, beta) * exp( kappa * cos(theta - pi * r) )
  }
}


par(mfrow=c(1, 1), mar = c(5, 5, 2, 1))

image(xs, ys, zs, col=rev(cont_colors), axes=FALSE, ann=FALSE)