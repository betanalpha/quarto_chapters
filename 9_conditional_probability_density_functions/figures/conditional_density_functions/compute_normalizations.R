# Partition Boundary
boundary <- function(x) {
  0.05 * x**3 - 0.2 * x +- 0.175
}

delta <- 1e-3
xs <- seq(-2.5, 2.5, delta)

plot(0, xlim=c(-2.5, 2.5), ylim=c(-2.5, 2.5), type="n")

for (d in c(-2, -1, 0, 1, 2, 3)) {
  ys <- sapply(xs, function(x) boundary(x) + d)
  lines(xs, ys)
}

# Discrete Partition Normalizations
delta <- 1e-2
xs <- seq(-6, +6, delta)
N <- length(xs)
ys <- xs

ds <- c(-2, -1, 0, 1, 2, 3)
norms <- c()

for (d in ds) {
  zs <- matrix(0, nrow=N, ncol=N)
  
  for (n in 1:N) {
    for (m in 1:N) {
      ymin <- boundary(xs[n]) + (d - 1)
      ymax <- boundary(xs[n]) + (d - 0)
      if (ymin < ys[m] & ys[m] < ymax) {
        zs[n, m] <- dnorm(xs[n], 0, 1) * dnorm(ys[m], 0, 1)
      }
    }
  }
  norms <- c(norms, sum(zs) * delta**2)
}

cat(sprintf("%.3f/%.3f,", ds, 1 / norms), "\n")

# Continuous Partition Normalizations
delta <- 1e-2
xs <- seq(-6, +6, delta)
N <- length(xs)

ds <- seq(-2.5, 3, 0.1)
norms <- c()

for (d in ds) {
  deltas <- rep(0, N - 1)
  zs <- rep(0, N - 1)
  
  for (n in 2:N) {
    y <- boundary(xs[n]) + d
    
    dx <- xs[n] - xs[n - 1]
    dy <- y - (boundary(xs[n - 1]) + d)
    deltas[n - 1] <- sqrt(dx**2 + dy**2)
    
    zs[n - 1] <- dnorm(xs[n], 0, 1) * dnorm(y, 0, 1)
  }
  norms <- c(norms, sum(zs * deltas))
}

cat(sprintf("%.3f/%.3f,", ds, 1 / norms), "\n")
