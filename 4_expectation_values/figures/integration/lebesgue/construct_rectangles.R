f <- function(x) { 
  1.5 * (sin(100 * pi / 180 * x) + 0.25 * x)
}

xs <- seq(-3, 3, 0.001)
fs <- sapply(xs, function(x) f(x))
plot(xs, fs, type="l")

z1 <- uniroot(function(x) f(x), c(-3, -1), maxiter = 10000)$root
z2 <- uniroot(function(x) f(x), c(-1, +1), maxiter = 10000)$root
z3 <- uniroot(function(x) f(x), c(+1, +3), maxiter = 10000)$root

abline(v=z1)
abline(v=z2)
abline(v=z3)

bins <- c(-3, z1, z2, z3, +3)
signs <- c(1, -1, 1, -1)

# Fine
delta <- 0.5

rects <- list()
n <- 1

for (b in 1:4) {
  xs <- seq(bins[b], bins[b + 1], 0.001)
  absfs <- sapply(xs, function(x) abs(f(x)))

  fmax <- max(absfs)
  xmax <- xs[which(absfs == fmax)]
            
  heights <- seq(delta, fmax, delta)
  h <- 0
  
  for (height in heights) {
    if (height < abs(f(bins[b]))) {
      l <- bins[b]
      u <- uniroot(function(x) abs(f(x)) - height, c(xmax, bins[b + 1]), maxiter = 10000)$root
    } else if (height < abs(f(bins[b + 1]))) {
      l <- uniroot(function(x) abs(f(x)) - height, c(bins[b], xmax), maxiter = 10000)$root
      u <- bins[b + 1]
    } else {
      l <- uniroot(function(x) abs(f(x)) - height, c(bins[b], xmax), maxiter = 10000)$root
      u <- uniroot(function(x) abs(f(x)) - height, c(xmax, bins[b + 1]), maxiter = 10000)$root
    }
    rects[[n]] <- c(l, u, h, signs[b] * height)
    
    h <- signs[b] * height
    n <- n + 1
  }
}

xs <- seq(-3, 3, 0.001)
fs <- sapply(xs, function(x) f(x))
plot(xs, fs, type="l")

for (r in rects) {
  rect(r[1], r[3], r[2], r[4], col="red", border=FALSE)
  cat(sprintf("%f/%f/%f/%f, ", r[1], r[3], r[2], r[4]))
}

# Finer
delta <- 0.25

rects <- list()
n <- 1

for (b in 1:4) {
  xs <- seq(bins[b], bins[b + 1], 0.001)
  absfs <- sapply(xs, function(x) abs(f(x)))
  
  fmax <- max(absfs)
  xmax <- xs[which(absfs == fmax)]
  
  heights <- seq(delta, fmax, delta)
  h <- 0
  
  for (height in heights) {
    if (height < abs(f(bins[b]))) {
      l <- bins[b]
      u <- uniroot(function(x) abs(f(x)) - height, c(xmax, bins[b + 1]), maxiter = 10000)$root
    } else if (height < abs(f(bins[b + 1]))) {
      l <- uniroot(function(x) abs(f(x)) - height, c(bins[b], xmax), maxiter = 10000)$root
      u <- bins[b + 1]
    } else {
      l <- uniroot(function(x) abs(f(x)) - height, c(bins[b], xmax), maxiter = 10000)$root
      u <- uniroot(function(x) abs(f(x)) - height, c(xmax, bins[b + 1]), maxiter = 10000)$root
    }
    rects[[n]] <- c(l, u, h, signs[b] * height)
    
    h <- signs[b] * height
    n <- n + 1
  }
}

xs <- seq(-3, 3, 0.001)
fs <- sapply(xs, function(x) f(x))
plot(xs, fs, type="l")

for (r in rects) {
  rect(r[1], r[3], r[2], r[4], col="red", border=FALSE)
  cat(sprintf("%f/%f/%f/%f, ", r[1], r[3], r[2], r[4]))
}

# Finest
delta <- 0.1

rects <- list()
n <- 1

for (b in 1:4) {
  xs <- seq(bins[b], bins[b + 1], 0.001)
  absfs <- sapply(xs, function(x) abs(f(x)))
  
  fmax <- max(absfs)
  xmax <- xs[which(absfs == fmax)]
  
  heights <- seq(delta, fmax, delta)
  h <- 0
  
  for (height in heights) {
    if (height < abs(f(bins[b]))) {
      l <- bins[b]
      u <- uniroot(function(x) abs(f(x)) - height, c(xmax, bins[b + 1]), maxiter = 10000)$root
    } else if (height < abs(f(bins[b + 1]))) {
      l <- uniroot(function(x) abs(f(x)) - height, c(bins[b], xmax), maxiter = 10000)$root
      u <- bins[b + 1]
    } else {
      l <- uniroot(function(x) abs(f(x)) - height, c(bins[b], xmax), maxiter = 10000)$root
      u <- uniroot(function(x) abs(f(x)) - height, c(xmax, bins[b + 1]), maxiter = 10000)$root
    }
    rects[[n]] <- c(l, u, h, signs[b] * height)
    
    h <- signs[b] * height
    n <- n + 1
  }
}

xs <- seq(-3, 3, 0.001)
fs <- sapply(xs, function(x) f(x))
plot(xs, fs, type="l")

for (r in rects) {
  rect(r[1], r[3], r[2], r[4], col="red", border=FALSE)
  cat(sprintf("%f/%f/%f/%f, ", r[1], r[3], r[2], r[4]))
}