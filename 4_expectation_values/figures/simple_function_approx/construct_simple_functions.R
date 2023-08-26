f <- function(x) { 
  3.6 * exp(-1 * (x - 0.65)**2 ) + 3.0 * exp(-1.75 * (x + 1.25)**2)
}

xs <- seq(-3, 0, 0.001)
fs <- sapply(xs, function(x) f(x))
fmax1 <- max(fs)
xmax1 <- xs[which(fs == fmax1)]

xs <- seq(-1, 1, 0.001)
fs <- sapply(xs, function(x) f(x))
fmin <- min(fs)
xmin <- xs[which(fs == fmin)]

xs <- seq(0, 3, 0.001)
fs <- sapply(xs, function(x) f(x))
fmax2 <- max(fs)
xmax2 <- xs[which(fs == fmax2)]

# Fine
xs <- seq(-3, 3, 0.1)
plot(xs, sapply(xs, function(x) f(x)), type="l")

delta <- 0.75
heights <- seq(0.0, 4, delta)

for (height in heights) {
  abline(h=height)
}

ascending <- c(rep(TRUE, sum(0 < heights & heights < fmax1)),
               rep(FALSE, sum(fmin - delta < heights & heights < fmax1)),
               rep(TRUE, sum(fmin < heights & heights < fmax2)),
               rep(FALSE, sum(0 < heights & heights < fmax2 )))

search_ls <- c(rep(-3, sum(0 < heights & heights < fmax1)),
               rep(xmax1, sum(fmin - delta < heights & heights < fmax1)),
               rep(xmin, sum(fmin < heights & heights < fmax2)),
               rep(xmax2, sum(0 < heights & heights < fmax2 )))

search_us <- c(rep(xmax1, sum(0 < heights & heights < fmax1)),
               rep(xmin, sum(fmin - delta < heights & heights < fmax1)),
               rep(xmax2, sum(fmin < heights & heights < fmax2)),
               rep(3, sum(0 < heights & heights < fmax2 )))

heights <- c(heights[0 < heights & heights < fmax1],
             rev(heights[fmin - delta < heights & heights < fmax1]),
             heights[fmin < heights & heights < fmax2],
             rev(heights[0 < heights & heights < fmax2]))

l <- -3
h <- 0

N <- length(heights)

# Skip root-finding around local minimum
skip <- 0
for (n in 2:(N - 1)) {
  if (heights[n] < heights[n - 1] & heights[n] < heights[n + 1]) 
    skip <- n
}
  
sf_lines <- list()

for (n in 1:N) {
  if (n != skip) {
    x <- uniroot(function(x) f(x) - heights[n], 
                 c(search_ls[n], search_us[n]), maxiter = 10000)$root
    if (ascending[n])
      sf_lines[[n]] <- c(l, x, h)
    else
      sf_lines[[n]] <- c(l, x, heights[n])
  }
  l <- x
  h <- heights[n]
}

sf_lines[[N + 1]] <- c(l, 3, 0)

xs <- seq(-3, 3, 0.1)
plot(xs, sapply(xs, function(x) f(x)), type="l")

for (line in sf_lines) {
  lines(c(line[1], line[2]), c(line[3], line[3]), col="red")
  cat(sprintf("%f/%f/%f, ", line[1], line[2], line[3]))
}

# Finer
xs <- seq(-3, 3, 0.1)
plot(xs, sapply(xs, function(x) f(x)), type="l")

delta <- 0.25
heights <- seq(0.0, 4, delta)

for (height in heights) {
  abline(h=height)
}

ascending <- c(rep(TRUE, sum(0 < heights & heights < fmax1)),
               rep(FALSE, sum(fmin - delta < heights & heights < fmax1)),
               rep(TRUE, sum(fmin < heights & heights < fmax2)),
               rep(FALSE, sum(0 < heights & heights < fmax2 )))

search_ls <- c(rep(-3, sum(0 < heights & heights < fmax1)),
               rep(xmax1, sum(fmin - delta < heights & heights < fmax1)),
               rep(xmin, sum(fmin < heights & heights < fmax2)),
               rep(xmax2, sum(0 < heights & heights < fmax2 )))

search_us <- c(rep(xmax1, sum(0 < heights & heights < fmax1)),
               rep(xmin, sum(fmin - delta < heights & heights < fmax1)),
               rep(xmax2, sum(fmin < heights & heights < fmax2)),
               rep(3, sum(0 < heights & heights < fmax2 )))

heights <- c(heights[0 < heights & heights < fmax1],
             rev(heights[fmin - delta < heights & heights < fmax1]),
             heights[fmin < heights & heights < fmax2],
             rev(heights[0 < heights & heights < fmax2]))

l <- -3
h <- 0

N <- length(heights)

# Skip root-finding around local minimum
skip <- 0
for (n in 2:(N - 1)) {
  if (heights[n] < heights[n - 1] & heights[n] < heights[n + 1]) 
    skip <- n
}

sf_lines <- list()

for (n in 1:N) {
  if (n != skip) {
    x <- uniroot(function(x) f(x) - heights[n], 
                 c(search_ls[n], search_us[n]), maxiter = 10000)$root
    if (ascending[n])
      sf_lines[[n]] <- c(l, x, h)
    else
      sf_lines[[n]] <- c(l, x, heights[n])
  }
  l <- x
  h <- heights[n]
}

sf_lines[[N + 1]] <- c(l, 3, 0)

xs <- seq(-3, 3, 0.1)
plot(xs, sapply(xs, function(x) f(x)), type="l")

for (line in sf_lines) {
  lines(c(line[1], line[2]), c(line[3], line[3]), col="red")
  cat(sprintf("%f/%f/%f, ", line[1], line[2], line[3]))
}


# Finest
xs <- seq(-3, 3, 0.1)
plot(xs, sapply(xs, function(x) f(x)), type="l")

delta <- 0.05
heights <- seq(0.0, 4, delta)

for (height in heights) {
  abline(h=height)
}

ascending <- c(rep(TRUE, sum(0 < heights & heights < fmax1)),
               rep(FALSE, sum(fmin - delta < heights & heights < fmax1)),
               rep(TRUE, sum(fmin < heights & heights < fmax2)),
               rep(FALSE, sum(0 < heights & heights < fmax2 )))

search_ls <- c(rep(-3, sum(0 < heights & heights < fmax1)),
               rep(xmax1, sum(fmin - delta < heights & heights < fmax1)),
               rep(xmin, sum(fmin < heights & heights < fmax2)),
               rep(xmax2, sum(0 < heights & heights < fmax2 )))

search_us <- c(rep(xmax1, sum(0 < heights & heights < fmax1)),
               rep(xmin, sum(fmin - delta < heights & heights < fmax1)),
               rep(xmax2, sum(fmin < heights & heights < fmax2)),
               rep(3, sum(0 < heights & heights < fmax2 )))

heights <- c(heights[0 < heights & heights < fmax1],
             rev(heights[fmin - delta < heights & heights < fmax1]),
             heights[fmin < heights & heights < fmax2],
             rev(heights[0 < heights & heights < fmax2]))

l <- -3
h <- 0

N <- length(heights)

# Skip root-finding around local minimum
skip <- 0
for (n in 2:(N - 1)) {
  if (heights[n] < heights[n - 1] & heights[n] < heights[n + 1]) 
    skip <- n
}

sf_lines <- list()

for (n in 1:N) {
  if (n != skip) {
    x <- uniroot(function(x) f(x) - heights[n], 
                 c(search_ls[n], search_us[n]), maxiter = 10000)$root
    if (ascending[n])
      sf_lines[[n]] <- c(l, x, h)
    else
      sf_lines[[n]] <- c(l, x, heights[n])
  }
  l <- x
  h <- heights[n]
}

sf_lines[[N + 1]] <- c(l, 3, 0)

xs <- seq(-3, 3, 0.1)
plot(xs, sapply(xs, function(x) f(x)), type="l")

for (line in sf_lines) {
  lines(c(line[1], line[2]), c(line[3], line[3]), col="red")
  cat(sprintf("%f/%f/%f, ", line[1], line[2], line[3]))
}