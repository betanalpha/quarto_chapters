################################################################################
#
# The code is copyright 2024 Michael Betancourt and licensed under the
# new BSD (3-clause) license:
#  https://opensource.org/licenses/BSD-3-Clause
#
# For more information see
#  https://github.com/betanalpha/mcmc_visualization_tools.
#
# Requires https://github.com/betanalpha/mcmc_diagnostics.
#
################################################################################

# Load required libraries
library(colormap)

# Graphic configuration
c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

c_light_teal <- c("#6B8E8E")
c_mid_teal <- c("#487575")
c_dark_teal <- c("#1D4F4F")

util <- new.env()
if (file.exists('mcmc_analysis_tools_rstan.R')) {
  source('mcmc_analysis_tools_rstan.R', local=util)
} else if (file.exists('mcmc_analysis_tools_other.R')) {
  source('mcmc_analysis_tools_other.R', local=util)
} else {
  stop(print0('mcmc_visualization_tools.R requires that ',
              'mcmc_analysis_tools_[rstan/other].R from ',
              'https://github.com/betanalpha/mcmc_diagnostics ',
              'is available.'))
}

################################################################################
# Utility Functions
################################################################################

# Emit stop if the two input arrays are not the same length.
# @param a First list
# @param a_name Name of first list
# @param b Second list
# @param b_name Name of second list
check_dimensions <- function(a, a_name, b, b_name) {
  if (length(a) != length(b)) {
    stop(sprintf('The arguments `%s` and `%s` are not the same length!',
                 a_name, b_name))
  }
}

# Check if the given expectand names exist in the samples object, removing any
# that do not.  Emit warning describing any missing names and emit stop if no
# expectand names are found.
# @param names A one-dimensional array of strings
# @param samples A named list
# @return One-dimensional array of valid expectand names
check_expectand_names <- function(names, samples) {
  all_names <- names(samples)

  bad_names <- setdiff(names, all_names)
  B <- length(bad_names)
  if (B > 0) {
    if (B == 1)
      warning(paste0(sprintf('The expectand name %s is not in the',
                             bad_names[1]),
                     ' `samples` object and will be ignored.'))
    else
      warning(paste0(sprintf('The expectand names %s are not in the',
                           paste(bad_names, collapse=", ")),
                     ' `samples` object and will be ignored.'))
  }

  good_names <- intersect(names, all_names)
  if (length(good_names) == 0) {
    stop('There are no valid expectand names.')
  }
  good_names
}

# Check how many values fall below bin_min and above bin_max and emit
# appropriate warming message.
# @param bin_min Lower threshold
# @param bin_max Upper threshold
# @param values A one-dimensional array of values to check
# @param name Value description
check_bin_containment <- function(bin_min, bin_max, values, name="value") {
  N <- length(values)

  N_low <- sum(values < bin_min)
  if (N_low > 0)
    if (N_low == 1)
      warning(sprintf('%i %s (%.1f%%) fell below the binning.',
                      N_low, name, 100 * N_low / N))
    else
      warning(sprintf('%i %ss (%.1f%%) fell below the binning.',
                      N_low, name, 100 * N_low / N))

  N_high <- sum(bin_max < values)
  if (N_high > 0)
    if (N_high == 1)
      warning(sprintf('%i %s (%.1f%%) fell above the binning.',
                      N_high, name, 100 * N_high / N))
    else
      warning(sprintf('%i %ss (%.1f%%) fell above the binning.',
                      N_high, name, 100 * N_high / N))
}

# Configure binning.  Any null arguments are automatically configured to match
# the behavior of the `values1` and `values2` arguments.  Additionally if the
# difference `bin_max - bin_min` does not divide `bin_delta` then `bin_min` and
# `bin_max` are tweaked to respectively smaller and larger values as needed.
# @param bin_min Lower threshold
# @param bin_max Upper threshold
# @param bin_delta Bin width
# @param values1 An array of values
# @param values1 An auxiliary array of values
# @return One-dimensional array of updated bin_min, bin_max, and bin_delta
#         values
configure_bins <- function(bin_min, bin_max, bin_delta,
                           values1, values2=NULL) {
  if (is.null(values2)) {
    # Adapt bin configuration to `values1`
    if (is.null(bin_min))
      bin_min <- min(values1)
    if (is.null(bin_max))
      bin_max <- max(values1)
  } else {
    # Adapt bin configuration to `values1` and `values2`
    if (is.null(bin_min))
      bin_min <- min(min(values1), min(values2))
    if (is.null(bin_max))
      bin_max <- max(max(values1), min(values2))
  }

  if (is.null(bin_delta))
    bin_delta <- (bin_max - bin_min) / 25

  # Tweak bin configuration so that `bin_delta`
  # evenly divides `bin_max - bin_min`
  N <- (bin_max - bin_min) / bin_delta
  excess <- N - floor(N)
  if (excess > 1e-15) {
    bin_min <- bin_min - 0.5 * bin_delta * excess
    bin_max <- bin_max + 0.5 * bin_delta * excess
  }

  c(bin_min, bin_max, bin_delta)
}

# Compute bin plotting.
# @param breaks Bin edges
# @return List of plotting indices and positions
configure_bin_plotting <- function(breaks) {
  B <- length(breaks) - 1
  idxs <- rep(1:B, each=2)
  xs <- sapply(1:length(idxs),
               function(b) if(b %% 2 == 1) breaks[idxs[b]]
                           else            breaks[idxs[b] + 1])
  list(idxs, xs)
}

################################################################################
# Data Visualizations
################################################################################

# Plot histogram outline.
# @param values Values that comprise the histogram
# @param bin_min Lower threshold
# @param bin_max Upper threshold
# @param bin_delta Bin width
# @param prob Boolean determining whether bin contents should be normalized so
#             that the histogram approximates a probability density function;
#             defaults to FALSE
# @param col Color of histogram; defaults to "black"
# @param add Boolean determining whether to add histogram outline to existing
#            plot or to create new axes; defaults to FALSE
# @param xlab Label for x-axis; defaults to empty string.
# @param main Plot title; defaults to empty string.
plot_line_hist <- function(values,
                           bin_min=NULL, bin_max=NULL, bin_delta=NULL,
                           prob=FALSE, col="black", add=FALSE,
                           xlab="", main="") {
  # Remove any NA values
  values <- values[!is.na(values)]

  # Construct binning configuration
  bin_config <- configure_bins(bin_min, bin_max, bin_delta, values)
  bin_min <- bin_config[1]
  bin_max <- bin_config[2]
  bin_delta <- bin_config[3]

  # Construct bins
  breaks <- seq(bin_min, bin_max, bin_delta)
  plot_config <- configure_bin_plotting(breaks)
  plot_idxs <- plot_config[[1]]
  plot_xs <- plot_config[[2]]

  # Check bin containment
  check_bin_containment(bin_min, bin_max, values)

  # Compute bin contents
  counts <- hist(values[bin_min <= values & values <= bin_max],
                 breaks=breaks, plot=FALSE)$counts

  ylab <- "Counts"
  if (prob) {
    ylab <- "Empirical Bin Probability / Bin Width"
    counts <- counts / (bin_delta * sum(counts))
  }

  if (add) {
    lines(plot_xs, counts[plot_idxs], col=col, lwd=2)
  } else {
    # Plot
    plot(plot_xs, counts[plot_idxs], main=main,
         type="l", col=col, lwd=2,
         xlab=xlab, xlim=c(bin_min, bin_max),
         ylab=ylab, ylim=c(0, 1.1 * max(counts)))
  }
}

# Plot the overlay of two line histograms.
# @param values1 Values that comprise the first histogram
# @param values2 Values that comprise the second histogram
# @param bin_min Lower threshold
# @param bin_max Upper threshold
# @param bin_delta Bin width
# @param prob Boolean determining whether bin contents should be normalized so
#             that the histogram approximates a probability density function;
#             defaults to FALSE
# @param xlab Label for x-axis; defaults to empty string.
# @param main Plot title; defaults to empty string.
# @param col1 Color of first histogram; defaults to "black"
# @param col2 Color of second histogram; defaults to c_mid_teal
plot_line_hists <- function(values1, values2,
                            bin_min=NULL, bin_max=NULL, bin_delta=NULL,
                            prob=FALSE, xlab="y", main="",
                            col1="black", col2=c_mid_teal) {
  # Remove any NA values
  values1 <- values1[!is.na(values1)]
  values2 <- values2[!is.na(values2)]

  # Construct binning configuration
  bin_config <- configure_bins(bin_min, bin_max, bin_delta,
                               values1, values2)
  bin_min <- bin_config[1]
  bin_max <- bin_config[2]
  bin_delta <- bin_config[3]

  # Construct bins
  breaks <- seq(bin_min, bin_max, bin_delta)
  plot_config <- configure_bin_plotting(breaks)
  plot_idxs <- plot_config[[1]]
  plot_xs <- plot_config[[2]]

  # Check bin containment
  check_bin_containment(bin_min, bin_max, values1)
  check_bin_containment(bin_min, bin_max, values2)

  # Compute bin contents
  counts1 <- hist(values1[bin_min <= values1 & values1 <= bin_max],
                  breaks=breaks, plot=FALSE)$counts
  counts2 <- hist(values2[bin_min <= values2 & values2 <= bin_max],
                  breaks=breaks, plot=FALSE)$counts

  ylab <- "Counts"
  if (prob) {
    ylab <- "Empirical Bin Probability / Bin Width"
    counts1 <- counts1 / (delta * sum(counts1))
    counts2 <- counts2 / (delta * sum(counts1))
  }

  # Plot
  ymax <- 1.1 * max(max(counts1), max(counts2))

  plot(plot_xs, counts1[plot_idxs], main=main,
       type="l", col=col1, lwd=2,
       xlab=xlab, xlim=c(bin_min, bin_max),
       ylab=ylab, ylim=c(0, ymax))
  lines(plot_xs, counts2[plot_idxs], col="white", lwd=4)
  lines(plot_xs, counts2[plot_idxs], col=col2, lwd=2)
}

################################################################################
# Pushforward Visualizations
################################################################################

# Overlay nested quantile intervals to visualize an ensemble of histograms.
# Individual quantiles are estimated as the average of the empirical quantiles
# across each Markov chain, a consistent quantile estimator for Markov chain
# Monte Carlo.
# @param samples A named list of two-dimensional arrays for
#                each expectand.  The first dimension of each element
#                indexes the Markov chains and the second dimension
#                indexes the sequential states within each Markov chain.
# @param val_name_prefix Prefix for the relevant variable names
# @param bin_min Lower threshold
# @param bin_max Upper threshold
# @param bin_delta Bin width
# @param baseline_values Baseline values for constructing a baseline histogram;
#                        defaults to NULL
# @param baseline_col Color for plotting baseline value; defaults to "black"
# @param xlab Label for x-axis; defaults to empty string
# @param display_ylim Plot limits for y-axis; defaults to NULL
# @param main Plot title; defaults to empty string
plot_hist_quantiles <- function(samples, val_name_prefix,
                                bin_min=NULL, bin_max=NULL, bin_delta=NULL,
                                baseline_values=NULL, baseline_col="black",
                                xlab="", display_ylim=NULL, main="") {
  # Construct relevant variable names and format corresponding values.
  # Order of the variables does not affect the shape of the histogram.
  names <- grep(paste0('^', val_name_prefix, '\\['),
                names(samples), value=TRUE)
  collapsed_values <- c(sapply(names, function(name) c(t(samples[[name]]),
                                                       recursive=TRUE)))

  # Construct binning configuration
  if (is.null(baseline_values))
    bin_config <- configure_bins(bin_min, bin_max, bin_delta,
                                 collapsed_values)
  else
    bin_config <- configure_bins(bin_min, bin_max, bin_delta,
                                 collapsed_values, baseline_values)
  bin_min <- bin_config[1]
  bin_max <- bin_config[2]
  bin_delta <- bin_config[3]

  # Construct bins
  breaks <- seq(bin_min, bin_max, bin_delta)
  plot_config <- configure_bin_plotting(breaks)
  plot_idxs <- plot_config[[1]]
  plot_xs <- plot_config[[2]]

  # Check bin containment
  check_bin_containment(bin_min, bin_max, collapsed_values,
                        "predictive value")
  if (!is.null(baseline_values))
    check_bin_containment(bin_min, bin_max, baseline_values,
                          "observed value")

  # Construct quantiles for bin contents
  B <- length(breaks) - 1
  probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

  counts <- rep(NA, B)
  quantiles <- matrix(0, nrow=9, ncol=B)

  bin_count <- function(x, b_low, b_high) {
    sum(b_low <= x & x < b_high)
  }

  bin_counters <- list()
  for (b in 1:B) {
    bin_counters[[b]] <-
      local({
        b_low <- breaks[b];
        b_high <- breaks[b + 1];
        function(x) bin_count(x, b_low, b_high)
      })

    if (!is.null(baseline_values))
      counts[b] = bin_counters[[b]](baseline_values)
  }

  bin_count_samples <-
    util$eval_expectand_pushforwards(samples,
                                     bin_counters,
                                     list('x'=array(names)))

  probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  quantiles <- sapply(bin_count_samples,
                      function(s)
                      util$ensemble_mcmc_quantile_est(s, probs))

  plot_quantiles <- do.call(cbind, lapply(plot_idxs,
                                          function(n) quantiles[1:9, n]))

  # Plot
  if (is.null(display_ylim)) {
    if (is.null(baseline_values)) {
      display_ylim <- c(0, max(quantiles[9,]))
    }
    else {
      display_ylim <- c(0, max(max(quantiles[9,]), max(counts)))
    }
  }

  plot(1, type="n", main=main,
       xlim=c(bin_min, bin_max), xlab=xlab,
       ylim=display_ylim, ylab="Counts")

  polygon(c(plot_xs, rev(plot_xs)),
          c(plot_quantiles[1,], rev(plot_quantiles[9,])),
          col = c_light, border = NA)
  polygon(c(plot_xs, rev(plot_xs)),
          c(plot_quantiles[2,], rev(plot_quantiles[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(plot_xs, rev(plot_xs)),
          c(plot_quantiles[3,], rev(plot_quantiles[7,])),
          col = c_mid, border = NA)
  polygon(c(plot_xs, rev(plot_xs)),
          c(plot_quantiles[4,], rev(plot_quantiles[6,])),
          col = c_mid_highlight, border = NA)
  for (b in 1:B) {
    idx1 <- 2 * b - 1
    idx2 <- 2 * b
    lines(plot_xs[idx1:idx2], plot_quantiles[5,idx1:idx2],
          col=c_dark, lwd=2)
  }

  if (!is.null(baseline_values)) {
    plot_counts <- do.call(cbind,
                           lapply(plot_idxs, function(n) counts[n]))

    lines(plot_xs, plot_counts, col="white", lty=1, lw=4)
    lines(plot_xs, plot_counts, col=baseline_col, lty=1, lw=2)
  }
}

# Overlay disconnected nested quantile intervals to visualize an ensemble of
# one-dimensional pushforward distributions.
# Individual quantiles are estimated as the average of the empirical quantiles
# across each Markov chain, a consistent quantile estimator for Markov chain
# Monte Carlo.
# @param samples A named list of two-dimensional arrays for
#                each expectand.  The first dimension of each element
#                indexes the Markov chains and the second dimension
#                indexes the sequential states within each Markov chain.
# @param names List of relevant variable names
# @param baseline_values Baseline values; defaults to NULL
# @param baseline_col Color for plotting baseline value; defaults to "black"
# @params residual Boolean value indicating whether to overlay quantiles and
#                  baseline values or plot their differences
# @param xlab Label for x-axis; defaults to empty string
# @param xticklabs Labels for x-axis tics; defaults to NULL
# @param ylab Label for y-axis; defaults to NULL
# @param display_ylim Plot limits for y-axis; defaults to NULL
# @param main Plot title; defaults to empty string
plot_disc_pushforward_quantiles <- function(samples, names,
                                            baseline_values=NULL,
                                            baseline_col="black",
                                            residual=FALSE,
                                            xlab="", xticklabs=NULL,
                                            ylab=NULL, display_ylim=NULL,
                                            main="") {
  # Check that baseline values are well-defined
  if (!is.null(baseline_values)) {
    if (length(baseline_values) != length(names)) {
      warning(paste0('The list of baseline values has the wrong',
                     ' dimension.  Baselines will not be plotted.'))
      baseline_values <- NULL
    }
  }

  # Check that names are in samples
  names <- check_expectand_names(names, samples)

  # Construct bins
  N <- length(names)
  bin_min <- 0.5
  bin_max <- N + 0.5
  bin_delta <- 1
  breaks <- seq(bin_min, bin_max, bin_delta)

  plot_config <- configure_bin_plotting(breaks)
  plot_idxs <- plot_config[[1]]
  plot_xs <- plot_config[[2]]

  # Construct marginal quantiles
  probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

  if(!is.null(baseline_values) & residual) {
    calc <- function(n) {
      util$ensemble_mcmc_quantile_est(samples[[names[n]]] -
                                      baseline_values[n],
                                      probs)
    }
    quantiles <- sapply(1:N, calc)
  } else {
    calc <- function(n) {
      util$ensemble_mcmc_quantile_est(samples[[names[n]]], probs)
    }
    quantiles <- sapply(1:N, calc)
  }

  plot_quantiles <- do.call(cbind, lapply(plot_idxs,
                                          function(n) quantiles[1:9, n]))

  # Plot
  if (is.null(display_ylim)) {
    if (is.null(baseline_values)) {
      display_ylim <- c(min(quantiles[1,]),
                        max(quantiles[9,]))
    }
    else {
      if (residual) {
        display_ylim <- c(min(c(0, quantiles[1,])),
                          max(c(0, quantiles[9,])))
      } else {
        display_ylim <- c(min(min(quantiles[1,]),
                              min(baseline_values)),
                          max(max(quantiles[9,]),
                              max(baseline_values)))
      }
    }
    delta <- 0.05 * (display_ylim[2] - display_ylim[1])
    display_ylim[1] <- display_ylim[1] - delta
    display_ylim[2] <- display_ylim[2] + delta
  }

  if (is.null(ylab)) {
    if (is.null(baseline_values) | !residual)
      ylab <- "Marginal Quantiles"
    else
      ylab <- "Marginal Quantiles - Baselines"
  }

  if (is.null(xticklabs)) {
    plot(1, type="n", main=main,
         xlim=c(bin_min, bin_max), xlab=xlab,
         ylim=display_ylim, ylab=ylab)
  } else {
    if (length(xticklabs) == N) {
      plot(1, type="n", main=main,
           xlim=c(bin_min, bin_max), xlab=xlab, xaxt="n",
           ylim=display_ylim, ylab=ylab)
      axis(1, at=1:N, labels=xticklabs)
    } else {
      warning(paste0('The list of x labels tick has the wrong',
                     ' dimension and baselines will not be plotted.'))
      plot(1, type="n", main=main,
           xlim=c(bin_min, bin_max), xlab=xlab,
           ylim=display_ylim, ylab=ylab)
    }
  }

  polygon(c(plot_xs, rev(plot_xs)),
          c(plot_quantiles[1,], rev(plot_quantiles[9,])),
          col = c_light, border = NA)
  polygon(c(plot_xs, rev(plot_xs)),
          c(plot_quantiles[2,], rev(plot_quantiles[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(plot_xs, rev(plot_xs)),
          c(plot_quantiles[3,], rev(plot_quantiles[7,])),
          col = c_mid, border = NA)
  polygon(c(plot_xs, rev(plot_xs)),
          c(plot_quantiles[4,], rev(plot_quantiles[6,])),
          col = c_mid_highlight, border = NA)
  for (n in 1:N) {
    idx1 <- 2 * n - 1
    idx2 <- 2 * n
    lines(plot_xs[idx1:idx2], plot_quantiles[5, idx1:idx2],
          col=c_dark, lwd=2)
  }

  if (!is.null(baseline_values)) {
    if (residual) {
      abline(h=0, col="#DDDDDD", lwd=2, lty=3)
    } else {
      for (n in 1:N) {
        idx1 <- 2 * n - 1
        idx2 <- 2 * n
        lines(plot_xs[idx1:idx2], rep(baseline_values[n], 2),
              col="white", lwd=4)
        lines(plot_xs[idx1:idx2], rep(baseline_values[n], 2),
              col=baseline_col, lwd=2)
      }
    }
  }
}

# Overlay connected nested quantile intervals to visualize an ensemble of
# one-dimensional pushforward distributions.
# Individual quantiles are estimated as the average of the empirical quantiles
# across each Markov chain, a consistent quantile estimator for Markov chain
# Monte Carlo.
# @param samples A named list of two-dimensional arrays for
#                each expectand.  The first dimension of each element
#                indexes the Markov chains and the second dimension
#                indexes the sequential states within each Markov chain.
# @param names List of relevant variable names
# @param plot_xs One-dimensional array of x-axis values
#                associated with each variable.
# @param baseline_values Baseline values; defaults to NULL
# @param baseline_col Color for plotting baseline value; defaults to "black"
# @params residual Boolean value indicating whether to overlay quantiles and
#                  baseline values or plot their differences
# @param xlab Label for x-axis; defaults to empty string
# @param display_xlim Plot limits for x-axis; defaults to NULL
# @param ylab Label for y-axis; defaults to NULL
# @param display_ylim Plot limits for y-axis; defaults to NULL
# @param main Plot title; defaults to empty string
plot_conn_pushforward_quantiles <- function(samples, names, plot_xs,
                                            baseline_values=NULL,
                                            baseline_col="black",
                                            residual=FALSE,
                                            xlab="", display_xlim=NULL,
                                            ylab=NULL, display_ylim=NULL,
                                            main="") {
  # Check dimensions
  check_dimensions(plot_xs, 'plot_xs', names, 'names')

  # Check that baseline values are well-defined
  if (!is.null(baseline_values)) {
    if (length(baseline_values) != length(names)) {
      warning(paste0('The list of baseline values has the wrong',
                     ' dimension.  Baselines will not be plotted.'))
      baseline_values <- NULL
    }
  }

  # Check that names are in samples
  names <- check_expectand_names(names, samples)

  # Construct quantiles for bin contents
  N <- length(names)
  probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

  if(!is.null(baseline_values) & residual) {
    calc <- function(n) {
      util$ensemble_mcmc_quantile_est(samples[[names[n]]] -
                                      baseline_values[n],
                                      probs)
    }
    plot_quantiles <- sapply(1:N, calc)
  } else {
    calc <- function(n) {
      util$ensemble_mcmc_quantile_est(samples[[names[n]]], probs)
    }
    plot_quantiles <- sapply(1:N, calc)
  }

  # Plot
  if (is.null(display_xlim))
    display_xlim <- range(plot_xs)

  if (is.null(display_ylim)) {
    if (is.null(baseline_values)) {
      display_ylim <- c(min(plot_quantiles[1,]),
                        max(plot_quantiles[9,]))
    }
    else {
      if (residual) {
        display_ylim <- c(min(c(0, plot_quantiles[1,])),
                          max(c(0, plot_quantiles[9,])))
      } else {
        display_ylim <- c(min(min(plot_quantiles[1,]),
                              min(baseline_values)),
                          max(max(plot_quantiles[9,]),
                              max(baseline_values)))
      }
    }
    delta <- 0.05 * (display_ylim[2] - display_ylim[1])
    display_ylim[1] <- display_ylim[1] - delta
    display_ylim[2] <- display_ylim[2] + delta
  }

  if (is.null(ylab)) {
    if (is.null(baseline_values) | !residual)
      ylab <- "Marginal Quantiles"
    else
      ylab <- "Marginal Quantiles - Baselines"
  }

  plot(1, type="n", main=main,
       xlim=display_xlim, xlab=xlab,
       ylim=display_ylim, ylab=ylab)

  polygon(c(plot_xs, rev(plot_xs)),
          c(plot_quantiles[1,], rev(plot_quantiles[9,])),
          col = c_light, border = NA)
  polygon(c(plot_xs, rev(plot_xs)),
          c(plot_quantiles[2,], rev(plot_quantiles[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(plot_xs, rev(plot_xs)),
          c(plot_quantiles[3,], rev(plot_quantiles[7,])),
          col = c_mid, border = NA)
  polygon(c(plot_xs, rev(plot_xs)),
          c(plot_quantiles[4,], rev(plot_quantiles[6,])),
          col = c_mid_highlight, border = NA)
  lines(plot_xs, plot_quantiles[5,], col=c_dark, lwd=2)

  if (!is.null(baseline_values)) {
    if (residual) {
      abline(h=0, col="#DDDDDD", lwd=2, lty=3)
    } else {
      lines(plot_xs, baseline_values, col="white", lwd=4)
      lines(plot_xs, baseline_values, col=baseline_col, lwd=2)
    }
  }
}

# Overlay an ensemble of function realizations to visualize a probability
# distribution over a space of one-dimensional functions.
# @param samples A named list of two-dimensional arrays for
#                each expectand.  The first dimension of each element
#                indexes the Markov chains and the second dimension
#                indexes the sequential states within each Markov chain.
# @param names List of relevant variable names
# @param plot_xs One-dimensional array of x-axis values
#                associated with each variable.
# @param N_plots Number of realizations to plot
# @param baseline_values Baseline values; defaults to NULL
# @param baseline_col Color for plotting baseline value; defaults to "black"
# @params residual Boolean value indicating whether to overlay quantiles and
#                  baseline values or plot their differences
# @param xlab Label for x-axis; defaults to empty string
# @param display_xlim Plot limits for x-axis; defaults to NULL
# @param ylab Label for y-axis; defaults to NULL
# @param display_ylim Plot limits for y-axis; defaults to NULL
# @param main Plot title; defaults to empty string
plot_realizations <- function(samples, names, plot_xs, N_plots=50,
                              baseline_values=NULL,
                              baseline_col="black",
                              residual=FALSE,
                              xlab="", display_xlim=NULL,
                              ylab=NULL, display_ylim=NULL,
                              main="") {
  # Check dimensions
  check_dimensions(plot_xs, 'plot_xs', names, 'names')

  # Check that baseline values are well-defined
  if (!is.null(baseline_values)) {
    if (length(baseline_values) != length(names)) {
      warning(paste0('The list of baseline values has the wrong',
                     ' dimension.  Baselines will not be plotted.'))
      baseline_values <- NULL
    }
  }

  # Check that names are in samples
  names <- check_expectand_names(names, samples)

  # Extract function values
  fs <- t(sapply(names, function(name)
                        c(samples[[name]], recursive=TRUE)))

  N <- dim(fs)[1]
  I <- dim(fs)[2]

  if (!is.null(baseline_values) & residual)
    for (i in 1:I)
      fs[,i] <- fs[,i] - baseline_values

  # Configure ensemble of function realizations
  J <- min(N_plots, I)
  plot_idx <- sapply(1:J, function(j) (I %/% J) * (j - 1) + 1)

  nom_colors <- c("#DCBCBC", "#C79999", "#B97C7C",
                  "#A25050", "#8F2727", "#7C0000")
  line_colors <- colormap(colormap=nom_colors, nshades=J)

  # Plot
  if (is.null(display_xlim))
    display_xlim <- range(plot_xs)

  if (is.null(display_ylim)) {
    if (is.null(baseline_values)) {
      display_ylim <- range(fs)
    }
    else {
      if (residual) {
        display_ylim <- range(c(0, fs))
      } else {
        display_ylim <- c(min(min(fs), min(baseline_values)),
                          max(max(fs), max(baseline_values)))
      }
    }
    delta <- 0.05 * (display_ylim[2] - display_ylim[1])
    display_ylim[1] <- display_ylim[1] - delta
    display_ylim[2] <- display_ylim[2] + delta
  }

  if (is.null(ylab)) {
    if (is.null(baseline_values) | !residual)
      ylab <- "Function Outputs"
    else
      ylab <- "Function Outputs - Baselines"
  }

  plot(1, type="n", main=main,
       xlim=display_xlim, xlab=xlab,
       ylim=display_ylim, ylab=ylab)

  for (j in 1:J) {
    f <- fs[, plot_idx[j]]
    lines(plot_xs, f, col=line_colors[j], lwd=3)
  }

  if (!is.null(baseline_values)) {
    if (residual) {
      abline(h=0, col="#DDDDDD", lwd=2, lty=3)
    } else {
      lines(plot_xs, baseline_values, col="white", lwd=4)
      lines(plot_xs, baseline_values, col=baseline_col, lwd=2)
    }
  }
}

# Overlay nested quantile intervals to visualize an ensemble of empirical means
# within the computed bins.
# Individual quantiles are estimated as the average of the empirical quantiles
# across each Markov chain, a consistent quantile estimator for Markov chain
# Monte Carlo.
# @param samples A named list of two-dimensional arrays for
#                each expectand.  The first dimension of each element
#                indexes the Markov chains and the second dimension
#                indexes the sequential states within each Markov chain.
# @param names List of relevant variable names
# @param obs_xs One-dimensional array of observed x-values on which to condition
# @param bin_min Lower threshold for conditioning
# @param bin_max Upper threshold for conditioning
# @param bin_delta Bin width for conditioning
# @param baseline_values Baseline values; defaults to NULL
# @param baseline_col Color for plotting baseline value; defaults to "black"
# @params residual Boolean value indicating whether to overlay quantiles and
#                  baseline values or plot their differences
# @param xlab Label for x-axis; defaults to empty string
# @param ylab Label for y-axis; defaults to NULL
# @param display_ylim Plot limits for y-axis; defaults to NULL
# @param main Plot title; defaults to empty string
plot_conditional_mean_quantiles <- function(samples, names, obs_xs,
                                            bin_min=NULL, bin_max=NULL,
                                            bin_delta=NULL,
                                            baseline_values=NULL,
                                            baseline_col="black",
                                            residual=FALSE,
                                            xlab="",
                                            ylab=NULL, display_ylim=NULL,
                                            main="") {
  # Check dimensions
  check_dimensions(obs_xs, 'obs_xs', names, 'names')

  # Check that baseline values are well-defined
  if (!is.null(baseline_values)) {
    if (length(baseline_values) != length(names)) {
      warning(paste0('The list of baseline values has the wrong',
                     ' dimension.  Baselines will not be plotted.'))
      baseline_values <- NULL
    }
  }

  # Check that names are in samples
  names <- check_expectand_names(names, samples)

  # Construct binning configuration
  bin_config <- configure_bins(bin_min, bin_max, bin_delta, obs_xs)
  bin_min <- bin_config[1]
  bin_max <- bin_config[2]
  bin_delta <- bin_config[3]

  # Construct bins
  breaks <- seq(bin_min, bin_max, bin_delta)
  plot_config <- configure_bin_plotting(breaks)
  plot_idxs <- plot_config[[1]]
  plot_xs <- plot_config[[2]]

  # Check bin containment
  check_bin_containment(bin_min, bin_max, obs_xs,
                        "conditioning value")

  # Construct quantiles for each conditional mean
  B <- length(breaks) - 1
  nonempty_bins <- sapply(1:B,
                          function(b)
                          sum(breaks[b] <= obs_xs & obs_xs < breaks[b + 1]) > 0)

  baseline_cond_means <- rep(NA, B)

  cond_mean <- function(y, x, b_low, b_high) {
    bin_idx <- which(b_low <= x & x < b_high)
    if (length(bin_idx))
      return(mean(y[bin_idx]))
    else
      return(mean(y))
  }

  if (!is.null(baseline_values)) {
    for (b in 1:B) {
      baseline_cond_means[b] <- cond_mean(baseline_values, obs_xs,
                                          breaks[b], breaks[b + 1])
    }
  }

  if (is.null(baseline_values) | !residual) {
    expectands <- list()
    for (b in 1:B) {
      expectands[[b]] <-
        local({
          x <- obs_xs;
          b_low <- breaks[b];
          b_high <- breaks[b + 1];
          function(y) cond_mean(y, x, b_low, b_high)
        })
    }
  } else {
    expectands <- list()
    for (b in 1:B) {
      expectands[[b]] <-
        local({
          x <- obs_xs;
          b_low <- breaks[b];
          b_high <- breaks[b + 1];
          baseline <- baseline_cond_means[b];
          function(y) cond_mean(y, x, b_low, b_high) - baseline
        })
    }
  }

  cond_mean_samples <-
    util$eval_expectand_pushforwards(samples,
                                     expectands,
                                     list('y'=array(names)))

  probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  mean_quantiles <- sapply(cond_mean_samples,
                           function(s)
                           util$ensemble_mcmc_quantile_est(s, probs))

  plot_quantiles <- do.call(cbind, lapply(plot_idxs,
                                          function(n) mean_quantiles[1:9, n]))

  # Plot
  if (is.null(display_ylim)) {
    if (is.null(baseline_values)) {
      display_ylim <- c(min(mean_quantiles[1, nonempty_bins]),
                        max(mean_quantiles[9, nonempty_bins]))
    }
    else {
      if (residual) {
        display_ylim <- c(min(c(0, mean_quantiles[1, nonempty_bins])),
                          max(c(0, mean_quantiles[9, nonempty_bins])))
      } else {
        display_ylim <- c(min(min(mean_quantiles[1, nonempty_bins]),
                              min(baseline_cond_means[nonempty_bins])),
                          max(max(mean_quantiles[9, nonempty_bins]),
                              max(baseline_cond_means[nonempty_bins])))
      }
    }
    delta <- 0.05 * (display_ylim[2] - display_ylim[1])
    display_ylim[1] <- display_ylim[1] - delta
    display_ylim[2] <- display_ylim[2] + delta
  }

  if (is.null(ylab)) {
    if (is.null(baseline_values) | !residual)
      ylab <- "Marginal Quantiles of Conditional Means"
    else
      ylab <- "Marginal Quantiles of Conditional Means - Baselines"
  }

  plot(1, type="n", main=main,
       xlim=c(bin_min, bin_max), xlab=xlab,
       ylim=display_ylim, ylab=ylab)

  polygon(c(plot_xs, rev(plot_xs)),
          c(plot_quantiles[1,], rev(plot_quantiles[9,])),
          col = c_light, border = NA)
  polygon(c(plot_xs, rev(plot_xs)),
          c(plot_quantiles[2,], rev(plot_quantiles[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(plot_xs, rev(plot_xs)),
          c(plot_quantiles[3,], rev(plot_quantiles[7,])),
          col = c_mid, border = NA)
  polygon(c(plot_xs, rev(plot_xs)),
          c(plot_quantiles[4,], rev(plot_quantiles[6,])),
          col = c_mid_highlight, border = NA)
  for (b in 1:B) {
    idx1 <- 2 * b - 1
    idx2 <- 2 * b
    if (nonempty_bins[b]) {
      lines(plot_xs[idx1:idx2], plot_quantiles[5, idx1:idx2],
            col=c_dark, lwd=2)
    } else {
      polygon(c(plot_xs[idx1:idx2], rev(plot_xs[idx1:idx2])),
              c(display_ylim[1], display_ylim[1],
                display_ylim[2], display_ylim[2]),
              col="#EEEEEE", border = NA)
    }
  }

  if (!is.null(baseline_values)) {
    if (residual) {
      abline(h=0, col="#DDDDDD", lwd=2, lty=3)
    } else {
      for (b in 1:B) {
        if (nonempty_bins[b]) {
          idx1 <- 2 * b - 1
          idx2 <- 2 * b
          lines(plot_xs[idx1:idx2], rep(baseline_cond_means[b], 2),
                col="white", lwd=4)
          lines(plot_xs[idx1:idx2], rep(baseline_cond_means[b], 2),
                col=baseline_col, lwd=2)
        }
      }
    }
  }
}


# Overlay nested quantile intervals to visualize an ensemble of empirical
# medians within the computed bins.
# Individual quantiles are estimated as the average of the empirical quantiles
# across each Markov chain, a consistent quantile estimator for Markov chain
# Monte Carlo.
# @param samples A named list of two-dimensional arrays for
#                each expectand.  The first dimension of each element
#                indexes the Markov chains and the second dimension
#                indexes the sequential states within each Markov chain.
# @param names List of relevant variable names
# @param obs_xs One-dimensional array of observed x-values on which to condition
# @param bin_min Lower threshold for conditioning
# @param bin_max Upper threshold for conditioning
# @param bin_delta Bin width for conditioning
# @param baseline_values Baseline values; defaults to NULL
# @param baseline_col Color for plotting baseline value; defaults to "black"
# @params residual Boolean value indicating whether to overlay quantiles and
#                  baseline values or plot their differences
# @param xlab Label for x-axis; defaults to empty string
# @param ylab Label for y-axis; defaults to NULL
# @param display_ylim Plot limits for y-axis; defaults to NULL
# @param main Plot title; defaults to empty string
plot_conditional_median_quantiles <- function(samples, names, obs_xs,
                                              bin_min=NULL, bin_max=NULL,
                                              bin_delta=NULL,
                                              baseline_values=NULL,
                                              baseline_col="black",
                                              residual=FALSE,
                                              xlab="",
                                              ylab=NULL, display_ylim=NULL,
                                              main="") {
  # Check dimensions
  check_dimensions(obs_xs, 'obs_xs', names, 'names')

  # Check that baseline values are well-defined
  if (!is.null(baseline_values)) {
    if (length(baseline_values) != length(names)) {
      warning(paste0('The list of baseline values has the wrong',
                     ' dimension.  Baselines will not be plotted.'))
      baseline_values <- NULL
    }
  }

  # Check that names are in samples
  names <- check_expectand_names(names, samples)

  # Construct binning configuration
  bin_config <- configure_bins(bin_min, bin_max, bin_delta, obs_xs)
  bin_min <- bin_config[1]
  bin_max <- bin_config[2]
  bin_delta <- bin_config[3]

  # Construct bins
  breaks <- seq(bin_min, bin_max, bin_delta)
  plot_config <- configure_bin_plotting(breaks)
  plot_idxs <- plot_config[[1]]
  plot_xs <- plot_config[[2]]

  # Check bin containment
  check_bin_containment(bin_min, bin_max, obs_xs,
                        "conditioning value")

  # Construct quantiles for each conditional mean
  B <- length(breaks) - 1
  nonempty_bins <- sapply(1:B,
                          function(b)
                          sum(breaks[b] <= obs_xs & obs_xs < breaks[b + 1]) > 0)

  baseline_cond_medians <- rep(NA, B)

  cond_median <- function(y, x, b_low, b_high) {
    bin_idx <- which(b_low <= x & x < b_high)
    if (length(bin_idx))
      return(median(y[bin_idx]))
    else
      return(median(y))
  }

  if (!is.null(baseline_values)) {
    for (b in 1:B) {
      baseline_cond_medians[b] <- cond_median(baseline_values, obs_xs,
                                              breaks[b], breaks[b + 1])
    }
  }

  if (is.null(baseline_values) | !residual) {
    expectands <- list()
    for (b in 1:B) {
      expectands[[b]] <-
        local({
          x <- obs_xs;
          b_low <- breaks[b];
          b_high <- breaks[b + 1];
          function(y) cond_median(y, x, b_low, b_high)
        })
    }
  } else {
    expectands <- list()
    for (b in 1:B) {
      expectands[[b]] <-
        local({
          x <- obs_xs;
          b_low <- breaks[b];
          b_high <- breaks[b + 1];
          baseline <- baseline_cond_medians[b];
          function(y) cond_median(y, x, b_low, b_high) - baseline
        })
    }
  }

  cond_median_samples <-
    util$eval_expectand_pushforwards(samples,
                                     expectands,
                                     list('y'=array(names)))

  probs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  median_quantiles <- sapply(cond_median_samples,
                             function(s)
                             util$ensemble_mcmc_quantile_est(s, probs))

  plot_quantiles <- do.call(cbind, lapply(plot_idxs,
                                          function(n) median_quantiles[1:9, n]))

  # Plot
  if (is.null(display_ylim)) {
    if (is.null(baseline_values)) {
      display_ylim <- c(min(median_quantiles[1, nonempty_bins]),
                        max(median_quantiles[9, nonempty_bins]))
    }
    else {
      if (residual) {
        display_ylim <- c(min(c(0, median_quantiles[1, nonempty_bins])),
                          max(c(0, median_quantiles[9, nonempty_bins])))
      } else {
        display_ylim <- c(min(min(median_quantiles[1, nonempty_bins]),
                              min(baseline_cond_medians[nonempty_bins])),
                          max(max(median_quantiles[9, nonempty_bins]),
                              max(baseline_cond_medians[nonempty_bins])))
      }
    }
    delta <- 0.05 * (display_ylim[2] - display_ylim[1])
    display_ylim[1] <- display_ylim[1] - delta
    display_ylim[2] <- display_ylim[2] + delta
  }

  if (is.null(ylab)) {
    if (is.null(baseline_values) | !residual)
      ylab <- "Marginal Quantiles of Conditional Medians"
    else
      ylab <- "Marginal Quantiles of Conditional Medians - Baselines"
  }

  plot(1, type="n", main=main,
       xlim=c(bin_min, bin_max), xlab=xlab,
       ylim=display_ylim, ylab=ylab)

  polygon(c(plot_xs, rev(plot_xs)),
          c(plot_quantiles[1,], rev(plot_quantiles[9,])),
          col = c_light, border = NA)
  polygon(c(plot_xs, rev(plot_xs)),
          c(plot_quantiles[2,], rev(plot_quantiles[8,])),
          col = c_light_highlight, border = NA)
  polygon(c(plot_xs, rev(plot_xs)),
          c(plot_quantiles[3,], rev(plot_quantiles[7,])),
          col = c_mid, border = NA)
  polygon(c(plot_xs, rev(plot_xs)),
          c(plot_quantiles[4,], rev(plot_quantiles[6,])),
          col = c_mid_highlight, border = NA)
  for (b in 1:B) {
    idx1 <- 2 * b - 1
    idx2 <- 2 * b
    if (nonempty_bins[b]) {
      lines(plot_xs[idx1:idx2], plot_quantiles[5, idx1:idx2],
            col=c_dark, lwd=2)
    } else {
      polygon(c(plot_xs[idx1:idx2], rev(plot_xs[idx1:idx2])),
              c(display_ylim[1], display_ylim[1],
                display_ylim[2], display_ylim[2]),
              col="#EEEEEE", border = NA)
    }
  }

  if (!is.null(baseline_values)) {
    if (residual) {
      abline(h=0, col="#DDDDDD", lwd=2, lty=3)
    } else {
      for (b in 1:B) {
        if (nonempty_bins[b]) {
          idx1 <- 2 * b - 1
          idx2 <- 2 * b
          lines(plot_xs[idx1:idx2], rep(baseline_cond_medians[b], 2),
                col="white", lwd=4)
          lines(plot_xs[idx1:idx2], rep(baseline_cond_medians[b], 2),
                col=baseline_col, lwd=2)
        }
      }
    }
  }
}
