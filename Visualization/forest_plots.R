time_varying_forest_plot <- function(
    data, x = "emmean", var, var_label,
    ll = "asymp.LCL", ul = "asymp.UCL",
    time_var = "tgroup", time_label = "Time Period",
    est_scale = "log", exponentiate = F) {
  require(ggplot2)
  require(forcats)
  # By default assumes data are on the log scale and will not exponentiate
  vert <- 0
  xlabel <- "log(Hazard)"
  if (est_scale == "exp" | (est_scale == "log" & exponentiate)) {
    vert <- 1
    xlabel <- "Hazard"
  }
  if (exponentiate) {
    data[, x] <- exp(data[, x])
    data[, ll] <- exp(data[, ll])
    data[, ul] <- exp(data[, ul])
  }
  # Base plot
  p <- ggplot(data, aes(
    x = !!sym(x), y = fct_rev(!!sym(time_var)),
    color = !!sym(time_var)
  )) +
    geom_point(size = 10) +
    geom_linerange(aes(xmin = !!sym(ll), xmax = !!sym(ul)), linewidth = 2) +
    theme_classic() +
    facet_grid(as.formula(paste0(var, "~."))) +
    ylab(time_label) +
    xlab(xlabel) +
    theme(legend.position = "none") +
    geom_vline(xintercept = vert, linetype = "dashed")
  return(p)
}
