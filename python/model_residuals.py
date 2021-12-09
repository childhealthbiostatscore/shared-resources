#-------------------------------------------------------------------------------
#
# This file contains several useful Python functions for statistical analysis.
#
# v 1.0
# Tim Vigers
# Last updated May 23, 2021
#
#-------------------------------------------------------------------------------

# Plot model residuals from a statmodels fit
def plot_resid(model):
    res = model.resid
    ax = sns.kdeplot(mdf.resid, fill = True)
    ax.set_title("KDE Plot of Model Residuals (Blue) and Normal Distribution (Black)")
    ax.set_xlabel("Residuals")