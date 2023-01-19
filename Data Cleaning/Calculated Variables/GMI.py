# Functions for calculating glucose metabolism indices.
# These functions assume that data is provided in the wide format with
# consistently labelled columns of the form "{measure prefix}_{timepoint}"
# All MINMOD code is based on Allen Downey's book "Modeling and Simulation in
# Python" (https://allendowney.github.io/ModSimPy/index.html).

# For testing
import pandas as pd
df = pd.read_csv("/Volumes/PEDS/RI Biostatistics Core/Shared/Shared Projects/" +
                 "Laura/Peds Endo/Petter Bjornstad/Data Harmonization/Data Clean/harmonized_dataset.csv",
                 low_memory=False)
df = pd.read_csv(
    "/Users/timvigers/Downloads/ModSimPy/chapters/glucose_insulin.csv")

# The following function is used for numerically solving the differential equations.
# Gb and Ib are the glucose and insulin measurements (respectively) at t=0 from the data.
# k1, k2, and k3 are parameters in the equations and will be chosen based on iterative least squares.


def minmod(t, initial_state, k1, k2, k3, I, Ib, Gb):
    G, X = initial_state
    dG = -k1 * (G - Gb) - X * G
    dX = k3 * (I(t) - Ib) - k2 * X
    return dG, dX


# Testing
import pandas as pd
import numpy as np
df = pd.read_csv(
    "/Users/timvigers/Dropbox/Work/ModSimPy/chapters/glucose_insulin.csv")

from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
initial_state = pd.Series({"G": 270, "X": 0}, name='state')
t_0 = df["time"].iloc[0]
t_end = df["time"].iloc[-1]
Gb = 92
Ib = 11
I = interp1d(df.time, df.insulin)
k1 = 0.02
k2 = 0.02
k3 = 1.5e-05
dt = 2
ts = np.array([i for i in range(t_0, t_end, dt)])


b = solve_ivp(fun=minmod, t_span=[t_0, t_end], t_eval=ts,
              y0=initial_state,
              args=[k1, k2, k3, I, Ib, Gb])

df["model glucose"] = pd.Series(b.y[0])
df["model insulin"] = pd.Series(b.y[1])

# Plot
from plotnine import ggplot, aes, geom_point, geom_line
ggplot(df) + aes(x="time", y="glucose") + \
    geom_point() + geom_line(aes(y="model glucose"))


def estimate_min_mod(df, glucose_prefix="glucose", insulin_prefix="insulin",
                     pre_baseline_prefix="minus", id_cols=["record_id", "co_enroll_id", "study", "visit", "procedure", "date"]):
    import pandas as pd
    from natsort import natsorted, ns
    from scipy.interpolate import interp1d
    # Find the glucose and insulin columns, sort them
    gluc_cols = df.filter(
        regex=(glucose_prefix + "_\\d{1,3}$")).columns.to_list()
    gluc_cols = natsorted(gluc_cols, alg=ns.IGNORECASE)
    neg_gluc = df.filter(regex=(
        glucose_prefix + "_" + pre_baseline_prefix + "_\\d{1,3}$")).columns.to_list()
    neg_gluc = natsorted(neg_gluc, alg=ns.IGNORECASE, reverse=True)
    gluc_cols = neg_gluc + gluc_cols

    ins_cols = df.filter(
        regex=(insulin_prefix + "_\\d{1,3}$")).columns.to_list()
    ins_cols = natsorted(ins_cols, alg=ns.IGNORECASE)
    neg_ins = df.filter(regex=(
        insulin_prefix + "_" + pre_baseline_prefix + "_\\d{1,3}$")).columns.to_list()
    neg_ins = natsorted(neg_ins, alg=ns.IGNORECASE, reverse=True)
    ins_cols = neg_ins + ins_cols
    # Convert to long
    gluc = df[id_cols + gluc_cols].copy()
    ins = df[id_cols + ins_cols].copy()
    # Here is the fun part! Solve ODEs
