# Functions for calculating glucose metabolism indices.
# These functions assume that data is provided in the wide format with
# consistently labelled columns of the form "{measure prefix}_{timepoint}"
# All MINMOD code is based on Allen Downey's book "Modeling and Simulation in
# Python" (https://allendowney.github.io/ModSimPy/index.html).

# The following function is used for numerically solving the differential equations.
# Gb and Ib are the glucose and insulin measurements (respectively) at t=0 from the data.
# k1, k2, and k3 are parameters in the equations and will be chosen based on iterative least squares.

from scipy.integrate import odeint


def vectorfield(w, t, p):
    """
    Defines the differential equations for the coupled spring-mass system.

    Arguments:
        w :  vector of the state variables:
                  w = [x1,y1,x2,y2]
        t :  time
        p :  vector of the parameters:
                  p = [m1,m2,k1,k2,L1,L2,b1,b2]
    """
    x1, y1, x2, y2 = w
    m1, m2, k1, k2, L1, L2, b1, b2 = p

    # Create f = (x1',y1',x2',y2'):
    f = [y1,
         (-b1 * y1 - k1 * (x1 - L1) + k2 * (x2 - x1 - L2)) / m1,
         y2,
         (-b2 * y2 - k2 * (x2 - x1 - L2)) / m2]
    return f


def minmod(s, t, p):
    """
    Defines the differential equations.
    Arguments:
        s :  vector of the state variables
        t :  time
        p :  vector of the parameters
    """
    G, X, I, Gb, Ib = s
    k1, k2, k3, = p
    dGdt = -k1 * (G - Gb) - X * G
    dXdt = k3 * (I(t) - Ib) - k2 * X
    return dGdt, dXdt


# Call the ODE solver.
wsol = odeint(minmod, s, t, args=(p,))


# Testing
import pandas as pd
import numpy as np
df = pd.read_csv("~/GitHub/ModSimPy/data/glucose_insulin.csv")


def solve_min_mod(data=df, k1=0.02, k2=0.02, k3=1.5e-05, dt=1):
    from scipy.integrate import solve_ivp
    from scipy.interpolate import interp1d
    # Pull starting parameters from data
    t_0 = data["time"].iloc[0]
    t_end = data["time"].iloc[-1]
    Gb = data["glucose"].iloc[0]
    Ib = data["insulin"].iloc[0]
    I = interp1d(data["time"], data["insulin"])
    ts = np.array([i for i in range(t_0, t_end, dt)])
    # Solve
    b = solve_ivp(fun=minmod, t_span=[t_0, t_end], t_eval=ts,
                  y0=[I, Ib, Gb],
                  args=[k1, k2, k3])


df["model glucose"] = pd.Series(b.y[0])
df["model insulin"] = pd.Series(b.y[1])

# Plot
from plotnine import ggplot, aes, geom_point, geom_line
ggplot(df) + aes(x="time", y="glucose") + \
    geom_point() + geom_line(aes(y="model glucose"))

# Optimize parameters by least squares


def error_func(initial_state, k1, k2, k3, data=df):
    # Solve ODI
    b = solve_ivp(fun=minmod, t_span=[t_0, t_end], t_eval=ts,
                  y0=[initial_state, I, Ib, Gb],
                  args=[initial_state, k1, k2, k3])
    # Calculate error
    data["model glucose"] = pd.Series(b.y[0])
    # compute the difference between the model
    # results and actual data
    errors = data["model glucose"] - data["glucose"]
    return errors


import scipy

t = scipy.optimize.leastsq(error_func, x0=[initial_state, k1, k2, k3],
                           args=[I, Ib, Gb])
