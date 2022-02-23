AUTHOR: Kristen Miller (kristen.miller@cuanschutz.edu)

Repeated measures on the same patient, but not all pts have multiple observations,
and among those who do, not that many observations.
Tried GLMM, but models did not converge.
GEE converges and accounts for correlation.

-- gee r package
-- MuMIn r package for multivariable model selection
-- Wald type tests using contrast statements

function for glmm.R: functions that fit univariate gee, GLMM, logistic regression
02_univaraite: calls the functions for glmm.R, fits models and outputs tables
03_multivariable: fits multivariable gee models using model selection procedures
