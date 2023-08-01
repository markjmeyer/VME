# VME
 Code to implement variational AICs on parameter-efficient variational mixed effects models is in the file vme.R. A description of the method is available on arXiv, [Meyer, Carter, & Malloy (2023)](https://arxiv.org/abs/2206.04012). Functions are described below. Folders contain code and data for the simulations and illustrations, additional details below and in the manuscript.

## Functions
Functions in the file vme.R include vme() and VAIC().

 ### parameter-efficient VME
 Build a parameter-efficient VME with general random effects structure using
```
vme(fixed, Z = NULL, rint.only = TRUE, fe = NULL, data, ID,
    ps.spline = list(df = 5, degree = 3, xi = 0.001),
    re.spline = list(df = NULL, degree = 3, xi = 0.001),
    priors = list(Ae = 0.01, Be = 0.01, Af = 0.01, Bf = 0.01, 
                  Au = 0.01, Bu = 0.01, sigB = 1000),
    tol = 1e-5, maxIter = 500, up = 100, dots = up/10,
    verbose = TRUE)
```
The argument fixed should be a formula containing at least the outcome and an intercept but can include additional fixed effects. For the random intercept only model, Z is left empty and rint.only = TRUE. For a general random effects structure, set rint.only = FALSE. User must provide a matrix for Z. Examples of possible random effect design matrices are in data_illustrations.R file found in the Illustrations folder. Functions to build these matrices are forth-coming.

Additional arguments detail the controls for random effects spline (re.spline) and any functional effets (ps.spline). Functional effects are specified through the argument fe as a one-sided formula. The priors arugment controls the priors for the model based on [Meyer, Carter, & Malloy (2023)](https://arxiv.org/abs/2206.04012). The remaining arguments include the tolerance (tol) for determining convergence and the total number of iterations to run (maxIter)--this may need to be increased if convergence is not achieved within maxIter iterations. The arugments up, dots, and verbose control console output for monitoring the model.

Returns an object of class vme. Class-based functions including summary(), coef(), and plot() can be called on objects of class vme.

 ### variational AIC
 Extract the VAIC of a VME using
```
VAIC(model)
```
 Takes an object, model, of class vme. VAIC is calculated via the function VAIC_vme() when vme() is called.

### Dependencies
Requires the following R packages: [Matrix](https://cran.r-project.org/web/packages/Matrix/index.html), splines (should be included in standard R release), [MASS](https://cran.r-project.org/web/packages/MASS/index.html), and [numbers](https://cran.r-project.org/web/packages/numbers/index.html).

## Simulations folder
This folder contains files that were used for the descrimination, bias, and MSE/MISE simulations described in [Meyer, Carter, & Malloy (2023)](https://arxiv.org/abs/2206.04012). The files sim_cycl_func.R, sim_peak_func.R, and sim_sigm_func.R run the random functions simulation while sim_risb.R runs the random intercept vs random slope simulation. To evaluate the results of the simulation and to replicate the tables in the manuscript, use the files functional_eval.R and insterslope_eval.R for the random function and random intercept vs random slope simulations, respectively.

## Illustrations folder
This folder contains a script and dataset for the data illustrations described in [Meyer, Carter, & Malloy (2023)](https://arxiv.org/abs/2206.04012). The script is titled data_illustrations.R. Data for the random intercept vs random slope illustration is in the accompanying file leadlong.txt. Data for the random function illustration is available in the [refund](https://cran.r-project.org/web/packages/refund/index.html) R package.
