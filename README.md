[![Travis-CI Build Status](https://travis-ci.org/veneficusnl/pmpp.svg?branch=master)](https://travis-ci.org/veneficusnl/pmpp)

# What is the Posterior Mean Panel Predictor?

## Accurate predictions with micro-panels

Micro-panels are longitudinal data sets that contain observations on multiple units at only a few points in time. Examples include the performance of start-up companies, developmental skills of small children or revenues and leverage of banks after significant regulatory changes. When working with micro-panels, it is challenging to build accurate predictive models, as the time series are too short to contain enough information on their own. 

Posterior Mean Panel Predictor (PMPP) takes an empirical-Bayes approach to computing forecasts with micro-panels. It uses cross-sectional information in the data to approximate the posterior mean of heterogeneous coefficients under a correlated random effects distribution. It has been shown to provide predictions of higher accuracy compared to the state-of-the-art methods for dynamic panel modelling. For more details, see the references in `pmpp()` function manual.

## Package features

The package allows for the following:

* Estimate the parameters of the PMPP model,
* Use the model to compute point forecasts,
* Compute prediction intervals with the Random-Window Block Bootstrap.

Additionally, the package exports a number of functions that can be used outside of the scope of PMPP modelling:

* `kde()` for computing a robust kernel density estimate;
* `kde2D()` for computing a robust 2-dimensional kernel density estimate;
* `create_fframe()` for adding time periods to a panel-structured data frame;
* `ssys_gmm()`, the suboptimal multi-step System-GMM estimator for AR(1) panel data model.

# How to use

The central function in the package is `pmpp()`. It estimates the model's coefficients and outputs an object of class `pmpp`. This class has the `plot` and `summary` methods, with the former plotting the distribution of individual-specific effects and the latter allowing to inspect model's coeffcients and fit measures. 

To compute predictions with the PMPP model, one needs to construct the forecast frame with `create_fframe()`. The forecast frame and the corresponding model object can be passed along to the `predict` method to obtain forecasts.

In order to calculate prediction intervals, the `pmpp_predinterval()` function can be used. This function, similarly to the `predict` method, takes the model object and the forecast frame as inputs. Be warned: bootstrapping of prediction interval might take time!

# Usage example

```
# Get data
data(EmplUK, package = "plm")
EmplUK <- dplyr::filter(EmplUK, year %in% c(1978, 1979, 1980, 1981, 1982))

# Run the model predicting employment
pmpp_model <- pmpp(dep_var = "emp", data = EmplUK)
summary(pmpp_model)

# Compute predictions for following three years
my_fframe <- create_fframe(EmplUK, 1983:1985)
prediction <- predict(pmpp_model, my_fframe)

# Compute prediction intervals
intervals <- pmpp_predinterval(pmpp_model, my_fframe, bootReps = 20, confidence = 0.95)
```
