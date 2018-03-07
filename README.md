# CaliCo (CALIbration COde): package for Bayesian calibration

## Description

This package provides a set of methods which allow to create a statistical model with the associated priors and then run the calibration on the wished parameters. Once the calibration terminated, the created object is used for predicting the behavior of the statistical model with other data.

```
devtools::install_github("mathieucarmassi/CaliCo")
```

## Use and example

See the [vignette](https://github.com/mathieucarmassi/CaliCo/blob/master/vignettes/CaliCo-introduction.Rmd). To build the vignettes on intallation, you need to run the following code:

```
devtools::install_github("mathieucarmassi/CaliCo", build_vignettes=TRUE)
```
