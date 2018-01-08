# CaliCo (CALIbration COde): package for Bayesian calibration

## Description

This package provides a set of methods which allow to create a statistical model with the associated priors and the calibration. Once the calibration terminated, the created object is used for predicting the behavior of the wanted code over other data.

```
devtools::install_github("mathieucarmassi/CaliCo")
```

## Use and example

See the package [vignette](https://github.com/mathieucarmassi/CaliCo/blob/master/vignettes/introduction.Rmd). To build the vignettes on intallation, you need the *ade4* package installed and to run the following code:

```
devtools::install_github("mathieucarmassi/CaliCo", build_vignettes=TRUE)
```
