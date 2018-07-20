# CaliCo 0.1.1

## Major changes

 - To establish a model with a Gaussian process that emulates the code, the user has now three possibilities. The first, he possesses a code written in R but no DOE. The second, he possesses a code written in R and wants to enforce its own DOE. Then, he has no code but possesses the DOE and the corresponding output of a code.
 - A pipe is defined to parametrize a statistical model.
 - The function plot for parametrized model takes now only two arguments and is simpler than in version 0.0.1
 - A fonction forecast replaces the function prediction which did not work very very. The new function forecast allows to predict, based on realized calibration, over a new data set the behaviour of the model.
 - Two functions chain and estimators are added to get access to the MCMC chains and the MAP and Mean a posteriori estimator of a calibrated object.
 - The ggplots given for the plot of the calibrated object has changed and now to graphics layout are given. The user is also free to load the ggplots into a variable and draw the graph he wants.
 - A function that improves calibration for the second and the fourth model have also been written and is called sequentialDesign.
 
## Bug fixes

 - The likelihoods written in the models are fixed
 - The R6 classes to create the statistical model have been simplified to have a more robust programmation
