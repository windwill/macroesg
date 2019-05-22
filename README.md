# macroesg
Macroeconomics based ESG

It includes calibration of a DSGE model, multifactor regression, and economic scenario generation.

# DSGE Model
The DSGE model is built using dynare and tested using Octave, an open source software that replicates Matlab.

The model is documented according to a research paper (to be updated later)

Dynare version 4.5.7 (https://www.dynare.org/download/)

Octave version 4.4.1 (https://ftp.gnu.org/gnu/octave/windows/octave-4.4.1-w64-installer.exe)

Data: us_data.m

Model: usb.mod

To run the model:

1. addpath c:\dynare\4.5.7\matlab
2. cd C:\dynare\us
3. dynare us

Folder "us" needs to include data file and model file

# Multifactor Regression
Regression models are used to build the relationships between macroeconomic factors and asset returns. R programs are used to test five model types: linear regression, GLM, CART, KNN, ANN. It also includes prediction of economic recession and repair correlation matrix for non-positive definite ones.

Data: inputmap.csv

Program: fundmapping.R

Check working directory in setwd("C:/dsge/r")

# ESG
