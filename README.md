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

Results are saved in us_results.mat in the dynare project folder. oo_ contains most of the useful information. To extract the results, you may load the mat in Octave and the output them to a csv file. Or you can do copy and paste if the data volume is small.

load('C:\dynare\us\us_results.mat');
csvwrite('C:\dynare\us\oo_vd.csv', oo_.variance_decomposition);


# Multifactor Regression
Regression models are used to build the relationships between macroeconomic factors and asset returns. R programs are used to test five model types: linear regression, GLM, CART, KNN, ANN. It also includes prediction of economic recession, repairing correlation matrix for non-positive definite ones, and Vector Autoregressive Model for macroeconomic factors

Data: inputmap.csv for multifactor regression; varinput.csv for VAR macroecnomic model

Program: fundmapping.R

Change working directory in setwd("C:/dsge/r")

# ESG
ESG using DSGE and multifactor regression. It also includes yield curve interpolation and extropolation, and bond return calculation based on bond yield curve and investment strategy.

Data
mapping.csv: linear regression model parameters;
normalchol.csv: Cholesky decomposition for expansion periods
recessionchol.csv: Cholesky decomposition for recession periods
termmix.csv: target term mix for bond fund
migration.csv: credit rating migration matrix, default rate, and recovery rate

Program: ESG.R

