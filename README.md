ts2sls
======

Two-sample two-stage least squares estimation

Currently, this is just a Stata ado file to calculate two-sample 2sls a la Inoue and Solon (2010).

For the future:
1. More numerically stable matrix calculations
2. Formatted output
3. Standard errors are not fully corrected in the case of multiple instruments
4. Perfect collinearity checking is broken

### Description 

ts2sls.ado
J. SHRADER
jgshrade@ucsd.edu
First Version: 2014/10/30

Stata program to calculate two-sample two-stage least squares (TS2SLS) estimate . Math is based on Inoue and Solon (2005), although variable names more closely follow the shorter version published as Inoue and Solon (2010).

Some of the code is taken more or less directly from Solomon Hsiang's ols_spatial_HAC.ado file.

If you find errors, please let me know.

#### Syntax

ts2sls y (x = z) [if] [in], group(group_var) [noconstant]

Where y is the outcome variable, x be the endogenous regressor, and z an exogenous instrument. I follow the notation of Inoue and Solon and call the data for estimating the reduced form (y as a function of z) "group 1" and the data for estimating the first stage (x as a function of z) "group 2".

Your data needs to be stacked, so if you are estimating from two different datasets, append them, then specify "group" appropriately. The dataset should look like this:
```
  Group  |    Y    |    X    |    Z   
---------+---------+---------+---------
    1    |    y    |    .    |    z_1
    2    |    .    |    x    |    z_2
```

### References 

Angrist, Joshua D., and Alan B. Krueger. "The effect of age at school entry on educational attainment: an application of instrumental variables with moments from two samples." Journal of the American Statistical Association 87, no. 418 (1992): 328-336.

Inoue, Atsushi, and Gary Solon. "Two-sample instrumental variables estimators." The Review of Economics and Statistics 92, no. 3 (2010): 557-561.

Inoue, Atsushi, and Gary Solon. "Two-Sample Instrumental Variables Estimators." NBER Working Paper (2005).
