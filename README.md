# ldr
Estimates the best linear dynamic range for dilution experiments.

## Motivation
There seems to be an overall lack of a decent, single step, automatic, determination of the linear dymanic range (LDR) for a chemical or molecular dilution series. This approach started out using the accepted standards for determination of LDR, except the case study this package was developed for did not conform to the standard S shaped dilution profile, and instead was affected by and mild rise in observed instruments response prior to flattening on the lower limit of detection. Therefor, an approach was taken to calculate all possible linear slopes and selected the best based on largest contiguous range given the fit passed the cut-offs for precision, accuracy, linearity and number of points rejected.

## Usage
```
>library(ldr)
>data.frame.fit <- linearFit(data.frame.org)

>str(data.frame.fit)

'data.frame':	N obs. of  12+ variables:
 $ fit_points      : int  4 3 7 11 4 7 11 3 6 6 ...
 $ fit_slope       : num  0.67504 0.64137 0.01474 -1.0613 0.00523 ...
 $ fit_intercept   : num  6.1 6.6 4.28 5.77 3.23 ...
 $ fit_rsquared    : num  0.9678 0.9636 0.0506 0.5762 0.0115 ...
 $ fit_cv_mean     : num  0.0623 0.1021 0.0618 0.1721 0.1719 ...
 $ fit_res_mean    : num  0.00965 0.00986 0.01225 0.0513 0.01751 ...
 $ fit_dil_min     : num  3.16e-01 1.00e-04 3.16e-05 1.00e-02 3.16e-03 1.00e-06 1.00e-03 3.16e-04 1.00e-01 3.00e-02 ...
 $ fit_dil_max     : num  1.00e+01 1.00e-03 1.00e-05 1.00e-01 1.00e-03 1.00e-03 3.16e-06 3.16e-03 3.16e+01 1.00e+01 ...
 $ fit_lloq        : int  642837 11426 13537 855382 1602 8684 28502 17280 27814 13670 ...
 $ fit_uloq        : int  5810478 50229 20813 22455915 1713 22025 36442 42639 1638195 859143 ...
 $ fit_fail        : logi  FALSE FALSE TRUE TRUE TRUE FALSE ...
 $ fit_fail_cat    : Factor w/ 17 levels "","3 or less conc",..: 1 1 4 2 4 1 4 4 1 1 ...

```

## Arguments
```
data      data.frame of dilution~response values
col_dil   column containing the variable for dilution concentration
col_res   column containing the variable for dilution response
col_uid   column containing the variable for unique id
rjct_pre  rejection criteria for min precision, default 0.25,
rjct_acc  rejection criteria for min accuracy, default 0.25,
rjct_rsq  rejection criteria for min linearity Pearson corrleation, default 0.90,
rjct_pts  rejection criteria for dilution points retained, default 0.75,
filter    filter dataset to top solution passing reject criteria, default TRUE
```

# Installation 
How do you install a package that’s sitting on GitHub?
First, you need to install the devtools package. You can do this from CRAN. Invoke R and then type
```
install.packages("devtools")
```
Load the devtools package.
```
library(devtools)
```
In most cases, you just use install_github("author/package"). 
```
install_github("jeffsocal/ldr")
```
