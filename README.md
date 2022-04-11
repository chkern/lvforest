# SC Forest for R

R package for "Survey Scale Forests: Estimating Valid Latent Variable Scores from Conditionally Causal Models"

We develop a survey scale forest (SC Forest) algorithm for the estimation of latent variable scores from conditionally causal models with one or more latent variables. SC Forest establishes conditional causality in confirmatory factor analysis (CFA) models with ordinal and/or numerical response variables. Through parametric model restrictions paired with a non-parametric tree-based machine learning approach, SC Forest estimates latent variables scores that fulfill the main criteria for construct validity.

### Installation

``` {.r}
if (!require("devtools")) install.packages("devtools")
devtools::install_github("chkern/scforest")
```

### Usage

``` {.r}
trained <- scforest.train(input = c("conf","dicho1","cat1","num1","cat2","rand1","rand2","rand3","rand4","rand5"), 
model <- 'LatVar1 =~ simuvar1 + beta2*simuvar2 + beta3*simuvar3 
LatVar2 =~ simuvar4 + beta5*simuvar5 + beta6*simuvar6
LatVar3 =~ simuvar7 + beta8*simuvar8 + beta9*simuvar9
LatVar4 =~ LatVar1 + lambda2*LatVar2 + lambda3*LatVar3 + delta*simuvar_effect',
ordered = c("simuvar1","simuvar2","simuvar3","simuvar4","simuvar5","simuvar6","simuvar7","simuvar8","simuvar9"), 
data = simu,
cutoff_rmsea = .03,
ctree_control = ctree_control(minbucket = 200, mtry = 1, testtype = "Teststatistic"))
```

``` {.r}
predicted <- scforest.predict(trained = trained, data = simu, idvar = "id", latvar = "LatVar4", exclude_unconf = T, stdscores = F)
```
