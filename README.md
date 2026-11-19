# FlowPredict

Predictive modelling of per-sample cytometry features (population frequencies, marker medians or any other numbers), in R:

- **Group prediction**: which small set of features best separates two groups, and how well? Lasso logistic regression and random forest, judged by repeated nested cross-validation (ROC AUC, balanced accuracy), with feature importance and the minimal panel of top-ranked features.
- **Immune age model**: how well do the features predict a number such as age? Elastic net with nested cross-validation, age gap corrected for regression to the mean, and a permutation check against chance.

Folds can keep each donor's samples together. Covariates can be added to every model (unpenalised).

The method follows the nested cross-validation design of Veiga et al. 2026 (*Nat Commun* 17:4670), Humblet-Baron et al. 2025 (*Alzheimer's Dement* 21:e70952) and Gemander et al. 2025 (*npj Vaccines* 10:140), and is the R counterpart of the predictive analysis and immune age model in [Honeychrome](https://github.com/salmansamson/honeychrome).

## Install

```r
install.packages( "remotes" )
remotes::install_github( "DrCytometer/FlowPredict", dependencies = TRUE )
```

Or source `inst/scripts/install_FlowPredict.R`, which also checks the installation.

## Use

```r
library( FlowPredict )
data <- read.feature.data( "features.csv", "metadata.csv", sample.column = "Samples",
                           transform = "log2" )
groups <- run.predictive( data, group.column = "Group", positive = "Case", negative = "Control" )
write.predictive.results( groups, "results/Case_vs_Control" )

age <- run.age.model( data, value.column = "Age", group.column = "Group",
                      model.groups = "Control" )
write.age.results( age, "results/immune_age", unit = "years" )
```

`inst/scripts/run_FlowPredict.R` is a complete, editable analysis script.

## Input

- **Feature CSV**: a `Samples` column and one numeric column per feature.
- **Metadata CSV**: the same `Samples` column and columns for groups, donors and covariates.

Use `transform = "log2"` for frequencies or counts (`log2(x + 0.01)`, as Honeychrome does), `"none"` for marker medians.

## Sample size

Group prediction needs at least 10 samples per group (30 recommended); the immune age model at least 20 samples with a value (40 recommended).
