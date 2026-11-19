# run_FlowPredict.R
#
# Predictive modelling of a feature table with FlowPredict:
#   1. group prediction (lasso logistic regression and random forest) for
#      each pair of groups listed in 'comparisons';
#   2. optionally, an immune age model for a numeric metadata column.
#
# Inputs
#   features.file : CSV, one row per sample: a 'Samples' column and any
#                   number of numeric feature columns.
#   metadata.file : CSV, one row per sample: the same 'Samples' column and
#                   the group, donor and covariate columns.
#
# Edit the settings below, then run the whole script (RStudio: "Source").
# Results (CSV tables, PDF and PNG figures) go to 'out.dir', one subfolder
# per analysis.


# ---- Settings ------------------------------------------------------------

# Input files and the name of their sample column
features.file <- "features.csv"
metadata.file <- "metadata.csv"
sample.column <- "Samples"

# "log2" for frequencies (%) or counts; "none" for marker medians or
# features that are already transformed or can be negative
transform <- "log2"

# Output folder
out.dir <- "FlowPredict_results"

# Group prediction: one line per comparison, giving the metadata column
# holding the groups, then the positive group, then the negative group.
# Set comparisons <- list() to skip group prediction.
comparisons <- list(
  c( "Group", "GroupA", "GroupB" )
)

# Optional: metadata column naming the donor of each sample, when a donor
# gave more than one sample (folds then keep a donor's samples together).
# NULL when every sample is from a different donor.
donor.column <- NULL

# Optional: metadata columns added to every model (e.g. c("Age", "Sex")).
# NULL for none.
covariates <- NULL

# Immune age model: TRUE to run it
run.age <- TRUE
age.column <- "Age"                 # numeric metadata column to predict
age.unit <- "years"                 # for axis labels; "" for none
age.group.column <- "Group"         # NULL to use every sample with a value
age.model.groups <- "GroupB"        # groups the model is trained on
age.permutations <- 100             # chance check shuffles (0 to skip)

# Cross-validation and reproducibility
outer.folds <- 5
repeats <- 5
inner.folds <- 3
seed <- 0


# ---- Analysis --------------------------------------------------------------

library( FlowPredict )

data <- read.feature.data( features.file, metadata.file, sample.column = sample.column,
                           transform = transform )
print( data )
dir.create( out.dir, showWarnings = FALSE, recursive = TRUE )

safe.name <- function( x ) gsub( "[^A-Za-z0-9]+", "_", x )

for ( comparison in comparisons ) {
  name <- paste0( safe.name( comparison[ 2 ] ), "_vs_", safe.name( comparison[ 3 ] ) )
  cat( "\n== Group prediction:", comparison[ 2 ], "vs", comparison[ 3 ], "\n" )
  result <- run.predictive( data, group.column = comparison[ 1 ],
                            positive = comparison[ 2 ], negative = comparison[ 3 ],
                            donor.column = donor.column, covariates = covariates,
                            outer.folds = outer.folds, repeats = repeats,
                            inner.folds = inner.folds, seed = seed )
  print( result )
  write.predictive.results( result, file.path( out.dir, name ), prefix = name )
  saveRDS( result, file.path( out.dir, name, paste0( name, "_result.rds" ) ) )
}

if ( run.age ) {
  name <- paste0( "immune_age_", safe.name( age.column ) )
  cat( "\n== Immune age model:", age.column, "\n" )
  age <- run.age.model( data, age.column, group.column = age.group.column,
                        model.groups = age.model.groups, donor.column = donor.column,
                        covariates = covariates, outer.folds = outer.folds,
                        repeats = repeats, inner.folds = inner.folds,
                        n.permutations = age.permutations, seed = seed )
  print( age )
  print( age$by.group )
  write.age.results( age, file.path( out.dir, name ), prefix = name, unit = age.unit )
  saveRDS( age, file.path( out.dir, name, paste0( name, "_result.rds" ) ) )
}

writeLines( capture.output( sessionInfo() ), file.path( out.dir, "sessionInfo.txt" ) )
cat( "\nDone. Results are in", normalizePath( out.dir ), "\n" )
