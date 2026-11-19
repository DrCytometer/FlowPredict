# install_FlowPredict.R
#
# Installs FlowPredict and its dependencies (glmnet, ranger, ggplot2) from
# GitHub and CRAN, then runs a short check on simulated data.
#
# Run once, in RStudio (open this file and click "Source") or from a
# terminal:  Rscript install_FlowPredict.R
#
# Needs R 4.1 or later and an internet connection.  No compiler is needed
# (FlowPredict is written in plain R; its dependencies install as CRAN
# binaries on Windows and macOS).

github.repo <- "DrCytometer/FlowPredict"

# To install from a downloaded copy instead of GitHub, put the path of the
# package folder or .tar.gz file here, e.g. "~/Downloads/FlowPredict_0.1.0.tar.gz"
local.copy <- ""


if ( getRversion() < "4.1.0" )
  stop( "FlowPredict needs R 4.1 or later; this is R ", getRversion(),
        ". Install a newer R from https://cran.r-project.org first." )

options( repos = c( CRAN = "https://cloud.r-project.org" ) )

if ( !requireNamespace( "remotes", quietly = TRUE ) )
  install.packages( "remotes" )

if ( nzchar( local.copy ) ) {
  remotes::install_local( path.expand( local.copy ), dependencies = TRUE, upgrade = "never",
                          force = TRUE )
} else {
  remotes::install_github( github.repo, dependencies = TRUE, upgrade = "never" )
}

# check the installation on a small simulated data set
library( FlowPredict )
set.seed( 1 )
n <- 40
group <- rep( c( "Case", "Control" ), each = n / 2 )
X <- matrix( rnorm( n * 10 ), n, 10, dimnames = list( paste0( "S", 1:n ), paste0( "Feature", 1:10 ) ) )
X[ group == "Case", 1 ] <- X[ group == "Case", 1 ] + 2
test.data <- feature.data(
  data.frame( Samples = rownames( X ), X, check.names = FALSE ),
  data.frame( Samples = rownames( X ), Group = group ) )
check <- run.predictive( test.data, "Group", "Case", "Control", repeats = 1, verbose = FALSE )

if ( min( check$performance$auc.mean ) > 0.8 ) {
  cat( "\nFlowPredict", as.character( packageVersion( "FlowPredict" ) ),
       "is installed and working.\n" )
} else {
  warning( "FlowPredict installed, but the check on simulated data gave an unexpected result. ",
           "Please report this with the output of sessionInfo()." )
}
