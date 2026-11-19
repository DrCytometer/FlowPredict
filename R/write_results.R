# write_results.R
#
# Writes the tables and figures of a predictive analysis or an immune age
# model to a folder.


#' Write the results of a predictive analysis
#'
#' @title Write Predictive Results
#'
#' @description
#' Writes the performance, importance, panel, ROC and per-sample prediction
#' tables as CSV files, and the ROC, importance and panel figures as PDF
#' and PNG, to \code{out.dir}.  File names start with \code{prefix}.
#'
#' @param result A \code{flowpredict.predictive} result.
#' @param out.dir Output folder (created if needed).
#' @param prefix File name prefix.
#' @param top Features shown in the importance figure.
#'
#' @return The paths written, invisibly.
#' @export
write.predictive.results <- function( result, out.dir, prefix = "predictive", top = 20 ) {
  dir.create( out.dir, showWarnings = FALSE, recursive = TRUE )
  s <- result$summary
  summary.df <- data.frame(
    item = c( "positive group", "negative group", paste( "samples:", names( s$n.per.group ) ),
              "features", "outer folds", "inner folds", "repeats", "donors kept together",
              "covariates", "minimal panel", "best panel AUC", "seed", "sample size note" ),
    value = c( s$positive, s$negative, s$n.per.group, s$n.features, s$outer.folds,
               s$inner.folds, s$repeats, if ( is.na( s$n.donors ) ) "no" else s$n.donors,
               if ( length( s$covariates ) ) paste( s$covariates, collapse = "; " ) else "none",
               s$minimal.panel, signif( s$best.panel.auc, 4 ), s$seed, s$message ) )
  tables <- list( summary = summary.df, performance = result$performance,
                  importance = result$importance, panel = result$panel, roc = result$roc,
                  predictions = result$predictions )
  figures <- list( roc = list( predictive.roc.plot( result ), 5.5, 6 ),
                   importance = list( predictive.importance.plot( result, top ), 11, 0.9 + 0.22 * top ),
                   panel = list( predictive.panel.plot( result ), 6, 4.5 ) )
  write.outputs( tables, figures, out.dir, prefix )
}


#' Write the results of an immune age model
#'
#' @title Write Immune Age Results
#'
#' @description
#' Writes the performance, predictions, importance, by-group and chance
#' tables as CSV files, and the prediction, gap, importance and chance
#' figures as PDF and PNG, to \code{out.dir}.
#'
#' @param result A \code{flowpredict.age} result.
#' @param out.dir Output folder (created if needed).
#' @param prefix File name prefix.
#' @param unit Optional unit of the variable, for axis labels.
#' @param top Features shown in the importance figure.
#'
#' @return The paths written, invisibly.
#' @export
write.age.results <- function( result, out.dir, prefix = "immune_age", unit = "", top = 20 ) {
  dir.create( out.dir, showWarnings = FALSE, recursive = TRUE )
  s <- result$summary
  summary.df <- data.frame(
    item = c( "variable", "samples in the model", "model groups", "features", "donors or blocks",
              "outer folds", "inner folds", "repeats", "chance shuffles", "covariates",
              "most chosen penalty (lambda / lasso share)", "seed", "sample size note" ),
    value = c( s$value.column, s$n,
               if ( length( s$model.groups ) ) paste( s$model.groups, collapse = "; " ) else "all",
               s$n.features, s$n.blocks, s$outer.folds, s$inner.folds, s$repeats,
               s$n.permutations,
               if ( length( s$covariates ) ) paste( s$covariates, collapse = "; " ) else "none",
               s$penalty, s$seed, s$message ) )
  tables <- list( summary = summary.df, performance = result$performance,
                  predictions = result$predictions, importance = result$importance,
                  by_group = result$by.group, chance = result$chance )
  figures <- list( predictions = list( age.prediction.plot( result, unit ), 6, 6 ),
                   gap = list( age.gap.plot( result, unit ), 5, 4.5 ),
                   importance = list( age.importance.plot( result, top ), 8, 0.9 + 0.22 * top ) )
  chance <- age.chance.plot( result )
  if ( !is.null( chance ) ) figures$chance <- list( chance, 5, 4 )
  write.outputs( tables, figures, out.dir, prefix )
}


write.outputs <- function( tables, figures, out.dir, prefix ) {
  paths <- character( 0 )
  for ( nm in names( tables ) ) {
    path <- file.path( out.dir, paste0( prefix, "_", nm, ".csv" ) )
    utils::write.csv( tables[[ nm ]], path, row.names = FALSE )
    paths <- c( paths, path )
  }
  for ( nm in names( figures ) ) {
    fig <- figures[[ nm ]]
    for ( ext in c( "pdf", "png" ) ) {
      path <- file.path( out.dir, paste0( prefix, "_", nm, ".", ext ) )
      ggplot2::ggsave( path, fig[[ 1 ]], width = fig[[ 2 ]], height = fig[[ 3 ]],
                       dpi = 200 )
      paths <- c( paths, path )
    }
  }
  invisible( paths )
}
