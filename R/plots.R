# plots.R
#
# ggplot2 figures of the predictive analysis and the immune age model.
# Each function returns a ggplot object; save with ggplot2::ggsave or
# write.predictive.results / write.age.results.


theme.flowpredict <- function() {
  ggplot2::theme_classic( base_size = 11 ) +
    ggplot2::theme( legend.position = "bottom", plot.title = ggplot2::element_text( face = "bold" ) )
}


#' ROC curves of the predictive analysis
#'
#' @title Predictive ROC Plot
#'
#' @description
#' Mean ROC curve of each model over repeats: the share of the positive
#' group detected (true positive rate) against the share of the negative
#' group wrongly called positive (false positive rate) as the probability
#' cut-off moves.  The diagonal is chance.
#'
#' @param result A \code{flowpredict.predictive} result.
#' @return A ggplot object.
#' @export
predictive.roc.plot <- function( result ) {
  roc <- result$roc
  perf <- result$performance
  lab <- stats::setNames( sprintf( "%s (AUC %.2f \u00b1 %.2f)", perf$model, perf$auc.mean,
                                   perf$auc.sd ), perf$model )
  roc$label <- factor( lab[ roc$model ], levels = lab )
  ggplot2::ggplot( roc, ggplot2::aes( fpr, tpr, colour = label ) ) +
    ggplot2::geom_abline( slope = 1, intercept = 0, linetype = 2, colour = "grey60" ) +
    ggplot2::geom_line( linewidth = 0.9 ) +
    ggplot2::scale_colour_brewer( palette = "Dark2", name = NULL ) +
    ggplot2::coord_equal( xlim = c( 0, 1 ), ylim = c( 0, 1 ) ) +
    ggplot2::labs( x = "False positive rate", y = "True positive rate",
                   title = paste( result$summary$positive, "vs", result$summary$negative ) ) +
    theme.flowpredict() +
    ggplot2::guides( colour = ggplot2::guide_legend( ncol = 1 ) )
}


#' Feature importance of the predictive analysis
#'
#' @title Predictive Importance Plot
#'
#' @description
#' The top features by lasso selection frequency (how often the lasso kept
#' the feature over all outer fits; features never kept are not shown) and
#' by random forest importance, side
#' by side, coloured by the group in which the feature is higher.
#'
#' @param result A \code{flowpredict.predictive} result.
#' @param top Number of features shown per panel.
#' @return A ggplot object.
#' @export
predictive.importance.plot <- function( result, top = 20 ) {
  imp <- result$importance
  s <- result$summary
  higher <- ifelse( imp$direction > 0, paste( "Higher in", s$positive ),
                    paste( "Higher in", s$negative ) )
  parts <- list()
  if ( "lasso.selection" %in% names( imp ) ) {
    o <- order( -imp$lasso.selection, -abs( imp$lasso.coef ) )
    o <- utils::head( o[ imp$lasso.selection[ o ] > 0 ], top )
    parts[[ 1 ]] <- data.frame( panel = "Lasso: selection frequency", feature = imp$feature[ o ],
                                value = imp$lasso.selection[ o ], higher = higher[ o ],
                                rank = seq_along( o ) )
  }
  if ( "forest.importance" %in% names( imp ) ) {
    o <- utils::head( order( -imp$forest.importance ), top )
    parts[[ length( parts ) + 1 ]] <- data.frame(
      panel = "Random forest: importance", feature = imp$feature[ o ],
      value = imp$forest.importance[ o ], higher = higher[ o ], rank = seq_along( o ) )
  }
  df <- do.call( rbind, parts )
  df$key <- factor( paste( df$panel, df$rank ), levels = rev( paste( df$panel, df$rank ) ) )
  labels <- stats::setNames( df$feature, as.character( df$key ) )
  fills <- stats::setNames( c( "#B2182B", "#2166AC" ),
                            paste( "Higher in", c( s$positive, s$negative ) ) )
  ggplot2::ggplot( df, ggplot2::aes( value, key, fill = higher ) ) +
    ggplot2::geom_col( width = 0.7 ) +
    ggplot2::facet_wrap( ~ panel, scales = "free", ncol = length( parts ) ) +
    ggplot2::scale_y_discrete( labels = labels ) +
    ggplot2::scale_fill_manual( values = fills, name = NULL ) +
    ggplot2::labs( x = NULL, y = NULL ) +
    theme.flowpredict() +
    ggplot2::theme( axis.text.y = ggplot2::element_text( size = 8 ) )
}


#' Panel size curve of the predictive analysis
#'
#' @title Predictive Panel Plot
#'
#' @description
#' Cross-validated AUC of a model using only the top N features (ranked
#' within each training split) for increasing N, with its SD over repeats.
#' The dashed line marks the minimal panel: the smallest N reaching 95\% of
#' the best mean AUC.
#'
#' @param result A \code{flowpredict.predictive} result.
#' @return A ggplot object.
#' @export
predictive.panel.plot <- function( result ) {
  p <- result$panel
  ggplot2::ggplot( p, ggplot2::aes( n.features, auc.mean ) ) +
    ggplot2::geom_ribbon( ggplot2::aes( ymin = pmax( 0, auc.mean - auc.sd ),
                                        ymax = pmin( 1, auc.mean + auc.sd ) ),
                          fill = "grey85" ) +
    ggplot2::geom_line( colour = "#1B9E77", linewidth = 0.9 ) +
    ggplot2::geom_point( colour = "#1B9E77", size = 1.8 ) +
    ggplot2::geom_vline( xintercept = result$summary$minimal.panel, linetype = 2,
                         colour = "grey40" ) +
    ggplot2::scale_x_log10() +
    ggplot2::labs( x = "Number of features (top ranked)", y = "Cross-validated AUC",
                   title = paste0( "Minimal panel: ", result$summary$minimal.panel,
                                   " feature(s)" ) ) +
    theme.flowpredict()
}


#' Predicted against actual values of the immune age model
#'
#' @title Immune Age Prediction Plot
#'
#' @description
#' Each sample's predicted value against its actual value.  Points on the
#' diagonal are predicted exactly; above it, predicted older than they
#' are.  Samples outside the model (predicted by the final model, not
#' cross-validated) are drawn as open circles.
#'
#' @param result A \code{flowpredict.age} result.
#' @param unit Optional unit for the axis labels (e.g. "weeks").
#' @return A ggplot object.
#' @export
age.prediction.plot <- function( result, unit = "" ) {
  pred <- result$predictions[ is.finite( result$predictions$actual ), , drop = FALSE ]
  pred$group[ pred$group == "" ] <- "Samples"
  pred$set <- ifelse( pred$cross.validated, "In the model (cross-validated)",
                      "Outside the model" )
  p <- result$performance
  u <- if ( nzchar( unit ) ) paste0( " (", unit, ")" ) else ""
  ggplot2::ggplot( pred, ggplot2::aes( actual, predicted, colour = group, shape = set ) ) +
    ggplot2::geom_abline( slope = 1, intercept = 0, linetype = 2, colour = "grey60" ) +
    ggplot2::geom_point( size = 2, alpha = 0.85 ) +
    ggplot2::scale_shape_manual( values = c( "In the model (cross-validated)" = 16,
                                             "Outside the model" = 1 ), name = NULL ) +
    ggplot2::scale_colour_brewer( palette = "Set2", name = NULL ) +
    ggplot2::labs( x = paste0( "Actual ", result$summary$value.column, u ),
                   y = paste0( "Predicted", u ),
                   title = sprintf( "MAE %.3g, R\u00b2 %.2f, Spearman %.2f", p$mae.mean,
                                    p$r2.mean, p$spearman.mean ) ) +
    theme.flowpredict() +
    ggplot2::guides( colour = ggplot2::guide_legend( ncol = 2 ),
                     shape = ggplot2::guide_legend( ncol = 1 ) )
}


#' Corrected age gap by group
#'
#' @title Immune Age Gap Plot
#'
#' @description
#' The corrected gap (predicted minus actual, with the regression to the
#' mean removed) of each sample, by group.  Above zero: predicted older
#' than the sample is.  Only groups whose samples have a value are shown.
#'
#' @param result A \code{flowpredict.age} result.
#' @param unit Optional unit for the axis label.
#' @return A ggplot object.
#' @export
age.gap.plot <- function( result, unit = "" ) {
  pred <- result$predictions[ is.finite( result$predictions$gap.corrected ), , drop = FALSE ]
  pred$group[ pred$group == "" ] <- "Samples"
  u <- if ( nzchar( unit ) ) paste0( " (", unit, ")" ) else ""
  ggplot2::ggplot( pred, ggplot2::aes( group, gap.corrected, fill = group ) ) +
    ggplot2::geom_hline( yintercept = 0, linetype = 2, colour = "grey60" ) +
    ggplot2::geom_boxplot( outlier.shape = NA, alpha = 0.6, width = 0.6 ) +
    ggplot2::geom_jitter( width = 0.15, height = 0, size = 1.5, colour = "grey30" ) +
    ggplot2::scale_fill_brewer( palette = "Set2", guide = "none" ) +
    ggplot2::labs( x = NULL, y = paste0( "Corrected gap", u ) ) +
    theme.flowpredict()
}


#' Feature importance of the immune age model
#'
#' @title Immune Age Importance Plot
#'
#' @description
#' The top features by selection frequency (share of outer fits keeping
#' the feature), with bars coloured by the sign of the mean coefficient:
#' rising or falling with the variable.
#'
#' @param result A \code{flowpredict.age} result.
#' @param top Number of features shown.
#' @return A ggplot object.
#' @export
age.importance.plot <- function( result, top = 20 ) {
  imp <- utils::head( result$importance, top )
  imp$key <- factor( seq_len( nrow( imp ) ), levels = rev( seq_len( nrow( imp ) ) ) )
  imp$trend <- ifelse( imp$coef > 0, "Rises with the variable", "Falls with the variable" )
  ggplot2::ggplot( imp, ggplot2::aes( selection, key, fill = trend ) ) +
    ggplot2::geom_col( width = 0.7 ) +
    ggplot2::scale_y_discrete( labels = stats::setNames( imp$feature, as.character( imp$key ) ) ) +
    ggplot2::scale_fill_manual( values = c( "Rises with the variable" = "#B2182B",
                                            "Falls with the variable" = "#2166AC" ),
                                name = NULL ) +
    ggplot2::labs( x = "Selection frequency", y = NULL ) +
    theme.flowpredict() +
    ggplot2::theme( axis.text.y = ggplot2::element_text( size = 8 ) )
}


#' Chance check of the immune age model
#'
#' @title Immune Age Chance Plot
#'
#' @description
#' Mean absolute error of the model refitted with the variable shuffled
#' (grey), against the real model's error (line).  A real error well to
#' the left of the shuffled ones means the prediction is not chance.
#'
#' @param result A \code{flowpredict.age} result.
#' @return A ggplot object, or NULL when the chance check was skipped.
#' @export
age.chance.plot <- function( result ) {
  if ( nrow( result$chance ) == 0 ) return( NULL )
  p <- result$performance
  ggplot2::ggplot( result$chance, ggplot2::aes( mae ) ) +
    ggplot2::geom_histogram( bins = 30, fill = "grey70", colour = "white" ) +
    ggplot2::geom_vline( xintercept = p$mae.mean, colour = "#B2182B", linewidth = 1 ) +
    ggplot2::labs( x = "Mean absolute error", y = "Shuffles",
                   title = sprintf( "Chance check: p = %.3g", p$p.chance ) ) +
    theme.flowpredict()
}
