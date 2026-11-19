# run_predictive.R
#
# Nested cross-validated prediction of group membership (two groups) from
# per-sample features: lasso logistic regression and random forest, with
# feature importance and the minimal panel (AUC of the top-N features).


MIN.PER.GROUP <- 10
RECOMMENDED.PER.GROUP <- 30
MIN.DONORS.PER.GROUP <- 4
MODEL.LABELS <- c( lasso = "Lasso logistic regression",
                   random.forest = "Random forest" )
RF.DEPTHS <- c( 2, 4, 0 )
RF.TREES <- 300
LASSO.N.LAMBDA <- 50
ROC.GRID <- seq( 0, 1, length.out = 101 )
PANEL.SIZES <- c( 1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200 )


#' Predict group membership from features
#'
#' @title Run Predictive Analysis
#'
#' @description
#' Asks "which small set of features best separates two groups, and how
#' well?".  This complements per-feature tests (which ask whether each
#' feature differs on its own): here the features are used together to
#' predict each sample's group, and the prediction is always judged on
#' samples the model has not seen.
#'
#' \strong{Models.}  Two models are fitted, each with the groups weighted
#' equally whatever their sizes:
#' \itemize{
#'   \item Lasso logistic regression: a logistic regression (the standard
#'     model for a yes/no outcome) with a penalty that pushes the weights
#'     of uninformative features to exactly zero, so the model keeps a few
#'     features.  The penalty strength is chosen by an inner
#'     cross-validation on the training samples only.
#'   \item Random forest: many decision trees, each grown on a random
#'     resample of the samples and features, voting together.  It can use
#'     interactions and thresholds the lasso cannot.  Its maximum tree
#'     depth (2, 4 or unlimited) is chosen by its out-of-bag AUC (each tree
#'     scored on the samples left out of its resample).
#' }
#'
#' \strong{Nested cross-validation.}  The samples are split into
#' \code{outer.folds} parts; each part in turn is held out, the model is
#' built (including penalty or depth choice, imputation and scaling) on the
#' rest, and the held-out samples are predicted.  Repeating with
#' \code{repeats} different random splits gives every sample
#' \code{repeats} predictions from models that never saw it.  With
#' \code{donor.column}, every split keeps all samples of one donor on the
#' same side, so a model is never tested on a donor it learnt from.
#'
#' \strong{Performance.}  ROC AUC: the probability that a randomly chosen
#' sample of the positive group gets a higher score than a randomly chosen
#' sample of the negative group (0.5 is chance, 1 is perfect separation).
#' Balanced accuracy: the mean of the share of each group classified
#' correctly at a probability cut-off of 0.5.  Both are given as mean and
#' SD over repeats.
#'
#' \strong{Importance.}  For the lasso, how often each feature is kept
#' (selection frequency, over all outer fits) and its mean standardised
#' coefficient; for the forest, the mean decrease in Gini impurity
#' (normalised to sum to 1).  Correlated features (e.g. a parent population
#' and its largest subset) share importance: the lasso tends to keep one of
#' them, so a feature it drops may still differ between groups.
#'
#' \strong{Minimal panel.}  Within each training split the features are
#' ranked by Welch t statistic, and a ridge logistic regression using only
#' the top N features is evaluated for increasing N.  The smallest N
#' reaching 95\% of the best mean AUC is the minimal panel.
#'
#' \strong{Covariates} (optional, not in Honeychrome): metadata columns such
#' as age or sex added to every model.  They are never penalised in the
#' lasso and are always included in the panel models, so the features are
#' judged on what they add beyond the covariates.  Samples missing a
#' covariate are dropped.
#'
#' At least 10 samples per group are required; below 30 the estimates are
#' noisy and the result carries a warning.
#'
#' @param data A \code{flowpredict.data} object from
#'   \code{\link{read.feature.data}}.
#' @param group.column Metadata column holding the groups.
#' @param positive,negative The two groups to compare; other samples are
#'   ignored.  Direction and coefficients are given for \code{positive}.
#' @param models Any of \code{"lasso"}, \code{"random.forest"}.
#' @param outer.folds,repeats,inner.folds Cross-validation sizes.
#' @param max.panel Largest panel evaluated.
#' @param donor.column Optional metadata column naming the donor of each
#'   sample, to keep a donor's samples in one fold.
#' @param covariates Optional metadata columns added to every model.
#' @param seed Random seed (results are reproducible for a given seed).
#' @param verbose Print progress.
#'
#' @return A list of class \code{flowpredict.predictive}:
#' \describe{
#'   \item{performance}{model, auc.mean, auc.sd, balanced.accuracy.mean,
#'     balanced.accuracy.sd.}
#'   \item{roc}{model, fpr, tpr: mean ROC curve over repeats.}
#'   \item{importance}{feature, direction (+1 higher in \code{positive}),
#'     lasso.selection, lasso.coef, forest.importance, forest.importance.sd
#'     (for the models run), sorted by lasso selection then forest
#'     importance.}
#'   \item{panel}{n.features, auc.mean, auc.sd.}
#'   \item{predictions}{sample, group, model, probability (mean
#'     out-of-fold probability of \code{positive}), probability.sd.}
#'   \item{summary}{list: positive, negative, n.per.group, n.features,
#'     minimal.panel, best.panel.auc, status, message, outer.folds,
#'     inner.folds, repeats, n.donors, covariates.}
#' }
#'
#' @export
run.predictive <- function( data, group.column, positive, negative,
                            models = c( "lasso", "random.forest" ),
                            outer.folds = 5, repeats = 5, inner.folds = 3,
                            max.panel = 50, donor.column = NULL,
                            covariates = NULL, seed = 0, verbose = TRUE ) {
  models <- match.arg( models, names( MODEL.LABELS ), several.ok = TRUE )
  meta <- data$metadata
  if ( !group.column %in% names( meta ) )
    stop( "No metadata column named '", group.column, "'.", call. = FALSE )
  groups <- as.character( meta[[ group.column ]] )
  for ( g in c( positive, negative ) )
    if ( !g %in% groups )
      stop( "Group '", g, "' not found in column '", group.column, "'. Groups: ",
            paste( sort( unique( stats::na.omit( groups ) ) ), collapse = ", " ),
            call. = FALSE )
  if ( !is.null( donor.column ) && !donor.column %in% names( meta ) )
    stop( "No metadata column named '", donor.column, "'.", call. = FALSE )

  keep <- groups %in% c( positive, negative ) & covariates.complete( meta, covariates )
  dropped.cov <- sum( groups %in% c( positive, negative ) ) - sum( keep )
  if ( dropped.cov > 0 && verbose )
    message( dropped.cov, " sample(s) missing a covariate were left out." )
  X <- data$features[ keep, , drop = FALSE ]
  y <- as.integer( groups[ keep ] == positive )
  samples <- rownames( X )
  n.per.group <- stats::setNames( c( sum( y == 1 ), sum( y == 0 ) ), c( positive, negative ) )
  status <- predictive.sample.size.status( n.per.group )
  if ( status$status == "block" ) stop( status$message, call. = FALSE )

  X <- X[ , usable.features( X ), drop = FALSE ]
  if ( ncol( X ) == 0 ) stop( "No usable features (all missing or constant).", call. = FALSE )
  features <- colnames( X )
  Z <- covariate.matrix( meta[ keep, , drop = FALSE ], covariates )
  n.cov <- if ( is.null( Z ) ) 0L else ncol( Z )
  XZ <- if ( n.cov ) cbind( X, Z ) else X
  feature.idx <- seq_along( features )

  donors <- NULL
  if ( !is.null( donor.column ) ) {
    donors <- as.character( meta[[ donor.column ]][ keep ] )
    donors[ is.na( donors ) ] <- samples[ is.na( donors ) ]
    if ( length( unique( donors ) ) == length( donors ) ) donors <- NULL
  }
  if ( is.null( donors ) ) {
    n.min <- min( n.per.group )
  } else {
    n.min <- min( length( unique( donors[ y == 1 ] ) ), length( unique( donors[ y == 0 ] ) ) )
    if ( n.min < MIN.DONORS.PER.GROUP )
      stop( "Keeping donors together needs at least ", MIN.DONORS.PER.GROUP,
            " donors in each group; the smaller group has ", n.min,
            ". Leave donor.column out to split samples alone.", call. = FALSE )
  }
  k.outer <- as.integer( min( outer.folds, n.min ) )
  k.inner <- if ( n.min > 3 ) as.integer( min( inner.folds, n.min - n.min %/% k.outer - 1 ) ) else 2L
  k.inner <- max( 2L, k.inner )

  sizes <- PANEL.SIZES[ PANEL.SIZES <= min( length( features ), max.panel ) ]
  upper <- min( length( features ), max.panel )
  if ( !upper %in% sizes ) sizes <- c( sizes, upper )

  oof <- lapply( stats::setNames( models, models ),
                 function( m ) matrix( NA_real_, repeats, length( y ) ) )
  panel.oof <- array( NA_real_, c( repeats, length( sizes ), length( y ) ) )
  lasso.nonzero <- numeric( length( features ) )
  lasso.coef.sum <- numeric( length( features ) )
  n.lasso.fits <- 0L
  forest.imp <- list()

  for ( r in seq_len( repeats ) ) {
    outer <- seeded( seed + r, if ( is.null( donors ) ) stratified.folds( y, k.outer )
                        else stratified.group.folds( y, donors, k.outer ) )
    for ( f in seq_len( k.outer ) ) {
      if ( verbose ) message( "Predictive: repeat ", r, "/", repeats, ", fold ", f, "/", k.outer )
      tr <- which( outer != f ); te <- which( outer == f )
      inner <- seeded( seed + 1000 + r,
                          if ( is.null( donors ) ) stratified.folds( y[ tr ], k.inner )
                          else stratified.group.folds( y[ tr ], donors[ tr ], k.inner ) )
      if ( "lasso" %in% models ) {
        fit <- fit.lasso( XZ[ tr, , drop = FALSE ], y[ tr ], inner, n.cov )
        oof$lasso[ r, te ] <- lasso.predict( fit, XZ[ te, , drop = FALSE ] )
        coef <- fit$coef[ feature.idx ]
        lasso.nonzero <- lasso.nonzero + ( coef != 0 )
        lasso.coef.sum <- lasso.coef.sum + coef
        n.lasso.fits <- n.lasso.fits + 1L
      }
      if ( "random.forest" %in% models ) {
        fit <- fit.forest( XZ[ tr, , drop = FALSE ], y[ tr ], seed + r )
        oof$random.forest[ r, te ] <- forest.predict( fit, XZ[ te, , drop = FALSE ] )
        forest.imp[[ length( forest.imp ) + 1 ]] <- fit$importance[ feature.idx ]
      }
      ranking <- order( -welch.abs.t( X[ tr, , drop = FALSE ], y[ tr ] ) )
      for ( si in seq_along( sizes ) ) {
        cols <- c( ranking[ seq_len( sizes[ si ] ) ],
                   if ( n.cov ) length( features ) + seq_len( n.cov ) )
        panel.oof[ r, si, te ] <- fit.predict.ridge( XZ[ tr, cols, drop = FALSE ], y[ tr ],
                                                     XZ[ te, cols, drop = FALSE ] )
      }
    }
  }

  perf <- list(); roc <- list(); preds <- list()
  for ( m in models ) {
    aucs <- apply( oof[[ m ]], 1, function( s ) roc.auc( y, s ) )
    baccs <- apply( oof[[ m ]], 1, function( s ) balanced.accuracy( y, as.integer( s >= 0.5 ) ) )
    tprs <- t( apply( oof[[ m ]], 1, function( s ) roc.on.grid( y, s, ROC.GRID ) ) )
    mean.tpr <- colMeans( matrix( tprs, nrow = repeats ) )
    mean.tpr[ 1 ] <- 0
    perf[[ m ]] <- data.frame( model = MODEL.LABELS[[ m ]],
                               auc.mean = mean( aucs ), auc.sd = sd.or.zero( aucs ),
                               balanced.accuracy.mean = mean( baccs ),
                               balanced.accuracy.sd = sd.or.zero( baccs ) )
    roc[[ m ]] <- data.frame( model = MODEL.LABELS[[ m ]], fpr = ROC.GRID, tpr = mean.tpr )
    preds[[ m ]] <- data.frame( sample = samples, group = ifelse( y == 1, positive, negative ),
                                model = MODEL.LABELS[[ m ]],
                                probability = colMeans( oof[[ m ]] ),
                                probability.sd = if ( repeats > 1 ) apply( oof[[ m ]], 2, stats::sd ) else 0 )
  }

  diff <- colMeans( X[ y == 1, , drop = FALSE ], na.rm = TRUE ) -
    colMeans( X[ y == 0, , drop = FALSE ], na.rm = TRUE )
  diff[ !is.finite( diff ) ] <- 0
  importance <- data.frame( feature = features, direction = sign( diff ),
                            stringsAsFactors = FALSE, check.names = FALSE )
  if ( "lasso" %in% models && n.lasso.fits ) {
    importance$lasso.selection <- lasso.nonzero / n.lasso.fits
    importance$lasso.coef <- lasso.coef.sum / n.lasso.fits
  }
  if ( "random.forest" %in% models && length( forest.imp ) ) {
    stack <- do.call( rbind, forest.imp )
    importance$forest.importance <- colMeans( stack )
    importance$forest.importance.sd <- if ( nrow( stack ) > 1 ) apply( stack, 2, stats::sd ) else 0
  }
  sort.cols <- intersect( c( "lasso.selection", "forest.importance" ), names( importance ) )
  ord <- do.call( order, c( lapply( sort.cols, function( cc ) -importance[[ cc ]] ),
                            list( seq_len( nrow( importance ) ) ) ) )
  importance <- importance[ ord, , drop = FALSE ]
  rownames( importance ) <- NULL

  panel.auc <- matrix( NA_real_, repeats, length( sizes ) )
  for ( r in seq_len( repeats ) ) for ( si in seq_along( sizes ) )
    panel.auc[ r, si ] <- roc.auc( y, panel.oof[ r, si, ] )
  panel <- data.frame( n.features = sizes, auc.mean = colMeans( panel.auc ),
                       auc.sd = if ( repeats > 1 ) apply( panel.auc, 2, stats::sd ) else 0 )
  best <- max( panel$auc.mean )
  minimal <- min( panel$n.features[ panel$auc.mean >= 0.95 * best ] )

  summary <- list( positive = positive, negative = negative, n.per.group = n.per.group,
                   n.features = length( features ), minimal.panel = minimal,
                   best.panel.auc = best, status = status$status, message = status$message,
                   outer.folds = k.outer, inner.folds = k.inner, repeats = repeats,
                   n.donors = if ( is.null( donors ) ) NA_integer_ else length( unique( donors ) ),
                   covariates = covariates, seed = seed )
  if ( verbose ) message( "Predictive: done (", length( features ), " features, minimal panel ",
                          minimal, ")." )
  structure( list( performance = do.call( rbind, unname( perf ) ),
                   roc = do.call( rbind, unname( roc ) ), importance = importance,
                   panel = panel, predictions = do.call( rbind, unname( preds ) ),
                   summary = summary ),
             class = "flowpredict.predictive" )
}


#' Sample size check for the predictive analysis
#'
#' @param n.per.group Named vector of samples per group.
#' @return list(status = "block" | "warn" | "ok", message).
#' @keywords internal
predictive.sample.size.status <- function( n.per.group ) {
  detail <- paste( names( n.per.group ), n.per.group, sep = "=", collapse = ", " )
  smallest <- min( n.per.group )
  if ( smallest < MIN.PER.GROUP )
    return( list( status = "block", message = paste0(
      "Predictive analysis needs at least ", MIN.PER.GROUP, " samples per group (", detail,
      "). With fewer, cross-validated performance estimates are too unstable to ",
      "interpret; use per-feature tests instead." ) ) )
  if ( smallest < RECOMMENDED.PER.GROUP )
    return( list( status = "warn", message = paste0(
      detail, ": results will be noisy. ", RECOMMENDED.PER.GROUP,
      " or more samples per group are recommended." ) ) )
  list( status = "ok", message = detail )
}


sd.or.zero <- function( v ) if ( length( v ) > 1 ) stats::sd( v ) else 0


class.weights <- function( y ) {
  n <- length( y )
  ifelse( y == 1, n / ( 2 * sum( y == 1 ) ), n / ( 2 * sum( y == 0 ) ) )
}


# Lasso logistic regression: penalty chosen by mean inner-fold AUC
# (largest penalty among ties); covariate columns (the last n.cov) are not
# penalised.
fit.lasso <- function( X, y, inner, n.cov = 0L ) {
  pf <- c( rep( 1, ncol( X ) - n.cov ), rep( 0, n.cov ) )
  pp <- fit.preprocess( X )
  Xs <- apply.preprocess( pp, X )
  w <- class.weights( y )
  path <- glmnet::glmnet( Xs, y, family = "binomial", alpha = 1, weights = w,
                          standardize = FALSE, nlambda = LASSO.N.LAMBDA,
                          penalty.factor = pf )
  lambda <- path$lambda
  fold.auc <- matrix( NA_real_, max( inner ), length( lambda ) )
  for ( k in sort( unique( inner ) ) ) {
    itr <- which( inner != k ); ite <- which( inner == k )
    if ( length( unique( y[ itr ] ) ) < 2 || length( unique( y[ ite ] ) ) < 2 ) next
    ipp <- fit.preprocess( X[ itr, , drop = FALSE ] )
    fit <- glmnet::glmnet( apply.preprocess( ipp, X[ itr, , drop = FALSE ] ), y[ itr ],
                           family = "binomial", alpha = 1, weights = class.weights( y[ itr ] ),
                           standardize = FALSE, lambda = lambda, penalty.factor = pf )
    s <- stats::predict( fit, apply.preprocess( ipp, X[ ite, , drop = FALSE ] ),
                         s = lambda, type = "response" )
    fold.auc[ k, ] <- apply( s, 2, function( v ) roc.auc( y[ ite ], v ) )
  }
  mean.auc <- colMeans( fold.auc, na.rm = TRUE )
  mean.auc[ !is.finite( mean.auc ) ] <- -Inf
  best <- which( mean.auc == max( mean.auc ) )[ 1 ]
  coef <- as.numeric( stats::coef( path, s = lambda[ best ] ) )[ -1 ]
  list( preprocess = pp, path = path, lambda = lambda[ best ], coef = coef )
}

lasso.predict <- function( fit, X ) {
  as.numeric( stats::predict( fit$path, apply.preprocess( fit$preprocess, X ),
                              s = fit$lambda, type = "response" ) )
}


# Random forest on median-imputed features, depth chosen by out-of-bag AUC.
fit.forest <- function( X, y, seed ) {
  pp <- fit.preprocess( X )
  Xi <- apply.preprocess( pp, X, scale = FALSE )
  colnames( Xi ) <- paste0( "f", seq_len( ncol( Xi ) ) )
  yf <- factor( y, levels = c( 0, 1 ) )
  cw <- c( length( y ) / ( 2 * sum( y == 0 ) ), length( y ) / ( 2 * sum( y == 1 ) ) )
  best <- NULL; best.auc <- -Inf
  for ( depth in RF.DEPTHS ) {
    rf <- ranger::ranger( x = Xi, y = yf, num.trees = RF.TREES, probability = TRUE,
                          mtry = max( 1L, floor( sqrt( ncol( Xi ) ) ) ),
                          max.depth = depth, class.weights = cw, importance = "impurity",
                          num.threads = worker.threads(), seed = seed )
    oob <- rf$predictions[ , "1" ]
    ok <- is.finite( oob )
    auc <- if ( length( unique( y[ ok ] ) ) == 2 ) roc.auc( y[ ok ], oob[ ok ] ) else -Inf
    if ( auc > best.auc ) { best <- rf; best.auc <- auc }
  }
  imp <- best$variable.importance
  imp <- if ( sum( imp ) > 0 ) imp / sum( imp ) else imp
  list( preprocess = pp, forest = best, importance = as.numeric( imp ) )
}

forest.predict <- function( fit, X ) {
  Xi <- apply.preprocess( fit$preprocess, X, scale = FALSE )
  colnames( Xi ) <- paste0( "f", seq_len( ncol( Xi ) ) )
  stats::predict( fit$forest, data = Xi, num.threads = worker.threads() )$predictions[ , "1" ]
}


# Ridge logistic regression with balanced class weights and inverse
# penalty C = 1 (minimises C * sum(w * log loss) + |beta|^2 / 2, intercept
# unpenalised), fitted by Newton's method; returns the predicted
# probability for the test rows.
fit.predict.ridge <- function( X, y, Xtest, C = 1 ) {
  pp <- fit.preprocess( X )
  A <- cbind( 1, apply.preprocess( pp, X ) )
  w <- class.weights( y )
  pen <- c( 0, rep( 1, ncol( A ) - 1 ) )
  beta <- numeric( ncol( A ) )
  for ( it in seq_len( 100 ) ) {
    p <- stats::plogis( drop( A %*% beta ) )
    grad <- C * drop( crossprod( A, w * ( p - y ) ) ) + pen * beta
    hess <- C * crossprod( A, A * ( w * p * ( 1 - p ) ) ) + diag( pen, ncol( A ) )
    step <- solve( hess, grad )
    beta <- beta - step
    if ( max( abs( step ) ) < 1e-8 ) break
  }
  stats::plogis( drop( cbind( 1, apply.preprocess( pp, Xtest ) ) %*% beta ) )
}


#' @export
print.flowpredict.predictive <- function( x, ... ) {
  s <- x$summary
  cat( "Predictive analysis:", s$positive, "vs", s$negative, "(",
       paste( names( s$n.per.group ), s$n.per.group, sep = " = ", collapse = ", " ), ")\n" )
  cat( s$n.features, "features;", s$repeats, "x", s$outer.folds, "fold nested CV",
       if ( !is.na( s$n.donors ) ) paste( "keeping", s$n.donors, "donors together" ) else "",
       "\n" )
  if ( s$status != "ok" ) cat( "Note:", s$message, "\n" )
  p <- x$performance
  for ( i in seq_len( nrow( p ) ) )
    cat( sprintf( "  %-27s AUC %.3f +/- %.3f, balanced accuracy %.3f +/- %.3f\n",
                  p$model[ i ], p$auc.mean[ i ], p$auc.sd[ i ],
                  p$balanced.accuracy.mean[ i ], p$balanced.accuracy.sd[ i ] ) )
  cat( "Minimal panel:", s$minimal.panel, "feature(s) (best panel AUC",
       sprintf( "%.3f", s$best.panel.auc ), ")\n" )
  invisible( x )
}
