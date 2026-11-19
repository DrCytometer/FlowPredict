# run_age_model.R
#
# Immune age model: nested cross-validated elastic net predicting a
# continuous variable (age, days after vaccination, any number) from
# per-sample features, with the age gap corrected for regression to the
# mean and a permutation check against chance.


MIN.AGE.SAMPLES <- 20
RECOMMENDED.AGE.SAMPLES <- 40
AGE.L1.RATIOS <- c( 0.2, 0.5, 0.8, 1.0 )
AGE.N.LAMBDA <- 30
AGE.MIN.OBSERVED <- 0.8


#' Predict a continuous variable (immune age) from features
#'
#' @title Run Immune Age Model
#'
#' @description
#' How well do the features predict a number such as age, which features
#' carry the prediction, and is each sample predicted older or younger than
#' it is?
#'
#' In plain terms, with age as the example: a model learns to guess each
#' sample's age from its immune profile alone, and is always judged on
#' samples it was not trained on.  If it guesses well, the immune system
#' changes with age in a consistent way; samples it places older or younger
#' than they are may be ahead or behind for their age.
#'
#' \strong{Model.}  Elastic net linear regression (Zou & Hastie 2005): a
#' penalised regression mixing the lasso, which keeps few features, and
#' ridge, which shares weight among correlated features such as a parent
#' population and its largest subset.  Features are median-imputed and
#' standardised inside each training split.  The penalty strength and the
#' lasso/ridge mix (0.2, 0.5, 0.8, 1) are chosen by an inner
#' cross-validation on the training samples.
#'
#' \strong{Nested cross-validation.}  As in \code{\link{run.predictive}}:
#' \code{outer.folds} folds, \code{repeats} times, so every sample is
#' predicted by models that never saw it.  Folds keep every donor whole
#' (\code{donor.column}), otherwise the same donor is in training and
#' testing and the model looks better than it is.
#'
#' \strong{Performance.}  Mean absolute error (MAE, in the variable's
#' units), R squared (share of the variance explained; 0 or below is no
#' better than predicting the mean) and Spearman correlation of predicted
#' with actual values, as mean and SD over repeats.  \code{mae.baseline} is
#' the error of predicting every sample at the mean.
#'
#' \strong{Gap.}  Predicted minus actual value.  Any such model pulls its
#' predictions towards the middle (the lowest values predicted too high,
#' the highest too low), so the raw gap follows the variable itself.
#' \code{gap.corrected} removes that: the gap's linear trend on the
#' variable is estimated inside each training split from its own
#' cross-validated predictions and subtracted (de Lange & Cole 2020,
#' NeuroImage Clin 26:102229).  Only the corrected gap should be compared
#' between samples or groups.
#'
#' \strong{Chance check.}  The whole procedure (one repeat) is rerun with
#' the variable shuffled among donors \code{n.permutations} times;
#' p = (1 + shuffles with a MAE at least as small) / (1 + shuffles).
#'
#' \strong{Other groups.}  Samples not in \code{model.groups} (e.g. a
#' disease group, or another tissue) are predicted by a final model fitted
#' on all model samples and marked \code{cross.validated = FALSE}.
#'
#' \strong{Covariates} (optional, not in Honeychrome): metadata columns
#' added to the model without penalty.  Note that the prediction then uses
#' them too, so it is no longer an immune-only age.
#'
#' At least 20 samples with a value are required; below 40 the result
#' carries a warning.  Features measured in fewer than 80\% of the model
#' samples are dropped.
#'
#' @param data A \code{flowpredict.data} object.
#' @param value.column Metadata column holding the numeric variable.
#' @param group.column,model.groups Optional: train only on samples whose
#'   \code{group.column} is in \code{model.groups}; other samples are
#'   predicted by the final model.  Both \code{NULL}: every sample with a
#'   value is in the model.
#' @param donor.column Optional donor column (folds keep donors whole).
#' @param covariates Optional metadata columns added without penalty.
#' @param outer.folds,repeats,inner.folds Cross-validation sizes.
#' @param n.permutations Shuffles for the chance check (0 to skip).
#' @param seed Random seed.
#' @param verbose Print progress.
#'
#' @return A list of class \code{flowpredict.age}:
#' \describe{
#'   \item{performance}{One row: mae.mean, mae.sd, r2.mean, r2.sd,
#'     spearman.mean, spearman.sd, mae.baseline, p.chance.}
#'   \item{predictions}{sample, group, actual, predicted, predicted.sd,
#'     gap, gap.corrected, cross.validated.}
#'   \item{importance}{feature, selection (share of outer fits keeping
#'     it), coef (mean standardised coefficient), direction.}
#'   \item{chance}{mae of each shuffle.}
#'   \item{by.group}{Per group: in.model, n, predicted.mean, predicted.sd,
#'     n.actual, mae, mean.gap, mean.gap.corrected, spearman.}
#'   \item{summary}{list: n, n.features, n.blocks, outer.folds,
#'     inner.folds, repeats, n.permutations, status, message, penalty.}
#' }
#'
#' @export
run.age.model <- function( data, value.column, group.column = NULL, model.groups = NULL,
                           donor.column = NULL, covariates = NULL,
                           outer.folds = 5, repeats = 5, inner.folds = 3,
                           n.permutations = 100, seed = 0, verbose = TRUE ) {
  meta <- data$metadata
  for ( cc in c( value.column, group.column, donor.column ) )
    if ( !cc %in% names( meta ) ) stop( "No metadata column named '", cc, "'.", call. = FALSE )
  value <- suppressWarnings( as.numeric( meta[[ value.column ]] ) )
  if ( all( is.na( value ) ) )
    stop( "Column '", value.column, "' has no numeric values.", call. = FALSE )
  groups <- if ( is.null( group.column ) ) rep( "", nrow( meta ) )
            else as.character( meta[[ group.column ]] )
  in.groups <- if ( is.null( model.groups ) ) rep( TRUE, nrow( meta ) ) else groups %in% model.groups
  cov.ok <- covariates.complete( meta, covariates )
  in.model <- in.groups & is.finite( value ) & cov.ok
  other <- !in.groups & cov.ok

  X <- data$features[ in.model, , drop = FALSE ]
  y <- value[ in.model ]
  n <- length( y )
  status <- age.sample.size.status( n )
  if ( status$status == "block" ) stop( status$message, call. = FALSE )
  blocks <- if ( is.null( donor.column ) ) rownames( X )
            else as.character( meta[[ donor.column ]][ in.model ] )
  blocks[ is.na( blocks ) ] <- rownames( X )[ is.na( blocks ) ]

  X <- X[ , usable.features( X, AGE.MIN.OBSERVED ), drop = FALSE ]
  if ( ncol( X ) == 0 )
    stop( "No usable features (each must be measured in at least ",
          round( 100 * AGE.MIN.OBSERVED ), "% of the samples and vary).", call. = FALSE )
  features <- colnames( X )
  Z <- covariate.matrix( meta, covariates )
  n.cov <- if ( is.null( Z ) ) 0L else ncol( Z )
  XZ <- if ( n.cov ) cbind( X, Z[ in.model, , drop = FALSE ] ) else X
  feature.idx <- seq_along( features )
  k.outer <- as.integer( max( 2, min( outer.folds, length( unique( blocks ) ) ) ) )

  pred <- matrix( NA_real_, repeats, n )
  corrected <- matrix( NA_real_, repeats, n )
  coefs <- list(); chosen <- character( 0 )
  # one stream of random numbers for all outer splits and shuffles
  seeded( seed, {
    for ( r in seq_len( repeats ) ) {
      folds <- shuffled.group.folds( blocks, k.outer )
      for ( f in seq_len( k.outer ) ) {
        if ( verbose ) message( "Immune age: repeat ", r, "/", repeats, ", fold ", f, "/", k.outer )
        tr <- which( folds != f ); te <- which( folds == f )
        fit <- fit.age( XZ[ tr, , drop = FALSE ], y[ tr ], blocks[ tr ], inner.folds, n.cov )
        pred[ r, te ] <- age.predict( fit, XZ[ te, , drop = FALSE ] )
        corrected[ r, te ] <- ( pred[ r, te ] - y[ te ] ) -
          ( fit$bias[ 1 ] + fit$bias[ 2 ] * y[ te ] )
        coefs[[ length( coefs ) + 1 ]] <- fit$coef[ feature.idx ]
        chosen <- c( chosen, sprintf( "%.3g / %g", fit$lambda, fit$alpha ) )
      }
    }
    scores <- do.call( rbind, lapply( seq_len( repeats ), function( r ) age.scores( y, pred[ r, ] ) ) )

    chance <- numeric( 0 )
    shared <- all( tapply( y, blocks, function( v ) length( unique( v ) ) ) <= 1 )
    for ( i in seq_len( n.permutations ) ) {
      if ( verbose ) message( "Immune age: chance check ", i, "/", n.permutations )
      if ( shared ) {
        uniq <- unique( blocks )
        vals <- y[ match( uniq, blocks ) ]
        yp <- vals[ sample.int( length( vals ) ) ][ match( blocks, uniq ) ]
      } else {
        yp <- y[ sample.int( n ) ]
      }
      pp <- rep( NA_real_, n )
      folds <- shuffled.group.folds( blocks, k.outer )
      for ( f in seq_len( k.outer ) ) {
        tr <- which( folds != f ); te <- which( folds == f )
        fit <- fit.age( XZ[ tr, , drop = FALSE ], yp[ tr ], blocks[ tr ], inner.folds, n.cov,
                        bias = FALSE )
        pp[ te ] <- age.predict( fit, XZ[ te, , drop = FALSE ] )
      }
      chance <- c( chance, age.scores( yp, pp )$mae )
    }
  } )
  mae.mean <- mean( scores$mae )
  p.chance <- if ( length( chance ) ) ( 1 + sum( chance <= mae.mean ) ) / ( 1 + length( chance ) )
              else NA_real_

  C <- do.call( rbind, coefs )
  importance <- data.frame( feature = features, selection = colMeans( C != 0 ),
                            coef = colMeans( C ), stringsAsFactors = FALSE )
  importance$direction <- sign( importance$coef )
  importance <- importance[ order( -importance$selection, -abs( importance$coef ) ), ]
  rownames( importance ) <- NULL

  samples <- rownames( X )
  mean.pred <- colMeans( pred )
  predictions <- data.frame(
    sample = samples, group = groups[ in.model ], actual = y, predicted = mean.pred,
    predicted.sd = if ( repeats > 1 ) apply( pred, 2, stats::sd ) else 0,
    gap = mean.pred - y, gap.corrected = colMeans( corrected ), cross.validated = TRUE,
    stringsAsFactors = FALSE )

  if ( any( other ) ) {
    if ( verbose ) message( "Immune age: predicting the samples outside the model" )
    fit <- fit.age( XZ, y, blocks, inner.folds, n.cov )
    Xo <- data$features[ other, features, drop = FALSE ]
    if ( n.cov ) Xo <- cbind( Xo, Z[ other, , drop = FALSE ] )
    po <- age.predict( fit, Xo )
    yo <- value[ other ]
    predictions <- rbind( predictions, data.frame(
      sample = rownames( data$features )[ other ], group = groups[ other ], actual = yo,
      predicted = po, predicted.sd = NA_real_, gap = po - yo,
      gap.corrected = ( po - yo ) - ( fit$bias[ 1 ] + fit$bias[ 2 ] * yo ),
      cross.validated = FALSE, stringsAsFactors = FALSE ) )
  }

  tab <- sort( table( chosen ), decreasing = TRUE )
  performance <- data.frame(
    mae.mean = mae.mean, mae.sd = sd.or.zero( scores$mae ),
    r2.mean = mean( scores$r2 ), r2.sd = sd.or.zero( scores$r2 ),
    spearman.mean = mean( scores$spearman ), spearman.sd = sd.or.zero( scores$spearman ),
    mae.baseline = mean( abs( y - mean( y ) ) ) * n / max( n - 1, 1 ),
    p.chance = p.chance )
  summary <- list( n = n, n.features = length( features ), n.blocks = length( unique( blocks ) ),
                   outer.folds = k.outer, inner.folds = inner.folds, repeats = repeats,
                   n.permutations = length( chance ), status = status$status,
                   message = status$message, value.column = value.column,
                   model.groups = model.groups, covariates = covariates, seed = seed,
                   penalty = if ( length( tab ) ) names( tab )[ 1 ] else "" )
  if ( verbose ) message( "Immune age: done (mean absolute error ", signif( mae.mean, 3 ), ")." )
  res <- structure( list( performance = performance, predictions = predictions,
                          importance = importance, chance = data.frame( mae = chance ),
                          summary = summary ),
                    class = "flowpredict.age" )
  res$by.group <- age.model.by.group( res )
  res
}


#' Sample size check for the immune age model
#'
#' @param n Samples with a value.
#' @return list(status = "block" | "warn" | "ok", message).
#' @keywords internal
age.sample.size.status <- function( n ) {
  if ( n < MIN.AGE.SAMPLES )
    return( list( status = "block", message = paste0(
      "An immune age model needs at least ", MIN.AGE.SAMPLES, " samples with a value (", n,
      " here): with fewer, its cross-validated error is too unstable to interpret." ) ) )
  if ( n < RECOMMENDED.AGE.SAMPLES )
    return( list( status = "warn", message = paste0(
      n, " samples: results will be noisy; ", RECOMMENDED.AGE.SAMPLES,
      " or more are recommended." ) ) )
  list( status = "ok", message = paste( n, "samples" ) )
}


age.scores <- function( y, pred ) {
  data.frame( mae = mean( abs( pred - y ) ),
              r2 = 1 - sum( ( pred - y )^2 ) / sum( ( y - mean( y ) )^2 ),
              spearman = suppressWarnings( stats::cor( pred, y, method = "spearman" ) ) )
}


# Elastic net with penalty and lasso/ridge mix chosen by grouped inner
# cross-validation (lowest mean squared error), and the linear bias of its
# predictions on the variable (intercept, slope), estimated from inner
# cross-validated predictions at the chosen penalty.
fit.age <- function( X, y, blocks, inner.folds, n.cov = 0L, bias = TRUE ) {
  pf <- c( rep( 1, ncol( X ) - n.cov ), rep( 0, n.cov ) )
  k <- as.integer( max( 2, min( inner.folds, length( unique( blocks ) ) ) ) )
  foldid <- balanced.group.folds( blocks, k )
  pp <- fit.preprocess( X )
  Xs <- apply.preprocess( pp, X )
  best <- NULL
  for ( a in AGE.L1.RATIOS ) {
    cv <- glmnet::cv.glmnet( Xs, y, family = "gaussian", alpha = a, foldid = foldid,
                             nlambda = AGE.N.LAMBDA, standardize = FALSE,
                             penalty.factor = pf, thresh = 1e-5 )
    err <- min( cv$cvm )
    if ( is.null( best ) || err < best$err ) best <- list( err = err, cv = cv, alpha = a )
  }
  lambda <- best$cv$lambda.min
  coef <- as.numeric( stats::coef( best$cv$glmnet.fit, s = lambda ) )[ -1 ]
  b <- c( 0, 0 )
  if ( bias ) {
    inner.pred <- rep( NA_real_, length( y ) )
    for ( f in sort( unique( foldid ) ) ) {
      itr <- which( foldid != f ); ite <- which( foldid == f )
      ipp <- fit.preprocess( X[ itr, , drop = FALSE ] )
      fit <- glmnet::glmnet( apply.preprocess( ipp, X[ itr, , drop = FALSE ] ), y[ itr ],
                             family = "gaussian", alpha = best$alpha,
                             lambda = best$cv$lambda, standardize = FALSE,
                             penalty.factor = pf, thresh = 1e-5 )
      inner.pred[ ite ] <- as.numeric( stats::predict(
        fit, apply.preprocess( ipp, X[ ite, , drop = FALSE ] ), s = lambda ) )
    }
    b <- as.numeric( stats::coef( stats::lm( I( inner.pred - y ) ~ y ) ) )
  }
  list( preprocess = pp, fit = best$cv$glmnet.fit, lambda = lambda, alpha = best$alpha,
        coef = coef, bias = b )
}

age.predict <- function( fit, X ) {
  as.numeric( stats::predict( fit$fit, apply.preprocess( fit$preprocess, X ), s = fit$lambda ) )
}


#' Immune age results by group
#'
#' @title Immune Age by Group
#'
#' @description
#' How the immune age model sees each group: the model's own samples
#' (cross-validated predictions) and every other group (predicted by the
#' final model), whether or not its samples have a value of the variable.
#'
#' @param result A \code{flowpredict.age} result.
#'
#' @return Data frame: group, in.model, n, predicted.mean, predicted.sd,
#'   n.actual, mae, mean.gap, mean.gap.corrected, spearman (NA with fewer
#'   than 3 samples with a value).
#'
#' @export
age.model.by.group <- function( result ) {
  pred <- result$predictions
  key <- paste( pred$group, pred$cross.validated, sep = "\r" )
  rows <- lapply( unique( key ), function( kk ) {
    part <- pred[ key == kk & is.finite( pred$predicted ), , drop = FALSE ]
    known <- part[ is.finite( part$actual ), , drop = FALSE ]
    nk <- nrow( known )
    data.frame(
      group = part$group[ 1 ], in.model = part$cross.validated[ 1 ], n = nrow( part ),
      predicted.mean = mean( part$predicted ),
      predicted.sd = if ( nrow( part ) > 1 ) stats::sd( part$predicted ) else NA_real_,
      n.actual = nk,
      mae = if ( nk ) mean( abs( known$predicted - known$actual ) ) else NA_real_,
      mean.gap = if ( nk ) mean( known$gap ) else NA_real_,
      mean.gap.corrected = if ( nk ) mean( known$gap.corrected ) else NA_real_,
      spearman = if ( nk >= 3 ) suppressWarnings(
        stats::cor( known$predicted, known$actual, method = "spearman" ) ) else NA_real_,
      stringsAsFactors = FALSE )
  } )
  out <- do.call( rbind, rows )
  out <- out[ order( !out$in.model ), , drop = FALSE ]
  rownames( out ) <- NULL
  out
}


#' @export
print.flowpredict.age <- function( x, ... ) {
  s <- x$summary; p <- x$performance
  cat( "Immune age model for '", s$value.column, "': ", s$n, " samples, ", s$n.features,
       " features; ", s$repeats, " x ", s$outer.folds, " fold nested CV\n", sep = "" )
  if ( s$status != "ok" ) cat( "Note:", s$message, "\n" )
  cat( sprintf( "  MAE %.3g +/- %.2g (predicting the mean: %.3g)\n", p$mae.mean, p$mae.sd,
                p$mae.baseline ) )
  cat( sprintf( "  R2 %.3f +/- %.3f, Spearman %.3f +/- %.3f\n", p$r2.mean, p$r2.sd,
                p$spearman.mean, p$spearman.sd ) )
  if ( is.finite( p$p.chance ) )
    cat( sprintf( "  Chance check: p = %.3g (%d shuffles)\n", p$p.chance, s$n.permutations ) )
  invisible( x )
}
