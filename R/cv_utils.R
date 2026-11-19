# cv_utils.R
#
# Shared pieces of the nested cross-validation: fold assignment (stratified,
# optionally keeping each donor's samples together), median imputation and
# standardisation fitted on training samples only, covariate matrices, and
# the performance measures (ROC AUC, ROC curve, balanced accuracy).


# Evaluate expr with the random seed set, restoring the caller's RNG state.
seeded <- function( seed, expr ) {
  had <- exists( ".Random.seed", envir = globalenv(), inherits = FALSE )
  if ( had ) old <- get( ".Random.seed", envir = globalenv() )
  on.exit( {
    if ( had ) assign( ".Random.seed", old, envir = globalenv() )
    else if ( exists( ".Random.seed", envir = globalenv(), inherits = FALSE ) )
      rm( ".Random.seed", envir = globalenv() )
  } )
  set.seed( seed )
  expr
}


# Stratified k folds: each class shuffled and dealt round the folds, so
# every fold has close to the overall class balance. Returns a fold number
# (1..k) per sample.
stratified.folds <- function( y, k ) {
  fold <- integer( length( y ) )
  start <- 0L
  for ( cl in sort( unique( y ) ) ) {
    idx <- which( y == cl )
    idx <- idx[ sample.int( length( idx ) ) ]
    fold[ idx ] <- ( ( seq_along( idx ) + start - 1L ) %% k ) + 1L
    start <- start + length( idx )
  }
  fold
}


# Stratified k folds keeping every donor whole. Donors are taken in random
# order (largest first) and each is placed in the fold whose class counts
# stay closest to the overall class balance, ties to the smallest fold.
stratified.group.folds <- function( y, donors, k ) {
  classes <- sort( unique( y ) )
  uniq <- unique( donors )
  uniq <- uniq[ sample.int( length( uniq ) ) ]
  counts <- t( vapply( uniq, function( d ) tabulate( match( y[ donors == d ], classes ),
                                                     length( classes ) ),
                       numeric( length( classes ) ) ) )
  if ( length( classes ) == 1 ) counts <- matrix( counts, ncol = 1 )
  ord <- order( -rowSums( counts ), seq_along( uniq ) )
  total <- colSums( counts )
  fold.counts <- matrix( 0, k, length( classes ) )
  donor.fold <- integer( length( uniq ) )
  for ( i in ord ) {
    cost <- vapply( seq_len( k ), function( f ) {
      trial <- fold.counts
      trial[ f, ] <- trial[ f, ] + counts[ i, ]
      share <- sweep( trial, 2, total, "/" )
      mean( apply( share, 2, stats::sd ) )
    }, numeric( 1 ) )
    best <- which( cost == min( cost ) )
    best <- best[ which.min( rowSums( fold.counts )[ best ] ) ]
    donor.fold[ i ] <- best
    fold.counts[ best, ] <- fold.counts[ best, ] + counts[ i, ]
  }
  donor.fold[ match( donors, uniq ) ]
}


# k folds keeping every block whole, blocks shuffled and dealt round.
shuffled.group.folds <- function( blocks, k ) {
  uniq <- unique( blocks )
  uniq <- uniq[ sample.int( length( uniq ) ) ]
  fold.of <- ( ( seq_along( uniq ) - 1L ) %% k ) + 1L
  fold.of[ match( blocks, uniq ) ]
}


# k folds keeping every block whole, without randomness: blocks from the
# largest down, each to the fold with the fewest samples so far.
balanced.group.folds <- function( blocks, k ) {
  tab <- table( blocks )
  uniq <- names( tab )[ order( -as.integer( tab ), names( tab ) ) ]
  size <- integer( k )
  fold.of <- integer( length( uniq ) )
  for ( i in seq_along( uniq ) ) {
    f <- which.min( size )
    fold.of[ i ] <- f
    size[ f ] <- size[ f ] + tab[[ uniq[ i ] ]]
  }
  fold.of[ match( blocks, uniq ) ]
}


# Median imputation and standardisation learnt from training rows.
# Columns empty in training are imputed with 0; constant columns are
# centred but not scaled.
fit.preprocess <- function( X ) {
  med <- apply( X, 2, stats::median, na.rm = TRUE )
  med[ !is.finite( med ) ] <- 0
  Xi <- impute.with( X, med )
  centre <- colMeans( Xi )
  sdev <- sqrt( colSums( sweep( Xi, 2, centre )^2 ) / nrow( Xi ) )
  sdev[ !is.finite( sdev ) | sdev < 1e-12 ] <- 1
  list( median = med, centre = centre, scale = sdev )
}

impute.with <- function( X, med ) {
  miss <- which( is.na( X ), arr.ind = TRUE )
  if ( nrow( miss ) ) X[ miss ] <- med[ miss[ , 2 ] ]
  X
}

apply.preprocess <- function( pp, X, scale = TRUE ) {
  X <- impute.with( X, pp$median )
  if ( scale ) X <- sweep( sweep( X, 2, pp$centre ), 2, pp$scale, "/" )
  X
}


# Covariate design matrix (no intercept; factors as indicator columns).
covariate.matrix <- function( metadata, covariates ) {
  if ( is.null( covariates ) || length( covariates ) == 0 ) return( NULL )
  missing.cols <- setdiff( covariates, names( metadata ) )
  if ( length( missing.cols ) )
    stop( "Covariate column(s) not in the metadata: ",
          paste( missing.cols, collapse = ", " ), call. = FALSE )
  df <- metadata[ , covariates, drop = FALSE ]
  for ( v in names( df ) ) if ( is.character( df[[ v ]] ) ) df[[ v ]] <- factor( df[[ v ]] )
  mf <- stats::model.frame( ~ ., df, na.action = stats::na.pass )
  M <- stats::model.matrix( ~ ., mf )[ , -1, drop = FALSE ]
  M <- M[ , apply( M, 2, function( v ) stats::sd( v, na.rm = TRUE ) > 1e-12 ), drop = FALSE ]
  rownames( M ) <- rownames( metadata )
  colnames( M ) <- paste0( "covariate: ", colnames( M ) )
  M
}


# Rows with every covariate present.
covariates.complete <- function( metadata, covariates ) {
  if ( is.null( covariates ) || length( covariates ) == 0 )
    return( rep( TRUE, nrow( metadata ) ) )
  stats::complete.cases( metadata[ , covariates, drop = FALSE ] )
}


# ROC AUC (Mann-Whitney), ties counted as half.
roc.auc <- function( y, score ) {
  ok <- is.finite( score )
  y <- y[ ok ]; score <- score[ ok ]
  n1 <- sum( y == 1 ); n0 <- sum( y == 0 )
  if ( n1 == 0 || n0 == 0 ) return( NA_real_ )
  r <- rank( score )
  ( sum( r[ y == 1 ] ) - n1 * ( n1 + 1 ) / 2 ) / ( n1 * n0 )
}


# ROC curve as (fpr, tpr), from (0, 0) to (1, 1).
roc.points <- function( y, score ) {
  thr <- sort( unique( score ), decreasing = TRUE )
  pos <- sum( y == 1 ); neg <- sum( y == 0 )
  tpr <- vapply( thr, function( t ) sum( score >= t & y == 1 ) / pos, numeric( 1 ) )
  fpr <- vapply( thr, function( t ) sum( score >= t & y == 0 ) / neg, numeric( 1 ) )
  list( fpr = c( 0, fpr ), tpr = c( 0, tpr ) )
}


# Mean ROC curve on a fixed false-positive grid (upper value at ties).
roc.on.grid <- function( y, score, grid ) {
  p <- roc.points( y, score )
  stats::approx( p$fpr, p$tpr, xout = grid, ties = max, rule = 2 )$y
}


balanced.accuracy <- function( y, predicted ) {
  mean( c( mean( predicted[ y == 1 ] == 1 ), mean( predicted[ y == 0 ] == 0 ) ) )
}


# |Welch t| per column, ignoring missing values; 0 where undefined.
welch.abs.t <- function( X, y ) {
  a <- X[ y == 1, , drop = FALSE ]; b <- X[ y == 0, , drop = FALSE ]
  na <- colSums( is.finite( a ) ); nb <- colSums( is.finite( b ) )
  ma <- colMeans( a, na.rm = TRUE ); mb <- colMeans( b, na.rm = TRUE )
  va <- apply( a, 2, stats::var, na.rm = TRUE )
  vb <- apply( b, 2, stats::var, na.rm = TRUE )
  t <- abs( ma - mb ) / sqrt( va / na + vb / nb )
  t[ !is.finite( t ) ] <- 0
  t
}


# Usable features: measured in at least min.observed of the rows and not
# constant.
usable.features <- function( X, min.observed = 0 ) {
  obs <- colMeans( !is.na( X ) )
  sdev <- apply( X, 2, stats::sd, na.rm = TRUE )
  sdev[ !is.finite( sdev ) ] <- 0
  ( obs > 0 ) & ( obs >= min.observed ) & ( sdev > 1e-12 )
}


worker.threads <- function() {
  max( 1L, ( parallel::detectCores( logical = TRUE ) %||% 2L ) - 1L )
}

`%||%` <- function( a, b ) if ( is.null( a ) || is.na( a ) ) b else a
