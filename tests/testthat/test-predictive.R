simulate.data <- function( n = 60, p = 40, signal = 1.5, seed = 1 ) {
  set.seed( seed )
  group <- rep( c( "A", "B" ), each = n / 2 )
  X <- matrix( stats::rnorm( n * p ), n, p, dimnames = list( paste0( "s", seq_len( n ) ),
                                                             paste0( "feature ", seq_len( p ) ) ) )
  X[ group == "A", 1:3 ] <- X[ group == "A", 1:3 ] + signal
  age <- 20 + 30 * stats::runif( n )
  X[ , 4 ] <- X[ , 4 ] + age / 10
  X[ , 5 ] <- X[ , 5 ] - age / 15
  X[ 2, 6 ] <- NA
  features <- data.frame( Samples = rownames( X ), X, check.names = FALSE )
  metadata <- data.frame( Samples = rownames( X ), Group = group, Age = age,
                          Donor = paste0( "d", rep( seq_len( n / 2 ), 2 ) ),
                          Sex = rep( c( "F", "M" ), length.out = n ) )
  feature.data( features, metadata )
}

test_that( "feature.data aligns samples and rejects bad input", {
  d <- simulate.data()
  expect_s3_class( d, "flowpredict.data" )
  expect_identical( rownames( d$features ), d$metadata$Samples )
  f <- data.frame( Samples = c( "a", "b", "c" ), x = 1:3 )
  m <- data.frame( Samples = c( "c", "a" ), g = c( "u", "v" ) )
  expect_message( d2 <- feature.data( f, m ), "dropped" )
  expect_identical( rownames( d2$features ), c( "a", "c" ) )
  expect_error( feature.data( data.frame( Samples = "a", x = "text" ), m ), "numeric" )
  expect_error( feature.data( f, m, sample.column = "ID" ), "no column" )
  expect_error( feature.data( data.frame( Samples = "a", x = -1 ), m, transform = "log2" ),
                "non-negative" )
} )

test_that( "run.predictive finds the signal and ranks the true features near the top", {
  d <- simulate.data()
  r <- run.predictive( d, "Group", "A", "B", repeats = 2, verbose = FALSE )
  expect_s3_class( r, "flowpredict.predictive" )
  expect_gt( min( r$performance$auc.mean ), 0.85 )
  expect_true( all( paste0( "feature ", 1:3 ) %in% r$importance$feature[ 1:5 ] ) )
  expect_true( all( r$importance$direction[ r$importance$feature %in% paste0( "feature ", 1:3 ) ] == 1 ) )
  expect_lte( r$summary$minimal.panel, 5 )
  expect_equal( nrow( r$predictions ), 2 * 60 )
} )

test_that( "run.predictive is near chance without signal and reproducible", {
  d <- simulate.data( signal = 0, seed = 2 )
  r1 <- run.predictive( d, "Group", "A", "B", repeats = 2, verbose = FALSE )
  r2 <- run.predictive( d, "Group", "A", "B", repeats = 2, verbose = FALSE )
  expect_lt( max( r1$performance$auc.mean ), 0.75 )
  expect_identical( r1$performance, r2$performance )
} )

test_that( "donor folds keep each donor on one side", {
  y <- rep( c( 1, 0 ), each = 20 )
  donors <- paste0( "d", rep( 1:20, 2 ) )
  folds <- seeded( 3, stratified.group.folds( y, donors, 5 ) )
  expect_true( all( tapply( folds, donors, function( f ) length( unique( f ) ) ) == 1 ) )
  d <- simulate.data()
  r <- run.predictive( d, "Group", "A", "B", donor.column = "Donor", repeats = 1,
                       verbose = FALSE )
  expect_equal( r$summary$n.donors, 30 )
} )

test_that( "too few samples are refused", {
  d <- simulate.data( n = 16 )
  expect_error( run.predictive( d, "Group", "A", "B", verbose = FALSE ), "at least 10" )
} )

test_that( "covariates are accepted and kept out of the importance table", {
  d <- simulate.data()
  r <- run.predictive( d, "Group", "A", "B", covariates = "Sex", repeats = 1,
                       models = "lasso", verbose = FALSE )
  expect_false( any( grepl( "covariate", r$importance$feature ) ) )
  expect_equal( nrow( r$importance ), 40 )
} )

test_that( "ridge panel model matches an unpenalised fit when the penalty is negligible", {
  set.seed( 4 )
  X <- matrix( stats::rnorm( 200 ), 100, 2 )
  y <- as.integer( X[ , 1 ] + stats::rnorm( 100 ) > 0 )
  p.ridge <- fit.predict.ridge( X, y, X, C = 1e8 )
  w <- class.weights( y )
  Xs <- scale( X ) * sqrt( 100 / 99 )
  p.glm <- suppressWarnings( stats::fitted( stats::glm( y ~ Xs, family = stats::binomial, weights = w ) ) )
  expect_equal( p.ridge, as.numeric( p.glm ), tolerance = 1e-5 )
} )

test_that( "roc.auc matches the Mann-Whitney statistic", {
  y <- c( 1, 1, 0, 0, 1, 0 )
  s <- c( 0.9, 0.4, 0.4, 0.1, 0.8, 0.3 )
  w <- suppressWarnings( stats::wilcox.test( s[ y == 1 ], s[ y == 0 ] )$statistic )
  expect_equal( roc.auc( y, s ), as.numeric( w ) / 9 )
} )
