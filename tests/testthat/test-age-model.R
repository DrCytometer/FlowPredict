simulate.age <- function( n = 60, p = 30, seed = 1 ) {
  set.seed( seed )
  age <- 20 + 40 * stats::runif( n )
  X <- matrix( stats::rnorm( n * p ), n, p, dimnames = list( paste0( "s", seq_len( n ) ),
                                                             paste0( "f", seq_len( p ) ) ) )
  X[ , 1 ] <- X[ , 1 ] + age / 8
  X[ , 2 ] <- X[ , 2 ] - age / 10
  group <- rep( c( "Control", "Case" ), length.out = n )
  features <- data.frame( Samples = rownames( X ), X, check.names = FALSE )
  metadata <- data.frame( Samples = rownames( X ), Group = group,
                          Age = ifelse( seq_len( n ) > n - 5, NA, age ) )
  feature.data( features, metadata )
}

test_that( "run.age.model predicts the variable and beats chance", {
  d <- simulate.age()
  a <- run.age.model( d, "Age", repeats = 2, n.permutations = 9, verbose = FALSE )
  expect_s3_class( a, "flowpredict.age" )
  expect_lt( a$performance$mae.mean, 0.6 * a$performance$mae.baseline )
  expect_gt( a$performance$spearman.mean, 0.7 )
  expect_equal( a$performance$p.chance, 0.1 )
  expect_true( all( c( "f1", "f2" ) %in% a$importance$feature[ 1:2 ] ) )
  expect_equal( sign( a$importance$coef[ a$importance$feature == "f2" ] ), -1 )
  expect_equal( a$summary$n, 55 )
} )

test_that( "corrected gap removes the trend of the gap on the variable", {
  d <- simulate.age( n = 80, seed = 3 )
  a <- run.age.model( d, "Age", repeats = 2, n.permutations = 0, verbose = FALSE )
  pr <- a$predictions
  raw <- stats::coef( stats::lm( gap ~ actual, pr ) )[ 2 ]
  corr <- stats::coef( stats::lm( gap.corrected ~ actual, pr ) )[ 2 ]
  expect_lt( raw, 0 )
  expect_lt( abs( corr ), abs( raw ) )
} )

test_that( "samples outside the model groups are predicted but not cross-validated", {
  d <- simulate.age()
  a <- run.age.model( d, "Age", group.column = "Group", model.groups = "Control",
                      repeats = 1, n.permutations = 0, verbose = FALSE )
  expect_true( all( a$predictions$group[ a$predictions$cross.validated ] == "Control" ) )
  expect_true( all( a$predictions$group[ !a$predictions$cross.validated ] == "Case" ) )
  expect_equal( nrow( a$by.group ), 2 )
  expect_true( a$by.group$in.model[ 1 ] )
} )

test_that( "too few samples are refused", {
  d <- simulate.age( n = 15 )
  expect_error( run.age.model( d, "Age", verbose = FALSE ), "at least 20" )
} )
