# read_feature_data.R
#
# Reads a feature table (one row per sample, one column per feature) and a
# metadata table (one row per sample: groups, donors, covariates) and lines
# them up by sample name.


#' Read a feature table and its sample metadata
#'
#' @title Read Feature Data
#'
#' @description
#' Reads two CSV files and matches their rows by sample name:
#' \itemize{
#'   \item the feature table: a column of sample names and any number of
#'     numeric feature columns (population frequencies, marker medians,
#'     counts or any other per-sample number);
#'   \item the metadata table: the same column of sample names and one
#'     column per sample property (group, donor, age, sex, batch).
#' }
#' Samples present in only one of the two files are dropped with a message.
#' Column names are kept exactly as written in the files (spaces, commas
#' and symbols included).
#'
#' Frequencies (percentages) are usually skewed: a population at 0.1\% and
#' one at 1\% differ more, biologically, than 50\% and 51\%.
#' \code{transform = "log2"} puts every feature on a log2 scale,
#' \code{log2(x + offset)}, as Honeychrome does for frequencies before
#' predictive modelling; the small offset keeps empty populations (0\%)
#' finite.  Use it for percentages or counts, not for marker medians, which
#' can be negative after unmixing.
#'
#' @param features.file Path to the feature CSV.
#' @param metadata.file Path to the metadata CSV.
#' @param sample.column Name of the sample column in both files.
#' @param transform \code{"none"} (features used as given) or \code{"log2"}.
#' @param offset Added before \code{log2}; 0.01 suits percentages.
#'
#' @return A list of class \code{flowpredict.data}: \code{features}, a
#'   numeric matrix (samples x features, sample names as row names), and
#'   \code{metadata}, a data frame in the same row order with the sample
#'   names in column \code{sample.column} and as row names.
#'
#' @export
read.feature.data <- function( features.file, metadata.file,
                               sample.column = "Samples",
                               transform = c( "none", "log2" ),
                               offset = 0.01 ) {
  transform <- match.arg( transform )
  features <- utils::read.csv( features.file, check.names = FALSE,
                               stringsAsFactors = FALSE,
                               na.strings = c( "", "NA", "NaN" ) )
  metadata <- utils::read.csv( metadata.file, check.names = FALSE,
                               stringsAsFactors = FALSE,
                               na.strings = c( "", "NA", "NaN" ) )
  feature.data( features, metadata, sample.column, transform, offset )
}


#' Build feature data from data frames already in R
#'
#' @title Feature Data From Data Frames
#'
#' @description
#' As \code{\link{read.feature.data}}, for a feature table and a metadata
#' table that are already data frames.
#'
#' @param features Data frame: a sample column and numeric feature columns.
#' @param metadata Data frame: a sample column and metadata columns.
#' @inheritParams read.feature.data
#'
#' @return A list of class \code{flowpredict.data}; see
#'   \code{\link{read.feature.data}}.
#'
#' @export
feature.data <- function( features, metadata, sample.column = "Samples",
                          transform = c( "none", "log2" ), offset = 0.01 ) {
  transform <- match.arg( transform )
  for ( tab in list( list( features, "feature" ), list( metadata, "metadata" ) ) ) {
    if ( !sample.column %in% names( tab[[ 1 ]] ) )
      stop( "The ", tab[[ 2 ]], " table has no column named '", sample.column,
            "'. Set sample.column to the name of the sample column.", call. = FALSE )
    s <- as.character( tab[[ 1 ]][[ sample.column ]] )
    if ( anyNA( s ) || any( s == "" ) )
      stop( "The ", tab[[ 2 ]], " table has empty sample names.", call. = FALSE )
    if ( anyDuplicated( s ) )
      stop( "Sample names repeated in the ", tab[[ 2 ]], " table: ",
            paste( utils::head( unique( s[ duplicated( s ) ] ), 5 ), collapse = ", " ),
            call. = FALSE )
  }

  f.samples <- as.character( features[[ sample.column ]] )
  m.samples <- as.character( metadata[[ sample.column ]] )
  only.f <- setdiff( f.samples, m.samples )
  only.m <- setdiff( m.samples, f.samples )
  if ( length( only.f ) )
    message( length( only.f ), " sample(s) in the feature table but not the metadata ",
             "were dropped: ", paste( utils::head( only.f, 5 ), collapse = ", " ),
             if ( length( only.f ) > 5 ) ", ..." )
  if ( length( only.m ) )
    message( length( only.m ), " sample(s) in the metadata but not the feature table ",
             "were dropped: ", paste( utils::head( only.m, 5 ), collapse = ", " ),
             if ( length( only.m ) > 5 ) ", ..." )
  samples <- f.samples[ f.samples %in% m.samples ]
  if ( length( samples ) == 0 )
    stop( "No sample names are shared by the two tables.", call. = FALSE )

  feat <- features[ match( samples, f.samples ), setdiff( names( features ), sample.column ),
                    drop = FALSE ]
  not.numeric <- names( feat )[ !vapply( feat, is.numeric, logical( 1 ) ) ]
  if ( length( not.numeric ) )
    stop( "Feature columns must be numeric; these are not: ",
          paste( utils::head( not.numeric, 10 ), collapse = ", " ),
          if ( length( not.numeric ) > 10 ) ", ...", ". Move metadata columns ",
          "to the metadata file.", call. = FALSE )
  X <- as.matrix( feat )
  storage.mode( X ) <- "double"
  rownames( X ) <- samples
  X[ !is.finite( X ) ] <- NA

  if ( transform == "log2" ) {
    if ( any( X < 0, na.rm = TRUE ) )
      stop( "transform = \"log2\" needs non-negative features; some are negative ",
            "(marker medians?). Use transform = \"none\".", call. = FALSE )
    X <- log2( X + offset )
  }

  meta <- metadata[ match( samples, m.samples ), , drop = FALSE ]
  rownames( meta ) <- samples

  structure( list( features = X, metadata = meta, sample.column = sample.column,
                   transform = transform, offset = offset ),
             class = "flowpredict.data" )
}


#' @export
print.flowpredict.data <- function( x, ... ) {
  cat( "FlowPredict data:", nrow( x$features ), "samples,", ncol( x$features ),
       "features", if ( x$transform == "log2" ) "(log2)" else "", "\n" )
  cat( "Metadata columns:", paste( names( x$metadata ), collapse = ", " ), "\n" )
  invisible( x )
}
