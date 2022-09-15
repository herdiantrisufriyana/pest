#' Transform test data for DI-VNN input
#'
#' Transform test data for deep-insight visible neural network input.
#'
#' @param test_data Tidy set, an ExpressionSet that includes testing set.
#' @param SGD1bit_fit Differential anayslis results, a dataframe of the 
#' differential analysis results that are already filtered for the significant 
#' predictors. This will be used for 1-bit stochastic gradient descent 
#' transformation.
#' @param similarity Feature similarity, a square matrix of numerics containing
#' feature-feature similarity measures.
#' @param mapping Feature three-dimensional mapping, a matrix of numerics with
#' rows for features and three columns for three dimensions where the features
#' are mapped onto.
#' @param ontology Ontology, a data frame with rows for ontologies and four
#' columns for source, target, similarity, and relation. Feature (source)-
#' ontology (target) relation should be annotated as 'feature', while ontology-
#' ontology relation should be annotated as 'is_a'. To differentiate between
#' feature and ontology names, a prefix of 'ONT:' precedes an ontology name. All
#' columns except similarity in \code{ontology} should be characters. Similarity
#' (a numeric) is a minimum threshold by which either features or ontologies
#' (source) belong to an ontology (target).
#' @param ranked Map ranking, a logical whether mapping uses rank instead the 
#' real values.
#' @param dims Dimension, an integer of pixelation for the square ontology map.
#' @param decreasing Order of angle, a logical of whether an optimal angle 
#' found by decreasing order of the number of channel. Decreasing means the 
#' angle that create the highest number of channels.
#' @param seed_num Random seed, an integer of seed when randomly taking one of 
#' multiple optimum angles with the same number of channels.
#'
#' @details
#' This function works for a DI-VNN pipeline with specific steps.
#' 
#' @return output TidySet, an ExpressionSet with three tables. Instance-feature
#' value data frame and outcome vector are compiled as a phenotype data frame
#' with rows for instances and columns for features and outcome. Instance-
#' feature value data frame and feature three-dimensional mapping matrix are
#' compiled as an expression matrix with rows for positions of features and
#' columns for instances. The mapping and similarity matrices and ontology data
#' frame are compiled as a feature data frame with rows for positions of
#' features and columns for feature names and ontological relations. For easier,
#' access the similarity matrix, ontomap four-dimensional array, ontotype list
#' of two-dimensional matrices, and ontology data frame are included in
#' experiment notes that can be called using Biobase function \code{notes}.
#' 

test_transformer=function(test_data
                          ,SGD1bit_fit
                          ,similarity
                          ,mapping
                          ,ontology
                          ,ranked=T
                          ,dims=7
                          ,decreasing=F
                          ,seed_num=33){
  
  # Get the outcome.
  outcome=
    test_data %>%
    select(outcome)
  
  # Get the features.
  raw=
    test_data %>%
    select(-outcome) %>%
    as.matrix() %>%
    .[rownames(outcome),,drop=F]
  
  # Filter the selected features.
  unnorm=
    raw %>%
    .[,rownames(SGD1bit_fit),drop=F]
  
  # Normalize the selected features quantile-to-quantile
  # using the averages based on the training set.
  norm=
    unnorm %>%
    t() %>%
    normalize.quantiles.use.target(SGD1bit_fit$AveExpr) %>%
    t() %>%
    `dimnames<-`(dimnames(unnorm))
  
  # Do 1-bit stochasti gradient descent transformation.
  cat('Transform features by 1-bit stochastic gradient descent (SGD)\n')
  value=
    norm %>%
    sweep(2,SGD1bit_fit$AveExpr,'-') %>%
    pbsapply(function(x){ifelse(x==0,0,ifelse(x>0,1,-1))}) %>%
    matrix(
      ncol=ncol(unnorm)
      ,byrow=F
      ,dimnames=dimnames(unnorm)
    ) %>%
    as.data.frame()
  
  # Compile into a tidy set.
  cat('Compile into a TidySet\n')
  TidySet.compile(
    value=value
    ,outcome=
      outcome %>%
      mutate(outcome=as.integer(outcome=='event')) %>%
      pull(outcome) %>%
      setNames(rownames(value))
    ,similarity=similarity
    ,mapping=mapping
    ,ontology=ontology
    ,ranked=ranked
    ,dims=dims
    ,decreasing=decreasing
    ,seed_num=seed_num
  )
  
}