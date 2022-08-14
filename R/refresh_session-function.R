#' Refresh keras backend session
#'
#' Refresh keras backend session for tensorflow 2.0.0.
#'
#' @details
#' None.
#' 
#' @return None.
#' 

refresh_session=function(){
  suppressWarnings(set.seed(33,sample.kind=sample.kind))
  library(tensorflow)
  library(keras)
  k_clear_session()
  tf$random$set_seed(33)
  gpu_options=
    tf$compat$v1$GPUOptions(
      allow_growth=TRUE,
      per_process_gpu_memory_fraction=1,
      visible_device_list='0'
    )
  config=tf$compat$v1$ConfigProto(gpu_options=gpu_options)
  session=tf$compat$v1$Session(config=config)
  graph=tf$compat$v1$get_default_graph
  tf$compat$v1$keras$backend$set_session(session)
}