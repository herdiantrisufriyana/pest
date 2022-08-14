#' DI-VNN calibration and evaluation for classification task in PROM study
#'
#' Calibrate and evaluate a DI-VNN model in a single pipeline for 
#' classification task in this PROM study.
#'
#' @param data Whole set, a list of ExpressionSet that includes training, 
#' calibration, and testing sets.
#' @param modeling DI-VNN modeling results, a list of Keras model object, 
#' early-stopping callback object, training index vector, validation index 
#' vector, a list of history of the Keras training, and a list of resampled 
#' evaluation results.
#' @param training Training strategy, a dataframe consisting details on the 
#' training strategy, e.g., internal validation method, etc.
#' @param batch_size Batch size, an integer of the number of instances used 
#' for each step or iteration.
#' @param title Title for log, a character of the title.
#' @param dir Path for saving the results, a character of the folder path.
#' @param file File name prefix, a character of the file name prefix.
#'
#' @details
#' This function specifically works for a DI-VNN pipeline in this PROM study.
#' 
#' @return Various results, a list of objects for the model (caret-like model 
#' object), calibrating model (caret object), and evaluation results (evalm 
#' object).
#' 

divnn_calibrator_evaluator=function(data
                                    ,modeling
                                    ,component
                                    ,batch_size=128
                                    ,title
                                    ,dir='.'
                                    ,file=NULL){
  
  # Print message for reproducible logging.
  cat(title,'\n')
  cat('Started:',as.character(now()),'\n')
  
  # Get predicted probability of DI-VNN.
  cat('Get predicted probability of DI-VNN\n')
  refresh_session()
  prediction=
    modeling$ontonet %>%
    
    # Predict the outcome of the calibration set.
    predict_generator(
      generator=
        generator.ontoarray(
          data[['calib']]
          ,seq(ncol(data[['calib']]))
          ,batch_size=batch_size
        )
      ,steps=ceiling(ncol(data[['calib']])/batch_size)
      ,verbose=1
    ) %>%
    
    # Give back the name of the DI-VNN ontologies.
    setNames(
      modeling$history$metrics %>%
        names() %>%
        .[str_detect(.,'val_')
          &str_detect(.,'_loss')&str_detect(.,'ONT|root')] %>%
        str_remove_all('val_|_loss')
    ) %>%
    
    # Set back number of instances, that is changed due to batch sizing.
    lapply(function(x){
      x[seq(ncol(data[['calib']])),,drop=F]
    }) %>%
    
    # Combine with the outcome data.
    c(list(outcome=data[['calib']]$outcome)) %>%
    as.data.frame() %>%
    select(root,outcome) %>%
    mutate(
      outcome=
        factor(
          ifelse(outcome==1,'event','nonevent')
          ,c('nonevent','event')
        )
    ) %>%
    
    # Mimic the standard input of evalm.
    rename(event=root) %>%
    mutate(nonevent=event) %>%
    select(nonevent,everything()) %>%
    
    # Give back the rownames of the calibration set.
    `rownames<-`(colnames(data[['calib']]))
  
  # Create an empty list to create a caret-like list.
  model=list()
  model$results=
    modeling$evaluation %>%
    lapply(as.data.frame) %>%
    do.call(rbind,.) %>%
    select(root_roc) %>%
    summarize(ROC=mean(root_roc),ROCSD=sd(root_roc))
  model$pred=
    prediction %>%
    rename(obs=outcome) %>%
    select(obs,everything())
  if(!is.null(file)) saveRDS(model,paste0(dir,'/',file,'.rds'))
  
  # Model calibration
  cat('Calibrate DI-VNN\n')
  calib_model=
    suppressWarnings(caret::train(
      outcome~.
      ,data=
        prediction %>%
        rownames_to_column(var='id') %>%
        select(id,outcome,event) %>%
        rename(prob=event) %>%
        column_to_rownames(var='id')
      ,method=component$calib_method
      ,weights=
        component$calib %>%
        rownames_to_column(var='id') %>%
        select(id,weight) %>%
        column_to_rownames(var='id') %>%
        pull(weight) %>%
        setNames(rownames(component$calib))
      ,metric='ROC'
      ,trControl=component$training_parameters$final_trControl
      ,tuneGrid=component$calib_tuning
      ,family=binomial(link='logit')
    ))
  if(!is.null(file)) saveRDS(calib_model,paste0(dir,'/calib_',file,'.rds'))
  
  # Model evaluation on validation sets
  eval_model=
    data %>%
    eval_divnn(
      modeling=modeling
      ,batch_size=batch_size
      ,calib_model=calib_model
    )
  if(!is.null(file)) saveRDS(eval_model,paste0(dir,'/eval_',file,'.rds'))
  
  # Clean up memory.
  gc()
  
  # Log off.
  cat('End:',as.character(now()))
  
  list(model=model,calib_model=calib_model,eval_model=eval_model)
  
}