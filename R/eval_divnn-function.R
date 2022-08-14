#' DI-VNN evaluation for classification task in PROM study
#'
#' Evaluate deep-insight visible neural network for classification task in this 
#' PROM study.
#'
#' @param data Whole set, a list of ExpressionSet that includes training, 
#' calibration, and testing sets.
#' @param modeling DI-VNN modeling results, a list of Keras model object, 
#' early-stopping callback object, training index vector, validation index 
#' vector, a list of history of the Keras training, and a list of resampled 
#' evaluation results.
#' @param batch_size Batch size, an integer of the number of instances used 
#' for each step or iteration.
#' @param calib_model Trained caret model for calibrating, a caret object of 
#' the output calibrating model of DI-VNN in \code{divnn_calibrator_evaluator} 
#' in which this function is used.
#'
#' @details
#' This function specifically works for a DI-VNN pipeline in this PROM study.
#' 
#' @return Evaluation results of DI-VNN, an evalm object of the results.
#' 

eval_divnn=function(data,modeling,batch_size=128,calib_model){
  
  # Create a list of character to call testing sets.
  list('int_val','test_set1','test_set2','test_set3') %>%
    lapply(function(x){
      
      # Evaluate non-calibrated DI-VNN for this subset.
      cat('Evaluate non-calibrated DI-VNN for set of:',x,'\n')
      modeling$ontonet %>%
        
        # Predict the outcome.
        predict_generator(
          generator=
            generator.ontoarray(
              data[[x]]
              ,seq(ncol(data[[x]]))
              ,batch_size=batch_size
            )
          ,steps=ceiling(ncol(data[[x]])/batch_size)
          ,verbose=1
        ) %>%
        
        # Give back the name of the DI-VNN ontologies.
        setNames(
          modeling$history$metrics %>%
            names() %>%
            .[str_detect(.,'val_')
              &str_detect(.,'_loss')
              &str_detect(.,'ONT|root')] %>%
            str_remove_all('val_|_loss')
        ) %>%
        
        # Set back number of instances, that is changed due to batch sizing.
        lapply(function(y){
          y[seq(ncol(data[[x]])),,drop=F]
        }) %>%
        
        # Combine with the outcome data.
        c(list(outcome=data[[x]]$outcome)) %>%
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
        mutate(nonevent=1-event) %>%
        select(nonevent,everything()) %>%
        `rownames<-`(colnames(data[[x]])) %>%
        select(nonevent,event,outcome) %>%
        mutate(obs=outcome) %>%
        select(-outcome) %>%
        
        # Evaluate the model.
        evalm(silent=T,showplots=F)
      
    }) %>%
    setNames(
      paste0('nocalib_',c('int_val','test_set1','test_set2','test_set3'))
    ) %>%
    c(list('int_val','test_set1','test_set2','test_set3') %>%
        lapply(function(x){
          
          # Evaluate calibrated DI-VNN for this subset.
          cat('Evaluate calibrated DI-VNN for set of:',x,'\n')
          modeling$ontonet %>%
            
            # Predict the outcome.
            predict_generator(
              generator=
                generator.ontoarray(
                  data[[x]]
                  ,seq(ncol(data[[x]]))
                  ,batch_size=batch_size
                )
              ,steps=ceiling(ncol(data[[x]])/batch_size)
              ,verbose=1
            ) %>%
            
            # Give back the name of the DI-VNN ontologies.
            setNames(
              modeling$history$metrics %>%
                names() %>%
                .[str_detect(.,'val_')
                  &str_detect(.,'_loss')
                  &str_detect(.,'ONT|root')] %>%
                str_remove_all('val_|_loss')
            ) %>%
            
            # Set back number of instances, that is changed due to batch sizing.
            lapply(function(y){
              y[seq(ncol(data[[x]])),,drop=F]
            }) %>%
            
            # Combine with the outcome data.
            c(list(outcome=data[[x]]$outcome)) %>%
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
            mutate(nonevent=1-event) %>%
            select(nonevent,everything()) %>%
            `rownames<-`(colnames(data[[x]])) %>%
            rownames_to_column(var='id') %>%
            select(id,outcome,event) %>%
            rename(prob=event) %>%
            column_to_rownames(var='id') %>%
            cbind(
              suppressWarnings(predict(calib_model,newdata=.,type='prob')) %>%
                rename_all(function(x)paste0(x,'2'))
            ) %>%
            select(-prob) %>%
            rename_all(function(x)str_remove_all(x,'2')) %>%
            mutate(obs=outcome) %>%
            select(-outcome) %>%
            select(nonevent,everything()) %>%
            
            # Evaluate the model.
            evalm(silent=T,showplots=F)
          
        }) %>%
        setNames(
          paste0('calib_',c('int_val','test_set1','test_set2','test_set3'))
        )
    )
}