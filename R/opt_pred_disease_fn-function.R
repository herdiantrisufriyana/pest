# Create a function to evaluate the models
opt_pred_disease_fn=function(X,Y,Z,Z2){
  if(X=='divnn'){
    K=Y[[X]]$eval_model[[paste0('calib_',Z2)]]
    Y[[X]]$eval_model=Y[[X]]$eval_model$calib_int_val
  }else{
    K=Z %>%
      lapply(X=1,Y=.,Z=Y[[X]],function(X,Y,Z){
        K=ExpressionSet(
          assayData=
            Y %>%
            select(-outcome) %>%
            t()
          ,phenoData=
            Y %>%
            select(outcome) %>%
            AnnotatedDataFrame()
        ) %>%
          .[rownames(Z$transformator$scaler),] %>%
          transformation(transformator_object=Z$transformator)
        
        K %>%
          .[Z$pc %>%
              colnames() %>%
              .[!.%in%c('outcome','weight')]
            ,] %>%
          exprs() %>%
          t() %>%
          as.data.frame() %>%
          rownames_to_column(var='id') %>%
          left_join(
            K %>%
              phenoData() %>%
              pData() %>%
              select(outcome) %>%
              rownames_to_column(var='id')
            ,by='id'
          ) %>%
          select(outcome,everything()) %>%
          column_to_rownames(var='id')
      }) %>%
      .[[1]] %>%
      cbind(predict(Y[[X]]$model,newdata=.,type='prob')) %>%
      rownames_to_column(var='id') %>%
      select(id,outcome,event) %>%
      rename(prob=event) %>%
      column_to_rownames(var='id') %>%
      cbind(
        suppressWarnings(
          predict(Y[[X]]$calib_model,newdata=.,type='prob')
        ) %>%
          rename_all(function(x)paste0(x,'2'))
      ) %>%
      select(-prob) %>%
      rename_all(function(x)str_remove_all(x,'2')) %>%
      mutate(obs=outcome) %>%
      select(-outcome) %>%
      select(nonevent,event,obs) %>%
      evalm(silent=T,showplots=F)
  }
  
  L=list(
    Y[[X]]$eval_model$optres$Group1 %>%
      setNames(paste0(colnames(.),'_train'))
    
    ,K %>%
      .$optres %>%
      .$Group1 %>%
      setNames(paste0(colnames(.),'_test'))
  ) %>%
    do.call(cbind,.)
  
  list(all_conclusions=L,test_eval=K)
}