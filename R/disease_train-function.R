# Create a function to train a non-DIVNN prediction model
disease_train=function(data
                       ,method
                       ,calib_method
                       ,tuning_trControl
                       ,final_trControl
                       ,tuningGrid
                       ,calib_tuningGrid
                       ,epv=20
                       ,sample.kind=sample.kind){
  
  component=list()
  
  component$data=data
  component$method=method
  component$calib_method=calib_method
  component$training_parameters=
    list(
      tuning_trControl=tuning_trControl
      ,final_trControl=final_trControl
    )
  component$tuning=tuningGrid
  component$calib_tuning=calib_tuningGrid
  
  suppressWarnings(set.seed(33,sample.kind=sample.kind))
  component$data=
    component$data %>%
    rownames_to_column(var='id') %>%
    mutate(
      partition=
        suppressWarnings(
          sapply(X=seq(nrow(.)),Y=outcome,function(X,Y){
            Z=createDataPartition(Y,p=0.8)[[1]]
            ifelse(X%in%Z,'train','calib')
          })
        )
    ) %>%
    select(id,partition,everything())
  
  component$class=
    component$data %>%
    select(outcome,partition) %>%
    table() %>%
    as.data.frame() %>%
    rename(n=Freq) %>%
    left_join(
      group_by(.,outcome) %>%
        summarize(subtotal=sum(n))
      ,by='outcome'
    ) %>%
    cbind(summarize(.,total=sum(n))) %>%
    mutate(
      weight=
        ifelse(
          outcome=='missing'
          ,NA
          ,1/(subtotal/total)*0.5
        )
    )
  
  component$data=
    component$data %>%
    left_join(
      component$class %>%
        select(partition,outcome,weight)
      ,by=c('partition','outcome')
    ) %>%
    select(partition,outcome,weight,everything()) %>%
    column_to_rownames(var='id')
  
  component$train=
    component$data %>%
    rownames_to_column(var='id') %>%
    filter(partition=='train' & outcome!='missing') %>%
    select(-partition) %>%
    column_to_rownames(var='id')
  
  component$calib=
    component$data %>%
    rownames_to_column(var='id') %>%
    filter(partition=='calib' & outcome!='missing') %>%
    select(-partition) %>%
    column_to_rownames(var='id')
  
  component$transformator=
    component$train %>%
    rownames_to_column(var='id') %>%
    select(-weight) %>%
    column_to_rownames(var='id') %>%
    select(-outcome) %>%
    rsdr(sample.kind=sample.kind) %>%
    transformator()
  
  component$pc=
    ExpressionSet(
      assayData=
        component$train %>%
        select(-outcome,-weight) %>%
        t()
      ,phenoData=
        component$train %>%
        select(outcome,weight) %>%
        AnnotatedDataFrame()
    ) %>%
    .[rownames(component$transformator$scaler),] %>%
    transformation(transformator_object=component$transformator)
  
  component$pc=
    component$pc %>%
    .[component$pc %>%
        preproc() %>%
        .$pve %>%
        sort(decreasing=T) %>%
        .[1:floor(min(table(component$train$outcome))/epv)] %>%
        names()
      ,] %>%
    exprs() %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var='id') %>%
    left_join(
      component$pc %>%
        phenoData() %>%
        pData() %>%
        select(outcome,weight) %>%
        rownames_to_column(var='id')
      ,by='id'
    ) %>%
    select(outcome,weight,everything()) %>%
    column_to_rownames(var='id')
  
  component$pc2=
    ExpressionSet(
      assayData=
        component$calib %>%
        select(-outcome,-weight) %>%
        t()
      ,phenoData=
        component$calib %>%
        select(outcome,weight) %>%
        AnnotatedDataFrame()
    ) %>%
    .[rownames(component$transformator$scaler),] %>%
    transformation(transformator_object=component$transformator)
  
  component$pc2=
    component$pc2 %>%
    .[component$pc %>%
        colnames() %>%
        .[!.%in%c('outcome','weight')]
      ,] %>%
    exprs() %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var='id') %>%
    left_join(
      component$pc2 %>%
        phenoData() %>%
        pData() %>%
        select(outcome,weight) %>%
        rownames_to_column(var='id')
      ,by='id'
    ) %>%
    select(outcome,weight,everything()) %>%
    column_to_rownames(var='id')
  
  data=
    component$pc %>%
    rownames_to_column(var='id') %>%
    select(-weight) %>%
    column_to_rownames(var='id')
  
  data2=
    component$pc2 %>%
    rownames_to_column(var='id') %>%
    select(-weight) %>%
    column_to_rownames(var='id')
  
  component$model=
    suppressWarnings(caret::train(
      outcome~.
      ,data=data
      ,method=component$method
      ,weights=
        component$train %>%
        rownames_to_column(var='id') %>%
        select(id,weight) %>%
        column_to_rownames(var='id') %>%
        pull(weight) %>%
        setNames(rownames(component$train))
      ,metric='ROC'
      ,trControl=component$training_parameters$tuning_trControl
      ,tuneGrid=component$tuning
    ))
  
  component$model=
    suppressWarnings(caret::train(
      outcome~.
      ,data=data
      ,method=component$method
      ,weights=
        component$train %>%
        rownames_to_column(var='id') %>%
        select(id,weight) %>%
        column_to_rownames(var='id') %>%
        pull(weight) %>%
        setNames(rownames(component$train))
      ,metric='ROC'
      ,trControl=component$training_parameters$final_trControl
      ,tuneGrid=component$model$bestTune
    ))
  
  component$calib_model=
    suppressWarnings(caret::train(
      outcome~.
      ,data=
        data2 %>%
        cbind(predict(component$model,newdata=.,type='prob')) %>%
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
  
  component$pc3=
    ExpressionSet(
      assayData=
        component$data[,colnames(component$train)] %>%
        rownames_to_column(var='id') %>%
        filter(outcome!='missing') %>%
        select(-weight) %>%
        column_to_rownames(var='id') %>%
        select(-outcome) %>%
        t()
      ,phenoData=
        component$data[,colnames(component$train)] %>%
        select(outcome,weight) %>%
        AnnotatedDataFrame()
    ) %>%
    .[rownames(component$transformator$scaler),] %>%
    transformation(transformator_object=component$transformator)
  
  component$pc3=
    component$pc3 %>%
    .[component$pc %>%
        colnames() %>%
        .[!.%in%c('outcome','weight')]
      ,] %>%
    exprs() %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var='id') %>%
    left_join(
      component$pc3 %>%
        phenoData() %>%
        pData() %>%
        select(outcome,weight) %>%
        rownames_to_column(var='id')
      ,by='id'
    ) %>%
    select(outcome,weight,everything()) %>%
    column_to_rownames(var='id')
  
  component$eval_model=
    component$pc3 %>%
    rownames_to_column(var='id') %>%
    filter(outcome!='missing') %>%
    select(-weight) %>%
    column_to_rownames(var='id') %>%
    cbind(predict(component$model,newdata=.,type='prob')) %>%
    rownames_to_column(var='id') %>%
    select(id,outcome,event) %>%
    rename(prob=event) %>%
    column_to_rownames(var='id') %>%
    cbind(
      suppressWarnings(predict(component$calib_model,newdata=.,type='prob')) %>%
        rename_all(function(x)paste0(x,'2'))
    ) %>%
    select(-prob) %>%
    rename_all(function(x)str_remove_all(x,'2')) %>%
    mutate(obs=outcome) %>%
    select(-outcome) %>%
    select(nonevent,everything()) %>%
    evalm(silent=T,showplots=F)
  
  rm(data,data2)
  
  component
}