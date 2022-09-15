# Create a function to prepare for a predictive modeling
model_component=function(data
                         ,method
                         ,calib_method
                         ,tuning_trControl
                         ,final_trControl
                         ,tuningGrid
                         ,calib_tuningGrid){
  
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
  
  component
}