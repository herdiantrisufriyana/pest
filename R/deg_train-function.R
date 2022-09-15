# Create a function to train a surrogate transcriptome model
deg_train=function(deg
                   ,method
                   ,calib_method
                   ,tuning_trControl
                   ,final_trControl
                   ,tuningGrid
                   ,path
                   ,GSE73685=GSE73685
                   ,outcome=outcome
                   ,features=features
                   ,sample.kind=sample.kind
                   ,rsdr=rsdr
                   ,transformator=transformator
                   ,transformation=transformation){
  
  # component=
  #   list(
  #     deg='BABAM2'
  #     ,method='glmnet'
  #     ,calib_method='gamLoess'
  #     ,training_parameters=training_parameters
  #     ,tuningGrid=expand.grid(alpha=seq(0,1,len=5),lambda=10^seq(-9,0,len=5))
  #     ,path='data/deg_models/'
  #   )
  
  component=list()
  
  component$deg=deg
  component$method=method
  component$calib_method=calib_method
  component$training_parameters=
    list(
      tuning_trControl=tuning_trControl
      ,final_trControl=final_trControl
    )
  component$tuning=tuningGrid
  component$path=path
  
  component$data=
    GSE73685 %>%
    phenoData() %>%
    pData() %>%
    rownames_to_column(var='id') %>%
    select(id,title) %>%
    separate(title,c('tissue_type','subject_id'),sep='_') %>%
    select(-subject_id) %>%
    right_join(
      features %>%
        rownames_to_column(var='id') %>%
        left_join(
          outcome[,component$deg,drop=F] %>%
            setNames('outcome') %>%
            rownames_to_column(var='id')
          ,by='id'
        ) %>%
        mutate(deg=component$deg) %>%
        select(id,outcome,everything())
      ,by='id'
    ) %>%
    mutate(
      outcome=
        case_when(
          outcome==0~'nonevent'
          ,outcome==1~'event'
          ,TRUE~'missing'
        )
    )
  
  component$class=
    component$data %>%
    select(outcome,tissue_type) %>%
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
        select(tissue_type,outcome,weight)
      ,by=c('tissue_type','outcome')
    ) %>%
    select(tissue_type,outcome,weight,everything()) %>%
    column_to_rownames(var='id')
  
  component$train=
    component$data %>%
    rownames_to_column(var='id') %>%
    filter(outcome!='missing')
  
  component$features=
    suppressWarnings(
      component$train %>%
        select(-id) %>%
        gather(feature,value,-tissue_type,-outcome,-weight) %>%
        group_by(feature,outcome,tissue_type) %>%
        summarize(sd=sd(value,na.rm=T)) %>%
        ungroup() %>%
        filter(!(is.na(sd) | sd==0)) %>%
        filter(
          feature %in%
            pull(
              group_by(.,feature,outcome) %>%
                summarize(n=n()) %>%
                spread(outcome,n) %>%
                # filter(event==7 & nonevent==7)
                filter(event>=1 & nonevent>=1)
              ,feature
            )
        ) %>%
        mutate(
          tissue_type=
            tissue_type %>%
            str_replace_all('\\s+','_') %>%
            str_to_lower()
        ) %>%
        # unite(strata,outcome,tissue_type,sep=':') %>%
        # spread(strata,sd)
        spread(outcome,sd)
    )
  
  component$train=
    component$train %>%
    mutate(outcome=factor(outcome,c('event','nonevent'))) %>%
    column_to_rownames(var='id') %>%
    select_at(
      c('tissue_type','outcome','weight'
        ,component$features$feature %>% .[!duplicated(.)])
    ) %>%
    .[,sapply(.,function(x)sum(is.na(x))==0)]
  
  component$features=
    component$features %>%
    filter(feature%in%colnames(component$train))
  
  component$train_tissue_type=
    component$train %>%
    select(outcome,tissue_type) %>%
    table() %>%
    as.data.frame() %>%
    spread(outcome,Freq) %>%
    filter(event>=3 & nonevent>=3) %>%
    pull(tissue_type) %>%
    str_replace_all('\\s+','_') %>%
    str_to_lower()
  
  if(length(component$features$feature %>% .[!duplicated(.)])>1
     & length(component$train_tissue_type)>=1
     & floor(min(table(component$train$outcome))/3)>=2){
    
    component$transformator=
      component$train %>%
      rownames_to_column(var='id') %>%
      select(-tissue_type,-weight) %>%
      column_to_rownames(var='id') %>%
      select(-outcome) %>%
      rsdr(sample.kind=sample.kind) %>%
      transformator()
    
    component$pc=
      ExpressionSet(
        assayData=
          component$train %>%
          select(-tissue_type,-outcome,-weight) %>%
          t()
        ,phenoData=
          component$train %>%
          select(tissue_type,outcome,weight) %>%
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
          .[1:floor(min(table(component$train$outcome))/3)] %>%
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
          select(tissue_type,outcome,weight) %>%
          rownames_to_column(var='id')
        ,by='id'
      ) %>%
      select(outcome,tissue_type,weight,everything()) %>%
      column_to_rownames(var='id')
    
    component$model=
      suppressWarnings(caret::train(
        outcome~.
        ,data=
          component$pc %>%
          rownames_to_column(var='id') %>%
          select(-tissue_type,-weight) %>%
          column_to_rownames(var='id')
        ,method=component$method
        ,weights=
          component$pc %>%
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
        ,data=
          component$pc %>%
          rownames_to_column(var='id') %>%
          select(-tissue_type,-weight) %>%
          column_to_rownames(var='id')
        ,method=component$method
        ,weights=
          component$pc %>%
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
      component$train_tissue_type %>%
      lapply(function(x){
        x=x %>%
          str_replace_all('_',' ') %>%
          str_to_title()
        
        y=component$pc %>%
          rownames_to_column(var='id') %>%
          filter(tissue_type==x) %>%
          select(-tissue_type,-weight) %>%
          column_to_rownames(var='id') %>%
          cbind(predict(component$model,newdata=.,type='prob')) %>%
          rownames_to_column(var='id') %>%
          select(id,outcome,event) %>%
          rename(prob=event) %>%
          column_to_rownames(var='id')
        
        if(sd(y$prob)==0){
          NULL
        }else{
          suppressWarnings(caret::train(
            outcome~.
            ,data=y
            ,method=component$calib_method
            ,weights=
              component$pc %>%
              rownames_to_column(var='id') %>%
              filter(tissue_type==x) %>%
              select(id,weight) %>%
              column_to_rownames(var='id') %>%
              lapply(X=1,Y=.,function(X,Y){
                Y %>%
                  pull(weight) %>%
                  setNames(rownames(Y))
              }) %>%
              .[[1]]
            ,metric='ROC'
            ,trControl=component$training_parameters$final_trControl
            ,family=binomial(link='logit')
          ))
        }
      }) %>%
      setNames(component$train_tissue_type)
    
    component$eval_model=
      component$train_tissue_type %>%
      lapply(function(x){
        if(is.null(component$calib_model[[x]])){
          NULL
        }else{
          x2=x %>%
            str_replace_all('_',' ') %>%
            str_to_title()
          
          component$pc %>%
            rownames_to_column(var='id') %>%
            filter(tissue_type==x2) %>%
            select(-tissue_type,-weight) %>%
            column_to_rownames(var='id') %>%
            cbind(predict(component$model,newdata=.,type='prob')) %>%
            rownames_to_column(var='id') %>%
            select(id,outcome,event) %>%
            rename(prob=event) %>%
            column_to_rownames(var='id') %>%
            cbind(
              suppressWarnings(
                predict(component$calib_model[[x]],newdata=.,type='prob')
              ) %>%
                rename_all(function(x)paste0(x,'2'))
            ) %>%
            select(-prob) %>%
            rename_all(function(x)str_remove_all(x,'2')) %>%
            mutate(obs=outcome) %>%
            select(-outcome) %>%
            evalm(silent=T,showplots=F)
        }
      }) %>%
      setNames(component$train_tissue_type)
    
    if(sum(component$eval_model %>% sapply(function(x)!is.null(x)))>0){
      
      saveRDS(component,paste0(component$path,component$deg,'.rds'))
      
      rm(component)
      gc()
      
      deg
    }else{
      rm(component)
      gc()
      
      NULL
    }
  }else{
    rm(component)
    gc()
    
    NULL
  }
  
}