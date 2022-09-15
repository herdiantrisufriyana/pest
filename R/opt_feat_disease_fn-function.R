# Create a function to get model weights
opt_feat_disease_fn=function(X,Y,Z){
  if(X=='divnn'){
    # Get the classification DI-VNN weights of the representation layers.
    Z$visualization$ontoarray %>%
      lapply(function(x){
        x$ontotype %>%
          lapply(X=seq(nrow(.)),Y=.,Z=x$output,function(X,Y,Z){
            K=Z[Y$x[X],Y$y[X],Y$z[X]]
            if(K==0){
              NULL
            }else{
              K
            }
          }) %>%
          setNames(x$ontotype$feature)
      }) %>%
      setNames(names(Z$visualization$ontoarray)) %>%
      unlist() %>%
      data.frame(output=.) %>%
      rownames_to_column(var='ontology') %>%
      rename(classification_output=output) %>%
      separate(ontology,c('ontology','predictor'),sep='\\.') %>%
      
      # Get the predictor position.
      left_join(
        Z$pre_object$output %>%
          fData() %>%
          rownames_to_column(var='xyz') %>%
          rename(predictor=feature) %>%
          select(xyz,predictor) %>%
          filter(!is.na(predictor)) %>%
          mutate(xyz=str_remove_all(xyz,'x')) %>%
          separate(xyz,c('dimension_1','yz'),sep='y') %>%
          separate(yz,c('dimension_2','dimension_3'),sep='z')
        ,by='predictor'
      ) %>%
      
      # Wrap up.
      arrange(ontology,dimension_3,dimension_2,dimension_1,predictor) %>%
      select(
        ontology
        ,dimension_1
        ,dimension_2
        ,dimension_3
        ,predictor
        ,everything()
      )
  }else{
    if(Y[[X]]$method=='glmnet'){
      # Get weights from the elastic net regression model.
      Z=coef(Y[[X]]$model$finalModel,Y[[X]]$model$bestTune$lambda) %>%
        as.matrix() %>%
        as.data.frame()
    }else{
      # Get weights from the tree-based models.
      Z=varImp(Y[[X]]$model) %>%
        .[[1]]
    }
    
    Z %>%
      rownames_to_column(var='term') %>%
      
      # Clean up.
      mutate(term=str_remove_all(term,'\\`|\\(|\\)')) %>%
      setNames(c('term','estimate')) %>%
      filter(estimate!=0) %>%
      arrange(desc(abs(estimate))) %>%
      filter(!term %in% c('Intercept')) %>%
      right_join(
        Y[[X]]$transformator$avg_rotm %>%
          .[,Y[[X]]$pc %>% colnames() %>% .[!.%in%c('outcome','weight')]] %>%
          as.data.frame() %>%
          rownames_to_column(var='tissue_feature') %>%
          separate(tissue_feature,c('tissue_type','feature'),sep=':') %>%
          gather(term,value,-feature,-tissue_type)
        ,by='term'
      ) %>%
      arrange(
        factor(tissue_type,c('maternal_blood',names(dea)))
        ,desc(estimate)
      )
  }
}