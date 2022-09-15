# DI-VNN feature selection and representation
divnn_pre_object=function(data,save=T,path=NULL){
  # Set input list
  input=list()
  
  # Use provider-wise, pre-calibrated, training set for modeling.
  input$data=data
  
  # Set output list
  output=list()
  
  # Get uncensored outcome.
  output$outcome=
    input$data %>%
    select(outcome)
  
  # Get features corresponding to uncensored outcome.
  output$raw=
    input$data %>%
    select(-outcome) %>%
    as.matrix()
  
  # Check if feature selection was conducted utilizing differential analysis
  # by moderated-t stats with Benjamini-Hochberg multiple testing correction.
  output$fit=
    t(output$raw) %>%
    normalize.quantiles() %>%
    `dimnames<-`(dimnames(t(output$raw))) %>%
    lmFit(model.matrix(~outcome,output$outcome)) %>%
    eBayes() %>%
    topTable(coef=2,nrow(.$coefficients),adjust.method='BH',sort.by='none') %>%
    filter(adj.P.Val<0.05)
  
  # Get the selected features as unnormalized ones.
  output$unnorm=
    output$raw %>%
    .[,rownames(output$fit)]
  
  # Normalize features quantile-to-quantile
  # using differential average of features as target.
  output$norm=
    output$unnorm %>%
    t() %>%
    normalize.quantiles.use.target(output$fit$AveExpr) %>%
    t() %>%
    `dimnames<-`(dimnames(output$unnorm))
  
  # Transform features by 1-bit stochastic gradient descent (SGD)
  # using the differential average.
  cat('Transform features by 1-bit stochastic gradient descent (SGD)\n')
  output$predictor_v=
    output$norm %>%
    sweep(2,output$fit$AveExpr,'-') %>%
    pbsapply(function(x){ifelse(x==0,0,ifelse(x>0,1,-1))}) %>%
    matrix(
      ncol=ncol(output$unnorm)
      ,byrow=F
      ,dimnames=dimnames(output$unnorm)
    ) %>%
    as.data.frame()
  
  # Get feature-wise mean of unnormalized features.
  output$predictor_m=
    output$unnorm %>%
    matrixStats::colMeans2()
  
  # Get feature-wise standard deviation (SD) of unnormalized features.
  output$predictor_s=
    output$unnorm %>%
    matrixStats::colSds()
  
  # Scale unnormalized features by the mean and SD.
  output$scaled=
    output$unnorm %>%
    sweep(2,output$predictor_m,'-') %>%
    sweep(2,output$predictor_s,'/')
  
  # Compute feature-feature Pearson correlation matrix
  # of unnormalized features.
  cat('Compute feature-feature Pearson correlation matrix\n')
  output$predictor_p=
    colnames(output$unnorm) %>%
    lapply(X=seq(length(.)-1),Y=.,function(X,Y){
      data.frame(
        predictor1=Y[X]
        ,predictor2=Y[(X+1):length(.)]
      )
    }) %>%
    do.call(rbind,.) %>%
    mutate(
      pearson=
        pbsapply(X=seq(nrow(.)),Y=.,Z=output$scaled,function(X,Y,Z){
          cor(
            Z[,Y$predictor1[X]]
            ,Z[,Y$predictor2[X]]
          )
        })
    ) %>%
    rbind(
      setNames(select(.,predictor2,predictor1,everything()),colnames(.))
      ,data.frame(
        predictor1=colnames(output$unnorm)
        ,predictor2=colnames(output$unnorm)
        ,pearson=1
      )
    ) %>%
    spread(predictor2,pearson) %>%
    column_to_rownames(var='predictor1') %>%
    as.matrix()
  
  # Conduct Barnes-Hut t-moderated stochastic neighbor embedding (t-SNE).
  cat('Conduct Barnes-Hut t-moderated stochastic neighbor embedding (t-SNE)\n')
  suppressWarnings(set.seed(33,sample.kind=sample.kind))
  output$predictor_tsne=
    output$predictor_p %>%
    Rtsne(dims=3,verbose=T,is_distance=T)
  
  # Construct clique-extracted ontology (CliXO) (choose your own OS).
  output$predictor_c=
    output$predictor_p %>%
    clixo(os='windows')
  
  # Compile input.
  input=list()
  
  input$value=output$predictor_v
  
  input$outcome=
    output$outcome %>%
    mutate(outcome=as.integer(outcome=='event')) %>%
    pull(outcome) %>%
    setNames(rownames(input$value))
  
  input$similarity=output$predictor_p
  
  input$mapping=
    output$predictor_tsne$Y %>%
    `rownames<-`(colnames(input$value))
  
  input$ontology=
    output$predictor_c %>%
    mutate(seq=seq(nrow(.))) %>%
    gather(key,value,-similarity,-relation,-seq) %>%
    mutate(
      value2=
        ifelse(
          str_detect(value,'CliXO:')
          ,str_remove_all(value,'CliXO:')
          ,NA
        ) %>%
        as.integer()
    ) %>%
    mutate(
      value=
        ifelse(
          str_detect(value,'CliXO:')
          ,paste0(
            'ONT:'
            ,str_pad(value2,str_count(max(value2,na.rm=T)),'left','0')
          )
          ,value
        )
    ) %>%
    select(-value2) %>%
    spread(key,value) %>%
    arrange(seq) %>%
    select(-seq) %>%
    select(source,target,everything())
  
  input$fit=output$fit
  
  # Compile into a TidySet.
  cat('Compile into a TidySet\n')
  output=
    TidySet.compile(
      value=input$value
      ,outcome=input$outcome
      ,similarity=input$similarity
      ,mapping=input$mapping
      ,ontology=input$ontology
      ,ranked=T
      ,dims=7
      ,decreasing=F
      ,seed_num=33
    )
  
  if(save){
    saveRDS(input,paste0(path,'_input.rds'))
    saveRDS(output,paste0(path,'_output.rds'))
  }
  
  list(input=input,output=output)
}