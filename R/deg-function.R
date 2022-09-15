# Create a function to define feature and outcome of a DEG
deg=function(the_dea,tissue,side=c('feature','outcome'),direction=F){
  eset=
    GSE73685 %>%
    .[,colnames(the_dea$eset)]
  
  if(side=='feature'){
    eset=
      eset %>%
      `assayData<-`(
        ExpressionSet(
          `dimnames<-`(
            normalize.quantiles.use.target(
              x=exprs(.)
              ,target=
                the_dea$result[rownames(.),'AveExpr1',drop=F] %>%
                pull(AveExpr1)
            )
            ,dimnames(.)
          )
        ) %>%
          assayData()
      )
  }
  
  eset %>%
    phenoData() %>%
    pData() %>%
    select(c('title',colnames(.) %>% .[str_detect(.,':ch1')])) %>%
    setNames(
      colnames(.) %>%
        str_remove_all(':ch1') %>%
        str_replace_all('\\s','_')
    ) %>%
    rownames_to_column(var='id') %>%
    filter(tissue_type==tissue) %>%
    select(subject_id,everything()) %>%
    arrange(subject_id) %>%
    lapply(X=1,Y=.,function(X,Y){
      GSE73685[,Y$id] %>%
        `colnames<-`(Y$subject_id) %>%
        exprs() %>%
        as.data.frame() %>%
        rownames_to_column(var='feature') %>%
        mutate(tissue_type=Y$tissue_type[1]) %>%
        gather(subject_id,value,-feature,-tissue_type)
    }) %>%
    .[[1]] %>%
    inner_join(
      the_dea$result %>%
        rownames_to_column(var='feature') %>%
        filter(adj.P.Val<0.05) %>%
        select(feature,AveExpr,AveExpr1,AveExpr2)
      ,by='feature'
    ) %>%
    mutate(type=side,regulation=direction) %>%
    mutate(
      outcome=
        ifelse(
          type=='feature'
          ,ifelse(
            AveExpr1<AveExpr
            ,ifelse(
              value<AveExpr
              ,0-(value-AveExpr)/(AveExpr1-AveExpr)
              ,0+(value-AveExpr)/(AveExpr2-AveExpr)
            )
            ,ifelse(
              value>AveExpr
              ,0-(value-AveExpr)/(AveExpr1-AveExpr)
              ,0+(value-AveExpr)/(AveExpr2-AveExpr)
            )
          )
          ,ifelse(
            regulation
            ,as.integer(value>=AveExpr1)
            ,ifelse(
              AveExpr1<AveExpr
              ,as.integer(value>=AveExpr)
              ,as.integer(value<AveExpr)
            )
          )
        )
    ) %>%
    select(-value,-type,-regulation,-AveExpr,-AveExpr1,-AveExpr2) %>%
    spread(feature,outcome)
}