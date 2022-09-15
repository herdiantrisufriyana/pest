# Create a function to derive the surrogate transcriptome
surrogate_deg=function(exprs_mat,the_dea_list1,component){
  
  # exprs_mat=
  #   GSE73685 %>%
  #   .[,phenoData(dea$amnion$eset) %>%
  #       pData() %>%
  #       filter(str_detect(str_to_lower(title),'maternal blood')) %>%
  #       rownames()] %>%
  #   exprs()
  # the_dea_list1=dea['amnion']
  # component=readRDS('data/deg_models/DIRAS2.rds')
  
  the_dea=the_dea_list1[[1]]
  the_tissue=names(the_dea_list1)[1]
  
  exprs_mat %>%
    dea_qn(the_dea$result) %>%
    .[rownames(component$transformator$avg_rotm),] %>%
    as.data.frame() %>%
    rownames_to_column(var='feature') %>%
    left_join(
      the_dea$result %>%
        rownames_to_column(var='feature') %>%
        filter(adj.P.Val<0.05) %>%
        select(feature,AveExpr,AveExpr1,AveExpr2)
      ,by='feature'
    ) %>%
    gather(id,value,-feature,-AveExpr,-AveExpr1,-AveExpr2) %>%
    mutate(
      norm_value=
        ifelse(
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
    ) %>%
    select(-value,-AveExpr,-AveExpr1,-AveExpr2) %>%
    spread(feature,norm_value) %>%
    column_to_rownames(var='id') %>%
    t() %>%
    ExpressionSet(
      phenoData=
        t(.) %>%
        as.data.frame() %>%
        mutate(outcome=NA) %>%
        select(outcome) %>%
        AnnotatedDataFrame()
    ) %>%
    transformation(transformator_object=component$transformator) %>%
    exprs() %>%
    t() %>%
    as.data.frame() %>%
    cbind(predict(component$model,newdata=.,type='prob')) %>%
    select(event) %>%
    rename(prob=event) %>%
    cbind(suppressWarnings(
      predict(component$calib_model[[the_tissue]],newdata=.,type='prob')
    )) %>%
    rename(calib_prob=event) %>%
    select(calib_prob) %>%
    setNames(paste0(the_tissue,':',component$deg))
}