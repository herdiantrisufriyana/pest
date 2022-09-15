# Create a function to preprocess the features
mb_predictor=function(exprs_mat,the_dea_list1){
  
  the_dea=the_dea_list1[[1]]
  the_tissue=names(the_dea_list1)[1]
  
  exprs_mat %>%
    dea_qn(the_dea$result) %>%
    as.data.frame() %>%
    rownames_to_column(var='feature') %>%
    left_join(
      the_dea$result %>%
        rownames_to_column(var='feature') %>%
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
    setNames(paste0(the_tissue,':',colnames(.)))
}