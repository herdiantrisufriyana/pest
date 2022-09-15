# Create a function to compare results of the evaluation
result_comparison=function(metric){
  list(
    all_GSE85307=opt_pred_disease1
    ,all_GSE86200=opt_pred_disease2
    ,all_GSE149437=opt_pred_disease3
    
    ,blood_GSE85307=opt_pred2_disease1
    ,blood_GSE86200=opt_pred2_disease2
    ,blood_GSE149437=opt_pred2_disease3
    
    ,repro_GSE85307=opt_pred3_disease1
    ,repro_GSE86200=opt_pred3_disease2
    ,repro_GSE149437=opt_pred3_disease3
  ) %>%
    lapply(X=names(.),Y=.,Z=metric,function(X,Y,Z){
      Y[[X]] %>%
        lapply(X=names(.),Y=.,Z=X,K=Z,function(X,Y,Z,K){
          Y[[X]]$all_conclusions[K,] %>%
            separate(CI_train,c('LB_train','UB_train'),sep='-') %>%
            separate(CI_test,c('LB_test','UB_test'),sep='-') %>%
            mutate_all(as.numeric) %>%
            gather(AUC,value) %>%
            separate(AUC,c('AUC','set'),sep='_') %>%
            spread(AUC,value) %>%
            mutate(model=X) %>%
            mutate(
              set=
                ifelse(
                  set=='train'
                  ,'GSE108497'
                  ,str_split_fixed(Z,'_',2)[,2]
                )
              ,tissue=str_split_fixed(Z,'_',2)[,1]
            )
        }) %>%
        do.call(rbind,.)
    }) %>%
    do.call(rbind,.) %>%
    filter(!duplicated(.)) %>%
    mutate(
      model=
        model %>%
        factor(c('pc_elnet','pc_rf','pc_gbm','divnn'))
      ,tissue=
        tissue %>%
        factor(c('all','blood','repro'))
      ,set=
        set %>%
        factor(c('GSE108497','GSE85307','GSE86200','GSE149437'))
    )
}