# Create a function to evaluate leaves of a biomarker tree
tree_eval=function(data,train_mod,train_eval){
  suppressWarnings(set.seed(33,sample.kind=sample.kind))
  train_mod %>%
    predict(newdata=data,type='prob') %>%
    rownames_to_column(var='id') %>%
    left_join(
      train_mod$finalModel$where %>%
        data.frame(node=.) %>%
        rownames_to_column(var='id') %>%
        left_join(
          data.frame(
            node=unique(.$node)
            ,leaf=LETTERS[seq(length(unique(.$node)))]
          )
          ,by='node'
        ) %>%
        right_join(
          train_mod %>%
            predict(type='prob') %>%
            rownames_to_column(var='id')
          ,by='id'
        ) %>%
        select(event,nonevent,node,leaf) %>%
        filter(!duplicated(.))
      ,by=c('event','nonevent')
    ) %>%
    mutate(th=train_eval$th[1]) %>%
    mutate(pred=ifelse(event>=th,'event','nonevent')) %>%
    cbind(select(data,outcome)) %>%
    lapply(X=seq(30),Y=.,function(X,Y){
      Y %>%
        .[sample(seq(nrow(.)),nrow(.),T),] %>%
        group_by(th,leaf,pred,outcome) %>%
        summarize(n=n()) %>%
        ungroup() %>%
        group_by(th,leaf,pred) %>%
        summarize(
          true_pred=sum(ifelse(outcome==pred,n,NA),na.rm=T)
          ,false_pred=sum(ifelse(outcome==pred,NA,n),na.rm=T)
          ,leaf_total=sum(n)
          ,leaf_ppv_npv=true_pred/leaf_total
        ) %>%
        ungroup() %>%
        group_by(th,pred) %>%
        mutate(
          pred_total=sum(leaf_total)
          ,leaf_p=leaf_total/pred_total
          ,pred_ppv_npv=sum(true_pred)/pred_total
        ) %>%
        ungroup() %>%
        mutate(
          total=sum(leaf_total)
          ,pred_p=sum(ifelse(pred=='event',leaf_total,NA),na.rm=T)/total
          ,outcome_p=sum(ifelse(pred=='event',true_pred,false_pred))/total
        ) %>%
        gather(metric,value,-th,-leaf,-pred) %>%
        arrange(th,leaf,pred) %>%
        mutate(b=X) %>%
        select(b,everything())
    }) %>%
    do.call(rbind,.) %>%
    group_by(th,leaf,pred,metric) %>%
    summarize(
      Score=mean(value)
      ,LB=mean(value)-qnorm(0.975)*sd(value)/sqrt(n())
      ,UB=mean(value)+qnorm(0.975)*sd(value)/sqrt(n())
    ) %>%
    mutate(
      metric=
        metric %>%
        factor(c(
          'true_pred'
          ,'false_pred'
          ,'leaf_total'
          ,'leaf_ppv_npv'
          ,'pred_total'
          ,'leaf_p'
          ,'pred_ppv_npv'
          ,'total'
          ,'pred_p'
          ,'outcome_p'
        ))
    ) %>%
    arrange(th,leaf,pred,metric)
}