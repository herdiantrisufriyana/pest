# Create a function to compare calibration of the models
calib_comparison=function(calib){
  list(
    all=disease_model
    ,blood=disease_model2
    ,repro=disease_model3
  ) %>%
    lapply(X=names(.),Y=.,Z=calib,function(X,Y,Z){
      Y[[X]] %>%
        lapply(X=names(.),Y=.,Z=X,K=Z,function(X,Y,Z,K){
          if(X=='divnn'){
            Y[[X]][[K]]$pred=
              Y[[X]][[K]]$pred %>%
              mutate(Resample='Resample01')
          }
          Y[[X]][[K]]$pred %>%
            mutate(pred=round(event*10)/10) %>%
            group_by(Resample,pred) %>%
            summarize(
              event=sum(obs=='event')
              ,nonevent=sum(obs=='nonevent')
              ,total=n()
              ,obs=mean(obs=='event',na.rm=T)
            ) %>%
            ungroup() %>%
            group_by(pred) %>%
            summarize(
              avg=mean(obs,na.rm=T)
              ,lb=mean(obs,na.rm=T)-qnorm(0.975)*sd(obs,na.rm=T)/sqrt(n())
              ,ub=mean(obs,na.rm=T)+qnorm(0.975)*sd(obs,na.rm=T)/sqrt(n())
              ,event=sum(event)
              ,nonevent=sum(nonevent)
              ,total=sum(total)
            ) %>%
            ungroup() %>%
            cbind(
              lm(avg~pred,data=.) %>%
                broom::tidy() %>%
                mutate(term=term %>% str_remove_all('[:punct:]') %>% str_to_lower()) %>%
                select(term,estimate,std.error) %>%
                mutate(term=ifelse(term=='intercept','i','s')) %>%
                rename(avg=estimate,se=std.error) %>%
                pivot_wider(names_from=term,values_from=c(avg,se)) %>%
                mutate_all(round,2)
            ) %>%
            mutate(fit=avg_i+avg_s*pred) %>%
            mutate(
              event=event/max(event)
              ,nonevent=nonevent/max(nonevent)
              ,event_from=-0.2-0.25
              ,event_to=-0.2-0.25+0.25*event
              ,nonevent_from=-0.2-0.25-0.25*nonevent
              ,nonevent_to=-0.2-0.25
            ) %>%
            mutate(model=X)
        }) %>%
        do.call(rbind,.) %>%
        mutate(tissue=X)
    }) %>%
    do.call(rbind,.) %>%
    mutate(
      model=
        model %>%
        factor(c('pc_elnet','pc_rf','pc_gbm','divnn'))
      ,tissue=
        tissue %>%
        factor(c('all','blood','repro'))
    )
}