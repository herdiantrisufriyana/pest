# Create a function to re-evaluate the biomarker tree leaves
re_evalm=function(data
                  ,train_mod
                  ,metric_list
                  ,ref_list
                  ,nonevent_name='nonevent'
                  ,event_name='event'){
  eval=
    train_mod %>%
    predict(type='prob') %>%
    cbind(select(train_mod$trainingData,.outcome)) %>%
    rename(outcome=.outcome) %>%
    rename_at(event_name,function(x)'prob') %>%
    select(prob,outcome) %>%
    cbind(
      data.frame(
        th=paste0('th_',str_pad(seq(0,100),'3','left','0'))
        ,value=seq(0,1,0.01)
      ) %>%
        spread(th,value)
    ) %>%
    gather(key,th,-prob,-outcome) %>%
    select(-key) %>%
    mutate(
      pred=factor(ifelse(prob>=th,event_name,nonevent_name),levels(outcome))
    ) %>%
    group_by(th) %>%
    summarize(
      tp=sum(outcome==event_name & pred==event_name)
      ,fn=sum(outcome==event_name & pred==nonevent_name)
      ,fp=sum(outcome==nonevent_name & pred==event_name)
      ,tn=sum(outcome==nonevent_name & pred==nonevent_name)
    ) %>%
    mutate(
      sens=tp/(tp+fn)
      ,spec=tn/(fp+tn)
      ,ppv=tp/(tp+fp)
      ,npv=tn/(tn+fn)
      ,nb=(tp-fp*(th/(1-th)))/(tp+fn+fp+tn)
      ,mcc=(tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    ) %>%
    gather(metric,value,-th)
  
  th_value=
    metric_list %>%
    lapply(X=seq(length(.)),Y=.,Z=ref_list,function(X,Y,Z){
      if(length(Y[[X]])>1){
        K=eval %>%
          group_by(th) %>%
          mutate(
            sd=
              ifelse(metric%in%Y[[X]],value,NA) %>%
              sd(na.rm=T)
          ) %>%
          ungroup() %>%
          mutate(
            th_value=
              ifelse(sd==min(sd,na.rm=T),th,NA) %>%
              median(na.rm=T) %>%
              round(2)
          ) %>%
          select(-sd)
      }else{
        K=eval %>%
          group_by(th) %>%
          mutate(
            metric_closest=
              abs(ifelse(metric%in%Y[[X]],value,NA)-Z[[X]]) %>%
              min(na.rm=T)
          ) %>%
          ungroup() %>%
          mutate(
            th_value=
              ifelse(metric_closest==min(metric_closest,na.rm=T),th,NA) %>%
              median(na.rm=T) %>%
              round(2)
          ) %>%
          select(-metric_closest)
      }
      
      K %>%
        select(th_value) %>%
        rename_at('th_value',function(x)paste0(x,X))
    }) %>%
    do.call(cbind,.)
  
  th_list=
    eval %>%
    cbind(th_value) %>%
    gather(th_type,th_value,-th,-metric,-value) %>%
    mutate(
      th_metrics=
        th_type %>%
        sapply(function(x){
          y=str_remove_all(x,'th_value') %>%
            as.numeric()
          paste0(
            paste0(metric_list[[y]],collapse=',')
            ,'|'
            ,paste0(ref_list[[y]],collapse=',')
          )
        })
      ,th_metrics=
        ifelse(
          str_detect(th_metrics,',')
          ,str_remove_all(th_metrics,'\\|')
          ,th_metrics
        )
    ) %>%
    mutate(th_diff=abs(th-th_value)) %>%
    group_by(th_type) %>%
    filter(th_diff==min(th_diff)) %>%
    ungroup() %>%
    select(-th_diff) %>%
    spread(metric,value) %>%
    select(th,th_metrics) %>%
    lapply(X=seq(nrow(.)),Y=.,function(X,Y)Y[X,,drop=T]) %>%
    lapply(as.character) %>%
    lapply(setNames,c('th','th_metrics'))
  
  th_list %>%
    lapply(function(x){
      y=train_mod %>%
        predict(newdata=data,type='prob') %>%
        cbind(select(data,outcome)) %>%
        rename_at(event_name,function(x)'prob') %>%
        select(prob,outcome) %>%
        mutate(
          pred=
            factor(
              ifelse(prob>=x['th'],event_name,nonevent_name)
              ,levels(outcome)
            )
        )
      
      suppressWarnings(set.seed(33,sample.kind=sample.kind))
      seq(30) %>%
        lapply(X=.,function(X){
          y %>%
            .[sample(seq(nrow(.)),nrow(.),T),] %>%
            summarize(
              b=X
              ,tp=sum(outcome==event_name & pred==event_name)
              ,fn=sum(outcome==event_name & pred==nonevent_name)
              ,fp=sum(outcome==nonevent_name & pred==event_name)
              ,tn=sum(outcome==nonevent_name & pred==nonevent_name)
            ) %>%
            mutate(
              sens=tp/(tp+fn)
              ,spec=tn/(fp+tn)
              ,ppv=tp/(tp+fp)
              ,npv=tn/(tn+fn)
              ,nb=
                (tp-fp*(as.numeric(x['th'])/(1-as.numeric(x['th']))))
              /(tp+fn+fp+tn)
              ,mcc=(tp*tn-fp*fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
            )
        }) %>%
        do.call(rbind,.) %>%
        gather(metric,value,-b) %>%
        mutate(value=ifelse(is.nan(value),NA,value)) %>%
        group_by(metric) %>%
        summarize(
          Score=
            mean(value,na.rm=T)
          ,LB=
            mean(value,na.rm=T)
          -qnorm(0.975)*sd(value,na.rm=T)/sqrt(sum(!is.na(value)))
          ,UB=
            mean(value,na.rm=T)
          +qnorm(0.975)*sd(value,na.rm=T)/sqrt(sum(!is.na(value)))
        ) %>%
        mutate(
          metric=
            factor(
              str_to_upper(metric)
              ,str_to_upper(
                c('tp','fn','fp','tn','sens','spec','ppv','npv','nb','mcc')
              )
            )
        ) %>%
        arrange(metric) %>%
        mutate(
          th=x['th']
          ,th_metrics=x['th_metrics']
        ) %>%
        select(th,th_metrics,everything())
    }) %>%
    do.call(rbind,.)
}