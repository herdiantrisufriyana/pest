# Create a function to compute AUROCs given a list
get_auroc=function(data,train_mod){
  data %>%
    pblapply(X=names(.),Y=.,Z=train_mod,function(X,Y,Z){
      Y[[X]] %>%
        lapply(X=names(.),Y=.,Z=Z[[X]],function(X,Y,Z){
          Y[[X]]=
            Y[[X]] %>%
            `colnames<-`(str_remove_all(colnames(.),'maternal_blood:'))
          Z[[X]] %>%
            predict(newdata=Y[[X]],type='prob') %>%
            cbind(select(Y[[X]],outcome)) %>%
            select(nonevent,event,outcome) %>%
            evalm(silent=T,showplots=F) %>%
            .$optres %>%
            .$Group1 %>%
            .['AUC-ROC',] %>%
            `rownames<-`(NULL) %>%
            mutate(grank=X)
        }) %>%
        do.call(rbind,.) %>%
        mutate(rank=X)
    }) %>%
    do.call(rbind,.) %>%
    unite(id,rank,grank,sep='.') %>%
    separate(CI,c('LB','UB'),sep='-') %>%
    mutate_at(c('LB','UB'),as.numeric) %>%
    left_join(
      train_mod %>%
        pblapply(X=names(.),Y=.,function(X,Y){
          Y[[X]] %>%
            lapply(X=names(.),Y=.,Z=X,function(X,Y,Z){
              K=Y[[X]] %>%
                .$finalModel %>%
                rpart.rules() %>%
                as.data.frame() %>%
                setNames(paste0('col',seq(ncol(.)))) %>%
                gather() %>%
                filter(
                  value!=''
                  & str_detect(value,'[:alpha:]')
                  & !value%in%c('when','is','to','null model')
                )
              if(nrow(K)==0){
                data.frame(
                  predictor=''
                  ,cand_predictor=
                    Y[[X]] %>%
                    .$trainingData %>%
                    colnames() %>%
                    .[.!='.outcome'] %>%
                    sort() %>%
                    paste0(collapse=',')
                  ,rank=Z
                  ,grank=X
                )
              }else{
                K %>%
                  pull(value) %>%
                  .[!duplicated(.)] %>%
                  setNames(NULL) %>%
                  data.frame(
                    predictor=.
                    ,cand_predictor=
                      Y[[X]] %>%
                      .$trainingData %>%
                      colnames() %>%
                      .[.!='.outcome'] %>%
                      sort() %>%
                      paste0(collapse=',')
                    ,rank=Z
                    ,grank=X
                  )
              }
            }) %>%
            do.call(rbind,.)
        }) %>%
        do.call(rbind,.) %>%
        group_by(rank,grank,cand_predictor) %>%
        summarize(predictors=paste0(sort(predictor),collapse=',')) %>%
        unite(id,rank,grank,sep='.')
      ,by='id'
    )
}