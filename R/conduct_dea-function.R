# Build a function to conduct differential expression analysis
conduct_dea=function(the_eset,nonevent='nonevent',paired=F,ignore_timepoint=F){
  
  model=list()
  
  model$eset=the_eset
  
  # Build the model and the null model
  if(paired){
    model$mod=model.matrix(
      ~outcome+subject_id
      ,data=
        model$eset %>%
        pData() %>%
        separate(title,c('tissue_type','subject_id'),sep='_')
    )
  }else if(ignore_timepoint){
    model$mod=model.matrix(~outcome+time_point,data=pData(model$eset))
  }else{
    model$mod=model.matrix(~outcome,data=pData(model$eset))
  }
  model$mod0=model.matrix(~1,data=pData(model$eset))
  
  # Computer surrogate variables
  model$svobj=sva(exprs(model$eset),model$mod,model$mod0,numSVmethod='leek')
  
  # Fit using limma with SV adjustment (removing heterogeneity)
  model$fit=lmFit(exprs(model$eset),cbind(model$mod,model$svobj$sv))

  # Compute Bayesian statistics
  model$eb=eBayes(model$fit)

  # Do multiple testing correction using BH method
  model$result=topTable(model$eb,coef=2,nrow(model$eset),adjust="BH",sort.by='none')

  model$result=
    model$result %>%
    rownames_to_column(var='feature') %>%
    left_join(
      model$eset %>%
        exprs() %>%
        t() %>%
        as.data.frame() %>%
        rownames_to_column(var='id') %>%
        left_join(
          model$eset %>%
            phenoData() %>%
            pData() %>%
            rownames_to_column(var='id') %>%
            select(id,outcome)
          ,by='id'
        ) %>%
        gather(feature,value,-id,-outcome) %>%
        group_by(feature,outcome) %>%
        summarize(avg=mean(value)) %>%
        ungroup() %>%
        mutate(
          outcome=ifelse(outcome==nonevent,'AveExpr1','AveExpr2')
        ) %>%
        spread(outcome,avg)
      ,by='feature'
    ) %>%
    arrange(factor(feature,rownames(model$result))) %>%
    column_to_rownames(var='feature')

  # Volcano plot
  model$vp=model$result  %>%
    rownames_to_column(var='hgnc_symbol') %>%
    mutate(
      deg=ifelse(adj.P.Val<0.05,paste0('adj. p-value<',0.05),'non-DEG'),
      deg=factor(deg,levels=c('non-DEG',paste0('adj. p-value<',0.05)))
    ) %>%
    mutate(
      minlog_p=-log10(P.Value),
      minlog_q=-log10(adj.P.Val)
    ) %>%
    ggplot(aes(x=logFC,y=minlog_p)) +
    geom_point(aes(color=deg)) +
    geom_text(
      aes(
        label=ifelse(abs(logFC)>2 & adj.P.Val<0.05,hgnc_symbol,NA),
        alpha=abs(logFC)*minlog_q,
        hjust=ifelse(logFC>0,0,1)
      ),
      vjust=0,size=2,na.rm=TRUE
    ) +
    geom_vline(xintercept=c(-2,0,2),linetype=3) +
    scale_x_continuous('logFC') +
    scale_y_continuous('-log10 p-value') +
    scale_color_manual('',values=c('#000000','#D55E00')) +
    scale_alpha(guide='none',range=c(-0.5,1)) +
    ggtitle(paste0(levels(the_eset$outcome),collapse=' vs. ')) +
    theme_minimal()

  print(model$vp)
  
  model
}
