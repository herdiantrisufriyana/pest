# Create a function to standardize features ignoring outlier
std_no_outlier=function(input,scaler){
  scaler=
    scaler %>%
    select(-outcome) %>%
    lapply(function(x){
      median=quantile(x,0.5)
      q1=quantile(x,0.25)
      q3=quantile(x,0.75)
      minq=q1-1.5*(q3-q1)
      maxq=q3+1.5*(q3-q1)
      outliers=x<minq | x>maxq
      avg=
        x %>%
        .[!outliers] %>%
        mean()
      
      std=
        x %>%
        .[!outliers] %>%
        sd()
      
      list(avg=avg,std=std)
    })
  
  input %>%
    lapply(X=seq(ncol(.)),Y=.,Z=scaler,function(X,Y,Z){
      if(X==1){
        Y[[X]]
      }else{
        (Y[[X]]-Z[[colnames(Y)[X]]]$avg)/Z[[colnames(Y)[X]]]$std
      }
    }) %>%
    as.data.frame() %>%
    `dimnames<-`(dimnames(input))
}