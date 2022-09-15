# Create a function for Q-Q normalization targeting controls
dea_qn=function(exprs_mat,the_dea_result){
  exprs_mat %>%
    normalize.quantiles.use.target(
      target=
        the_dea_result[rownames(.),'AveExpr1',drop=F] %>%
        pull(AveExpr1)
    ) %>%
    `dimnames<-`(dimnames(exprs_mat))
}