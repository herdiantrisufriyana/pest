# Build a function to list suspected outliers by RLE
suspected_outliers=function(the_assay_data_unnorm){
  sweep(
    exprs(the_assay_data_unnorm),1,
    the_assay_data_unnorm %>%
      exprs() %>%
      rowMedians()
  ) %>%
    as.data.frame() %>%
    rownames_to_column(var='probe_id') %>%
    gather(sample_id,RLE,-probe_id) %>%
    mutate(sample_id=reorder(sample_id,RLE,mean)) %>%
    left_join(
      group_by(.,sample_id) %>%
        summarise(
          median=quantile(RLE,0.5),
          q1=quantile(RLE,0.25),
          q3=quantile(RLE,0.75)
        ) %>%
        ungroup() %>%
        mutate(
          minq=q1-1.5*(q3-q1),
          maxq=q3+1.5*(q3-q1)
        ),
      by='sample_id'
    ) %>%
    cbind(
      summarise(
        .,
        all_median=quantile(RLE,0.5),
        all_q1=quantile(RLE,0.25),
        all_q3=quantile(RLE,0.75)
      ) %>%
        ungroup() %>%
        mutate(
          all_minq=all_q1-1.5*(all_q3-all_q1),
          all_maxq=all_q3+1.5*(all_q3-all_q1)
        )
    ) %>%
    mutate(
      suspect_outlier=
        ifelse(
          median>all_maxq |
            median<all_minq |
            maxq>all_maxq |
            minq<all_minq,
          'yes','no'
        )
    )
}
