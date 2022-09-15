# Build a function to download annotation data
download_annotation=function(the_platform,the_mart){
  getBM(mart=the_mart,attributes=c(
    the_platform,
    'ensembl_gene_id','entrezgene_id','hgnc_symbol'
  )) %>%
    apply(2,as.character) %>%
    apply(2,trimws) %>%
    apply(2,function(x) gsub('^$|^ $',NA,x)) %>%
    as.data.frame(stringsAsFactors=FALSE) %>%
    setNames(c('probe_id','ensembl_gene_id','entrezgene_id','hgnc_symbol')) %>%
    filter(!(is.na(probe_id))) %>%
    group_by(probe_id) %>%
    summarise_all(function(x){
      ifelse(
        length(unique(x))>1,
        paste(unique(x),collapse='///'),
        x
      )
    }) %>%
    ungroup() %>%
    filter(
      !(is.na(ensembl_gene_id)) &
        !(is.na(entrezgene_id)) &
        !(is.na(hgnc_symbol))
    ) %>%
    filter(
      !(str_detect(ensembl_gene_id,'///')) &
        !(str_detect(entrezgene_id,'///')) &
        !(str_detect(hgnc_symbol,'///'))
    ) %>%
    column_to_rownames(var='probe_id')
}