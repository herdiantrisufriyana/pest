# Create a function to show the shortest path
show_shortest_path=function(interaction_data,title=NULL){
  combination=
    interaction_data$inputs %>%
    expand.grid(.,.,stringsAsFactors=F) %>%
    filter(Var1!=Var2) %>%
    mutate(
      Var3=
        seq(nrow(.)) %>%
        sapply(X=.,Y=Var1,Z=Var2,function(X,Y,Z){
          c(Y[X],Z[X]) %>%
            sort() %>%
            paste(collapse=',')
        })
    ) %>%
    filter(!duplicated(Var3))
  
  shortest_path=
    combination %>%
    lapply(X=seq(nrow(.)),Y=.,function(X,Y){
      Z=suppressWarnings(
        interaction_data$nodes %>%
          rename_all(str_remove_all,'#') %>%
          select(-identifier) %>%
          rename(name=node,x=x_position,y=y_position) %>%
          right_join(
            interaction_data$nodes %>%
              rename_all(str_remove_all,'#') %>%
              select(node,x_position,y_position) %>%
              rename(node2=node,xend=x_position,yend=y_position) %>%
              right_join(
                interaction_data$edges %>%
                  select(-node1_string_id,-node2_string_id) %>%
                  `colnames<-`(str_remove_all(colnames(.),'#')) %>%
                  gather(relation,score,-node1,-node2,-combined_score) %>%
                  filter(score!=0) %>%
                  right_join(
                    graph_from_data_frame(.) %>%
                      induced_subgraph(
                        shortest_paths(.,Y$Var1[X],Y$Var2[X]) %>%
                          .$vpath %>%
                          unlist() %>%
                          names()
                      ) %>%
                      get.edgelist() %>%
                      as.data.frame() %>%
                      setNames(paste0('node',1:2))
                    ,by=paste0('node',1:2)
                  ) %>%
                  filter(!duplicated(.)) %>%
                  rename(name=node1)
                ,by='node2'
              )
            ,by='name'
          ) %>%
          select(-node2) %>%
          select(x,y,name,xend,yend,everything())
      )
      
      K=Z %>%
        ggplot(aes(x=x,y=y,xend=xend,yend=yend)) +
        geom_edges(
          aes(color=relation)
          ,size=1
          ,show.legend=F
        ) +
        geom_nodes(
          aes(color=color)
          ,size=5
          ,show.legend=F
        ) +
        geom_nodelabel_repel(
          aes(label=name,fill=name%in%interaction_data$inputs)
          ,label.padding=unit(0.15,'line')
          ,alpha=0.75
          ,size=3
          ,show.legend=F
        ) +
        facet_wrap(~relation) +
        scale_y_reverse() +
        scale_fill_manual(values=c('white','yellow')) +
        theme_blank() +
        ggtitle(title)
      
      print(K)
      Z
    })
  
  list(combination=combination,shortest_path=shortest_path)
}