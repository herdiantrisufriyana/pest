# Build a function to use only the common genes among all sets
take_common_genes=function(the_eset,clist){
  the_eset[Reduce(intersect,lapply(clist,rownames)),]
}