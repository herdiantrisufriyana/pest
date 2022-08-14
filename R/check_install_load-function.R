#' Check, install, and load a package
#'
#' Install/change a package with a specific version, and load it.
#'
#' @param pkg Package, a character of the name.
#' @param version Version, a character of the version.
#' @param repo Repository, a character of the repository which are cran, bioc, 
#' devtools, and github.
#' @param load Load, a logical whether the package is loaded.
#'
#' @details
#' None.
#' 
#' @return None.
#' 

check_install_load=function(pkg
                            ,version
                            ,repo=c('cran','bioc','devtools','github')
                            ,load=T){
  
  # Check if this package is already installed; otherwise, set install true
  if(pkg%in%installed.packages()[,1]){
    
    # Set install true if only the version is not the expected one
    if(packageVersion(pkg)==version){
      is.install=FALSE
    }else{
      is.install=TRUE
    }
    
  }else{
    is.install=TRUE
  }
  
  # Install the oackage depending on the repo
  if(is.install){
    if(repo=='cran'){
      install.packages(pkg)
    }else if(repo=='bioc'){
      BiocManager::install(pkg,update=F)
    }else if(repo=='devtools'){
      devtools::install_version(pkg,version,upgrade=F)
    }else{
      devtools::install_github(paste0(repo,'/',pkg),upgrade=F)
    }
  }
  
  # If warranted, load the library
  if(load)library(pkg,character.only=T)
}