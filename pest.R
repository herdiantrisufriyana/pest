################################################################################
################################################################################
# Set up programming environment
################################################################################
################################################################################

#####{Stop codes if R version is false}
r_version='4.0.2'
if(paste0(version$major,'.',version$minor)!=r_version){
  stop(paste0('This code uses R version of ',r_version))
}
################################################################################

#####{Ensure the zip files exist}
if(file.exists('pest_files.zip')){
  if(!file.exists('renv')) unzip('pest_files.zip')
}else{
  if(!file.exists('renv')){
    stop(
      paste0(
        'Please manually dowload the zip files (17.47 GB) from this link '
        ,'(https://drive.google.com/file/d/12sIGr0ys07WyMNDdZ55hmhTzF7PRoBA4/view?usp=sharing) '
        ,'Put the zip file the same folder with this pest.Rmd and pest.R'
        ,r_version
      )
    )
  }
}
################################################################################

#####{Set up reproducible environment, include=FALSE}
renv::activate()
renv::restore(prompt=F)
reticulate::install_miniconda(update=F,force=T)
################################################################################

#####{Set sample kind, include=FALSE}
# sample.kind=NULL # if using R 3.5 or earlier
sample.kind='Rounding' # if using R 3.6 or later
################################################################################

#####{Set to run or not run very heavy computations, include=FALSE}
# Many computations were very heavy;
# thus, we provided the RDS files as substitutes and load only ones
# that can be ran in most computers.
# Set to TRUE if you want to run the very heavy computations.
run_heavy_computation=F
load_outlier_comp=F
load_trained_divnn=F
load_emulators=F
show_best_tree=F
################################################################################

#####{Recommended memory limit, include=FALSE}
# Minimum 1 GB free memory should be allocated
memory.limit(size=1000000000)

# Minimum connection size
Sys.setenv('VROOM_CONNECTION_SIZE'=128000*1024)
################################################################################

#####{Load packages, include=FALSE}
library(tidyverse)
library(dslabs)
library(GEOquery)
library(Biobase)
library(biomaRt)
library(oligo)
library(limma)
library(illuminaio)
library(parallelDist)
library(WGCNA)
library(sva)
library(pbapply)
library(lubridate)
library(rsdr)
library(caret)
library(glmnet)
library(Rborist)
library(gbm)
library(gam)
library(MLeval)
library(doParallel)
library(preprocessCore)
library(matrixStats)
library(Rtsne)
library(clixo)
library(divnn)
library(igraph)
library(reticulate)
library(ggnetwork)
library(rpart.plot)
library(knitr)
library(kableExtra)

options(dplyr.summarise.inform=F)
dslabs::ds_theme_set()
select=dplyr::select
summarize=dplyr::summarize
rename=dplyr::rename
slice=dplyr::slice
if(load_trained_divnn){
  reticulate::use_condaenv('./renv/python/condaenvs/renv-python',required=T)
}
renv::use_python(name='./renv/python/condaenvs/renv-python')
################################################################################






################################################################################
################################################################################
################################################################################
# Data load
################################################################################
################################################################################
################################################################################

#####{Create empty lists for DEA-related variables, include=FALSE}
pdata=list()
edata=list()
fdata=list()
rdata=list()
cdata=list()
odata=list()
ldata=list()
gdata=list()
################################################################################

################################################################################
## Phenotype data
################################################################################

#####{Load GSE108497 phenotype, include=FALSE}
# Expression set with the processed transcript expression
pdata$GSE108497$eset=
  suppressWarnings(suppressMessages(
    getGEO(
      filename='data/GSE108497/GSE108497_series_matrix.txt.gz'
      ,getGPL=F
      ,parseCharacteristics=F
    )
  ))

# Metadata of phenotype variables
pdata$GSE108497$phenolabel=
  suppressWarnings(suppressMessages(
    read_csv('data/GSE108497/GSE108497_varMetadata.csv') %>%
      mutate(code=str_to_lower(code))
  ))

# Phenotype data
pdata$GSE108497$phenotype=
  pdata$GSE108497$eset %>%
  phenoData() %>%
  pData() %>%
  select(c('title',colnames(.) %>% .[str_detect(.,':ch1')])) %>%
  setNames(colnames(.) %>% str_remove_all(':ch1'))
################################################################################

#####{Load GSE73685 phenotype, include=FALSE}
# Expression set with the processed transcript expression
pdata$GSE73685$eset=
  suppressWarnings(suppressMessages(
    getGEO(
      filename='data/GSE73685/GSE73685_series_matrix.txt.gz'
      ,getGPL=F
      ,parseCharacteristics=F
    )
  ))
################################################################################

#####{Load GSE85307 phenotype, include=FALSE}
# Expression set with the processed transcript expression
pdata$GSE85307$eset=
  suppressWarnings(suppressMessages(
    getGEO(
      filename='data/GSE85307/GSE85307-GPL6244_series_matrix.txt.gz'
      ,getGPL=F
      ,parseCharacteristics=F
    )
  ))
################################################################################

#####{Load GSE86200 phenotype, include=FALSE}
# Expression set with the processed transcript expression
pdata$GSE86200$eset=
  suppressWarnings(suppressMessages(
    getGEO(
      filename='data/GSE86200/GSE86200_series_matrix.txt.gz'
      ,getGPL=F
      ,parseCharacteristics=F
    )
  ))
################################################################################

#####{Load GSE149437 phenotype, include=FALSE}
# Expression set with the processed transcript expression
pdata$GSE149437$eset=
  suppressWarnings(suppressMessages(
    getGEO(
      filename='data/GSE149437/GSE149437_series_matrix.txt.gz'
      ,getGPL=F
      ,parseCharacteristics=F
    )
  ))
################################################################################

#####{Load GSE177477 phenotype, include=FALSE}
# Expression set with the processed transcript expression
pdata$GSE177477$eset=
  suppressWarnings(suppressMessages(
    getGEO(
      filename='data/GSE177477/GSE177477_series_matrix.txt.gz'
      ,getGPL=F
      ,parseCharacteristics=F
    )
  ))

# Add phenotype data from the publication
pdata$GSE177477$eset=
  pdata$GSE177477$eset %>%
  `pData<-`(
    phenoData(.) %>%
      pData() %>%
      mutate(`id:ch1`=sapply(title,function(x)str_split_fixed(x,'\\s+',2)[1])) %>%
      rownames_to_column(var='gsm') %>%
      left_join(
        suppressWarnings(suppressMessages(
          read_csv('data/GSE177477/GSE177477_phenotype.csv')
        ))
        ,by='id:ch1'
      ) %>%
      column_to_rownames(var='gsm')
  ) %>%
  `varMetadata<-`(
    varMetadata(.) %>%
      rownames_to_column(var='variable') %>%
      select(-labelDescription) %>%
      left_join(
        suppressWarnings(suppressMessages(
          read_csv('data/GSE177477/GSE177477_varMetadata.csv')
        ))
        ,by='variable'
      ) %>%
      column_to_rownames(var='variable')
  )
################################################################################

#####{Load GSE192902 phenotype, include=FALSE}
# Expression set with the processed transcript expression
pdata$GSE192902$eset=
  suppressWarnings(suppressMessages(
    getGEO(
      filename='data/GSE192902/GSE192902_series_matrix.txt.gz'
      ,getGPL=F
      ,parseCharacteristics=F
    )
  ))
################################################################################

#####{Summarize periods of all the datasets, include=FALSE}
data_period=
  rbind(
    
    pdata$GSE108497$eset %>%
      phenoData() %>%
      pData() %>%
      select(submission_date) %>%
      mutate(submission_date=mdy(submission_date)) %>%
      gather() %>%
      mutate(dataset='GSE108497')
    
    ,pdata$GSE73685$eset %>%
      phenoData() %>%
      pData() %>%
      select(submission_date) %>%
      mutate(submission_date=mdy(submission_date)) %>%
      gather() %>%
      mutate(dataset='GSE73685')
    
    ,pdata$GSE85307$eset %>%
      phenoData() %>%
      pData() %>%
      select(submission_date,`microarray_hybridization_date:ch1`) %>%
      mutate(
        submission_date=
          mdy(submission_date)
        ,`microarray_hybridization_date:ch1`=
          mdy(`microarray_hybridization_date:ch1`)
      ) %>%
      gather() %>%
      mutate(dataset='GSE85307')
    
    ,pdata$GSE86200$eset %>%
      phenoData() %>%
      pData() %>%
      select(submission_date,`microarray hybridization date:ch1`) %>%
      mutate(
        submission_date=
          mdy(submission_date)
        ,`microarray hybridization date:ch1`=
          mdy(`microarray hybridization date:ch1`)
      ) %>%
      gather() %>%
      mutate(dataset='GSE86200')
    
    ,pdata$GSE149437$eset %>%
      phenoData() %>%
      pData() %>%
      select(submission_date,`sample_date:ch1`) %>%
      mutate(
        submission_date=
          mdy(submission_date)
        ,`sample_date:ch1`=
          dmy(`sample_date:ch1`)
      ) %>%
      gather() %>%
      mutate(dataset='GSE149437')
    
    ,pdata$GSE177477$eset %>%
      phenoData() %>%
      pData() %>%
      select(submission_date) %>%
      mutate(submission_date=mdy(submission_date)) %>%
      gather() %>%
      mutate(dataset='GSE177477') %>%
      rbind(
        data.frame(
          key='sample_date'
          ,value=ymd(c('2020/04/01','2020/10/31'))
          ,dataset='GSE177477'
        )
      )
    
    ,pdata$GSE192902$eset %>%
      phenoData() %>%
      pData() %>%
      select(submission_date) %>%
      mutate(submission_date=mdy(submission_date)) %>%
      gather() %>%
      mutate(dataset='GSE192902') %>%
      rbind(
        data.frame(
          key='guidelines_date'
          ,value=ymd(c('2020/06/01','2020/06/01'))
          ,dataset='GSE192902'
        )
      )
  ) %>%
  group_by(dataset,key) %>%
  summarize(
    from=min(value)
    ,to=max(value)
  ) %>%
  mutate(
    dataset=
      dataset %>%
      factor(names(pdata) %>% .[c(2,1,3:length(.))] %>% rev())
    ,key=
      key %>%
      str_remove_all('\\:ch1') %>%
      str_replace_all('_',' ') %>%
      str_to_sentence() %>%
      str_remove_all(' date') %>%
      paste0(
        case_when(
          .=='Submission'~'*'
          ,.=='Guidelines'~'†'
          ,TRUE~''
        )
      ) %>%
      factor(.[!duplicated(.)] %>% rev())
    ,dataset_group=
      case_when(
        dataset%in%c('GSE73685')~'1. Derivation dataset'
        ,dataset%in%c('GSE108497')~'2. Development dataset'
        ,dataset%in%c('GSE85307','GSE86200','GSE149437')
        ~'3. Replication datasets'
        ,dataset%in%c('GSE177477')~'4. COVID-19 dataset'
        ,dataset%in%c('GSE192902')
        ~'5. Preeclampsia gene set during the pandemic'
        ,TRUE~''
      )
  )
################################################################################






################################################################################
## Expression data
################################################################################

#####{Load GSE108497 expression, include=FALSE}
# Load the raw data
if(run_heavy_computation){
  edata$GSE1084970$raw=
    suppressWarnings(suppressMessages(
      read_tsv('data/GSE108497/GSE108497_non-normalized_data.txt')
    )) %>%
    lapply(X=1,Y=.,function(X,Y){
      
      Z=pdata$GSE108497$eset %>%
        phenoData() %>%
        pData() %>%
        select(description) %>%
        rename(key=description) %>%
        rownames_to_column(var='id') %>%
        mutate(id=paste0('sample_',id))
      
      K=Y %>%
        colnames() %>%
        data.frame(key=.) %>%
        left_join(Z,by='key') %>%
        arrange(factor(key,colnames(Y))) %>%
        mutate(key=ifelse(is.na(id),key,id)) %>%
        pull(key)
      
      gather(Y,key,value,-ID_REF) %>%
        left_join(Z,by='key') %>%
        mutate(key=ifelse(is.na(id),key,id)) %>%
        select(-id) %>%
        spread(key,value) %>%
        .[,K]
      
    }) %>%
    .[[1]] %>%
    lapply(X=1,Y=.,function(X,Y){
      Z='data/GSE108497/GSE108497_non-normalized_data_modified.txt'
      write_tsv(Y,Z)
      Z
    }) %>%
    .[[1]] %>%
    read.ilmn(probeid="ID_REF",expr="sample_",verbose=F)
}

# Background correction without normalization
if(run_heavy_computation){
  edata$GSE108497$unnorm=
    edata$GSE108497$raw %>%
    nec() %>%
    .$E %>%
    ExpressionSet()
  saveRDS(edata$GSE108497$unnorm,'data/GSE108497/GSE108497_unnorm.rds')
}else if(load_outlier_comp){
  edata$GSE108497$unnorm=readRDS('data/GSE108497/GSE108497_unnorm.rds')
}

# Background correction with normalization
if(run_heavy_computation){
  edata$GSE108497$norm=
    edata$GSE108497$raw %>%
    neqc() %>%
    .$E %>%
    ExpressionSet()
  saveRDS(edata$GSE108497$norm,'data/GSE108497/GSE108497_norm.rds')
}else{
  edata$GSE108497$norm=readRDS('data/GSE108497/GSE108497_norm.rds')
}
################################################################################

#####{Load GSE73685 expression, include=FALSE}
# Load the raw data
if(run_heavy_computation){
  edata$GSE73685$raw=
    list.files('data/GSE73685/raw',full.names=T) %>%
    read.celfiles(sampleNames=colnames(pdata$GSE73685$eset),verbose=F)
}

# Background correction without normalization
if(run_heavy_computation){
  edata$GSE73685$unnorm=
    edata$GSE73685$raw %>%
    rma(background=T,normalize=F)
  saveRDS(edata$GSE73685$unnorm,'data/GSE73685/GSE73685_unnorm.rds')
}else if(load_outlier_comp){
  edata$GSE73685$unnorm=readRDS('data/GSE73685/GSE73685_unnorm.rds')
}

# Background correction with normalization
if(run_heavy_computation){
  edata$GSE73685$norm=
    edata$GSE73685$raw %>%
    rma(background=T,normalize=T)
  saveRDS(edata$GSE73685$norm,'data/GSE73685/GSE73685_norm.rds')
}else{
  edata$GSE73685$norm=readRDS('data/GSE73685/GSE73685_norm.rds')
}
################################################################################

#####{Load GSE85307 expression, include=FALSE}
# Load the raw data
if(run_heavy_computation){
  edata$GSE85307$raw=
    list.files('data/GSE85307/raw',full.names=T) %>%
    read.celfiles(sampleNames=colnames(pdata$GSE85307$eset),verbose=F)
}

# Background correction without normalization
if(run_heavy_computation){
  edata$GSE85307$unnorm=
    edata$GSE85307$raw %>%
    rma(background=T,normalize=F)
  saveRDS(edata$GSE85307$unnorm,'data/GSE85307/GSE85307_unnorm.rds')
}else if(load_outlier_comp){
  edata$GSE85307$unnorm=readRDS('data/GSE85307/GSE85307_unnorm.rds')
}

# Background correction with normalization
if(run_heavy_computation){
  edata$GSE85307$norm=
    edata$GSE85307$raw %>%
    rma(background=T,normalize=T)
  saveRDS(edata$GSE85307$norm,'data/GSE85307/GSE85307_norm.rds')
}else{
  edata$GSE85307$norm=readRDS('data/GSE85307/GSE85307_norm.rds')
}
################################################################################

#####{Load GSE86200 expression, include=FALSE}
# Load the raw data
if(run_heavy_computation){
  edata$GSE86200$raw=
    list.files('data/GSE86200/raw',full.names=T) %>%
    read.idat(bgxfile='data/GSE86200/HumanHT-12_V4_0_R2_15002873_B.bgx') %>%
    lapply(X=1,Y=.,function(X,Y){
      recol=function(x){
        x %>%
          `colnames<-`(
            colnames(.) %>%
              str_remove_all('data/GSE86200/raw/') %>%
              str_split_fixed('_',2) %>%
              .[,1]
          ) 
        
      }
      Y$E=recol(Y$E)
      Y$other$NumBeads=recol(Y$other$NumBeads)
      Y$other$STDEV=recol(Y$other$STDEV)
      Y
    }) %>%
    .[[1]]
}

# Background correction without normalization
if(run_heavy_computation){
  edata$GSE86200$unnorm=
    edata$GSE86200$raw %>%
    nec() %>%
    .$E %>%
    .[rownames(.) %>% .[!duplicated(.)],] %>%
    as.data.frame() %>%
    rownames_to_column(var='Array_Address_Id') %>%
    mutate(seq=seq(nrow(.))) %>%
    left_join(
      edata$GSE86200$raw$genes %>%
        select(Array_Address_Id,Probe_Id) %>%
        filter(!duplicated(.)) %>%
        mutate(Array_Address_Id=as.character(Array_Address_Id))
      ,by='Array_Address_Id'
    ) %>%
    arrange(seq) %>%
    select(-seq,-Array_Address_Id) %>%
    select(Probe_Id,everything()) %>%
    column_to_rownames(var='Probe_Id') %>%
    as.matrix() %>%
    ExpressionSet()
  saveRDS(edata$GSE86200$unnorm,'data/GSE86200/GSE86200_unnorm.rds')
}else if(load_outlier_comp){
  edata$GSE86200$unnorm=readRDS('data/GSE86200/GSE86200_unnorm.rds')
}

# Background correction with normalization
if(run_heavy_computation){
  edata$GSE86200$norm=
    edata$GSE86200$raw %>%
    neqc() %>%
    .$E %>%
    .[rownames(.) %>% .[!duplicated(.)],] %>%
    as.data.frame() %>%
    rownames_to_column(var='Array_Address_Id') %>%
    mutate(seq=seq(nrow(.))) %>%
    left_join(
      edata$GSE86200$raw$genes %>%
        select(Array_Address_Id,Probe_Id) %>%
        filter(!duplicated(.)) %>%
        mutate(Array_Address_Id=as.character(Array_Address_Id))
      ,by='Array_Address_Id'
    ) %>%
    arrange(seq) %>%
    select(-seq,-Array_Address_Id) %>%
    select(Probe_Id,everything()) %>%
    column_to_rownames(var='Probe_Id') %>%
    as.matrix() %>%
    ExpressionSet()
  saveRDS(edata$GSE86200$norm,'data/GSE86200/GSE86200_norm.rds')
}else{
  edata$GSE86200$norm=readRDS('data/GSE86200/GSE86200_norm.rds')
}
################################################################################

#####{Create calling list to load GSE149437 expression, include=FALSE}
# Load the raw data
suppressWarnings(set.seed(33,sample.kind=sample.kind))
GSE149437_calling_list=
  pdata$GSE149437$eset %>%
  phenoData() %>%
  pData() %>%
  select(c('title',colnames(.) %>% .[str_detect(.,':ch1')])) %>%
  setNames(
    colnames(.) %>%
      str_remove_all(':ch1') %>%
      str_replace_all('\\s','_')
  ) %>%
  mutate(
    blood_draw_at_gestation_week=as.numeric(gestational_age)
    ,time_point=case_when(
      between(ceiling(blood_draw_at_gestation_week),0,15)~'<16 weeks'
      ,between(ceiling(blood_draw_at_gestation_week),16,23)~'16-23 weeks'
      ,between(ceiling(blood_draw_at_gestation_week),24,31)~'24-31 weeks'
      ,between(ceiling(blood_draw_at_gestation_week),32,40)~'32-40 weeks'
      ,TRUE~''
    )
  ) %>%
  rename(outcome=group) %>%
  select(title,outcome,time_point) %>%
  rownames_to_column(var='id') %>%
  group_by(outcome,time_point) %>%
  mutate(subset=sample(rep(1:10,ceiling(n()/10))[1:n()],n(),F)) %>%
  ungroup() %>%
  unite(filename,id,title,sep='_') %>%
  mutate(filename=paste0('data/GSE149437/raw/',filename,'.CEL.gz')) %>%
  lapply(X=1:10,Y=.,function(X,Y){
    Y %>%
      filter(subset==X) %>%
      pull(filename)
  })
################################################################################

#####{Load GSE149437 expression, include=FALSE}
# Load the raw data
if(run_heavy_computation){
  edata$GSE149437$raw=
    list.files('data/GSE149437/raw',full.names=T) %>%
    read.celfiles(sampleNames=colnames(pdata$GSE149437$eset),verbose=F)
  
  gc()
}

if(run_heavy_computation){
  # Chunk the raw data
  edata$GSE149437$raw=
    GSE149437_calling_list %>%
    pblapply(X=seq(length(.)),Y=.,Z=edata$GSE149437$raw,function(X,Y,Z){
      K=Z %>%
        .[,Y[[X]] %>%
            str_remove_all('data\\/GSE149437\\/raw/|\\.CEL\\.gz') %>%
            str_split_fixed('_',2) %>%
            .[,1]]
      
      K %>%
        rma(background=T,normalize=F) %>%
        saveRDS(paste0(
          'data/GSE149437/GSE149437_unnorm_'
          ,str_pad(X,2,'left','0')
          ,'.rds'
        ))
      
      K %>%
        rma(background=T,normalize=T) %>%
        saveRDS(paste0(
          'data/GSE149437/GSE149437_norm_'
          ,str_pad(X,2,'left','0')
          ,'.rds'
        ))
      
      if(X==1){
        Z
      }else{
        NULL
      }
    }) %>%
    .[[1]]
}else{
  # Background correction without normalization (per chunk)
  edata$GSE149437$unnorm=
    seq(length(GSE149437_calling_list)) %>%
    pblapply(function(x){
      readRDS(paste0(
        'data/GSE149437/GSE149437_unnorm_'
        ,str_pad(x,2,'left','0')
        ,'.rds'
      ))
    })
  
  # Background correction with normalization (per chunk)
  edata$GSE149437$norm=
    seq(length(GSE149437_calling_list)) %>%
    pblapply(function(x){
      readRDS(paste0(
        'data/GSE149437/GSE149437_norm_'
        ,str_pad(x,2,'left','0')
        ,'.rds'
      ))
    })
  
  # Background correcting
  # Calculating Expression
  # Background correcting
  # Normalizing
  # Calculating Expression
  #   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=57m 33s
  #   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=03s  
  #   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=03s
  
  # Background correcting
  # Calculating Expression
  # Background correcting
  # Normalizing
  # Calculating Expression
  #   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=01h 05m 53s
  #   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=03s  
  #   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=03s
  
  # Background correction without normalization (combining chunks)
  edata$GSE149437$unnorm=
    edata$GSE149437$unnorm %>%
    lapply(X=1,Y=.,function(X,Y){
      for(i in seq(length(Y))){
        if(i==1){
          Z=Y[[i]]
        }else{
          Z=a4Base::combineTwoExpressionSet(Z,Y[[i]])
        }
      }
      Z
    }) %>%
    .[[1]] %>%
    .[,colnames(pdata$GSE149437$eset)]
  
  # Background correction with normalization (combining chunks)
  edata$GSE149437$norm=
    edata$GSE149437$norm %>%
    lapply(X=1,Y=.,function(X,Y){
      for(i in seq(length(Y))){
        if(i==1){
          Z=Y[[i]]
        }else{
          Z=a4Base::combineTwoExpressionSet(Z,Y[[i]])
        }
      }
      Z
    }) %>%
    .[[1]] %>%
    .[,colnames(pdata$GSE149437$eset)]
}
################################################################################

#####{Load GSE177477 expression, include=FALSE}
# Load the raw data
if(run_heavy_computation){
  edata$GSE177477$raw=
    list.files('data/GSE177477/raw',full.names=T) %>%
    .[!str_detect(.,'chp')] %>%
    read.celfiles(sampleNames=colnames(pdata$GSE177477$eset),verbose=F)
}

# Background correction without normalization
if(run_heavy_computation){
  edata$GSE177477$unnorm=
    edata$GSE177477$raw %>%
    rma(background=T,normalize=F)
  saveRDS(edata$GSE177477$unnorm,'data/GSE177477/GSE177477_unnorm.rds')
}else if(load_outlier_comp){
  edata$GSE177477$unnorm=readRDS('data/GSE177477/GSE177477_unnorm.rds')
}

# Background correction with normalization
if(run_heavy_computation){
  edata$GSE177477$norm=
    edata$GSE177477$raw %>%
    rma(background=T,normalize=T)
  saveRDS(edata$GSE177477$norm,'data/GSE177477/GSE177477_norm.rds')
}else{
  edata$GSE177477$norm=readRDS('data/GSE177477/GSE177477_norm.rds')
}
################################################################################






################################################################################
## Feature data
################################################################################

#####{Download annotation database, eval=FALSE, include=FALSE}
# mart=useMart("ensembl")
# saveRDS(mart,file='data/mart.rds')
# ensembl=useDataset("hsapiens_gene_ensembl",mart)
# saveRDS(ensembl,file='data/ensembl.rds')
################################################################################

#####{Load the saved database, include=FALSE}
mart=readRDS(file='data/mart.rds') # Downloaded on 2021-07-27
ensembl=readRDS(file='data/ensembl.rds') # Downloaded on 2021-07-27
listAttributes(ensembl)
################################################################################

#####{Download annotation data, include=FALSE}
source('R/download_annotation-function.R')

# GPL10558, downloaded on 2021-07-27
# GPL10558=
#   download_annotation(
#     the_platform='illumina_humanht_12_v4',
#     the_mart=ensembl
#   )
# saveRDS(GPL10558,'data/GPL10558.rds')
GPL10558=readRDS('data/GPL10558.rds')

# GPL6244, downloaded on 2021-07-27
# GPL6244=
#   download_annotation(
#     the_platform='affy_hugene_1_0_st_v1',
#     the_mart=ensembl
#   )
# saveRDS(GPL6244,'data/GPL6244.rds')
GPL6244=readRDS('data/GPL6244.rds')

# GPL28460, download on 2021-09-12
# GPL28460=
#   download_annotation(
#     the_platform='affy_hta_2_0',
#     the_mart=ensembl
#   )
# saveRDS(GPL28460,'data/GPL28460.rds')
GPL28460=
  readRDS('data/GPL28460.rds') %>%
  rownames_to_column(var='id') %>%
  mutate(id=paste0(id,'.1')) %>%
  column_to_rownames(var='id')

# GPL23159, download on 2022-02-11
# GPL23159=list()
# GPL23159$refseq_mrna=
#   download_annotation(
#     the_platform='refseq_mrna',
#     the_mart=ensembl
#   )
# GPL23159$ensembl_transcript_id=
#   download_annotation(
#     the_platform='ensembl_transcript_id',
#     the_mart=ensembl
#   )
# GPL23159$soft_file=
#   suppressWarnings(suppressMessages(
#     read_tsv('data/GSE177477/GSE177477_family.soft.gz',skip=137) %>%
#       setNames(
#         colnames(.) %>%
#           str_remove_all('\\s+|[:punct:]') %>%
#           str_to_lower()
#       )
#   )) %>%
#   filter(
#     probesetid
#     %in% intersect(
#       .$probesetid
#       ,rownames(edata$GSE177477$norm)
#     )
#   ) %>%
#   right_join(
#     lapply(X=1,Y=.,function(X,Y){
#         Y$spotid10 %>%
#           pblapply(function(x)str_split(x,' /// ')[[1]]) %>%
#           pblapply(function(x)x[str_count(x,' // ')==8]) %>%
#           pblapply(sapply,function(x)str_split(x,' // ')[[1]]) %>%
#           pblapply(`dimnames<-`,NULL) %>%
#           pblapply(t) %>%
#           pblapply(as.data.frame) %>%
#           pblapply(X=seq(length(.)),Y=Y$probesetid,Z=.,function(X,Y,Z){
#             Z[[X]] %>%
#               mutate(probesetid=Y[X]) %>%
#               select(probesetid,everything())
#           }) %>%
#           do.call(rbind,.)
#       }) %>%
#       .[[1]]
#     ,by='probesetid'
#   )
# GPL23159=
#   GPL23159$soft_file %>%
#   left_join(
#     GPL23159$refseq_mrna %>%
#       rownames_to_column(var='V1') %>%
#       mutate(V2='RefSeq') %>%
#       rbind(
#         GPL23159$ensembl_transcript_id %>%
#           rownames_to_column(var='V1') %>%
#           mutate(V2='ENSEMBL')
#       )
#     ,by=c('V1','V2')
#   ) %>%
#   select(probesetid,ensembl_gene_id,entrezgene_id,hgnc_symbol) %>%
#   filter(!duplicated(.)) %>%
#   gather(annotation,value,-probesetid) %>%
#   filter(!is.na(value)) %>%
#   group_by(probesetid,annotation) %>%
#   summarise(value=paste0(value,collapse='///')) %>%
#   spread(annotation,value) %>%
#   right_join(
#     edata$GSE177477$norm %>%
#       rownames() %>%
#       data.frame(probesetid=.)
#     ,by='probesetid'
#   ) %>%
#   filter(
#     !(is.na(ensembl_gene_id)) &
#       !(is.na(entrezgene_id)) &
#       !(is.na(hgnc_symbol))
#   ) %>%
#   filter(
#     !(str_detect(ensembl_gene_id,'///')) &
#       !(str_detect(entrezgene_id,'///')) &
#       !(str_detect(hgnc_symbol,'///'))
#   ) %>%
#   column_to_rownames(var='probesetid')
# saveRDS(GPL23159,'data/GPL23159.rds')
GPL23159=readRDS('data/GPL23159.rds')
################################################################################

#####{Take common transcripts of the dataset and the annotation, include=FALSE}
fdata$GSE108497=
  GPL10558 %>%
  .[intersect(rownames(edata$GSE108497$norm),rownames(.)),]

fdata$GSE73685=
  GPL6244 %>%
  .[intersect(rownames(edata$GSE73685$norm),rownames(.)),]

fdata$GSE85307=
  GPL6244 %>%
  .[intersect(rownames(edata$GSE85307$norm),rownames(.)),]

fdata$GSE86200=
  GPL10558 %>%
  .[intersect(rownames(edata$GSE86200$norm),rownames(.)),]

fdata$GSE149437=
  GPL28460 %>%
  .[intersect(rownames(edata$GSE149437$norm),rownames(.)),]

fdata$GSE177477=
  GPL23159 %>%
  .[intersect(rownames(edata$GSE177477$norm),rownames(.)),]
################################################################################






################################################################################
################################################################################
################################################################################
# Data preprocessing
################################################################################
################################################################################
################################################################################

################################################################################
## Removal of technical outliers
################################################################################

#####{Identify outliers of unnormalized expression by RLE, include=FALSE}
source('R/suspected_outliers-function.R')

if(load_outlier_comp){
  rdata$GSE108497=
    edata$GSE108497$unnorm %>%
    .[intersect(rownames(GPL10558),rownames(.)),] %>%
    suspected_outliers()
  
  rdata$GSE73685=
    edata$GSE73685$unnorm %>%
    .[intersect(rownames(GPL6244),rownames(.)),] %>%
    suspected_outliers()
  
  rdata$GSE85307=
    edata$GSE85307$unnorm %>%
    .[intersect(rownames(GPL6244),rownames(.)),] %>%
    suspected_outliers()
  
  rdata$GSE86200=
    edata$GSE86200$unnorm %>%
    .[intersect(rownames(GPL10558),rownames(.)),] %>%
    suspected_outliers()
  
  rdata$GSE149437=
    edata$GSE149437$unnorm %>%
    .[intersect(rownames(GPL28460),rownames(.)),] %>%
    suspected_outliers()
  
  rdata$GSE177477=
    edata$GSE177477$unnorm %>%
    .[intersect(rownames(GPL23159),rownames(.)),] %>%
    suspected_outliers()
}
################################################################################

#####{Identify outliers of normalized expression by HClust, include=FALSE}
if(load_outlier_comp){
  cdata$GSE108497=
    edata$GSE108497$norm %>%
    .[intersect(rownames(GPL10558),rownames(.)),] %>%
    exprs() %>%
    t() %>%
    parallelDist(method='manhattan') %>%
    hclust(method='complete')
  
  cdata$GSE73685=
    edata$GSE73685$norm %>%
    .[intersect(rownames(GPL6244),rownames(.)),] %>%
    exprs() %>%
    t() %>%
    parallelDist(method='manhattan') %>%
    hclust(method='complete')
  
  cdata$GSE85307=
    edata$GSE85307$norm %>%
    .[intersect(rownames(GPL6244),rownames(.)),] %>%
    exprs() %>%
    t() %>%
    parallelDist(method='manhattan') %>%
    hclust(method='complete')
  
  cdata$GSE86200=
    edata$GSE86200$norm %>%
    .[intersect(rownames(GPL10558),rownames(.)),] %>%
    exprs() %>%
    t() %>%
    parallelDist(method='manhattan') %>%
    hclust(method='complete')
  
  cdata$GSE149437=
    edata$GSE149437$norm %>%
    .[intersect(rownames(GPL28460),rownames(.)),] %>%
    exprs() %>%
    t() %>%
    parallelDist(method='manhattan') %>%
    hclust(method='complete')
  
  cdata$GSE177477=
    edata$GSE177477$norm %>%
    .[intersect(rownames(GPL23159),rownames(.)),] %>%
    exprs() %>%
    t() %>%
    parallelDist(method='manhattan') %>%
    hclust(method='complete')
}
################################################################################

#####{HClust outliers, eval=FALSE, fig.height=10, fig.width=20, include=FALSE}
if(load_outlier_comp){
  plot(cdata$GSE108497,cex=0.6); rect.hclust(cdata$GSE108497,k=2,border=2:3)
  cutree(cdata$GSE108497,k=2) %>% table() # 2
  
  plot(cdata$GSE73685,cex=0.6); rect.hclust(cdata$GSE73685,k=3,border=2:4)
  cutree(cdata$GSE73685,k=3) %>% table() # 3
  
  plot(cdata$GSE85307,cex=0.6); rect.hclust(cdata$GSE85307,k=3,border=2:5)
  cutree(cdata$GSE85307,k=3) %>% table() # 3
  
  plot(cdata$GSE86200,cex=0.6); rect.hclust(cdata$GSE86200,k=4,border=2:5)
  cutree(cdata$GSE86200,k=4) %>% table() # 4
  
  plot(cdata$GSE149437,cex=0.6); rect.hclust(cdata$GSE149437,k=3,border=2:4)
  cutree(cdata$GSE149437,k=3) %>% table() # 3
  
  plot(cdata$GSE177477,cex=0.6); rect.hclust(cdata$GSE177477,k=2,border=2:4)
  cutree(cdata$GSE177477,k=2) %>% table() # 2
}
################################################################################

#####{Identify the outliers by both RLE and HClust, include=FALSE}
if(load_outlier_comp){
  odata$GSE108497=
    cutree(cdata$GSE108497,k=2) %>%
    .[. %in% 2] %>%
    names() %>%
    intersect(
      rdata$GSE108497 %>%
        filter(suspect_outlier=='yes') %>%
        .$sample_id %>%
        .[!duplicated(.)]
    )
  
  odata$GSE73685=
    cutree(cdata$GSE73685,k=3) %>%
    .[. %in% 3] %>%
    names() %>%
    intersect(
      rdata$GSE73685 %>%
        filter(suspect_outlier=='yes') %>%
        .$sample_id %>%
        .[!duplicated(.)]
    )
  
  odata$GSE85307=
    cutree(cdata$GSE85307,k=3) %>%
    .[. %in% 3] %>%
    names() %>%
    intersect(
      rdata$GSE85307 %>%
        filter(suspect_outlier=='yes') %>%
        .$sample_id %>%
        .[!duplicated(.)]
    )
  
  odata$GSE86200=
    cutree(cdata$GSE86200,k=4) %>%
    .[. %in% 4] %>%
    names() %>%
    intersect(
      rdata$GSE86200 %>%
        filter(suspect_outlier=='yes') %>%
        .$sample_id %>%
        .[!duplicated(.)]
    )
  
  odata$GSE149437=
    cutree(cdata$GSE149437,k=3) %>%
    .[. %in% 3] %>%
    names() %>%
    intersect(
      rdata$GSE149437 %>%
        filter(suspect_outlier=='yes') %>%
        .$sample_id %>%
        .[!duplicated(.)]
    )
  
  odata$GSE177477=
    cutree(cdata$GSE177477,k=2) %>%
    .[. %in% 2] %>%
    names() %>%
    intersect(
      rdata$GSE177477 %>%
        filter(suspect_outlier=='yes') %>%
        .$sample_id %>%
        .[!duplicated(.)]
    )
  
  saveRDS(odata,'data/odata.rds')
}else{
  odata=readRDS('data/odata.rds')
}
################################################################################






################################################################################
## Removal of low-expressed probe sets
################################################################################

#####{Determine transcript-wise medians of expression, include=FALSE}
ldata$GSE108497$row_median=
  edata$GSE108497$norm %>%
  .[rownames(.) %>% intersect(rownames(GPL10558))
    ,colnames(.) %>% .[!.%in%odata$GSE108497]] %>%
  exprs() %>%
  rowMedians() %>%
  setNames(intersect(rownames(GPL10558),rownames(edata$GSE108497$norm)))

ldata$GSE73685$row_median=
  edata$GSE73685$norm %>%
  .[rownames(.) %>% intersect(rownames(GPL6244))
    ,colnames(.) %>% .[!.%in%odata$GSE73685]] %>%
  exprs() %>%
  rowMedians() %>%
  setNames(intersect(rownames(GPL6244),rownames(edata$GSE73685$norm)))

ldata$GSE85307$row_median=
  edata$GSE85307$norm %>%
  .[rownames(.) %>% intersect(rownames(GPL6244))
    ,colnames(.) %>% .[!.%in%odata$GSE85307]] %>%
  exprs() %>%
  rowMedians() %>%
  setNames(intersect(rownames(GPL6244),rownames(edata$GSE85307$norm)))

ldata$GSE86200$row_median=
  edata$GSE86200$norm %>%
  .[rownames(.) %>% intersect(rownames(GPL10558))
    ,colnames(.) %>% .[!.%in%odata$GSE86200]] %>%
  exprs() %>%
  rowMedians() %>%
  setNames(intersect(rownames(GPL10558),rownames(edata$GSE86200$norm)))

ldata$GSE149437$row_median=
  edata$GSE149437$norm %>%
  .[rownames(.) %>% intersect(rownames(GPL28460))
    ,colnames(.) %>% .[!.%in%odata$GSE149437]] %>%
  exprs() %>%
  rowMedians() %>%
  setNames(intersect(rownames(GPL28460),rownames(edata$GSE149437$norm)))

ldata$GSE177477$row_median=
  edata$GSE177477$norm %>%
  .[rownames(.) %>% intersect(rownames(GPL23159))
    ,colnames(.) %>% .[!.%in%odata$GSE177477]] %>%
  exprs() %>%
  rowMedians() %>%
  setNames(intersect(rownames(GPL23159),rownames(edata$GSE177477$norm)))
################################################################################

#####{Identify the low-expressed transcripts, include=FALSE}
ldata$GSE108497$probe=
  ldata$GSE108497$row_median %>%
  .[. <= quantile(.,0.2)] %>%
  names()

ldata$GSE73685$probe=
  ldata$GSE73685$row_median %>%
  .[. <= quantile(.,0.2)] %>%
  names()

ldata$GSE85307$probe=
  ldata$GSE85307$row_median %>%
  .[. <= quantile(.,0.2)] %>%
  names()

ldata$GSE86200$probe=
  ldata$GSE86200$row_median %>%
  .[. <= quantile(.,0.2)] %>%
  names()

ldata$GSE149437$probe=
  ldata$GSE149437$row_median %>%
  .[. <= quantile(.,0.2)] %>%
  names()

ldata$GSE177477$probe=
  ldata$GSE177477$row_median %>%
  .[. <= quantile(.,0.2)] %>%
  names()
################################################################################






################################################################################
## Summarization from probe sets to genes
################################################################################

#####{Summarize to genes but no outliers and low-exp transcripts, include=FALSE}
# GSE108497
gdata$GSE108497$filtered=
  edata$GSE108497$norm %>%
  .[rownames(.) %>%
      .[.%in%intersect(.,rownames(GPL10558))
        & !.%in%ldata$GSE108497$probe]
    ,colnames(.) %>% .[!.%in%odata$GSE108497]]

gdata$GSE108497$hgnc=
  collapseRows(
    datET=exprs(gdata$GSE108497$filtered)
    ,rowID=rownames(gdata$GSE108497$filtered)
    ,rowGroup=
      fdata$GSE108497 %>%
      .[rownames(gdata$GSE108497$filtered),'hgnc_symbol',drop=F] %>%
      arrange(factor(hgnc_symbol,rownames(gdata$GSE108497$filtered))) %>%
      pull(hgnc_symbol)
  )

# GSE73685
gdata$GSE73685$filtered=
  edata$GSE73685$norm %>%
  .[rownames(.) %>%
      .[.%in%intersect(.,rownames(GPL6244))
        & !.%in%ldata$GSE73685$probe]
    ,colnames(.) %>% .[!.%in%odata$GSE73685]]

gdata$GSE73685$hgnc=
  collapseRows(
    datET=exprs(gdata$GSE73685$filtered)
    ,rowID=rownames(gdata$GSE73685$filtered)
    ,rowGroup=
      fdata$GSE73685 %>%
      .[rownames(gdata$GSE73685$filtered),'hgnc_symbol',drop=F] %>%
      arrange(factor(hgnc_symbol,rownames(gdata$GSE73685$filtered))) %>%
      pull(hgnc_symbol)
  )

# GSE85307
gdata$GSE85307$filtered=
  edata$GSE85307$norm %>%
  .[rownames(.) %>%
      .[.%in%intersect(.,rownames(GPL6244))
        & !.%in%ldata$GSE85307$probe]
    ,colnames(.) %>% .[!.%in%odata$GSE85307]]

gdata$GSE85307$hgnc=
  collapseRows(
    datET=exprs(gdata$GSE85307$filtered)
    ,rowID=rownames(gdata$GSE85307$filtered)
    ,rowGroup=
      fdata$GSE85307 %>%
      .[rownames(gdata$GSE85307$filtered),'hgnc_symbol',drop=F] %>%
      arrange(factor(hgnc_symbol,rownames(gdata$GSE85307$filtered))) %>%
      pull(hgnc_symbol)
  )

# GSE86200
gdata$GSE86200$filtered=
  edata$GSE86200$norm %>%
  .[rownames(.) %>%
      .[.%in%intersect(.,rownames(GPL10558))
        & !.%in%ldata$GSE86200$probe]
    ,colnames(.) %>% .[!.%in%odata$GSE86200]]

gdata$GSE86200$hgnc=
  collapseRows(
    datET=exprs(gdata$GSE86200$filtered)
    ,rowID=rownames(gdata$GSE86200$filtered)
    ,rowGroup=
      fdata$GSE86200 %>%
      .[rownames(gdata$GSE86200$filtered),'hgnc_symbol',drop=F] %>%
      arrange(factor(hgnc_symbol,rownames(gdata$GSE86200$filtered))) %>%
      pull(hgnc_symbol)
  )

# GSE149437
gdata$GSE149437$filtered=
  edata$GSE149437$norm %>%
  .[rownames(.) %>%
      .[.%in%intersect(.,rownames(GPL28460))
        & !.%in%ldata$GSE149437$probe]
    ,colnames(.) %>% .[!.%in%odata$GSE149437]]

gdata$GSE149437$hgnc=
  collapseRows(
    datET=exprs(gdata$GSE149437$filtered)
    ,rowID=rownames(gdata$GSE149437$filtered)
    ,rowGroup=
      fdata$GSE149437 %>%
      .[rownames(gdata$GSE149437$filtered),'hgnc_symbol',drop=F] %>%
      arrange(factor(hgnc_symbol,rownames(gdata$GSE149437$filtered))) %>%
      pull(hgnc_symbol)
  )

# GSE177477
gdata$GSE177477$filtered=
  edata$GSE177477$norm %>%
  .[rownames(.) %>%
      .[.%in%intersect(.,rownames(GPL23159))
        & !.%in%ldata$GSE177477$probe]
    ,colnames(.) %>% .[!.%in%odata$GSE177477]]

gdata$GSE177477$hgnc=
  collapseRows(
    datET=exprs(gdata$GSE177477$filtered)
    ,rowID=rownames(gdata$GSE177477$filtered)
    ,rowGroup=
      fdata$GSE177477 %>%
      .[rownames(gdata$GSE177477$filtered),'hgnc_symbol',drop=F] %>%
      arrange(factor(hgnc_symbol,rownames(gdata$GSE177477$filtered))) %>%
      pull(hgnc_symbol)
  )
################################################################################






################################################################################
################################################################################
################################################################################
# Selecting the common genes among all the microarray platforms
################################################################################
################################################################################
################################################################################

#####{Compile into new expression sets, include=FALSE}
# Compile each of the expression sets
GSE108497=
  ExpressionSet(
    assayData=
      gdata$GSE108497$hgnc$datETcollapsed
    ,phenoData=
      pdata$GSE108497$eset %>%
      .[,colnames(gdata$GSE108497$hgnc$datETcollapsed)] %>%
      phenoData() %>%
      pData() %>%
      AnnotatedDataFrame()
    ,featureData=
      gdata$GSE108497$hgnc$group2row %>%
      as.data.frame() %>%
      setNames(c('hgnc_symbol','probe_id')) %>%
      left_join(
        fdata$GSE108497 %>%
          rownames_to_column(var='probe_id')
        ,by=c('hgnc_symbol','probe_id')
      ) %>%
      column_to_rownames(var='hgnc_symbol') %>% 
      AnnotatedDataFrame()
    ,protocolData=
      pdata$GSE108497$eset %>%
      .[,colnames(gdata$GSE108497$hgnc$datETcollapsed)] %>%
      protocolData() %>%
      pData() %>%
      AnnotatedDataFrame()
    ,experimentData=experimentData(pdata$GSE108497$eset)
    ,annotation=annotation(pdata$GSE108497$eset)
  )

GSE73685=
  ExpressionSet(
    assayData=
      gdata$GSE73685$hgnc$datETcollapsed
    ,phenoData=
      pdata$GSE73685$eset %>%
      .[,colnames(gdata$GSE73685$hgnc$datETcollapsed)] %>%
      phenoData() %>%
      pData() %>%
      AnnotatedDataFrame()
    ,featureData=
      gdata$GSE73685$hgnc$group2row %>%
      as.data.frame() %>%
      setNames(c('hgnc_symbol','probe_id')) %>%
      left_join(
        fdata$GSE73685 %>%
          rownames_to_column(var='probe_id')
        ,by=c('hgnc_symbol','probe_id')
      ) %>%
      column_to_rownames(var='hgnc_symbol') %>% 
      AnnotatedDataFrame()
    ,protocolData=
      pdata$GSE73685$eset %>%
      .[,colnames(gdata$GSE73685$hgnc$datETcollapsed)] %>%
      protocolData() %>%
      pData() %>%
      AnnotatedDataFrame()
    ,experimentData=experimentData(pdata$GSE73685$eset)
    ,annotation=annotation(pdata$GSE73685$eset)
  )

GSE85307=
  ExpressionSet(
    assayData=
      gdata$GSE85307$hgnc$datETcollapsed
    ,phenoData=
      pdata$GSE85307$eset %>%
      .[,colnames(gdata$GSE85307$hgnc$datETcollapsed)] %>%
      phenoData() %>%
      pData() %>%
      AnnotatedDataFrame()
    ,featureData=
      gdata$GSE85307$hgnc$group2row %>%
      as.data.frame() %>%
      setNames(c('hgnc_symbol','probe_id')) %>%
      left_join(
        fdata$GSE85307 %>%
          rownames_to_column(var='probe_id')
        ,by=c('hgnc_symbol','probe_id')
      ) %>%
      column_to_rownames(var='hgnc_symbol') %>% 
      AnnotatedDataFrame()
    ,protocolData=
      pdata$GSE85307$eset %>%
      .[,colnames(gdata$GSE85307$hgnc$datETcollapsed)] %>%
      protocolData() %>%
      pData() %>%
      AnnotatedDataFrame()
    ,experimentData=experimentData(pdata$GSE85307$eset)
    ,annotation=annotation(pdata$GSE85307$eset)
  )

GSE86200=
  ExpressionSet(
    assayData=
      gdata$GSE86200$hgnc$datETcollapsed
    ,phenoData=
      pdata$GSE86200$eset %>%
      .[,colnames(gdata$GSE86200$hgnc$datETcollapsed)] %>%
      phenoData() %>%
      pData() %>%
      AnnotatedDataFrame()
    ,featureData=
      gdata$GSE86200$hgnc$group2row %>%
      as.data.frame() %>%
      setNames(c('hgnc_symbol','probe_id')) %>%
      left_join(
        fdata$GSE86200 %>%
          rownames_to_column(var='probe_id')
        ,by=c('hgnc_symbol','probe_id')
      ) %>%
      column_to_rownames(var='hgnc_symbol') %>% 
      AnnotatedDataFrame()
    ,protocolData=
      pdata$GSE86200$eset %>%
      .[,colnames(gdata$GSE86200$hgnc$datETcollapsed)] %>%
      protocolData() %>%
      pData() %>%
      AnnotatedDataFrame()
    ,experimentData=experimentData(pdata$GSE86200$eset)
    ,annotation=annotation(pdata$GSE86200$eset)
  )

GSE149437=
  ExpressionSet(
    assayData=
      gdata$GSE149437$hgnc$datETcollapsed
    ,phenoData=
      pdata$GSE149437$eset %>%
      .[,colnames(gdata$GSE149437$hgnc$datETcollapsed)] %>%
      phenoData() %>%
      pData() %>%
      AnnotatedDataFrame()
    ,featureData=
      gdata$GSE149437$hgnc$group2row %>%
      as.data.frame() %>%
      setNames(c('hgnc_symbol','probe_id')) %>%
      left_join(
        fdata$GSE149437 %>%
          rownames_to_column(var='probe_id')
        ,by=c('hgnc_symbol','probe_id')
      ) %>%
      column_to_rownames(var='hgnc_symbol') %>% 
      AnnotatedDataFrame()
    ,protocolData=
      pdata$GSE149437$eset %>%
      .[,colnames(gdata$GSE149437$hgnc$datETcollapsed)] %>%
      protocolData() %>%
      pData() %>%
      AnnotatedDataFrame()
    ,experimentData=experimentData(pdata$GSE149437$eset)
    ,annotation=annotation(pdata$GSE149437$eset)
  )

GSE177477=
  ExpressionSet(
    assayData=
      gdata$GSE177477$hgnc$datETcollapsed
    ,phenoData=
      pdata$GSE177477$eset %>%
      .[,colnames(gdata$GSE177477$hgnc$datETcollapsed)] %>%
      phenoData() %>%
      pData() %>%
      AnnotatedDataFrame()
    ,featureData=
      gdata$GSE177477$hgnc$group2row %>%
      as.data.frame() %>%
      setNames(c('hgnc_symbol','probe_id')) %>%
      left_join(
        fdata$GSE177477 %>%
          rownames_to_column(var='probe_id')
        ,by=c('hgnc_symbol','probe_id')
      ) %>%
      column_to_rownames(var='hgnc_symbol') %>% 
      AnnotatedDataFrame()
    ,protocolData=
      pdata$GSE177477$eset %>%
      .[,colnames(gdata$GSE177477$hgnc$datETcollapsed)] %>%
      protocolData() %>%
      pData() %>%
      AnnotatedDataFrame()
    ,experimentData=experimentData(pdata$GSE177477$eset)
    ,annotation=annotation(pdata$GSE177477$eset)
  )

# Take common genes among the expression sets, except that of COVID-19
source('R/take_common_genes-function.R')
GSE108497=
  GSE108497 %>%
  take_common_genes(list(GSE108497,GSE73685,GSE85307,GSE86200,GSE149437))

GSE73685=
  GSE73685 %>%
  take_common_genes(list(GSE108497,GSE73685,GSE85307,GSE86200,GSE149437))

GSE85307=
  GSE85307 %>%
  take_common_genes(list(GSE108497,GSE73685,GSE85307,GSE86200,GSE149437))

GSE86200=
  GSE86200 %>%
  take_common_genes(list(GSE108497,GSE73685,GSE85307,GSE86200,GSE149437))

GSE149437=
  GSE149437 %>%
  take_common_genes(list(GSE108497,GSE73685,GSE85307,GSE86200,GSE149437))

GSE177477=
  GSE177477 %>%
  take_common_genes(
    list(GSE108497,GSE73685,GSE85307,GSE86200,GSE149437,GSE177477)
  )
################################################################################






################################################################################
################################################################################
################################################################################
#	Derivation of maternal-fetal interface transcriptome from maternal blood
################################################################################
################################################################################
################################################################################

#####{Create a list of tissue-pairwise expression sets, include=FALSE}
surrogate=
  GSE73685 %>%
  phenoData() %>%
  pData() %>%
  select(c('title',colnames(.) %>% .[str_detect(.,':ch1')])) %>%
  setNames(
    colnames(.) %>%
      str_remove_all(':ch1') %>%
      str_replace_all('\\s','_')
  ) %>%
  rownames_to_column(var='id') %>%
  lapply(
    X=.$tissue_type %>%
      .[!duplicated(.)] %>%
      .[str_to_lower(.)!='maternal blood'] %>%
      sort()
    ,Y=.
    ,function(X,Y){
      group_by(.,subject_id) %>%
        filter(
          str_to_lower(tissue_type)==str_to_lower(X)
          | str_to_lower(tissue_type)=='maternal blood'
        ) %>%
        filter(n()==2) %>%
        ungroup() %>%
        select(id,subject_id,tissue_type) %>%
        mutate(
          tissue_type=
            tissue_type %>%
            str_to_lower() %>%
            str_replace_all('\\s+','_')
        ) %>%
        spread(tissue_type,id) %>%
        column_to_rownames(var='subject_id')
    }) %>%
  setNames(
    GSE73685 %>%
      phenoData() %>%
      pData() %>%
      select(c('title',colnames(.) %>% .[str_detect(.,':ch1')])) %>%
      setNames(
        colnames(.) %>%
          str_remove_all(':ch1') %>%
          str_replace_all('\\s','_')
      ) %>%
      .$tissue_type %>%
      .[!duplicated(.)] %>%
      .[str_to_lower(.)!='maternal blood'] %>%
      sort() %>%
      str_to_lower() %>%
      str_replace_all('\\s+','_')
  ) %>%
  lapply(function(x){
    GSE73685 %>%
      .[,c(x[[1]],x[[2]])] %>%
      `phenoData<-`(
        phenoData(.) %>%
          pData() %>%
          rownames_to_column(var='id') %>%
          left_join(
            x %>%
              rownames_to_column(var='subject_id') %>%
              gather(tissue_type,id,-subject_id) %>%
              select(-subject_id)
            ,by='id'
          ) %>%
          mutate(
            tissue_type=
              tissue_type %>%
              factor(c('maternal_blood',unique(.) %>% .[.!='maternal_blood']))
          ) %>%
          select(tissue_type,everything()) %>%
          rename(outcome=tissue_type) %>%
          column_to_rownames(var='id') %>%
          AnnotatedDataFrame()
      )
  })
################################################################################

#####{Differential expression analysis per tissue pair, include=FALSE}
source('R/conduct_dea-function.R')

dea=
  surrogate %>%
  lapply(conduct_dea,nonevent='maternal_blood',paired=T)
################################################################################

#####{Create a function to define feature and outcome of a DEG, include=FALSE}
source('R/deg-function.R')
################################################################################

#####{Compute outcomes of individual DEGs, include=FALSE}
# Compute the individual DEGs
outcome=
  dea %>%
  lapply(X=names(.),Y=.,function(X,Y){
    Y[[X]] %>%
      deg(
        X %>%
          str_replace_all('_',' ') %>%
          str_to_title()
        ,'outcome'
        ,direction=F
      )
  }) %>%
  setNames(names(dea)) %>%
  lapply(X=1,Y=.,function(X,Y){
    Z=Y %>%
      lapply(colnames) %>%
      unlist() %>%
      unique() %>%
      .[!.%in%c('tissue_type','subject_id')]
    
    Y %>%
      lapply(X=names(.),Y=.,Z=Z,function(X,Y,Z){
        K=expand.grid(subject_id=Y[[X]]$subject_id,feature=Z)
        
        Y[[X]] %>%
          gather(feature,value,-tissue_type,-subject_id) %>%
          right_join(K,by=c('subject_id','feature')) %>%
          spread(feature,value) %>%
          .[,c('tissue_type','subject_id',Z)] %>%
          filter(!is.na(tissue_type))
      })
  }) %>%
  .[[1]] %>%
  do.call(rbind,.) %>%
  unite(title,tissue_type,subject_id,sep='_')

# Join the individual DEGs to the expression set
outcome=
  GSE73685 %>%
  `phenoData<-`(
    phenoData(.) %>%
      pData() %>%
      rownames_to_column(var='id') %>%
      left_join(outcome,by='title') %>%
      column_to_rownames(var='id') %>%
      .[,colnames(.) %>% .[!.%in%colnames(phenoData(GSE73685))]] %>%
      AnnotatedDataFrame()
  ) %>%
  pData()

# Filter only genes with 80% of minor-outcome sizes, greater than or equal to 20
outcome=
  outcome %>%
  .[
    ,pull(
      outcome %>%
        rownames_to_column(var='id') %>%
        gather(feature,value,-id) %>%
        group_by(feature,value) %>%
        summarize(n=n()) %>%
        ungroup() %>%
        mutate(
          value=
            case_when(
              value==0~'nonevent'
              ,value==1~'event'
              ,TRUE~'missing'
            )
        ) %>%
        spread(value,n) %>%
        mutate(minimum=ifelse(event<=nonevent,event,nonevent)) %>%
        filter(floor(0.8*minimum)>=20) %>%
        arrange(desc(minimum))
      ,feature
    )
  ]
################################################################################

#####{Compute features of individual DEGs, include=FALSE}
# Compute the individual DEGs
features=
  dea %>%
  lapply(deg,'Maternal Blood','feature') %>%
  lapply(X=1,Y=.,function(X,Y){
    Z=Y %>%
      lapply(colnames) %>%
      unlist() %>%
      unique() %>%
      .[!.%in%c('tissue_type','subject_id')]
    
    Y %>%
      lapply(X=names(.),Y=.,Z=Z,function(X,Y,Z){
        K=expand.grid(subject_id=Y[[X]]$subject_id,feature=Z)
        
        Y[[X]] %>%
          gather(feature,value,-tissue_type,-subject_id) %>%
          right_join(K,by=c('subject_id','feature')) %>%
          spread(feature,value) %>%
          .[,c('tissue_type','subject_id',Z)] %>%
          filter(!is.na(tissue_type)) %>%
          mutate(
            tissue_type=
              X %>%
              str_replace_all('_',' ') %>%
              str_to_title()
          )
      })
  }) %>%
  .[[1]] %>%
  do.call(rbind,.) %>%
  unite(title,tissue_type,subject_id,sep='_')

# Join the individual DEGs to the expression set
features=
  GSE73685 %>%
  `phenoData<-`(
    phenoData(.) %>%
      pData() %>%
      rownames_to_column(var='id') %>%
      left_join(features,by='title') %>%
      column_to_rownames(var='id') %>%
      .[,colnames(.) %>% .[!.%in%colnames(phenoData(GSE73685))]] %>%
      AnnotatedDataFrame()
  ) %>%
  pData()
################################################################################

#####{Determine training parameters, include=FALSE}
# Create an empty list to save training parameters.
training_parameters=list()

# Define classification tuning to apply 5-fold cross validation.
training_parameters$tuning_trControl=
  trainControl(
    method='cv'
    ,number=5
    ,summaryFunction=twoClassSummary
    ,classProbs=T
    ,savePredictions=T
    ,allowParallel=T
  )

# Define classification training to apply 30-time bootstrapping.
training_parameters$final_trControl=
  trainControl(
    method='boot'
    ,number=30
    ,summaryFunction=twoClassSummary
    ,classProbs=T
    ,savePredictions=T
    ,allowParallel=T
  )
################################################################################

#####{Create functions for resampled dimensional reduction, include=FALSE}
# Create a function for fitting the resampled dimensional reduction weights
rsdr=function(data,sd_cutoff=0,state=33,sample.kind=sample.kind){
  
  set.seed(state,sample.kind=sample.kind)
  idx_list=
    round(seq(1,nrow(data),len=10)) %>%
    lapply(X=2:length(.),Y=.,function(X,Y){
      seq(ifelse(X==2,1,Y[X-1]+1),Y[X])
    })
  idx_list=
    idx_list%>%
    lapply(X=seq(10),Y=.,function(X,Y){
      Y[-X] %>% unlist
    })
  dimr_func=prcomp
  rotm_name='rotation'
  prefix='PC'
  
  dimr_exec=function(idx_list,data){
    rs_data=
      data %>%
      .[idx_list,,drop=F]
    
    rs_nzv=
      rs_data %>%
      gather() %>%
      group_by(key) %>%
      summarize(sd_value=sd(value,na.rm=T),.groups='drop') %>%
      filter(sd_value>sd_cutoff)
    
    rs_data=rs_data[,rs_nzv$key,drop=F]
    avg=summarize_all(rs_data,mean) %>% gather()
    avg=setNames(avg$value,avg$key)
    std=summarize_all(rs_data, sd) %>% gather()
    std=setNames(std$value,std$key)
    
    dmr=
      as.matrix(rs_data) %>%
      sweep(2,avg,'-') %>%
      sweep(2,std,'/') %>%
      dimr_func()
    
    dmr$rotm=
      dmr[[rotm_name]] %>%
      `dimnames<-`(list(names(avg),paste0(prefix,seq(ncol(.)))))
    
    list(avg=avg,std=std,dmr=dmr["rotm"])
  }
  
  rsdr_object=
    idx_list %>%
    pblapply(
      FUN=dimr_exec
      ,data
    )
  
  gc()
  
  range_dim=
    rsdr_object %>%
    sapply(function(x)dim(x$dmr$rotm)[2]) %>%
    range() %>%
    .[!is.na(.)] %>%
    unique()
  
  list(
    rsdr_object=rsdr_object
    ,input_dim=colnames(data), 
    rs_method='CV'
    ,rs_number='10'
    ,dr_method='PCA', 
    sd_cutoff=sd_cutoff
    ,state=state
    ,min_output_dim=range_dim[1]
    ,max_output_dim=
      range_dim %>%
      .[length(.)]
    ,dimr_func=dimr_func
    ,rotm_name=rotm_name
    ,prefix=prefix
  )
}

# Create a function to build inputs for dimension transformation
transformator=function(rsdr_object){
  avg_rotm=
    rsdr_object %>%
    composition() %>%
    as.matrix()
  
  scaler=
    rsdr_object$rsdr_object %>%
    lapply(X=seq(length(.)),Y=.,function(X,Y){
      data.frame(
        input_dim=names(Y[[X]]$avg)
        ,avg=Y[[X]]$avg
        ,std=Y[[X]]$std
      ) %>%
        `rownames<-`(NULL)
    }) %>%
    do.call(rbind, .)
  
  scaler=
    scaler %>%
    group_by(input_dim) %>%
    summarize_all(function(x)
      mean(x,na.rm=T)
    ) %>%
    ungroup() %>%
    column_to_rownames(var='input_dim')
  
  list(avg_rotm=avg_rotm,scaler=scaler)
}

# Create a function for dimension transformation
transformation=function(tidy_set,transformator_object){
  
  avg_rotm=transformator_object$avg_rotm
  scaler=transformator_object$scaler
  
  scaled_data=
    tidy_set %>%
    exprs() %>%
    t() %>%
    as.data.frame() %>%
    cbind(
      pData(tidy_set) %>%
        mutate(outcome=as.numeric(outcome=='event'))
    ) %>%
    select(outcome,everything()) %>%
    .[,rownames(avg_rotm)] %>%
    sweep(2,scaler$avg,'-') %>%
    sweep(2,scaler$std,'/') %>%
    as.matrix()
  
  converted_data=scaled_data %*% avg_rotm
  converted_data=as.data.frame(converted_data)
  
  top_n=ncol(converted_data)
  
  dr_var=apply(converted_data,2,var)
  
  pve=dr_var/sum(dr_var)
  
  top_n_pve_dimensions=
    pve %>%
    sort(decreasing=T) %>%
    names() %>%
    .[seq(top_n)]
  
  converted_data=converted_data[,top_n_pve_dimensions]
  
  ExpressionSet(
    assayData=
      converted_data %>% t()
    ,phenoData=
      tidy_set %>% 
      pData() %>%
      cbind(
        tidy_set %>%
          exprs() %>%
          t() %>%
          as.data.frame()
      ) %>% 
      AnnotatedDataFrame(
        varMetadata=
          tidy_set %>%
          phenoData() %>% 
          varMetadata() %>%
          rbind(
            tidy_set %>%
              rownames() %>%
              data.frame(rowname = ., labelDescription = NA) %>%
              column_to_rownames(var = "rowname"))
      )
    ,featureData=
      avg_rotm %>%
      as.data.frame() %>%
      .[rownames(tidy_set)[rownames(tidy_set) %in% rownames(.)]
        ,colnames(converted_data)[colnames(converted_data) %in% colnames(.)]
        ,drop=F] %>%
      t() %>%
      as.data.frame() %>% 
      AnnotatedDataFrame(
        varMetadata=
          tidy_set %>%
          fData()
      ) %>%
      `varMetadata<-`(
        varMetadata(.) %>%
          select(labelDescription,everything())
      )
    ,experimentData=
      tidy_set %>%
      `preproc<-`(
        preproc(.) %>%
          c(list(pve=pve))
      ) %>% 
      experimentData()
    ,annotation=tidy_set %>% annotation()
    ,protocolData=tidy_set %>% protocolData())
}
################################################################################

#####{Create a function to train a surrogate transcriptome model, include=FALSE}
source('R/deg_train-function.R')
################################################################################

#####{Train surrogate transcriptome models for the genes, include=FALSE}
if(run_heavy_computation){
  cat('Conduct DEG modeling\n')
  cat('Started:',as.character(now()),'\n')
  
  dir.create('data/deg_models/',showWarnings=F)
  
  # Define the method
  the_method='glmnet'
  
  # Create multiple CPU clusters for parallel computing.
  cl=detectCores()-2
  cl=makePSOCKcluster(cl)
  registerDoParallel(cl)
  
  # Close the CPU clusters and clean up memory at exit.
  on.exit(stopCluster(cl))
  on.exit(registerDoSEQ())
  on.exit(rm(cl))
  on.exit(rm(the_method))
  on.exit(gc())
  
  # Bring these packages to the CPU clusters.
  if(the_method=='glmnet'){
    clusterEvalQ(cl,{
      library('parallel')
      library('doParallel')
      library('tidyverse')
      library('Biobase')
      library('pbapply')
      library('lubridate')
      library('rsdr')
      library('caret')
      library('glmnet')
      library('gam')
      library('MLeval')
    })
  }else if(the_method=='Rborist'){
    clusterEvalQ(cl,{
      library('parallel')
      library('doParallel')
      library('tidyverse')
      library('Biobase')
      library('pbapply')
      library('lubridate')
      library('rsdr')
      library('caret')
      library('Rborist')
      library('gam')
      library('MLeval')
    })
  }else if(the_method=='gbm'){
    clusterEvalQ(cl,{
      library('parallel')
      library('doParallel')
      library('tidyverse')
      library('Biobase')
      library('pbapply')
      library('lubridate')
      library('rsdr')
      library('caret')
      library('gbm')
      library('gam')
      library('MLeval')
    })
  }
  
  deg_models=
    colnames(outcome) %>%
    pblapply(
      deg_train
      ,method='glmnet'
      ,calib_method='gamLoess'
      ,tuning_trControl=training_parameters$tuning_trControl
      ,final_trControl=training_parameters$final_trControl
      ,tuningGrid=expand.grid(alpha=seq(0,1,len=5),lambda=10^seq(-9,0,len=5))
      ,path='data/deg_models/'
      ,GSE73685=GSE73685
      ,outcome=outcome
      ,features=features
      ,sample.kind=sample.kind
      ,rsdr=rsdr
      ,transformator=transformator
      ,transformation=transformation
      ,cl=cl
    )
  
  deg_models %>%
    unlist() %>%
    saveRDS('data/deg_models.rds')
  
  cat('End:',as.character(now()))
  
  readRDS('data/deg_models.rds') %>%
    pblapply(function(x){
      readRDS(paste0('data/deg_models/',x,'.rds')) %>%
        .$eval_model %>%
        lapply(X=names(.),Y=.,function(X,Y){
          Y[[X]]$optres$Group1[c('AUC-ROC','MCC'),] %>%
            rownames_to_column(var='metric') %>%
            mutate(tissue_type=X) %>%
            rename_all(str_to_lower)
        }) %>%
        do.call(rbind,.) %>%
        mutate(deg=x)
    }) %>%
    do.call(rbind,.) %>%
    saveRDS('data/opt_pred_deg.rds')
  
  opt_pred_deg=
    readRDS('data/opt_pred_deg.rds') %>%
    filter(metric=='MCC') %>%
    select(-metric,-ci) %>%
    rename(mcc=score) %>%
    lapply(X=1,Y=.,function(X,Y){
      lapply(X=Y$tissue_type %>% .[!duplicated(.)],Y=.,function(X,Y){
        Y %>%
          filter(tissue_type==X) %>%
          select(-tissue_type)
      }) %>%
        setNames(Y$tissue_type %>% .[!duplicated(.)])
    }) %>%
    .[[1]]
  
  readRDS('data/deg_models.rds') %>%
    pblapply(function(x){
      
      y=readRDS(paste0('data/deg_models/',x,'.rds')) %>%
        .[c('transformator','model')]
      
      # Get weights from the elastic net regression model.
      z=coef(y$model$finalModel,y$model$bestTune$lambda) %>%
        as.matrix() %>%
        as.data.frame() %>%
        rownames_to_column(var='term') %>%
        
        # Clean up.
        mutate(term=str_remove_all(term,'\\`|\\(|\\)')) %>%
        setNames(c('term','estimate')) %>%
        filter(estimate!=0) %>%
        arrange(desc(abs(estimate)))
      
      k=z %>%
        filter(!term %in% c('Intercept')) %>%
        column_to_rownames(var='term')
      
      y$transformator$avg_rotm[,rownames(k),drop=F] %>%
        dimnames() %>%
        setNames(c('feature','dimension'))
      
    }) %>%
    setNames(readRDS('data/deg_models.rds')) %>%
    saveRDS('data/opt_feat_deg.rds')
  
  opt_feat_deg=
    readRDS('data/opt_feat_deg.rds') %>%
    lapply(function(x)x$feature) %>%
    Reduce(union,x=.) %>%
    data.frame(feature=.)
  
}else{
  cat(readRDS('data/log.rds')[['deg_models']])
  deg_models=readRDS('data/deg_models.rds')
  opt_pred_deg=
    readRDS('data/opt_pred_deg.rds') %>%
    filter(metric=='MCC') %>%
    select(-metric,-ci) %>%
    rename(mcc=score) %>%
    lapply(X=1,Y=.,function(X,Y){
      lapply(X=Y$tissue_type %>% .[!duplicated(.)],Y=.,function(X,Y){
        Y %>%
          filter(tissue_type==X) %>%
          select(-tissue_type)
      }) %>%
        setNames(Y$tissue_type %>% .[!duplicated(.)])
    }) %>%
    .[[1]]
  
  opt_feat_deg=
    readRDS('data/opt_feat_deg.rds') %>%
    lapply(function(x)x$feature) %>%
    Reduce(union,x=.) %>%
    data.frame(feature=.)
}
################################################################################






################################################################################
################################################################################
################################################################################
# Predictive modeling
################################################################################
################################################################################
################################################################################

################################################################################
## Construct the development and validation datasets
################################################################################

#####{Create a function for Q-Q normalization targeting controls, include=FALSE}
source('R/dea_qn-function.R')
################################################################################

#####{Create a function to derive the surrogate transcriptome, include=FALSE}
source('R/surrogate_deg-function.R')
################################################################################

#####{Create surrogate transcriptomes of GSE108497, include=FALSE}
if(run_heavy_computation){
  cat('Create surrogate transcriptomes of GSE108497\n')
  cat('Started:',as.character(now()),'\n')
  
  GSE108497_surrogate=
    opt_pred_deg %>%
    lapply(X=names(.),Y=.,Z=exprs(GSE108497),function(X,Y,Z){
      cat('Create surrogate transcriptome of',str_remove_all(X,'_'),'\n')
      L=Y[[X]] %>%
        pblapply(X=seq(nrow(.)),Y=.,Z=X,K=Z,function(X,Y,Z,K){
          L=surrogate_deg(
            K
            ,dea[Z]
            ,readRDS(paste0('data/deg_models/',Y$deg[X],'.rds'))
          )
          L[[1]]=(1+L[[1]]*Y$mcc[X])/2
          L
        }) %>%
        do.call(cbind,.)
      cat('\n')
      L
    }) %>%
    do.call(cbind,.)
  
  cat('End:',as.character(now()))
  saveRDS(GSE108497_surrogate,'data/GSE108497_surrogate.rds')
}else{
  cat(readRDS('data/log.rds')[['GSE108497_surrogate']])
  GSE108497_surrogate=readRDS('data/GSE108497_surrogate.rds')
}  
################################################################################

#####{Create surrogate transcriptomes of GSE85307, include=FALSE}
if(run_heavy_computation){
  cat('Create surrogate transcriptomes of GSE85307\n')
  cat('Started:',as.character(now()),'\n')
  
  GSE85307_surrogate=
    opt_pred_deg %>%
    lapply(X=names(.),Y=.,Z=exprs(GSE85307),function(X,Y,Z){
      cat('Create surrogate transcriptome of',str_remove_all(X,'_'),'\n')
      L=Y[[X]] %>%
        pblapply(X=seq(nrow(.)),Y=.,Z=X,K=Z,function(X,Y,Z,K){
          L=surrogate_deg(
            K
            ,dea[Z]
            ,readRDS(paste0('data/deg_models/',Y$deg[X],'.rds'))
          )
          L[[1]]=(1+L[[1]]*Y$mcc[X])/2
          L
        }) %>%
        do.call(cbind,.)
      cat('\n')
      L
    }) %>%
    do.call(cbind,.)
  
  cat('End:',as.character(now()))
  saveRDS(GSE85307_surrogate,'data/GSE85307_surrogate.rds')
}else{
  cat(readRDS('data/log.rds')[['GSE85307_surrogate']])
  GSE85307_surrogate=readRDS('data/GSE85307_surrogate.rds')
}  
################################################################################

#####{Create surrogate transcriptomes of GSE86200, include=FALSE}
if(run_heavy_computation){
  cat('Create surrogate transcriptomes of GSE86200\n')
  cat('Started:',as.character(now()),'\n')
  
  GSE86200_surrogate=
    opt_pred_deg %>%
    lapply(X=names(.),Y=.,Z=exprs(GSE86200),function(X,Y,Z){
      cat('Create surrogate transcriptome of',str_remove_all(X,'_'),'\n')
      L=Y[[X]] %>%
        pblapply(X=seq(nrow(.)),Y=.,Z=X,K=Z,function(X,Y,Z,K){
          L=surrogate_deg(
            K
            ,dea[Z]
            ,readRDS(paste0('data/deg_models/',Y$deg[X],'.rds'))
          )
          L[[1]]=(1+L[[1]]*Y$mcc[X])/2
          L
        }) %>%
        do.call(cbind,.)
      cat('\n')
      L
    }) %>%
    do.call(cbind,.)
  
  cat('End:',as.character(now()))
  saveRDS(GSE86200_surrogate,'data/GSE86200_surrogate.rds')
}else{
  cat(readRDS('data/log.rds')[['GSE86200_surrogate']])
  GSE86200_surrogate=readRDS('data/GSE86200_surrogate.rds')
}  
################################################################################

#####{Create surrogate transcriptomes of GSE149437, include=FALSE}
if(run_heavy_computation){
  cat('Create surrogate transcriptomes of GSE149437\n')
  cat('Started:',as.character(now()),'\n')
  
  GSE149437_surrogate=
    opt_pred_deg %>%
    lapply(X=names(.),Y=.,Z=exprs(GSE149437),function(X,Y,Z){
      cat('Create surrogate transcriptome of',str_remove_all(X,'_'),'\n')
      L=Y[[X]] %>%
        pblapply(X=seq(nrow(.)),Y=.,Z=X,K=Z,function(X,Y,Z,K){
          L=surrogate_deg(
            K
            ,dea[Z]
            ,readRDS(paste0('data/deg_models/',Y$deg[X],'.rds'))
          )
          L[[1]]=(1+L[[1]]*Y$mcc[X])/2
          L
        }) %>%
        do.call(cbind,.)
      cat('\n')
      L
    }) %>%
    do.call(cbind,.)
  
  cat('End:',as.character(now()))
  saveRDS(GSE149437_surrogate,'data/GSE149437_surrogate.rds')
}else{
  cat(readRDS('data/log.rds')[['GSE149437_surrogate']])
  GSE149437_surrogate=readRDS('data/GSE149437_surrogate.rds')
}  
################################################################################

#####{Conduct DEA for maternal blood genes of GSE108497, include=FALSE}
# Determine the true outcome from the phenotype data
GSE108497_phenotype=
  GSE108497 %>%
  phenoData() %>%
  pData() %>%
  select(c('title',colnames(.) %>% .[str_detect(.,':ch1')])) %>%
  setNames(colnames(.) %>% str_remove_all(':ch1')) %>%
  filter(
    apl==0
    & lac==0
    & fd==0
    & nnd==0
  ) %>%
  mutate(early=as.integer(as.numeric(ga_at_end_of_pregnancy)<34)) %>%
  mutate(time_point=ifelse(is.na(time_point),'non_pregnant',time_point)) %>%
  filter(
    time_point
    %in%c('<16 weeks','16-23 weeks','24-31 weeks','32-40 weeks')
  ) %>%
  filter(
    (early==0 & pe==0 & sga==0)
    |
      (early==1 & pe==1 & sga==0)
    |
      (early==0 & pe==1 & sga==0)
    |
      (early==0 & pe==0 & sga==1)
  ) %>%
  mutate(
    diagnosis=
      case_when(
        (early==0 & pe==0 & sga==0)~'normotensive-AGA'
        ,(early==1 & pe==1 & sga==0)~'EOPE'
        ,(early==0 & pe==1 & sga==0)~'LOPE'
        ,(early==0 & pe==0 & sga==1)~'FGR-SGA'
        ,TRUE~'other'
      )
  ) %>%
  mutate(
    outcome=
      ifelse(
        diagnosis%in%c('EOPE','LOPE')
        ,'event'
        ,'nonevent'
      ) %>%
      factor(c('event','nonevent'))
  ) %>%
  select(outcome,diagnosis,time_point,everything())

# Conduct the DEA
GSE108497_maternal_blood=
  GSE108497 %>%
  .[rownames(.) %>% .[.%in%opt_feat_deg$feature]
    ,rownames(GSE108497_phenotype)] %>%
  `phenoData<-`(
    GSE108497_phenotype %>%
      AnnotatedDataFrame()
  ) %>%
  conduct_dea(ignore_timepoint=T)
################################################################################

################################################################################

#####{Get the average expressions of GSE108497 for all Q-Q norm, include=FALSE}
GSE108497_AveExpr1=
  GSE108497_maternal_blood$result %>%
  rownames_to_column(var='feature') %>%
  filter(adj.P.Val<0.05 & feature%in%opt_feat_deg$feature) %>%
  select(feature,AveExpr1) %>%
  mutate(feature=paste0('maternal_blood',':',feature)) %>%
  rbind(
    dea %>%
      lapply(X=names(.),Y=.,function(X,Y){
        Y[[X]]$result %>%
          .[opt_pred_deg[[X]]$deg,] %>%
          rownames_to_column(var='feature') %>%
          mutate(feature=paste0(X,':',feature)) %>%
          select(feature,AveExpr1)
      }) %>%
      do.call(rbind,.)
  ) %>%
  column_to_rownames(var='feature')
################################################################################

#####{Get the features of GSE108497, include=FALSE}
GSE108497_features=
  GSE108497 %>%
  .[,rownames(GSE108497_phenotype)] %>%
  .[GSE108497_maternal_blood$result %>%
      rownames_to_column(var='feature') %>%
      filter(adj.P.Val<0.05 & feature%in%opt_feat_deg$feature) %>%
      pull(feature),] %>%
  exprs() %>%
  t() %>%
  as.data.frame() %>%
  rename_all(function(x)paste0('maternal_blood',':',x)) %>%
  cbind(
    GSE108497_surrogate %>%
      .[rownames(GSE108497_phenotype),]
  )
################################################################################

#####{Construct the development dataset, include=FALSE}
train_set=
  GSE108497_phenotype %>%
  select(outcome) %>%
  rownames_to_column(var='id') %>%
  left_join(
    GSE108497_features %>%
      t() %>%
      dea_qn(the_dea_result=GSE108497_AveExpr1) %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column(var='id')
    ,by='id'
  ) %>%
  column_to_rownames(var='id')
################################################################################

#####{DEA of all the features of the development dataset, include=FALSE}
GSE108497_all_tissue=
  ExpressionSet(
    assayData=
      train_set %>%
      select(-outcome) %>%
      t()
    ,phenoData=
      train_set %>%
      select(outcome) %>%
      rownames_to_column(var='id') %>%
      left_join(
        GSE108497_phenotype %>%
          select(time_point) %>%
          rownames_to_column(var='id')
        ,by='id'
      ) %>%
      column_to_rownames(var='id') %>%
      .[rownames(train_set),] %>%
      AnnotatedDataFrame()
  ) %>%
  conduct_dea(ignore_timepoint=T)
################################################################################

#####{DEA of blood transcripts of the development dataset, include=FALSE}
GSE108497_sel_mat_blood=
  ExpressionSet(
    assayData=
      train_set %>%
      select(-outcome) %>%
      .[,str_detect(colnames(.),'maternal_blood')] %>%
      t()
    ,phenoData=
      train_set %>%
      select(outcome) %>%
      rownames_to_column(var='id') %>%
      left_join(
        GSE108497_phenotype %>%
          select(time_point) %>%
          rownames_to_column(var='id')
        ,by='id'
      ) %>%
      column_to_rownames(var='id') %>%
      .[rownames(train_set),] %>%
      AnnotatedDataFrame()
  ) %>%
  conduct_dea(ignore_timepoint=T)
################################################################################

#####{DEA of surrogate transcripts of the development dataset, include=FALSE}
GSE108497_repro_tissue=
  ExpressionSet(
    assayData=
      train_set %>%
      select(-outcome) %>%
      .[,!str_detect(colnames(.),'maternal_blood')] %>%
      t()
    ,phenoData=
      train_set %>%
      select(outcome) %>%
      rownames_to_column(var='id') %>%
      left_join(
        GSE108497_phenotype %>%
          select(time_point) %>%
          rownames_to_column(var='id')
        ,by='id'
      ) %>%
      column_to_rownames(var='id') %>%
      .[rownames(train_set),] %>%
      AnnotatedDataFrame()
  ) %>%
  conduct_dea(ignore_timepoint=T)
################################################################################

#####{Conduct DEA for maternal blood genes of GSE85307, include=FALSE}
# Determine the true outcome from the phenotype data
GSE85307_phenotype=
  GSE85307 %>%
  phenoData() %>%
  pData() %>%
  select(c('title',colnames(.) %>% .[str_detect(.,':ch1')])) %>%
  setNames(
    colnames(.) %>%
      str_remove_all(':ch1') %>%
      str_replace_all('\\s','_')
  ) %>%
  mutate(
    blood_draw_at_gestation_week=
      as.numeric(blood_draw_at_gestation_week)
    ,time_point=case_when(
      between(ceiling(blood_draw_at_gestation_week),0,15)~'<16 weeks'
      ,between(ceiling(blood_draw_at_gestation_week),16,23)~'16-23 weeks'
      ,between(ceiling(blood_draw_at_gestation_week),24,31)~'24-31 weeks'
      ,between(ceiling(blood_draw_at_gestation_week),32,40)~'32-40 weeks'
      ,TRUE~''
    )
    ,early=ifelse(gestation_weeks<=34,'Early_Onset','Late_Onset')
  ) %>%
  unite(diagnosis,early,pregnancy_condition) %>%
  # filter(`vitamind_25(0h)_ng/ml_whole_blood_baseline`<30) %>%
  mutate(
    outcome=
      ifelse(
        diagnosis%in%c('Early_Onset_Preeclampsia','Late_Onset_Preeclampsia')
        ,'event'
        ,'nonevent'
      ) %>%
      factor(c('event','nonevent'))
  ) %>%
  select(outcome,diagnosis,time_point,everything())

# Conduct the DEA
GSE85307_maternal_blood=
  GSE85307 %>%
  .[rownames(.) %>% .[.%in%opt_feat_deg$feature]
    ,rownames(GSE85307_phenotype)] %>%
  `phenoData<-`(
    GSE85307_phenotype %>%
      AnnotatedDataFrame()
  ) %>%
  conduct_dea(ignore_timepoint=T)
################################################################################

#####{Get the features of GSE85307, include=FALSE}
GSE85307_features=
  GSE85307 %>%
  .[,rownames(GSE85307_phenotype)] %>%
  .[GSE108497_maternal_blood$result %>%
      rownames_to_column(var='feature') %>%
      filter(adj.P.Val<0.05 & feature%in%opt_feat_deg$feature) %>%
      pull(feature),] %>%
  exprs() %>%
  t() %>%
  as.data.frame() %>%
  rename_all(function(x)paste0('maternal_blood',':',x)) %>%
  cbind(
    GSE85307_surrogate %>%
      .[rownames(GSE85307_phenotype),]
  )
################################################################################

#####{Construct the replication dataset 1, include=FALSE}
test_set1=
  GSE85307_phenotype %>%
  select(outcome) %>%
  rownames_to_column(var='id') %>%
  left_join(
    GSE85307_features %>%
      t() %>%
      dea_qn(the_dea_result=GSE108497_AveExpr1) %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column(var='id')
    ,by='id'
  ) %>%
  column_to_rownames(var='id')
################################################################################

#####{DEA of all the features of the replication dataset 1, include=FALSE}
GSE85307_all_tissue=
  ExpressionSet(
    assayData=
      test_set1 %>%
      select(-outcome) %>%
      t()
    ,phenoData=
      test_set1 %>%
      select(outcome) %>%
      rownames_to_column(var='id') %>%
      left_join(
        GSE85307_phenotype %>%
          select(time_point) %>%
          rownames_to_column(var='id')
        ,by='id'
      ) %>%
      column_to_rownames(var='id') %>%
      .[rownames(test_set1),] %>%
      AnnotatedDataFrame()
  ) %>%
  conduct_dea(ignore_timepoint=T)
################################################################################

#####{DEA of blood transcripts of the replication dataset 1, include=FALSE}
GSE85307_sel_mat_blood=
  ExpressionSet(
    assayData=
      test_set1 %>%
      select(-outcome) %>%
      .[,str_detect(colnames(.),'maternal_blood')] %>%
      t()
    ,phenoData=
      test_set1 %>%
      select(outcome) %>%
      rownames_to_column(var='id') %>%
      left_join(
        GSE85307_phenotype %>%
          select(time_point) %>%
          rownames_to_column(var='id')
        ,by='id'
      ) %>%
      column_to_rownames(var='id') %>%
      .[rownames(test_set1),] %>%
      AnnotatedDataFrame()
  ) %>%
  conduct_dea(ignore_timepoint=T)
################################################################################

#####{DEA of surrogate transcripts of the replication dataset 1, include=FALSE}
GSE85307_repro_tissue=
  ExpressionSet(
    assayData=
      test_set1 %>%
      select(-outcome) %>%
      .[,!str_detect(colnames(.),'maternal_blood')] %>%
      t()
    ,phenoData=
      test_set1 %>%
      select(outcome) %>%
      rownames_to_column(var='id') %>%
      left_join(
        GSE85307_phenotype %>%
          select(time_point) %>%
          rownames_to_column(var='id')
        ,by='id'
      ) %>%
      column_to_rownames(var='id') %>%
      .[rownames(test_set1),] %>%
      AnnotatedDataFrame()
  ) %>%
  conduct_dea(ignore_timepoint=T)
################################################################################

#####{Conduct DEA for maternal blood genes of GSE86200, include=FALSE}
# Determine the true outcome from the phenotype data
GSE86200_phenotype=
  GSE86200 %>%
  phenoData() %>%
  pData() %>%
  select(c('title',colnames(.) %>% .[str_detect(.,':ch1')])) %>%
  setNames(
    colnames(.) %>%
      str_remove_all(':ch1') %>%
      str_replace_all('\\s','_')
  ) %>%
  mutate(
    blood_draw_at_gestation_week=
      floor(as.numeric(gestational_day_at_enrollment)/7)
    ,blood_draw_at_gestation_week=
      ifelse(str_detect(title,'32-38'),34,blood_draw_at_gestation_week)
    ,time_point=case_when(
      between(ceiling(blood_draw_at_gestation_week),0,15)~'<16 weeks'
      ,between(ceiling(blood_draw_at_gestation_week),16,23)~'16-23 weeks'
      ,between(ceiling(blood_draw_at_gestation_week),24,31)~'24-31 weeks'
      ,between(ceiling(blood_draw_at_gestation_week),32,40)~'32-40 weeks'
      ,TRUE~''
    )
    ,trt=
      ifelse(
        str_detect(treatment_group,'Intervention')
        ,'Treated','Untreated'
      )
  ) %>%
  unite(diagnosis,pregnancy_condition,trt) %>%
  mutate(
    outcome=
      ifelse(
        diagnosis%in%c('Preeclampsia_Untreated','Preeclampsia_Treated')
        ,'event'
        ,'nonevent'
      ) %>%
      factor(c('event','nonevent'))
  ) %>%
  select(outcome,diagnosis,time_point,everything())

# Conduct the DEA
GSE86200_maternal_blood=
  GSE86200 %>%
  .[rownames(.) %>% .[.%in%opt_feat_deg$feature]
    ,rownames(GSE86200_phenotype)] %>%
  `phenoData<-`(
    GSE86200_phenotype %>%
      AnnotatedDataFrame()
  ) %>%
  conduct_dea(ignore_timepoint=T)
################################################################################

#####{Get the features of GSE86200, include=FALSE}
GSE86200_features=
  GSE86200 %>%
  .[,rownames(GSE86200_phenotype)] %>%
  .[GSE108497_maternal_blood$result %>%
      rownames_to_column(var='feature') %>%
      filter(adj.P.Val<0.05 & feature%in%opt_feat_deg$feature) %>%
      pull(feature),] %>%
  exprs() %>%
  t() %>%
  as.data.frame() %>%
  rename_all(function(x)paste0('maternal_blood',':',x)) %>%
  cbind(
    GSE86200_surrogate %>%
      .[rownames(GSE86200_phenotype),]
  )
################################################################################

#####{Construct the replication dataset 2, include=FALSE}
test_set2=
  GSE86200_phenotype %>%
  select(outcome) %>%
  rownames_to_column(var='id') %>%
  left_join(
    GSE86200_features %>%
      t() %>%
      dea_qn(the_dea_result=GSE108497_AveExpr1) %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column(var='id')
    ,by='id'
  ) %>%
  column_to_rownames(var='id')
################################################################################

#####{DEA of all the features of the replication dataset 2, include=FALSE}
GSE86200_all_tissue=
  ExpressionSet(
    assayData=
      test_set2 %>%
      select(-outcome) %>%
      t()
    ,phenoData=
      test_set2 %>%
      select(outcome) %>%
      rownames_to_column(var='id') %>%
      left_join(
        GSE86200_phenotype %>%
          select(time_point) %>%
          rownames_to_column(var='id')
        ,by='id'
      ) %>%
      column_to_rownames(var='id') %>%
      .[rownames(test_set2),] %>%
      AnnotatedDataFrame()
  ) %>%
  conduct_dea(ignore_timepoint=T)
################################################################################

#####{DEA of blood transcripts of the replication dataset 2, include=FALSE}
GSE86200_sel_mat_blood=
  ExpressionSet(
    assayData=
      test_set2 %>%
      select(-outcome) %>%
      .[,str_detect(colnames(.),'maternal_blood')] %>%
      t()
    ,phenoData=
      test_set2 %>%
      select(outcome) %>%
      rownames_to_column(var='id') %>%
      left_join(
        GSE86200_phenotype %>%
          select(time_point) %>%
          rownames_to_column(var='id')
        ,by='id'
      ) %>%
      column_to_rownames(var='id') %>%
      .[rownames(test_set2),] %>%
      AnnotatedDataFrame()
  ) %>%
  conduct_dea(ignore_timepoint=T)
################################################################################

#####{DEA of surrogate transcripts of the replication dataset 2, include=FALSE}
GSE86200_repro_tissue=
  ExpressionSet(
    assayData=
      test_set2 %>%
      select(-outcome) %>%
      .[,!str_detect(colnames(.),'maternal_blood')] %>%
      t()
    ,phenoData=
      test_set2 %>%
      select(outcome) %>%
      rownames_to_column(var='id') %>%
      left_join(
        GSE86200_phenotype %>%
          select(time_point) %>%
          rownames_to_column(var='id')
        ,by='id'
      ) %>%
      column_to_rownames(var='id') %>%
      .[rownames(test_set2),] %>%
      AnnotatedDataFrame()
  ) %>%
  conduct_dea(ignore_timepoint=T)
################################################################################

#####{Conduct DEA for maternal blood genes of GSE149437, include=FALSE}
# Determine the true outcome from the phenotype data
GSE149437_phenotype=
  GSE149437 %>%
  phenoData() %>%
  pData() %>%
  select(c('title',colnames(.) %>% .[str_detect(.,':ch1')])) %>%
  setNames(
    colnames(.) %>%
      str_remove_all(':ch1') %>%
      str_replace_all('\\s','_')
  ) %>%
  mutate(
    blood_draw_at_gestation_week=as.numeric(gestational_age)
    ,time_point=case_when(
      between(ceiling(blood_draw_at_gestation_week),0,15)~'<16 weeks'
      ,between(ceiling(blood_draw_at_gestation_week),16,23)~'16-23 weeks'
      ,between(ceiling(blood_draw_at_gestation_week),24,31)~'24-31 weeks'
      ,between(ceiling(blood_draw_at_gestation_week),32,40)~'32-40 weeks'
      ,TRUE~''
    )
  ) %>%
  rename(diagnosis=group) %>%
  # filter(diagnosis%in%c('Early_Preeclampsia','Control')) %>%
  mutate(
    outcome=
      ifelse(
        diagnosis%in%c('Early_Preeclampsia')
        ,'event'
        ,'nonevent'
      ) %>%
      factor(c('event','nonevent'))
  ) %>%
  select(outcome,diagnosis,time_point,everything())

# Conduct the DEA
GSE149437_maternal_blood=
  GSE149437 %>%
  .[rownames(.) %>% .[.%in%opt_feat_deg$feature]
    ,rownames(GSE149437_phenotype)] %>%
  `phenoData<-`(
    GSE149437_phenotype %>%
      AnnotatedDataFrame()
  ) %>%
  conduct_dea(ignore_timepoint=T)
################################################################################

#####{Get the features of GSE149437, include=FALSE}
GSE149437_features=
  GSE149437 %>%
  .[,rownames(GSE149437_phenotype)] %>%
  .[GSE108497_maternal_blood$result %>%
      rownames_to_column(var='feature') %>%
      filter(adj.P.Val<0.05 & feature%in%opt_feat_deg$feature) %>%
      pull(feature),] %>%
  exprs() %>%
  t() %>%
  as.data.frame() %>%
  rename_all(function(x)paste0('maternal_blood',':',x)) %>%
  cbind(
    GSE149437_surrogate %>%
      .[rownames(GSE149437_phenotype),]
  )
################################################################################

#####{Construct the replication dataset 3, include=FALSE}
test_set3=
  GSE149437_phenotype %>%
  select(outcome) %>%
  rownames_to_column(var='id') %>%
  left_join(
    GSE149437_features %>%
      t() %>%
      dea_qn(the_dea_result=GSE108497_AveExpr1) %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column(var='id')
    ,by='id'
  ) %>%
  column_to_rownames(var='id')
################################################################################

#####{DEA of all the features of the replication dataset 3, include=FALSE}
GSE149437_all_tissue=
  ExpressionSet(
    assayData=
      test_set3 %>%
      select(-outcome) %>%
      t()
    ,phenoData=
      test_set3 %>%
      select(outcome) %>%
      rownames_to_column(var='id') %>%
      left_join(
        GSE149437_phenotype %>%
          select(time_point) %>%
          rownames_to_column(var='id')
        ,by='id'
      ) %>%
      column_to_rownames(var='id') %>%
      .[rownames(test_set3),] %>%
      AnnotatedDataFrame()
  ) %>%
  conduct_dea(ignore_timepoint=T)
################################################################################

#####{DEA of blood transcripts of the replication dataset 3, include=FALSE}
GSE149437_sel_mat_blood=
  ExpressionSet(
    assayData=
      test_set3 %>%
      select(-outcome) %>%
      .[,str_detect(colnames(.),'maternal_blood')] %>%
      t()
    ,phenoData=
      test_set3 %>%
      select(outcome) %>%
      rownames_to_column(var='id') %>%
      left_join(
        GSE149437_phenotype %>%
          select(time_point) %>%
          rownames_to_column(var='id')
        ,by='id'
      ) %>%
      column_to_rownames(var='id') %>%
      .[rownames(test_set3),] %>%
      AnnotatedDataFrame()
  ) %>%
  conduct_dea(ignore_timepoint=T)
################################################################################

#####{DEA of surrogate transcripts of the replication dataset 3, include=FALSE}
GSE149437_repro_tissue=
  ExpressionSet(
    assayData=
      test_set3 %>%
      select(-outcome) %>%
      .[,!str_detect(colnames(.),'maternal_blood')] %>%
      t()
    ,phenoData=
      test_set3 %>%
      select(outcome) %>%
      rownames_to_column(var='id') %>%
      left_join(
        GSE149437_phenotype %>%
          select(time_point) %>%
          rownames_to_column(var='id')
        ,by='id'
      ) %>%
      column_to_rownames(var='id') %>%
      .[rownames(test_set3),] %>%
      AnnotatedDataFrame()
  ) %>%
  conduct_dea(ignore_timepoint=T)
################################################################################

#####{Determine the true outcome from the phenotype of GSE177477, include=FALSE}
GSE177477_phenotype=
  GSE177477 %>%
  phenoData() %>%
  pData() %>%
  select(c('title',colnames(.) %>% .[str_detect(.,':ch1')])) %>%
  setNames(
    colnames(.) %>%
      str_remove_all(':ch1') %>%
      str_replace_all('\\s','_')
  ) %>%
  separate(
    comorbids
    ,c('asthma','cad','caf','t2dm','htn','ihd','sp_cabg')
    ,sep='\\|'
  ) %>%
  mutate_at(
    c('asthma','cad','caf','t2dm','htn','ihd','sp_cabg')
    ,function(x) ifelse(x=='',0,ifelse(is.na(x),0,1))
  ) %>%
  mutate(
    severity=
      factor(
        ifelse(
          is.na(severity)
          ,'uninfected'
          ,severity
        )
        ,c('uninfected','mild','severe')
      )
    ,sex=
      factor(
        ifelse(
          is.na(sex)
          ,'U'
          ,sex
        )
        ,c('U','F','M')
      )
  ) %>%
  mutate(
    outcome=
      ifelse(
        severity%in%c('severe')
        ,'event'
        ,'nonevent'
      ) %>%
      factor(c('event','nonevent'))
  ) %>%
  select(outcome,severity,sex,everything())
################################################################################

#####{Explore DEGs of all the datasets, include=FALSE}
degs_all_datasets=
  list(
    GSE108497=GSE108497_maternal_blood
    ,GSE85307=GSE85307_maternal_blood
    ,GSE86200=GSE86200_maternal_blood
    ,GSE149437=GSE149437_maternal_blood
    ,`GSE108497 and GSE149437`=
      GSE108497_maternal_blood$result %>%
      rownames_to_column(var='gene') %>%
      filter(adj.P.Val<0.05) %>%
      filter(
        gene%in%rownames(filter(
          GSE149437_maternal_blood$result
          ,adj.P.Val<0.05
        ))
      ) %>%
      column_to_rownames(var='gene') %>%
      list(result=.)
  ) %>%
  lapply(X=names(.),Y=.,function(X,Y){
    Y[[X]]$result %>%
      mutate(
        regulation=
          ifelse(
            adj.P.Val>=0.05
            ,'non-DEG'
            ,ifelse(
              logFC>=0
              ,ifelse(logFC>2,'strong up','weak up')
              ,ifelse(logFC<(-2),'strong down','weak down')
            )
          )
      ) %>%
      rownames_to_column(var='deg') %>%
      filter(regulation!='non-DEG') %>%
      mutate(set=X)
  }) %>%
  setNames(
    c('GSE108497'
      ,'GSE85307'
      ,'GSE86200'
      ,'GSE149437'
      ,'GSE108497 and GSE149437')
  )
################################################################################






################################################################################
## Model training
################################################################################

#####{Create a function to train a non-DIVNN prediction model, include=FALSE}
source('R/disease_train-function.R')
################################################################################

#####{Create a function to infer clique-extracted ontology, include=FALSE}
clixo=function(similarity
               ,alpha=0.01
               ,beta=0.5
               ,feature_name='feature'
               ,onto_prefix='CliXO'
               ,os='windows'){
  
  if(os=='windows'){
    clixo_path='clixo_0.3/clixo_0.3_windows'
  }else{
    clixo_path='clixo_0.3'
  }
  
  system(paste(c(
    'git clone https://github.com/herdiantrisufriyana/clixo_0.3/'
  ),collapse=' '))
  on.exit({
    system('rm -r clixo_0.3')
  })
  
  similarity %>%
    as.data.frame() %>%
    rownames_to_column(var='source') %>%
    gather(target,similarity,-source) %>%
    .[lower.tri(similarity) %>%
        as.logical()
      ,] %>%
    arrange(desc(similarity)) %>%
    write_tsv(paste0('clixo_0.3/input.tsv',collapse=''),col_names=F)
  
  if(os=='windows'){
    filecon=file('clixo_0.3/ontology.cx')
    on.exit({
      system('rm -r clixo_0.3')
      close(filecon)
    })
    
    suppressWarnings(system(paste(c(
      paste0(clixo_path,'/clixo')
      ,'clixo_0.3/input.tsv'
      ,alpha
      ,beta
      ,feature_name
    ),collapse=' '),intern=T)) %>%
      paste(collapse='\n') %>% 
      writeLines(filecon)
  }else{
    system(paste(c(
      'bash -c'
      ,paste0('\"',clixo_path,'/clixo')
      ,'clixo_0.3/input.tsv'
      ,alpha
      ,beta
      ,feature_name
      ,'> clixo_0.3/ontology.cx\"'
    ),collapse=' '))
  }
  
  cx=
    suppressMessages(suppressWarnings(read_tsv(
      'clixo_0.3/ontology.cx'
      ,col_names=
        c('target','source','relation','similarity')
      ,col_types=
        list(col_character(),col_character(),col_character(),col_double())
      ,skip=
        sum(
          1
          ,suppressMessages(suppressWarnings(read_tsv(
            'clixo_0.3/ontology.cx'
          ))) %>%
            setNames("column") %>%
            pull(column) %>%
            str_detect('#') %>% 
            which() %>%
            max()
        )
    ))) %>%
    mutate(
      relation=ifelse(relation=='default','is_a',feature_name)
      ,source=ifelse(relation=='is_a',paste0('CliXO:',source),source)
      ,target=paste0(onto_prefix,':',target)
    ) %>%
    select(source,target,similarity,relation)
  
  cx$target=
    paste0(
      onto_prefix
      ,':'
      ,str_pad(
        as.numeric(gsub(paste0(onto_prefix,':'),'',cx$target))
        ,str_count(max(as.numeric(gsub(paste0(onto_prefix,':'),'',cx$target))))
        ,'left'
        ,'0'
      )
    )
  
  cx$source[cx$relation=='is_a']=
    paste0(
      onto_prefix
      ,':'
      ,str_pad(
        as.numeric(gsub(
          paste0(onto_prefix,':')
          ,''
          ,cx$source[cx$relation=='is_a']
        ))
        ,str_count(max(as.numeric(gsub(
          paste0(onto_prefix,':'),'',cx$source[cx$relation=='is_a']
        ))))
        ,'left'
        ,'0'
      )
    )
  
  cx
}
################################################################################

#####{Create a function to prepare for a predictive modeling, include=FALSE}
source('R/model_component-function.R')
################################################################################

#####{DI-VNN feature selection and representation, include=FALSE}
source('R/divnn_pre_object-function.R')
################################################################################

#####{Create a function to create training function given lambda, include=FALSE}
source('R/trainer_generator-function.R')
################################################################################

#####{Create a function to transform test data for DI-VNN input, include=FALSE}
source('R/test_transformer-function.R')
################################################################################

#####{Build a function for DI-VNN evaluation, include=FALSE}
source('R/eval_divnn-function.R')
################################################################################

#####{Build a function for DI-VNN calibration and evaluation, include=FALSE}
source('R/divnn_calibrator_evaluator-function.R')
################################################################################

#####{Create a function to refresh keras backend session, include=FALSE}
source('R/refresh_session-function.R')
if(load_trained_divnn) refresh_session()
################################################################################

#####{Create empty lists for non-DIVNN prediction models, include=FALSE}
disease_model=list()
disease_model2=list()
disease_model3=list()
################################################################################

#####{Train PC-ENR using all the features, include=FALSE}
if(run_heavy_computation){
  disease_model$pc_elnet=
    train_set %>%
    .[,c('outcome'
         ,GSE108497_all_tissue$result %>%
           rownames_to_column(var='feature') %>%
           filter(adj.P.Val<0.05) %>%
           pull(feature))] %>%
    disease_train(
      method='glmnet'
      ,calib_method='gamLoess'
      ,tuning_trControl=training_parameters$tuning_trControl
      ,final_trControl=training_parameters$final_trControl
      ,tuningGrid=expand.grid(alpha=seq(0,1,len=5),lambda=10^seq(-9,0,len=5))
      ,calib_tuningGrid=
        expand.grid(span=seq(0.15,0.65,len=10),degree=seq(0,1,len=2))
      ,epv=10
      ,sample.kind=sample.kind
    )
  
  saveRDS(disease_model$pc_elnet,'data/pc_elnet.rds')
}else{
  cat(readRDS('data/log.rds')[['pc_elnet']])
  disease_model$pc_elnet=readRDS('data/pc_elnet.rds')
}
################################################################################

#####{Train PC-ENR using blood transcripts, include=FALSE}
if(run_heavy_computation){
  disease_model2$pc_elnet=
    train_set %>%
    .[,c('outcome'
         ,GSE108497_sel_mat_blood$result %>%
           rownames_to_column(var='feature') %>%
           filter(adj.P.Val<0.05) %>%
           pull(feature))] %>%
    disease_train(
      method='glmnet'
      ,calib_method='gamLoess'
      ,tuning_trControl=training_parameters$tuning_trControl
      ,final_trControl=training_parameters$final_trControl
      ,tuningGrid=expand.grid(alpha=seq(0,1,len=5),lambda=10^seq(-9,0,len=5))
      ,calib_tuningGrid=
        expand.grid(span=seq(0.15,0.65,len=10),degree=seq(0,1,len=2))
      ,epv=10
      ,sample.kind=sample.kind
    )
  
  saveRDS(disease_model2$pc_elnet,'data/pc_elnet2.rds')
}else{
  cat(readRDS('data/log.rds')[['pc_elnet2']])
  disease_model2$pc_elnet=readRDS('data/pc_elnet2.rds')
}
################################################################################

#####{Train PC-ENR using surrogate transcripts, include=FALSE}
if(run_heavy_computation){
  disease_model3$pc_elnet=
    train_set %>%
    .[,c('outcome'
         ,GSE108497_repro_tissue$result %>%
           rownames_to_column(var='feature') %>%
           filter(adj.P.Val<0.05) %>%
           pull(feature))] %>%
    disease_train(
      method='glmnet'
      ,calib_method='gamLoess'
      ,tuning_trControl=training_parameters$tuning_trControl
      ,final_trControl=training_parameters$final_trControl
      ,tuningGrid=expand.grid(alpha=seq(0,1,len=5),lambda=10^seq(-9,0,len=5))
      ,calib_tuningGrid=
        expand.grid(span=seq(0.15,0.65,len=10),degree=seq(0,1,len=2))
      ,epv=10
      ,sample.kind=sample.kind
    )
  # |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=02s
  
  saveRDS(disease_model3$pc_elnet,'data/pc_elnet3.rds')
}else{
  cat(readRDS('data/log.rds')[['pc_elnet3']])
  disease_model3$pc_elnet=readRDS('data/pc_elnet3.rds')
}
################################################################################

#####{Train PC-RF using all the features, include=FALSE}
if(run_heavy_computation){
  disease_model$pc_rf=
    train_set %>%
    .[,c('outcome'
         ,GSE108497_all_tissue$result %>%
           rownames_to_column(var='feature') %>%
           filter(adj.P.Val<0.05) %>%
           pull(feature))] %>%
    disease_train(
      method='Rborist'
      ,calib_method='gamLoess'
      ,tuning_trControl=training_parameters$tuning_trControl
      ,final_trControl=training_parameters$final_trControl
      ,tuningGrid=
        expand.grid(predFixed=seq(5,45,len=5),minNode=seq(20,100,len=5))
      ,calib_tuningGrid=
        expand.grid(span=seq(0.15,0.65,len=10),degree=seq(0,1,len=2))[-1,]
      ,epv=5
      ,sample.kind=sample.kind
    )
  
  saveRDS(disease_model$pc_rf,'data/pc_rf.rds')
}else{
  cat(readRDS('data/log.rds')[['pc_rf']])
  disease_model$pc_rf=readRDS('data/pc_rf.rds')
}
################################################################################

#####{Train PC-RF using blood transcripts, include=FALSE}
if(run_heavy_computation){
  disease_model2$pc_rf=
    train_set %>%
    .[,c('outcome'
         ,GSE108497_sel_mat_blood$result %>%
           rownames_to_column(var='feature') %>%
           filter(adj.P.Val<0.05) %>%
           pull(feature))] %>%
    disease_train(
      method='Rborist'
      ,calib_method='gamLoess'
      ,tuning_trControl=training_parameters$tuning_trControl
      ,final_trControl=training_parameters$final_trControl
      ,tuningGrid=expand.grid(predFixed=seq(5,45,len=5),minNode=seq(20,100,len=5))
      ,calib_tuningGrid=
        expand.grid(span=seq(0.15,0.65,len=10),degree=seq(0,1,len=2))[-1,]
      ,epv=5
      ,sample.kind=sample.kind
    )
  
  saveRDS(disease_model2$pc_rf,'data/pc_rf2.rds')
}else{
  cat(readRDS('data/log.rds')[['pc_rf2']])
  disease_model2$pc_rf=readRDS('data/pc_rf2.rds')
}
################################################################################

#####{Train PC-RF using surrogate transcripts, include=FALSE}
if(run_heavy_computation){
  disease_model3$pc_rf=
    train_set %>%
    .[,c('outcome'
         ,GSE108497_repro_tissue$result %>%
           rownames_to_column(var='feature') %>%
           filter(adj.P.Val<0.05) %>%
           pull(feature))] %>%
    disease_train(
      method='Rborist'
      ,calib_method='gamLoess'
      ,tuning_trControl=training_parameters$tuning_trControl
      ,final_trControl=training_parameters$final_trControl
      ,tuningGrid=
        expand.grid(predFixed=seq(5,45,len=5),minNode=seq(20,100,len=5))
      ,calib_tuningGrid=
        expand.grid(span=seq(0.15,0.65,len=10),degree=seq(0,1,len=2))[-1,]
      ,epv=5
      ,sample.kind=sample.kind
    )
  
  saveRDS(disease_model3$pc_rf,'data/pc_rf3.rds')
}else{
  cat(readRDS('data/log.rds')[['pc_rf3']])
  disease_model3$pc_rf=readRDS('data/pc_rf3.rds')
}
################################################################################

#####{Train PC-GBM using all the features, include=FALSE}
if(run_heavy_computation){
  disease_model$pc_gbm=
    train_set %>%
    .[,c('outcome'
         ,GSE108497_all_tissue$result %>%
           rownames_to_column(var='feature') %>%
           filter(adj.P.Val<0.05) %>%
           pull(feature))] %>%
    disease_train(
      method='gbm'
      ,calib_method='gamLoess'
      ,tuning_trControl=training_parameters$tuning_trControl
      ,final_trControl=training_parameters$final_trControl
      ,tuningGrid=
        expand.grid(
          interaction.depth=1
          ,shrinkage=seq(0.0005,0.05,len=25)
          ,n.minobsinnode=20
        ) %>%
        mutate(n.trees=seq(0,2500,len=26) %>% .[-1] %>% rev()) %>%
        select(n.trees,everything())
      ,calib_tuningGrid=
        expand.grid(span=seq(0.15,0.65,len=10),degree=seq(0,1,len=2))
      ,epv=10
      ,sample.kind=sample.kind
    )
  
  saveRDS(disease_model$pc_gbm,'data/pc_gbm.rds')
}else{
  cat(readRDS('data/log.rds')[['pc_gbm']])
  disease_model$pc_gbm=readRDS('data/pc_gbm.rds')
}
################################################################################

#####{Train PC-GBM using blood transcripts, include=FALSE}
if(run_heavy_computation){
  disease_model2$pc_gbm=
    train_set %>%
    .[,c('outcome'
         ,GSE108497_sel_mat_blood$result %>%
           rownames_to_column(var='feature') %>%
           filter(adj.P.Val<0.05) %>%
           pull(feature))] %>%
    disease_train(
      method='gbm'
      ,calib_method='gamLoess'
      ,tuning_trControl=training_parameters$tuning_trControl
      ,final_trControl=training_parameters$final_trControl
      ,tuningGrid=
        expand.grid(
          interaction.depth=1
          ,shrinkage=seq(0.0005,0.05,len=25)
          ,n.minobsinnode=20
        ) %>%
        mutate(n.trees=seq(0,2500,len=26) %>% .[-1] %>% rev()) %>%
        select(n.trees,everything())
      ,calib_tuningGrid=
        expand.grid(span=seq(0.15,0.65,len=10),degree=seq(0,1,len=2))
      ,epv=10
      ,sample.kind=sample.kind
    )
  
  saveRDS(disease_model2$pc_gbm,'data/pc_gbm2.rds')
}else{
  cat(readRDS('data/log.rds')[['pc_gbm2']])
  disease_model2$pc_gbm=readRDS('data/pc_gbm2.rds')
}
################################################################################

#####{Train PC-GBM using surrogate transcripts, include=FALSE}
if(run_heavy_computation){
  disease_model3$pc_gbm=
    train_set %>%
    .[,c('outcome'
         ,GSE108497_repro_tissue$result %>%
           rownames_to_column(var='feature') %>%
           filter(adj.P.Val<0.05) %>%
           pull(feature))] %>%
    disease_train(
      method='gbm'
      ,calib_method='gamLoess'
      ,tuning_trControl=training_parameters$tuning_trControl
      ,final_trControl=training_parameters$final_trControl
      ,tuningGrid=
        expand.grid(
          interaction.depth=1
          ,shrinkage=seq(0.0005,0.05,len=25)
          ,n.minobsinnode=20
        ) %>%
        mutate(n.trees=seq(0,2500,len=26) %>% .[-1] %>% rev()) %>%
        select(n.trees,everything())
      ,calib_tuningGrid=
        expand.grid(span=seq(0.15,0.65,len=10),degree=seq(0,1,len=2))
      ,epv=10
      ,sample.kind=sample.kind
    )
  
  saveRDS(disease_model3$pc_gbm,'data/pc_gbm3.rds')
}else{
  cat(readRDS('data/log.rds')[['pc_gbm3']])
  disease_model3$pc_gbm=readRDS('data/pc_gbm3.rds')
}
################################################################################

#####{Create empty lists for DIVNN prediction models, include=FALSE}
disease_divnn=list()
disease_divnn2=list()
disease_divnn3=list()
lambda=10^seq(-10,-1,len=10)
################################################################################

#####{Prepare DIVNN using all the features, include=FALSE}
disease_divnn$component=
  model_component(
    data=
      train_set %>%
      .[,c('outcome'
           ,GSE108497_all_tissue$result %>%
             rownames_to_column(var='feature') %>%
             filter(adj.P.Val<0.05) %>%
             pull(feature))]
    ,method='glmnet'
    ,calib_method='gamLoess'
    ,tuning_trControl=training_parameters$tuning_trControl
    ,final_trControl=training_parameters$final_trControl
    ,tuningGrid=expand.grid(alpha=seq(0,1,len=5),lambda=10^seq(-9,0,len=5))
    ,calib_tuningGrid=
      expand.grid(span=seq(0.15,0.65,len=10),degree=seq(0,1,len=2))
  )

if(run_heavy_computation){
  disease_divnn$pre_object=
    divnn_pre_object(
      data=
        disease_divnn$component$train %>%
        rownames_to_column(var='id') %>%
        select(-weight) %>%
        column_to_rownames(var='id')
      ,save=T
      ,path='data/disease_divnn'
    )
}else{
  cat(readRDS('data/log.rds')[['disease_divnn_pre_object']])
  disease_divnn$pre_object=
    list(
      input=readRDS('data/disease_divnn_input.rds')
      ,output=readRDS('data/disease_divnn_output.rds')
    )
}
################################################################################

#####{Hyperparameter tuning for a DI-VNN using all the features, include=FALSE}
if(run_heavy_computation){
  cat('Conduct hyperparameter tuning for DI-VNN model\n')
  cat('Started:',as.character(now()),'\n')
  
  surrogate_model=
    trainer_generator(
      disease_divnn$pre_object$output
      ,class_weight=
        disease_divnn$component$class %>%
        filter(partition=='train') %>%
        mutate(outcome=as.integer(outcome=='event')) %>%
        arrange(outcome) %>%
        lapply(X=1,Y=.,function(X,Y){
          setNames(Y$weight,Y$outcome)
        }) %>%
        .[[1]]
      ,epochs=5
      ,patience=round(5/2)
      ,batch_size=8
      ,warm_up=0.05
      ,lr=2^-6
      ,min_lr=2^-6/8
      ,tuning_mode=T
      ,verbose=1
    )
  
  disease_divnn$tuning_divnn=list()
  i=1
  for(j in seq(i,length(lambda))){
    suppressWarnings(set.seed(33,sample.kind=sample.kind))
    disease_divnn$tuning_divnn[[j]]=surrogate_model(lambda[j])
  }
  
  rm(i,j,surrogate_model)
  cat('End:',as.character(now()))
  saveRDS(disease_divnn$tuning_divnn,'data/tuning_divnn.rds')
}else{
  cat(readRDS('data/log.rds')[['tuning_divnn']])
  disease_divnn$tuning_divnn=readRDS('data/tuning_divnn.rds')
}
################################################################################

#####{Conduct modeling for a DI-VNN using all the features, include=FALSE}
if(run_heavy_computation){
  cat('Conduct modeling for DI-VNN\n')
  cat('Started:',as.character(now()),'\n')
  
  surrogate_model=
    trainer_generator(
      disease_divnn$pre_object$output
      ,path='data/ontonet'
      ,class_weight=
        disease_divnn$component$class %>%
        filter(partition=='train') %>%
        mutate(outcome=as.integer(outcome=='event')) %>%
        arrange(outcome) %>%
        lapply(X=1,Y=.,function(X,Y){
          setNames(Y$weight,Y$outcome)
        }) %>%
        .[[1]]
      ,epochs=500
      ,patience=100
      ,batch_size=8
      ,warm_up=0.05
      ,lr=2^-6
      ,min_lr=2^-6/8
      ,tuning_mode=F
      ,checkpoint=T
      ,verbose=1
    )
  suppressWarnings(set.seed(33,sample.kind=sample.kind))
  disease_divnn$modeling_divnn=surrogate_model(lambda[8])
  
  rm(surrogate_model)
  cat('End:',as.character(now()))
  save_model_weights_hdf5(
    disease_divnn$modeling_divnn$ontonet
    ,'data/ontonet.h5'
  )
  saveRDS(disease_divnn$modeling_divnn,'data/modeling_divnn.rds')
}else if(load_trained_divnn){
  cat(readRDS('data/log.rds')[['modeling_divnn']])
  disease_divnn$modeling_divnn=readRDS('data/modeling_divnn.rds')
  refresh_session()
  disease_divnn$modeling_divnn$ontonet=
    disease_divnn$pre_object$output %>%
    generator.ontonet(l2_norm=lambda[8]) %>%
    load_model_weights_hdf5('data/ontonet.h5') %>%
    compile(
      optimizer=optimizer_sgd(
        lr=disease_divnn$modeling_divnn$history$metrics$lr %>%
          .[which.max(
            disease_divnn$modeling_divnn$history$metrics$val_root_roc
          )]
        ,momentum=0.9
      )
      ,loss='mean_squared_error'
      ,loss_weights=c(rep(
        0.3/(0.3*(length(.$outputs)-1)+1),length(.$outputs)-1)
        ,1/(0.3*(length(.$outputs)-1)+1)
      )
      ,metrics=c(
        tf$keras$metrics$AUC(name='roc')
        ,tf$keras$metrics$TruePositives(name='tp')
        ,tf$keras$metrics$FalseNegatives(name='fn')
        ,tf$keras$metrics$FalsePositives(name='fp')
        ,tf$keras$metrics$TrueNegatives(name='tn')
      )
    )
}
################################################################################

#####{Prepare pre-objects of the datasets using all the features, include=FALSE}
if(run_heavy_computation){
  disease_divnn$other_outputs=
    list(
      calib=
        disease_divnn$component$calib %>%
        rownames_to_column(var='id') %>%
        select(-weight) %>%
        column_to_rownames(var='id')
      ,int_val=
        disease_divnn$component$data %>%
        .[,colnames(disease_divnn$component$train)] %>%
        rownames_to_column(var='id') %>%
        select(-weight) %>%
        column_to_rownames(var='id')
      ,test_set1=test_set1
      ,test_set2=test_set2
      ,test_set3=test_set3
    ) %>%
    lapply(X=names(.),Y=.,function(X,Y){
      cat('Create tidy set for',X,'\n')
      Z=Y[[X]] %>%
        test_transformer(
          SGD1bit_fit=disease_divnn$pre_object$input$fit
          ,similarity=disease_divnn$pre_object$input$similarity
          ,mapping=disease_divnn$pre_object$input$mapping
          ,ontology=disease_divnn$pre_object$input$ontology
          ,ranked=T
          ,dims=7
          ,decreasing=F
          ,seed_num=33
        )
      cat('\n')
      Z
    }) %>%
    setNames(c('calib','int_val','test_set1','test_set2','test_set3'))
  
  saveRDS(disease_divnn$other_outputs,'data/disease_divnn_other_outputs.rds')
}else{
  cat(readRDS('data/log.rds')[['disease_divnn_other_outputs']])
  disease_divnn$other_outputs=
    readRDS('data/disease_divnn_other_outputs.rds')
}
################################################################################

#####{Calibrate and evaluate DI-VNN using all the features, include=FALSE}
if(run_heavy_computation){
  disease_model$divnn=
    divnn_calibrator_evaluator(
      data=disease_divnn$other_outputs
      ,modeling=disease_divnn$modeling_divnn
      ,component=disease_divnn$component
      ,batch_size=8
      ,title='Calibrate and evaluate DI-VNN'
      ,dir='data'
      ,file='divnn'
    )
}else{
  cat(readRDS('data/log.rds')[['divnn']])
  disease_model$divnn=
    list(
      model=readRDS('data/divnn.rds')
      ,calib_model=readRDS('data/calib_divnn.rds')
      ,eval_model=readRDS('data/eval_divnn.rds')
    )
}
################################################################################

#####{Visualization tables using all the features, include=FALSE}
if(run_heavy_computation){
  disease_divnn$visualization=list()
  
  cat('Get ontonet visualisation table.\n')
  disease_divnn$visualization$ontonet=
    disease_divnn$pre_object$output %>%
    viz.ontonet(
      feature=F
      ,eval.results=disease_divnn$modeling_divnn$evaluation
      ,eval.metric='roc'
      ,eval.pal=c('#E64B35FF','#00A087FF')
    )
  
  cat('Get ontoarray visualisation table.\n')
  disease_divnn$visualization$ontoarray=
    disease_divnn$pre_object$output %>%
    viz.ontoarray(disease_divnn$modeling_divnn$ontonet,batch_size=8,verbose=T)
  
  saveRDS(disease_divnn$visualization,'data/visualization.rds')
}else{
  cat(readRDS('data/log.rds')[['visualization']])
  disease_divnn$visualization=readRDS('data/visualization.rds')
}
################################################################################

#####{Prepare DIVNN using blood transcripts, include=FALSE}
disease_divnn2$component=
  model_component(
    data=
      train_set %>%
      .[,c('outcome'
           ,GSE108497_sel_mat_blood$result %>%
             rownames_to_column(var='feature') %>%
             filter(adj.P.Val<0.05) %>%
             pull(feature))]
    ,method='glmnet'
    ,calib_method='gamLoess'
    ,tuning_trControl=training_parameters$tuning_trControl
    ,final_trControl=training_parameters$final_trControl
    ,tuningGrid=expand.grid(alpha=seq(0,1,len=5),lambda=10^seq(-9,0,len=5))
    ,calib_tuningGrid=
      expand.grid(span=seq(0.15,0.65,len=10),degree=seq(0,1,len=2))
  )

if(run_heavy_computation){
  disease_divnn2$pre_object=
    divnn_pre_object(
      data=
        disease_divnn2$component$train %>%
        rownames_to_column(var='id') %>%
        select(-weight) %>%
        column_to_rownames(var='id')
      ,save=T
      ,path='data/disease_divnn2'
    )
}else{
  cat(readRDS('data/log.rds')[['disease_divnn2_pre_object']])
  disease_divnn2$pre_object=
    list(
      input=readRDS('data/disease_divnn2_input.rds')
      ,output=readRDS('data/disease_divnn2_output.rds')
    )
}
################################################################################

#####{Hyperparameter tuning for a DI-VNN using blood transcripts, include=FALSE}
if(run_heavy_computation){
  cat('Conduct hyperparameter tuning for DI-VNN 2 model\n')
  cat('Started:',as.character(now()),'\n')
  
  surrogate_model=
    trainer_generator(
      disease_divnn2$pre_object$output
      ,class_weight=
        disease_divnn2$component$class %>%
        filter(partition=='train') %>%
        mutate(outcome=as.integer(outcome=='event')) %>%
        arrange(outcome) %>%
        lapply(X=1,Y=.,function(X,Y){
          setNames(Y$weight,Y$outcome)
        }) %>%
        .[[1]]
      ,epochs=5
      ,patience=round(5/2)
      ,batch_size=8
      ,warm_up=0.05
      ,lr=2^-6
      ,min_lr=2^-6/8
      ,tuning_mode=T
      ,verbose=1
    )
  
  disease_divnn2$tuning_divnn=list()
  i=1
  for(j in seq(i,length(lambda))){
    suppressWarnings(set.seed(33,sample.kind=sample.kind))
    disease_divnn2$tuning_divnn[[j]]=surrogate_model(lambda[j])
  }
  
  rm(i,j,surrogate_model)
  cat('End:',as.character(now()))
  saveRDS(disease_divnn2$tuning_divnn,'data/tuning_divnn2.rds')
}else{
  cat(readRDS('data/log.rds')[['tuning_divnn2']])
  disease_divnn2$tuning_divnn=readRDS('data/tuning_divnn2.rds')
}
################################################################################

#####{Conduct modeling for a DI-VNN using blood transcripts, include=FALSE}
if(run_heavy_computation){
  cat('Conduct modeling for DI-VNN 2\n')
  cat('Started:',as.character(now()),'\n')
  
  surrogate_model=
    trainer_generator(
      disease_divnn2$pre_object$output
      ,path='data/ontonet2'
      ,class_weight=
        disease_divnn2$component$class %>%
        filter(partition=='train') %>%
        mutate(outcome=as.integer(outcome=='event')) %>%
        arrange(outcome) %>%
        lapply(X=1,Y=.,function(X,Y){
          setNames(Y$weight,Y$outcome)
        }) %>%
        .[[1]]
      ,epochs=500
      ,patience=100
      ,batch_size=8
      ,warm_up=0.05
      ,lr=2^-6
      ,min_lr=2^-6/8
      ,tuning_mode=F
      ,checkpoint=T
      ,verbose=1
    )
  suppressWarnings(set.seed(33,sample.kind=sample.kind))
  disease_divnn2$modeling_divnn=surrogate_model(lambda[10])
  
  rm(surrogate_model)
  cat('End:',as.character(now()))
  save_model_weights_hdf5(
    disease_divnn2$modeling_divnn$ontonet
    ,'data/ontonet2.h5'
  )
  saveRDS(disease_divnn2$modeling_divnn,'data/modeling_divnn2.rds')
}else if(load_trained_divnn){
  cat(readRDS('data/log.rds')[['modeling_divnn2']])
  disease_divnn2$modeling_divnn=readRDS('data/modeling_divnn2.rds')
  refresh_session()
  disease_divnn2$modeling_divnn$ontonet=
    disease_divnn2$pre_object$output %>%
    generator.ontonet(l2_norm=lambda[10]) %>%
    load_model_weights_hdf5('data/ontonet2.h5') %>%
    compile(
      optimizer=optimizer_sgd(
        lr=disease_divnn2$modeling_divnn$history$metrics$lr %>%
          .[which.max(
            disease_divnn2$modeling_divnn$history$metrics$val_root_roc
          )]
        ,momentum=0.9
      )
      ,loss='mean_squared_error'
      ,loss_weights=c(rep(
        0.3/(0.3*(length(.$outputs)-1)+1),length(.$outputs)-1)
        ,1/(0.3*(length(.$outputs)-1)+1)
      )
      ,metrics=c(
        tf$keras$metrics$AUC(name='roc')
        ,tf$keras$metrics$TruePositives(name='tp')
        ,tf$keras$metrics$FalseNegatives(name='fn')
        ,tf$keras$metrics$FalsePositives(name='fp')
        ,tf$keras$metrics$TrueNegatives(name='tn')
      )
    )
}
################################################################################

#####{Prepare pre-objects of the datasets using blood transcripts,include=FALSE}
if(run_heavy_computation){
  disease_divnn2$other_outputs=
    list(
      calib=
        disease_divnn2$component$calib %>%
        rownames_to_column(var='id') %>%
        select(-weight) %>%
        column_to_rownames(var='id')
      ,int_val=
        disease_divnn2$component$data %>%
        .[,colnames(disease_divnn2$component$train)] %>%
        rownames_to_column(var='id') %>%
        select(-weight) %>%
        column_to_rownames(var='id')
      ,test_set1=test_set1
      ,test_set2=test_set2
      ,test_set3=test_set3
    ) %>%
    lapply(X=names(.),Y=.,function(X,Y){
      cat('Create tidy set for',X,'\n')
      Z=Y[[X]] %>%
        test_transformer(
          SGD1bit_fit=disease_divnn2$pre_object$input$fit
          ,similarity=disease_divnn2$pre_object$input$similarity
          ,mapping=disease_divnn2$pre_object$input$mapping
          ,ontology=disease_divnn2$pre_object$input$ontology
          ,ranked=T
          ,dims=7
          ,decreasing=F
          ,seed_num=33
        )
      cat('\n')
      Z
    }) %>%
    setNames(c('calib','int_val','test_set1','test_set2','test_set3'))
  
  saveRDS(disease_divnn2$other_outputs,'data/disease_divnn2_other_outputs.rds')
}else{
  cat(readRDS('data/log.rds')[['disease_divnn2_other_outputs']])
  disease_divnn2$other_outputs=
    readRDS('data/disease_divnn2_other_outputs.rds')
}
################################################################################

#####{Calibrate and evaluate DI-VNN using blood transcripts, include=FALSE}
if(run_heavy_computation){
  disease_model2$divnn=
    divnn_calibrator_evaluator(
      data=disease_divnn2$other_outputs
      ,modeling=disease_divnn2$modeling_divnn
      ,component=disease_divnn2$component
      ,batch_size=8
      ,title='Calibrate and evaluate DI-VNN 2'
      ,dir='data'
      ,file='divnn2'
    )
}else{
  cat(readRDS('data/log.rds')[['divnn2']])
  disease_model2$divnn=
    list(
      model=readRDS('data/divnn2.rds')
      ,calib_model=readRDS('data/calib_divnn2.rds')
      ,eval_model=readRDS('data/eval_divnn2.rds')
    )
}
################################################################################

#####{Visualization tables using blood transcripts, include=FALSE}
if(run_heavy_computation){
  disease_divnn2$visualization=list()
  
  cat('Get ontonet 2 visualisation table.\n')
  disease_divnn2$visualization$ontonet=
    disease_divnn2$pre_object$output %>%
    viz.ontonet(
      feature=F
      ,eval.results=disease_divnn2$modeling_divnn$evaluation
      ,eval.metric='roc'
      ,eval.pal=c('#E64B35FF','#00A087FF')
    )
  
  cat('Get ontoarray 2 visualisation table.\n')
  disease_divnn2$visualization$ontoarray=
    disease_divnn2$pre_object$output %>%
    viz.ontoarray(disease_divnn2$modeling_divnn$ontonet,batch_size=8,verbose=T)
  
  saveRDS(disease_divnn2$visualization,'data/visualization2.rds')
}else{
  cat(readRDS('data/log.rds')[['visualization2']])
  disease_divnn2$visualization=readRDS('data/visualization2.rds')
}
################################################################################

#####{Prepare DIVNN using surrogate transcripts, include=FALSE}
disease_divnn3$component=
  model_component(
    data=
      train_set %>%
      .[,c('outcome'
           ,GSE108497_repro_tissue$result %>%
             rownames_to_column(var='feature') %>%
             filter(adj.P.Val<0.05) %>%
             pull(feature))]
    ,method='glmnet'
    ,calib_method='gamLoess'
    ,tuning_trControl=training_parameters$tuning_trControl
    ,final_trControl=training_parameters$final_trControl
    ,tuningGrid=expand.grid(alpha=seq(0,1,len=5),lambda=10^seq(-9,0,len=5))
    ,calib_tuningGrid=
      expand.grid(span=seq(0.15,0.65,len=10),degree=seq(0,1,len=2))
  )

if(run_heavy_computation){
  disease_divnn3$pre_object=
    divnn_pre_object(
      data=
        disease_divnn3$component$train %>%
        rownames_to_column(var='id') %>%
        select(-weight) %>%
        column_to_rownames(var='id')
      ,save=T
      ,path='data/disease_divnn3'
    )
}else{
  cat(readRDS('data/log.rds')[['disease_divnn3_pre_object']])
  disease_divnn3$pre_object=
    list(
      input=readRDS('data/disease_divnn3_input.rds')
      ,output=readRDS('data/disease_divnn3_output.rds')
    )
}
################################################################################

#####{Hyperparameter tuning for DIVNN using surrogate transcript, include=FALSE}
if(run_heavy_computation){
  cat('Conduct hyperparameter tuning for DI-VNN 3 model\n')
  cat('Started:',as.character(now()),'\n')
  
  surrogate_model=
    trainer_generator(
      disease_divnn3$pre_object$output
      ,class_weight=
        disease_divnn3$component$class %>%
        filter(partition=='train') %>%
        mutate(outcome=as.integer(outcome=='event')) %>%
        arrange(outcome) %>%
        lapply(X=1,Y=.,function(X,Y){
          setNames(Y$weight,Y$outcome)
        }) %>%
        .[[1]]
      ,epochs=5
      ,patience=round(5/2)
      ,batch_size=8
      ,warm_up=0.05
      ,lr=2^-6
      ,min_lr=2^-6/8
      ,tuning_mode=T
      ,verbose=1
    )
  
  disease_divnn3$tuning_divnn=list()
  i=1
  for(j in seq(i,length(lambda))){
    suppressWarnings(set.seed(33,sample.kind=sample.kind))
    disease_divnn3$tuning_divnn[[j]]=surrogate_model(lambda[j])
  }
  
  rm(i,j,surrogate_model)
  cat('End:',as.character(now()))
  saveRDS(disease_divnn3$tuning_divnn,'data/tuning_divnn3.rds')
}else{
  cat(readRDS('data/log.rds')[['tuning_divnn3']])
  disease_divnn3$tuning_divnn=readRDS('data/tuning_divnn3.rds')
}
################################################################################

#####{Conduct modeling for a DI-VNN using surrogate transcripts, include=FALSE}
if(run_heavy_computation){
  cat('Conduct modeling for DI-VNN 3\n')
  cat('Started:',as.character(now()),'\n')
  
  surrogate_model=
    trainer_generator(
      disease_divnn3$pre_object$output
      ,path='data/ontonet3'
      ,class_weight=
        disease_divnn3$component$class %>%
        filter(partition=='train') %>%
        mutate(outcome=as.integer(outcome=='event')) %>%
        arrange(outcome) %>%
        lapply(X=1,Y=.,function(X,Y){
          setNames(Y$weight,Y$outcome)
        }) %>%
        .[[1]]
      ,epochs=500
      ,patience=100
      ,batch_size=8
      ,warm_up=0.05
      ,lr=2^-6
      ,min_lr=2^-6/8
      ,tuning_mode=F
      ,checkpoint=T
      ,verbose=1
    )
  suppressWarnings(set.seed(33,sample.kind=sample.kind))
  disease_divnn3$modeling_divnn=surrogate_model(lambda[8])
  
  rm(surrogate_model)
  cat('End:',as.character(now()))
  save_model_weights_hdf5(
    disease_divnn3$modeling_divnn$ontonet
    ,'data/ontonet3.h5'
  )
  saveRDS(disease_divnn3$modeling_divnn,'data/modeling_divnn3.rds')
}else if(load_trained_divnn){
  cat(readRDS('data/log.rds')[['modeling_divnn3']])
  disease_divnn3$modeling_divnn=readRDS('data/modeling_divnn3.rds')
  refresh_session()
  disease_divnn3$modeling_divnn$ontonet=
    disease_divnn3$pre_object$output %>%
    generator.ontonet(l2_norm=lambda[8]) %>%
    load_model_weights_hdf5('data/ontonet3.h5') %>%
    compile(
      optimizer=optimizer_sgd(
        lr=disease_divnn3$modeling_divnn$history$metrics$lr %>%
          .[which.max(
            disease_divnn3$modeling_divnn$history$metrics$val_root_roc
          )]
        ,momentum=0.9
      )
      ,loss='mean_squared_error'
      ,loss_weights=c(rep(
        0.3/(0.3*(length(.$outputs)-1)+1),length(.$outputs)-1)
        ,1/(0.3*(length(.$outputs)-1)+1)
      )
      ,metrics=c(
        tf$keras$metrics$AUC(name='roc')
        ,tf$keras$metrics$TruePositives(name='tp')
        ,tf$keras$metrics$FalseNegatives(name='fn')
        ,tf$keras$metrics$FalsePositives(name='fp')
        ,tf$keras$metrics$TrueNegatives(name='tn')
      )
    )
}
################################################################################

#####{Prepare pre-objects using surrogate transcripts, include=FALSE}
if(run_heavy_computation){
  disease_divnn3$other_outputs=
    list(
      calib=
        disease_divnn3$component$calib %>%
        rownames_to_column(var='id') %>%
        select(-weight) %>%
        column_to_rownames(var='id')
      ,int_val=
        disease_divnn3$component$data %>%
        .[,colnames(disease_divnn3$component$train)] %>%
        rownames_to_column(var='id') %>%
        select(-weight) %>%
        column_to_rownames(var='id')
      ,test_set1=test_set1
      ,test_set2=test_set2
      ,test_set3=test_set3
    ) %>%
    lapply(X=names(.),Y=.,function(X,Y){
      cat('Create tidy set for',X,'\n')
      Z=Y[[X]] %>%
        test_transformer(
          SGD1bit_fit=disease_divnn3$pre_object$input$fit
          ,similarity=disease_divnn3$pre_object$input$similarity
          ,mapping=disease_divnn3$pre_object$input$mapping
          ,ontology=disease_divnn3$pre_object$input$ontology
          ,ranked=T
          ,dims=7
          ,decreasing=F
          ,seed_num=33
        )
      cat('\n')
      Z
    }) %>%
    setNames(c('calib','int_val','test_set1','test_set2','test_set3'))
  
  saveRDS(disease_divnn3$other_outputs,'data/disease_divnn3_other_outputs.rds')
}else{
  cat(readRDS('data/log.rds')[['disease_divnn3_other_outputs']])
  disease_divnn3$other_outputs=
    readRDS('data/disease_divnn3_other_outputs.rds')
}
################################################################################

#####{Calibrate and evaluate DI-VNN using surrogate transcripts, include=FALSE}
if(run_heavy_computation){
  disease_model3$divnn=
    divnn_calibrator_evaluator(
      data=disease_divnn3$other_outputs
      ,modeling=disease_divnn3$modeling_divnn
      ,component=disease_divnn3$component
      ,batch_size=8
      ,title='Calibrate and evaluate DI-VNN 3'
      ,dir='data'
      ,file='divnn3'
    )
}else{
  cat(readRDS('data/log.rds')[['divnn3']])
  disease_model3$divnn=
    list(
      model=readRDS('data/divnn3.rds')
      ,calib_model=readRDS('data/calib_divnn3.rds')
      ,eval_model=readRDS('data/eval_divnn3.rds')
    )
}
################################################################################

#####{Visualization tables using surrogate transcripts, include=FALSE}
if(run_heavy_computation){
  disease_divnn3$visualization=list()
  
  cat('Get ontonet 3 visualisation table.\n')
  disease_divnn3$visualization$ontonet=
    disease_divnn3$pre_object$output %>%
    viz.ontonet(
      feature=F
      ,eval.results=disease_divnn3$modeling_divnn$evaluation
      ,eval.metric='roc'
      ,eval.pal=c('#E64B35FF','#00A087FF')
    )
  
  cat('Get ontoarray 3 visualisation table.\n')
  disease_divnn3$visualization$ontoarray=
    disease_divnn3$pre_object$output %>%
    viz.ontoarray(disease_divnn3$modeling_divnn$ontonet,batch_size=8,verbose=T)
  
  saveRDS(disease_divnn3$visualization,'data/visualization3.rds')
}else{
  cat(readRDS('data/log.rds')[['visualization3']])
  disease_divnn3$visualization=readRDS('data/visualization3.rds')
}
################################################################################






################################################################################
## Model evaluation
################################################################################

#####{Create a function to get model weights, include=FALSE}
source('R/opt_feat_disease_fn-function.R')
################################################################################

#####{Get weights of the models using all the features, include=FALSE}
opt_feat_disease=
  disease_model %>%
  lapply(X=names(.),Y=.,Z=disease_divnn,opt_feat_disease_fn) %>%
  setNames(names(disease_model))
################################################################################

#####{Get weights of the models using blood transcripts, include=FALSE}
opt_feat_disease2=
  disease_model2 %>%
  lapply(X=names(.),Y=.,Z=disease_divnn2,opt_feat_disease_fn) %>%
  setNames(names(disease_model))
################################################################################

#####{Get weights of the models using surrogate features, include=FALSE}
opt_feat_disease3=
  disease_model3 %>%
  lapply(X=names(.),Y=.,Z=disease_divnn3,opt_feat_disease_fn) %>%
  setNames(names(disease_model))
################################################################################

#####{Create a function to evaluate the models, include=FALSE}
source('R/opt_pred_disease_fn-function.R')
################################################################################

#####{Evaluate using all the features in the replication 1, include=FALSE}
opt_pred_disease1=
  disease_model %>%
  lapply(X=names(.),Y=.,Z=test_set1,Z2='test_set1',opt_pred_disease_fn) %>%
  setNames(names(disease_model))
################################################################################

#####{Evaluate using blood transcripts in the replication 1, include=FALSE}
opt_pred2_disease1=
  disease_model2 %>%
  lapply(X=names(.),Y=.,Z=test_set1,Z2='test_set1',opt_pred_disease_fn) %>%
  setNames(names(disease_model))
################################################################################

#####{Evaluate using surrogate transcripts in the replication 1, include=FALSE}
opt_pred3_disease1=
  disease_model3 %>%
  lapply(X=names(.),Y=.,Z=test_set1,Z2='test_set1',opt_pred_disease_fn) %>%
  setNames(names(disease_model))
################################################################################

#####{Evaluate using all the features in the replication 2, include=FALSE}
opt_pred_disease2=
  disease_model %>%
  lapply(X=names(.),Y=.,Z=test_set2,Z2='test_set2',opt_pred_disease_fn) %>%
  setNames(names(disease_model))
################################################################################

#####{Evaluate using blood transcripts in the replication 2, include=FALSE}
opt_pred2_disease2=
  disease_model2 %>%
  lapply(X=names(.),Y=.,Z=test_set2,Z2='test_set2',opt_pred_disease_fn) %>%
  setNames(names(disease_model))
################################################################################

#####{Evaluate using surrogate transcripts in the replication 2, include=FALSE}
opt_pred3_disease2=
  disease_model3 %>%
  lapply(X=names(.),Y=.,Z=test_set2,Z2='test_set2',opt_pred_disease_fn) %>%
  setNames(names(disease_model))
################################################################################

#####{Evaluate using all the features in the replication 3, include=FALSE}
opt_pred_disease3=
  disease_model %>%
  lapply(X=names(.),Y=.,Z=test_set3,Z2='test_set3',opt_pred_disease_fn) %>%
  setNames(names(disease_model))
################################################################################

#####{Evaluate using blood transcripts in the replication 3, include=FALSE}
opt_pred2_disease3=
  disease_model2 %>%
  lapply(X=names(.),Y=.,Z=test_set3,Z2='test_set3',opt_pred_disease_fn) %>%
  setNames(names(disease_model))
################################################################################

#####{Evaluate using surrogate transcripts in the replication 3, include=FALSE}
opt_pred3_disease3=
  disease_model3 %>%
  lapply(X=names(.),Y=.,Z=test_set3,Z2='test_set3',opt_pred_disease_fn) %>%
  setNames(names(disease_model))
################################################################################

#####{Create a function to compare results of the evaluation, include=FALSE}
source('R/result_comparison-function.R')
################################################################################

#####{Compare the AUC-ROCs, echo=FALSE, fig.height=5.34309, fig.width=7.48031}
auc_comparison=result_comparison('AUC-ROC')
################################################################################

#####{Compare the PPVs, include=FALSE}
ppv_comparison=result_comparison('PREC')
################################################################################

#####{Compare the NPVs, include=FALSE}
npv_comparison=result_comparison('NPV')
################################################################################






################################################################################
################################################################################
################################################################################
# Biomarker emulation
################################################################################
################################################################################
################################################################################

#####{Get PC weights from the best model, include=FALSE}
if(run_heavy_computation){
  bestmod_raw_weights=
    disease_model3$pc_gbm$model %>%
    varImp() %>%
    .[[1]] %>%
    rownames_to_column(var='term') %>%
    
    # Clean up.
    filter(Overall!=0) %>%
    arrange(desc(Overall)) %>%
    
    left_join(
      disease_model3$pc_gbm$transformator$avg_rotm %>%
        as.data.frame() %>%
        rownames_to_column(var='feature') %>%
        gather(term,weight,-feature)
      ,by='term'
    ) %>%
    select(-Overall) %>%
    separate(feature,c('tissue','gene'),sep=':') %>%
    mutate(w_direction=ifelse(weight>=0,'positive','negative')) %>%
    group_by(term,tissue,w_direction) %>%
    arrange(term,tissue,w_direction,desc(abs(weight))) %>%
    ungroup() %>%
    left_join(
      pull(.,gene) %>%
        .[!duplicated(.)] %>%
        pblapply(function(x){
          y=readRDS(paste0('data/deg_models/',x,'.rds'))
          
          coef(y$model$finalModel,y$model$bestTune$lambda) %>%
            as.matrix() %>%
            as.data.frame() %>%
            rownames_to_column(var='term') %>%
            
            # Clean up.
            mutate(term=str_remove_all(term,'\\`|\\(|\\)')) %>%
            setNames(c('term','estimate')) %>%
            filter(estimate!=0) %>%
            arrange(desc(abs(estimate))) %>%
            filter(term!='Intercept') %>%
            left_join(
              y$transformator$avg_rotm %>%
                t() %>%
                as.data.frame() %>%
                rownames_to_column(var='term')
              ,by='term'
            ) %>%
            gather(blood_t,bt_weight,-term,-estimate) %>%
            rename(surro_pc=term,spc_weight=estimate) %>%
            mutate(gene=x) %>%
            select(gene,blood_t,bt_weight,surro_pc,spc_weight)
        }) %>%
        do.call(rbind,.)
      ,by='gene'
    ) %>%
    left_join(
      opt_pred_deg %>%
        do.call(rbind,.) %>%
        rownames_to_column(var='tissue') %>%
        separate(tissue,c('tissue','seq'),sep='\\.') %>%
        select(-seq) %>%
        rename(gene=deg)
      ,by=c('tissue','gene')
    )
  
  bestmod_raw_weights %>%
    pull(tissue) %>%
    .[!duplicated(.)] %>%
    sapply(function(x){
      bestmod_raw_weights %>%
        filter(tissue==x) %>%
        saveRDS(paste0('data/bestmod_raw_weights_',x,'.rds'))
    })
}else if(show_best_tree){
  cat(readRDS('data/log.rds')[['bestmod_raw_weights']])
  bestmod_raw_weights=
    list.files('data/',pattern='bestmod_raw_weights_',full.names=T) %>%
    lapply(readRDS) %>%
    do.call(rbind,.) %>%
    arrange(term,tissue)
}
################################################################################

#####{Get PC weights from the best model, include=FALSE}
if(show_best_tree){
  bestmod_weights=
    bestmod_raw_weights %>%
    group_by(term,tissue,gene,weight,w_direction,mcc,blood_t) %>%
    summarize(accu_w=sum(bt_weight*spc_weight,na.rm=T),n=n()) %>%
    ungroup() %>%
    group_by(tissue,blood_t) %>%
    summarize(accu_w=sum(accu_w*mcc*weight,na.rm=T),n=sum(n,na.rm=T)) %>%
    ungroup() %>%
    mutate(avg_accu_w=accu_w/n) %>%
    arrange(desc(abs(avg_accu_w))) %>%
    group_by(tissue) %>%
    mutate(rank=seq(n())) %>%
    ungroup()
}

# 5 biomarkers is considerably low-cost
# 20 is minimum rank per tissue to include all of them
# after filtering top 5 in any tissues
if(run_heavy_computation & show_best_tree){
  bestmod_list=
    seq(20) %>%
    pblapply(function(x){
      seq(5) %>%
        lapply(X=.,Y=x,function(X,Y){
          bestmod_raw_weights %>%
            left_join(bestmod_weights,by=c('tissue','blood_t')) %>%
            filter(rank%in%seq(Y)) %>%
            select(tissue,gene,blood_t,avg_accu_w,rank) %>%
            filter(!duplicated(.)) %>%
            group_by(blood_t) %>%
            mutate(n_tissue=length(tissue %>% .[!duplicated(.)])) %>%
            arrange(desc(n_tissue),tissue,rank) %>%
            ungroup() %>%
            left_join(
              select(.,blood_t) %>%
                filter(!duplicated(.)) %>%
                mutate(global_rank=seq(n()))
              ,by='blood_t'
            ) %>%
            filter(global_rank<=X) %>%
            mutate(grank=X)
        }) %>%
        do.call(rbind,.) %>%
        mutate(rank=x)
    }) %>%
    do.call(rbind,.)
  saveRDS(bestmod_list,'data/bestmod_list.rds')
}else{
  cat(readRDS('data/log.rds')[['bestmod_list']])
  bestmod_list=readRDS('data/bestmod_list.rds')
}
################################################################################

#####{Create a function to preprocess the features, include=FALSE}
source('R/mb_predictor-function.R')
################################################################################

#####{Create a function to standardize features ignoring outlier, include=FALSE}
source('R/std_no_outlier-function.R')
################################################################################

#####{Create a function to re-evaluate the biomarker tree leaves, include=FALSE}
source('R/re_evalm-function.R')
################################################################################

#####{Determine metric and reference to choose thresholds}
th_selection=list(c('sens','spec'))
th_values=list(NULL)
################################################################################

#####{Get features to train biomarker tree and evaluate it, include=FALSE}
if(run_heavy_computation){
  # Preprocess the features per combination in the development dataset
  bestmod_dev=
    bestmod_list %>%
    pull(rank) %>%
    .[!duplicated(.)] %>%
    sort() %>%
    pblapply(X=.
             ,Y=GSE108497_phenotype %>%
               select(outcome) %>%
               rownames_to_column(var='id') %>%
               left_join(
                 GSE108497 %>%
                   exprs() %>%
                   mb_predictor(
                     list(maternal_blood=GSE108497_maternal_blood)
                   ) %>%
                   rownames_to_column(var='id')
                 ,by='id'
               ) %>%
               column_to_rownames(var='id')
             ,function(X,Y){
               bestmod_list %>%
                 filter(rank==X) %>%
                 pull(grank) %>%
                 .[!duplicated(.)] %>%
                 sort() %>%
                 lapply(X=.,Y=X,Z=Y,function(X,Y,Z){
                   Z %>%
                     select_at(c(
                       'outcome'
                       ,paste0(
                         'maternal_blood:'
                         ,bestmod_list$blood_t[bestmod_list$grank==X & bestmod_list$rank==Y]
                       )
                     ))
                 }) %>%
                 setNames(paste0(
                   'grank'
                   ,bestmod_list %>%
                     filter(rank==X) %>%
                     pull(grank) %>%
                     .[!duplicated(.)] %>%
                     sort()
                 ))
             }) %>%
    setNames(paste0(
      'rank'
      ,bestmod_list %>%
        pull(rank) %>%
        .[!duplicated(.)] %>%
        sort()
    ))
  
  saveRDS(bestmod_dev,'data/bestmod_dev.rds')
  
  # Create multiple CPU clusters for parallel computing.
  cl=detectCores()-2
  cl=makePSOCKcluster(cl)
  registerDoParallel(cl)
  
  # Close the CPU clusters and clean up memory at exit.
  on.exit(stopCluster(cl))
  on.exit(registerDoSEQ())
  on.exit(rm(cl))
  on.exit(rm(the_method))
  on.exit(gc())
  
  # Bring these packages to the CPU clusters.
  clusterEvalQ(cl,{
    library('parallel')
    library('doParallel')
    library('tidyverse')
    library('pbapply')
    library('caret')
    library('rpart')
  })
  
  # Train the biomarker tree models
  bestmod_mod=
    bestmod_dev %>%
    pblapply(X=.
             ,training_parameters=training_parameters
             ,sample.kind=sample.kind
             ,std_no_outlier=std_no_outlier
             ,cl=cl
             ,function(X,training_parameters,sample.kind,std_no_outlier){
               X %>%
                 lapply(function(x){
                   suppressWarnings(set.seed(33,sample.kind=sample.kind))
                   suppressWarnings(
                     x %>%
                       std_no_outlier(x) %>%
                       `colnames<-`(
                         str_remove_all(colnames(.),'maternal_blood:')
                       ) %>%
                       caret::train(
                         outcome~.
                         ,data=.
                         ,method='rpart'
                         ,weights=
                           x %>%
                           `colnames<-`(
                             str_remove_all(colnames(.),'maternal_blood:')
                           ) %>%
                           pull(outcome) %>%
                           table() %>%
                           as.data.frame() %>%
                           setNames(c('outcome','Freq')) %>%
                           mutate(total=sum(Freq)) %>%
                           mutate(weight=1/(Freq/total)*0.5) %>%
                           right_join(
                             x %>%
                               rownames_to_column(var='id')
                             ,by='outcome'
                           ) %>%
                           pull(weight) %>%
                           setNames(rownames(x))
                         ,metric='ROC'
                         ,trControl=training_parameters$final_trControl
                         ,tuneLength=10
                         ,control=
                           rpart.control(
                             maxdepth=
                               sum(str_detect(colnames(x),'maternal_blood'))
                           )
                         ,parms=list(split='gini')
                       )
                   )
                 })
             })
  
  bestmod_mod %>%
    lapply(X=names(.),Y=.,function(X,Y){
      saveRDS(Y[[X]],paste0('data/bestmod_mod/bestmod_mod_',X,'.rds'))
    })
  
  # Evaluate the biomarker tree models using development dataset
  bestmod_edev=
    bestmod_mod %>%
    pblapply(function(x){
      x %>%
        lapply(function(x){
          x$trainingData %>%
            rename(outcome=.outcome) %>%
            re_evalm(x,th_selection,th_values) %>%
            column_to_rownames(var='metric')
        })
    })
  
  saveRDS(bestmod_edev,'data/bestmod_edev.rds')
}else if(load_emulators){
  cat(readRDS('data/log.rds')[['bestmod_dev']],'\n')
  cat(readRDS('data/log.rds')[['bestmod_mod']],'\n')
  cat(readRDS('data/log.rds')[['bestmod_edev']])
  bestmod_dev=readRDS('data/bestmod_dev.rds')
  bestmod_mod=
    list.files('data/bestmod_mod/',pattern='bestmod_mod_',full.names=T) %>%
    .[order(as.numeric(str_remove_all(.,'data/bestmod_mod/bestmod_mod_|rank|\\.rds')))] %>%
    lapply(readRDS) %>%
    setNames(
      list.files('data/bestmod_mod/',pattern='bestmod_mod_',full.names=T) %>%
        .[order(as.numeric(
          str_remove_all(.,'data/bestmod_mod/bestmod_mod_|rank|\\.rds')
        ))] %>%
        str_remove_all('data/bestmod_mod/bestmod_mod_|\\.rds')
    )
  bestmod_edev=readRDS('data/bestmod_edev.rds')
}
################################################################################

#####{Get features of the replication dataset 3 and evaluate it, include=FALSE}
if(run_heavy_computation){
  bestmod_rep3=
    bestmod_list %>%
    pull(rank) %>%
    .[!duplicated(.)] %>%
    sort() %>%
    pblapply(X=.
             ,Y=GSE149437_phenotype %>%
               select(outcome) %>%
               rownames_to_column(var='id') %>%
               left_join(
                 GSE149437 %>%
                   exprs() %>%
                   mb_predictor(list(maternal_blood=GSE108497_maternal_blood)) %>%
                   rownames_to_column(var='id')
                 ,by='id'
               ) %>%
               column_to_rownames(var='id')
             ,function(X,Y){
               bestmod_list %>%
                 filter(rank==X) %>%
                 pull(grank) %>%
                 .[!duplicated(.)] %>%
                 sort() %>%
                 lapply(X=.,Y=X,Z=Y,function(X,Y,Z){
                   Z %>%
                     select_at(c(
                       'outcome'
                       ,paste0(
                         'maternal_blood:'
                         ,bestmod_list$blood_t[bestmod_list$grank==X & bestmod_list$rank==Y]
                       )
                     ))
                 }) %>%
                 setNames(paste0(
                   'grank'
                   ,bestmod_list %>%
                     filter(rank==X) %>%
                     pull(grank) %>%
                     .[!duplicated(.)] %>%
                     sort()
                 ))
             }) %>%
    setNames(paste0(
      'rank'
      ,bestmod_list %>%
        pull(rank) %>%
        .[!duplicated(.)] %>%
        sort()
    ))
  
  saveRDS(bestmod_rep3,'data/bestmod_rep3.rds')
  
  bestmod_erep3=
    bestmod_mod %>%
    pblapply(X=names(.),Y=.,function(X,Y){
      Y[[X]] %>%
        lapply(X=names(.),Y=.,Z=X,function(X,Y,Z){
          bestmod_rep3[[Z]][[X]] %>%
            std_no_outlier(bestmod_dev[[Z]][[X]]) %>%
            `colnames<-`(
              str_remove_all(colnames(.),'maternal_blood:')
            ) %>%
            re_evalm(Y[[X]],th_selection,th_values) %>%
            column_to_rownames(var='metric')
        }) %>%
        setNames(names(Y[[X]]))
    }) %>%
    setNames(names(bestmod_mod))
  
  saveRDS(bestmod_erep3,'data/bestmod_erep3.rds')
}else if(load_emulators){
  cat(readRDS('data/log.rds')[['bestmod_rep3']],'\n')
  cat(readRDS('data/log.rds')[['bestmod_erep3']])
  bestmod_rep3=readRDS('data/bestmod_rep3.rds')
  bestmod_erep3=readRDS('data/bestmod_erep3.rds')
}
################################################################################

#####{Get features to train logFC2 DEG and evaluate it, include=FALSE}
# Preprocess the features of logFC2 DEG in the development dataset
logFC2_dev=
  GSE108497_phenotype %>%
  select(outcome) %>%
  rownames_to_column(var='id') %>%
  left_join(
    GSE108497 %>%
      exprs() %>%
      mb_predictor(list(maternal_blood=GSE108497_maternal_blood)) %>%
      rownames_to_column(var='id')
    ,by='id'
  ) %>%
  column_to_rownames(var='id') %>%
  select_at(c(
    'outcome'
    ,paste0(
      'maternal_blood:'
      ,GSE108497_maternal_blood$result %>%
        filter(adj.P.Val<0.05 & abs(logFC)>2) %>%
        rownames()
    )
  ))

# Preprocess the features of logFC2 DEG in the replication dataset 3
logFC2_rep3=
  GSE149437_phenotype %>%
  select(outcome) %>%
  rownames_to_column(var='id') %>%
  left_join(
    GSE149437 %>%
      exprs() %>%
      mb_predictor(list(maternal_blood=GSE108497_maternal_blood)) %>%
      rownames_to_column(var='id')
    ,by='id'
  ) %>%
  column_to_rownames(var='id') %>%
  select_at(c(
    'outcome'
    ,paste0(
      'maternal_blood:'
      ,GSE108497_maternal_blood$result %>%
        filter(adj.P.Val<0.05 & abs(logFC)>2) %>%
        rownames()
    )
  ))

# Train the logFC2-DEG tree model
logFC2_mod=
  logFC2_dev %>%
  lapply(X=1,Y=.,function(X,Y){
    suppressWarnings(set.seed(33,sample.kind=sample.kind))
    suppressWarnings(
      Y %>%
        std_no_outlier(Y) %>%
        `colnames<-`(
          str_remove_all(colnames(.),'maternal_blood:')
        ) %>%
        caret::train(
          outcome~.
          ,data=.
          ,method='rpart'
          ,weights=
            Y %>%
            `colnames<-`(
              str_remove_all(colnames(.),'maternal_blood:')
            ) %>%
            pull(outcome) %>%
            table() %>%
            as.data.frame() %>%
            setNames(c('outcome','Freq')) %>%
            mutate(total=sum(Freq)) %>%
            mutate(weight=1/(Freq/total)*0.5) %>%
            right_join(
              Y %>%
                rownames_to_column(var='id')
              ,by='outcome'
            ) %>%
            pull(weight) %>%
            setNames(rownames(Y))
          ,metric='ROC'
          ,trControl=training_parameters$final_trControl
          ,tuneLength=10 
          ,control=
            rpart.control(
              maxdepth=
                sum(str_detect(colnames(Y),'maternal_blood'))
            )
          ,parms=list(split='gini')
        )
    )
  }) %>%
  .[[1]]

# Evaluate the logFC2-DEG tree model using development dataset
logFC2_edev=
  logFC2_dev %>%
  std_no_outlier(logFC2_dev) %>%
  `colnames<-`(
    str_remove_all(colnames(.),'maternal_blood:')
  ) %>%
  re_evalm(logFC2_mod,th_selection,th_values) %>%
  column_to_rownames(var='metric')

# Evaluate the logFC2-DEG tree model using replication dataset 3
logFC2_erep3=
  logFC2_rep3 %>%
  std_no_outlier(logFC2_dev) %>%
  `colnames<-`(
    str_remove_all(colnames(.),'maternal_blood:')
  ) %>%
  re_evalm(logFC2_mod,th_selection,th_values) %>%
  column_to_rownames(var='metric')
################################################################################

#####{Get features to train devrep3DEG and evaluate it, include=FALSE}
# Get the possible combinations of devrep3DEG
devrep3DEG_list=
  GSE108497_maternal_blood$result %>%
  filter(adj.P.Val<0.05) %>%
  rownames() %>%
  intersect(
    GSE149437_maternal_blood$result %>%
      filter(adj.P.Val<0.05) %>%
      rownames()
  ) %>%
  lapply(X=seq(2),Y=.,function(X,Y){
    if(X==1){
      as.list(Y)
    }else{
      rep(list(Y),X) %>%
        expand.grid() %>%
        lapply(X=seq(nrow(.)),Y=.,Z=X,function(X,Y,Z){
          Y=Y %>%
            mutate_all(as.character)
          K=c()
          for(i in seq(Z)){
            K[i]=Y[[i]][X]
          }
          K %>%
            sort() %>%
            data.frame(value=.) %>%
            mutate(key=paste0('Var',seq(Z))) %>%
            spread(key,value)
        }) %>%
        do.call(rbind,.) %>%
        filter(!duplicated(.)) %>%
        mutate(seq=seq(nrow(.))) %>%
        gather(key,value,-seq) %>%
        arrange(seq,key) %>%
        group_by(seq) %>%
        filter(!duplicated(value)) %>%
        filter(n()==X) %>%
        lapply(X=.$seq %>% .[!duplicated(.)],Y=.,function(X,Y){
          Y %>%
            filter(seq==X) %>%
            pull(value)
        })
    }
  }) %>%
  lapply(X=1,Y=.,function(X,Y){
    for(i in seq(length(Y))){
      if(i==1){
        Z=Y[[i]]
      }else{
        Z=Z %>%
          c(Y[[i]])
      }
    }
    Z
  }) %>%
  .[[1]]

if(run_heavy_computation){
  # Preprocess the features of devrep3DEG in the development dataset
  devrep3DEG_dev=
    GSE108497_phenotype %>%
    select(outcome) %>%
    rownames_to_column(var='id') %>%
    left_join(
      GSE108497 %>%
        exprs() %>%
        mb_predictor(list(maternal_blood=GSE108497_maternal_blood)) %>%
        rownames_to_column(var='id')
      ,by='id'
    ) %>%
    column_to_rownames(var='id') %>%
    pblapply(X=seq(length(devrep3DEG_list)),Y=devrep3DEG_list,Z=.,function(X,Y,Z){
      Z %>%
        select_at(c(
          'outcome'
          ,paste0(
            'maternal_blood:'
            ,Y[[X]]
          )
        ))
    })
  
  saveRDS(devrep3DEG_dev,'data/devrep3DEG_dev.rds')
  
  # Preprocess the features of devrep3DEG in the replication dataset 3
  devrep3DEG_rep3=
    GSE149437_phenotype %>%
    select(outcome) %>%
    rownames_to_column(var='id') %>%
    left_join(
      GSE149437 %>%
        exprs() %>%
        mb_predictor(list(maternal_blood=GSE108497_maternal_blood)) %>%
        rownames_to_column(var='id')
      ,by='id'
    ) %>%
    column_to_rownames(var='id') %>%
    pblapply(X=seq(length(devrep3DEG_list)),Y=devrep3DEG_list,Z=.,function(X,Y,Z){
      Z %>%
        select_at(c(
          'outcome'
          ,paste0(
            'maternal_blood:'
            ,Y[[X]]
          )
        ))
    })
  
  saveRDS(devrep3DEG_rep3,'data/devrep3DEG_rep3.rds')
  
  # Create multiple CPU clusters for parallel computing.
  cl=detectCores()-2
  cl=makePSOCKcluster(cl)
  registerDoParallel(cl)
  
  # Close the CPU clusters and clean up memory at exit.
  on.exit(stopCluster(cl))
  on.exit(registerDoSEQ())
  on.exit(rm(cl))
  on.exit(rm(the_method))
  on.exit(gc())
  
  # Bring these packages to the CPU clusters.
  clusterEvalQ(cl,{
    library('parallel')
    library('doParallel')
    library('tidyverse')
    library('pbapply')
    library('caret')
    library('rpart')
  })
  
  # Train the devrep3DEG tree model
  devrep3DEG_mod=
    devrep3DEG_dev %>%
    pblapply(X=.
             ,training_parameters=training_parameters
             ,sample.kind=sample.kind
             ,std_no_outlier=std_no_outlier
             ,cl=cl
             ,function(X,training_parameters,sample.kind,std_no_outlier){
               suppressWarnings(set.seed(33,sample.kind=sample.kind))
               suppressWarnings(
                 X %>%
                   std_no_outlier(X) %>%
                   `colnames<-`(
                     str_remove_all(colnames(.),'maternal_blood:')
                   ) %>%
                   caret::train(
                     outcome~.
                     ,data=.
                     ,method='rpart'
                     ,weights=
                       X %>%
                       `colnames<-`(
                         str_remove_all(colnames(.),'maternal_blood:')
                       ) %>%
                       pull(outcome) %>%
                       table() %>%
                       as.data.frame() %>%
                       setNames(c('outcome','Freq')) %>%
                       mutate(total=sum(Freq)) %>%
                       mutate(weight=1/(Freq/total)*0.5) %>%
                       right_join(
                         X %>%
                           rownames_to_column(var='id')
                         ,by='outcome'
                       ) %>%
                       pull(weight) %>%
                       setNames(rownames(X))
                     ,metric='ROC'
                     ,trControl=training_parameters$final_trControl
                     ,tuneLength=10 
                     ,control=
                       rpart.control(
                         maxdepth=
                           sum(str_detect(colnames(X),'maternal_blood'))
                       )
                     ,parms=list(split='gini')
                   )
               )
             })
  
  devrep3DEG_mod %>%
    lapply(X=seq(length(.)),Y=.,function(X,Y){
      saveRDS(Y[[X]],paste0('data/devrep3DEG_mod/devrep3DEG_mod_',X,'.rds'))
    })
  
  # Evaluate the devrep3DEG tree model using development dataset
  devrep3DEG_edev=
    devrep3DEG_mod %>%
    pblapply(X=seq(length(.)),Y=.,function(X,Y){
      devrep3DEG_dev[[X]] %>%
        std_no_outlier(devrep3DEG_dev[[X]]) %>%
        `colnames<-`(
          str_remove_all(colnames(.),'maternal_blood:')
        ) %>%
        re_evalm(Y[[X]],th_selection,th_values) %>%
        column_to_rownames(var='metric')
    })
  
  saveRDS(devrep3DEG_edev,'data/devrep3DEG_edev.rds')
  
  # Evaluate the devrep3DEG tree model using replication dataset 3
  devrep3DEG_erep3=
    devrep3DEG_mod %>%
    pblapply(X=seq(length(.)),Y=.,function(X,Y){
      devrep3DEG_rep3[[X]] %>%
        std_no_outlier(devrep3DEG_dev[[X]]) %>%
        `colnames<-`(
          str_remove_all(colnames(.),'maternal_blood:')
        ) %>%
        re_evalm(Y[[X]],th_selection,th_values) %>%
        column_to_rownames(var='metric')
    })
  
  saveRDS(devrep3DEG_erep3,'data/devrep3DEG_erep3.rds')
}else if(load_emulators){
  cat(readRDS('data/log.rds')[['devrep3DEG_dev']],'\n')
  cat(readRDS('data/log.rds')[['devrep3DEG_rep3']],'\n')
  cat(readRDS('data/log.rds')[['devrep3DEG_mod']],'\n')
  cat(readRDS('data/log.rds')[['devrep3DEG_edev']],'\n')
  cat(readRDS('data/log.rds')[['devrep3DEG_erep3']])
  devrep3DEG_dev=readRDS('data/devrep3DEG_dev.rds')
  devrep3DEG_rep3=readRDS('data/devrep3DEG_rep3.rds')
  devrep3DEG_mod=
    list.files('data/devrep3DEG_mod/',pattern='devrep3DEG_mod_',full.names=T) %>%
    .[order(as.numeric(
      str_remove_all(.,'data/devrep3DEG_mod/devrep3DEG_mod_|rank|\\.rds')
    ))] %>%
    lapply(readRDS)
  devrep3DEG_edev=readRDS('data/devrep3DEG_edev.rds')
  devrep3DEG_erep3=readRDS('data/devrep3DEG_erep3.rds')
}
################################################################################

#####{Get features to train devnorep3DEG and evaluate it, include=FALSE}
# Get the possible combinations of devnorep3DEG
devnorep3DEG_list=
  GSE108497_maternal_blood$result %>%
  filter(adj.P.Val<0.05) %>%
  rownames() %>%
  setdiff(
    GSE149437_maternal_blood$result %>%
      filter(adj.P.Val<0.05) %>%
      rownames()
  ) %>%
  lapply(X=seq(1),Y=.,function(X,Y){
    if(X==1){
      as.list(Y)
    }else{
      Z=rep(list(Y),X) %>%
        expand.grid() %>%
        pblapply(X=seq(nrow(.)),Y=.,Z=X,function(X,Y,Z){
          Y=Y %>%
            mutate_all(as.character)
          K=c()
          for(i in seq(Z)){
            K[i]=Y[[i]][X]
          }
          K %>%
            sort() %>%
            data.frame(value=.) %>%
            mutate(key=paste0('Var',seq(Z))) %>%
            spread(key,value)
        }) %>%
        do.call(rbind,.) %>%
        filter(!duplicated(.)) %>%
        mutate(seq=seq(nrow(.))) %>%
        gather(key,value,-seq) %>%
        arrange(seq,key) %>%
        group_by(seq) %>%
        filter(!duplicated(value)) %>%
        filter(n()==X) %>%
        lapply(X=.$seq %>% .[!duplicated(.)],Y=.,function(X,Y){
          Y %>%
            filter(seq==X) %>%
            pull(value)
        })
    }
  }) %>%
  lapply(X=1,Y=.,function(X,Y){
    for(i in seq(length(Y))){
      if(i==1){
        Z=Y[[i]]
      }else{
        Z=Z %>%
          c(Y[[i]])
      }
    }
    Z
  }) %>%
  .[[1]]

if(run_heavy_computation){
  # Preprocess the features of devnorep3DEG in the development dataset
  devnorep3DEG_dev=
    GSE108497_phenotype %>%
    select(outcome) %>%
    rownames_to_column(var='id') %>%
    left_join(
      GSE108497 %>%
        exprs() %>%
        mb_predictor(list(maternal_blood=GSE108497_maternal_blood)) %>%
        # mutate_all(function(x)ifelse(x>=0,1,0)) %>%
        rownames_to_column(var='id')
      ,by='id'
    ) %>%
    column_to_rownames(var='id') %>%
    pblapply(X=seq(length(devnorep3DEG_list)),Y=devnorep3DEG_list,Z=.,function(X,Y,Z){
      Z %>%
        select_at(c(
          'outcome'
          ,paste0(
            'maternal_blood:'
            ,Y[[X]]
          )
        ))
    })
  
  saveRDS(devnorep3DEG_dev,'data/devnorep3DEG_dev.rds')
  
  # Preprocess the features of devnorep3DEG in the replication dataset 3
  devnorep3DEG_rep3=
    GSE149437_phenotype %>%
    select(outcome) %>%
    rownames_to_column(var='id') %>%
    left_join(
      GSE149437 %>%
        exprs() %>%
        mb_predictor(list(maternal_blood=GSE108497_maternal_blood)) %>%
        # mutate_all(function(x)ifelse(x>=0,1,0)) %>%
        rownames_to_column(var='id')
      ,by='id'
    ) %>%
    column_to_rownames(var='id') %>%
    pblapply(X=seq(length(devnorep3DEG_list)),Y=devnorep3DEG_list,Z=.,function(X,Y,Z){
      Z %>%
        select_at(c(
          'outcome'
          ,paste0(
            'maternal_blood:'
            ,Y[[X]]
          )
        ))
    })
  
  saveRDS(devnorep3DEG_rep3,'data/devnorep3DEG_rep3.rds')
  
  # Create multiple CPU clusters for parallel computing.
  cl=detectCores()-2
  cl=makePSOCKcluster(cl)
  registerDoParallel(cl)
  
  # Close the CPU clusters and clean up memory at exit.
  on.exit(stopCluster(cl))
  on.exit(registerDoSEQ())
  on.exit(rm(cl))
  on.exit(rm(the_method))
  on.exit(gc())
  
  # Bring these packages to the CPU clusters.
  clusterEvalQ(cl,{
    library('parallel')
    library('doParallel')
    library('tidyverse')
    library('pbapply')
    library('caret')
    library('rpart')
  })
  
  # Train the devnorep3DEG tree model
  devnorep3DEG_mod=
    devnorep3DEG_dev %>%
    pblapply(X=.
             ,training_parameters=training_parameters
             ,sample.kind=sample.kind
             ,std_no_outlier=std_no_outlier
             ,cl=cl
             ,function(X,training_parameters,sample.kind,std_no_outlier){
               suppressWarnings(set.seed(33,sample.kind=sample.kind))
               suppressWarnings(
                 X %>%
                   std_no_outlier(X) %>%
                   `colnames<-`(
                     str_remove_all(colnames(.),'maternal_blood:')
                   ) %>%
                   caret::train(
                     outcome~.
                     ,data=.
                     ,method='rpart'
                     ,weights=
                       X %>%
                       `colnames<-`(
                         str_remove_all(colnames(.),'maternal_blood:')
                       ) %>%
                       pull(outcome) %>%
                       table() %>%
                       as.data.frame() %>%
                       setNames(c('outcome','Freq')) %>%
                       mutate(total=sum(Freq)) %>%
                       mutate(weight=1/(Freq/total)*0.5) %>%
                       right_join(
                         X %>%
                           rownames_to_column(var='id')
                         ,by='outcome'
                       ) %>%
                       pull(weight) %>%
                       setNames(rownames(X))
                     ,metric='ROC'
                     ,trControl=training_parameters$final_trControl
                     ,tuneLength=10 
                     ,control=
                       rpart.control(
                         maxdepth=
                           sum(str_detect(colnames(X),'maternal_blood'))
                       )
                     ,parms=list(split='gini')
                   )
               )
             })
  
  devnorep3DEG_mod %>%
    lapply(X=seq(length(.)),Y=.,function(X,Y){
      saveRDS(Y[[X]],paste0('data/devnorep3DEG_mod/devnorep3DEG_mod_',X,'.rds'))
    })
  
  # Evaluate the devnorep3DEG tree model using development dataset
  devnorep3DEG_edev=
    devnorep3DEG_mod %>%
    pblapply(X=seq(length(.)),Y=.,function(X,Y){
      devnorep3DEG_dev[[X]] %>%
        std_no_outlier(devnorep3DEG_dev[[X]]) %>%
        `colnames<-`(
          str_remove_all(colnames(.),'maternal_blood:')
        ) %>%
        re_evalm(Y[[X]],th_selection,th_values) %>%
        column_to_rownames(var='metric')
    })
  
  saveRDS(devnorep3DEG_edev,'data/devnorep3DEG_edev.rds')
  
  # Evaluate the devnorep3DEG tree model using replication dataset 3
  devnorep3DEG_erep3=
    devnorep3DEG_mod %>%
    pblapply(X=seq(length(.)),Y=.,function(X,Y){
      devnorep3DEG_rep3[[X]] %>%
        std_no_outlier(devnorep3DEG_dev[[X]]) %>%
        `colnames<-`(
          str_remove_all(colnames(.),'maternal_blood:')
        ) %>%
        re_evalm(Y[[X]],th_selection,th_values) %>%
        column_to_rownames(var='metric')
    })
  
  saveRDS(devnorep3DEG_erep3,'data/devnorep3DEG_erep3.rds')
}else if(load_emulators){
  cat(readRDS('data/log.rds')[['devnorep3DEG_dev']],'\n')
  cat(readRDS('data/log.rds')[['devnorep3DEG_rep3']],'\n')
  cat(readRDS('data/log.rds')[['devnorep3DEG_mod']],'\n')
  cat(readRDS('data/log.rds')[['devnorep3DEG_edev']],'\n')
  cat(readRDS('data/log.rds')[['devnorep3DEG_erep3']])
  devnorep3DEG_dev=readRDS('data/devnorep3DEG_dev.rds')
  devnorep3DEG_rep3=readRDS('data/devnorep3DEG_rep3.rds')
  devnorep3DEG_mod=
    list.files('data/devnorep3DEG_mod/',pattern='devnorep3DEG_mod_',full.names=T) %>%
    .[order(as.numeric(
      str_remove_all(.,'data/devnorep3DEG_mod/devnorep3DEG_mod_|rank|\\.rds')
    ))] %>%
    lapply(readRDS)
  devnorep3DEG_edev=readRDS('data/devnorep3DEG_edev.rds')
  devnorep3DEG_erep3=readRDS('data/devnorep3DEG_erep3.rds')
}
################################################################################

#####{Get features to train recentDEG and evaluate it, include=FALSE}
# Get the possible combinations of recentDEG
recentDEG_list=
  c('CAMK2G'
    ,'DERA'
    ,'FAM46A'
    ,'KIAA1109'
    ,'LRRC58'
    ,'MYLIP'
    ,'NDUFV3'
    ,'NMRK1'
    ,'PI4KA'
    ,'PRTFDC1'
    ,'PYGO2'
    ,'RNF149'
    ,'TFIP11'
    ,'TRIM21'
    ,'USB1'
    ,'YWHAQP5'
    ,'Y_RNA'
  ) %>%
  # "FAM46A"  "PI4KA"   "PRTFDC1" "USB1"    "YWHAQP5" "Y_RNA"
  intersect(rownames(GSE108497_maternal_blood$result)) %>%
  # "MYLIP"
  intersect(rownames(GSE177477)) %>%
  lapply(X=seq(2),Y=.,function(X,Y){
    if(X==1){
      as.list(Y)
    }else{
      rep(list(Y),X) %>%
        expand.grid() %>%
        lapply(X=seq(nrow(.)),Y=.,Z=X,function(X,Y,Z){
          Y=Y %>%
            mutate_all(as.character)
          K=c()
          for(i in seq(Z)){
            K[i]=Y[[i]][X]
          }
          K %>%
            sort() %>%
            data.frame(value=.) %>%
            mutate(key=paste0('Var',seq(Z))) %>%
            spread(key,value)
        }) %>%
        do.call(rbind,.) %>%
        filter(!duplicated(.)) %>%
        mutate(seq=seq(nrow(.))) %>%
        gather(key,value,-seq) %>%
        arrange(seq,key) %>%
        group_by(seq) %>%
        filter(!duplicated(value)) %>%
        filter(n()==X) %>%
        lapply(X=.$seq %>% .[!duplicated(.)],Y=.,function(X,Y){
          Y %>%
            filter(seq==X) %>%
            pull(value)
        })
    }
  }) %>%
  lapply(X=1,Y=.,function(X,Y){
    for(i in seq(length(Y))){
      if(i==1){
        Z=Y[[i]]
      }else{
        Z=Z %>%
          c(Y[[i]])
      }
    }
    Z
  }) %>%
  .[[1]]

if(run_heavy_computation){
  # Preprocess the features of recentDEG in the development dataset
  recentDEG_dev=
    GSE108497_phenotype %>%
    select(outcome) %>%
    rownames_to_column(var='id') %>%
    left_join(
      GSE108497 %>%
        exprs() %>%
        mb_predictor(list(maternal_blood=GSE108497_maternal_blood)) %>%
        # mutate_all(function(x)ifelse(x>=0,1,0)) %>%
        rownames_to_column(var='id')
      ,by='id'
    ) %>%
    column_to_rownames(var='id') %>%
    pblapply(X=seq(length(recentDEG_list)),Y=recentDEG_list,Z=.,function(X,Y,Z){
      Z %>%
        select_at(c(
          'outcome'
          ,paste0(
            'maternal_blood:'
            ,Y[[X]]
          )
        ))
    })
  
  saveRDS(recentDEG_dev,'data/recentDEG_dev.rds')
  
  # Preprocess the features of recentDEG in the replication dataset 3
  recentDEG_rep3=
    GSE149437_phenotype %>%
    select(outcome) %>%
    rownames_to_column(var='id') %>%
    left_join(
      GSE149437 %>%
        exprs() %>%
        mb_predictor(list(maternal_blood=GSE108497_maternal_blood)) %>%
        # mutate_all(function(x)ifelse(x>=0,1,0)) %>%
        rownames_to_column(var='id')
      ,by='id'
    ) %>%
    column_to_rownames(var='id') %>%
    pblapply(X=seq(length(recentDEG_list)),Y=recentDEG_list,Z=.,function(X,Y,Z){
      Z %>%
        select_at(c(
          'outcome'
          ,paste0(
            'maternal_blood:'
            ,Y[[X]]
          )
        ))
    })
  
  saveRDS(recentDEG_rep3,'data/recentDEG_rep3.rds')
  
  # Create multiple CPU clusters for parallel computing.
  cl=detectCores()-2
  cl=makePSOCKcluster(cl)
  registerDoParallel(cl)
  
  # Close the CPU clusters and clean up memory at exit.
  on.exit(stopCluster(cl))
  on.exit(registerDoSEQ())
  on.exit(rm(cl))
  on.exit(rm(the_method))
  on.exit(gc())
  
  # Bring these packages to the CPU clusters.
  clusterEvalQ(cl,{
    library('parallel')
    library('doParallel')
    library('tidyverse')
    library('pbapply')
    library('caret')
    library('rpart')
  })
  
  # Train the recentDEG tree model
  recentDEG_mod=
    recentDEG_dev %>%
    pblapply(X=.
             ,training_parameters=training_parameters
             ,sample.kind=sample.kind
             ,std_no_outlier=std_no_outlier
             ,cl=cl
             ,function(X,training_parameters,sample.kind,std_no_outlier){
               suppressWarnings(set.seed(33,sample.kind=sample.kind))
               suppressWarnings(
                 X %>%
                   std_no_outlier(X) %>%
                   `colnames<-`(
                     str_remove_all(colnames(.),'maternal_blood:')
                   ) %>%
                   caret::train(
                     outcome~.
                     ,data=.
                     ,method='rpart'
                     ,weights=
                       X %>%
                       `colnames<-`(
                         str_remove_all(colnames(.),'maternal_blood:')
                       ) %>%
                       pull(outcome) %>%
                       table() %>%
                       as.data.frame() %>%
                       setNames(c('outcome','Freq')) %>%
                       mutate(total=sum(Freq)) %>%
                       mutate(weight=1/(Freq/total)*0.5) %>%
                       right_join(
                         X %>%
                           rownames_to_column(var='id')
                         ,by='outcome'
                       ) %>%
                       pull(weight) %>%
                       setNames(rownames(X))
                     ,metric='ROC'
                     ,trControl=training_parameters$final_trControl
                     ,tuneLength=10 
                     ,control=
                       rpart.control(
                         maxdepth=
                           sum(str_detect(colnames(X),'maternal_blood'))
                       )
                     ,parms=list(split='gini')
                   )
               )
             })
  
  recentDEG_mod %>%
    lapply(X=seq(length(.)),Y=.,function(X,Y){
      saveRDS(Y[[X]],paste0('data/recentDEG_mod/recentDEG_mod_',X,'.rds'))
    })
  
  # Evaluate the recentDEG tree model using development dataset
  recentDEG_edev=
    recentDEG_mod %>%
    pblapply(X=seq(length(.)),Y=.,function(X,Y){
      recentDEG_dev[[X]] %>%
        std_no_outlier(recentDEG_dev[[X]]) %>%
        `colnames<-`(
          str_remove_all(colnames(.),'maternal_blood:')
        ) %>%
        re_evalm(Y[[X]],th_selection,th_values) %>%
        column_to_rownames(var='metric')
    })
  
  saveRDS(recentDEG_edev,'data/recentDEG_edev.rds')
  
  # Evaluate the recentDEG tree model using replication dataset 3
  recentDEG_erep3=
    recentDEG_mod %>%
    pblapply(X=seq(length(.)),Y=.,function(X,Y){
      recentDEG_rep3[[X]] %>%
        std_no_outlier(recentDEG_dev[[X]]) %>%
        `colnames<-`(
          str_remove_all(colnames(.),'maternal_blood:')
        ) %>%
        re_evalm(Y[[X]],th_selection,th_values) %>%
        column_to_rownames(var='metric')
    })
  
  saveRDS(recentDEG_erep3,'data/recentDEG_erep3.rds')
}else if(load_emulators){
  cat(readRDS('data/log.rds')[['recentDEG_dev']],'\n')
  cat(readRDS('data/log.rds')[['recentDEG_rep3']],'\n')
  cat(readRDS('data/log.rds')[['recentDEG_mod']],'\n')
  cat(readRDS('data/log.rds')[['recentDEG_edev']],'\n')
  cat(readRDS('data/log.rds')[['recentDEG_erep3']])
  recentDEG_dev=readRDS('data/recentDEG_dev.rds')
  recentDEG_rep3=readRDS('data/recentDEG_rep3.rds')
  recentDEG_mod=
    list.files('data/recentDEG_mod/',pattern='recentDEG_mod_',full.names=T) %>%
    .[order(as.numeric(
      str_remove_all(.,'data/recentDEG_mod/recentDEG_mod_|rank|\\.rds')
    ))] %>%
    lapply(readRDS)
  recentDEG_edev=readRDS('data/recentDEG_edev.rds')
  recentDEG_erep3=readRDS('data/recentDEG_erep3.rds')
}
################################################################################

#####{Create a function to compute AUROCs given a list, include=FALSE}
source('R/get_auroc-function.R')
################################################################################

#####{Compute AUROCs for all the methods of predictor discovery, include=FALSE}
if(run_heavy_computation){
  emul_aurocs=list()
  
  emul_aurocs$bestmod_mod$train=
    bestmod_dev %>%
    get_auroc(bestmod_mod)
  
  emul_aurocs$bestmod_mod$test=
    bestmod_rep3 %>%
    get_auroc(bestmod_mod)
  
  emul_aurocs$logFC2_mod$train=
    list(ext=list(reme=logFC2_dev)) %>%
    get_auroc(list(ext=list(reme=logFC2_mod)))
  
  emul_aurocs$logFC2_mod$test=
    list(ext=list(reme=logFC2_rep3)) %>%
    get_auroc(list(ext=list(reme=logFC2_mod)))
  
  emul_aurocs$devrep3DEG_mod$train=
    devrep3DEG_dev %>%
    lapply(X=seq(length(.)),Y=.,function(X,Y){
      list(both=Y[[X]])
    }) %>%
    setNames(paste0('mod',str_pad(seq(length(.)),3,'left','0'))) %>%
    get_auroc(
      devrep3DEG_mod %>%
        lapply(X=seq(length(.)),Y=.,function(X,Y){
          list(both=Y[[X]])
        }) %>%
        setNames(paste0('mod',str_pad(seq(length(.)),3,'left','0')))
    )
  
  emul_aurocs$devrep3DEG_mod$test=
    devrep3DEG_rep3 %>%
    lapply(X=seq(length(.)),Y=.,function(X,Y){
      list(both=Y[[X]])
    }) %>%
    setNames(paste0('mod',str_pad(seq(length(.)),3,'left','0'))) %>%
    get_auroc(
      devrep3DEG_mod %>%
        lapply(X=seq(length(.)),Y=.,function(X,Y){
          list(both=Y[[X]])
        }) %>%
        setNames(paste0('mod',str_pad(seq(length(.)),3,'left','0')))
    )
  
  emul_aurocs$devnorep3DEG_mod$train=
    devnorep3DEG_dev %>%
    lapply(X=seq(length(.)),Y=.,function(X,Y){
      list(not_in=Y[[X]])
    }) %>%
    setNames(paste0('mod',str_pad(seq(length(.)),3,'left','0'))) %>%
    get_auroc(
      devnorep3DEG_mod %>%
        lapply(X=seq(length(.)),Y=.,function(X,Y){
          list(not_in=Y[[X]])
        }) %>%
        setNames(paste0('mod',str_pad(seq(length(.)),3,'left','0')))
    )
  
  emul_aurocs$devnorep3DEG_mod$test=
    devnorep3DEG_rep3 %>%
    lapply(X=seq(length(.)),Y=.,function(X,Y){
      list(not_in=Y[[X]])
    }) %>%
    setNames(paste0('mod',str_pad(seq(length(.)),3,'left','0'))) %>%
    get_auroc(
      devnorep3DEG_mod %>%
        lapply(X=seq(length(.)),Y=.,function(X,Y){
          list(not_in=Y[[X]])
        }) %>%
        setNames(paste0('mod',str_pad(seq(length(.)),3,'left','0')))
    )
  
  emul_aurocs$recentDEG_mod$train=
    recentDEG_dev %>%
    lapply(X=seq(length(.)),Y=.,function(X,Y){
      list(recent=Y[[X]])
    }) %>%
    setNames(paste0('mod',str_pad(seq(length(.)),3,'left','0'))) %>%
    get_auroc(
      recentDEG_mod %>%
        lapply(X=seq(length(.)),Y=.,function(X,Y){
          list(recent=Y[[X]])
        }) %>%
        setNames(paste0('mod',str_pad(seq(length(.)),3,'left','0')))
    )
  
  emul_aurocs$recentDEG_mod$test=
    recentDEG_rep3 %>%
    lapply(X=seq(length(.)),Y=.,function(X,Y){
      list(recent=Y[[X]])
    }) %>%
    setNames(paste0('mod',str_pad(seq(length(.)),3,'left','0'))) %>%
    get_auroc(
      recentDEG_mod %>%
        lapply(X=seq(length(.)),Y=.,function(X,Y){
          list(recent=Y[[X]])
        }) %>%
        setNames(paste0('mod',str_pad(seq(length(.)),3,'left','0')))
    )
  
  saveRDS(emul_aurocs,'data/emul_aurocs.rds')
}else{
  cat(readRDS('data/log.rds')[['emul_aurocs']])
  emul_aurocs=readRDS('data/emul_aurocs.rds')
}
################################################################################

#####{Summarize the AUROC comparison for the predictor discovery, include=FALSE}
if(run_heavy_computation){
  sum_aurocs=
    emul_aurocs %>%
    lapply(X=names(.),Y=.,function(X,Y){
      Y[[X]] %>%
        lapply(X=names(.),Y=.,function(X,Y){
          Y[[X]] %>%
            mutate(set=X)
        }) %>%
        do.call(rbind,.) %>%
        mutate(mod_list=X)
    }) %>%
    do.call(rbind,.) %>%
    pivot_wider(names_from=c('set'),values_from=c('Score','LB','UB')) %>%
    lapply(X=1,Y=.,function(X,Y){
      
      Z=GSE177477_phenotype %>%
        select(outcome) %>%
        rownames_to_column(var='id') %>%
        left_join(
          GSE177477 %>%
            exprs() %>%
            mb_predictor(list(maternal_blood=GSE108497_maternal_blood)) %>%
            rownames_to_column(var='id')
          ,by='id'
        ) %>%
        left_join(
          GSE177477 %>%
            phenoData() %>%
            pData() %>%
            rownames_to_column(var='gsm') %>%
            select(c('gsm',colnames(.) %>% .[str_detect(.,':ch1')])) %>%
            setNames(
              colnames(.) %>%
                str_remove_all(':ch1') %>%
                str_replace_all('\\s','_')
            ) %>%
            select(-id) %>%
            rename(id=gsm) %>%
            select(id,symptons,severity)
          ,by='id'
        ) %>%
        mutate(
          outcome=
            factor(ifelse(
              !symptons%in%c('Healthy controls')
              ,'event','nonevent'
            ),levels(outcome))
        ) %>%
        select(-symptons,-severity) %>%
        `colnames<-`(str_remove_all(colnames(.),'maternal_blood:')) %>%
        column_to_rownames(var='id')
      
      K=list(
        bestmod_mod=
          bestmod_mod
        ,logFC2_mod=
          list(ext=list(reme=logFC2_mod))
        ,devrep3DEG_mod=
          devrep3DEG_mod %>%
          lapply(X=seq(length(.)),Y=.,function(X,Y){
            list(both=Y[[X]])
          }) %>%
          setNames(paste0('mod',str_pad(seq(length(.)),3,'left','0')))
        ,devnorep3DEG_mod=
          devnorep3DEG_mod %>%
          lapply(X=seq(length(.)),Y=.,function(X,Y){
            list(not_in=Y[[X]])
          }) %>%
          setNames(paste0('mod',str_pad(seq(length(.)),3,'left','0')))
        ,recentDEG_mod=
          recentDEG_mod %>%
          lapply(X=seq(length(.)),Y=.,function(X,Y){
            list(recent=Y[[X]])
          }) %>%
          setNames(paste0('mod',str_pad(seq(length(.)),3,'left','0')))
      )
      
      Y=Y %>%
        mutate(
          biom_avail=
            cand_predictor %>%
            sapply(function(x){
              y=x %>%
                str_split(',') %>%
                .[[1]]
              all(y%in%colnames(Z))
            })
        )
      
      pblapply(X=seq(nrow(Y)),Y=Y,Z=Z,K=K,function(X,Y,Z,K){
        L=Y[X,,drop=F] %>%
          separate(id,c('level1','level2'),sep='\\.')
        if(L$biom_avail){
          M=K %>%
            .[[L$mod_list]] %>%
            .[[L$level1]] %>%
            .[[L$level2]] %>%
            predict(newdata=Z,type='prob') %>%
            cbind(select(Z,outcome)) %>%
            select(nonevent,event,outcome) %>%
            evalm(silent=T,showplots=F) %>%
            .$optres %>%
            .$Group1 %>%
            .['AUC-ROC',] %>%
            `rownames<-`(NULL) %>%
            separate(CI,c('LB','UB'),sep='-') %>%
            mutate_at(c('LB','UB'),as.numeric) %>%
            rename_all(function(x)paste0(x,'_covid'))
        }else{
          M=data.frame(
            Score_covid=NA
            ,LB_covid=NA
            ,UB_covid=NA
          )
        }
        cbind(L,M)
      }) %>%
        do.call(rbind,.)
      
    }) %>%
    .[[1]]
  
  saveRDS(sum_aurocs,'data/sum_aurocs.rds')
}else{
  cat(readRDS('data/log.rds')[['sum_aurocs']])
  sum_aurocs=readRDS('data/sum_aurocs.rds')
}
################################################################################

#####{Assess the eligibility of the emulated biomarkers, include=FALSE}
eligible_biomarkers=
  sum_aurocs %>%
  mutate(
    pred_avail=predictors!=''
    ,replicable=Score_test>=LB_train & Score_test<=UB_train
    ,covidable=
      case_when(
        Score_covid>=UB_train ~ 'upper'
        ,Score_covid>=LB_train & Score_covid<=UB_train ~ 'on_train'
        ,Score_covid<=LB_train ~ 'lower'
        ,TRUE ~ ''
      )
    ,covid_signif=LB_covid>=0.5
    ,quality=
      ifelse(
        biom_avail
        & pred_avail
        & replicable
        & covidable=='lower'
        & !is.na(covid_signif)
        & !covid_signif
        ,'high','low'
      ) %>%
      factor(c('low','high'))
  )
################################################################################






################################################################################
################################################################################
################################################################################
# Justification of the biological relevance
################################################################################
################################################################################
################################################################################

################################################################################
## Best biomarkers
################################################################################

#####{Get the best tree from the best method of discovery, include=FALSE}
if(load_emulators){
  best_each_method=
    eligible_biomarkers %>%
    filter(quality=='high' & mod_list=='bestmod_mod') %>%
    arrange(desc(LB_train),desc(Score_test),LB_covid,level1,level2) %>%
    slice(1) %>%
    unite(id,level1,level2,sep='.') %>%
    pull(id) %>%
    lapply(function(x){
      y=str_split(x,'\\.')[[1]]
      list(
        mod=bestmod_mod[[y[1]]][[y[2]]]
        ,eval=bestmod_edev[[y[1]]][[y[2]]]
        ,train=bestmod_dev[[y[1]]][[y[2]]]
        ,test=bestmod_rep3[[y[1]]][[y[2]]]
        ,covid=
          GSE177477_phenotype %>%
          select(outcome) %>%
          rownames_to_column(var='id') %>%
          left_join(
            GSE177477 %>%
              exprs() %>%
              mb_predictor(list(maternal_blood=GSE108497_maternal_blood)) %>%
              rownames_to_column(var='id')
            ,by='id'
          ) %>%
          left_join(
            GSE177477 %>%
              phenoData() %>%
              pData() %>%
              rownames_to_column(var='gsm') %>%
              select(c('gsm',colnames(.) %>% .[str_detect(.,':ch1')])) %>%
              setNames(
                colnames(.) %>%
                  str_remove_all(':ch1') %>%
                  str_replace_all('\\s','_')
              ) %>%
              select(-id) %>%
              rename(id=gsm) %>%
              select(id,symptons,severity)
            ,by='id'
          ) %>%
          mutate(
            outcome=
              factor(ifelse(
                !symptons%in%c('Healthy controls')
                ,'event','nonevent'
              ),levels(outcome))
          ) %>%
          select(-symptons,-severity) %>%
          column_to_rownames(var='id') %>%
          select_at(colnames(bestmod_dev[[y[1]]][[y[2]]]))
      )
    }) %>%
    .[[1]]
  
  saveRDS(best_each_method,'data/best_each_method.rds')
}else{
  best_each_method=readRDS('data/best_each_method.rds')
}
################################################################################

#####{Candidate and selected predictors of biomarker trees, include=FALSE}
if(load_emulators){
  bestmod_biomarkers=
    bestmod_dev %>%
    pblapply(function(x){
      x %>%
        lapply(function(x){
          x[rownames(x)%in%rownames(
            GSE108497 %>%
              phenoData() %>%
              pData() %>%
              select(c('title',colnames(.) %>% .[str_detect(.,':ch1')])) %>%
              setNames(colnames(.) %>% str_remove_all(':ch1')) %>%
              filter(
                apl==0
                & lac==0
                & fd==0
                & nnd==0
              ) %>%
              mutate(early=as.integer(as.numeric(ga_at_end_of_pregnancy)<34)) %>%
              mutate(time_point=ifelse(is.na(time_point),'non_pregnant',time_point)) %>%
              select(sle,early,pe,sga,time_point) %>%
              filter(
                (early==0 & pe==0 & sga==0)
                |
                  (early==1 & pe==1 & sga==0)
                |
                  (early==0 & pe==1 & sga==0)
                |
                  (early==0 & pe==0 & sga==1)
              ) %>%
              filter(!(early==1 & pe==1 & sga==0))
          ),]
        })
    }) %>%
    get_auroc(bestmod_mod)
  
  saveRDS(bestmod_biomarkers,'data/bestmod_biomarkers.rds')
}else{
  bestmod_biomarkers=readRDS('data/bestmod_biomarkers.rds')
}
################################################################################

#####{The best biomarkers, echo=FALSE, fig.height=3.46457, fig.width=3.46457}
if(load_emulators){
  saveRDS(bestmod_mod$rank2$grank3,'data/best_biomarkers_mod.rds')
}else{
  best_biomarkers_mod=readRDS('data/best_biomarkers_mod.rds')
}
################################################################################

#####{Show averages and SDs to scale the biomarkers, eval=FALSE, include=FALSE}
if(load_emulators){
  saveRDS(bestmod_dev$rank2$grank3,'data/best_biomarkers_dev.rds')
}else{
  best_biomarkers_dev=readRDS('data/best_biomarkers_dev.rds')
}
################################################################################






################################################################################
## miRNA
################################################################################

#####{Load the DIANA-mited miRNA data for PE, include=FALSE}
mirna_pe_fgrpe=
  read_tsv('data/diana_mited_20220302.tsv') %>%
  select(
    -S.No
    ,-Project_ID
    ,-Collection
    ,-Tissue_subregion
    ,-Cell_Line
    ,-Organism
    ,-Gender
    ,-Tissue_definition
  ) %>%
  gather(
    miRNA
    ,log2RPM
    ,-Sample_ID
    ,-Tissue_or_organ_of_origin
    ,-Disease
    ,-Health_state
  ) %>%
  mutate(
    Health_state=
      ifelse(
        Health_state=='disease'
        ,'event'
        ,'nonevent'
      ) %>%
      factor(c('nonevent','event'))
  )
################################################################################

#####{Load the GeneCards miRNA ID targeting the predictor genes, include=FALSE}
mirna_target=
  read_csv('data/mirna.csv') %>%
  separate(mirna,c('miRNA','id'),sep=' \\(') %>%
  mutate(
    id=str_remove_all(id,'\\)')
    ,gene=factor(gene,unique(gene))
  ) %>%
  left_join(mirna_pe_fgrpe,by='miRNA')
################################################################################

#####{Filter PE as outcome of the miRNA of interests, include=FALSE}
mirna_target_pe=
  mirna_target %>%
  filter(
    !(Health_state=='event'
      & Disease=='Preeclampsia, intrauterine growth restriction')
    | is.na(Health_state)
  ) %>%
  mutate(
    Disease=
      ifelse(
        Health_state=='event'
        ,Disease
        ,NA
      ) %>%
      .[!is.na(.)] %>%
      .[!duplicated(.)]
  )
################################################################################

#####{Conduct regression analyses of the miRNAs on PE, include=FALSE}
mirna_target_OR_pe=
  mirna_target_pe %>%
  filter(!is.na(Sample_ID)) %>%
  group_by(gene,miRNA,Disease,Health_state) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  spread(Health_state,n,fill=0) %>%
  mutate(
    total=nonevent+event
    ,p=event/total
  ) %>%
  filter(!(p==0|p==1)) %>%
  left_join(mirna_target_pe,by=c('gene','miRNA','Disease')) %>%
  group_by(gene,miRNA,Disease) %>%
  do(broom::tidy(suppressWarnings(
    glm(Health_state~log2RPM,data=.,family=binomial(link='logit'))
  ))) %>%
  filter(term=='log2RPM') %>%
  mutate(
    LB=round(exp(estimate-std.error),2)
    ,UB=round(exp(estimate+std.error),2)
    ,OR=round(exp(estimate),2)
  ) %>%
  select(gene,miRNA,Disease,OR,LB,UB) %>%
  ungroup()
################################################################################

#####{Determine total effect of miRNAs on PE per predictor, include=FALSE}
mirna_target_reg_pe=
  mirna_target_pe %>%
  mutate(
    availability=
      ifelse(
        miRNA
        %in%pull(
          mirna_target_pe %>%
            filter(!is.na(Sample_ID)) %>%
            group_by(gene,miRNA,Disease,Health_state) %>%
            summarize(n=n()) %>%
            ungroup() %>%
            spread(Health_state,n,fill=0) %>%
            mutate(
              total=nonevent+event
              ,p=event/total
            ) %>%
            filter(!(p==0|p==1))
          ,miRNA
        )
        ,'available'
        ,'not_available'
      )
  ) %>%
  select(gene,miRNA,availability) %>%
  filter(!duplicated(.)) %>%
  group_by(gene,availability) %>%
  summarize(
    n=n()
    ,unavailable_miRNAs=
      paste0(
        ifelse(availability=='not_available',miRNA,NA) %>%
          .[!is.na(.)]
        ,collapse='|'
      )
  ) %>%
  ungroup() %>%
  group_by(gene) %>%
  mutate(
    unavailable_miRNAs=
      paste0(
        unavailable_miRNAs %>%
          .[.!='']
        ,collapse='|'
      )
  ) %>%
  spread(availability,n,fill=0) %>%
  mutate(
    total=available+not_available
    ,p=available/total
    ,availability=
      paste0(
        available
        ,'/'
        ,total
        ,' ('
        ,round(p*100,ifelse(total<10,0,ifelse(total<100,1,2)))
        ,'%)'
      )
  ) %>%
  select(gene,availability,unavailable_miRNAs) %>%
  left_join(
    mirna_target_OR_pe %>%
      mutate(
        significance=
          ifelse(!is.na(OR) & (LB>1 | UB<1),'significant','not_significant')
      ) %>%
      group_by(gene,significance) %>%
      summarize(n=n()) %>%
      spread(significance,n,fill=0) %>%
      mutate(
        total=significant+not_significant
        ,p=significant/total
        ,significance=
          paste0(
            significant
            ,'/'
            ,total
            ,' ('
            ,round(p*100,ifelse(total<10,0,ifelse(total<100,1,2)))
            ,'%)'
          )
      ) %>%
      select(gene,significance)
    ,by='gene'
  ) %>%
  left_join(
    mirna_target_OR_pe %>%
      filter(!is.na(OR) & (LB>1 | UB<1)) %>%
      mutate(
        regulation=ifelse(UB<1,'downregulated','upregulated')
      ) %>%
      group_by(gene,regulation) %>%
      summarize(n=n()) %>%
      spread(regulation,n,fill=0) %>%
      mutate(
        total=downregulated+upregulated
        ,p=downregulated/total
        ,downregulation=
          paste0(
            downregulated
            ,'/'
            ,total
            ,' ('
            ,round(p*100,ifelse(total<10,0,ifelse(total<100,1,2)))
            ,'%)'
          )
      ) %>%
      select(gene,downregulation)
    ,by='gene'
  )

suppressWarnings(set.seed(33,sample.kind=sample.kind))
mirna_target_reg_pe=
  mirna_target_reg_pe %>%
  left_join(
    mirna_target_OR_pe %>%
      filter(!is.na(OR) & (LB>1 | UB<1)) %>%
      mutate(
        regulation=ifelse(UB<1,'downregulated','upregulated')
      ) %>%
      lapply(X=seq(100),Y=.,function(X,Y){
        Y %>%
          mutate(
            regulation=
              regulation %>%
              sample(length(.),F)
          ) %>%
          group_by(gene,regulation) %>%
          summarize(n=n()) %>%
          ungroup() %>%
          group_by(gene) %>%
          mutate(
            total=sum(n)
            ,p_rand=n/total
          ) %>%
          mutate(b=X)
      }) %>%
      do.call(rbind,.) %>%
      select(b,gene,regulation,p_rand) %>%
      spread(regulation,p_rand,fill=0) %>%
      mutate(
        d_rand=downregulated-upregulated
      ) %>%
      select(-downregulated,-upregulated) %>%
      left_join(
        mirna_target_OR_pe %>%
          filter(!is.na(OR) & (LB>1 | UB<1)) %>%
          mutate(
            regulation=ifelse(UB<1,'downregulated','upregulated')
          ) %>%
          group_by(gene,regulation) %>%
          summarize(n=n()) %>%
          ungroup() %>%
          group_by(gene) %>%
          mutate(
            total=sum(n)
            ,p=n/total
          ) %>%
          select(gene,regulation,p) %>%
          spread(regulation,p,fill=0) %>%
          mutate(
            d=downregulated-upregulated
          ) %>%
          select(-downregulated,-upregulated)
        ,by=c('gene')
      ) %>%
      group_by(gene) %>%
      summarize(p.value=mean(abs(d)>abs(d_rand)))
    ,by='gene'
  ) %>%
  left_join(
    mirna_target_OR_pe %>%
      filter(!is.na(OR) & (LB>1 | UB<1)) %>%
      mutate_at(c('OR','LB','UB'),format,digits=2,scientific=T)
    ,by='gene'
  )
################################################################################

#####{Filter PE-FGR as outcome of the miRNA of interests, include=FALSE}
mirna_target_fgrpe=
  mirna_target %>%
  filter(
    !(Health_state=='event'
      & Disease=='Preeclampsia')
    | is.na(Health_state)
  ) %>%
  mutate(
    Disease=
      ifelse(
        Health_state=='event'
        ,Disease
        ,NA
      ) %>%
      .[!is.na(.)] %>%
      .[!duplicated(.)]
  )
################################################################################

#####{Conduct regression analyses of the miRNAs on PE-FGR, include=FALSE}
mirna_target_OR_fgrpe=
  mirna_target_fgrpe %>%
  filter(!is.na(Sample_ID)) %>%
  group_by(gene,miRNA,Disease,Health_state) %>%
  summarize(n=n()) %>%
  ungroup() %>%
  spread(Health_state,n,fill=0) %>%
  mutate(
    total=nonevent+event
    ,p=event/total
  ) %>%
  filter(!(p==0|p==1)) %>%
  left_join(mirna_target_fgrpe,by=c('gene','miRNA','Disease')) %>%
  group_by(gene,miRNA,Disease) %>%
  do(broom::tidy(suppressWarnings(
    glm(Health_state~log2RPM,data=.,family=binomial(link='logit'))
  ))) %>%
  filter(term=='log2RPM') %>%
  mutate(
    LB=round(exp(estimate-std.error),2)
    ,UB=round(exp(estimate+std.error),2)
    ,OR=round(exp(estimate),2)
  ) %>%
  select(gene,miRNA,Disease,OR,LB,UB) %>%
  ungroup()
################################################################################

#####{Determine total effect of miRNAs on PE-FGR per predictor, include=FALSE}
mirna_target_reg_fgrpe=
  mirna_target_fgrpe %>%
  mutate(
    availability=
      ifelse(
        miRNA
        %in%pull(
          mirna_target_fgrpe %>%
            filter(!is.na(Sample_ID)) %>%
            group_by(gene,miRNA,Disease,Health_state) %>%
            summarize(n=n()) %>%
            ungroup() %>%
            spread(Health_state,n,fill=0) %>%
            mutate(
              total=nonevent+event
              ,p=event/total
            ) %>%
            filter(!(p==0|p==1))
          ,miRNA
        )
        ,'available'
        ,'not_available'
      )
  ) %>%
  select(gene,miRNA,availability) %>%
  filter(!duplicated(.)) %>%
  group_by(gene,availability) %>%
  summarize(
    n=n()
    ,unavailable_miRNAs=
      paste0(
        ifelse(availability=='not_available',miRNA,NA) %>%
          .[!is.na(.)]
        ,collapse='|'
      )
  ) %>%
  ungroup() %>%
  group_by(gene) %>%
  mutate(
    unavailable_miRNAs=
      paste0(
        unavailable_miRNAs %>%
          .[.!='']
        ,collapse='|'
      )
  ) %>%
  spread(availability,n,fill=0) %>%
  mutate(
    total=available+not_available
    ,p=available/total
    ,availability=
      paste0(
        available
        ,'/'
        ,total
        ,' ('
        ,round(p*100,ifelse(total<10,0,ifelse(total<100,1,2)))
        ,'%)'
      )
  ) %>%
  select(gene,availability,unavailable_miRNAs) %>%
  left_join(
    mirna_target_OR_fgrpe %>%
      mutate(
        significance=
          ifelse(!is.na(OR) & (LB>1 | UB<1),'significant','not_significant')
      ) %>%
      group_by(gene,significance) %>%
      summarize(n=n()) %>%
      spread(significance,n,fill=0) %>%
      mutate(
        total=significant+not_significant
        ,p=significant/total
        ,significance=
          paste0(
            significant
            ,'/'
            ,total
            ,' ('
            ,round(p*100,ifelse(total<10,0,ifelse(total<100,1,2)))
            ,'%)'
          )
      ) %>%
      select(gene,significance)
    ,by='gene'
  ) %>%
  left_join(
    mirna_target_OR_fgrpe %>%
      filter(!is.na(OR) & (LB>1 | UB<1)) %>%
      mutate(
        regulation=ifelse(UB<1,'downregulated','upregulated')
      ) %>%
      group_by(gene,regulation) %>%
      summarize(n=n()) %>%
      spread(regulation,n,fill=0) %>%
      mutate(
        total=downregulated+upregulated
        ,p=downregulated/total
        ,downregulation=
          paste0(
            downregulated
            ,'/'
            ,total
            ,' ('
            ,round(p*100,ifelse(total<10,0,ifelse(total<100,1,2)))
            ,'%)'
          )
      ) %>%
      select(gene,downregulation)
    ,by='gene'
  )

suppressWarnings(set.seed(33,sample.kind=sample.kind))
mirna_target_reg_fgrpe=
  mirna_target_reg_fgrpe %>%
  left_join(
    mirna_target_OR_fgrpe %>%
      filter(!is.na(OR) & (LB>1 | UB<1)) %>%
      mutate(
        regulation=ifelse(UB<1,'downregulated','upregulated')
      ) %>%
      lapply(X=seq(100),Y=.,function(X,Y){
        suppressWarnings(set.seed(33+X,sample.kind=sample.kind))
        Y %>%
          mutate(
            regulation=
              regulation %>%
              sample(length(.),F)
          ) %>%
          group_by(gene,regulation) %>%
          summarize(n=n()) %>%
          ungroup() %>%
          group_by(gene) %>%
          mutate(
            total=sum(n)
            ,p_rand=n/total
          ) %>%
          mutate(b=X)
      }) %>%
      do.call(rbind,.) %>%
      select(b,gene,regulation,p_rand) %>%
      spread(regulation,p_rand,fill=0) %>%
      mutate(
        d_rand=downregulated-upregulated
      ) %>%
      select(-downregulated,-upregulated) %>%
      left_join(
        mirna_target_OR_fgrpe %>%
          filter(!is.na(OR) & (LB>1 | UB<1)) %>%
          mutate(
            regulation=ifelse(UB<1,'downregulated','upregulated')
          ) %>%
          group_by(gene,regulation) %>%
          summarize(n=n()) %>%
          ungroup() %>%
          group_by(gene) %>%
          mutate(
            total=sum(n)
            ,p=n/total
          ) %>%
          select(gene,regulation,p) %>%
          spread(regulation,p,fill=0) %>%
          mutate(
            d=downregulated-upregulated
          ) %>%
          select(-downregulated,-upregulated)
        ,by=c('gene')
      ) %>%
      group_by(gene) %>%
      summarize(p.value=mean(abs(d)>abs(d_rand)))
    ,by='gene'
  ) %>%
  left_join(
    mirna_target_OR_fgrpe %>%
      filter(!is.na(OR) & (LB>1 | UB<1)) %>%
      mutate_at(c('OR','LB','UB'),format,digits=2,scientific=T)
    ,by='gene'
  )
################################################################################






################################################################################
## PPI
################################################################################

#####{Create a function to show the shortest path, include=FALSE}
source('R/show_shortest_path-function.R')
################################################################################

#####{Load the STRING PPI data of the biomarkers and surrogates, include=FALSE}
ppi_string=
  list(
    blood_20220303=c('ITGA5','P2RX7','IRF6')
    ,blcord_20220303=c('ITGA5','FANCI','TSEN15')
    ,bldeci_20220303=c('ITGA5','TPX2','WIPF3','ARID2')
    ,blfund_20220303=c('ITGA5','ARID2')
    ,blamni_20220303=c('P2RX7','INSM1')
    ,blamni2_20220303=c('IRF6','ALS2CL')
    ,cordam_20220304=c('SELV','INSM1','ALS2CL','TMEM38B')
    ,bldeciplac_20220314=
      c('ITGA5','TPX2','WIPF3','ARID2'
        ,'IRF6','TMEM38B','PLOD2'
      )
    ,shopath_all_20220304=
      c('ITGA5'
        ,'P2RX7'
        ,'IRF6'
        ,'FANCI'
        ,'TSEN15'
        ,'TPX2'
        ,'WIPF3'
        ,'ARID2'
        ,'INSM1'
        ,'ALS2CL'
        ,'SELV'
        ,'TMEM38B'
      )
    ,shopath_blood_20220304=
      c('ITGA5'
        ,'P2RX7'
        ,'IRF6'
      )
    ,shopath_bldecifund_20220304=
      c('ITGA5'
        ,'TPX2'
        ,'WIPF3'
        ,'ARID2'
      )
    ,shopath_blcordamplac_20220304=
      c('ITGA5'
        ,'P2RX7'
        ,'IRF6'
        ,'INSM1'
        ,'ALS2CL'
        ,'SELV'
        ,'TMEM38B'
      )
    ,shopath_bldeci_20220304=
      c('ITGA5'
        ,'TPX2'
        ,'WIPF3'
        ,'ARID2'
      )
    ,shopath_blfund_20220304=
      c('ITGA5'
        ,'ARID2'
      )
    ,shopath_blcord_20220310=
      c('ITGA5'
        ,'FANCI'
        ,'TSEN15'
      )
    ,shopath_blamni_20220310=
      c('P2RX7'
        ,'IRF6'
        ,'INSM1'
        ,'ALS2CL'
      )
    ,shopath_amniplac_20220310=
      c('INSM1'
        ,'ALS2CL'
        ,'TMEM38B'
      )
    ,shopath_amnicord_20220310=
      c('INSM1'
        ,'ALS2CL'
        ,'SELV'
      )
    ,shopath_cordplac_20220310=
      c('SELV'
        ,'TMEM38B'
      )
    ,shopath_bldeciplac_20220314=
      c('ITGA5','TPX2','WIPF3','ARID2'
        ,'IRF6','TMEM38B','PLOD2'
      )
  )

ppi_string=
  ppi_string %>%
  lapply(X=names(.),Y=.,function(X,Y){
    list(
      inputs=Y[[X]]
      ,nodes=read_tsv(paste0('data/',X,'_string_network_coordinates.tsv'))
      ,edges=read_tsv(paste0('data/',X,'_string_interactions.tsv'))
      ,pathways=read_tsv(paste0('data/',X,'_enrichment.all.tsv'))
      ,pubmeds=read_tsv(paste0('data/',X,'_enrichment.PMID.tsv'))
    )
  }) %>%
  setNames(names(ppi_string))
################################################################################

#####{Identify and show the shortest paths, include=FALSE}
ppi_string=
  ppi_string %>%
  pblapply(X=names(.),Y=.,function(X,Y){
    Y[[X]] %>%
      .[!names(.)%in%c('combination','shortest_path')] %>%
      c(show_shortest_path(Y[[X]],X))
  }) %>%
  setNames(names(ppi_string))
################################################################################

#####{Select the pathways overrepresented by the shortest paths, include=FALSE}
if(run_heavy_computation){
  ppi_pathways=
    ppi_string %>%
    lapply(X=names(.),Y=.,function(X,Y){
      Z=Y[[X]]$pubmeds %>%
        mutate(`#category`='PubMed') %>%
        rename(`term ID`=`#term ID`) %>%
        select(`#category`,everything()) %>%
        list() %>%
        c(list(Y[[X]]$pathways)) %>%
        .[2:1] %>%
        do.call(rbind,.) %>%
        `colnames<-`(
          colnames(.) %>%
            str_remove_all('[:punct:]') %>%
            str_to_lower() %>%
            str_replace_all('\\s+','_')
        ) %>%
        mutate(id=X) %>%
        select(id,everything())
      
      K=Y[[X]]$combination %>%
        lapply(X=seq(nrow(.)),Y=.,Z=Y[[X]]$shortest_path,function(X,Y,Z){
          K=Z[[X]] %>%
            select(name,x,y) %>%
            rename(node1=name) %>%
            filter(!duplicated(.)) %>%
            left_join(
              Z[[X]] %>%
                select(name,xend,yend) %>%
                rename(node2=name,x=xend,y=yend) %>%
                filter(!duplicated(.))
              ,by=c('x','y')
            ) %>%
            select(-x,-y) %>%
            graph_from_data_frame() %>%
            shortest_paths(Y$Var1[X],Y$Var2[X]) %>%
            .$vpath %>%
            unlist() %>%
            names()
          
          if(length(K)<=2){
            paste0(K,collapse=',')
          }else{
            seq(2,length(K)) %>%
              lapply(function(x){
                seq(1,length(K)-x+1,1) %>%
                  lapply(function(y){
                    K[seq(y,y+x-1)] %>%
                      paste0(collapse=',')
                  })
              }) %>%
              unlist()
          }
        }) %>%
        unlist() %>%
        sapply(function(x){
          paste0(sort(str_split(x,',')[[1]]),collapse=',')
        }) %>%
        .[!duplicated(.)]
      
      L=K %>%
        pblapply(function(x){
          Z %>%
            mutate(
              matched=
                matching_proteins_in_your_network_labels %>%
                sapply(function(x){
                  paste0(sort(str_split(x,',')[[1]]),collapse=',')
                })
              ,matched=matched==x
            ) %>%
            select(matched) %>%
            setNames(x)
        }) %>%
        do.call(cbind,.)
      
      cbind(Z,L) %>%
        pivot_longer(
          seq(ncol(Z)+1,ncol(Z)+ncol(L))
          ,names_to='gene'
          ,values_to='matched'
        ) %>%
        group_by_at(seq(1,ncol(Z))) %>%
        summarize(
          genes=
            ifelse(matched,gene,NA) %>%
            .[!is.na(.)] %>%
            sort() %>%
            paste0(collapse='|')
          ,matched=sum(matched)
        ) %>%
        ungroup() %>%
        arrange(factor(term_id,Z$term_id)) %>%
        mutate(
          background_gene_class=
            ifelse(
              matched==1
              & background_gene_count>=15 & background_gene_count<=200
              ,1
              ,ifelse(
                matched==1
                & background_gene_count>=10 & background_gene_count<=500
                ,2
                ,3
              )
            )
        ) %>%
        filter(matched==1) %>%
        group_by(genes) %>%
        filter(background_gene_class==min(background_gene_class)) %>%
        ungroup() %>%
        mutate(
          genes=factor(genes,unique(genes))
          ,category=factor(category,unique(category))
        ) %>%
        group_by(genes,category) %>%
        arrange(desc(strength),desc(observed_gene_count)) %>%
        filter(
          strength==max(strength)
          & observed_gene_count==max(observed_gene_count)
        ) %>%
        ungroup() %>%
        mutate_at(c('genes','category'),as.character) %>%
        select(
          id,genes,observed_gene_count,term_id,category
          ,background_gene_class
          ,everything()
        )
    }) %>%
    setNames(names(ppi_string))
  
  saveRDS(ppi_pathways,'data/ppi_pathways.rds')
}else{
  cat(readRDS('data/log.rds')[['ppi_pathways']])
  ppi_pathways=readRDS('data/ppi_pathways.rds')
}
################################################################################






################################################################################
## PTM
################################################################################

#####{Load PTM and protein data, include=FALSE}
ptm_genes=read_csv('data/ptm.csv')
ptm_effects=
  read_csv('data/ptm_effect.csv') %>%
  # https://en.wikipedia.org/wiki/ADP-ribosylation
  rbind(
    data.frame(
      ptm='ADP-ribosylation'
      ,effect=
        c('cell signaling'
          ,'DNA repair'
          ,'transcription' # gene regulation
          ,'apoptosis')
      ,covalent_attachment='small chemical group'
    )
  )
protein_exprs=read_csv('data/protein.csv')
################################################################################
