##keeping loading scripts in this repo
library(amlresistancenetworks)
library(dplyr)





#' getCytokine data
#' gets cytokine sensitivity data and metadata in one call
#' @import dplyr
#' @import readxl
getCytokineSensData<-function(){
  library(dplyr)
  syn<-synapseLogin()
  metadata<-readxl::read_xlsx(syn$get('syn22175391')$path)%>%
    mutate(sample=stringr::str_replace(`Sample name (DMS)`,'PTRC_Ex15_',''))%>%
    dplyr::select(sample,CellType,TimePoint,Treatment)%>%
    subset(!is.na(sample))
  
  ##there is an error in the metadata, fixing here
  metadata[grep('M13_PAR_02',metadata$sample),'TimePoint']<-0
  
  metadata$CellType<-sapply(metadata$CellType,function(x) {
    switch(x,Late_M='Late MOLM-13',Late_MR='Late MOLM-13 Tr Resistant',
           M='MOLM-13',MR='MOLM-13 Tr Resistant')})
  
  #  ifelse(x=='Late_M','Late MOLM-13',ifelse(x=='Late_MR','Late MOLM-13 Tr Resistant',x))})
  metadata$Treatment<-unlist(sapply(metadata$Treatment,function(x){
    switch(x,M='MCP-1',T='Trametinib',TW='Trametinib Withdrawn',none="none",`T+M`='Trametinib+MCP-1')}))
  
  ##this is experiment 2
  pdat<-read.csv2(syn$get('syn22862628')$path,sep='\t')%>%
    tidyr::pivot_longer(-c(Protein,Gene),names_to='tsamp',values_to='LogRatio')%>%
    mutate(sample=stringr::str_replace(tsamp,'X',''))%>%
    dplyr::select(-tsamp)%>%
    left_join(metadata)%>%
    subset(sample!='Peptide')
  
  pdat$Batch='Experiment 2'
  
  phdat<-read.csv2(syn$get('syn22862617')$path,sep='\t')%>%
    tidyr::pivot_longer(-c(Protein,Gene,site,Peptide),
                        names_to='tsamp',values_to='LogRatio')%>%
    mutate(sample=stringr::str_replace(tsamp,'X',''))%>%
    dplyr::select(-tsamp)%>%
    left_join(metadata)
  phdat$Batch='Experiment 2'
  
  pdat2<-read.csv2(syn$get('syn22173207')$path,sep='\t')%>%
    tidyr::pivot_longer(-c(Protein,Gene),names_to='tsamp',values_to='LogRatio')%>%
    mutate(sample=stringr::str_replace(tsamp,'X',''))%>%
    dplyr::select(-tsamp)%>%
    left_join(metadata)
  pdat2$Batch='Experiment 1'
  
  phdat2<-read.csv2(syn$get('syn22173201')$path,sep='\t')%>%
    tidyr::pivot_longer(-c(Protein,Gene,site,Peptide),
                        names_to='tsamp',values_to='LogRatio')%>%
    mutate(sample=stringr::str_replace(tsamp,'X',''))%>%
    dplyr::select(-tsamp)%>%
    left_join(metadata)
  phdat2$Batch='Experiment 1'
  
  fullp<-rbind(pdat,pdat2)
  p.synid='syn22986326'
  #phrows=syn$tableQuery(paste0('select * from ',ph.synid))
  #syn$delete(phrows)
  #syn$store(p.synid,fullp)
  synTableUpdate(fullp,p.synid)
  
  fullph<-rbind(phdat,phdat2)
  ph.synid='syn22986341'
  
  #phrows=syn$tableQuery(paste0('select * from ',ph.synid))
  #syn$delete(phrows)
  synTableUpdate(fullph,ph.synid)
  
  #synTableStore(rbind(pdat2,pdat),'Cytokine-induced Drug Sensitivity Proteomics')
  #synTableStore(rbind(phdat2,phdat),'Cytokine-induced Drug Sensitivity Phospho-proteomics')
  
}


updateCytokinePhospho<-function(){
  library(dplyr)
  syn<-synapseLogin()
  metadata<-readxl::read_xlsx(syn$get('syn22175391')$path)%>%
    mutate(sample=stringr::str_replace(`Sample name (DMS)`,'PTRC_Ex15_',''))%>%
    dplyr::select(sample,CellType,TimePoint,Treatment)%>%
    subset(!is.na(sample))
  
  metadata$CellType<-sapply(metadata$CellType,function(x) {
    switch(x,Late_M='Late MOLM-13',Late_MR='Late MOLM-13 Tr Resistant',
           M='MOLM-13',MR='MOLM-13 Tr Resistant')})
  
  #  ifelse(x=='Late_M','Late MOLM-13',ifelse(x=='Late_MR','Late MOLM-13 Tr Resistant',x))})
  metadata$Treatment<-unlist(sapply(metadata$Treatment,function(x){
    switch(x,M='MCP-1',T='Trametinib',TW='Trametinib Withdrawn',none="none",`T+M`='Trametinib+MCP-1')}))
  
  ##there is an error in the metadata, fixing here
  metadata[grep('M13_PAR_02',metadata$sample),'TimePoint']<-0
  
  phdat<-read.csv2(syn$get('syn24305135')$path,sep='\t')%>%
    tidyr::pivot_longer(-c(Protein,Gene,site,Peptide),
                        names_to='tsamp',values_to='LogRatio')%>%
    mutate(sample=stringr::str_replace(tsamp,'X',''))%>%
    dplyr::select(-tsamp)%>%
    left_join(metadata)
  
  phdat$Batch='Experiment 2'
  
  phdat2<-read.csv2(syn$get('syn24375434')$path,sep='\t')%>% ##TBD
    tidyr::pivot_longer(-c(Protein,Gene,site,Peptide),
                        names_to='tsamp',values_to='LogRatio')%>%
    mutate(sample=stringr::str_replace(tsamp,'X',''))%>%
    dplyr::select(-tsamp)%>%
    left_join(metadata)
  phdat2$Batch='Experiment 1'
  
  synid='syn24389738'
  #phrows=syn$tableQuery(paste0('select * from ',synid))
  #syn$delete(phrows)
  #syn$store(synid,rbind(phdat2,phdat))
  synTableUpdate(rbind(phdat2,phdat),synid)
  #synTableStore(rbind(phdat2,phdat),'Cytokine-induced Drug Sensitivity Phosphoproteomics Unnormalized')
  
  
}

