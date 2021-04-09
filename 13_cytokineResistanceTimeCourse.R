##process cytokine data
library(amlresistancenetworks)
library(dplyr)


###load all the data

phosData<-querySynapseTable('syn22986341')%>%subset(!is.nan(LogRatio))%>%
  mutate(Gene=unlist(Gene))%>%
  mutate(site=unlist(site))

uncorrectedPhosData<-querySynapseTable('syn24389738')%>%subset(!is.nan(LogRatio))%>%
  mutate(Gene=unlist(Gene))%>%
  mutate(site=unlist(site))

clinvars<-phosData%>%
  dplyr::select(Sample='sample',CellType,TimePoint,Treatment)%>%
  distinct()


##what are we doing again?
summary<-phosData%>%
  dplyr::select(sample,CellType,TimePoint,Treatment)%>%
  distinct()%>%
  mutate(conditionName=stringr::str_c(CellType,TimePoint,Treatment,sep='_'))

print(summary)

phosMat<-phosData%>%dplyr::select(sample,site,LogRatio)%>%
  tidyr::pivot_wider(values_from=LogRatio,names_from=sample,
                     values_fn=list(LogRatio=mean),values_fill=list(LogRatio=0.0))%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames('site')

fullPhosMat<-uncorrectedPhosData%>%dplyr::select(sample,site,LogRatio)%>%
  tidyr::pivot_wider(values_from=LogRatio,names_from=sample,
                     values_fn=list(LogRatio=mean),values_fill=list(LogRatio=0.0))%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames('site')


##
#' @param dat.table
plotAllData<-function(dat.table){
  library(ggfortify)
  met<-dat.table%>%dplyr::select(sample,CellType,TimePoint,Treatment)%>%
    distinct()
  #%>%
  #  tibble::column_to_rownames('sample')
    
  mat<-dat.table%>%dplyr::select(Gene,LogRatio,sample)%>%
    distinct()%>%
    mutate(LogRatio=as.numeric(LogRatio))%>%
    tidyr::pivot_wider(names_from='sample',values_from='LogRatio',values_fn=list(LogRatio=function(x) mean(x,na.rm=T)),values_fill=list(LogRatio=0))%>%
  tibble::remove_rownames()%>%
    tibble::column_to_rownames('Gene')
  
  autoplot(prcomp(t(mat)),data=met,colour='Treatment',shape='CellType')
 
}

##plot kinase activity
plotKinDat<-function(kindat,prefix='all'){
  library(pheatmap)
  ##create matrix of kinase scores
  mat <-kindat%>%
    ungroup()%>%
    tidyr::pivot_wider(-c(meanNKINscore,numSubstr),
                                              values_from=meanLFC,
                                                names_from=Sample,
                                                values_fn=list(meanLFC=mean))%>%
    tibble::remove_rownames()%>%
    tibble::column_to_rownames('Kinase')
  
  kinAts<-kindat%>%
    ungroup()%>%
    dplyr::select(Kinase,numSubstr)%>%
    distinct()%>%
    group_by(Kinase)%>%
    summarize(substrates=mean(numSubstr))%>%
    tibble::remove_rownames()%>%
    tibble::column_to_rownames('Kinase')
  
  sampAts<-phosData%>%
    dplyr::select(sample,CellType,TimePoint,Treatment)%>%
    distinct()%>%
    tibble::remove_rownames()%>%
    tibble::column_to_rownames('sample')
  
  sampAts$TimePoint=as.factor(sampAts$TimePoint)
  
  vars=names(sort(apply(mat,1,var),decreasing=T)[1:150])
 
  pheatmap(mat[vars,],cellwidth = 8,cellheight=8,clustering_distance_cols = 'correlation',
          clustering_distance_rows = 'correlation',
          annotation_row = kinAts,annotation_col=sampAts,
          file=paste0(prefix,'cytokineKinaseHeatmap.pdf'),height=20,width=8) 
}


####Show kinase activity
kindat<-mapPhosphoToKinase(dplyr::rename(phosData,Sample='sample', LogFoldChange='LogRatio'))
parental<-mapPhosphoToKinase(dplyr::rename(filter(phosData,CellType=='MOLM-13'),Sample='sample', LogFoldChange='LogRatio'))
uncorrectedKinDat<-mapPhosphoToKinase(dplyr::rename(uncorrectedPhosData,Sample='sample', LogFoldChange='LogRatio'))

plotKinDat(kindat)
plotKinDat(parental,'molm13')
plotKinDat(uncorrectedKinDat,'uncorrected')

plots=list(plotAllData(protData),plotAllData(phosData),plotAllData(uncorrectedPhosData))
cowplot::plot_grid(plotlist=plots,labels=c("Bulk Proteomics",'Phosphoprotomics','uncorrected Phospho'),nrow=2)
ggsave('pcaOfSamples.png')

l

#' plot all the KSEA 
#' @param condList
#' @return data frame
doAllKSEAplots<-function(condList,pdat=phosData){
  
  gene.to.site<-dplyr::select(pdat,Gene,site,Peptide)%>%distinct()%>%
    dplyr::mutate(residue=stringr::str_replace(site,paste0(Gene,'-'),''))%>%
    dplyr::mutate(residue=stringr::str_replace_all(residue,"([STY])", ";\\1"))%>%
    dplyr::mutate(residue=stringr::str_replace(residue,"^;", ""))%>%
    dplyr::mutate(residue=stringr::str_replace_all(residue,"([sty])", ""))
  
  full.df<-purrr::map_df(names(condList),.f=function(clName){ 
    condList[[clName]]%>%
      tibble::rownames_to_column('site')%>%
      left_join(gene.to.site)%>%
      dplyr::select(Gene,Peptide,residue,value='logFC',p_adj='adj.P.Val')%>%
      amlresistancenetworks::computeKSEA(.,prefix=clName,0.05)%>%
      mutate(Condition=clName)%>%
      as.data.frame()
  })
  return(full.df)
  
}



   

#' build networks from data frame
#' @param data.res
#' @param gene.col
#' @param weight.col
#' @param condition.col
#' @return network list?
runNetworksFromDF<-function(data,gene.col='Kinase.Gene',
                              weight.col='aveSubstrateLog2FC',
                              condition.col='Condition',extra.col=c('Substrate.Gene','Source','log2FC'),
                              signif=0.05){
  res = data%>%
   # dplyr::select(cond=condition.col,value=weight.col,Gene=gene.col,p.value)%>%
    mutate(signif=p.value<signif)%>%
      dplyr::select(c(condition.col,weight.col,gene.col,'signif',extra.col))%>%distinct()%>%
    dplyr::rename(cond=condition.col,value=weight.col,Gene=gene.col)%>%
    group_by(cond)%>%
    dplyr::select(c('cond','Gene','value',extra.col,'signif'))%>%
    group_map(~ amlresistancenetworks::computeProteinNetwork(.x),keep=TRUE)
  return(res)
}

#####now do various comparisons
timeCoursePhospho<-function(){
  
  latePhos<-list(lateTram_vs_lateCombo=limmaTwoFactorDEAnalysis(phosMat,
                                                          filter(summary,conditionName=='Late MOLM-13_0_Trametinib')$sample,
                                                          filter(summary,conditionName=='Late MOLM-13_0_Trametinib+MCP-1')$sample),
               m13_vs_lateTram=limmaTwoFactorDEAnalysis(phosMat,
                                                filter(summary,conditionName=='MOLM-13_0_none')$sample,
                                                filter(summary,conditionName=='Late MOLM-13_0_Trametinib')$sample),
               m13_vs_lateCombo=limmaTwoFactorDEAnalysis(phosMat,
                                                              filter(summary,conditionName=='MOLM-13_0_none')$sample,
                                                                filter(summary,conditionName=='Late MOLM-13_0_Trametinib+MCP-1')$sample))
 
  earlyLatePhos<-list(early_60m_combo_vs_late_combo=limmaTwoFactorDEAnalysis(phosMat,
                                                                         filter(summary,conditionName=='MOLM-13_60_Trametinib+MCP-1')$sample,
                                   filter(summary,conditionName=='Late MOLM-13_0_Trametinib+MCP-1')$sample),
                      early_5m_combo_vs_late_combo=limmaTwoFactorDEAnalysis(phosMat,
                                                                         filter(summary,conditionName=='MOLM-13_5_Trametinib+MCP-1')$sample,
                                   filter(summary,conditionName=='Late MOLM-13_0_Trametinib+MCP-1')$sample),
                      early_both_combo_vs_late_combo=limmaTwoFactorDEAnalysis(phosMat,
                                                                               filter(summary,conditionName%in%c('MOLM-13_5_Trametinib+MCP-1','MOLM-13_60_Trametinib+MCP-1'))$sample,
                                                                              filter(summary,conditionName=='Late MOLM-13_0_Trametinib+MCP-1')$sample))
                      
 
  ph3<-doAllKSEAplots(latePhos)
  resdf<-do.call(rbind,lapply(names(lateProt),function(x) data.frame(lateProt[[x]],Condition=x)))
 
  earlyLateP<-plotConditionsInFlow(earlyLateProt,title='Bulk Proteomics in early vs late',0.05)
  ggsave('earlyLateProt.png',earlyLateP,width=11,height=6)
 # r2<-doAllGOplots(earlyLateProt)              
  
  earlyLatePh<-plotConditionsInFlow(earlyLatePhos,title='Phosphoproteomics in late',0.05)
  ggsave('earlyLatePhos.png',earlyLatePh,width=11,height=6)
  
  earlyLatePhresdf<-do.call(rbind,lapply(names(earlyLatePhos),function(x) data.frame(earlyLatePhos[[x]],Condition=x)))
  
  ph4<-doAllKSEAplots(earlyLatePhos)
 # lateNets<-runNetworksFromDF(ph4)
  
  #lateHeatmap<-kseaZscoreHeatmap(list(ph3,ph4),'lateResistantTreatmentKSEAzscoreHeatmap.pdf')
  return(list(ph3,ph4))

}




##this looks at the early resistant cell lines
compareTramTreatmentCombos<-function(){
  #' here we get various treatments of MCP1 and tram at different time points

  m13Phos<-list(molm13_vs_tram_5min=limmaTwoFactorDEAnalysis(phosMat,
                                                             filter(summary,conditionName=='MOLM-13_0_none')$sample,
                                                             filter(summary,conditionName=='MOLM-13_5_Trametinib')$sample),
                molm13_vs_tram_60min=limmaTwoFactorDEAnalysis(phosMat,
                                                              filter(summary,conditionName=='MOLM-13_0_none')$sample,
                                                              filter(summary,conditionName=='MOLM-13_60_Trametinib')$sample),
                molm13_vs_MCP1_5min=limmaTwoFactorDEAnalysis(phosMat,
                                                             filter(summary,conditionName=='MOLM-13_0_none')$sample,
                                                             filter(summary,conditionName=='MOLM-13_5_MCP-1')$sample),
                molm13_vs_MCP1_60min=limmaTwoFactorDEAnalysis(phosMat,
                                                              filter(summary,conditionName=='MOLM-13_0_none')$sample,
                                                              filter(summary,conditionName=='MOLM-13_60_MCP-1')$sample),
                molm13_vs_MCP1_tram_5min=limmaTwoFactorDEAnalysis(phosMat,
                                                                  filter(summary,conditionName=='MOLM-13_0_none')$sample,
                                                                  filter(summary,conditionName=='MOLM-13_5_Trametinib+MCP-1')$sample),
                molm13_vs_MCP1_tram_60min=limmaTwoFactorDEAnalysis(phosMat,
                                                                   filter(summary,conditionName=='MOLM-13_0_none')$sample,
                                                                   filter(summary,conditionName=='MOLM-13_60_Trametinib+MCP-1')$sample))

 
  p3<-doAllKSEAplots(m13Phos)
  

  tramMCPPhos<-list(tram_vs_mcp1tram_5min=limmaTwoFactorDEAnalysis(phosMat,
                                                                     filter(summary,conditionName=='MOLM-13_5_Trametinib')$sample,
                                                                     filter(summary,conditionName=='MOLM-13_5_Trametinib+MCP-1')$sample),
                      tram_vs_mcp1tram_60min=limmaTwoFactorDEAnalysis(phosMat,
                                                                      filter(summary,conditionName=='MOLM-13_60_Trametinib')$sample,
                                                                      filter(summary,conditionName=='MOLM-13_60_Trametinib+MCP-1')$sample),
                      mcp1_vs_mcp1tram_5min=limmaTwoFactorDEAnalysis(phosMat,
                                                                     filter(summary,conditionName=='MOLM-13_5_MCP-1')$sample,
                                                                     filter(summary,conditionName=='MOLM-13_5_Trametinib+MCP-1')$sample),
                      mcp1_vs_mcp1tram_60min=limmaTwoFactorDEAnalysis(phosMat,
                                                                      filter(summary,conditionName=='MOLM-13_60_MCP-1')$sample,
                                                                      filter(summary,conditionName=='MOLM-13_60_Trametinib+MCP-1')$sample),
                    late_tram_vs_mcp1tram_5min=manualDEAnalysis(phosMat,
                                                                   filter(summary,conditionName=='Late MOLM-13_0_Trametinib')$sample,
                                                                   filter(summary,conditionName=='MOLM-13_5_Trametinib+MCP-1')$sample),
                    late_tram_vs_mcp1tram_60min=manualDEAnalysis(phosMat,
                                                                    filter(summary,conditionName=='Late MOLM-13_0_Trametinib')$sample,
                                                                    filter(summary,conditionName=='MOLM-13_60_Trametinib+MCP-1')$sample),
                    late_mcp1_vs_mcp1tram_5min=manualDEAnalysis(phosMat,
                                                                   filter(summary,conditionName=='Late MOLM-13_0_Trametinib+MCP-1')$sample,
                                                                   filter(summary,conditionName=='MOLM-13_5_Trametinib+MCP-1')$sample),
                    late_mcp1_vs_mcp1tram_60min=manualDEAnalysis(phosMat,
                                                                    filter(summary,conditionName=='Late MOLM-13_0_Trametinib+MCP-1')$sample,
                                                                    filter(summary,conditionName=='MOLM-13_60_Trametinib+MCP-1')$sample)
                    )
  
  p5<-doAllKSEAplots(tramMCPPhos)
 # tramMCP=runNetworksFromDF(ph3)
  #lateHeatmap<-kseaZscoreHeatmap(list(p3,p5),'earlyTreatmentKSEAzscoreHeatmap.pdf')
  
  ph3<-plotConditionsInFlow(tramMCPPhos,title='Phospho effects of tram vs combo',0.05)
  ggsave('phoscomboVsMCPTramIndiv.png',ph3,width=11,height=6)
  return(list(p3,p5))
}


compareTramInResistCells<-function(){
  #' how does trametinib affect resistant cells?
##des MCP1 loo like resistant cells?

  mcp1ResistPhos<-list(resist_vs_mcp1_5min=limmaTwoFactorDEAnalysis(phosMat,
                                                                filter(summary,conditionName=='MOLM-13 Tr Resistant_0_none')$sample,
                                                                filter(summary,conditionName=='MOLM-13 Tr Resistant_5_MCP-1')$sample),
                   resist_vs_mcp1_60min=limmaTwoFactorDEAnalysis(phosMat,
                                                                 filter(summary,conditionName=='MOLM-13 Tr Resistant_0_none')$sample,
                                                                 filter(summary,conditionName=='MOLM-13 Tr Resistant_60_MCP-1')$sample))
  
  p3<-plotConditionsInFlow(mcp1ResistPhos,title='Effects of tram/mcp-1 in resistant',0.05)
  ggsave('mcp1PhosInResistantCells.png',p3,width=11,height=6)
  ph3<-doAllKSEAplots(mcp1ResistPhos)
  #mcp1Resistnetworks<-runNetworksFromDF(ph3)
}

##now pool time points for bulk differences

#
##now let's use the heatmap code!
tcvals<-compareTramTreatmentCombos()
othervals<-doLateComparisons()
