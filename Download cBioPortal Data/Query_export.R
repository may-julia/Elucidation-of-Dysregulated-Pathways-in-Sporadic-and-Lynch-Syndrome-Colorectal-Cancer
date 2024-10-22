queryExport<-function(query_input,molecular_P){
  clinical_set<-query_input
  molecular_set<-molecular_P
  
  #MSI
  has_msi<-colnames(clinical_set)[grep(".*msi.*|.*mantis.*|molecular.*subtype|^Subtype",colnames(clinical_set),ignore.case = TRUE)]
  if (length(has_msi)==0){
    print("no MSI data")
    MSIvec<-rep(NA,nrow(clinical_set))
    clinical_set<-cbind(clinical_set,MSIvec)
  }
  else if ("MSI_STATUS"%in%has_msi){
    MSIvec<-clinical_set$MSI_STATUS
    MSIvec[which(MSIvec=="MSI-high"|MSIvec=="MSI-H"|MSIvec=="MSI")]<-"MSI"
    MSIvec[which(MSIvec=="MSI-stable"|MSIvec=="MSI-S"|MSIvec=="MSS")]<-"MSS"
    MSIvec[which(MSIvec!="MSI-stable"&MSIvec!="MSI-S"&MSIvec!="MSS"&MSIvec!="MSI-high"&MSIvec!="MSI-H"&MSIvec!="MSI")]<-NA
    
  }
  else if ("MSI_TYPE"%in%has_msi) {
    MSIvec<-clinical_set$MSI_TYPE
    MSIvec[which(MSIvec=="Instable")]<-"MSI"
    MSIvec[which(MSIvec=="Stable")]<-"MSS"
    MSIvec[which(MSIvec!="Stable"&MSIvec!="Instable")]<-NA
  }
  else if ("MOLECULAR_SUBTYPE"%in%has_msi) {
    MSIvec<-clinical_set$MOLECULAR_SUBTYPE
    MSIvec[which(MSIvec=="MSI"|MSIvec=="MSI-H")]<-"MSI"
    MSIvec[which(MSIvec=="MSS"|MSIvec=="MSI-S")]<-"MSS"
    MSIvec[which(MSIvec!="MSI"&MSIvec!="MSI-H"&MSIvec!="MSS"&MSIvec!="MSI-S")]<-NA
  }
  else if ("SUBTYPE"%in%has_msi) {
    SUBTYPE<-clinical_set$SUBTYPE
    MANTIS<-clinical_set$MSI_SCORE_MANTIS
    SENSOR<-clinical_set$MSI_SENSOR_SCORE
    MSIvec<-rep(NA,length(SUBTYPE))
    MSIvec[which(SUBTYPE=="COAD_MSI"|SUBTYPE=="READ_MSI"|MANTIS>=0.4|SENSOR>=3.5)]<-"MSI"
    MSIvec[which(MANTIS<0.4&SENSOR<3.5)]<-"MSS"
  }else{print("no MSI data")
    MSIvec<-rep(NA,nrow(clinical_set))
    clinical_set<-cbind(clinical_set,MSIvec)}
  print(unique(MSIvec))
  ground_truth<-read_csv("~/Analysis/2024/1. Download cBioportal data/allmuts2_crc_new.csv")
  muts2<-mutationData(cbio,molecularProfileIds = molecular_set[which(molecular_set$molecularAlterationType=="MUTATION_EXTENDED"),]$molecularProfileId, clinical_set$sampleId,entrezGeneIds =c(4292,4436,5395,4072,2956))#Add msh6
  
  
  #KRAS-BRAF
  muts<-mutationData(cbio,molecularProfileIds = molecular_set[which(molecular_set$molecularAlterationType=="MUTATION_EXTENDED"),]$molecularProfileId, clinical_set$sampleId,entrezGeneIds =c(3845,673))
  mutsx<-mutationData(cbio,molecularProfileIds = molecular_set[which(molecular_set$molecularAlterationType=="MUTATION_EXTENDED"),]$molecularProfileId, clinical_set$sampleId,entrezGeneIds = 1:200000)
  
  muts1a<-muts[[1]]
  muts2a<-muts2[[1]]
  
  if(length(muts[[1]])==0){
    cat("no KRAS or BRAF mutation data in",inf2$name[study])
    KBvec<-rep(NA,nrow(clinical_set))
    KRASvec<-rep(NA,nrow(clinical_set))
    BRAFvec<-rep(NA,nrow(clinical_set))
    KB_var_vec<-rep(NA,nrow(clinical_set))
    KB_mut_vec<-rep(NA,nrow(clinical_set))
    KRAS_var_vec<-rep(NA,nrow(clinical_set))
    KRAS_mut_vec<-rep(NA,nrow(clinical_set))
    BRAF_var_vec<-rep(NA,nrow(clinical_set))
    BRAF_mut_vec<-rep(NA,nrow(clinical_set))
    
  }
  else if (length(muts[[1]])>0){
    mutstibblex<-mutsx[[1]]
    genex<-tstrsplit(mutstibblex$keyword," ")[[1]]
    length(unique(genex))
    pt_with_any_variant<-unique(mutstibblex$patientId)
    
    mutstibble<-muts[[1]]
    m2<-mutstibble[grep("BRAF V600.*|KRAS.*",mutstibble$keyword),]
    m2_KRAS<-mutstibble[grep("KRAS.*",mutstibble$keyword),]
    m2_BRAF<-mutstibble[grep("BRAF V600.*",mutstibble$keyword),]
    
    pt_with_molecular_data<-unique(mutstibble$patientId)
    pt_KRAS_BRAF_pos<-unique(m2$patientId)
    pt_KRAS_pos<-unique(m2_KRAS$patientId)
    pt_BRAF_pos<-unique(m2_BRAF$patientId)
    
    m3<-mutstibblex[grep("BRAF V600.*|KRAS.*",mutstibblex$keyword,invert=TRUE),]
    pt_KRAS_BRAF_neg<-unique(m3$patientId)[which(!unique(m3$patientId)%in%pt_KRAS_BRAF_pos)]
    
    m3_KRAS<-mutstibblex[grep("KRAS.*",mutstibblex$keyword,invert=TRUE),]
    pt_KRAS_neg<-unique(m3_KRAS$patientId)[which(!unique(m3_KRAS$patientId)%in%pt_KRAS_pos)]
    
    m3_BRAF<-mutstibblex[grep("BRAF V600.*",mutstibblex$keyword,invert=TRUE),]
    pt_BRAF_neg<-unique(m3_BRAF$patientId)[which(!unique(m3_BRAF$patientId)%in%pt_BRAF_pos)]
    
    KBvec<-rep(NA,nrow(clinical_set))
    KBvec[which(clinical_set$patientId%in%pt_KRAS_BRAF_pos)]<-"variant_positive"
    KBvec[which(clinical_set$patientId%in%pt_KRAS_BRAF_neg)]<-"variant_negative"
    
    KRASvec<-rep(NA,nrow(clinical_set))
    KRASvec[which(clinical_set$patientId%in%pt_KRAS_pos)]<-"KRAS_variant_positive"
    KRASvec[which(clinical_set$patientId%in%pt_KRAS_neg)]<-"KRAS_variant_negative"
    
    BRAFvec<-rep(NA,nrow(clinical_set))
    BRAFvec[which(clinical_set$patientId%in%pt_BRAF_pos)]<-"BRAF_V600E_variant_positive"
    BRAFvec[which(clinical_set$patientId%in%pt_BRAF_neg)]<-"BRAF_V600E_variant_negative"
    
  }
  #MMR-mut
  
  
  if(length(muts2[[1]])==0|!"keyword"%in%colnames(muts2a)){
    cat("no MMR mutation data in",inf2$name[study])
    MMR_var_vec<-rep(NA,nrow(clinical_set))
    MMR_mut_vec<-rep(NA,nrow(clinical_set))
    
    KB_var_vec<-rep(NA,nrow(clinical_set))
    KB_mut_vec<-rep(NA,nrow(clinical_set))
    
    MLH1_var_vec<-rep(NA,nrow(clinical_set))
    MSH2_var_vec<-rep(NA,nrow(clinical_set))
    PMS2_var_vec<-rep(NA,nrow(clinical_set))
    EPCAM_var_vec<-rep(NA,nrow(clinical_set))
    MSH6_var_vec<-rep(NA,nrow(clinical_set))
    
    MLH1_mut_vec<-rep(NA,nrow(clinical_set))
    MSH2_mut_vec<-rep(NA,nrow(clinical_set))
    PMS2_mut_vec<-rep(NA,nrow(clinical_set))
    EPCAM_mut_vec<-rep(NA,nrow(clinical_set))
    MSH6_mut_vec<-rep(NA,nrow(clinical_set))
    
    KRAS_mut_vec<-rep(NA,nrow(clinical_set))
    BRAF_mut_vec<-rep(NA,nrow(clinical_set))
    
    KRASvec<-rep(NA,nrow(clinical_set))
    BRAFvec<-rep(NA,nrow(clinical_set))
    
    KRAS_var_vec<-rep(NA,nrow(clinical_set))
    BRAF_var_vec<-rep(NA,nrow(clinical_set))
    
    KRAS_mut_vec<-rep(NA,nrow(clinical_set))
    BRAF_mut_vec<-rep(NA,nrow(clinical_set))
    
    
    
  }
  else if (length(muts2[[1]])>0&"keyword"%in%colnames(muts2a)){
    mutstibblex<-mutsx[[1]]
    genex<-tstrsplit(mutstibblex$keyword," ")[[1]]
    length(unique(genex))
    pt_with_any_variant<-unique(mutstibblex$patientId)
    
    #variants in MMR
    muts2tibble<-muts2[[1]]
    m2.2<-muts2tibble[grep("MLH1.*|MSH2.*|PMS2.*|EPCAM.*|MSH6.*",muts2tibble$keyword),]
    m2.2b<-muts2tibble[grep("MLH1.*|MSH2.*|PMS2.*|EPCAM.*|MSH6.*|KRAS.*|BRAF.*",muts2tibble$keyword),]
    if (length(muts[[1]])>0){
    m2.2c<-mutstibble[grep("KRAS.*|BRAF.*",mutstibble$keyword),]
    }
    
    gene<-tstrsplit(m2.2$keyword," ")[[1]]
    gene2<-paste("GI=",gene,sep="")
    forvcf<-m2.2[,c("chr","startPosition","referenceAllele","variantAllele")]
    nrowvcf<-nrow(m2.2)
    m2.2$newID<-paste(m2.2$chr, m2.2$startPosition,m2.2$referenceAllele,m2.2$variantAllele,sep="_")

    #repeat for kras and braf
    if (length(muts[[1]])>0){
      gene<-tstrsplit(m2.2c$keyword," ")[[1]]
      gene2<-paste("GI=",gene,sep="")
      forvcf<-m2.2c[,c("chr","startPosition","referenceAllele","variantAllele")]
     nrowvcf<-nrow(m2.2c)
      m2.2c$newID<-paste(m2.2c$chr, m2.2c$startPosition,m2.2c$referenceAllele,m2.2c$variantAllele,sep="_")
  }
    #get MMR mutations
    pt_with_molecular_data2<-unique(muts2tibble$patientId)
    pt_MMR_var_pos<-unique(m2.2$patientId)
    if(length(muts[[1]])>0){
      pt_KB_var_pos<-unique(m2.2c$patientId)}
    else{
      pt_KB_var_pos<-NA
    }
    
    m3.2<-mutstibblex[grep("MLH1.*|MSH2.*|PMS2.*|EPCAM.*|MSH6.*",mutstibblex$keyword,invert = TRUE),]
    pt_MMR_var_neg<-unique(m3.2$patientId)[which(!unique(m3.2$patientId)%in%pt_MMR_var_pos)]
    MMR_var_vec<-rep(NA,nrow(clinical_set))
    MMR_var_vec[which(clinical_set$patientId%in%pt_MMR_var_pos)]<-"MMR_variant_positive"
    MMR_var_vec[which(clinical_set$patientId%in%pt_MMR_var_neg)]<-"MMR_variant_negative"
    
    m3.2c<-mutstibblex[grep("KRAS.*|BRAF.8",mutstibblex$keyword,invert = TRUE),]
    pt_KB_var_neg<-unique(m3.2c$patientId)[which(!unique(m3.2c$patientId)%in%pt_KB_var_pos)]
    KB_var_vec<-rep(NA,nrow(clinical_set))
    KB_var_vec[which(clinical_set$patientId%in%pt_KB_var_pos)]<-"KB_variant_positive"
    KB_var_vec[which(clinical_set$patientId%in%pt_KB_var_neg)]<-"KB_variant_negative"
    
    ##MMR mutations
    m2.3<-m2.2[m2.2$newID%in%ground_truth$newID,]
    pt_MMR_mut_pos<-unique(m2.3$patientId)
    #pt_MMR_mut_neg<-unique(m3.2$patientId)[which(!unique(m3.2$patientId)%in%pt_MMR_var_pos)]
    MMR_mut_vec<-rep(NA,nrow(clinical_set))
    MMR_mut_vec[which(clinical_set$patientId%in%pt_MMR_mut_pos)]<-"MMR_mutation_positive"
    MMR_mut_vec[which((!clinical_set$patientId%in%pt_MMR_mut_pos)&clinical_set$patientId%in%pt_with_any_variant)]<-"MMR_mutation_negative"
    
    #KRAS BRAF mutations
    if(length(muts[[1]])>0){
    m2.3c<-m2.2c[m2.2c$newID%in%ground_truth$newID,]
    pt_KB_mut_pos<-unique(m2.3c$patientId)
    KB_mut_vec<-rep(NA,nrow(clinical_set))
    KB_mut_vec[which(clinical_set$patientId%in%pt_KB_mut_pos)]<-"KB_mutation_positive"
    KB_mut_vec[which((!clinical_set$patientId%in%pt_KB_mut_pos)&clinical_set$patientId%in%pt_with_any_variant)]<-"KB_mutation_negative"
    }else{
      pt_KB_mut_pos<-NA
      KB_mut_vec<-rep(NA,nrow(clinical_set))
    }
    
    #Individual genes with variants
    m2.2_MLH1<-m2.2[grep("MLH1.*",m2.2$keyword),]
    m2.2_MSH2<-m2.2[grep("MSH2.*",m2.2$keyword),]
    m2.2_PMS2<-m2.2[grep("PMS2.*",m2.2$keyword),]
    m2.2_EPCAM<-m2.2[grep("EPCAM.*",m2.2$keyword),]
    m2.2_MSH6<-m2.2[grep("MSH6.*",m2.2$keyword),]
    
    if(length(muts[[1]])>0){
    m2.2c_KRAS<-m2.2c[grep("KRAS.*",m2.2c$keyword),]
    m2.2c_BRAF<-m2.2c[grep("BRAF.*",m2.2c$keyword),]
    }
    
    pt_MLH1_var_pos<-unique(m2.2_MLH1$patientId)
    pt_MSH2_var_pos<-unique(m2.2_MSH2$patientId)
    pt_PMS2_var_pos<-unique(m2.2_PMS2$patientId)
    pt_EPCAM_var_pos<-unique(m2.2_EPCAM$patientId)
    pt_MSH6_var_pos<-unique(m2.2_MSH6$patientId)
   
    if(length(muts[[1]])>0){
    pt_KRAS_var_pos<-unique(m2.2c_KRAS$patientId)
    pt_BRAF_var_pos<-unique(m2.2c_BRAF$patientId)
    }
    
    m3.2<-mutstibblex[grep("MLH1.*|MSH2.*|PMS2.*|EPCAM.*|MSH6.*",mutstibblex$keyword,invert=TRUE),]
    if(length(muts[[1]])>0){
    m3.2c<-mutstibblex[grep("KRAS.*|BRAF.*",mutstibblex$keyword,invert=TRUE),]
    }
    
    m3.2_MLH1<-mutstibblex[grep("MLH1.*",mutstibblex$keyword,invert=TRUE),]#all variant data that does not contain mlh1
    pt_MLH1_var_neg<-unique(m3.2_MLH1$patientId)[which(!unique(m3.2$patientId)%in%pt_MLH1_var_pos)]
    m3.2_MSH2<-mutstibblex[grep("MSH2.*",mutstibblex$keyword,invert=TRUE),]
    pt_MSH2_var_neg<-unique(m3.2_MSH2$patientId)[which(!unique(m3.2$patientId)%in%pt_MSH2_var_pos)]
    m3.2_PMS2<-mutstibblex[grep("PMS2.*",mutstibblex$keyword,invert=TRUE),]
    pt_PMS2_var_neg<-unique(m3.2_PMS2$patientId)[which(!unique(m3.2$patientId)%in%pt_PMS2_var_pos)]
    m3.2_EPCAM<-mutstibblex[grep("EPCAM.*",mutstibblex$keyword,invert=TRUE),]
    pt_EPCAM_var_neg<-unique(m3.2_EPCAM$patientId)[which(!unique(m3.2$patientId)%in%pt_EPCAM_var_pos)]
    m3.2_MSH6<-mutstibblex[grep("MSH6.*",mutstibblex$keyword,invert=TRUE),]
    pt_MSH6_var_neg<-unique(m3.2_MSH6$patientId)[which(!unique(m3.2$patientId)%in%pt_MSH6_var_pos)]
    
    if(length(muts[[1]])>0){
    m3.2c_KRAS<-mutstibblex[grep("KRAS.*",mutstibblex$keyword,invert=TRUE),]
    pt_KRAS_var_neg<-unique(m3.2c_KRAS$patientId)[which(!unique(m3.2c$patientId)%in%pt_KRAS_var_pos)]
    m3.2c_BRAF<-mutstibblex[grep("BRAF.*",mutstibblex$keyword,invert=TRUE),]
    pt_BRAF_var_neg<-unique(m3.2c_BRAF$patientId)[which(!unique(m3.2c$patientId)%in%pt_BRAF_var_pos)]
    }
    
    MLH1_var_vec<-rep(NA,nrow(clinical_set))
    MLH1_var_vec[which(clinical_set$patientId%in%pt_MLH1_var_pos)]<-"MLH1_variant_positive"
    MLH1_var_vec[which(clinical_set$patientId%in%pt_MLH1_var_neg)]<-"MLH1_variant_negative"
    
    MSH2_var_vec<-rep(NA,nrow(clinical_set))
    MSH2_var_vec[which(clinical_set$patientId%in%pt_MSH2_var_pos)]<-"MSH2_variant_positive"
    MSH2_var_vec[which(clinical_set$patientId%in%pt_MSH2_var_neg)]<-"MSH2_variant_negative"
    
    PMS2_var_vec<-rep(NA,nrow(clinical_set))
    PMS2_var_vec[which(clinical_set$patientId%in%pt_PMS2_var_pos)]<-"PMS2_variant_positive"
    PMS2_var_vec[which(clinical_set$patientId%in%pt_PMS2_var_neg)]<-"PMS2_variant_negative"
    
    EPCAM_var_vec<-rep(NA,nrow(clinical_set))
    EPCAM_var_vec[which(clinical_set$patientId%in%pt_EPCAM_var_pos)]<-"EPCAM_variant_positive"
    EPCAM_var_vec[which(clinical_set$patientId%in%pt_EPCAM_var_neg)]<-"EPCAM_variant_negative"
    
    MSH6_var_vec<-rep(NA,nrow(clinical_set))
    MSH6_var_vec[which(clinical_set$patientId%in%pt_MSH6_var_pos)]<-"MSH6_variant_positive"
    MSH6_var_vec[which(clinical_set$patientId%in%pt_MSH6_var_neg)]<-"MSH6_variant_negative"
    
    if(length(muts[[1]])>0){
    KRAS_var_vec<-rep(NA,nrow(clinical_set))
    KRAS_var_vec[which(clinical_set$patientId%in%pt_KRAS_var_pos)]<-"KRAS_variant_positive"
    KRAS_var_vec[which(clinical_set$patientId%in%pt_KRAS_var_neg)]<-"KRAS_variant_negative"
    
    BRAF_var_vec<-rep(NA,nrow(clinical_set))
    BRAF_var_vec[which(clinical_set$patientId%in%pt_BRAF_var_pos)]<-"BRAF_variant_positive"
    BRAF_var_vec[which(clinical_set$patientId%in%pt_BRAF_var_neg)]<-"BRAF_variant_negative"
    }else{
      KRAS_var_vec<-rep(NA,nrow(clinical_set))
      BRAF_var_vec<-rep(NA,nrow(clinical_set))
    }
    
    #Individual genes with mutations
    m2.3_MLH1<-m2.3[grep("MLH1.*",m2.3$keyword),]
    m2.3_MSH2<-m2.3[grep("MSH2.*",m2.3$keyword),]
    m2.3_PMS2<-m2.3[grep("PMS2.*",m2.3$keyword),]
    m2.3_EPCAM<-m2.3[grep("EPCAM.*",m2.3$keyword),]
    m2.3_MSH6<-m2.3[grep("MSH6.*",m2.3$keyword),]
    
    if(length(muts[[1]])>0){
    m2.3c_KRAS<-m2.3c[grep("KRAS.*",m2.3c$keyword),]
    m2.3c_BRAF<-m2.3c[grep("BRAF.*",m2.3c$keyword),]
    }
    
    pt_MLH1_mut_pos<-unique(m2.3_MLH1$patientId)
    pt_MSH2_mut_pos<-unique(m2.3_MSH2$patientId)
    pt_PMS2_mut_pos<-unique(m2.3_PMS2$patientId)
    pt_EPCAM_mut_pos<-unique(m2.3_EPCAM$patientId)
    pt_MSH6_mut_pos<-unique(m2.3_MSH6$patientId)
    
    if(length(muts[[1]])>0){
    pt_KRAS_mut_pos<-unique(m2.3c_KRAS$patientId)
    pt_BRAF_mut_pos<-unique(m2.3c_BRAF$patientId)
    }
    
    MLH1_mut_vec<-rep(NA,nrow(clinical_set))
    MLH1_mut_vec[which(clinical_set$patientId%in%pt_MLH1_mut_pos)]<-"MLH1_mutation_positive"
    MLH1_mut_vec[which((!clinical_set$patientId%in%pt_MLH1_mut_pos)&clinical_set$patientId%in%pt_with_any_variant)]<-"MLH1_mutation_negative" 
    
    MSH2_mut_vec<-rep(NA,nrow(clinical_set))
    MSH2_mut_vec[which(clinical_set$patientId%in%pt_MSH2_mut_pos)]<-"MSH2_mutation_positive"
    MSH2_mut_vec[which((!clinical_set$patientId%in%pt_MSH2_mut_pos)&clinical_set$patientId%in%pt_with_any_variant)]<-"MSH2_mutation_negative"
    
    PMS2_mut_vec<-rep(NA,nrow(clinical_set))
    PMS2_mut_vec[which(clinical_set$patientId%in%pt_PMS2_mut_pos)]<-"PMS2_mutation_positive"
    PMS2_mut_vec[which((!clinical_set$patientId%in%pt_PMS2_mut_pos)&clinical_set$patientId%in%pt_with_any_variant)]<-"PMS2_mutation_negative"
    
    EPCAM_mut_vec<-rep(NA,nrow(clinical_set))
    EPCAM_mut_vec[which(clinical_set$patientId%in%pt_EPCAM_mut_pos)]<-"EPCAM_mutation_positive"
    EPCAM_mut_vec[which((!clinical_set$patientId%in%pt_EPCAM_mut_pos)&clinical_set$patientId%in%pt_with_any_variant)]<-"EPCAM_mutation_negative"
    
    MSH6_mut_vec<-rep(NA,nrow(clinical_set))
    MSH6_mut_vec[which(clinical_set$patientId%in%pt_MSH6_mut_pos)]<-"MSH6_mutation_positive"
    MSH6_mut_vec[which((!clinical_set$patientId%in%pt_MSH6_mut_pos)&clinical_set$patientId%in%pt_with_any_variant)]<-"MSH6_mutation_negative"
    if(length(muts[[1]])>0){
    KRAS_mut_vec<-rep(NA,nrow(clinical_set))
    KRAS_mut_vec[which(clinical_set$patientId%in%pt_KRAS_mut_pos)]<-"KRAS_mutation_positive"
    KRAS_mut_vec[which((!clinical_set$patientId%in%pt_KRAS_mut_pos)&clinical_set$patientId%in%pt_with_any_variant)]<-"KRAS_mutation_negative"
    
    BRAF_mut_vec<-rep(NA,nrow(clinical_set))
    BRAF_mut_vec[which(clinical_set$patientId%in%pt_BRAF_mut_pos)]<-"BRAF_mutation_positive"
    BRAF_mut_vec[which((!clinical_set$patientId%in%pt_BRAF_mut_pos)&clinical_set$patientId%in%pt_with_any_variant)]<-"BRAF_mutation_negative"
    }else{
      KRAS_mut_vec<-rep(NA,nrow(clinical_set))
      BRAF_mut_vec<-rep(NA,nrow(clinical_set))
    }
  }
  augmented_data<-cbind(clinical_set,
                        MSIvec,
                        KBvec,KRASvec,BRAFvec,
                        KB_var_vec,KRAS_var_vec,BRAF_var_vec,
                        KB_mut_vec,KRAS_mut_vec,BRAF_mut_vec,
                        MMR_var_vec,MLH1_var_vec,MSH2_var_vec,PMS2_var_vec,EPCAM_var_vec,MSH6_var_vec,
                        MMR_mut_vec,MLH1_mut_vec,MSH2_mut_vec,PMS2_mut_vec,EPCAM_mut_vec,MSH6_mut_vec)
  return(augmented_data)
}
