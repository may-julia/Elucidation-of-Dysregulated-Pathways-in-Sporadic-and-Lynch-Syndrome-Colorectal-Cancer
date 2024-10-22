# 1. Install packages ####
if (!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("cBioPortalData")
install.packages("combinat")
install.packages("xlsx")

# 2. Load Libs ####
library(cBioPortalData)
library(xlsx)
library(combinat)
library(tidyr)
library(data.table)

# 3. query export function ####
source("~/Analysis/2024/1. Download cBioportal data/Query_export.R")

# 4. Initialise API ####
cbio <- cBioPortal()
cbio
cbiotags<-tags(cbio)
types<-cbio$getAllCancerTypesUsingGET()
searchOps(cbio,"Studies")
studies<-cbio$getAllStudiesUsingGET()

parsedResponse <- httr::content(studies)
cat("Number of elements in the response:", length(parsedResponse))

studies<-getStudies(cbio)
studies
cat("Answer 1: There are", nrow(studies), "studies in cBioPortal")
cat("Answer 2: The studies spans", length(unique(studies$cancerTypeId)), "cancer types")
cat("Answer 3: There are ", sum(studies$allSampleCount), "samples in cBioPortal")
cat("Answer 4: The study with the most samples is", studies[which.max(studies$allSampleCount), "name"][[1]])
ls("package:cBioPortalData")


# 5. Get studies (eg. CRC) ####
inf<-getStudies(cbio,buildReport = TRUE)
unique(inf$cancerTypeId)
inf[which(inf$cancerTypeId=="coadread"),]
crc1<-inf[grep("rect|colo[n,r]",inf$name,ignore.case = TRUE),]$name # change to include more cancer types 
crc2<-crc1[-c(intersect(grep("TCGA",crc1),grep("Nature|Firehose",crc1)))] # exclude the over lapping studies 
crc2 <- "Colorectal Adenocarcinoma (TCGA, PanCancer Atlas)" 

inf2<-inf[which(inf$name%in%crc2),]
inf2$description
inf2$studyId
inf2$name
inf2$pmid

# 6. Collect sample, clinical, molecular profile data for all studies in lists ####
sample_list<-list()
clinical_list<-list()
molprof_list<-list()
samlist<-list()

for (name in inf2$studyId){
  sample_list[[name]]<-allSamples(cbio,name)
  clinical_list[[name]]<-clinicalData(cbio,name) 
  molprof_list[[name]]<-molecularProfiles(cbio,name)
  samlist[[name]]<-getSampleInfo(cbio,name,projection = "DETAILED")
}

# 7. Export loop ####
for (study in 1:length(clinical_list)){
  name<-inf2$name[study]
  print(name)
  molprof<-molprof_list[[study]]  
  workdata<-clinical_list[[study]]
  aug_dat<-queryExport(workdata,molprof)
  aug_dat2 <- data.frame(lapply(aug_dat, function(x) {
    gsub(",", ";", x)
  }))
  
  write.xlsx2(aug_dat,file.path("~/Analysis/2024/1. Download cBioportal data",paste(name,"_clin_aug_01_11.xlsx")),quote = FALSE,sep="\t",row.names=FALSE)
  
}
