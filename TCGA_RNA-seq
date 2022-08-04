library("TCGAbiolinks")
library("SummarizedExperiment")
library("data.table")
library("tidyverse")

#Download TCGA-PAAD RNA-seq data and transform the dataframe
GDCprojects<-getGDCprojects()
head(GDCprojects[c("project_id", "name")])
TCGAbiolinks:::getProjectSummary("TCGA-PAAD")
TCGA_PAAD_RNASeq_hiseq <- getLinkedOmicsData(
  project = "TCGA-PAAD",
  dataset = "RNAseq (HiSeq, Gene level)"
)
t(TCGA_PAAD_RNASeq_hiseq)->TCGA_PAAD_RNASeq_hiseq
TCGA_PAAD_RNASeq_hiseq <- as.data.frame(as.matrix(x=TCGA_PAAD_RNASeq_hiseq))
names(TCGA_PAAD_RNASeq_hiseq) <- as.matrix(TCGA_PAAD_RNASeq_hiseq[1, ])
TCGA_PAAD_RNASeq_hiseq <- TCGA_PAAD_RNASeq_hiseq[-1, ]
TCGA_PAAD_RNASeq_hiseq[] <- lapply(TCGA_PAAD_RNASeq_hiseq, function(x) type.convert(as.character(x)))
rownames_to_column(TCGA_PAAD_RNASeq_hiseq,var = "Sample ID") -> TCGA_PAAD_RNASeq_hiseq

#Extract patient info
query_TCGA = GDCquery(
  project = "TCGA-PAAD",
  data.category = "Transcriptome Profiling", 
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor"), data.type = "Gene Expression Quantification")

GDCdownload(query = query_TCGA)
tcga_data = GDCprepare(query_TCGA)

tcga_data@colData$patient
Patient_ID<-tcga_data@colData$patient
Patient_ID <- as.data.frame(as.matrix(x=Patient_ID))
colnames(Patient_ID) <- c("Sample ID")

tcga_data@colData$vital_status
Vital_status<-tcga_data@colData$vital_status

tcga_data$deceased = tcga_data$vital_status == "Dead"
tcga_data$deceased 

tcga_data$days_to_death
tcga_data$days_to_last_follow_up

tcga_data$overall_survival = ifelse(tcga_data$deceased,
                                    tcga_data$days_to_death,
                                    tcga_data$days_to_last_follow_up)

tcga_data$overall_survival
Overall_survival<-tcga_data@colData$overall_survival
Overall_survival <- as.data.frame(as.matrix(x=Overall_survival))
colnames(Overall_survival) <- c("Survival time")


tcga_data$deceased <- as.integer(tcga_data$deceased == "TRUE")
deceased<-tcga_data$deceased
deceased <- as.data.frame(as.matrix(x=deceased))
colnames(deceased) <- c("Survival event")

tcga_data@colData$'paper_Pathologist Reviewed Tumor Cellularity'
Tumor_cellularity<-tcga_data@colData$'paper_Pathologist Reviewed Tumor Cellularity'
Tumor_cellularity <- as.data.frame(as.matrix(x=Tumor_cellularity))
colnames(Tumor_cellularity) <- c("Filter_A")

combined <- cbind(Patient_ID, Overall_survival, deceased, Tumor_cellularity)
combined[order(combined$'Sample ID'),] ->combined

#Combine RNA-seq data and patient info
TCGA_PAAD_RNASeq_hiseq$`Sample ID` <- NULL
combined <- cbind(combined, TCGA_PAAD_RNASeq_hiseq)
combined %>% drop_na() -> combined

#Export data as .csv file
write.csv(combined,file="combined.csv", row.names = FALSE)
