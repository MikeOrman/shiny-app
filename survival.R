# Compute DFS benefit between primary subtype patients stratified by expression of the gene in the summary table
# Input = summary table from Genomic and Transcriptomic analysis
clinical.sample <- read.table("TCGA clinical sample.txt", header = TRUE, sep = "\t", check.names = FALSE)
clinical.patient <- read.table("TCGA clinical patient.txt", header = TRUE, sep = "\t", check.names = FALSE)
library(survival)
survival <- function(input){
  last.col <- ncol(input)+1
  # For each gene in the summary table:
  for (k in 1:nrow(input)){
    gene <- input$`Hugo Symbol`[k]
    # Obtain sample names for primary subtype
    sample.names <- colnames(TCGA.coloss)[2:ncol(TCGA.coloss)]
    # Create data frame of rows = samples | columns = Patient ID, Group, Time, Status
    data <- data.frame()
    for (i in 1:length(sample.names)){
      index <- which(sample.names[i] == clinical.sample$SAMPLE_ID)
      sample <- clinical.sample$PATIENT_ID[index]
      index2 <- which(sample == clinical.patient$PATIENT_ID)
      data[i,1] <- clinical.patient$PATIENT_ID[index2]
      if (clinical.patient$DFS_MONTHS[index2] != "[Not Available]") {data[i,2] <- clinical.patient$DFS_MONTHS[index2]}
      if (clinical.patient$DFS_STATUS[index2] != "[Not Available]") {data[i,3] <- clinical.patient$DFS_STATUS[index2]}
    }
    for (i in 1:length(sample.names)){
      index <- which(sample.names[i] == colnames(TCGA.mRNA))
      index2 <- which(gene == TCGA.mRNA$Hugo_Symbol)
      data[i,4] <- TCGA.mRNA[index2,index]
    }
    data <- na.omit(data)
    colnames(data) <- c("Patient ID", "DFS Months", "DFS Status", "Expression")
    # Rename status column
    for (i in 1:nrow(data)){
      if (data[i,3] == "0:DiseaseFree") {data[i,3] <- 0}
      if (data[i,3] == "1:Recurred/Progressed") {data[i,3] <- 1}
    }
    # Order dataset by gene expression. Subset into high and low expressors
    data <- data[order(data$Expression, decreasing = TRUE),]
    # Stratify patients into n=2 quantiles.
    quantiles <- quantile(data$Expression, 0.50)
    high.expressing <- data[data$Expression >= quantiles[1],]
    high.expressing[,5] <- "High Expressors"
    colnames(high.expressing)[5] <- "Group"
    low.expressing <- data[data$Expression <= quantiles[1],]
    low.expressing[,5] <- "Low Expressors"
    colnames(low.expressing)[5] <- "Group"
    # Prepare dataframe for survival analysis
    high.expressing <- high.expressing[,c(1, 5, 2, 3)]
    low.expressing <- low.expressing[,c(1, 5, 2, 3)]
    surv_data <- rbind(high.expressing, low.expressing)
    surv_data$`DFS Months` <- as.numeric(surv_data$`DFS Months`)
    surv_data$`DFS Status` <- as.numeric(surv_data$`DFS Status`)
    # Compute survival statistics
    surv_object <- Surv(time = surv_data$`DFS Months`, event = surv_data$`DFS Status`)
    test <- survdiff(surv_object ~ Group, data = surv_data, rho = 0)
    chisq <- test$chisq
    pval <- pchisq(chisq, length(test$n)-1, lower.tail = FALSE)
    input[k,last.col] <- pval
  }
  colnames(input)[last.col] <- "Survival pval"
  return(input)
}