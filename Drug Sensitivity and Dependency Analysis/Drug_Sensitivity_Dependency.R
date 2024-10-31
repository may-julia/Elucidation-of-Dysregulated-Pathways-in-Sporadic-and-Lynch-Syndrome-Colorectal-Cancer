setwd("~/Analysis/2024/R_drugs")

library(dplyr)
library(data.table)
library(readr)
library(tidyr)
library(zoo)
library(purrr)
library(broom)
library(readxl)
library(grid)
library(ggplot2)

# Load the CCLE data ####
if (file.exists('CCLE_mrna_processed_20Q4.csv')) {
  cat('\n Loading mRNA transcription data for the cell lines \n')
  CCLEmrna <- fread('CCLE_mrna_processed_20Q4.csv')
  sampleInfo <- fread('Achilles_sample_info_20Q4.csv')
  pharos <- fread('pharos_PI3K.csv')
} else {
  cat('\n Processing the CCLE mRNA expression data \n')
  CCLEmrna <- fread('CCLE_expression.csv')
  
  #Clean up the data
  new_names <- gsub('\\(\\d+\\)', '', names(CCLEmrna))
  new_names <- trimws(new_names)
  unique_indices <- !duplicated(new_names)
  CCLEmrna <- CCLEmrna[, unique_indices, with = FALSE]
  setnames(CCLEmrna, new_names[unique_indices])
  
  #Filter for PI3K_AKT genes
  if (file.exists('pharos_PI3K.csv')) {
    cat('\n Loading Pharos data \n')
    pharos <- fread('pharos_PI3K.csv')
  } else {
    cat('\n Pharos data not found \n')
  }
  selected_columns <- c(1, which(names(CCLEmrna) %in% pharos$HugoSymbol))
  CCLEmrna <- CCLEmrna[, ..selected_columns, with = FALSE]
  
  #Change the cell line ID to those of the common names from the DepMap IDs
  if (file.exists('Achilles_sample_info_20Q4.csv')) {
    cat('\n Loading DepMap IDs data \n')
    sampleInfo <- fread('Achilles_sample_info_20Q4.csv')
  } else {
    cat('\n Sample info not found \n')
  }
  setnames(CCLEmrna, old = names(CCLEmrna)[1], new = names(sampleInfo)[1])
  CCLEmrna <- merge(sampleInfo[, 1:2, with = FALSE], CCLEmrna, by = names(sampleInfo)[1])
  CCLEmrna <- CCLEmrna[, -1, with = FALSE]
  setnames(CCLEmrna, old = names(CCLEmrna)[1], new = "cell_line")
  
  #Save the processed data to csv
  fwrite(CCLEmrna, 'CCLE_mrna_processed_20Q4.csv')
}

# Process the achilles data to that it matches the CCLE data  ####
if (file.exists('Achilles_gene_effect_CRISPR_processed_20Q4.csv')) {
  cat('\n Loading crispr transcription data for the cell lines \n')
  crispr <- fread('Achilles_gene_effect_CRISPR_processed_20Q4.csv')
} else {
  cat('\n Processing the crispr data \n')
  crispr <- fread('Achilles_gene_effect_CRISPR_20Q4.csv')
  
  #Clean up the data
  new_names <- gsub('\\(\\d+\\)', '', names(crispr))
  new_names <- trimws(new_names)
  unique_indices <- !duplicated(new_names)
  crispr <- crispr[, unique_indices, with = FALSE]
  setnames(crispr, new_names[unique_indices])
  
  #Filter for PI3K_AKT genes
  selected_columns <- c(1, which(names(crispr) %in% pharos$HugoSymbol))
  crispr <- crispr[, ..selected_columns, with = FALSE]
  
  #Change the cell line ID to those of the common names from the DepMap IDs
  setnames(crispr, old = names(crispr)[1], new = names(sampleInfo)[1])
  crispr <- merge(sampleInfo[, 1:2, with = FALSE], crispr, by = names(sampleInfo)[1])
  crispr <- crispr[, -1, with = FALSE]
  setnames(crispr, old = names(crispr)[1], new = "cell_line")
  
  #Save the processed data to csv
  fwrite(crispr, 'Achilles_gene_effect_CRISPR_processed_20Q4.csv')
}

# Get the clustering data and the missing data and fill the missing data 
setnames(sampleInfo, old = names(sampleInfo)[2], new = "cell_line")
setnames(sampleInfo, old = names(sampleInfo)[18], new = "disease")
achillesCCLEmatched <- inner_join(sampleInfo %>% select(disease, cell_line), crispr, by = "cell_line")
# Fill missing values using linear interpolation for specific columns
achillesCCLEmatched <- achillesCCLEmatched %>%
  mutate(across(all_of(names(achillesCCLEmatched)[3:ncol(achillesCCLEmatched)]), ~ na.approx(.x, na.rm = FALSE)))

mostEssential <- achillesCCLEmatched %>%
  select(-c(disease, cell_line)) %>%
  mutate(across(everything(), ~ . < -0.50))

mostEssentialSummary <- mostEssential %>%
  summarise(across(everything(), sum)) %>%
  pivot_longer(cols = everything(), names_to = "HugoSymbol", values_to = "dependentCellLines")

# Perform an inner join with the pharos table
setnames(pharos, old = names(pharos)[3], new = "HugoSymbol")
pharos_selected <- pharos %>% select(2, 3)
mostEssential <- inner_join(pharos_selected, mostEssentialSummary, by = "HugoSymbol")

# Sort the results based on the number of dependent cell lines
mostEssential <- mostEssential %>%
  arrange(desc(dependentCellLines))

# load and process the CCLE mutation data #####
if (file.exists('CCLE_mutations_processed_20Q4.RData')) {
  cat('\n Loading processed mutations data \n')
  CCLEmutations <- fread('CCLE_mutations_processed_20Q4.csv')
} else {
  cat('\n Processing the CCLE mutations data \n')
  CCLEmutations <- fread('CCLE_mutations_20Q4.csv')
  
  #Filter for PI3K_AKT genes
  CCLEmutations <- CCLEmutations %>%
    filter(Hugo_Symbol %in% pharos$HugoSymbol)
  
  # Add missing Tumor_Sample_Barcode column if necessary
  if (!'Tumor_Sample_Barcode' %in% colnames(CCLEmutations)) {
    CCLEmutations <- CCLEmutations %>% 
      mutate(Tumor_Sample_Barcode = DepMap_ID)
  }
  
  source("~/Analysis/2024/R_drugs/processMAF.R")
  CCLEmutations <- processMAF(CCLEmutations)
  
  setnames(CCLEmutations, old = names(CCLEmutations)[1], new = "DepMap_ID")
  CCLEmutations <- CCLEmutations %>%
    inner_join(sampleInfo %>% select(DepMap_ID, cell_line), by = 'DepMap_ID') %>%
    select(DepMap_ID, cell_line, everything())

  CCLEmutations <- CCLEmutations %>% 
    filter(cell_line %in% crispr$cell_line)
  
  # make the mutation data the sample length as the crispr data get
  # the missing genes and add them to the mutation data using NaN
  existing_genes <- colnames(CCLEmutations)
  missing_genes <- setdiff(colnames(crispr), existing_genes)
  for (gene in missing_genes) {
    CCLEmutations[[gene]] <- NA
  }
  
  # Ensure the mutation data has the same rows as the CRISPR data
  CCLEmutations <- CCLEmutations %>%
    right_join(crispr %>% select(cell_line), by = "cell_line")
  
  #Save the processed data to csv
  fwrite(CCLEmutations, 'CCLE_mrna_processed_20Q4.RData')
}


# Find the correlation between mRNA to CRISPR gene in the data ########

# Check if the .RData file already exists
if (file.exists("crispr_mrna_corr_results.RData")) {
  # If it does, load the file
  load("crispr_mrna_corr_results.RData")
} else {
  # Merge tables on 'cell_line'
  mergedData <- inner_join(achillesCCLEmatched, CCLEmrna, by = "cell_line")
  names(mergedData) <- sub("\\.x$", "_achillesCCLEmatched", names(mergedData))
  names(mergedData) <- sub("\\.y$", "_CCLEmrna", names(mergedData))

  # Convert to factor
  if (is.character(mergedData$disease)) {
    mergedData$disease <- as.factor(mergedData$disease)
  }
  
  # Get the list of genes (removing the _CCLEmrna suffix)
  genes <- unique(c(names(achillesCCLEmatched)[3:ncol(achillesCCLEmatched)], names(CCLEmrna)[2:ncol(CCLEmrna)]))
  
  # Get the list of cancer types
  cancerTypes <- levels(mergedData$disease)
  
  # Specify the number of rows
  num_rows <- length(genes) * length(cancerTypes)
  
  # Initialize row counter and results storage
  row <- 1
  corr_results <- data.frame(CancerType = character(), Gene = character(), R2 = numeric(), P_Value = numeric(), stringsAsFactors = FALSE)
  
  # Loop through each cancer type
  for (jj in seq_along(cancerTypes)) {
    cat(sprintf("\nRunning analysis for cancer study number %d of %d\n", jj, length(cancerTypes)))
    
    # Get the current cancer type
    cancerType <- cancerTypes[jj]
    
    # Loop through each gene
    for (ii in seq_along(genes)) {
      gene <- genes[ii]
      
      if (ii %% 200 == 0) {
        cat(sprintf("\nRunning analysis for gene number %d of %d, in Cancer Study %d of %d\n", ii, length(genes), jj, length(cancerTypes)))
      }
      
      # Subset the data for the current cancer type
      subData <- mergedData %>% filter(disease == cancerType)
      
      if (nrow(subData) > 0) {
        # Compute the correlation
        tryCatch({
          gene_mrna <- paste0(gene, "_CCLEmrna")
          gene_achilles <- paste0(gene, "_achillesCCLEmatched")
          
          if (gene_mrna %in% names(subData) && gene_achilles %in% names(subData)) {
            cor_test <- cor.test(subData[[gene_mrna]], subData[[gene_achilles]], method = "pearson", use = "complete.obs")
            
            # Store the result
            corr_results[row, ] <- c(cancerType, gene, cor_test$estimate^2, cor_test$p.value)
            row <- row + 1
          } else {
            cat(sprintf("Skipping gene %s for cancer type %s: data not found.\n", gene, cancerType))
          }
        }, error = function(e) {
          # Handle the error gracefully
          cat(sprintf("Error processing gene %s for cancer type %s: %s\n", gene, cancerType, e$message))
        })
      } else {
        cat(sprintf("No data for cancer type %s.\n", cancerType))
      }
    }
  }
  
  # Remove unused rows
  if (row > 1) {
    corr_results <- corr_results[1:(row - 1), ]
  } else {
    cat("No valid correlations found.\n")
  }
  
  # Convert columns to appropriate types
  corr_results <- as.data.frame(corr_results)
  corr_results$R2 <- as.numeric(as.character(corr_results$R2))
  corr_results$P_Value <- as.numeric(as.character(corr_results$P_Value))
  
  # Sort the table by p-value
  corr_results <- corr_results %>% arrange(P_Value)
  
  # Save the table to an .RData file
  save(corr_results, file = "crispr_mrna_corr_results.RData")
  
}

# Relatioship Between Developmental Stage and Drug Response 
# load the GDSC drug response data 
gdsc1 <- read_excel("GDSC1_fitted_dose_response_27Oct23.xlsx")
gdsc2 <- read_excel("GDSC2_fitted_dose_response_27Oct23.xlsx")
gdscDoseResponse <- rbind(gdsc1, gdsc2)

# Convert PATHWAY_NAME to a factor 
gdscDoseResponse$PATHWAY_NAME <- as.factor(gdscDoseResponse$PATHWAY_NAME)

# Filter rows where PATHWAY_NAME includes "PI3K"
gdscDoseResponse <- gdscDoseResponse %>% filter(grepl("PI3K", as.character(PATHWAY_NAME)))

# Remove dashes from CELL_LINE_NAME
gdscDoseResponse$CELL_LINE_NAME <- gsub("-", "", gdscDoseResponse$CELL_LINE_NAME)
crispr$cell_line <- gsub("-", "", crispr$cell_line)

# Rename the 5th column to "cell_line"
names(gdscDoseResponse)[5] <- "cell_line"

# Return only the GDSC cell lines that have CRISPR data
gdscDoseResponse <- gdscDoseResponse %>% filter(cell_line %in% crispr$cell_line)

# Convert PATHWAY_NAME to a factor 
gdscDoseResponse$PATHWAY_NAME <- as.factor(as.character(gdscDoseResponse$PATHWAY_NAME))

# Compare the Response of the Cell Lines for the Achille Dataset ####
# are cell line with higher dependence score on a particular gene significantly more responpsive to a particular drug 

library(dplyr)
library(readxl)
library(parallel)

# Rename the column TargetDevelopmentLevel to DevStage
colnames(pharos)[which(colnames(pharos) == "Target Development Level")] <- "DevStage"

library(readxl)
library(dplyr)
library(parallel)
library(openxlsx)

# Initialize an empty data frame for results
crisprGdscResults <- data.frame()

# Function to perform the analysis for each drug
perform_analysis <- function(curDrug, validCrispr, crispr, gdscDoseResponse) {
  appendTable <- data.frame(DRUG_NAME = character(), HugoSymbol = character(), meanHigherRank = numeric(), 
                            meanLowerRank = numeric(), tStat = numeric(), lowerBound = numeric(), 
                            upperBound = numeric(), pValue = numeric(), numHighRank = numeric(), 
                            numLowRank = numeric(), stringsAsFactors = FALSE)
  
  for (jj in 1:ncol(validCrispr)) {
    geneName <- colnames(validCrispr)[jj]
    
    if (jj %% 200 == 0) {
      cat(sprintf("Running analysis for drug %s and gene #%d: %s\n", curDrug, jj, geneName))
    }
    
    curDepCellLines <- crispr$cell_line[validCrispr[[jj]] < -0.50]
    
    if (length(curDepCellLines) < 6 || length(curDepCellLines) > nrow(validCrispr) - 6) {
      next
    }
    
    curGDSC <- gdscDoseResponse %>% filter(DRUG_NAME == curDrug)
    locInGDSC <- curGDSC$cell_line %in% curDepCellLines
    
    numInstances <- table(locInGDSC)
    if (any(numInstances < 6)) {
      next
    }
    
    if (sum(locInGDSC) < 6 || sum(!locInGDSC) < 6) {
      next
    }
    
    ttest_result <- t.test(
      curGDSC$Z_SCORE[locInGDSC],
      curGDSC$Z_SCORE[!locInGDSC],
      var.equal = FALSE
    )
    
    meanHighRank <- mean(curGDSC$Z_SCORE[locInGDSC])
    meanLowerRank <- mean(curGDSC$Z_SCORE[!locInGDSC])
    
    appendTable <- rbind(appendTable, data.frame(DRUG_NAME = curDrug, HugoSymbol = geneName, 
                                                 meanHigherRank = meanHighRank, meanLowerRank = meanLowerRank, 
                                                 tStat = ttest_result$statistic, lowerBound = ttest_result$conf.int[1], 
                                                 upperBound = ttest_result$conf.int[2], pValue = ttest_result$p.value, 
                                                 numHighRank = sum(locInGDSC), numLowRank = sum(!locInGDSC)))
  }
  
  cat("Structure of appendTable for drug:", curDrug, "\n")
  str(appendTable)
  
  cat(sprintf("Completed analysis for drug: %s\n", curDrug))
  return(appendTable)
}


tryCatch({
  crisprGdscResults <- read_excel('Supplementary_File_4.xlsx', sheet = 'Between Cell Line Response')
  cat("Loaded existing data from 'Supplementary File 4.xlsx'\n")
}, error = function(e) {
  cat("Error loading file or file does not exist. Proceeding with analysis...\n")
  
  validCrispr <- crispr %>% select(-cell_line)
  lowerBound <- colSums(validCrispr < -0.5) >= 10
  upperBound <- colSums(validCrispr < -0.5) <= (nrow(validCrispr) - 10)
  validGenes <- lowerBound & upperBound
  validCrispr <- validCrispr %>% select(which(validGenes))
  
  gdscDoseResponse$DRUG_NAME <- as.factor(gdscDoseResponse$DRUG_NAME)
  gdscDrugs <- as.character(unique(gdscDoseResponse$DRUG_NAME))
  
  cat("Number of drugs to analyze: ", length(gdscDrugs), "\n")
  
  results_list <- lapply(gdscDrugs, function(curDrug) {
    cat("Analyzing drug: ", curDrug, "\n")
    result <- perform_analysis(curDrug, validCrispr, crispr, gdscDoseResponse)
    if (nrow(result) > 0) {
      cat("Result for drug", curDrug, ":\n")
      print(head(result))
    }
    return(result)
  })
  
  crisprGdscResults <- do.call(rbind, results_list)
  
  if (nrow(crisprGdscResults) == 0) {
    stop("Error: The results from the computation are empty.")
  }
  
  crisprGdscResults <- as.data.frame(crisprGdscResults)
  
  cat("Structure of crisprGdscResults before filtering:\n")
  str(crisprGdscResults)
  
  crisprGdscResults <- crisprGdscResults %>% 
    filter(pValue < 0.05) %>%
    mutate(FDR = p.adjust(pValue, method = "BH"))
  
  if (nrow(crisprGdscResults) == 0) {
    stop("Error: No significant results after filtering by pValue < 0.05.")
  }
  
  cat("Intermediate crisprGdscResults after filtering:\n")
  print(head(crisprGdscResults))
  
  # Join with pharos to add DevStage
  crisprGdscResults <- inner_join(crisprGdscResults, pharos %>% select(HugoSymbol, DevStage), by = "HugoSymbol")
  
  # Get unique DRUG_NAME, PUTATIVE_TARGET, and PATHWAY_NAME from gdscDoseResponse
  gdscTargets <- gdscDoseResponse %>% 
    select(DRUG_NAME, PUTATIVE_TARGET, PATHWAY_NAME) %>% 
    distinct()
  
  # Join with crisprGdscResults to add PUTATIVE_TARGET and PATHWAY_NAME
  crisprGdscResults <- inner_join(crisprGdscResults, gdscTargets, by = "DRUG_NAME")
  
  # Reorder columns to move PUTATIVE_TARGET and PATHWAY_NAME after DRUG_NAME
  crisprGdscResults <- crisprGdscResults %>%
    select(DRUG_NAME, PUTATIVE_TARGET, PATHWAY_NAME, everything())
  
  # Sort the results by FDR in ascending order
  crisprGdscResults <- crisprGdscResults %>%
    arrange(FDR)
  
  write.xlsx(crisprGdscResults, "Supplementary_File_4.xlsx", sheetName = "Between Cell Line Response")
  
  print("Analysis and saving complete.")
})

cat("Final structure of crisprGdscResults:\n")
str(crisprGdscResults)

# Create a table with genes along the column and drug along the rows #####
# Convert the drug name to a factor
crisprGdscResults$DRUG_NAME <- as.factor(crisprGdscResults$DRUG_NAME)

# Create the preallocated table
unique_drugs <- unique(crisprGdscResults$DRUG_NAME)
unique_genes <- unique(crisprGdscResults$HugoSymbol)
drugTable <- matrix(0, nrow = length(unique_genes), ncol = length(unique_drugs),
                    dimnames = list(unique_genes, as.character(unique_drugs)))
drugTable <- as.data.frame(drugTable)

# Add the gene name to the table
drugTable$HugoSymbol <- rownames(drugTable)
drugTable <- drugTable[, c("HugoSymbol", setdiff(names(drugTable), "HugoSymbol"))]


# Add the t-values to the table
for (ii in 2:ncol(drugTable)) {
  curDrug <- colnames(drugTable)[ii]
  curData <- subset(crisprGdscResults, DRUG_NAME == curDrug)
  gene_indices <- match(curData$HugoSymbol, drugTable$HugoSymbol)
  drugTable[gene_indices, curDrug] <- curData$tStat
}

# Remove rows with all missing data 
drugTable <- drugTable[rowSums(drugTable[, -1] == 0) != (ncol(drugTable) - 1), ]
drugTable <- drugTable[, colSums(drugTable == 0) != nrow(drugTable)]

# Display the first 10 rows and columns
print(drugTable[1:10, 1:10])

# Produce a heatmap with correlation coefficient data
library(gplots)
heatmap.2(as.matrix(drugTable[, -1]), trace = "none", col = bluered(256), 
          dendrogram = "both", scale = "row", margins = c(5, 10),
          Colv = TRUE, Rowv = TRUE, labCol = colnames(drugTable)[-1], 
          labRow = drugTable$HugoSymbol)

# Perform a Chi-Square test for Crispr Response  ######
# Add the developmental stage to the data
if (!"DevStage" %in% names(crisprGdscResults)) {
  crisprGdscResults <- merge(pharos[, c("HugoSymbol", "DevStage")], crisprGdscResults, by = "HugoSymbol")
}

# Convert DevStage to a factor (categorical in R)
if (is.character(crisprGdscResults$DevStage)) {
  crisprGdscResults$DevStage <- as.factor(crisprGdscResults$DevStage)
}

# Sort the table by the FDR column in ascending order
crisprGdscResults <- crisprGdscResults[order(crisprGdscResults$FDR), ]

# Get the most significant results based on numHighRank and numLowRank
locTop <- crisprGdscResults$numHighRank > 20 & crisprGdscResults$numLowRank > 20

# Return only those instances
top9Results <- crisprGdscResults[locTop, ]

# top 9 results to be plotted
top9Results <- top9Results[1:9, ]

# Define colors for CRISPR (switched)
crisprColor <- c("#FC4E2A", "#2B83BA")  # Red for "Lower", Blue for "Higher"

# Function to create the scatter and box plots
colourBoxPlot <- function(plotData, groups, colors, includeScatter, showYAxisLabel) {
  group_labels <- ifelse(groups, "Higher", "Lower")
  p <- ggplot(data.frame(plotData, group_labels), aes(x = factor(group_labels), y = plotData)) +
    geom_boxplot(aes(fill = factor(group_labels)), width = 0.3, outlier.shape = 19, outlier.size = 1.5, outlier.colour = "black") +
    scale_fill_manual(values = colors) +
    scale_colour_manual(values = colors) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "bold", size = 8),
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      legend.position = "none"
    ) +
    labs(
      x = "Dependency"
    )
  if (showYAxisLabel) {
    p <- p + labs(y = "IC_50 Z-score")
  } else {
    p <- p + theme(axis.title.y = element_blank())
  }
  if (includeScatter) {
    p <- p + geom_jitter(width = 0.2, aes(colour = factor(group_labels)), alpha = 0.5, size = 1)
  }
  return(p)
}

# Letters for annotating plots
theLetters <- letters[1:9]

# Create a list to store the plots
plots <- list()

# Loop through each row of top9Results to plot individual panels
for (ii in 1:nrow(top9Results)) {
  # Get the GDSC data for the current drug
  curGDSC <- gdscDoseResponse[gdscDoseResponse$DRUG_NAME == top9Results$DRUG_NAME[ii], ]
  
  # Get the cell lines highly dependent on the current genes
  curDepCellLines <- crispr$cell_line[crispr[[top9Results$HugoSymbol[ii]]] < -0.50]
  
  # Get the location of those cell lines in the current GDSC data
  locInGDSC <- curGDSC$cell_line %in% curDepCellLines
  
  # Fill outliers in Z_SCORE with nearest mean
  dReponses <- curGDSC$Z_SCORE
  dReponses[is.na(dReponses)] <- ave(dReponses, FUN = function(x) mean(x, na.rm = TRUE))
  
  # Determine whether to show y-axis label
  showYAxisLabel <- ii %% 3 == 1
  
  # Plot the box plots
  includeScatter <- TRUE
  p <- colourBoxPlot(dReponses, locInGDSC, crisprColor, includeScatter, showYAxisLabel) +
    labs(title = paste(top9Results$DRUG_NAME[ii], ":", top9Results$HugoSymbol[ii])) 
  
  plots[[ii]] <- p
}

# Function to add annotations to the grid
add_annotations <- function(grid_plots, annotations, ncol, nrow) {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow, ncol, heights = unit(rep(1, nrow), "null"), widths = unit(rep(1, ncol), "null"))))
  for (i in seq_along(grid_plots)) {
    row <- (i - 1) %/% ncol + 1
    col <- (i - 1) %% ncol + 1
    pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
    print(grid_plots[[i]], newpage = FALSE)
    grid.text(annotations[i], x = unit(1, "npc") - unit(9, "cm"), y = unit(1, "npc") - unit(0.5, "cm"), 
              just = c("left", "top"), gp = gpar(fontsize = 12, fontface = "bold"))
    popViewport()
  }
}

# Call the function to add annotations and arrange the plots
add_annotations(plots, theLetters, ncol = 3, nrow = 3)

# Compare the Response of the Cell Lines for the Mutations Dataset ######

gdscDrugs <- factor(unique(gdscDoseResponse$DRUG_NAME))

# Initialize results dataframe
CCLEmutGdscResults <- data.frame()

# Function to perform the analysis for each drug
perform_analysis <- function(curDrug, validCCLEmuts, CCLEmutations, gdscDoseResponse) {
  appendTable <- data.frame(DRUG_NAME = character(),
                            geneName = character(),
                            meanHighRank = numeric(),
                            meanLowerRank = numeric(),
                            tStat = numeric(),
                            lowerBound = numeric(),
                            upperBound = numeric(),
                            pValue = numeric(),
                            numHighRank = numeric(),
                            numLowRank = numeric(),
                            stringsAsFactors = FALSE)
  
  for (jj in 1:ncol(validCCLEmuts)) {
    geneName <- colnames(validCCLEmuts)[jj]
    
    if (jj %% 200 == 0) {
      cat(sprintf("Running analysis for drug %s and gene #%d: %s\n",
                  curDrug, jj, geneName))
    }
    
    curMutCellLines <- CCLEmutations$cell_line[!is.na(validCCLEmuts[[jj]])]
    
    if (length(curMutCellLines) < 10 || length(curMutCellLines) > nrow(validCCLEmuts) - 10) {
      next
    }
    
    curGDSC <- gdscDoseResponse %>% filter(DRUG_NAME == curDrug)
    locInGDSC <- curGDSC$cell_line %in% curMutCellLines
    
    numInstances <- table(locInGDSC)
    if (any(numInstances < 8)) {
      next
    }
    
    higherResponse <- curGDSC$Z_SCORE[locInGDSC]
    lowerResponse <- curGDSC$Z_SCORE[!locInGDSC]
    
    if (length(higherResponse) < 2 || length(lowerResponse) < 2) {
      next
    }
    
    ttest_result <- t.test(higherResponse, lowerResponse, var.equal = FALSE)
    
    meanHighRank <- mean(higherResponse)
    meanLowerRank <- mean(lowerResponse)
    
    appendTable <- rbind(appendTable, data.frame(DRUG_NAME = curDrug,
                                                 geneName = geneName,
                                                 meanHighRank = meanHighRank,
                                                 meanLowerRank = meanLowerRank,
                                                 tStat = ttest_result$statistic,
                                                 lowerBound = ttest_result$conf.int[1],
                                                 upperBound = ttest_result$conf.int[2],
                                                 pValue = ttest_result$p.value,
                                                 numHighRank = sum(locInGDSC),
                                                 numLowRank = sum(!locInGDSC)))
  }
  
  return(appendTable)
}


tryCatch({
  CCLEmutGdscResults <- read_excel('Supplementary_File_4.xlsx', sheet = 'Between Cell Mut Response')
  cat("Loaded existing data from 'Supplementary File 4.xlsx'\n")
}, error = function(e) {
  cat("Error loading file or file does not exist. Proceeding with analysis...\n")
  
  validCCLEmuts <- CCLEmutations[, 3:ncol(CCLEmutations)]
  
  mutFreq <- colSums(!sapply(validCCLEmuts, is.na))
  validGenes <- mutFreq >= 10
  validCCLEmuts <- validCCLEmuts[, validGenes]
  
  gdscDoseResponse$DRUG_NAME <- as.factor(gdscDoseResponse$DRUG_NAME)
  gdscDrugs <- factor(unique(gdscDoseResponse$DRUG_NAME))
  
  results_list <- lapply(gdscDrugs, function(curDrug) {
    cat("Analyzing drug: ", curDrug, "\n")
    result <- perform_analysis(curDrug, validCCLEmuts, CCLEmutations, gdscDoseResponse)
    if (nrow(result) > 0) {
      cat("Result for drug", curDrug, ":\n")
      print(head(result))
    }
    return(result)
  })
  
  CCLEmutGdscResults <- do.call(rbind, results_list)
  
  if (nrow(CCLEmutGdscResults) == 0) {
    stop("Error: The results from the computation are empty.")
  }
  
  CCLEmutGdscResults <- as.data.frame(CCLEmutGdscResults)
  
  cat("Structure of CCLEmutGdscResults before filtering:\n")
  str(CCLEmutGdscResults)
})

cat('\nNow cleaning up the t-test data\n')

# Rename columns
colnames(CCLEmutGdscResults)[2:10] <- c('HugoSymbol', 'meanMutated', 'meanNotMutated', 'tStat', 
                                        'lowerBound', 'upperBound', 'pValue', 'numMutated', 'numNotMutated')

# Add FDR column
CCLEmutGdscResults$FDR <- p.adjust(CCLEmutGdscResults$pValue, method = "BH")

# Filter for significant results
CCLEmutGdscResults <- CCLEmutGdscResults %>% filter(pValue < 0.05)

# Remove empty rows
CCLEmutGdscResults <- CCLEmutGdscResults[!is.na(CCLEmutGdscResults$HugoSymbol) & CCLEmutGdscResults$HugoSymbol != "", ]

# Add developmental stage to the table
CCLEmutGdscResults <- inner_join(pharos %>% select(HugoSymbol, DevStage), CCLEmutGdscResults, by = "HugoSymbol")

# Add drug targets to the table
gdscTargets <- gdscDoseResponse %>% distinct(DRUG_NAME, PUTATIVE_TARGET, PATHWAY_NAME)
CCLEmutGdscResults <- inner_join(CCLEmutGdscResults, gdscTargets, by = "DRUG_NAME")

# Reorder columns
CCLEmutGdscResults <- CCLEmutGdscResults %>%
  select(DRUG_NAME, PUTATIVE_TARGET, PATHWAY_NAME, everything())

# Sort the table by FDR
CCLEmutGdscResults <- CCLEmutGdscResults %>% arrange(FDR)

cat('\nSaving the data to Excel\n')

library(dplyr)
library(readxl)
library(openxlsx)

# Save the results to the supplementary file
write.xlsx(CCLEmutGdscResults, "Supplementary_File_4.1.xlsx", sheetName = "Between Cell Mut Response")

# Convert columns to factors
CCLEmutGdscResults$HugoSymbol <- as.factor(CCLEmutGdscResults$HugoSymbol)
CCLEmutGdscResults$PUTATIVE_TARGET <- as.factor(CCLEmutGdscResults$PUTATIVE_TARGET)

# Create a table with genes along the column and drug along the rows ######
# Here is the mutation data
data <- CCLEmutGdscResults

# Convert the drug name to categorical
if (is.character(data$DRUG_NAME)) {
  data$DRUG_NAME <- as.factor(data$DRUG_NAME)
}

# Here is the preallocated table
unique_drug_names <- unique(data$DRUG_NAME)
unique_hugo_symbols <- unique(data$HugoSymbol)
drugTable <- data.frame(matrix(nrow = length(unique_hugo_symbols), 
                               ncol = length(unique_drug_names)))
colnames(drugTable) <- as.character(unique_drug_names)
drugTable$HugoSymbol <- unique_hugo_symbols

# Initialize the table with zeros
drugTable[is.na(drugTable)] <- 0

# Add the t-values to the table
for (ii in 1:length(unique_drug_names)) {
  
  # Get the current drug data
  curData <- data %>% filter(DRUG_NAME == unique_drug_names[ii])
  
  # Get the indices of the intersecting Hugo symbols
  intersect_indices <- match(curData$HugoSymbol, drugTable$HugoSymbol)
  
  # Add the drug data to the drug table
  drugTable[intersect_indices, as.character(unique_drug_names[ii])] <- curData$tStat
}

# Remove rows with all missing data
locMissing <- drugTable[, 1:(ncol(drugTable) - 1)] == 0

# Remove rows and columns with all missing data
drugTable <- drugTable[rowSums(locMissing) != ncol(locMissing), ]
drugTable <- drugTable[, colSums(locMissing) != nrow(locMissing)]

# Display the first 10 rows and columns of the table
print(drugTable[1:10, 1:10])

library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)

# Sort the table by the FDR column in ascending order
CCLEmutGdscResults <- CCLEmutGdscResults[order(CCLEmutGdscResults$FDR), ]

# Get the most significant results based on numMutated and numNotMutated
locTop <- CCLEmutGdscResults$numMutated > 10 & 
  CCLEmutGdscResults$numNotMutated > 10

# Return only those instances
top9Results <- CCLEmutGdscResults[locTop, ]

# Now get the top 9 results to be plotted
top9Results <- top9Results[1:9, ]

# Get the valid CCLE mutations (excluding the first two columns)
validCCLEmuts <- CCLEmutations[, 3:ncol(CCLEmutations)]

# Find the mutation frequency of each gene
mutFreq <- colSums(!is.na(validCCLEmuts))

# Get the genes that have at least 10 cell lines mutated or at least 10 cell lines not mutated
validGenes <- mutFreq >= 10

# Subset validCCLEmuts to include only these genes
validCCLEmuts <- validCCLEmuts[, validGenes]


# Define colors for mutations
CCLEmutationsColor <- c("#FC4E2A", "#2B83BA")  # Red for "Mutated", Blue for "Unmutated"

# Function to create the scatter and box plots
colourBoxPlot <- function(plotData, groups, colors, includeScatter, showYAxisLabel) {
  group_labels <- ifelse(groups, "Mutated", "Unmutated")
  p <- ggplot(data.frame(plotData, group_labels), aes(x = factor(group_labels), y = plotData)) +
    geom_boxplot(aes(fill = factor(group_labels)), width = 0.3, outlier.shape = 19, outlier.size = 1.5, outlier.colour = "black") +
    scale_fill_manual(values = colors) +
    scale_colour_manual(values = colors) +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "bold", size = 8),
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      legend.position = "none"
    ) +
    labs(
      x = "Mutation Status"
    )
  if (showYAxisLabel) {
    p <- p + labs(y = "IC_50 Z-score")
  } else {
    p <- p + theme(axis.title.y = element_blank())
  }
  if (includeScatter) {
    p <- p + geom_jitter(width = 0.2, aes(colour = factor(group_labels)), alpha = 0.5, size = 1)
  }
  return(p)
}

# Letters for annotating plots
theLetters <- letters[1:9]

# Create a list to store the plots
plots <- list()

# Loop through each row of top9Results to plot individual panels
for (ii in 1:nrow(top9Results)) {
  # Get the GDSC data for the current drug
  curGDSC <- gdscDoseResponse[gdscDoseResponse$DRUG_NAME == top9Results$DRUG_NAME[ii], ]
  
  # Get the cell lines highly dependent on the current genes
  curMutCellLines <- CCLEmutations$cell_line[!is.na(validCCLEmuts[[top9Results$HugoSymbol[ii]]])]
  
  # Get the location of those cell lines in the current GDSC data
  locInGDSC <- curGDSC$cell_line %in% curMutCellLines
  
  # Fill outliers in Z_SCORE with nearest mean
  dReponses <- curGDSC$Z_SCORE
  dReponses[is.na(dReponses)] <- ave(dReponses, FUN = function(x) mean(x, na.rm = TRUE))
  
  # Determine whether to show y-axis label
  showYAxisLabel <- ii %% 3 == 1
  
  # Plot the box plots
  includeScatter <- TRUE
  p <- colourBoxPlot(dReponses, locInGDSC, CCLEmutationsColor, includeScatter, showYAxisLabel) +
    labs(title = paste(top9Results$DRUG_NAME[ii], ":", top9Results$HugoSymbol[ii], ":", top9Results$DevStage[ii])) 
  
  plots[[ii]] <- p
}

# Function to add annotations to the grid
add_annotations <- function(grid_plots, annotations, ncol, nrow) {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow, ncol, heights = unit(rep(1, nrow), "null"), widths = unit(rep(1, ncol), "null"))))
  for (i in seq_along(grid_plots)) {
    row <- (i - 1) %/% ncol + 1
    col <- (i - 1) %% ncol + 1
    pushViewport(viewport(layout.pos.row = row, layout.pos.col = col))
    print(grid_plots[[i]], newpage = FALSE)
    grid.text(annotations[i], x = unit(0.95, "npc"), y = unit(0.95, "npc"), 
              just = c("left", "top"), gp = gpar(fontsize = 12, fontface = "bold"))
    popViewport()
  }
}

# Call the function to add annotations and arrange the plots
add_annotations(plots, theLetters, ncol = 3, nrow = 3)

# Add a legend to the figure
createLegendInternal <- function(x, y, labels, colors, title, fontSizes, rectAndTextBox) {
  legend <- legendGrob(labels = labels, 
                       pch = 16, 
                       gp = gpar(col = colors, fontsize = fontSizes[1]), 
                       title = title, title_gp = gpar(fontsize = fontSizes[2]))
  grid.draw(legend, vp = viewport(x = x, y = y, just = c("right", "top")))
}

createLegendInternal(0.86, 0.90, c("Mutated", "Unmutated"), CCLEmutationsColor, "Mutations", c(11, 10), c(0.1, 0.12))

