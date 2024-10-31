library(readxl)
library(dplyr)
library(readr)

# Load LLS and LSp sample ids
TCGACOAD <- read_excel("~/Analysis/2024/CRC/TCGACOAD.xlsx")
# Load mutation data
mutation_data <- read.delim("~/Analysis/2024/CRC/coadread_tcga_pan_can_atlas_2018/data_mutations.txt", header = TRUE, sep = "\t")

likely_LS_samples <- na.omit(TCGACOAD$Likely_Lynch) 
likely_LS_samples <- likely_LS_samples[likely_LS_samples != "TCGA-AA-A010-01"]
likely_sporadic_samples <- na.omit(TCGACOAD$Likely_Sporadic)

mutation_data <- mutation_data %>%
  mutate(variant = paste(Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2, sep = "_"))

num_likely_LS <- length(intersect(unique(mutation_data$Tumor_Sample_Barcode), likely_LS_samples))
num_likely_sporadic <- length(intersect(unique(mutation_data$Tumor_Sample_Barcode), likely_sporadic_samples))

all_samples <- c(likely_LS_samples, likely_sporadic_samples)

# Load allmuts2 which contains all pathogenic variants in all CRC cbioportal studies
allmuts2 <- read_delim("~/Analysis/2024/1. Download cBioportal data/allmuts2_crc_new.csv", 
                       delim = ",", escape_double = FALSE, trim_ws = TRUE)

# Filter for pathogenic variants 
mutation_data <- mutation_data %>%
  filter(variant %in% allmuts2$newID)

# Subsetting mutation_data
mutation_data_likely_LS <- mutation_data[mutation_data$Tumor_Sample_Barcode %in% likely_LS_samples, ]
mutation_data_likely_sporadic <- mutation_data[mutation_data$Tumor_Sample_Barcode %in% likely_sporadic_samples, ]

length(unique(mutation_data_likely_LS$Tumor_Sample_Barcode))
length(unique(mutation_data_likely_sporadic$Tumor_Sample_Barcode))

# number of somatically altered pathogenic genes per sample ####
# Group by sample and calculate the number of unique genes per sample
sample_gene_counts <- mutation_data_likely_LS %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(Unique_Genes = length(unique(Hugo_Symbol)))

# Calculate the median number of unique genes mutated per sample
median_genes_mutated <- median(sample_gene_counts$Unique_Genes)

# Calculate Q1 and Q3 for the IQR
Q1 <- quantile(sample_gene_counts$Unique_Genes, 0.25)
Q3 <- quantile(sample_gene_counts$Unique_Genes, 0.75)

IQR <- Q3 - Q1

# Print the results
print(paste("Median number of unique genes mutated per sample:", median_genes_mutated))
print(paste("Q1 of unique genes mutated per sample:", round(Q1, 2)))
print(paste("Q3 of unique genes mutated per sample:", round(Q3, 2)))
print(paste("IQR of unique genes mutated per sample:", round(IQR, 2)))

# Likely Sporadic 
sample_gene_counts_sporadic <- mutation_data_likely_sporadic %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(Unique_Genes = length(unique(Hugo_Symbol)), .groups = 'drop')

likely_sporadic_df <- data.frame(Tumor_Sample_Barcode = likely_sporadic_samples)

median_genes_mutated_sporadic <- median(sample_gene_counts_sporadic$Unique_Genes)
Q1_sporadic <- quantile(sample_gene_counts_sporadic$Unique_Genes, 0.25)
Q3_sporadic <- quantile(sample_gene_counts_sporadic$Unique_Genes, 0.75)

IQR_sporadic <- Q3_sporadic - Q1_sporadic

gene_mutation_range_sporadic <- range(sample_gene_counts_sporadic$Unique_Genes)

# Print the results
print(paste("Median number of unique genes mutated per sample:", median_genes_mutated_sporadic))
print(paste("Q1 of unique genes mutated per sample:", round(Q1_sporadic, 2)))
print(paste("Q3 of unique genes mutated per sample:", round(Q3_sporadic, 2)))
print(paste("IQR of unique genes mutated per sample:", round(IQR_sporadic, 2)))
print(paste("Range of unique genes mutated per sample: Min =", gene_mutation_range_sporadic[1], "Max =", gene_mutation_range_sporadic[2]))

# number of somatically altered pathogenic variants per sample ####
#Likely LS
sample_variant_counts <- mutation_data_likely_LS %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(Unique_Variants = length(unique(variant)))

median_variant_mutated <- median(sample_variant_counts$Unique_Variants)
print(paste("Median number of unique variants mutated per sample:", median_variant_mutated))

quartiles <- quantile(sample_variant_counts$Unique_Variants, probs = c(0.25, 0.75))
iqr_variants <- paste("IQR,", quartiles[1], "-", quartiles[2])
iqr_variants

# Likely Sporadic 
sample_variant_counts_sporadic <- mutation_data_likely_sporadic %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(Unique_Variants = length(unique(variant)), .groups = 'drop')

median_variant_sporadic <- median(sample_variant_counts_sporadic$Unique_Variants)

print(paste("Median number of variants mutated per sample:", median_variant_sporadic))

quartiles_sporadic <- quantile(sample_variant_counts_sporadic$Unique_Variants, probs = c(0.25, 0.75))
iqr_variants_sporadic <- paste("IQR,", quartiles_sporadic[1], "-", quartiles_sporadic[2])
iqr_variants_sporadic

group1 <- sample_variant_counts$Unique_Variants
group2 <- sample_variant_counts_sporadic$Unique_Variants

# Perform the Mann-Whitney U Test
test_result <- wilcox.test(group1, group2, alternative = "two.sided")

print(test_result$p.value)

formatted_p_value <- format(test_result$p.value, scientific = FALSE)

# Get the number of samples with a mutation in each gene ####
gene_mutation_counts_likely_LS <- mutation_data_likely_LS %>%
  group_by(Hugo_Symbol) %>%                  # Group data by gene symbol
  summarise(Unique_Samples = n_distinct(Tumor_Sample_Barcode)) %>%  # Count unique sample IDs
  arrange(desc(Unique_Samples))              # Arrange the data in descending order of counts

gene_mutation_counts_likely_sporadic <- mutation_data_likely_sporadic %>%
  group_by(Hugo_Symbol) %>%                  # Group data by gene symbol
  summarise(Unique_Samples = n_distinct(Tumor_Sample_Barcode)) %>%  # Count unique sample IDs
  arrange(desc(Unique_Samples))              # Arrange the data in descending order of counts

One_third_gene_mutation_counts_likely_LS <- gene_mutation_counts_likely_LS[gene_mutation_counts_likely_LS$Unique_Samples >= 4, ]

# Merge the data frames on Hugo_Symbol
merged_mutation_counts <- merge(
  gene_mutation_counts_likely_LS, 
  gene_mutation_counts_likely_sporadic, 
  by = "Hugo_Symbol", 
  suffixes = c("_likely_LS", "_likely_Sporadic"), 
  all = TRUE  
)

# Replace NA values with 0 in the merged data frame
merged_mutation_counts[is.na(merged_mutation_counts)] <- 0

# Initialize a data frame to store the results
fisher_results <- data.frame(Hugo_Symbol = character(), p_value = numeric(), stringsAsFactors = FALSE)

# Loop through each gene to perform Fisher's test
for (i in 1:nrow(merged_mutation_counts)) {
  gene <- merged_mutation_counts$Hugo_Symbol[i]
  count_LS <- merged_mutation_counts$Unique_Samples_likely_LS[i]
  count_Sporadic <- merged_mutation_counts$Unique_Samples_likely_Sporadic[i]
  
  # Correctly construct the contingency table
  matrix <- matrix(
    c(count_LS, num_likely_LS - count_LS,
      count_Sporadic, num_likely_sporadic - count_Sporadic),
    nrow = 2,                      
    byrow = TRUE,                
    dimnames = list(
      c("Mutated", "Not Mutated"), 
      c("LS", "Sporadic")
    )
  )
  
  # Perform Fisher's Exact Test
  test <- fisher.test(matrix)
  
  # Append the results
  fisher_results <- rbind(fisher_results, data.frame(Hugo_Symbol = gene, p_value = test$p.value))
}

# Benjamini Hochberg procedure
fisher_results$adj_p_value <- p.adjust(fisher_results$p_value, method = "BH")

# Filter results to include only significant genes after adjustment
significant_genes <- fisher_results[fisher_results$adj_p_value <= 0.05, ]

# Add counts and frequencies
# Merge significant genes with mutation counts for Likely Lynch
significant_genes <- merge(
  significant_genes, 
  gene_mutation_counts_likely_LS[, c("Hugo_Symbol", "Unique_Samples")],
  by = "Hugo_Symbol",
  all.x = TRUE,
  suffixes = c("", "_LS")
)

# Merge significant genes with mutation counts for Likely Sporadic
significant_genes <- merge(
  significant_genes, 
  gene_mutation_counts_likely_sporadic[, c("Hugo_Symbol", "Unique_Samples")],
  by = "Hugo_Symbol",
  all.x = TRUE,
  suffixes = c("", "_Sporadic")
)

# Replace NA values in the new columns with 0
significant_genes[is.na(significant_genes)] <- 0

# Rename columns for clarity
names(significant_genes)[names(significant_genes) == "Unique_Samples"] <- "Num_Samples_Likely_LS"
names(significant_genes)[names(significant_genes) == "Unique_Samples_Sporadic"] <- "Num_Samples_Likely_Sporadic"

print(significant_genes)

significant_genes$Freq_Likely_LS <- significant_genes$Num_Samples_Likely_LS / num_likely_LS
significant_genes$Freq_Likely_Sporadic <- significant_genes$Num_Samples_Likely_Sporadic / num_likely_sporadic

# Enrichment 
# Calculate the log2-based ratio
significant_genes$Log2 <- log2(significant_genes$Freq_Likely_LS /significant_genes$Freq_Likely_Sporadic)

significant_genes <- significant_genes %>%
  mutate(Enriched_in = case_when(
    Log2 > 0 ~ "Likely Lynch",
    Log2 <= 0 ~ "Likely Sporadic",
    TRUE ~ ""
  ))

significant_genes$p_value <- format(significant_genes$p_value, scientific = FALSE)
significant_genes$adj_p_value <- format(significant_genes$adj_p_value, scientific = FALSE)

# Filter out genes found in less than 1/3 of lynch samples
filtered_significant_genes <- significant_genes[significant_genes$Num_Samples_Likely_LS >= (num_likely_LS*(1/3)), ]

# Filter for only oncogenes and TSGs
oncokb <- read.delim("~/Analysis/2024/Mutation and Copy Number Alteration Analyses/cancerGeneList.tsv", header = TRUE, sep = "\t")

oncokb <- oncokb %>%
  filter(!(Is.Oncogene == "No" & Is.Tumor.Suppressor.Gene == "No"))

filtered_significant_genes <- semi_join(filtered_significant_genes, oncokb, by = c("Hugo_Symbol" = "Hugo.Symbol"))

file_path <- "~/Analysis/2024/Mutation and Copy Number Alteration Analyses/New_August_24/CRC_fishers.csv"

# Save the data frame as a CSV file
write.csv(significant_genes, file = file_path, row.names = FALSE)


file_path <- "~/Analysis/2024/Mutation and Copy Number Alteration Analyses/New_August_24/filtered_CRC_fishers.csv"
write.csv(filtered_significant_genes, file = file_path, row.names = FALSE)

#Enrichr
# install.packages("enrichR")

library(enrichR)
websiteLive <- getOption("enrichR.live")
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes   
}

if (websiteLive) dbs <- listEnrichrDbs()
if (websiteLive) head(dbs)
head <- head(dbs)

dbs <- c("KEGG_2021_Human")
#CRC likely ls
if (websiteLive) {
  enriched <- enrichr(filtered_significant_genes$Hugo_Symbol, dbs)
}

if (websiteLive) enriched[["KEGG_2021_Human"]]

if (websiteLive) {
  plotEnrich(enriched[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
}

png("~/Analysis/2024/Mutation and Copy Number Alteration Analyses/New_August_24/enrichr_CRC.png") 

if (websiteLive) {
  plotEnrich(enriched[[1]], showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value")
} 

dev.off()


