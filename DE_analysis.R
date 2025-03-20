# Loading Required Libraries
library(edgeR)
library(ggplot2)
library(knitr)

# Parameters (replace these with actual values or parameterize your script)
params <- list()
params$TableName <- "C:/Users/kejia/OneDrive/Desktop/Kejia_research/MLthings/withoutbaseline/NS3astar/S3C9CdenovoMBX/R3forDE.csv"
params$TestCols <- "2"
params$ControlCols <- "1"
params$saveName <- "C:/Users/kejia/OneDrive/Desktop/Kejia_research/MLthings/withoutbaseline/NS3astar/S3C9CdenovoMBX/R3DE.txt"

# Load the data

TestCols <- as.numeric(unlist(strsplit(params$TestCols,',')))
ControlCols <- as.numeric(unlist(strsplit(params$ControlCols, ',')))

# Load and prepare the dataset
counts <- read.table(params$TableName, header = TRUE, sep = ",", quote = "\"", check.names = FALSE)
print(ncol(counts))  # This should now correctly show the number of columns
print(length(TestCols) + length(ControlCols))
print(ncol(counts)-1)

if(ncol(counts)-1 != length(TestCols) + length(ControlCols)){
  stop("Number of control and test columns does not match the dataset!")
}

group <- factor(c(rep("control", length(ControlCols)), rep("test", length(TestCols))))

# Assign sequence to row names
row.names(counts) <- counts[,1]
counts <- counts[,-1]

# Differential Expression Analysis using edgeR
counts[is.na(counts)] <- 0
RG <- DGEList(counts = counts, group = group)
RG <- calcNormFactors(RG)

tmm_data <- cpm(RG, normalized.lib.sizes = TRUE)
final_table <- data.frame(sequence = rownames(RG$counts),tmm_data)
write.table(final_table, params$saveName, quote = FALSE, sep = "\t", row.names = FALSE)

