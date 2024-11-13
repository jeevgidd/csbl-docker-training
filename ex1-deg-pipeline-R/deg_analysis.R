# Load required libraries
library(limma)
library(edgeR)
library(ggplot2)
library(dplyr)

# Read data
counts <- read.csv("data/counts.csv", row.names=1)
metadata <- read.csv("data/metadata.csv", row.names=1)

# Create DGEList object
dge <- DGEList(counts = counts)

# Filter low expression genes
keep <- filterByExpr(dge, group=metadata$condition)
dge <- dge[keep, ]

# Calculate normalization factors
dge <- calcNormFactors(dge)

# Create design matrix
design <- model.matrix(~condition, data=metadata)

# Convert counts to log2-CPM with associated weights
v <- voom(dge, design, plot=TRUE)

# Fit linear model
fit <- lmFit(v, design)
fit <- eBayes(fit)

# Get results
results <- topTable(fit, coef=2, n=Inf)
results$gene_id <- rownames(results)

# Add regulation status column
results$regulation <- "NS"  # Non-significant
results$regulation[results$adj.P.Val < 0.05 & results$logFC > 1] <- "UP"
results$regulation[results$adj.P.Val < 0.05 & results$logFC < -1] <- "DOWN"

# Convert regulation to factor with specific order
results$regulation <- factor(results$regulation, levels = c("UP", "DOWN", "NS"))

ggplot(results, aes(x=logFC, y=-log10(adj.P.Val), color=regulation)) +
  geom_point(alpha=0.6, size=2) +
  theme_minimal(base_size = 12) +
  scale_color_manual(values=c("UP"="red", "DOWN"="blue", "NS"="grey")) +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="darkgrey") +
  geom_vline(xintercept = c(-1, 1), linetype="dashed", color="darkgrey") +
  labs(title="Differential Expression Analysis",
       subtitle="Red: Upregulated, Blue: Downregulated (FDR < 0.05, |log2FC| > 1)",
       x="log2 Fold Change",
       y="-log10 adjusted p-value",
       color="Regulation") +
  theme(
    plot.title = element_text(size=16, face="bold"),
    plot.subtitle = element_text(size=10),
    legend.position = "right",
    legend.title = element_text(size=12),
    legend.text = element_text(size=10),
    axis.title = element_text(size=12),
    axis.text = element_text(size=10)
  )


# Create enhanced volcano plot
dir.create("results", showWarnings = FALSE)
pdf("results/volcano_plot.pdf", width=10, height=8)
ggplot(results, aes(x=logFC, y=-log10(adj.P.Val), color=regulation)) +
  geom_point(alpha=0.6, size=2) +
  theme_minimal(base_size = 12) +
  scale_color_manual(values=c("UP"="red", "DOWN"="blue", "NS"="grey")) +
  geom_hline(yintercept = -log10(0.05), linetype="dashed", color="darkgrey") +
  geom_vline(xintercept = c(-1, 1), linetype="dashed", color="darkgrey") +
  labs(title="Differential Expression Analysis",
       subtitle="Red: Upregulated, Blue: Downregulated (FDR < 0.05, |log2FC| > 1)",
       x="log2 Fold Change",
       y="-log10 adjusted p-value",
       color="Regulation") +
  theme(
    plot.title = element_text(size=16, face="bold"),
    plot.subtitle = element_text(size=10),
    legend.position = "right",
    legend.title = element_text(size=12),
    legend.text = element_text(size=10),
    axis.title = element_text(size=12),
    axis.text = element_text(size=10)
  )
dev.off()


# Save results
write.csv(results, "results/deg_results.csv")

# Print summary
cat("\nDifferential Expression Analysis Summary:\n")
cat("Total genes analyzed:", nrow(results), "\n")
cat("Significantly DE genes (FDR < 0.05, |log2FC| > 1):", sum(results$regulation != "NS"), "\n")
cat("  - Upregulated:", sum(results$regulation == "UP"), "\n")
cat("  - Downregulated:", sum(results$regulation == "DOWN"), "\n")