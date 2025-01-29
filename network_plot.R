# Load the required libraries
library(org.Hs.eg.db)
library(clusterProfiler)
library(pathview)
library(ggplot2)
library(cowplot)
library(enrichplot)
library(showtext)
library(ggraph)
library(Cairo)
library(grid)

# Font for the title
font_add("Libre Baskerville Bold", "libre-baskerville.bold.ttf")
font_add("TimesNewRomanPSMT", "C:/WINDOWS/Fonts/times.ttf")
showtext_auto()


# Read the file
gene_data <- read.csv("file.csv", header = T)
gene <- gene_data$Entrez_ID

# Kegg enrichment analysis
kegg_res <- enrichKEGG(
  gene = gene, 
  organism = 'hsa', 
  pAdjustMethod = 'BH', 
  pvalueCutoff = 0.05)

head(kegg_res)

# Similarity matrix
kegg_sim <- pairwise_termsim(kegg_res)

# Set the gene readability
kegg_gene <- setReadable(kegg_sim, 'org.Hs.eg.db', 'ENTREZID') #@result

# Select the fold change for the network
fold_change <- gene_data$logFC
names(fold_change) <- gene_data$GENES

# Set the overlap value
options(ggrepel.max.overlaps = Inf)

# Specify the important categories for visualization
category <- c("NF-kappa B signaling pathway", "Cytokine-cytokine receptor interaction", 
              "Toll-like receptor signaling pathway", "TNF signaling pathway",
              "NOD-like receptor signaling pathway")

# Create a cnetplot for the enriched pathways
plot <- cnetplot(
  kegg_gene,
  showCategory = category,  # Number of pathways to display
  circular = F,
  color.params = list(foldChange = fold_change, category = "#F70D1A", edge = T),
  cex.params = list(category_label = 0.99, gene_label = 1),
  layout = 'circle',
  shadowtext = "all"
  )

# Upgrade with ggplot2 functions
plot_2 <- plot + 
  ggtitle("Pathway Dynamics in AML relapse: A KEGG Perspective",
          subtitle = "Acute myelogenous leukemia, also called AML, is one of the most common type of leukemia in adults. \nUnraveling the Complex Molecular Pathways Underlying AML Relapse: \nBridging Genomic Alterations to Clinical Outcomes Through \nKey Players Like 'NFKBIA' and 'IL6'.") +
  scale_color_viridis_c(option = "viridis", name = 'LogFC')+ 
  theme(
    plot.title = element_text(family = "Libre Baskerville Bold",
                              size = 30,
                              face = "bold", 
                              colour = "black"),
    plot.background = element_rect(fill = NA, colour = NA),
    panel.background = element_blank(),
    text = element_text(color = "black", size = 10),
    legend.position = "right",
    legend.text = element_text(family = "TimesNewRomanPSMT", 
                               size = 14),
    plot.subtitle = element_text(family = "TimesNewRomanPSMT",
                                 size = 18),
    plot.margin = margin(10,15,15,15)
  )+
  scale_edge_alpha(range = c(1, 1))+
  guides(size = "none")

# Gradient background
gradient <- rasterGrob(
  colorRampPalette(c("#D8BFD8", "#F0F8FF"))(256),
  width = unit(1, "npc"), height = unit(1, "npc"),
  interpolate = TRUE
)

grid.newpage()
grid.draw(gradient)
print(plot_2, newpage = FALSE)
