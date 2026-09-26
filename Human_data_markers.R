library(Seurat)
library(ggplot2)

human_data <- readRDS(/Users/brendakarumbo/Desktop/PhD/Project/IM_samples/stacas.rds)

# ---- 1. Cell type order (matches your reference plot) --------------------
# Types in the reference that aren't in your object (cDC1, cDC2, ASDC,
# Plasmablast, NK_CD56bright) are dropped automatically; your "DC" sits
# where cDC1/cDC2 were.
celltype_order <- c(
  "HSPC", "GMP", "Neutrophil Precursor", "Neutrophil",
  "CD14 Mono", "CD16 Mono",
  "DC", "cDC1", "cDC2", "pDC", "ASDC",
  "B naive", "B intermediate", "B memory",
  "Plasmablast", "Plasma cell precursor", "Plasma cell",
  "CD4 Naive", "CD4 TCM", "CD4 TEM", "CD4 CTL", "Treg", "CD4 Proliferating",
  "CD8 Naive", "CD8 TCM","CD8 TEM","CD8 Proliferating",
  "MAIT", "gdT", "dnT", "ILC",
  "NK", "NK_CD56bright", "NK Proliferating",
  "Platelet", "Eryth"
)

# Change "celltype.l2" to whatever metadata column holds your labels
Idents(human_data) <- "celltypes"

present <- intersect(celltype_order, levels(Idents(human_data)))
missing <- setdiff(levels(Idents(human_data)), celltype_order)
if (length(missing)) message("Not in order list (added at end): ",
                             paste(missing, collapse = ", "))
Idents(human_data) <- factor(Idents(human_data), levels = c(present, missing))

# ---- 2. Marker genes, grouped in the same order --------------------------
markers <- list(
  "HSPC"           = c("CD34", "SPINK2", "CRHBP"),
  "Neutrophils"    = c("DEFA3","ELANE", "CAMP", "FCGR3B","CSF3R"),
  "Monocytes"      = c("CD14","CX3CR1",
                       "C1QA", "MS4A7"),
  "DC"             = c("CLEC9A", "FCER1A", "PLD4"),
  "pDC"            = c("SIGLEC6","LILRA4","IL3RA"),
  "B cells"        = c( "IGHM","TCL1A", "NIBAN3","MS4A1", "CD79A", "BANK1"),
  "Plasmablast"    = c("TNFRSF17","JCHAIN", "MZB1"),
  "CD4 T"          = c( "CD4", "GZMK", "RGS1"),
  "Proliferating"  = c("MCM3", "CCNB2", "STMN1", "CDK6"),
  "CD8 T"          = c("CD8A", "CD8B", "CCL5", "CST7"),
  "MAIT / gdT"     = c("KLRB1", "TRDC", "TRGC1", "TRGC2", "TRAV1-2"),
  "NK / cytotoxic" = c("NKG7", "GNLY", "PCLAF", "TYMS"),
  "Platelets"      = c("PPBP","F13A1", "PF4"),
  "Erythroid"      = c("HBB", "SLC25A37", "ALAS2")
)

# ---- 3. Plot --------------------------------------------------------------
p <- DotPlot(human_data, features = markers, cluster.idents = FALSE,
             dot.scale = 5, col.min = -1, col.max = 2.5) +
  scale_colour_gradient(low = "grey85", high = "#B2182B",
                        name = "Average\nexpression") +
  scale_y_discrete(limits = rev(c(present, missing))) +   # first type at top
  guides(size = guide_legend(title = "Percent\nexpressed")) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x   = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                 face = "italic", size = 12),
    axis.text.y   = element_text(size = 13, colour = "black"),
    axis.line     = element_blank(),
    panel.border  = element_rect(colour = "grey60", fill = NA, linewidth = 0.4),
    panel.spacing = unit(1, "mm"),
    strip.background = element_blank(),
    strip.text.x  = element_text(angle = 90, hjust = 0, size = 12,
                                 face = "bold"),
    strip.clip    = "off",
    legend.title  = element_text(size = 13),
    legend.text   = element_text(size = 12)
  )

print(p)
ggsave("dotplot_clean.pdf", p, width = 18, height = 12)
ggsave("dotplot_clean.png", p, width = 18, height = 12, dpi = 300)




