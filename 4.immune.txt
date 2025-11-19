library(IOBR)
library(ltc)
names(palettes)
pltc(ltc("paloma"))

estimate <- deconvo_tme(eset = eset, method = "estimate")
input_estimate <- merge(estimate, pdata_cmoic, by = "ID")
input_estimate$CMOIC <- factor(input_estimate$CMOIC, levels = c("CS1", "CS2", "CS3", "CS4", "CS5"))
head(input_estimate)

comparision <- combn(unique(as.character(input_estimate$CMOIC)), 2, simplify = FALSE)

p1 <- ggviolin(
  input_estimate, "CMOIC", "StromalScore_estimate",
  fill = "CMOIC",
  palette = ggsci::pal_lancet()(6),
  add = "boxplot", add.params = list(fill = "white")
) +
  stat_compare_means(label.y = 3400, size = 10) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 24),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 28)
  )

p2 <- ggviolin(
  input_estimate, "CMOIC", "ImmuneScore_estimate",
  fill = "CMOIC",
  palette = ggsci::pal_lancet()(6),
  add = "boxplot", add.params = list(fill = "white")
) +
  stat_compare_means(label.y = 4800, size = 10) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 24),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 28)
  )

p3 <- ggviolin(
  input_estimate, "CMOIC", "ESTIMATEScore_estimate",
  fill = "CMOIC",
  palette = ggsci::pal_lancet()(6),
  add = "boxplot", add.params = list(fill = "white")
) +
  stat_compare_means(label.y = 6800, size = 10) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 24),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 28)
  )

p3 | p2 | p1

cell <- deconvo_tme(eset = eset, method = "cibersort", arrays = TRUE, perm = 1000, absolute.mode = TRUE)
head(cell)

input_cell <- merge(input_estimate, cell, by = "ID")
input_cell <- surv.info %>%
  rownames_to_column(var = "ID") %>%
  dplyr::select(ID, Immune.Subtype) %>%
  merge(., input_cell, by = "ID") %>%
  arrange(CMOIC) %>%
  dplyr::select(-TumorPurity_estimate)
head(input_cell)

dat.cell.p <- data.frame()
for (i in colnames(input_cell[, 7:28])) {
  aa <- input_cell[, c(i, "CMOIC")]
  colnames(aa)[1] <- "var"
  aa <- kruskal.test(var ~ CMOIC, data = aa)
  aa <- data.frame(var = i, p = aa$p.value)
  dat.cell.p <- rbind(dat.cell.p, aa)
}
dat.cell.p$p.label <- ifelse(
  dat.cell.p$p < 0.0001, "****",
  ifelse(
    dat.cell.p$p < 0.001, "***",
    ifelse(
      dat.cell.p$p < 0.01, "**",
      ifelse(dat.cell.p$p < 0.05, "*", "")
    )
  )
)
cellname <- str_split(dat.cell.p$var, "_CI") %>% do.call(rbind, .) %>% .[, 1]
dat.cell.p$cell <- paste0(cellname, dat.cell.p$p.label)

input_cellb <- input_cell
colnames(input_cellb)[7:28] <- dat.cell.p$cell
input_cellb <- input_cellb %>%
  dplyr::select(
    T_cells_CD4_naive, `T_cells_CD4_memory_resting****`,
    `T_cells_CD4_memory_activated****`, `T_cells_CD8****`,
    `T_cells_follicular_helper****`, `T_cells_regulatory_(Tregs)**`,
    T_cells_gamma_delta,
    `B_cells_naive**`, `B_cells_memory**`,
    Plasma_cells, `NK_cells_resting**`, NK_cells_activated,
    `Macrophages_M0****`, `Macrophages_M1****`, `Macrophages_M2****`,
    `Dendritic_cells_resting****`, `Dendritic_cells_activated****`,
    `Mast_cells_resting****`, `Mast_cells_activated*`,
    `Monocytes****`, Neutrophils, Eosinophils
  )

input_cell2 <- t(scale(input_cellb))
dim(input_cell2)

library(ComplexHeatmap)
ha <- HeatmapAnnotation(
  CS = input_cell$CMOIC,
  ImmuneSubtype = input_cell$Immune.Subtype,
  ESTIMATEScore = input_cell$ESTIMATEScore_estimate,
  ImmuneScore = input_cell$ImmuneScore_estimate,
  StromalScore = input_cell$StromalScore_estimate,
  col = list(
    CS = c(
      "CS1" = ggsci::pal_lancet()(6)[1], "CS2" = ggsci::pal_lancet()(6)[2],
      "CS3" = ggsci::pal_lancet()(6)[3], "CS4" = ggsci::pal_lancet()(6)[4],
      "CS5" = ggsci::pal_lancet()(6)[5]
    ),
    ImmuneSubtype = c(
      "C1" = ltc::palettes$paloma[1], "C2" = ltc::palettes$paloma[2],
      "C3" = ltc::palettes$paloma[3], "C4" = ltc::palettes$paloma[4],
      "C6" = ltc::palettes$paloma[5]
    )
  ),
  na_col = "white",
  show_legend = rep(TRUE, 6),
  annotation_height = unit(rep(4, 6), "mm"),
  annotation_legend_param = list(
    CS = list(title = "CS"),
    ImmuneSubtype = list(title = "Immune Subtype"),
    ESTIMATEScore = list(title = "ESTIMATE Score"),
    ImmuneScore = list(title = "Immune Score"),
    StromalScore = list(title = "Stromal Score")
  ),
  border = TRUE
)

col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("#3288BD", "white", "#D53E4F"))
ht <- Heatmap(
  input_cell2,
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  border = TRUE,
  top_annotation = ha,
  heatmap_legend_param = list(title = "value", legend_height = unit(4, "cm")),
  column_split = input$CMOIC,
  row_split = c(rep("Lymphoid", 12), rep("Myeloid", 10)),
  show_column_dend = FALSE,
  show_row_dend = FALSE,
  column_title_gp = gpar(fontsize = 18, fontface = "bold"),
  clustering_method_columns = "single",
  clustering_distance_columns = "euclidean",
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "single",
  column_dend_height = unit(10, "mm")
)
print(ht)

colnames(input_cell)

ggstatsplot::ggbetweenstats(
  data = input_cell,
  x = CMOIC,
  y = Monocytes_CIBERSORT,
  messages = FALSE
) +
  ylab(i) +
  xlab(NULL) +
  theme(
    plot.title = element_text(size = 24),
    axis.text.x = element_text(size = 18),
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text.y = element_text(size = 18)
  )

input_expr <- combine_pd_eset(eset = eset, pdata = pdata_cmoic, scale = TRUE)
ggstatsplot::ggbetweenstats(
  data = input_expr,
  x = CMOIC,
  y = CTLA4,
  messages = FALSE
) +
  ylab(i) +
  xlab(NULL) +
  theme(
    plot.title = element_text(size = 24),
    axis.text.x = element_text(size = 18),
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.text.y = element_text(size = 18)
  )

library(ComplexHeatmap)
library(circlize)
library(ChAMPdata)
library(data.table)
library(genefu)
data("probe.features")

heatmap.BlWtRd <- c("#1F66AC", "#75AFD3", "grey90", "#FAB99B", "#B2192B")

immunomodulator <- read.table(
  "r0.data/immunomodulator.txt",
  sep = "\t", row.names = 1, check.names = FALSE,
  stringsAsFactors = FALSE, header = TRUE
)
expr <- eset
is.element(rownames(immunomodulator), rownames(expr))

sinfo <- pdata_cmoic %>% column_to_rownames(var = "ID")
colnames(sinfo) <- "subtype"
expr <- expr[rownames(immunomodulator), ] %>% na.omit()

meth <- fread(
  "r0.data/TCGA.BRCA.sampleMap_HumanMethylation450.gz",
  sep = "\t", check.names = FALSE, stringsAsFactors = FALSE,
  header = TRUE, data.table = FALSE
)
rownames(meth) <- meth$sample
meth <- meth[, -1]
meth <- as.data.frame(na.omit(meth))

probeOfInterest <- probe.features[which(probe.features$gene %in% rownames(immunomodulator)), ]
probeOfInterest <- probeOfInterest[intersect(rownames(probeOfInterest), rownames(meth)), ]
is.element(rownames(immunomodulator), probeOfInterest$gene)

meth <- meth[rownames(probeOfInterest), ]
meth$gene <- probeOfInterest$gene
meth <- as.data.frame(apply(
  meth[, setdiff(colnames(meth), "gene")], 2,
  function(x) tapply(x, INDEX = factor(meth$gene), FUN = median, na.rm = TRUE)
))
meth <- meth %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var = "ID") %>%
  mutate(ID = str_sub(ID, 1, 12)) %>%
  distinct(ID, .keep_all = TRUE) %>%
  column_to_rownames(var = "ID") %>%
  t() %>% as.data.frame()

cna <- read.table(
  "r0.data/TCGA.BRCA.sampleMap_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz",
  sep = "\t", row.names = 1, check.names = FALSE,
  stringsAsFactors = FALSE, header = TRUE
)
cna$gene <- sapply(strsplit(rownames(cna), "|", fixed = TRUE), "[", 1)
cna <- cna[!duplicated(cna$gene), ]
cna <- cna[, setdiff(colnames(cna), "gene")]
is.element(rownames(immunomodulator), rownames(cna))

cna <- cna[intersect(rownames(cna), rownames(immunomodulator)), ]
cna <- cna %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var = "ID") %>%
  mutate(ID = str_sub(ID, 1, 12)) %>%
  distinct(ID, .keep_all = TRUE) %>%
  column_to_rownames(var = "ID") %>%
  t() %>% as.data.frame()
cna[cna > 1] <- 1
cna[cna < -1] <- -1

comsam <- Reduce(intersect, list(colnames(expr), colnames(meth), colnames(cna)), rownames(sinfo))
sinfo <- sinfo[comsam, , drop = FALSE]
expr <- expr[, comsam]
meth <- meth[, comsam]
cna <- cna[, comsam]
immunomodulator <- immunomodulator[rownames(expr), ]

ICI.genes <- rownames(immunomodulator)
dat.ICI <- eset[ICI.genes, ] %>% na.omit() %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var = "ID") %>%
  merge(pdata_cmoic, ., by = "ID")
dim(dat.ICI)

input_ICI <- dat.ICI %>% arrange(CMOIC)

dat.ICI.p <- data.frame()
for (i in colnames(input_ICI[, 3:ncol(input_ICI)])) {
  aa <- input_ICI[, c(i, "CMOIC")]
  colnames(aa)[1] <- "var"
  aa <- kruskal.test(var ~ CMOIC, data = aa)
  aa <- data.frame(var = i, p = aa$p.value)
  dat.ICI.p <- rbind(dat.ICI.p, aa)
}
dat.ICI.p$p.label <- ifelse(
  dat.ICI.p$p < 0.0001, "****",
  ifelse(
    dat.ICI.p$p < 0.001, "***",
    ifelse(
      dat.ICI.p$p < 0.01, "**",
      ifelse(dat.ICI.p$p < 0.05, "*", "ns")
    )
  )
)

(n.subt <- length(unique(sinfo$subtype)))
subt <- unique(sinfo$subtype)

expMat <- as.data.frame(t(expr[rownames(immunomodulator), ]))
expMat$subtype <- sinfo[rownames(expMat), "subtype"]
expMat <- as.data.frame(t(apply(
  expMat[, setdiff(colnames(expMat), "subtype")], 2,
  function(x) tapply(
    x,
    INDEX = factor(expMat$subtype),
    FUN = median,
    na.rm = TRUE
  )
)))

corExpMeth <- ampFreq <- delFreq <- as.data.frame(matrix(
  NA,
  nrow = nrow(immunomodulator),
  ncol = n.subt,
  dimnames = list(
    rownames(immunomodulator),
    unique(sinfo$subtype)
  )
))

for (i in rownames(immunomodulator)) {
  if (!is.element(i, rownames(expr)) | !is.element(i, rownames(meth))) {
    corExpMeth[i, ] <- NA
  } else {
    for (j in subt) {
      sam <- rownames(sinfo[which(sinfo$subtype == j), , drop = FALSE])
      expr.subset <- as.numeric(expr[i, sam])
      meth.subset <- as.numeric(meth[i, sam])
      ct <- cor.test(expr.subset, meth.subset, method = "spearman")
      corExpMeth[i, j] <- ct$estimate
    }
  }
}

for (i in rownames(immunomodulator)) {
  if (!is.element(i, rownames(cna))) {
    ampFreq[i, ] <- NA
    delFreq[i, ] <- NA
  } else {
    ampFreqInAll <- sum(as.numeric(cna[i, ]) == 1) / ncol(cna)
    delFreqInAll <- sum(as.numeric(cna[i, ]) == -1) / ncol(cna)
    for (j in subt) {
      sam <- rownames(sinfo[which(sinfo$subtype == j), , drop = FALSE])
      cna.subset <- cna[, sam]
      ampFreqInSubt <- sum(as.numeric(cna.subset[i, ]) == 1) / length(sam)
      delFreqInSubt <- sum(as.numeric(cna.subset[i, ]) == -1) / length(sam)

      ampFreqInDiff <- ampFreqInSubt - ampFreqInAll
      delFreqInDiff <- delFreqInSubt - delFreqInAll

      ampFreq[i, j] <- ampFreqInDiff
      delFreq[i, j] <- delFreqInDiff
    }
  }
}

annCol <- data.frame(
  subtype = subt,
  row.names = subt
)
annCol <- annCol[order(annCol$subtype), , drop = FALSE]
annColors <- list()
annColors[["subtype"]] <- c(
  "CS1" = ggsci::pal_lancet()(6)[1],
  "CS2" = ggsci::pal_lancet()(6)[2],
  "CS3" = ggsci::pal_lancet()(6)[3],
  "CS4" = ggsci::pal_lancet()(6)[4],
  "CS5" = ggsci::pal_lancet()(6)[5]
)
top_anno <- HeatmapAnnotation(
  df = annCol,
  col = annColors,
  gp = gpar(col = "grey80"),
  simple_anno_size = unit(3.5, "mm"),
  show_legend = FALSE,
  show_annotation_name = FALSE,
  border = FALSE
)

annRow <- immunomodulator
annRow[which(annRow$Category == "Co-stimulator"), "Category"] <- "Co-stm"
annRow[which(annRow$Category == "Co-inhibitor"), "Category"] <- "Co-ihb"
annRow[which(annRow$Category == "Cell adhesion"), "Category"] <- "Cell\nadhesion"
annRow[which(annRow$Category == "Antigen presentation"), "Category"] <- "Antigen\npresentation"
annRow$Category <- factor(
  annRow$Category,
  levels = c("Co-stm", "Co-ihb", "Ligand", "Receptor", "Cell\nadhesion", "Antigen\npresentation", "Other")
)
annRow$ICI <- factor(annRow$ICI, levels = c("Stimulatory", "Inhibitory", "N/A"))
annRowColors <- list("ICI" = c("Stimulatory" = "#E59E02", "Inhibitory" = "black", "N/A" = "#888888"))
left_anno <- HeatmapAnnotation(
  df = annRow[, "ICI", drop = FALSE],
  which = "row",
  gp = gpar(col = "grey80"),
  col = annRowColors,
  simple_anno_size = unit(3.5, "mm"),
  show_annotation_name = FALSE,
  border = FALSE
)
pvalue_col_fun <- colorRamp2(c(0, 5, 10), c("#f7bf95", "white", "#fe8ca1"))
right_anno <- rowAnnotation(
  pvalue = anno_simple(-log10(dat.ICI.p$p), col = pvalue_col_fun,
                       pch = dat.ICI.p$p.label)
)

col_expr <- colorRamp2(seq(min(na.omit(expMat)), max(na.omit(expMat)), length.out = 5), heatmap.BlWtRd)
hm.expr <- Heatmap(
  matrix = as.matrix(expMat),
  col = col_expr,
  border = NA,
  rect_gp = gpar(col = "grey80"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10),
  show_column_names = FALSE,
  column_names_side = "top",
  row_split = annRow$Category,
  top_annotation = top_anno,
  left_annotation = left_anno,
  name = "mRNA\nExpression",
  width = ncol(expMat) * unit(4, "mm"),
  height = nrow(expMat) * unit(3.5, "mm")
)

col_corExprMeth <- colorRamp2(seq(min(na.omit(corExpMeth)), max(na.omit(corExpMeth)), length.out = 5), heatmap.BlWtRd)
hm.corExprMeth <- Heatmap(
  matrix = as.matrix(corExpMeth),
  col = col_corExprMeth,
  border = NA,
  rect_gp = gpar(col = "grey80"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10),
  show_column_names = FALSE,
  column_names_side = "top",
  row_split = annRow$Category,
  row_title = NULL,
  top_annotation = top_anno,
  name = "Expression\nvs. Methylation",
  width = ncol(expMat) * unit(4, "mm"),
  height = nrow(expMat) * unit(3.5, "mm")
)

col_ampFreq <- colorRamp2(seq(min(na.omit(ampFreq)), max(na.omit(ampFreq)), length.out = 5), heatmap.BlWtRd)
hm.ampFreq <- Heatmap(
  matrix = as.matrix(ampFreq),
  col = col_ampFreq,
  border = NA,
  rect_gp = gpar(col = "grey80"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10),
  show_column_names = FALSE,
  column_names_side = "top",
  row_split = annRow$Category,
  row_title = NULL,
  top_annotation = top_anno,
  name = "Amplification\nFrequency",
  width = ncol(expMat) * unit(4, "mm"),
  height = nrow(expMat) * unit(3.5, "mm")
)

col_delFreq <- colorRamp2(seq(min(na.omit(delFreq)), max(na.omit(delFreq)), length.out = 5), heatmap.BlWtRd)
hm.delFreq <- Heatmap(
  matrix = as.matrix(delFreq),
  col = col_delFreq,
  border = NA,
  rect_gp = gpar(col = "grey70"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10),
  show_column_names = FALSE,
  column_names_side = "top",
  row_split = annRow$Category,
  row_title = NULL,
  top_annotation = top_anno,
  right_annotation = right_anno,
  name = "Deletion\nFrequency",
  width = ncol(expMat) * unit(4, "mm"),
  height = nrow(expMat) * unit(3.5, "mm")
)

lgd_pvalue <- Legend(
  title = "p-value", col_fun = pvalue_col_fun, at = c(0, 5, 10),
  labels = c("-log10(0)", "-log10(5)", "-log10(10)")
)
lgd_sig <- Legend(
  pch = c("****", "***", "**", "*", "ns"), type = "points",
  labels = c("< 0.0001", "< 0.001", "< 0.01", "< 0.05", "> 0.05")
)

draw(
  hm.expr + hm.corExprMeth + hm.ampFreq + hm.delFreq,
  heatmap_legend_side = "bottom",
  annotation_legend_list = list(lgd_pvalue, lgd_sig)
)

input_dat <- input_cell[, c("CMOIC", "Immune.Subtype")]

data_summary <- input_dat %>%
  group_by(CMOIC, Immune.Subtype) %>%
  summarise(Count = n()) %>%
  ungroup()

total_counts <- data_summary %>%
  group_by(CMOIC) %>%
  summarise(Total = sum(Count))

data_summary <- data_summary %>%
  left_join(total_counts, by = "CMOIC") %>%
  mutate(Percentage = Count / Total * 100)

ggplot(data_summary, aes(x = CMOIC, y = Percentage, fill = Immune.Subtype)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(11, "Set3")[1:5]) +
  geom_text(aes(label = round(Percentage, 1)), position = position_stack(vjust = 0.5), size = 8) +
  theme_minimal() +
  theme(
    legend.text = element_text(size = 22),
    legend.title = element_text(size = 22),
    legend.position = "bottom",
    axis.text.x = element_text(size = 24, colour = ggsci::pal_lancet()(5)),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 28)
  ) +
  labs(y = "Percentage", fill = "Immune Subtype")

load("../_Immunity_breast_TCGA/Immunity/Scores_83_Signatures.rda")
input <- Scores_83_Signatures %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "ID") %>%
  dplyr::filter(
    ID == "Tcell_receptors_score" | ID == "Bcell_receptors_score" |
      ID == "STAT1_score" | ID == "IFNG_score_21050467" | ID == "MHC1_21978456" |
      ID == "HER2_Immune_PCA_18006808" | ID == "ICR_SCORE" | ID == "APM1" |
      ID == "Troester_WoundSig_19887484" | ID == "MHC2_21978456"
  ) %>%
  column_to_rownames(var = "ID")
input[, 1:3]
rownames(input) <- c(
  "APM", "BCR", "STAT1", "IFNgamma", "TCR", "MHCI", "MHCII",
  "Wound_healing", "HER2", "ICR"
)
input_immune <- input %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var = "ID") %>%
  merge(pdata_cmoic, ., by = "ID")

data_summary <- input_immune %>%
  pivot_longer(., 3:ncol(input_immune), names_to = "Pathway", values_to = "Score") %>%
  group_by(Pathway, CMOIC) %>%
  summarise(mean_score = mean(Score, na.rm = TRUE)) %>%
  ungroup()
data_summary$Pathway <- factor(data_summary$Pathway, levels = sort(unique(data_summary$Pathway), decreasing = TRUE))

head(data_summary)
ggplot(data_summary, aes(x = CMOIC, y = Pathway, fill = mean_score)) +
  geom_tile(color = "white", size = 0.8) +
  geom_text(aes(label = round(mean_score, 2)), size = 4, color = "black") +
  scale_fill_gradientn(
    colors = c(RColorBrewer::brewer.pal(11, "Set3")[1], "white", RColorBrewer::brewer.pal(11, "Set3")[2]),
    values = scales::rescale(c(min(data_summary$mean_score), 0, 0.5, max(data_summary$mean_score))),
    limits = c(min(data_summary$mean_score), max(data_summary$mean_score)),
    breaks = c(min(data_summary$mean_score), max(data_summary$mean_score)),
    labels = c("Lower score", "Higher score"),
    name = "Mean score"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14, colour = ggsci::pal_lancet()(5)),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(size = 18, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  labs(title = "Immune signatures", fill = "Mean score") +
  guides(fill = guide_colorbar(barwidth = 10, barheight = 1, title.position = "top", title.hjust = 0.5))
