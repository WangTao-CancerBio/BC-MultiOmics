
library(Seurat)
library(edgeR)
library(vcd)
library(pheatmap)
library(scRNAtoolVis)
library(infercnv)


(load("CombSeurat_GSE161529.rda"))
head(SampleStats)
SamplesComb


Pat <- c("N_PM0019_Total","N_0093_total","N_PM0095_Total","N_PM0233_Total",
         "HER2_AH0308","TN-0126","ER_MH0025","ER_MH0032",
         "HER2_MH0176","ER_MH0042","HER2_PM0337","TN_B1_MH4031"
)
SamplesComb.tar <- SamplesComb[!str_detect(SamplesComb,"^N_.*")]
SamplesComb.tar <- SamplesComb.tar[!str_detect(SamplesComb.tar,"_LN$")]
SamplesComb.tar

dimUsed <- 30
Anchors <- FindIntegrationAnchors(object.list=CombSeurat[SamplesComb.tar], dims=1:dimUsed,
                                  anchor.features=1000, scale=TRUE, k.anchor=5, k.filter=30,
                                  k.score=90, max.features=100)
so <- IntegrateData(anchorset=Anchors, dims=1:dimUsed, k.weight=100)
DefaultAssay(so) <- "integrated"
so <- ScaleData(so, verbose=FALSE)
so <- RunPCA(so, npcs=dimUsed, verbose=FALSE)
so <- RunTSNE(so, dims=1:dimUsed, seed.use=2018)
so <- RunUMAP(so, dims=1:dimUsed, seed.use=2018)
TumLN <- so


resolution <- 0.7

TumLN <- FindNeighbors(TumLN, dims=1:dimUsed, verbose=FALSE)
TumLN <- FindClusters(TumLN, resolution=resolution, verbose=FALSE)

DimPlot(TumLN, group.by=c("seurat_clusters"),label = T,repel = T,reduction = "tsne")+ NoLegend()
ggsave("sc.pdf",width = 6,height = 6)


srat <- TumLN



table(srat$seurat_clusters)
DefaultAssay(srat) <- "RNA"
Idents(srat) <- "seurat_clusters"

marker_genes.All <- c(c("EPCAM","KRT19"),
                      c("COL1A1","PDGFRA"),
                      c("VWF","CD36"),
                      c("CD3E","CD8A","CD4","IL7R"),
                      c("MS4A1", "CD79A","CD19","JCHAIN"),
                      c("ITGAX","ITGAM","CD68","CD14","CD86","CD163"),
                      c("GNLY", "NKG7", "HAVR2C","NCAM1"),
                      c("RGS5","CSPG4","PDGFRB","MCAM"),
                      c("CD14", "FCGR3A", "LYZ"),
                      c("FCER1A", "PPBP"),
                      c("HBB","HBA2"),
                      c("NKG7","TPSAB1"),
                      c("PTPRC")
)
scCustomize::Stacked_VlnPlot(srat,features = marker_genes.All,group.by = "seurat_clusters")

srat@meta.data$cell <- plyr::mapvalues(
          from = c(
                    0,3,4,2,6,8,13,14,
                    7,15,
                    12,
                    5,18,19,
                    1,9,17,
                    11,
                    10,16
          ),
          to = c(rep("Epithelial cells",8),
                 rep("Fibroblasts",2),
                 rep("Ednothelial cells",1),
                 rep("Macrophages",3),
                 rep("T cells",3),
                 rep("Plasma cells",1),
                 rep("Mast cells",2)
          ),
          x = srat@meta.data$seurat_clusters)
table(srat$cell)
DimPlot(srat, group.by=c("seurat_clusters","cell"),label.box = T,label = T,repel = F,reduction = "tsne")+ NoLegend()
Idents(srat) <- "cell"
all_markers <- FindAllMarkers(object = srat)
clustGroup <- table(srat$cell) %>% names() %>% as.numeric() %>% na.omit()
for(i in clustGroup){
          dat.markers <- subset(all_markers,cluster==i) %>%
                    arrange(desc(avg_log2FC))
          write.csv(dat.markers,paste0("cellMarkers_",i,".csv"))
}

DimPlot(srat, group.by=c("cell","group"),label = T,repel = T)

marker_genes.all <- c("EPCAM","KRT19",
                      "COL1A1",
                      "CD36",
                      c("CD3E","CD8A","CD4","IL7R"),
                      c("GATA2"),
                      "CD68","CD14",
                      c("JCHAIN","MZB1")
)
col.P2 <- c(RColorBrewer::brewer.pal(9,"Set1"),RColorBrewer::brewer.pal(12,"Set3"))
srat$cell <- factor(srat$cell,levels = c("Epithelial cells","Fibroblasts","Ednothelial cells",
                                         "T cells","Mast cells","Macrophages",
                                         "Plasma cells"))
scCustomize::Stacked_VlnPlot(srat,features = marker_genes.all,group.by = "cell",colors_use = col.P2,x_lab_rotate = 45)
ggsave(file="cellMarkers.pdf",height = 10, width = 10)


DefaultAssay(srat) <- "RNA"

DimPlot(
  srat,
  group.by = c("cell"),
  label = TRUE,
  raster = FALSE,
  label.box = TRUE,
  repel = TRUE,
  cols = col.P2,
  reduction = "tsne"
) + NoLegend() + labs(title = "Cell types")

ggsave("celltypes.png", width = 8, height = 8)

DimPlot(
  srat,
  group.by = c("seurat_clusters"),
  label = TRUE,
  raster = FALSE,
  label.box = TRUE,
  repel = TRUE,
  cols = c(col.P2, "darkblue"),
  reduction = "tsne"
) + NoLegend() + labs(title = "Clusters")

ggsave("seurat_clusters.png", width = 8, height = 8)

DefaultAssay(srat) <- "RNA"

scCustomize::Plot_Density_Custom(
  seurat_object = srat,
  features = marker_genes.all,
  viridis_palette = "plasma",
  reduction = "tsne"
) +
  theme(
    title = element_text(size = 14, face = "bold"),
    axis.line.x = element_line(linewidth = 0.5),
    axis.line.y = element_line(linewidth = 0.5)
  )

ggsave("markers.png", width = 17, height = 11)

library(infercnv)
table(srat$cell)

srat.cnv <- subset(
  srat,
  cell == "Epithelial cells" |
    cell == "T cells" |
    cell == "Macrophages"
)

dim(srat.cnv)
srat.cnv$epi <- as.character(srat.cnv$cell)
table(srat.cnv$epi)

dat.interCNV.matrix <- as.matrix(srat.cnv@assays$RNA@counts)

library(AnnoProbe)
geneInfor <- annoGene(rownames(dat.interCNV.matrix), "SYMBOL", "human")
geneInfor <- geneInfor[with(geneInfor, order(chr, start)), c(1, 4:6)]
geneInfor <- geneInfor[!duplicated(geneInfor[, 1]), ]
geneInfor <- geneInfor %>%
  mutate(sort = as.numeric(str_extract(chr, "\\d.*"))) %>%
  arrange(sort) %>%
  select(-sort)

write.table(
  geneInfor,
  file = "geneInfor.txt",
  sep = "\t",
  quote = FALSE,
  col.names = FALSE,
  row.names = FALSE
)

dat.interCNV.anno <- data.frame(
  colnames(srat.cnv),
  srat.cnv$epi
) %>%
  remove_rownames()

write.table(
  dat.interCNV.anno,
  "dat.interCNV.anno.txt",
  sep = "\t",
  quote = FALSE,
  col.names = FALSE,
  row.names = FALSE
)

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = dat.interCNV.matrix,
  annotations_file = "dat.interCNV.anno.txt",
  delim = "\t",
  gene_order_file = "geneInfor.txt",
  ref_group_names = c("T cells", "Macrophages")
)

rm(dat.interCNV.matrix)
options(scipen = 100)
gc()
options(bitmapType = "Xlib")

infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = "interCNV",
  cluster_by_groups = TRUE,
  denoise = TRUE,
  num_threads = 60,
  HMM = FALSE
)

options(bitmapType = "cairo")

infercnv::plot_cnv(
  infercnv_obj,
  plot_chr_scale = TRUE,
  output_filename = "better_plot",
  output_format = "pdf",
  custom_color_pal = color.palette(
    c("#8DD3C7", "white", "#BC80BD"),
    c(2, 2)
  )
)

infer_CNV_obj <- readRDS("run.final.infercnv_obj")
expr <- infer_CNV_obj@expr.data
data_cnv <- as.data.frame(expr)

tmp1 <- expr[, infer_CNV_obj@reference_grouped_cell_indices$`T cells`]
tmp2 <- expr[, infer_CNV_obj@reference_grouped_cell_indices$`Macrophages`]
tmp <- cbind(tmp1, tmp2)
down <- mean(rowMeans(tmp)) - 2 * mean(apply(tmp, 1, sd))
up <- mean(rowMeans(tmp)) + 2 * mean(apply(tmp, 1, sd))
oneCopy <- up - down
a1 <- down - 2 * oneCopy
a2 <- down - 1 * oneCopy
a3 <- up + 1 * oneCopy
a4 <- up + 2 * oneCopy

cnv_score_table <- infer_CNV_obj@expr.data
cnv_score_mat <- as.matrix(cnv_score_table)

cnv_score_table[cnv_score_mat > 0 & cnv_score_mat < a2] <- "A"
cnv_score_table[cnv_score_mat >= a2 & cnv_score_mat < down] <- "B"
cnv_score_table[cnv_score_mat >= down & cnv_score_mat < up] <- "C"
cnv_score_table[cnv_score_mat >= up & cnv_score_mat <= a3] <- "D"
cnv_score_table[cnv_score_mat > a3 & cnv_score_mat <= a4] <- "E"
cnv_score_table[cnv_score_mat > a4] <- "F"

cnv_score_table_pts <- cnv_score_mat

cnv_score_table_pts[cnv_score_table == "A"] <- 2
cnv_score_table_pts[cnv_score_table == "B"] <- 1
cnv_score_table_pts[cnv_score_table == "C"] <- 0
cnv_score_table_pts[cnv_score_table == "D"] <- 1
cnv_score_table_pts[cnv_score_table == "E"] <- 2
cnv_score_table_pts[cnv_score_table == "F"] <- 2

cell_scores_CNV <- as.data.frame(colSums(cnv_score_table_pts))
colnames(cell_scores_CNV) <- "cnv_score"

meta <- srat.cnv@meta.data[,c("cell","seurat_clusters")] %>%
        rownames_to_column(var = "CB")
set.seed(20210418)
kmeans.result <- kmeans(t(infercnv_obj@expr.data), 10)
kmeans_df <- data.frame(kmeans_class=kmeans.result$cluster)
kmeans_df$kmeans_class <- paste0("epi_",kmeans_df$kmeans_class)
kmeans_df <- kmeans_df %>% 
        rownames_to_column(var = "CB") %>%
        inner_join(meta,by="CB")
cell_scores_CNV <- cell_scores_CNV %>%
        rownames_to_column(var = "CB") %>%
        merge(.,kmeans_df,by="CB")
cell_scores_CNV$cell <- as.character(cell_scores_CNV$cell)
cell_scores_CNV$kmeans_class2 <- as.character(cell_scores_CNV$kmeans_class)
cell_scores_CNV$kmeans_class2[cell_scores_CNV$cell=="T cells" | cell_scores_CNV$cell=="Macrophages"] <- cell_scores_CNV$cell[cell_scores_CNV$cell=="T cells" | cell_scores_CNV$cell=="Macrophages"]
table(cell_scores_CNV$kmeans_class,cell_scores_CNV$kmeans_class2)
ggviolin(cell_scores_CNV, "kmeans_class2", "cnv_score", fill = "kmeans_class2", add = "boxplot",
         palette = ggsci::pal_frontiers()(12),
         order = c("T cells", "Macrophages", paste0("epi_", 1:10))) +
        labs(x = NULL, y = "CNV score (log10)") +
        geom_hline(yintercept = median(cell_scores_CNV$cnv_score[cell_scores_CNV$kmeans_class2 == "T cells"]),
                   linetype = "dashed") +
        NoLegend() +
        scale_y_continuous(
                trans = scales::log10_trans(),
                breaks = scales::trans_breaks("log10", function(x) 10**x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(1, 10000)
        ) +
        coord_cartesian(ylim = c(1, 10000)) +
        theme(axis.text.x = element_text(size = 16,angle = 45,hjust = 1))
ggsave("CNV_level.pdf",width = 8,height = 6)
rm(infer_CNV_obj,infercnv_obj,obs,TumLN)
tumor.id <- c("epi_1","epi_2","epi_3","epi_4","epi_5","epi_6","epi_7","epi_8","epi_9")
cell_scores_CNV$tumor <- ifelse(cell_scores_CNV$kmeans_class2 %in% tumor.id,"tumor","normal")
table(cell_scores_CNV$tumor,cell_scores_CNV$cell)
tumor.id <- cell_scores_CNV$CB[cell_scores_CNV$tumor=="tumor"]
srat.tumor <- srat.cnv[,tumor.id]
DimPlot(srat.tumor, group.by=c("cell"),label.box = TRUE,label = TRUE,repel = FALSE,reduction = "tsne")+ NoLegend()
dim(srat.cnv);dim(srat.tumor)
table(srat.tumor$cell,srat.tumor$epi)
save(SamplesComb.tar,Anchors,resolution,so,srat,srat.cnv,srat.tumor,file = "sratDATA.rda")
rm(Pat,SamplesComb.tar,Anchors,TumLN,resolution,so)
rm(CombSeurat)
rm(srat,srat.cnv)
expr.sc <- as.matrix(srat.tumor@assays$RNA$data)
sc.ntp.pred <- runNTP(expr       = expr.sc,
                       templates  = marker.up$templates,
                       scaleFlag  = TRUE,
                       centerFlag = TRUE,
                       doPlot     = TRUE,
                       fig.name   = "NTP_HEATMAP_FOR_sc") 
rm(expr.sc)
sc.ntp.pred.P <- sc.ntp.pred$ntp.res %>% dplyr::filter(FDR<0.05)
srat.tumor.clust <- srat.tumor[,rownames(sc.ntp.pred.P)]
srat.tumor.clust <- AddMetaData(srat.tumor.clust, sc.ntp.pred.P$prediction, col.name ="CS")
srat.tumor.clust$CS
DimPlot(srat.tumor.clust, group.by=c("CS"),label.box = TRUE,label = TRUE,repel = FALSE,reduction = "tsne",cols=ggsci::pal_lancet()(6))+ NoLegend()
ggsave("cs.png",width=8,height=8)
DimPlot(srat.tumor.clust, group.by=c("CS"),label = TRUE,raster=FALSE,label.box = TRUE,repel = TRUE,cols=ggsci::pal_lancet()(6),label.size = 8,reduction = 'tsne')+NoLegend()+labs(title = "Cell subtype")
ggsave("seurat_cs.png",width = 8,height = 8)
DimPlot(srat.tumor.clust, group.by=c("CS"),repel = TRUE,reduction = "tsne",split.by="CS",ncol=5,
        cols = col.P2)+NoLegend()+labs(title = "")+
        theme(axis.line = element_blank(),
              axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank())
ggsave("scsubtype.pdf",width = 20,height = 4)
Idents(srat.tumor.clust) <- "CS"
srat.markers <- FindAllMarkers(object = srat.tumor.clust)
jjVolcano(diffData = srat.markers,legend.position = c(.9,.1))+xlab("Cell subtype")
ggsave("volMarkers.png",width = 8,height = 6)
