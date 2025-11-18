library("MOVICS")

mo.data <- dat.MOVICS[1:6]

count <- dat.MOVICS$count
tpm <- dat.MOVICS$mRNA.expr
maf <- dat.MOVICS$maf
segment <- dat.MOVICS$segment
surv.info <- dat.MOVICS$clin.info

dat.MOVICS$mRNA.expr <- dat.MOVICS$mRNA.expr %>% .[rowSums(.) > 1, ]
elite.mRNA <- getElites(
  dat = dat.MOVICS$mRNA.expr,
  method = "sd",
  elite.pct = 0.035
)
mo.data$mRNA.expr <- elite.mRNA$elite.dat
mo.data$mRNA.expr <- na.omit(mo.data$mRNA.expr)
dim(mo.data$mRNA.expr)

elite.miRNA <- getElites(
  dat = dat.MOVICS$miRNA.expr,
  method = "sd",
  elite.pct = 0.4
)
mo.data$miRNA.expr <- elite.miRNA$elite.dat
mo.data$miRNA.expr <- na.omit(mo.data$miRNA.expr)
dim(mo.data$miRNA.expr)

elite.lncRNA <- getElites(
  dat = dat.MOVICS$lncRNA.expr,
  method = "sd",
  elite.pct = 0.1
)
mo.data$lncRNA.expr <- elite.lncRNA$elite.dat
mo.data$lncRNA.expr <- na.omit(mo.data$lncRNA.expr)
dim(mo.data$lncRNA.expr)

elite.meth <- getElites(
  dat = dat.MOVICS$meth.beta,
  method = "sd",
  elite.pct = 0.06
)
mo.data$meth.beta <- elite.meth$elite.dat
mo.data$meth.beta <- na.omit(mo.data$meth.beta)
dim(mo.data$meth.beta)

elite.CNV <- getElites(
  dat = dat.MOVICS$CNV.GISTIC,
  method = "sd",
  elite.pct = 0.04
)
mo.data$CNV.GISTIC <- elite.CNV$elite.dat
mo.data$CNV.GISTIC <- na.omit(mo.data$CNV.GISTIC)
dim(mo.data$CNV.GISTIC)

elite.mut <- getElites(
  dat = dat.MOVICS$mut.status,
  method = "freq",
  elite.num = 30,
  elite.pct = 0.1
)
mo.data$mut.status <- elite.mut$elite.dat
dim(mo.data$mut.status)

optk.num <- getClustNum(
  data = mo.data,
  is.binary = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
  try.N.clust = 2:8,
  fig.path = "r1.MOVICS/",
  fig.name = "CLUSTER NUMBER OF TCGA-BRCA"
)

dat.cluster <- getMOIC(
  data = mo.data,
  methodslist = list(
    "iClusterBayes",
    "SNF",
    "PINSPlus",
    "NEMO",
    "COCA",
    "LRAcluster",
    "ConsensusClustering",
    "IntNMF",
    "CIMLR",
    "MoCluster"
  ),
  N.clust = 5,
  type = c("gaussian", "gaussian", "gaussian", "gaussian", "gaussian", "binomial")
)

dat.cmoic <- getConsensusMOIC(
  moic.res.list = dat.cluster,
  fig.name = "CONSENSUS HEATMAP",
  distance = "euclidean",
  linkage = "ward.D"
)

getSilhouette(
  sil = dat.cmoic$sil,
  fig.name = "SILHOUETTE",
  height = 5.5,
  width = 5
)

dat.surv <- compSurv(
  moic.res = dat.cmoic,
  surv.info = surv.info,
  convt.time = "y",
  surv.median.line = "h",
  xyrs.est = c(5, 10),
  fig.name = "KAPLAN-MEIER CURVE OF CONSENSUSMOIC"
)
print(dat.surv)

mo.data$meth.beta <- log2(mo.data$meth.beta / (1 - mo.data$meth.beta))

plotdata <- getStdiz(
  data = mo.data,
  halfwidth = c(2, 2, 2, 2, 2, NA),
  centerFlag = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE),
  scaleFlag = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)
)

feat <- dat.cluster$MoCluster$feat.res
feat1 <- feat[which(feat$dataset == "mRNA.expr"), ][1:10, "feature"]
feat2 <- feat[which(feat$dataset == "miRNA.expr"), ][1:10, "feature"]
feat3 <- feat[which(feat$dataset == "lncRNA.expr"), ][1:10, "feature"]
feat4 <- feat[which(feat$dataset == "meth.beta"), ][1:10, "feature"]
feat5 <- feat[which(feat$dataset == "CNV.GISTIC"), ][1:10, "feature"]
feat6 <- feat[which(feat$dataset == "mut.status"), ][1:10, "feature"]
annRow <- list(feat1, feat2, feat3, feat4, feat5, feat6)

mRNA.col <- c(col_a[1], "white", col_a[2])
miRNA.col <- c("#6699CC", "white", "#FF3C38")
lncRNA.col <- c("#6699CC", "white", "#FF3C38")
meth.col <- c("#0074FE", "#96EBF9", "#FEE900", "#F00003")
CNV.col <- c("#8FB2C9", "white", "#F0A1A8")
mut.col <- c("grey90", "black")
col.list <- list(mRNA.col, miRNA.col, lncRNA.col, meth.col, CNV.col, mut.col)

getMoHeatmap(
  data = plotdata,
  row.title = c("mRNA", "miRNA", "lncRNA", "Methylation", "CNV", "Mutation"),
  is.binary = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
  legend.name = c("mRNA.TPM", "miRNA.TPM", "lncRNA.TPM", "M value", "CNV", "Mutated"),
  clust.res = dat.cmoic$clust.res,
  clust.dend = NULL,
  show.rownames = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
  show.colnames = FALSE,
  color = col.list,
  annCol = NULL,
  annColors = NULL,
  width = 10,
  height = 5,
  fig.name = "COMPREHENSIVE HEATMAP"
)

annCol <- surv.info[, c("PAM50", "Stage", "Age"), drop = FALSE]

annColors <- list(
  Age = circlize::colorRamp2(
    breaks = c(
      min(annCol$Age, na.rm = TRUE),
      median(annCol$Age, na.rm = TRUE),
      max(annCol$Age, na.rm = TRUE)
    ),
    colors = c("#0000AA", "#555555", "#AAAA00")
  ),
  PAM50 = c(
    Basal = "blue",
    Her2 = "red",
    LumA = "yellow",
    LumB = "green",
    Normal = "black"
  ),
  Stage = c(
    "Stage I" = "green",
    "Stage II" = "blue",
    "Stage III" = "red",
    "Stage IV" = "yellow"
  )
)

getMoHeatmap(
  data = plotdata,
  row.title = c("mRNA", "miRNA", "lncRNA", "Methylation", "CNV", "Mutation"),
  is.binary = c(FALSE, FALSE, FALSE, FALSE, TRUE),
  legend.name = c("mRNA.TPM", "miRNA.TPM", "lncRNA.TPM", "M value", "CNV", "Mutated"),
  clust.res = dat.cmoic$clust.res,
  clust.dend = NULL,
  show.rownames = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
  show.colnames = FALSE,
  show.row.dend = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
  annRow = annRow,
  color = col.list,
  annCol = annCol,
  annColors = annColors,
  width = 10,
  height = 5,
  fig.name = "COMPREHENSIVE HEATMAP OF CONSENSUSMOIC"
)

dat.surv <- compSurv(
  moic.res = dat.cmoic,
  surv.info = surv.info,
  convt.time = "m",
  surv.median.line = "h",
  xyrs.est = c(5, 10),
  fig.name = "KAPLAN-MEIER CURVE OF CONSENSUSMOIC"
)
print(dat.surv)

com.clin <- compClinvar(
  moic.res = dat.cmoic,
  var2comp = surv.info,
  strata = "Subtype",
  factorVars = c("PAM50", "Stage", "fustat"),
  nonnormalVars = "futime",
  exactVars = "Stage",
  doWord = TRUE,
  tab.name = "SUMMARIZATION OF CLINICAL FEATURES"
)
print(com.clin$compTab)

com.mut <- compMut(
  moic.res = dat.cmoic,
  mut.matrix = dat.MOVICS$mut.status,
  doWord = TRUE,
  doPlot = TRUE,
  freq.cutoff = 0.05,
  p.adj.cutoff = 0.05,
  innerclust = TRUE,
  annCol = annCol,
  annColors = annColors,
  width = 6,
  height = 2,
  fig.name = "ONCOPRINT FOR SIGNIFICANT MUTATIONS",
  tab.name = "INDEPENDENT TEST BETWEEN SUBTYPE AND MUTATION"
)
print(com.mut)

head(maf)
maf2 <- maf %>%
  mutate(Tumor_Sample_Barcode = str_sub(Tumor_Sample_Barcode, 1, 12))

com.tmb <- compTMB(
  moic.res = dat.cmoic,
  maf = maf2,
  rmDup = TRUE,
  rmFLAGS = FALSE,
  exome.size = 38,
  test.method = "nonparametric",
  fig.name = "DISTRIBUTION OF TMB AND TITV"
)
head(com.tmb$TMB.dat)

colnames(segment) <- c("sample", "chrom", "start", "end", "value")
head(segment)

com.fga <- compFGA(
  moic.res = dat.cmoic,
  segment = segment,
  iscopynumber = TRUE,
  cnathreshold = 0.3,
  test.method = "nonparametric",
  fig.name = "BARPLOT OF FGA"
)
head(com.fga$summary)

com.drug <- compDrugsen(
  moic.res = dat.cmoic,
  norm.expr = mo.data[[1]][, dat.cmoic$clust.res$samID],
  drugs = c("Cisplatin", "Paclitaxel"),
  tissueType = "breast",
  test.method = "nonparametric",
  prefix = "BOXVIOLIN OF ESTIMATED IC50"
)
head(com.drug$Cisplatin)

surv.info$Stage <- factor(
  surv.info$Stage,
  levels = c("Stage I", "Stage II", "Stage III", "Stage IV")
)

com.agree <- compAgree(
  moic.res = dat.cmoic,
  subt2comp = surv.info[, c("PAM50", "Stage")],
  doPlot = TRUE,
  box.width = 0.2,
  fig.name = "AGREEMENT OF CONSENSUSMOIC WITH PAM50 AND PSTAGE"
)
print(com.agree)

runDEA(
  dea.method = "deseq2",
  expr = count,
  moic.res = dat.cmoic,
  prefix = "TCGA-BRCA"
)

marker.up <- runMarker(
  moic.res = dat.cmoic,
  dea.method = "deseq2",
  prefix = "TCGA-BRCA",
  dat.path = getwd(),
  res.path = getwd(),
  p.cutoff = 0.05,
  p.adj.cutoff = 0.05,
  dirct = "up",
  n.marker = 100,
  doplot = TRUE,
  norm.expr = fpkm,
  annCol = annCol,
  annColors = annColors,
  show_rownames = FALSE,
  fig.name = "UPREGULATED BIOMARKER HEATMAP"
)
head(marker.up$templates)

marker.dn <- runMarker(
  moic.res = cmoic.brca,
  dea.method = "limma",
  prefix = "TCGA-BRCA",
  dirct = "down",
  n.marker = 50,
  doplot = TRUE,
  norm.expr = fpkm,
  annCol = annCol,
  annColors = annColors,
  fig.name = "DOWNREGULATED BIOMARKER HEATMAP"
)

MSIGDB.FILE <- system.file(
  "extdata",
  "c5.bp.v7.1.symbols.xls",
  package = "MOVICS",
  mustWork = TRUE
)

gsea.up <- runGSEA(
  moic.res = cmoic.brca,
  dea.method = "deseq",
  prefix = "TCGA-BRCA",
  dat.path = getwd(),
  res.path = getwd(),
  msigdb.path = MSIGDB.FILE,
  norm.expr = fpkm,
  dirct = "up",
  p.cutoff = 0.05,
  p.adj.cutoff = 0.25,
  gsva.method = "gsva",
  norm.method = "mean",
  fig.name = "UPREGULATED PATHWAY HEATMAP"
)

print(gsea.up$gsea.list$CS1[1:6, 3:6])
head(round(gsea.up$grouped.es, 3))

gsea.dn <- runGSEA(
  moic.res = cmoic.brca,
  dea.method = "deseq2",
  prefix = "TCGA-BRCA",
  msigdb.path = MSIGDB.FILE,
  norm.expr = fpkm,
  dirct = "down",
  p.cutoff = 0.05,
  p.adj.cutoff = 0.25,
  gsva.method = "gsva",
  norm.method = "mean",
  fig.name = "DOWNREGULATED PATHWAY HEATMAP"
)

GSET.FILE <- system.file(
  "extdata",
  "gene sets of interest.gmt",
  package = "MOVICS",
  mustWork = TRUE
)

gsva.res <- runGSVA(
  moic.res = cmoic.brca,
  norm.expr = fpkm,
  gset.gmt.path = GSET.FILE,
  gsva.method = "gsva",
  annCol = annCol,
  annColors = annColors,
  fig.path = getwd(),
  fig.name = "GENE SETS OF INTEREST HEATMAP",
  height = 5,
  width = 8
)

print(gsva.res$raw.es[1:3, 1:3])
print(gsva.res$scaled.es[1:3, 1:3])
