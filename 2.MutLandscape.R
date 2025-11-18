library(ComplexHeatmap)
library(RColorBrewer)
library(deconstructSigs)
library(BSgenome.Hsapiens.UCSC.hg38)
library(NMF)
require(maftools)

load("TCGA-BRCA_maf.rda")

maf <- as.data.frame(maf[, c(1:2, 4:13, 16)])
rmSilence <- FALSE

if (rmSilence) {
  maf <- as.data.frame(maf[which(maf$Variant_Type == "SNP" & maf$Variant_Classification != "Silent"), ])
} else {
  maf <- as.data.frame(maf[which(maf$Variant_Type == "SNP"), ])
}

snp.count <- mut.to.sigs.input(
  mut.ref  = maf,
  sample.id = "Tumor_Sample_Barcode",
  chr       = "Chromosome",
  pos       = "Start_Position",
  ref       = "Reference_Allele",
  alt       = "Tumor_Seq_Allele2",
  bsg       = BSgenome.Hsapiens.UCSC.hg38
)

cut.off <- 0.06
mut.wt <- data.frame()
sigs.out.list <- list()
index <- 1

cosmic.sbs <- read.table(
  "COSMIC_SBS.txt",
  header    = TRUE,
  sep       = "\t",
  row.names = 1
) |>
  t() |>
  as.data.frame()

for (sample in rownames(snp.count)) {
  cat(paste0(sample, " starts and ", length(rownames(snp.count)) - index, " samples remain to be analyzed!\n"))
  tmp <- whichSignatures(
    tumor.ref        = snp.count,
    signatures.ref   = cosmic.sbs,
    sample.id        = sample,
    contexts.needed  = TRUE,
    tri.counts.method = "exome2genome",
    signature.cutoff = cut.off
  )
  index <- index + 1
  sigs.out.list[[sample]] <- tmp
  tmp <- data.frame(c(tmp$weights, unknown = tmp$unknown), row.names = sample)
  mut.wt <- rbind.data.frame(mut.wt, tmp)
}

mutsig <- mut.wt |>
  rownames_to_column(var = "id") |>
  mutate(id = str_sub(id, 1, 12)) |>
  distinct(id, .keep_all = TRUE) |>
  column_to_rownames(var = "id")

cna.region <- read.table(
  "20160128-BRCA-all_lesions.conf_99.txt",
  sep             = "\t",
  row.names       = 1,
  check.names     = FALSE,
  stringsAsFactors = FALSE,
  header          = TRUE
)
colnames(cna.region) <- str_sub(colnames(cna.region), 1, 12)

cna.gene <- read.table(
  "20160128-BRCA-all_thresholded.by_genes.txt",
  sep             = "\t",
  row.names       = 1,
  check.names     = FALSE,
  stringsAsFactors = FALSE,
  header          = TRUE
)
colnames(cna.gene) <- str_sub(colnames(cna.gene), 1, 12)

mut <- dat.MOVICS$mut.status

tmb <- com.tmb$TMB.dat[, c(1, 4)] |>
  remove_rownames() |>
  column_to_rownames(var = "samID")

subt <- dat.cmoic$clust.res |>
  dplyr::select(-samID)
subt$clust <- paste0("CS", subt$clust)

coPatID.genom <- Reduce(
  intersect,
  list(
    rownames(mutsig),
    colnames(mut),
    rownames(tmb),
    rownames(subt)
  )
)

mutsig <- mutsig[coPatID.genom, ]
mut    <- mut[, coPatID.genom]
tmb    <- tmb[coPatID.genom, ] |> as.data.frame()
rownames(tmb) <- coPatID.genom
colnames(tmb) <- "log10TMB"

subt <- subt[coPatID.genom, ] |> as.data.frame()
rownames(subt) <- coPatID.genom
colnames(subt) <- "CMOIC"

clust.col <- ggsci::pal_lancet()(9)[1:5]
blue <- "#5bc0eb"
red  <- "#f25f5c"

ha <- HeatmapAnnotation(
  Group = mutsig$CMOIC,
  col   = list(
    Group = c(
      CS1 = clust.col[1],
      CS2 = clust.col[2],
      CS3 = clust.col[3],
      CS4 = clust.col[4],
      CS5 = clust.col[5]
    )
  ),
  annotation_name_gp    = gpar(fontsize = 25, col = "red"),
  na_col                = "white",
  show_legend           = rep(TRUE, 5),
  annotation_height     = unit(rep(4, 6), "mm"),
  annotation_legend_param = list(
    Group = list(title = "Group")
  ),
  border = TRUE
)

col_fun <- circlize::colorRamp2(c(0, 1), c("white", red))

pheatmap(
  t(mutsig[, -ncol(mutsig)]),
  col             = col_fun,
  top_annotation  = ha,
  show_rownames   = TRUE,
  show_colnames   = FALSE,
  cluster_rows    = TRUE,
  cluster_cols    = FALSE,
  show_column_dend = FALSE,
  show_row_dend    = FALSE
)

mutsig.per <- mutsig
mutsig.per[mutsig.per != 0] <- 1
mutsig.per <- data.frame(
  name = colnames(mutsig.per),
  per  = colSums(mutsig.per) / nrow(mutsig.per)
)

mutsig <- mutsig[, c("SBS2", "SBS13", "SBS7b", "SBS7d")]
mutsig$APOBEC <- mutsig$SBS2 + mutsig$SBS13
mutsig$CMOIC  <- subt[rownames(mutsig), "CMOIC"]
mutsig <- mutsig[order(mutsig$CMOIC, -mutsig$APOBEC, decreasing = FALSE), ]

mutgene <- com.mut$`Gene (Mutated)`

onco.input <- mut[mutgene, rownames(mutsig)]
onco.input[onco.input == 1]      <- "Mutated"
onco.input[onco.input != "Mutated"] <- ""

alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(
      x, y,
      w - unit(0.5, "mm"),
      h - unit(0.5, "mm"),
      gp = gpar(fill = "#dcddde", col = "#dcddde")
    )
  },
  Mutated = function(x, y, w, h) {
    grid.rect(
      x, y,
      w - unit(0.5, "mm"),
      h - unit(0.5, "mm"),
      gp = gpar(fill = "#A60000", col = "#A60000")
    )
  }
)

col <- c(Mutated = "#A60000")

my_ann <- subt[rownames(mutsig), , drop = FALSE] |> na.omit()

my_annotation <- HeatmapAnnotation(
  df  = my_ann,
  col = list(
    CMOIC = c(
      CS1 = clust.col[1],
      CS2 = clust.col[2],
      CS3 = clust.col[3],
      CS4 = clust.col[4],
      CS5 = clust.col[5]
    )
  )
)

top_anno <- anno_barplot(
  as.numeric(tmb[rownames(mutsig), "log10TMB"]),
  border = FALSE,
  gp     = gpar(fill = "#3379B4", border = NA, lty = "blank"),
  height = unit(2.5, "cm")
)

tmp <- mutsig[, c("SBS2", "SBS13", "SBS7b", "SBS7d")]
tmp$Others <- 1 - rowSums(tmp)

top_anno2 <- anno_barplot(
  as.matrix(tmp),
  border = FALSE,
  gp     = gpar(
    fill   = c(brewer.pal(6, "Paired")[c(2, 1, 6, 5)], "grey90"),
    border = NA,
    lty    = "blank"
  ),
  height = unit(2, "cm")
)

tmp <- as.data.frame(t(mut[mutgene, rownames(mutsig)]))
mut.order <- names(sort(colSums(tmp), decreasing = TRUE))
tmp$CMOIC <- subt[rownames(tmp), "CMOIC"]

dat.mut <- data.frame()

for (i in colnames(tmp[, -ncol(tmp)])) {
  aa <- tmp[, c(i, "CMOIC")]
  colnames(aa)[1] <- "var"
  aa <- kruskal.test(var ~ CMOIC, data = aa)
  aa <- data.frame(var = i, p = aa$p.value)
  dat.mut <- rbind(dat.mut, aa)
}

dat.mut$p.label <- ifelse(
  dat.mut$p < 0.0001, "****",
  ifelse(
    dat.mut$p < 0.001, "***",
    ifelse(
      dat.mut$p < 0.01, "**",
      ifelse(dat.mut$p < 0.05, "*", "")
    )
  )
)

dat.mut$var2 <- paste0(dat.mut$var, dat.mut$p.label)
rownames(onco.input) <- dat.mut$var2

pct <- NULL

for (i in mut.order) {
  tmp1 <- tmp[, c(i, "CMOIC")]
  x <- as.data.frame.array(table(tmp1[, 1], tmp1$CMOIC))
  tmp1 <- x[2, ] / apply(x, 2, sum)
  tmp1 <- tmp1 / sum(tmp1)
  pct  <- rbind.data.frame(pct, tmp1)
}

rownames(pct) <- mut.order

right_anno <- anno_barplot(
  as.matrix(pct),
  which  = "row",
  border = FALSE,
  gp     = gpar(fill = clust.col, border = NA, lty = "blank"),
  bar_width = 0.6,
  width     = unit(1.8, "cm"),
  height    = unit(1, "cm")
)

op1 <- oncoPrint(
  onco.input[mut.order, rownames(my_ann)],
  alter_fun           = alter_fun,
  col                 = col,
  bottom_annotation   = NULL,
  top_annotation      = c(
    HeatmapAnnotation(TMB = top_anno),
    my_annotation,
    HeatmapAnnotation(MutSig = top_anno2)
  ),
  column_order        = rownames(my_ann),
  right_annotation    = rowAnnotation(PCT = right_anno),
  show_pct            = TRUE,
  column_title        = "",
  show_heatmap_legend = TRUE,
  column_split        = my_ann$CMOIC,
  column_title_gp     = gpar(fontsize = 8),
  row_names_gp        = gpar(fontsize = 8),
  column_names_gp     = gpar(fontsize = 8)
)

op1

cna <- cna.region[1:(nrow(cna.region) / 2), c(1, 8, 9:(ncol(cna.region) - 1))]
rownames(cna) <- paste0(gsub(" ", "", cna$Descriptor), "-", substr(rownames(cna), 1, 3))
cna.modified <- cna[1:nrow(cna), 3:ncol(cna)]

onco.input2 <- cna.modified[, rownames(mutsig)]

tmp1 <- onco.input2[1:28, ]
tmp1[tmp1 == 1] <- "Gain"
tmp1[tmp1 == 2] <- "Gain"
tmp1[tmp1 == 0] <- ""

tmp2 <- onco.input2[29:70, ]
tmp2[tmp2 == 1] <- "Loss"
tmp2[tmp2 == 2] <- "Loss"
tmp2[tmp2 == 0] <- ""
onco.input2 <- rbind.data.frame(tmp1, tmp2)

alter_fun2 <- list(
  background = function(x, y, w, h) {
    grid.rect(
      x, y,
      w - unit(0.5, "mm"),
      h - unit(0.5, "mm"),
      gp = gpar(fill = "#dcddde", col = "#dcddde")
    )
  },
  Gain = function(x, y, w, h) {
    grid.rect(
      x, y,
      w - unit(0.5, "mm"),
      h - unit(0.5, "mm"),
      gp = gpar(fill = red, col = red)
    )
  },
  Loss = function(x, y, w, h) {
    grid.rect(
      x, y,
      w - unit(0.5, "mm"),
      h - unit(0.5, "mm"),
      gp = gpar(fill = blue, col = blue)
    )
  }
)

col2 <- c(Gain = red, Loss = blue)

lesion.order <- rownames(onco.input2)
tmp <- as.data.frame(t(cna.modified[lesion.order, rownames(mutsig)]))
tmp[tmp > 0] <- 1
tmp$CMOIC <- as.character(subt[rownames(tmp), "CMOIC"])

dat.chisq.amp <- dat.chisq.del <- data.frame()

for (i in colnames(tmp[, -ncol(tmp)])) {
  aa <- tmp[, c(i, "CMOIC")]
  colnames(aa)[1] <- "var"
  aa <- kruskal.test(var ~ CMOIC, data = aa)
  aa <- data.frame(var = i, p = aa$p.value)
  if (str_detect(i, "Amp")) {
    dat.chisq.amp <- rbind(dat.chisq.amp, aa)
  } else {
    dat.chisq.del <- rbind(dat.chisq.del, aa)
  }
}

colnames(dat.chisq.amp)[2] <- colnames(dat.chisq.del)[2] <- "p"

dat.chisq.amp$p.label <- ifelse(
  dat.chisq.amp$p < 0.0001, "****",
  ifelse(
    dat.chisq.amp$p < 0.001, "***",
    ifelse(
      dat.chisq.amp$p < 0.01, "**",
      ifelse(dat.chisq.amp$p < 0.05, "*", "")
    )
  )
)
dat.chisq.amp$var2 <- paste0(dat.chisq.amp$var, dat.chisq.amp$p.label)

dat.chisq.del$p.label <- ifelse(
  dat.chisq.del$p < 0.0001, "****",
  ifelse(
    dat.chisq.del$p < 0.001, "***",
    ifelse(
      dat.chisq.del$p < 0.01, "**",
      ifelse(dat.chisq.del$p < 0.05, "*", "")
    )
  )
)
dat.chisq.del$var2 <- paste0(dat.chisq.del$var, dat.chisq.del$p.label)

dat.chisq.amp <- dat.chisq.amp |> arrange(p)
dat.chisq.del <- dat.chisq.del |> arrange(p)

lesion.order <- c(dat.chisq.amp$var[1:5], dat.chisq.del$var[1:5])
lesion.order <- c(
  "3q26.32-Amp", "6q21-Amp", "8q24.21-Amp",
  "10p15.1-Amp", "20q13.2-Amp",
  "5q11.2-Del", "5q21.3-Del", "8p23.2-Del",
  "14q24.1-Del", "19p13.3-Del"
)

lesion.order.label <- paste0(lesion.order, "****")
tmp <- tmp[, c(lesion.order, "CMOIC")]

pct <- NULL

for (i in lesion.order) {
  tmp1 <- tmp[, c(i, "CMOIC")]
  x <- as.data.frame.array(table(tmp1[, 1], tmp1$CMOIC))
  tmp1 <- x[2, ] / apply(x, 2, sum)
  tmp1 <- tmp1 / sum(tmp1)
  pct  <- rbind.data.frame(pct, tmp1)
}

rownames(pct) <- lesion.order

right_anno2 <- anno_barplot(
  as.matrix(pct),
  which  = "row",
  border = FALSE,
  gp     = gpar(fill = clust.col, border = NA, lty = "blank"),
  bar_width = 0.6,
  width     = unit(1.8, "cm"),
  height    = unit(1, "cm")
)

op2 <- oncoPrint(
  onco.input2[lesion.order.label, rownames(my_ann)],
  alter_fun           = alter_fun2,
  col                 = col2,
  bottom_annotation   = NULL,
  top_annotation      = NULL,
  column_order        = rownames(my_ann),
  right_annotation    = rowAnnotation(PCT = right_anno2),
  row_order           = lesion.order.label,
  show_pct            = TRUE,
  column_title        = "",
  show_heatmap_legend = TRUE,
  column_split        = my_ann$CMOIC,
  column_title_gp     = gpar(fontsize = 8),
  row_names_gp        = gpar(fontsize = 8),
  column_names_gp     = gpar(fontsize = 8)
)

op2

load("humanGTF.rda")

cytoband <- read.table("cytoBand.txt", sep = "\t", header = FALSE)
colnames(cytoband) <- NULL
cytoband <- cytoband[, 1:4]

dat.humanGTF <- humanGTF[c(4:6, 1)]
colnames(dat.humanGTF) <- NULL

cytoTOsymbol <- function(cyto, cytoband, dat.humanGTF) {
  suppressMessages(library(bedtoolsr))
  suppressMessages(library(AnnoProbe))
  cyto.chr <- str_extract(cyto, "\\d")
  cyto2 <- paste0(
    str_extract(cyto, "[qp]"),
    strsplit(cyto, split = str_extract(cyto, "[qp]"))[[1]][2]
  )
  cytoband.tar <- cytoband[str_detect(cytoband[, 4], cyto2), ]
  output.bedtools <- bt.intersect(a = dat.humanGTF, b = cytoband.tar)
  output.bedtools <- output.bedtools[output.bedtools[, 1] == paste0("chr", cyto.chr), ] |>
    dplyr::arrange(V4) |>
    remove_rownames()
  output.bedtools
}

amp.genes.chr <- strsplit(dat.chisq.amp[1, 1], split = "-")[[1]][1]
id.amp <- amp.genes.chr
heat.amp <- cytoTOsymbol(id.amp, cytoband, dat.humanGTF)

cna <- cna.gene
cna <- cna[c("MYC", "PVT1", "CCDC26", "GSDMC"), rownames(mutsig)]
onco.input3 <- cna

tmp <- onco.input3 |>
  t() |>
  as.data.frame() |>
  mutate(CMOIC = mutsig$CMOIC)

dat.mut2 <- data.frame()

for (i in colnames(tmp[, -ncol(tmp)])) {
  aa <- tmp[, c(i, "CMOIC")]
  colnames(aa)[1] <- "var"
  aa <- kruskal.test(var ~ CMOIC, data = aa)
  aa <- data.frame(var = i, p = aa$p.value)
  dat.mut2 <- rbind(dat.mut2, aa)
}

colnames(dat.mut2)[2] <- "p"

dat.mut2$p.label <- ifelse(
  dat.mut2$p < 0.0001, "****",
  ifelse(
    dat.mut2$p < 0.001, "***",
    ifelse(
      dat.mut2$p < 0.01, "**",
      ifelse(dat.mut2$p < 0.05, "*", "")
    )
  )
)
dat.mut2$var2 <- paste0(dat.mut2$var, dat.mut2$p.label)
rownames(onco.input3) <- dat.mut2$var2

onco.input3[onco.input3 == 1]  <- "Gain"
onco.input3[onco.input3 == 2]  <- "High_balanced_gain"
onco.input3[onco.input3 == 0]  <- ""
onco.input3[onco.input3 == -1] <- "Loss"
onco.input3[onco.input3 == -2] <- "High_balanced_loss"

alter_fun3 <- list(
  background = function(x, y, w, h) {
    grid.rect(
      x, y,
      w - unit(0.5, "mm"),
      h - unit(0.5, "mm"),
      gp = gpar(fill = "#dcddde", col = "#dcddde")
    )
  },
  Gain = function(x, y, w, h) {
    grid.rect(
      x, y,
      w - unit(0.5, "mm"),
      h - unit(0.5, "mm"),
      gp = gpar(fill = brewer.pal(6, "Paired")[5], col = brewer.pal(6, "Paired")[5])
    )
  },
  High_balanced_gain = function(x, y, w, h) {
    grid.rect(
      x, y,
      w - unit(0.5, "mm"),
      h - unit(0.5, "mm"),
      gp = gpar(fill = brewer.pal(6, "Paired")[6], col = brewer.pal(6, "Paired")[6])
    )
  },
  Loss = function(x, y, w, h) {
    grid.rect(
      x, y,
      w - unit(0.5, "mm"),
      h - unit(0.5, "mm"),
      gp = gpar(fill = brewer.pal(6, "Paired")[1], col = brewer.pal(6, "Paired")[1])
    )
  },
  High_balanced_loss = function(x, y, w, h) {
    grid.rect(
      x, y,
      w - unit(0.5, "mm"),
      h - unit(0.5, "mm"),
      gp = gpar(fill = brewer.pal(6, "Paired")[2], col = brewer.pal(6, "Paired")[2])
    )
  }
)

col3 <- c(
  Gain               = brewer.pal(6, "Paired")[5],
  High_balanced_gain = brewer.pal(6, "Paired")[6],
  Loss               = brewer.pal(6, "Paired")[1],
  High_balanced_loss = brewer.pal(6, "Paired")[2]
)

op3 <- oncoPrint(
  onco.input3[, rownames(my_ann)],
  alter_fun           = alter_fun3,
  col                 = col3,
  bottom_annotation   = NULL,
  top_annotation      = NULL,
  column_order        = rownames(my_ann),
  right_annotation    = NULL,
  show_pct            = TRUE,
  column_title        = "",
  show_heatmap_legend = TRUE,
  column_split        = my_ann$CMOIC,
  column_title_gp     = gpar(fontsize = 8),
  row_names_gp        = gpar(fontsize = 8),
  column_names_gp     = gpar(fontsize = 8)
)

op3

del.genes.chr <- strsplit(dat.chisq.del[1, 1], split = "-")[[1]][1]
id.amp <- del.genes.chr
heat.amp <- cytoTOsymbol(id.amp, cytoband, dat.humanGTF)

cna <- cna.gene
cna <- cna[c("ITGA1", "DDX4", "GPBP1", "RAB3C"), rownames(mutsig)]
onco.input3 <- cna

tmp <- onco.input3 |>
  t() |>
  as.data.frame() |>
  mutate(CMOIC = mutsig$CMOIC)

dat.mut2 <- data.frame()

for (i in colnames(tmp[, -ncol(tmp)])) {
  aa <- tmp[, c(i, "CMOIC")]
  colnames(aa)[1] <- "var"
  aa <- kruskal.test(var ~ CMOIC, data = aa)
  aa <- data.frame(var = i, p = aa$p.value)
  dat.mut2 <- rbind(dat.mut2, aa)
}

colnames(dat.mut2)[2] <- "p"

dat.mut2$p.label <- ifelse(
  dat.mut2$p < 0.0001, "****",
  ifelse(
    dat.mut2$p < 0.001, "***",
    ifelse(
      dat.mut2$p < 0.01, "**",
      ifelse(dat.mut2$p < 0.05, "*", "")
    )
  )
)
dat.mut2$var2 <- paste0(dat.mut2$var, dat.mut2$p.label)
rownames(onco.input3) <- dat.mut2$var2

onco.input3[onco.input3 == 1]  <- "Gain"
onco.input3[onco.input3 == 2]  <- "High_balanced_gain"
onco.input3[onco.input3 == 0]  <- ""
onco.input3[onco.input3 == -1] <- "Loss"
onco.input3[onco.input3 == -2] <- "High_balanced_loss"

alter_fun3 <- list(
  background = function(x, y, w, h) {
    grid.rect(
      x, y,
      w - unit(0.5, "mm"),
      h - unit(0.5, "mm"),
      gp = gpar(fill = "#dcddde", col = "#dcddde")
    )
  },
  Gain = function(x, y, w, h) {
    grid.rect(
      x, y,
      w - unit(0.5, "mm"),
      h - unit(0.5, "mm"),
      gp = gpar(fill = brewer.pal(6, "Paired")[5], col = brewer.pal(6, "Paired")[5])
    )
  },
  High_balanced_gain = function(x, y, w, h) {
    grid.rect(
      x, y,
      w - unit(0.5, "mm"),
      h - unit(0.5, "mm"),
      gp = gpar(fill = brewer.pal(6, "Paired")[6], col = brewer.pal(6, "Paired")[6])
    )
  },
  Loss = function(x, y, w, h) {
    grid.rect(
      x, y,
      w - unit(0.5, "mm"),
      h - unit(0.5, "mm"),
      gp = gpar(fill = brewer.pal(6, "Paired")[1], col = brewer.pal(6, "Paired")[1])
    )
  },
  High_balanced_loss = function(x, y, w, h) {
    grid.rect(
      x, y,
      w - unit(0.5, "mm"),
      h - unit(0.5, "mm"),
      gp = gpar(fill = brewer.pal(6, "Paired")[2], col = brewer.pal(6, "Paired")[2])
    )
  }
)

col3 <- c(
  Gain               = brewer.pal(6, "Paired")[5],
  High_balanced_gain = brewer.pal(6, "Paired")[6],
  Loss               = brewer.pal(6, "Paired")[1],
  High_balanced_loss = brewer.pal(6, "Paired")[2]
)

op4 <- oncoPrint(
  onco.input3[, rownames(my_ann)],
  alter_fun           = alter_fun3,
  col                 = col3,
  bottom_annotation   = NULL,
  top_annotation      = NULL,
  column_order        = rownames(my_ann),
  right_annotation    = NULL,
  show_pct            = TRUE,
  column_title        = "",
  show_heatmap_legend = TRUE,
  column_split        = my_ann$CMOIC,
  column_title_gp     = gpar(fontsize = 8),
  row_names_gp        = gpar(fontsize = 8),
  column_names_gp     = gpar(fontsize = 8)
)

op4

lgd.mutsig <- Legend(
  labels    = c("SBS2", "SBS13", "SBS7b", "SBS7d", "others"),
  title     = "MutSig",
  legend_gp = gpar(fill = c(brewer.pal(6, "Paired")[c(2, 1, 6, 5)], "grey90"))
)

lgd.cna.region <- Legend(
  labels    = c("Gain", "Loss"),
  title     = "CNA (arm-level)",
  legend_gp = gpar(fill = c(red, blue))
)

lgd.cna.gene <- Legend(
  labels    = c("Gain", "High_balanced_gain", "Loss", "High_balanced_loss"),
  title     = "CNA (gene-level)",
  legend_gp = gpar(fill = brewer.pal(6, "Paired")[c(5, 6, 1, 2)])
)

lgd.cna.gene.amp <- Legend(labels = amp.genes.chr)
lgd.cna.gene.del <- Legend(labels = del.genes.chr)

lgd_list <- list(
  lgd.mutsig,
  lgd.cna.region,
  lgd.cna.gene,
  lgd.cna.gene.amp,
  lgd.cna.gene.del
)

draw(op1 %v% op2 %v% op3 %v% op4, annotation_legend_list = lgd_list)

ylim1 <- boxplot.stats(com.tmb$TMB.dat$log10TMB)$stats[c(1, 5)]

com.tmb$TMB.dat |>
  ggplot(aes(Subtype, log10TMB, fill = Subtype)) +
  geom_boxplot(notch = TRUE, alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = Subtype), alpha = 0.5) +
  coord_cartesian(ylim = ylim1 * 1.15) +
  scale_fill_manual(values = clust.col) +
  scale_color_manual(values = clust.col) +
  stat_compare_means(size = 10, label = "p.format", label.y = 0.8) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title      = element_text(size = 24, colour = "red"),
    axis.title.x    = element_blank(),
    axis.text.x     = element_text(size = 20),
    axis.text.y     = element_text(size = 20),
    title           = element_text(size = 24)
  )

ylim1 <- boxplot.stats(mutsig$Signature.1)$stats[c(1, 5)]

mutsig2 <- mutsig |>
  dplyr::select(CMOIC, everything()) |>
  melt()

mutsig2$value[mutsig2$value == 0] <- NA
mutsig2 <- na.omit(mutsig2)

mutsig2 |>
  ggplot(aes(CMOIC, value, fill = CMOIC)) +
  geom_boxplot(notch = FALSE, alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = CMOIC), alpha = 0.5) +
  scale_fill_manual(values = clust.col) +
  scale_color_manual(values = clust.col) +
  stat_compare_means(size = 10, label = "p.format", label.y = 0.6) +
  facet_wrap(~variable) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title      = element_text(size = 24, colour = "red"),
    axis.title.x    = element_blank(),
    axis.text.x     = element_text(size = 20),
    axis.text.y     = element_text(size = 20),
    title           = element_text(size = 24)
  )

subset(mutsig2, variable == "APOBEC") |>
  ggplot(aes(CMOIC, value, fill = CMOIC)) +
  geom_boxplot(notch = FALSE, alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = CMOIC), alpha = 0.5) +
  scale_fill_manual(values = clust.col) +
  scale_color_manual(values = clust.col) +
  ylab("APOBEC") +
  stat_compare_means(size = 10, label = "p.format", label.y = 0.97) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title      = element_text(size = 24, colour = "red"),
    axis.title.x    = element_blank(),
    axis.text.x     = element_text(size = 20),
    axis.text.y     = element_text(size = 20),
    title           = element_text(size = 24)
  )

mutsig3 <- subset(mutsig2, variable == "APOBEC" | variable == "SBS2" | variable == "SBS13")
mutsig3$variable <- factor(mutsig3$variable, levels = c("APOBEC", "SBS2", "SBS13"))

mutsig3 |>
  ggplot(aes(CMOIC, value, fill = CMOIC)) +
  geom_boxplot(notch = FALSE, alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = CMOIC), alpha = 0.5) +
  scale_fill_manual(values = clust.col) +
  scale_color_manual(values = clust.col) +
  ylab("Mutational Signature Load") +
  stat_compare_means(size = 10, label = "p.format") +
  facet_wrap(~variable, scales = "free_y") +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title      = element_text(size = 24, colour = "red"),
    strip.text      = element_text(size = 28),
    axis.title.x    = element_blank(),
    axis.text.x     = element_text(size = 20),
    axis.text.y     = element_text(size = 20),
    title           = element_text(size = 24)
  )

com.mut2 <- com.mut
colnames(com.mut2)[1] <- "gene"
com.mut2$pvalue <- as.numeric(com.mut2$pvalue)

com.mut2 <- com.mut2 |>
  dplyr::filter(pvalue < 0.05) |>
  dplyr::select(-TMB, -padj, -pvalue) |>
  pivot_longer(., 2:6, names_to = "CS", values_to = "mut") |>
  mutate(
    percentage = str_extract(mut, "\\(.*?\\)"),
    percentage = str_replace_all(percentage, "[\\(\\)\\s]", ""),
    percentage = str_replace(percentage, "%", ""),
    percentage = as.numeric(percentage)
  )

gene_order <- c(
  "TP53", "PIK3CA", "TTN", "CDH1", "GATA3",
  "SPTA1", "MAP3K1", "SYNE1", "DST", "MUC17", "ZFHX4"
)

com.mut2 <- com.mut2 |>
  mutate(
    gene = factor(gene, levels = gene_order),
    CS   = factor(CS, levels = c("CS1", "CS2", "CS3", "CS4", "CS5"))
  )

color_schemes <- list(
  CS1 = ggsci::pal_lancet()(9)[1],
  CS2 = ggsci::pal_lancet()(9)[2],
  CS3 = ggsci::pal_lancet()(9)[3],
  CS4 = ggsci::pal_lancet()(9)[4],
  CS5 = ggsci::pal_lancet()(9)[5]
)

mut.output <- list()

for (i in 1:5) {
  input <- com.mut2[com.mut2$CS == paste0("CS", i), ]
  p <- ggplot(input, aes(gene, CS)) +
    geom_tile(aes(fill = percentage), color = "grey", size = 0.8, width = 1, height = 1) +
    geom_text(aes(label = paste0(percentage, "%")), size = 3, color = "black") +
    scale_fill_gradient(low = "white", high = color_schemes[[paste0("CS", i)]]) +
    theme_minimal() +
    theme(
      axis.title.x  = element_blank(),
      axis.ticks.x  = element_blank(),
      axis.title.y  = element_blank(),
      axis.text.x   = element_blank(),
      axis.text.y   = element_text(size = 8),
      plot.margin   = unit(c(0.01, 0.01, 0.01, 0.01), "cm"),
      legend.position = "none",
      panel.spacing = unit(0.01, "lines")
    )
  if (i == 5) {
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  }
  mut.output[[i]] <- p
}

pp <- plot_grid(
  plotlist    = mut.output,
  ncol        = 1,
  align       = "v",
  rel_heights = c(rep(1, 4), 1.5)
)

print(pp)
