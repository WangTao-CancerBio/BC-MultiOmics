library(PharmacoGx)
library(parallel)
library(ggrepel)
library(ggthemes)
library(gridExtra)
library(cowplot)
library(ISOpureR) 
library(impute) 
library(pRRophetic) 
library(SimDesign) 
library(cowplot) 
display.progress = function (index, totalN, breakN=20) {
        if ( index %% ceiling(totalN/breakN)  ==0  ) {
                cat(paste(round(index*100/totalN), "% ", sep=""))
        }
} 
(load("../../TCGA project/dataTCGAbiolinks/TCGA-BRCA-TPM_Symbol.rda"))
expr <- log2(tpm+1)
normsam <- colnames(expr[,which(substr(colnames(expr),14,15) == "11")])
tumosam <- colnames(expr[,which(substr(colnames(expr),14,15) == "01")])
normexpr <- as.matrix(expr[,normsam])
tumoexpr <- as.matrix(expr[,tumosam])
runpure <- F
if(runpure) {
        set.seed(123)
        ISOpureS1model <- ISOpure.step1.CPE(tumoexpr, normexpr)
        set.seed(456);
        ISOpureS2model <- ISOpure.step2.PPE(tumoexpr,normexpr,ISOpureS1model)
        pure.tumoexpr <- ISOpureS2model$cc_cancerprofiles
}
if(!runpure) {
        pure.tumoexpr <- tumoexpr
}
pure.tumoexpr <- eset
(load("../_CA/1_BRCA_prognosis_YB/r0.data/drug/trainDATA.rda"))
dim(trainExpr);dim(trainPtype)
keepgene <- apply(pure.tumoexpr, 1, mad) > 0.5
testExpr <- log2(pure.tumoexpr[keepgene,] + 1)
comgene <- intersect(rownames(trainExpr),rownames(testExpr)) 
trainExpr <- as.matrix(trainExpr[comgene,])
testExpr <- as.matrix(testExpr[comgene,])
dim(trainExpr);dim(testExpr)
outTab <- NULL
for (i in 1:ncol(trainPtype)) { 
        display.progress(index = i,totalN = ncol(trainPtype))
        d <- colnames(trainPtype)[i]
        tmp <- log2(as.vector(trainPtype[,d]) + 0.00001)
        ptypeOut <- quiet(calcPhenotype(trainingExprData = trainExpr,
                                        trainingPtype = tmp,
                                        testExprData = testExpr,
                                        powerTransformPhenotype = F,
                                        selection = 1))
        ptypeOut <- 2^ptypeOut - 0.00001
        outTab <- rbind.data.frame(outTab,ptypeOut)
}
dimnames(outTab) <- list(colnames(trainPtype),colnames(testExpr))
ctrp.pred.auc <- outTab
keepgene <- apply(ccl.expr, 1, mad) > 0.5
trainExpr <- log2(ccl.expr[keepgene,] + 1)
colnames(trainExpr) <- sapply(strsplit(colnames(trainExpr),"_",fixed = True),"[",1)
trainPtype <- as.data.frame(prism.auc.knn)
rownames(trainPtype) <- prism.ccl.anno[rownames(trainPtype),"cell_line_display_name"]
comccl <- intersect(rownames(trainPtype),colnames(trainExpr))
trainExpr <- trainExpr[,comccl]
trainPtype <- trainPtype[comccl,]
keepgene <- apply(pure.tumoexpr, 1, mad) > 0.5
testExpr <- log2(pure.tumoexpr[keepgene,] + 1)
comgene <- intersect(rownames(trainExpr),rownames(testExpr))
trainExpr <- as.matrix(trainExpr[comgene,])
testExpr <- as.matrix(testExpr[comgene,])
outTab <- None

drug.cut <- 0.05
candidate.ctrp <- candidate.prism <- list()
for(j in 1:5){
    top.pps <- pdata_cmoic$ID[pdata_cmoic$CMOIC == paste0("CS", j)]
    bot.pps <- pdata_cmoic$ID[pdata_cmoic$CMOIC != paste0("CS", j)]
    ctrp.log2fc <- c()
    for (i in 1:nrow(ctrp.pred.auc)) {
        display.progress(index = i, totalN = nrow(ctrp.pred.auc))
        d <- rownames(ctrp.pred.auc)[i]
        a <- mean(as.numeric(ctrp.pred.auc[d, top.pps]))
        b <- mean(as.numeric(ctrp.pred.auc[d, bot.pps]))
        fc <- b / a
        log2fc <- log2(fc)
        names(log2fc) <- d
        ctrp.log2fc <- c(ctrp.log2fc, log2fc)
    }
    candidate.ctrp[[j]] <- ctrp.log2fc[ctrp.log2fc > drug.cut]

    prism.log2fc <- c()
    for (i in 1:nrow(prism.pred.auc)) {
        display.progress(index = i, totalN = nrow(prism.pred.auc))
        d <- rownames(prism.pred.auc)[i]
        a <- mean(as.numeric(prism.pred.auc[d, top.pps]))
        b <- mean(as.numeric(prism.pred.auc[d, bot.pps]))
        fc <- b / a
        log2fc <- log2(fc)
        names(log2fc) <- d
        prism.log2fc <- c(prism.log2fc, log2fc)
    }
    candidate.prism[[j]] <- prism.log2fc[prism.log2fc > drug.cut]
}

candidate.prism <- lapply(candidate.prism, function(x) {
    names(x) <- gsub("\\s*\\([^\\)]+\\)", "", names(x))
    return(x)
})

top_n <- 3
top_candidates.ctrp <- lapply(candidate.ctrp, function(x) {
    x <- sort(x, decreasing = TRUE)
    head(x, top_n)
})
drug.cut <- 0.05
candidate.ctrp <- candidate.prism <- list()
for(j in 1:5){
        top.pps <- pdata_cmoic$ID[pdata_cmoic$CMOIC==paste0("CS",j)]
        bot.pps <- pdata_cmoic$ID[pdata_cmoic$CMOIC!=paste0("CS",j)]
        ctrp.log2fc <- c()
        for (i in 1:nrow(ctrp.pred.auc)) {
                display.progress(index = i,totalN = nrow(ctrp.pred.auc))
                d <- rownames(ctrp.pred.auc)[i]
                a <- mean(as.numeric(ctrp.pred.auc[d,top.pps]))
                b <- mean(as.numeric(ctrp.pred.auc[d,bot.pps]))
                fc <- b/a
                log2fc <- log2(fc); names(log2fc) <- d
                ctrp.log2fc <- c(ctrp.log2fc,log2fc)
        }
        candidate.ctrp[[j]] <- ctrp.log2fc[ctrp.log2fc > drug.cut]

        prism.log2fc <- c()
        for (i in 1:nrow(prism.pred.auc)) {
                display.progress(index = i,totalN = nrow(prism.pred.auc))
                d <- rownames(prism.pred.auc)[i]
                a <- mean(as.numeric(prism.pred.auc[d,top.pps]))
                b <- mean(as.numeric(prism.pred.auc[d,bot.pps]))
                fc <- b/a
                log2fc <- log2(fc); names(log2fc) <- d
                prism.log2fc <- c(prism.log2fc,log2fc)
        }
        candidate.prism[[j]] <- prism.log2fc[prism.log2fc > drug.cut]
}
candidate.prism <- lapply(candidate.prism, function(x) {
        names(x) <- gsub("\\s*\\([^\\)]+\\)", "", names(x))
        return(x)
})
top_n <- 3
top_candidates.ctrp <- lapply(candidate.ctrp, function(x) {
        x <- sort(x, decreasing = TRUE)
        head(x, top_n)
})
top_candidates.prism <- lapply(candidate.prism, function(x) {
        x <- sort(x, decreasing = TRUE)
        head(x, top_n)
})
ctrp.candidate <- unlist(top_candidates.ctrp)
prism.candidate <- unlist(top_candidates.prism)
ctrp.candidate.C2 <- top_candidates.ctrp[[2]] %>% names(.)
prism.candidate.C2 <- top_candidates.prism[[2]] %>% names(.)
