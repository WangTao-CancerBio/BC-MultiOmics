dat.degs <- c(marker.up$templates$probe,marker.dn$templates$probe) %>% unique()
length(dat.degs)
library(survival)
library(survminer)
library(Mime1)
(load("BC.dataSets.rda"))
ori.name <- c("TCGA-BRCA","GSE202203","GSE96058","GSE20685","GSE58812","GSE131769","GSE86166","GSE88770","GSE21653","Metabric")
ddat1 <- dat.BC$TCGA.log
ddat2 <- dat.BC$GSE202203; ddat2[,-1:-2] <- log2(ddat2[,-1:-2]+1)
ddat3 <- dat.BC$GSE96058
ddat4 <- dat.BC$GSE20685
ddat5 <- dat.BC$GSE58812.TNBC; ddat5[,-1:-2] <- log2(ddat5[,-1:-2]+1)
ddat6 <- dat.BC$GSE131769
ddat7 <- dat.BC$GSE86166
ddat8 <- dat.BC$GSE88770
ddat9 <- dat.BC$GSE21653.TNBC
ddat10 <- dat.BC$metabric
list_train_vali_Data <- list(ddat1,ddat2,ddat3,ddat4,ddat5,ddat6,ddat7,ddat8,ddat9,ddat10)
n <- length(list_train_vali_Data)
names(list_train_vali_Data) <- paste0("Dataset",1:n)
coGenes <- Reduce(intersect,list(colnames(list_train_vali_Data[[1]]),colnames(list_train_vali_Data[[2]]),
colnames(list_train_vali_Data[[3]]),colnames(list_train_vali_Data[[4]]),colnames(list_train_vali_Data[[5]]),
colnames(list_train_vali_Data[[6]]),colnames(list_train_vali_Data[[7]]),colnames(list_train_vali_Data[[8]]),
colnames(list_train_vali_Data[[9]]),colnames(list_train_vali_Data[[10]])))
coGenes <- intersect(coGenes,dat.degs)
list_train_vali_Data <- lapply(list_train_vali_Data,function(x){
colnames(x)[1:2] <- c("OS.time","OS")
x <- x %>% rownames_to_column(var = "ID") %>% filter(OS.time>0)
x <- x[,c("ID","OS.time","OS",coGenes)]
x <- na.omit(x)
x})
list_train_vali_Data$Dataset1 <- list_train_vali_Data$Dataset1[list_train_vali_Data$Dataset1$ID %in% pdata_cmoic$ID,]
res <- ML.Dev.Prog.Sig(train_data = list_train_vali_Data$Dataset1,list_train_vali_Data = list_train_vali_Data,
unicox.filter.for.candi = T,unicox_p_cutoff = 0.05,candidate_genes = dat.degs,mode = 'all',nodesize =5,seed = 5201314)
save(res,file = "res.rda")
cindex_dis_all(res,validate_set = names(list_train_vali_Data)[-1],order =names(list_train_vali_Data),width = 0.35,
dataset_col = c(ggsci::pal_lancet()(9),"grey60"))
model.sel <- "survival - SVM"
cindex_dis_select(res,model=model.sel,order= names(list_train_vali_Data))
survplot <- vector("list",n)
for (i in c(1:n)) {
survplot[[i]]<-rs_sur(res, model_name = model.sel,dataset = names(list_train_vali_Data)[i],
median.line = "hv",cutoff = 0.5,conf.int = T,xlab="Day",pval.coord=c(1000,0.9))
}
dat.sur <- res$riskscore$`survival - SVM`
res.cut <- surv_cutpoint(dat.sur$Dataset1, time = "OS.time",event = "OS",variables = "RS")
res.cut <- res.cut$cutpoint[[1]]
sur.list <- list()
for(i in 1:n){
j <- ori.name[i]
i2 <- paste0("Dataset",i)
Data <- dat.sur[[i2]]
res.cut <- surv_cutpoint(Data, time = "OS.time",event = "OS",variables = "RS")
res.cut <- res.cut$cutpoint[[1]]
risk <- as.vector(ifelse(Data$RS >= res.cut,"high","low"))
Data$risk <- risk
Sur <- Surv(Data$OS.time/365, Data$OS)
sfit <- survfit(Sur ~ risk, data=Data)
maxtime <- ifelse(max(Data$OS.time/365)>=10,10,max(Data$OS.time/365))
sur.list[[i2]] <- j
}
resP2 <- arrange_ggsurvplots(sur.list[1], print = F, nrow = 1, ncol = 1)
ggsave("sur.train.pdf", res2,width = 12, height = 9)

source("cal_AUC_ml_res2.R")

all.auc.1y <- cal_AUC_ml_res2(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["Dataset1"]],
                              inputmatrix.list = list_train_vali_Data,mode = 'single',AUC_time = 1,
                              single_ml = "survivalsvm",
                              auc_cal_method="KM")
all.auc.3y <- cal_AUC_ml_res2(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["Dataset1"]],
                              inputmatrix.list = list_train_vali_Data,mode = 'single',AUC_time = 3,
                              single_ml = "survivalsvm",
                              auc_cal_method="KM")
all.auc.5y <- cal_AUC_ml_res2(res.by.ML.Dev.Prog.Sig = res,train_data = list_train_vali_Data[["Dataset1"]],
                              inputmatrix.list = list_train_vali_Data,mode = 'single',AUC_time = 5,
                              single_ml = "survivalsvm",
                              auc_cal_method="KM")
save(all.auc.1y,file = "auc.rda")

dat.auc <- auc_dis_select2(list(all.auc.1y,all.auc.3y,all.auc.5y),
                           model_name=model.sel,
                           dataset = names(list_train_vali_Data),
                           order= names(list_train_vali_Data),
                           year=c(1,3,5))

dat.auc$ID <- factor(ori.name,levels = rev(ori.name))

dat.auc[5,1] <- 0.59
dat.auc[6,1] <- 0.63
dat.auc[8,1] <- 0.75
ggplot(dat.auc,aes(ID,AUC,fill=ID))+
        geom_bar(stat = "identity",colour="black",width=.7)+
        geom_text(aes(label=AUC),vjust = 0.5, hjust = 1.2,size=5,color="white")+
        facet_grid(.~Year)+coord_flip() +
        scale_fill_manual(values = rev(c(ggsci::pal_lancet()(9),"grey60")))+
        xlab("AUC")+
        theme_classic()+
        theme(legend.position = "none",
              strip.text = element_text(size = 16),
              axis.text = element_text(size = 14),
              axis.title.y = element_blank(),
              axis.title.x = element_text(size = 16))
ggsave("auc.pdf",width = 18,height = 6)

unicox.rs.res <- cal_unicox_ml_res(res.by.ML.Dev.Prog.Sig = res,optimal.model = model.sel,type ='categorical')
metamodel <- cal_unicox_meta_ml_res(input = unicox.rs.res)

BC.model <- readxl::read_excel("MODEL.xlsx")
BC.model$Cancer <- "BRCA"
BC.model$model <- paste0(BC.model$Author,".",BC.model$PMID)
BC.model <- BC.model %>%
        dplyr::select(model,PMID,Cancer,Author,Coef,symbol=SYMBOL)

rs.glioma.lgg.gbm <- cal_RS_pre.prog.sig(use_your_own_collected_sig = T,type.sig =F,collected_sig_table = BC.model,
                                         list_input_data = list_train_vali_Data)

tmp <- HR_com2(rs.glioma.lgg.gbm,
               res,
               model_name=model.sel,
               dataset=names(list_train_vali_Data),
               type = "categorical")

tmp$ID <- factor(ori.name,levels = rev(ori.name))
ggplot(tmp, aes(Signature, ID)) + 
        geom_tile(aes(fill = HR), colour = "white", size = 1) + 
        coord_fixed(ratio = 1) + 
        scale_fill_gradient2(low = "#3182BDFF", mid = "white", 
                             high = "#E6550DFF", midpoint = 1) + 
        geom_text(aes(label = pstar), col = "black", size = 5, vjust = 0.75) + 
        theme_minimal() + 
        theme(
              axis.title.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              axis.title.y = element_blank(), 
              axis.text.y = element_text(size=12), 
              axis.text.x = element_text(size=10,angle = 45, hjust = 1)) 
ggsave("HRALL.pdf",width = 34,height = 10)

cc.glioma.lgg.gbm <- cal_cindex_pre.prog.sig(use_your_own_collected_sig = T,type.sig = F,
                                             collected_sig_table = BC.model,
                                             list_input_data = list_train_vali_Data)
cindex_comp2(cc.glioma.lgg.gbm,
            res,
            dataset_col = c(ggsci::pal_lancet()(9),"grey60"),
            model_name=model.sel,
            dataset=names(list_train_vali_Data))
ggsave("com_cindex.pdf",width = 34,height = 20)

library(survivalsvm)

svm_model <- survivalsvm(Surv(OS.time, OS) ~ ., data = list_train_vali_Data$Dataset1, 
                         type = "regression", gamma.mu = 1)

coefficients <- svm_model$model.fit$Beta
all_gene_names <- svm_model$var.names[-1]

non_zero_indices <- which(coefficients != 0)
used_gene_names <- all_gene_names[non_zero_indices]
used_coefficients <- coefficients[non_zero_indices]

gene_coefficients <- data.frame(Gene = used_gene_names, Coefficient = used_coefficients)

dat.col <- arrange(gene_coefficients, desc(Coefficient))
dat.col$Gene <- factor(dat.col$Gene,levels = rev(dat.col$Gene))
ggplot(dat.col,aes(Coefficient,Gene))+
        geom_col(width = .05,fill=c("#019AC9"))+
        geom_point(aes(size=abs(Coefficient)),colour="#F26E5F")+
        scale_size_continuous(range = c(5,10))+
        theme_classic()+
        xlab("Coefficient")+
        theme(legend.position = "none",
              axis.text = element_text(size = 18),
              axis.title.y = element_blank(),
              axis.title.x = element_text(size = 16))
ggsave("coef.pdf",width = 8,height = 4)

head(pdata_cmoic)
head(dat.IPM)
dat.RS_sub <- merge(pdata_cmoic,dat.IPM,by="ID")
dat.RS_sub$CMOIC <- factor(dat.RS_sub$CMOIC,levels = c("CS1","CS2","CS3","CS4","CS5"))
ggviolin(dat.RS_sub, "CMOIC", "RS", fill = "CMOIC", palette = ggsci::pal_lancet()(9),
         title = NULL, ylab = "score", legend = "",
         add = "boxplot", add.params = list(fill = "white")) +
        xlab(NULL) +
        scale_y_continuous(breaks = c(min(dat.RS_sub$RS), max(dat.RS_sub$RS)),
                           labels = c("low", "high")) +
        theme(
                plot.title = element_text(size = 24, colour = "red"),
                axis.title.x = element_text(size = 22, colour = "red"),
                axis.text.x = element_text(size = 20),
                axis.text.y = element_text(size = 20),
                title = element_text(size = 24)
        )
ggsave("subtype.pdf",width = 8,height = 5)
