library(IOBR)

crRegulators <- read.csv("chromatin remodelling.csv")
load("DNArepairGenesets.rda")
DNArepairGenesets$`Chromatin remodelling` <- crRegulators$regulators

eset <- dat.TPM
pdata_cmoic <- mutsig %>%
        rownames_to_column(var = "ID") %>%
        dplyr::select(ID, CMOIC)

sig_tme_DNA <- calculate_sig_score(pdata           = NULL,
                                   eset            = eset,
                                   signature       = DNArepairGenesets,
                                   method          = "ssGSEA",
                                   mini_gene_count = 2)

sig_tme_DNA <- t(column_to_rownames(sig_tme_DNA, var = "ID"))
sig_tme_DNA[1:5, 1:3]

input_DNA <- combine_pd_eset(eset = sig_tme_DNA, pdata = pdata_cmoic, scale = TRUE)

p_DNA <- sig_heatmap(input         = input_DNA, 
                     features      = rownames(sig_tme_DNA)[-c(4, 9)],
                     group         = "CMOIC", 
                     palette_group = "jama", 
                     palette       = 6,
                     path          = "result")

data_summary <- input_DNA[, -c(6, 11)] %>%
        pivot_longer(., 3:10, names_to = "Pathway", values_to = "Score") %>%
        group_by(Pathway, CMOIC) %>%
        summarise(mean_score = mean(Score, na.rm = TRUE)) %>%
        ungroup()

ggplot(data_summary, aes(x = CMOIC, y = Pathway, fill = mean_score)) +
        geom_tile(color = "white", size = 0.8) +
        geom_text(aes(label = round(mean_score, 2)), size = 4, color = "black") +
        scale_fill_gradientn(colors = c("#3288BD", "white", "#D53E4F"), 
                             limits = c(min(data_summary$mean_score), max(data_summary$mean_score)), 
                             breaks = c(min(data_summary$mean_score), 0, max(data_summary$mean_score)),
                             labels = c("Lower score", "0", "Higher score"),
                             name = "Mean score") +
        theme_minimal() +
        theme(axis.text.x = element_text(size = 12, colour = ggsci::pal_lancet()(9)[1:5]),
              axis.text.y = element_text(size = 12),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              plot.title = element_text(size = 14, face = "bold"),
              legend.position = "bottom",
              legend.title = element_text(size = 10),
              legend.text = element_text(size = 8)) +
        labs(title = "DNA Repair Signaling", fill = "Mean score") +
        guides(fill = guide_colorbar(barwidth = 10, barheight = 1, title.position = "top", title.hjust = 0.5))

load("GS.score.rda")
GS.score <- GS.score[, -1]
GS.score <- GS.score %>%
        rownames_to_column(var = "ID") %>%
        merge(pdata_cmoic, ., by = "ID")

selected_columns <- GS.score[, -1] %>%
        select(CMOIC, Intratumor.Heterogeneity, SNV.Neoantigens, 
               Nonsilent.Mutation.Rate, Aneuploidy.Score, 
               Homologous.Recombination.Defects, LOH_n_seg, 
               LOH_frac_altered)

normalize <- function(x) (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

normalized_columns <- selected_columns %>%
        mutate(across(-CMOIC, normalize))

group_means <- normalized_columns %>%
        group_by(CMOIC) %>%
        summarise(across(everything(), mean, na.rm = TRUE))

melted_means <- melt(group_means, id.vars = "CMOIC")
melted_means$variable <- gsub("_", " ", melted_means$variable)
melted_means$variable <- gsub("\\.", " ", melted_means$variable)
melted_means$variable <- factor(melted_means$variable, levels = sort(unique(melted_means$variable), decreasing = TRUE))

ggplot(melted_means, aes(x = CMOIC, y = variable, fill = value)) +
        geom_tile(color = "white") +
        scale_fill_gradient(low = "white", high = "#D53E4F", name = "Normalized Mean Value") +
        geom_text(aes(label = round(value, 2)), color = "black", size = 4) +
        theme_minimal() +
        labs(title = "DNA Damage Measures") +
        theme_minimal() +
        theme(axis.text.x = element_text(size = 12, colour = ggsci::pal_lancet()(9)[1:5]),
              axis.text.y = element_text(size = 12),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              plot.title = element_text(size = 14, face = "bold"),
              legend.position = "bottom",
              legend.title = element_text(size = 10),
              legend.text = element_text(size = 8)) +
        guides(fill = guide_colorbar(barwidth = 10, barheight = 1, title.position = "top", title.hjust = 0.5))


breastGenesets <- read.csv("breastSigs.csv", header = T)
breastGenesets <- split(breastGenesets, breastGenesets[,1])
length(breastGenesets)
breastGenesets <- lapply(breastGenesets, function(x) x[,-1])

(load("breastSigs.rda"))

sig_tme_BC <- calculate_sig_score(pdata           = NULL,
                                  eset            = eset,
                                  signature       = breastGenesets,
                                  method          = "ssGSEA",
                                  mini_gene_count = 2)

sig_tme_BC <- t(column_to_rownames(sig_tme_BC, var = "ID"))
sig_tme_BC[1:5, 1:3]

input_BC <- combine_pd_eset(eset = sig_tme_BC, pdata = pdata_cmoic, scale = T)

p_BC <- sig_heatmap(input         = input_BC, 
                    features      = rownames(sig_tme_BC),
                    group         = "CMOIC", 
                    palette_group = "jama", 
                    palette       = 6,
                    path          = "result" )

data_summary <- input_BC %>%
        pivot_longer(.,3:ncol(input_BC), names_to = "Pathway", values_to = "Score") %>%
        group_by(Pathway, CMOIC) %>%
        summarise(mean_score = mean(Score, na.rm = TRUE)) %>%
        ungroup()
data_summary$Pathway <- factor(data_summary$Pathway, levels = sort(unique(data_summary$Pathway),decreasing = T)

ggplot(data_summary, aes(x = CMOIC, y = Pathway, fill = mean_score)) +
        geom_tile(color = "white", size = 0.8) +
        geom_text(aes(label = round(mean_score, 2)), size = 4, color = "black") +
        scale_fill_gradientn(colors = c("#3288BD", "white", "#D53E4F")) +
        theme_minimal() +
        labs(title = "Breast cancer signatures")

(load("OncoPathGenesets.rda"))

sig_tme_onco <- calculate_sig_score(pdata=NULL,eset=eset,signature=OncoPathGenesets,method="ssGSEA",mini_gene_count=2)
sig_tme_onco <- t(column_to_rownames(sig_tme_onco, var = "ID"))
input_onco <- combine_pd_eset(eset = sig_tme_onco, pdata = pdata_cmoic, scale = T)

data_summary <- input_onco %>%
        pivot_longer(.,3:ncol(input_onco), names_to = "Pathway", values_to = "Score") %>%
        group_by(Pathway, CMOIC) %>%
        summarise(mean_score = mean(Score, na.rm = TRUE)) %>%
        ungroup()
data_summary$Pathway <- factor(data_summary$Pathway, levels = sort(unique(data_summary$Pathway),decreasing = T)

ggplot(data_summary, aes(x = CMOIC, y = Pathway, fill = mean_score)) +
        geom_tile(color = "white", size = 0.8) +
        geom_text(aes(label = round(mean_score, 2)), size = 4, color = "black") +
        scale_fill_gradientn(colors = c("#3288BD", "white", "#D53E4F")) +
        theme_minimal() +
        labs(title = "Oncogenic signaling")

dat.regulon <- read.table("BRCA_TCGA-regulons.txt",sep="\t",header=T)
TF.regulon <- unique(dat.regulon$TF)

input <- eset[TF.regulon,] %>% na.omit() %>% t() %>% as.data.frame() %>%
        rownames_to_column("ID") %>% merge(.,pdata_cmoic,by="ID") %>%
        dplyr::select(ID,CMOIC,everything()) %>% dplyr::arrange(CMOIC) %>%
        column_to_rownames("ID")

dat.p <- data.frame()
for(i in colnames(input[,-1])){
        aa <- input[,c(i,"CMOIC")]
        colnames(aa)[1] <- "var"
        aa <- kruskal.test(var~CMOIC,data=aa)
        aa <- data.frame(var=i,p=aa$p.value)
        dat.p <- rbind(dat.p,aa)
}
dat.p <- dat.p %>% arrange(p)

input2 <- input[,dat.p$var[1:30]]
input2 <- t(scale(input2))

col_fun <- colorRamp2(c(-1,0,1),c("#3288BD","white","#D53E4F"))

Heatmap(input2, col=col_fun, cluster_rows=T, cluster_columns=F)

sig_tme_metabolism <- calculate_sig_score(pdata=NULL,eset=eset,signature=signature_metabolism,method="ssGSEA",mini_gene_count=2)
sig_tme_metabolism <- t(column_to_rownames(sig_tme_metabolism,"ID"))
input_metabolism <- combine_pd_eset(eset=sig_tme_metabolism,pdata=pdata_cmoic,scale=T)

dat.metabolism.p <- data.frame()
for(i in colnames(input_metabolism[,-1:-2])){
        aa <- input_metabolism[,c(i,"CMOIC")]
        colnames(aa)[1] <- "var"
        aa <- kruskal.test(var~CMOIC,data=aa)
        aa <- data.frame(var=i,p=aa$p.value)
        dat.metabolism.p <- rbind(dat.metabolism.p,aa)
}
dat.metabolism.p <- dat.metabolism.p %>% arrange(p)

input_metabolism2 <- input_metabolism[,dat.metabolism.p$var[1:30]]
input_metabolism2 <- t(scale(input_metabolism2))

Heatmap(input_metabolism2, col=col_fun, cluster_rows=T, cluster_columns=F)
