####This is a code to reproduce statistical analysis of proteomics data from Lobov et al., 2022 "Crenigacestat (LY3039478) inhibits osteogenic differentiation of human valve interstitial cells and seems promising for treatment of cardiovascular calcification"
#Openinig_the_data
#protein_expression_data
#setwd("")

dat <- data.frame(read.csv("proteins_cre.csv"))


dat$Accession <- sub("\\|.*", "", dat$Accession)                   # Extract first three characters

#qualitative analysis
dat1 <- dat[,c(3,50:52, 56:58, 68:72, 76:84, 88:93)]
head(dat1)
head(dat1)

str(dat1)
rownames(dat1) <- dat1[,1]
dat1 <- dat1[,-1]
head(dat1)

## Qualititative analysis
#Extraction group-specific proteins
Cre <- dat1[which(rowMeans(!is.na(dat1[,c(1:3, 9:11, 18:20)])) >= 9/9), ]
DMSO <- dat1[which(rowMeans(!is.na(dat1[,c(4:6, 12:14, 21:23)])) >= 9/9), ]
Dif <- dat1[which(rowMeans(!is.na(dat1[,c(7:8, 15:17, 24:26)])) >= 9/9), ]

#Venn diagram
library(VennDiagram)
library(RColorBrewer)
myCol1 <- brewer.pal(3, "Pastel2")

#This code will draw diagram to working directory
venn.diagram(
  x = list(rownames(Cre), rownames(DMSO), rownames(Dif)), 
  category.names = c("Cre" , "DMSO" , "Dif"),
  filename = '#venn_diagramm.png',
  resolution = 600,
  fill = myCol1,
  output=F
)

#Names of proteins specific for each group
library(gplots)
v.table2 <- venn(list(rownames(Cre), rownames(DMSO), rownames(Dif)))
print(v.table2)

#sample_info
library(readxl)
fact <- data.frame(read_excel("sample_info.xlsx"))

rownames(fact) <- fact[,1]
fact <- fact[,-1]
fact$Type <- as.factor(fact$Type)
fact$Type

fact$Donor <- as.factor(fact$Donor)
fact$Donor



#Cre_vs_DMSO+contr
dat2 <- dat[,c(3,50:52, 56:67, 70:72, 76:78, 82:84, 88:90)]
head(dat2)
colnames(dat2) <- sub("Area.", "", colnames(dat2))                   # Extract first three characters
head(dat2)

str(dat2)
rownames(dat2) <- dat2[,1]
dat2 <- dat2[,-1]
head(dat2)


#Removing rows with a lot of missing values
dat2 <- dat2[which(rowMeans(!is.na(dat2)) >= 0.85), ]
mean(complete.cases(dat2))
colSums(is.na(dat2))

#knn imputation of missng values
library(impute)
tdat <- t(dat2)
dat_knn1 <- impute.knn(tdat, k = 5)
dat_knn <- t(dat_knn1$data)
mean(complete.cases(dat_knn))


#Normalization and data QC
library(RColorBrewer)
#tiff('Raw_dat.tiff', units="in", width=16, height=8, res=600, compression = 'lzw')
pal <- brewer.pal(n = 9, name = "Set1")
cols <- pal[fact$Type]
boxplot(dat_knn, outline = FALSE, col = cols, main = "Raw data")
legend("topright", levels(fact$Type), fill = pal, bty = "n", xpd = T)
#dev.off()
colSums(dat_knn)
#log transformation
dat_log <- log2(dat_knn+1)
head(dat_log)
mean(complete.cases(dat_log))
boxplot(dat_log, outline = FALSE, col = cols, main = "Log-transformed data")
legend("topright", levels(fact$Type), fill = pal, bty = "n", xpd = T)

#Quantile normalization
library(limma)
dat_norm <- normalizeQuantiles(dat_log)
head(dat_norm)
#tiff('Norm_dat.tiff', units="in", width=16, height=8, res=300, compression = 'lzw')
boxplot(dat_norm, col = cols, main = "Normalized data")
legend("topright", levels(fact$Type), fill = pal, bty = "n", xpd = T)
#dev.off()
colSums(is.na(dat_norm))
mean(complete.cases(dat_norm))


#Removing batch effect by Combat
library(sva)
pheno <- fact
edata <- dat_norm
colnames(pheno)[2] <- "batch"
batch = pheno$batch
mod <- model.matrix(~ Type, data = pheno)

# non-parametric adjustment, mean-only version
combat_edata2 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=FALSE, mean.only=TRUE)


## Clusterisation
t_dat1 <- t(combat_edata2)

library(mixOmics)
#PLS-DA 
ordination.optimum.splsda <- splsda(t_dat1, fact$Type, ncomp = 3, keepX = c(15,15,15))
selectVar(ordination.optimum.splsda, comp=1)
selectVar(ordination.optimum.splsda, comp=2)
selectVar(ordination.optimum.splsda, comp=3)

#tiff('PLSDA_cre_contr.tiff', units="in", width=16, height=8, res=600, compression = 'lzw')
layout(matrix(c(1, 2, 3, 3, 3, 3), 2, 3))
plotLoadings(ordination.optimum.splsda, comp = 1, size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", contrib = "max", legend = FALSE, col.ties="black", ndisplay = 15)
plotLoadings(ordination.optimum.splsda, comp = 2, size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", contrib = "max",ndisplay = 15,  legend = FALSE, col.ties="black")
plotIndiv(ordination.optimum.splsda, ind.names = F, ellipse = T, style = "graphics", abline = TRUE, cex = 1.5, size.axis = 1.2, size.xlabel = 1.5, size.ylabel = 1.5, title = "PLS-DA ordination", size.title = 1.5, legend=TRUE)
#dev.off()
layout(1,1)

dat_pca <- pca(t(as.data.frame(combat_edata2)), ncomp = 3, center = TRUE)
#tiff('PCA_type_all.tiff', units="in", width=16, height=8, res=600, compression = 'lzw')
plotIndiv(dat_pca, comp = c(1, 2), ind.names = T, 
          group = fact$Type, legend = TRUE, ellipse = TRUE,
          title = 'PCA')
#dev.off()


#Dif_Cre_vs_dif+DMSO. Cre-ref
Xc <- model.matrix(~ 0+fact$Type)
Xc
colnames(Xc) <- c("Control","Cre","DMSO")

fitc1 <- lmFit(dat_norm, design = Xc, method = "robust", maxit = 10000)


desC <- makeContrasts((Cre-(DMSO+Control)/2),
                      levels=fact$Type)
fitc2 <- contrasts.fit(fitc1, desC)


# Empirical Bayes statistics
efit_c <- eBayes(fitc2)

# Dif_expr_table
topTable(efit_c, coef = 1)
full_list_c <- topTable(efit_c, coef = 1, number = length(dat_norm))
#write.csv(full_list_c,'Dif_expr_DMSO+contr_vs_Cre_cor.csv')
head(full_list_c)

vulk_dat <- data.frame(read_xlsx("dif_expre_supl.xlsx", sheet = 2 )) #dif_expre_supl.xlsx is table with dif. expr. results and gene names added
head(vulk_dat)
str(vulk_dat)
vulk_dat$logFC <- as.numeric(vulk_dat$logFC)
vulk_dat$adj.P.Val <- as.numeric(vulk_dat$adj.P.Val)
head(vulk_dat)
rownames(vulk_dat) = make.names(vulk_dat[,2], unique=TRUE)

#Vulcano
library(EnhancedVolcano)

#tiff('Vulcano_cre_vsDMSOContr_cor.tiff', units="in", width=11, height=8, res=300, compression = 'lzw')
EnhancedVolcano(vulk_dat,
                lab = rownames(vulk_dat),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 1,
             # xlim = c(-3, 5),
              #  ylim = c(0, 10),
                title ="Cre versus DMSO+Contr",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)
#dev.off()

#Manually check top differentialy expressed proteins
boxplot(combat_edata2[c("Q9Y3Z3"),] ~ Type, data = fact,
        varwidth = TRUE, log = "y", las = 1)
boxplot(dat_knn[c("P09601"),] ~ Type, data = fact,
        varwidth = TRUE, log = "y", las = 1)


#Pathway enrichment analysis by rWikiPathways according to http://www.bioconductor.org/packages/devel/bioc/vignettes/rWikiPathways/inst/doc/Pathway-Analysis.html
library(rWikiPathways)

load.libs <- c(
  "DOSE",
  "GO.db",
  "GSEABase",
  "org.Hs.eg.db",
  "clusterProfiler",
  "dplyr",
  "tidyr",
  "ggplot2",
  "stringr",
  "RColorBrewer",
  "rWikiPathways",
  "RCy3")
library(pacman)

head(vulk_dat)
vulk_dat <- vulk_dat[,-1]

up.genes <- vulk_dat[vulk_dat$logFC > 1 & vulk_dat$adj.P.Val < 0.05, 1] 
dn.genes <- vulk_dat[vulk_dat$logFC < -1 & vulk_dat$adj.P.Val < 0.05, 1]
bkgd.genes <- vulk_dat[,1]

up.genes.entrez <- clusterProfiler::bitr(up.genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
dn.genes.entrez <- bitr(dn.genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
bkgd.genes.entrez <- bitr(bkgd.genes,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)

egobp <- clusterProfiler::enrichGO(
  gene     = up.genes.entrez[[2]],
  universe = bkgd.genes.entrez[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
  readable = TRUE)

head(egobp,10)
barplot(egobp, showCategory = 20)
dotplot(egobp, showCategory = 20)

egobp <- clusterProfiler::enrichGO(
  gene     = dn.genes.entrez[[2]],
  universe = bkgd.genes.entrez[[2]],
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  pAdjustMethod = "fdr",
  pvalueCutoff = 0.05, #p.adjust cutoff (https://github.com/GuangchuangYu/clusterProfiler/issues/104)
  readable = TRUE)

head(egobp,20)
barplot(egobp, showCategory = 30)
dotplot(egobp, showCategory = 20)

#tiff('enrichment_plot.tiff', units="in", width=8, height=6, res=300, compression = 'lzw')
dotplot(egobp, showCategory = 20)
#dev.off()
