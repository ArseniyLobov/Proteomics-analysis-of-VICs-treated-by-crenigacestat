###This is a code to reproduce statistical analysis of proteomics data from Lobov et al., 2021 "Crenigacestat (LY3039478) inhibits osteogenic differentiation of human valve interstitial cells and seems promising for treatment of cardiovascular calcification"

##Openinig_the_data
#We reccomed to set the woring directory to make easy to reproduce the code
#setwd("your directory")
dat <- data.frame(read.csv("data.csv"))
dat$Accession <- sub("\\|.*", "", dat$Accession)                   

#exctracting the data used in the Lobov et al., 2021
dat1 <- dat[,c(3,50:52, 56:58, 68:72, 76:84, 88:93)]
head(dat1)
colnames(dat1) <- sub("Area.", "", colnames(dat1))                   
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

## Quantitative analysis
#Removing rows with a lot of missing values
dat2 <- dat1[which(rowMeans(!is.na(dat1)) >= 0.85), ]
mean(complete.cases(dat2))
colSums(is.na(dat2))

#knn imputation of missng values
library(impute)
tdat <- t(dat2)
dat_knn1 <- impute.knn(tdat, k = 5)
dat_knn <- t(dat_knn1$data)
mean(complete.cases(dat_knn))

#Opening the factor matrix
library(readxl)
fact <- data.frame(read_excel("Legend.xlsx"))

rownames(fact) <- fact[,1]
fact <- fact[,-1]
fact$Type <- as.factor(fact$Type)
fact$Type

fact$Donor <- as.factor(fact$Donor)
fact$Donor


#Normalization and data QC
library(RColorBrewer)
#tiff('Raw_dat.tiff', units="in", width=16, height=8, res=600, compression = 'lzw')
pal <- brewer.pal(n = 9, name = "Set1")
cols <- pal[fact$Type]
boxplot(dat_knn, outline = FALSE, col = cols, main = "Raw data")
legend("topright", levels(fact$Type), fill = pal, bty = "n", xpd = T)
#dev.off()
colSums(dat_knn)
#log transformation (погуглите зачем это делать). В чем разница с тем, что было? почему так удобно? зачем прибавлять 1 к dat_knn?
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


#Differential expression analysis: Cre_vs_dif+DMSO. Cre-ref
Xc <- model.matrix(~ 0+fact$Type)
Xc
colnames(Xc) <- c("Cre","Dif","DMSO")

fitc1 <- lmFit(dat_norm, design = Xc, method = "robust", maxit = 10000)

desC <- makeContrasts((Cre-(DMSO+Dif)/2),
              levels=fact$Type)
fitc2 <- contrasts.fit(fitc1, desC)

# Empirical Bayes statistics
efit_c <- eBayes(fitc2)

# Dif_expr_table
topTable(efit_c, coef = 1)
full_list_c <- topTable(efit_c, coef = 1, number = length(dat_norm))
#write.csv(full_list_c,'Dif_expr_DMSO+Dif_vs_Cre.csv')
head(full_list_c)


#Vulcano
library(EnhancedVolcano)

#tiff('Vulcano_cre_vsDMSOdiff.tiff', units="in", width=11, height=8, res=300, compression = 'lzw')
EnhancedVolcano(full_list_c,
                lab = rownames(full_list_c),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 1,
                #xlim = c(-2.5, 4),
                #ylim = c(0, 10),
                title ="Cre versus DMSO+Diff",
                labSize = 4.0,
                boxedLabels = F,
                colAlpha = 1)
#dev.off()



#Manually check top differentially expressed proteins
boxplot(combat_edata2[c("P09601"),] ~ Type, data = fact,
        varwidth = TRUE, log = "y", las = 1)
boxplot(dat_knn[c("P09601"),] ~ Type, data = fact,
        varwidth = TRUE, log = "y", las = 1)


boxplot(combat_edata2[c("P29536"),] ~ Type, data = fact,
        varwidth = TRUE, log = "y", las = 1)
boxplot(dat_knn[c("P29536"),] ~ Type, data = fact,
        varwidth = TRUE, log = "y", las = 1)


## Ordination by sPLS-DA
t_dat1 <- t(combat_edata2)

library(mixOmics)
#PLS-DA 
ordination.optimum.splsda <- splsda(t_dat1, fact$Type, ncomp = 3, keepX = c(15,15,15))
selectVar(ordination.optimum.splsda, comp=1)
selectVar(ordination.optimum.splsda, comp=2)
selectVar(ordination.optimum.splsda, comp=3)

#tiff('PLSDA_cre_dif.tiff', units="in", width=16, height=8, res=600, compression = 'lzw')
layout(matrix(c(1, 2, 3, 3, 3, 3), 2, 3))
plotLoadings(ordination.optimum.splsda, comp = 1, size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", contrib = "max", legend = FALSE, col.ties="black", ndisplay = 15)
plotLoadings(ordination.optimum.splsda, comp = 2, size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", contrib = "max",ndisplay = 15,  legend = FALSE, col.ties="black")
plotIndiv(ordination.optimum.splsda, ind.names = F, ellipse = T, style = "graphics", abline = TRUE, cex = 1.5, size.axis = 1.2, size.xlabel = 1.5, size.ylabel = 1.5, title = "PLS-DA ordination", size.title = 1.5, legend=TRUE)
#dev.off()
layout(1,1)