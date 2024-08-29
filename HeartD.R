
## 1.1 DEG1 ####
setwd("/Users/emperor/Documents/DataR/HeartD/GSE20681/")
GSE20681_series_matrix <- read.table("GSE20681_series_matrix.txt.gz",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE)
GPL4133 <- read.table("GPL4133-12599.txt",sep = "\t",comment.char = "#", stringsAsFactors = F,header = T, fill=TRUE, quote="")
GPL4133 = data.frame(GPL4133[,c(1,10)])

trim <- function( x ) {
  gsub(" ", "", x)
}
GPL4133$ID = trim(GPL4133$ID)

Exp <- merge(GSE20681_series_matrix,GPL4133,by.x = "ID_REF",by.y = "ID",all=FALSE) # 45015
Exp = Exp[,c(200,2:199)]
colnames(Exp)[1] = "Symbol"
# Exp = data.frame(Exp[-grep("///",Exp$"Symbol"),]) #去除一个探针对应多个Symbol的行(去一对多???
Exp <- Exp[Exp$Symbol != "" ,]  # 3471
# Exp = na.omit(Exp)
meanfun <- function(x) {
  x1 <- data.frame(unique(x[,1]))
  colnames(x1) <- c("Symbol")
  for (i in 2:ncol(x)){
    x2 <- data.frame(tapply(x[,i],x[,1],mean))
    x2[,2] <- rownames(x2)
    colnames(x2) <- c(colnames(x)[i], "Symbol")
    x1 <- merge(x1,x2,by.x = "Symbol",by.y = "Symbol",all=FALSE)
  }
  return(x1)
}
Exp <- meanfun(Exp) # 19749
write.table(Exp,"Exp1.txt",row.names = F,quote = F,sep="\t")








## 在全基因全样本表达谱里筛??? Case ??? Control 样本
sample_sheet <- read.table("./sample1.txt",sep = "\t",comment.char = "#", stringsAsFactors = F,header = T, fill=TRUE, quote="")

Case_sample_type = data.frame(sample_sheet[grep("case",sample_sheet$Group),])
Case_sample_type = Case_sample_type[1:99,1]

Control_sample_type = data.frame(sample_sheet[grep("control",sample_sheet$Group),])
Control_sample_type = Control_sample_type[1:99,1]

numData = Exp
rownames(numData)=numData[,1]
numData = numData[,-1]

# 筛选疾病样???
Casesample_nindex=match(Case_sample_type,colnames(numData))
CasesampleMatrix = numData[,Casesample_nindex];
CasesampleMatrix = data.frame(Symbol=rownames(CasesampleMatrix),CasesampleMatrix)
write.table(CasesampleMatrix, 'CasesampleMatrix.txt',sep = '\t', quote=F, row.names=F)

# 筛选正常样???
Controlsample_nindex=match(Control_sample_type,colnames(numData))
ControlsampleMatrix = numData[,Controlsample_nindex];
ControlsampleMatrix = data.frame(Symbol=rownames(ControlsampleMatrix),ControlsampleMatrix)
write.table(ControlsampleMatrix, 'ControlsampleMatrix.txt',sep = '\t', quote=F, row.names=F)


## 合并
DEG_EXP = merge(CasesampleMatrix,ControlsampleMatrix,by.x = "Symbol",by.y = "Symbol",all=FALSE)
write.table(DEG_EXP, 'DEG_EXP.txt',sep = '\t', quote=F, row.names=F)



## DEG (P.Value < 0.05)
library(limma)
# 把ExprData_result.txt编辑成rep("Control",5),rep("asthma",5))形式
exp = DEG_EXP
for (i in 2:ncol(exp)){
  exp[,i] = as.numeric(exp[,i])
}
# 开始差异分???
rownames(exp)<-exp[,1]
exp<-exp[,-1]
exp <- 16*exp
# exp <- log2(exp)
# exp = exp[,1:12]
exp<-as.matrix(exp)

samps<-factor(c(rep("case",ncol(CasesampleMatrix)-1),rep("control",ncol(ControlsampleMatrix)-1)))
design <- model.matrix(~0+samps);
colnames(design) <- c("case","control")
fit <- lmFit(exp,design)
cont.matrix<-makeContrasts(case-control,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
final<-topTable(fit2, coef=1, number=dim(exp)[1], adjust.method="BH", sort.by="B", resort.by="M")

DEG = data.frame(DEG=rownames(final),final)
DEG_sort_p <- DEG[DEG$P.Value < 0.05 & abs(DEG$logFC)>1.5  ,] # 339
write.table(DEG,"DegData_limma1.txt",quote=FALSE,sep="\t",row.names=F)
write.table(DEG_sort_p, 'DEG_sort_p1.txt',sep = '\t', quote=F, row.names=F)

DE_Gene = data.frame(DEG_sort_p[,1])
colnames(DE_Gene) = "DEG"
write.table(DE_Gene, 'DE_Gene1.txt',sep = '\t', quote=F, row.names=F)

a = DE_Gene

## 1.1.1 热图  ####
EXPexp <- read.table('Exp1.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
DEG_sort_p <- read.table('DEG_sort_p1.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
consensusClass <- read.table('./sample1.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40

# 注释的文件，记得行名是TCGA.ZB.A96A.01A这些
annotation <- data.frame(cluster = consensusClass[, 3])
rownames(annotation) <- consensusClass[, 1]

library(gplots)
EXP.data <- EXPexp
da <- DEG_sort_p

Deg_exprData <- EXP.data[EXP.data$Symbol %in% da$DEG, ]
rownames(Deg_exprData) <- Deg_exprData[,1]
Deg_exprData <- Deg_exprData[,-1]

# 排序，注释文件的排序ok，merge到exp上，用排序后的exp
Deg_exprData3 <- Deg_exprData
Deg_exprData3 <- t(Deg_exprData3)
Deg_exprData3 <- data.frame(id=rownames(Deg_exprData3), Deg_exprData3)
# Deg_exprData3P <- cbind(annotation, Deg_exprData3)
Deg_exprData3P <- merge(data.frame(id=rownames(annotation), annotation), Deg_exprData3, by = "id", all = F)

# 矩阵排序
Deg_exprData3P <- Deg_exprData3P[order(Deg_exprData3P$cluster, decreasing = T), ]
Deg_exprData3Q <- Deg_exprData3P[, -2] # 去掉cluster1 ，2 
rownames(Deg_exprData3Q) <- Deg_exprData3Q[, 1]
Deg_exprData3Q <- Deg_exprData3Q[, -1]
Deg_exprData3QQ <- scale(Deg_exprData3Q)
Deg_exprData3Q <- data.frame(Deg_exprData3Q)
library(pheatmap)
pheatmap(t(Deg_exprData3QQ),cluster_rows = T,
         show_rownames = F,
         cluster_cols = F,
         show_colnames = F, 
         annotation_col = annotation,
         scale = "row"
         ) # column 




## 1.1.1-火山图 ####
library(ggplot2)
library(Cairo)
DEG <- read.table("DegData_limma1.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE,check.names=F)
title <- "mRNA"
filePath <- "Volcano picture of DEMs.pdf"

data <- DEG
# 设置颜色域
data$threshold <- as.factor(ifelse(data$P.Value < 0.05 & abs(data$logFC) >= 1.5,ifelse(data$logFC > 1.5 ,'Up','Down'),'Not'))
# 说明：注意logFC的参数，正常用1.5，若实验倍数差异不明显，可使用0.5（后面参数同样修改）
##Construct the plot object
# with legend
# 生成文件
# Cairo(file="Volcan1.pdf", type="png",units="in",bg="white",width=5.5, height=5, pointsize=12, dpi=300)
plot2 = ggplot(data=data, aes(x=logFC, y =-log(P.Value), colour=threshold,fill=threshold)) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_point(alpha=0.4, size=1.2) +
  xlim(c(-5, 3)) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept=c(-1.5, 1.5),lty=4,col="grey",lwd=0.6)+
  geom_hline(yintercept = -log(0.05),lty=4,col="grey",lwd=0.6)+
  theme(legend.position="right",
        panel.grid=element_blank(),
        legend.title = element_blank(),
        legend.text= element_text(face="bold", color="black",family = "Times", size=8),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold", color="black", size=12),
        axis.text.y = element_text(face="bold",  color="black", size=12),
        axis.title.x = element_text(face="bold", color="black", size=12),
        axis.title.y = element_text(face="bold",color="black", size=12))+
  labs(x="log2 (fold change)",y="-log (p_value)",title= title)
plot2
ggsave(filePath, plot2, width = 5.15, height = 4.65)



# 教程：https://cloud.tencent.com/developer/article/2441324
# ComBat ####
library(limma)
library(sva)
library(ggplot2)
library(FactoMineR) # PCA函数
library(factoextra) # fviz_pca_ind函数
Exp1 <- read.table("./Exp1.txt",sep = "\t",comment.char = "#", stringsAsFactors = F,header = T, fill=TRUE, quote="")
Exp2 <- read.table("../GSE42148/Exp2.txt",sep = "\t",comment.char = "#", stringsAsFactors = F,header = T, fill=TRUE, quote="")
xxgene <- Exp1[Exp1$Symbol %in% Exp2$Symbol, ]
Exp1x <-  Exp1[Exp1$Symbol %in% xxgene$Symbol, ]
Exp2x <-  Exp2[Exp2$Symbol %in% xxgene$Symbol, ]
Expxxx <- merge(Exp1x, Exp2x, by = "Symbol", all = F)
sample1 <- read.table("./sample1.txt",sep = "\t",comment.char = "#", stringsAsFactors = F,header = T, fill=TRUE, quote="")
sample1$batch  <- "1"
sample2 <- read.table("../GSE42148/sample2.txt",sep = "\t",comment.char = "#", stringsAsFactors = F,header = T, fill=TRUE, quote="")
sample2$batch  <- "2"
sample <- rbind(sample1, sample2)


# PCA分析--------------------------------------------------
# PCA图时要求是行名时样本名，列名时探针名，因此此时需要转换
exp1 <- Expxxx
rownames(exp1) <- exp1[, 1]
exp1 <- exp1[, -1]
exp2 <- t(exp1)
exp2 <- data.frame(exp2)
exp2[1:4, 1:4]
#              1007_s_at  1053_at   117_at   121_at
# GSM71019.CEL 10.115170 5.345168 6.348024 8.901739
# GSM71020.CEL  8.628044 5.063598 6.663625 9.439977
# GSM71021.CEL  8.779235 5.113116 6.465892 9.540738
# GSM71022.CEL  9.248569 5.179410 6.116422 9.254368


# 设置批次信息
batch <- sample$batch # 批次
# 设置生物学分类，告诉函数不要把生物学差异整没了 
sample$Group <- factor(sample$Group, levels = c("control", "case"))
mod <- model.matrix(~as.factor(Group), data=sample)

# x需要注意基因要相同，且基因在列名中，整个矩阵是数字
expr_combat <- ComBat(dat = exp1, batch = batch, mod = mod,par.prior=TRUE, ref.batch=1)
expr_limma <- removeBatchEffect(exp1,batch=batch,batch2=NULL,
                                covariates=NULL,design= mod)

expr_limma <- data.frame(expr_limma)
exp3 <- t(expr_limma)
exp3 <- data.frame(exp3)
exp3[1:4, 1:4]
#              A1BG     A1CF    A2BP1    A2LD1
# GSM518885 8.239918 5.145095 5.086594 7.719400
# GSM518886 8.576959 5.175852 5.095249 8.084000
# GSM518887 8.532465 5.146992 5.097216 7.945976
# GSM518888 8.566251 5.150919 5.163170 8.337744

# 添加分组信息
# ac <- data.frame(
#   row.names = rownames(exp), 
#   Group = pheno$cancer)
ac <- sample[, c(1,3)]


# 绘图
dat.pca <- PCA(exp3, graph = FALSE)
p.pca <- fviz_pca_ind(dat.pca,
                      # 只显示点而不显示文本，默认都显示
                      geom.ind = "point",  #c( "point", "text" ) / "point",
                      # geom.ind = 'text',
                      # 设定分类种类
                      col.ind = sample$batch, # ac$Group，color by groups
                      # 设定颜色
                      palette = "jco", # jco/Dark2
                      # 添加椭圆
                      addEllipses = TRUE, # Concentration ellipses
                      # 添加图例标题
                      legend.title = "Groups") +
  # ggtitle(this_title) + 
  theme(plot.title = element_text(size=16, hjust = 0.5))
p.pca


# 单独
# Exp1x, sample1
expSingle1 <- Expxxx
rownames(expSingle1) <- expSingle1[, 1]
expSingle1 <- expSingle1[, -1]
expSingle11 <- t(expSingle1)
expSingle11 <- data.frame(expSingle11)
expSingle11[1:4, 1:4]

ac <- sample

# 绘图
dat.pca <- PCA(expSingle11, graph = FALSE)
p.pca <- fviz_pca_ind(dat.pca,
                      # 只显示点而不显示文本，默认都显示
                      geom.ind = "point",  #c( "point", "text" ) / "point",
                      # geom.ind = 'text',
                      # 设定分类种类
                      col.ind = ac$batch, # ac$Group，color by groups
                      # 设定颜色
                      palette = "jco", # jco/Dark2
                      # 添加椭圆
                      addEllipses = TRUE, # Concentration ellipses
                      # 添加图例标题
                      legend.title = "Groups") +
  # ggtitle(this_title) + 
  theme(plot.title = element_text(size=16, hjust = 0.5))
p.pca



## 1.2 DEG2 ####
setwd("/Users/emperor/Documents/DataR/HeartD/GSE42148/")
GSE42148_series_matrix <- read.table("GSE42148_series_matrix.txt.gz",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE)
GPL13607 <- read.table("GPL13607-20416.1.txt",sep = "\t",comment.char = "#", stringsAsFactors = F,header = T, fill=TRUE, quote="")
GPL13607 = data.frame(GPL13607[,c(1,6)])

trim <- function( x ) {
  gsub(" ", "", x)
}
GPL13607$ID = trim(GPL13607$ID)

Exp <- merge(GSE42148_series_matrix,GPL13607,by.x = "ID_REF",by.y = "ID",all=FALSE) # 45015
Exp = Exp[,c(26,2:25)]
colnames(Exp)[1] = "Symbol"
# Exp = data.frame(Exp[-grep("///",Exp$"Symbol"),]) #去除一个探针对应多个Symbol的行(去一对多???
Exp <- Exp[Exp$Symbol != "" ,]  # 3471
# Exp = na.omit(Exp)
meanfun <- function(x) {
  x1 <- data.frame(unique(x[,1]))
  colnames(x1) <- c("Symbol")
  for (i in 2:ncol(x)){
    x2 <- data.frame(tapply(x[,i],x[,1],mean))
    x2[,2] <- rownames(x2)
    colnames(x2) <- c(colnames(x)[i], "Symbol")
    x1 <- merge(x1,x2,by.x = "Symbol",by.y = "Symbol",all=FALSE)
  }
  return(x1)
}
Exp <- meanfun(Exp) # 34949
write.table(Exp,"Exp2.txt",row.names = F,quote = F,sep="\t")

## 在全基因全样本表达谱里筛??? Case ??? Control 样本
sample_sheet <- read.table("./sample2.txt",sep = "\t",comment.char = "#", stringsAsFactors = F,header = T, fill=TRUE, quote="")

Case_sample_type = data.frame(sample_sheet[grep("case",sample_sheet$Group),])
Case_sample_type = Case_sample_type[1:13,1]

Control_sample_type = data.frame(sample_sheet[grep("control",sample_sheet$Group),])
Control_sample_type = Control_sample_type[1:11,1]

numData = Exp
rownames(numData)=numData[,1]
numData = numData[,-1]

# 筛选疾病样???
Casesample_nindex=match(Case_sample_type,colnames(numData))
CasesampleMatrix = numData[,Casesample_nindex];
CasesampleMatrix = data.frame(Symbol=rownames(CasesampleMatrix),CasesampleMatrix)
write.table(CasesampleMatrix, 'CasesampleMatrix.txt',sep = '\t', quote=F, row.names=F)

# 筛选正常样???
Controlsample_nindex=match(Control_sample_type,colnames(numData))
ControlsampleMatrix = numData[,Controlsample_nindex];
ControlsampleMatrix = data.frame(Symbol=rownames(ControlsampleMatrix),ControlsampleMatrix)
write.table(ControlsampleMatrix, 'ControlsampleMatrix.txt',sep = '\t', quote=F, row.names=F)


## 合并
DEG_EXP = merge(CasesampleMatrix,ControlsampleMatrix,by.x = "Symbol",by.y = "Symbol",all=FALSE)
write.table(DEG_EXP, 'DEG_EXP.txt',sep = '\t', quote=F, row.names=F)



## DEG (P.Value < 0.05)
library(limma)
# 把ExprData_result.txt编辑成rep("Control",5),rep("asthma",5))形式
exp = DEG_EXP
for (i in 2:ncol(exp)){
  exp[,i] = as.numeric(exp[,i])
}
# 开始差异分???
rownames(exp)<-exp[,1]
exp<-exp[,-1]
exp <- 2.8*exp
# exp <- log2(exp)
# exp = exp[,1:12]
exp<-as.matrix(exp)

samps<-factor(c(rep("case",ncol(CasesampleMatrix)-1),rep("control",ncol(ControlsampleMatrix)-1)))
design <- model.matrix(~0+samps);
colnames(design) <- c("case","control")
fit <- lmFit(exp,design)
cont.matrix<-makeContrasts(case-control,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
final<-topTable(fit2, coef=1, number=dim(exp)[1], adjust.method="BH", sort.by="B", resort.by="M")

DEG = data.frame(DEG=rownames(final),final)
DEG_sort_p <- DEG[DEG$P.Value < 0.05 & abs(DEG$logFC)>1.5  ,] # 1424
write.table(DEG,"DegData_limma2.txt",quote=FALSE,sep="\t",row.names=F)
write.table(DEG_sort_p, 'DEG_sort_p2.txt',sep = '\t', quote=F, row.names=F)

DE_Gene = data.frame(DEG_sort_p[,1])
colnames(DE_Gene) = "DEG"
write.table(DE_Gene, 'DE_Gene2.txt',sep = '\t', quote=F, row.names=F)

b= DE_Gene



## 1.2.1 热图 ####
EXPexp <- read.table('Exp2.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
DEG_sort_p <- read.table('DEG_sort_p2.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
consensusClass <- read.table('./sample2.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40

# 注释的文件，记得行名是TCGA.ZB.A96A.01A这些
annotation <- data.frame(cluster = consensusClass[, 3])
rownames(annotation) <- consensusClass[, 1]

library(gplots)
EXP.data <- EXPexp
da <- DEG_sort_p

Deg_exprData <- EXP.data[EXP.data$Symbol %in% da$DEG, ]
rownames(Deg_exprData) <- Deg_exprData[,1]
Deg_exprData <- Deg_exprData[,-1]

# 排序，注释文件的排序ok，merge到exp上，用排序后的exp
Deg_exprData3 <- Deg_exprData
Deg_exprData3 <- t(Deg_exprData3)
Deg_exprData3 <- data.frame(id=rownames(Deg_exprData3), Deg_exprData3)
Deg_exprData3P <- cbind(annotation, Deg_exprData3)
# Deg_exprData3P <- merge(data.frame(id=rownames(annotation), annotation), Deg_exprData3, by = "id", all = F)

# 矩阵排序
# Deg_exprData3P <- Deg_exprData3P[order(Deg_exprData3P$cluster, decreasing = T), ]
Deg_exprData3Q <- Deg_exprData3P[, -1] # 去掉cluster1 ，2 
rownames(Deg_exprData3Q) <- Deg_exprData3Q[, 1]
Deg_exprData3Q <- Deg_exprData3Q[, -1]

library(pheatmap)
pheatmap(t(Deg_exprData3Q),cluster_rows = T, # Deg_exprData3Q[, 500:1000]
         show_rownames = F,
         cluster_cols = F,
         show_colnames = F, 
         annotation_col = annotation,
         scale = "row"
) # column 



CFPgene <- read.table("../CFPgene.txt",sep = "\t",comment.char = "#", stringsAsFactors = F,header = T, fill=TRUE, quote="")
xxx <- merge(CFPgene,DEG_sort_p, by.x = "Symbol", by.y = "DEG", all = F )

# 2 WGCNA ####
setwd("/Users/emperor/Documents/DataR/HeartD/2WGCNA")
# setwd("/Users/emperor/Documents/DataR/HeartD/2WGCNA.revised/")

# 用GSE20681矩阵-样本量多
eXP <- read.table("../GSE20681/Exp1.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE,check.names=F)
# sur.clini.OK <- readRDS("./WGCNA.All/SurvivalInput.clini.OK.Rds")

library(WGCNA)

# 构建GSE42148数据库中2031个差异+铜、铁死亡pyrotosis等基因的所有共表达谱
DEG_sort_p1 <- read.table("../GSE20681//DEG_sort_p1.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE,check.names=F)
gene1 <- data.frame(symbol=DEG_sort_p1[, 1])

DEG_sort_p2 <- read.table("../GSE42148/DEG_sort_p2.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE,check.names=F)
gene2 <- data.frame(symbol=DEG_sort_p2[, 1])
CFPgene <- read.table("../CFPgene.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE,check.names=F)
gene3 <- data.frame(symbol=CFPgene[, 1])
gene <- rbind(gene1,gene2, gene3)
gene <- unique(gene) # 2232
DEGsAllSampleMatrix =eXP[eXP$Symbol %in% gene$symbol, ]
# DEGsAllSampleMatrix <- t(DEGsAllSampleMatrix)
# DEGsAllSampleMatrix <- data.frame(gene=rownames(DEGsAllSampleMatrix), DEGsAllSampleMatrix)

write.table(DEGsAllSampleMatrix, './DEGsAllSampleMatrix.revised.txt',sep = '\t', quote=F, row.names=F)

## 开始共表达分析
# source("https://bioconductor.org/biocLite.R")
# biocLite(c("AnnotationDbi", "impute","GO.db", "preprocessCore"))
# site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# install.packages(c("WGCNA", "stringr", "reshape2"), repos=site)

library(WGCNA)
library(dynamicTreeCut)
library(fastcluster)

## 输入表达矩阵
# exprMat <- "DEGsCaseSampleMatrix.txt"
exprMat <- "DEGsAllSampleMatrix.revised.txt"
type = "unsigned"
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)



# 对二元变量，如样本性状信息计算相关性时，
# 或基因表达严重依赖于疾病状态时，需设置下面参数
maxPOutliers = ifelse(corType=="pearson",1,0.05)

# 关联样品性状的二元变量时，设置
robustY = ifelse(corType=="pearson",T,F)

##导入数据##
dataExpr <- read.table(exprMat, sep='\t',  header=T, 
                       quote="", comment="", check.names=F)
# 去除重复值，去除重复名字
# dataExpr <- dplyr::distinct(dataExpr, Symbol, .keep_all = T)
rownames(dataExpr) <- dataExpr[, 1]
dataExpr <- dataExpr[, -1]

dim(dataExpr)  # 1237  198
head(dataExpr)[,1:8]
# TCGA.5L.AAT0.01A

## 数据筛选
## 筛选中位绝对偏差前75%的基因，至少MAD大于0.01
## 筛选后会降低运算量，也会失去部分信息
## 也可不做筛选，使MAD大于0即可
# apply函数只能用于处理矩阵类型的数据，也就是说所有的数据必须是同一类型。因此要使用apply函数的话，需要将数据类型转换成矩阵类型。
# apply函数一般有三个参数，第一个参数代表矩阵对象，第二个参数代表要操作矩阵的维度，1表示对行进行处理，2表示对列进行处理。第三个参数就是处理数据的函数。apply会分别一行或一列处理该矩阵的数据。
# MAD（Median absolute deviation, 中位数绝对偏差）是单变量数据集中样本差异性的稳健度量。mad是一个健壮的统计量，对于数据集中异常值的处理比标准差更具有弹性，可以大大减少异常值对于数据集的影响。
# m.mad <- apply(dataExpr,1,mad)
# dataExprVar <- dataExpr[which(m.mad > 
# max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]

# 转换为样品在行，基因在列的矩阵
# dataExpr <- as.data.frame(t(dataExprVar))

dataExpr <- as.data.frame(t(dataExpr))
for (i in 1:ncol(dataExpr)){
  dataExpr[,i] = as.numeric(dataExpr[,i])
}

## 检测缺失值(剔除基因)
gsg = goodSamplesGenes(dataExpr, verbose = 3);
gsg$allOK
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(dataExpr)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}
## Removing genes: XBP1, KREMEN1, PTPRCAP, ANKRD36BP1, TOMM6, SDHAP2, TNFRSF6B, PIGY, SLC7A5P2, MIR4461


## 查看是否有离群样品
sampleTree = hclust(dist(dataExpr), method = "average");
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))#这是为什么？
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
# 保存图片 Sample clustering to detect outliers



abline(h = 8, col = "red");#
clust = cutreeStatic(sampleTree, cutHeight = 100, minSize = 5)
table(clust)
# clust
# 0   1   2   3 
# 20 786 154   5   

# keepSamples = (clust==1) ## 剔除离群样本
# dataExpr = dataExpr[keepSamples, ]
nGenes = ncol(dataExpr)  # 1237
nSamples = nrow(dataExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(dataExpr, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.95;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=-0.9,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
abline(h=0.1,col="red")
# 保存图片 Scale independence & Mean connectivity

# 根据图表选择 power=20
power = 12

## 一般不使用推荐power
# power = sft$powerEstimate
# power
# [1] 1



## 经验power (无满足条件的power时选用)
# 无向网络在power小于15或有向网络power小于30内，没有一个power值可以使
# 无标度网络图谱结构R^2达到0.8，平均连接度较高如在100以上，可能是由于
# 部分样品与其他样品差别太大。这可能由批次效应、样品异质性或实验条件对
# 表达影响太大等造成。可以通过绘制样品聚类查看分组信息和有无异常样品。
# 如果这确实是由有意义的生物变化引起的，也可以使用下面的经验power值。
if (is.na(power)){
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                               ifelse(type == "unsigned", 6, 12))
                 )
  )
}  # power=20




### 网络构建
## 一步法网络构建：One-step network construction and module detection##
# power: 上一步计算的软阈值
# maxBlockSize: 计算机能处理的最大模块的基因数量 (默认5000)；
#  4G内存电脑可处理8000-10000个，16G内存电脑可以处理2万个，32G内存电脑可
#  以处理3万个
#  计算资源允许的情况下最好放在一个block里面。
# corType: pearson or bicor
# numericLabels: 返回数字而不是颜色作为模块的名字，后面可以再转换为颜色
# saveTOMs：最耗费时间的计算，存储起来，供后续使用
# mergeCutHeight: 合并模块的阈值，越大模块越少
# 必须做数值化
net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize =28,
                       reassignThreshold = 0, mergeCutHeight = 0.15,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType, 
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       saveTOMFileBase = "expr_d5_TOM",
                       verbose = 3)


# 根据模块中基因数目的多少，降序排列，依次编号为 `1-最大模块数`。
# **0 (grey)**表示**未**分入任何模块的基因。 
table(net$colors)
# 0    1    2    3    4 
# 1426  193   94   62   58   

# 0   1   2   3 
# 967 122 109  39 
length(net$colors)  # 1237
length(table(net$colors))  # 4
mergedColors = labels2colors(net$colors)
table(mergedColors)
mergedColors

# blue     brown      grey turquoise     
# 109        39      967       122        
write.table(table(mergedColors), 'Number_color.txt',sep = '\t', quote=F, row.names=F)

## 层级聚类树展示各个模块
## 灰色的为**未分类**到模块的基因。
# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
# Plot the dendrogram and the module colors underneath
# 如果对结果不满意，还可以recutBlockwiseTrees，节省计算时间
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# 保存图片 Cluster Dendrogram

# 解决聚类图和cut阈值聚类图 ####
# 合并阈值之后，再做一张这个图，借用别人的代码和自己的合并一起，弄清输入文件net$dendrograms[[1]]就是别人的geneTree，做出cut图
merge_modules.my = mergeCloseModules(dataExpr, moduleColors,
                                     cutHeight = 0.02,  # 关键，决定新旧模块数量情况
                                     verbose = 3)  # multiExpr[[1]]$data
# 合并后的颜色：
mergedColors.my = merge_modules.my$colors;
# 新模块的特征向量基因：
mergedMEs.my = merge_modules.my$newMEs;
plotDendroAndColors(net$dendrograms[[1]], cbind(moduleColors, mergedColors.my),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# 心得，这里的net$dendrograms[[1]]就是geneTree
# adjacency = adjacency(dataExpr, power = power)
# TOM = TOMsimilarity(adjacency) 
# dissTOM = 1-TOM
# geneTree = hclust(as.dist(dissTOM), method = "average");
length(unique(moduleColors))
length(unique(mergedColors.my))


# 解决METree，新旧模块对比-注意使用旧的模块颜色moduleColors，还是新的合并后的模块颜色mergedColors.my ####
# 新模块内容：
MEList = moduleEigengenes(dataExpr, 
                          colors = moduleColors  # 旧模块颜色 moduleColors
)  # dataExpr = multiExpr[[1]]$data
MEs = MEList$eigengenes
# 计算根据模块特征向量基因计算模块相异度：
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");

# 新模块情况看看
plot(METree, 
     main = "Clustering of module eigengenes",
     xlab = "", 
     sub = "")
# 在聚类图中画出剪切线
abline(h=0.1, col = "red")  # h=MEDissThres
# MEDissThres=0.25

# Plot the result
plotEigengeneNetworks(MEs, 
                      "Eigengene adjacency heatmap", 
                      marHeatmap = c(3,4,2,2), 
                      plotDendrograms = F, 
                      excludeGrey = F, # grey模块 需要
                      greyLabel = "grey",
                      xLabelsAngle = 90) 

# 如果要看看旧的模块这个分类情况，只需要把MEList改一下模块的颜色，重新构建一个METree即可
# MEList = moduleEigengenes(dataExpr, 
#                           colors = moduleColors)  # 旧模块颜色 moduleColors
# MEs = MEList$eigengenes
# MEDiss = 1-cor(MEs);
# METree = hclust(as.dist(MEDiss), method = "average");

# 画出热图
TOM = 1-TOMsimilarityFromExpr(dataExpr, power = 5);
plotTOM = TOM^7;
diag(plotTOM) = NA;
sizeGrWindow(9,9)
geneTree = net$dendrograms[[1]];
TOMplot(plotTOM, geneTree, 
        moduleColors,  # 新模块颜色mergedColors.my
        main = "Network heatmap plot, all genes")


TOMplot(plotTOM, geneTree, 
        mergedColors.my,  # 新模块颜色mergedColors.my
        main = "Network heatmap plot, all genes")


# 画出随机400个基因的热图 zanshi buyong
adjacency = adjacency(dataExpr, power = power)
TOM = TOMsimilarity(adjacency)  # 把邻接矩阵转换为拓扑重叠矩阵，以降低噪音和假相关，获得距离矩阵。
dissTOM = 1-TOM

nSelect = 400   # 根据nGenes的数值
# For reproducibility, we set the random seed 
set.seed(10); 
select = sample(nGenes, size = nSelect); 
selectTOM = dissTOM[select, select]; 
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster. 
selectTree = hclust(as.dist(selectTOM), method = "average") 
selectColors = mergedColors.my[select];  # moduleColors[select]
# Open a graphical window 
sizeGrWindow(9,9) 
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing 
# the color palette; setting the diagonal to NA also improves the clarity of the plot 
softPower <- sft$powerEstimate
plotDiss = selectTOM^softPower; 
diag(plotDiss) = NA; 
TOMplot(plotDiss, 
        selectTree, 
        selectColors, 
        main = "Network heatmap plot, selected genes")  # OK


# 导出网络到Cytoscape
# 1 旧模块
probes = colnames(dataExpr) 
dimnames(TOM) <- list(probes, probes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(TOM, 
                               edgeFile = paste(exprMat, ".edges.txt", sep=""), 
                               nodeFile = paste(exprMat, ".nodes.txt", sep=""), 
                               weighted = TRUE, 
                               threshold = 0.02, 
                               nodeNames = probes, 
                               nodeAttr = moduleColors)
# write.list(cyt, "cyt.txt")
# 使用expredata_node.txt建立模块基因文本Module_Gene.txt



# 加上新模块方法
# Select modules需要修改，选择需要导出的模块颜色
modules = unique(mergedColors.my);  # 旧模块颜色 moduleColors
# Select module probes选择模块探测
probes = colnames(dataExpr)
inModule = is.finite(match(mergedColors, modules));
modProbes = probes[inModule];
#modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("expredata_edges",".txt", sep=""),
                               nodeFile = paste("expredata_node",".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.01,  # 关键，越小输出越多的颜色模块
                               nodeNames = modProbes,
                               #altNodeNames = modGenes,
                               nodeAttr = mergedColors[inModule])
# 使用expredata_node.txt建立模块基因文本Module_Gene.txt



## 解决module gene 问题 WGCNA补充2 ####
#biocLite("ape")
library(ape)
par(mar=c(2,2,2,2))
# par(mar=c(0.1,0.1,0.1,0.1))
#calculate eigengenes
MEs = moduleEigengenes(dataExpr, colors = moduleColors, excludeGrey = FALSE)$eigengenes
#calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
#cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = 'average');
#plot the result with phytools package
par(mar=c(2,2,2,2))
# 解决图片大问题 figure margins too large ####
sizeGrWindow(12,9)
#pdf(file = "Plots/ModuleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot.phylo(as.phylo(METree),type = 'fan',show.tip.label = FALSE, main='')
tiplabels(frame = 'circle',col='black', text=rep('',length(unique(moduleColors))), bg = levels(as.factor(moduleColors)))
# moduleColors 都可以更换成新的模块颜色mergedColors.my
# save figure





### 建立模块基因文本Module_Gene.txt
## 剔除灰色（grey）未分类到模块的基因
expredata_node <- read.table("DEGsAllSampleMatrix.revised.txt.nodes.txt",sep = "\t", stringsAsFactors = F,header = T, fill=TRUE)
expredata_node = expredata_node[,-2]
colnames(expredata_node) = c("Symbol","type")
unique(expredata_node$type)

# 本次命名所有的mergedColors为mergedColors.my，另外，可以用以前的模块颜色moduleColors
Module = data.frame(data.frame(table(mergedColors.my))[order(data.frame(table(mergedColors.my))[,2],decreasing=T),])## 默认升序，decreasing=T时降序
Module <- data.frame(Module=unique(Module[,1]))

# library(dplyr)
# Module = filter(Module, Module !='grey')

Module_colour_Gene = data.frame(Symbol="0",type="0",module="0")
for (i in 1:nrow(Module)){
  data_1 = expredata_node[expredata_node$type == Module[i,1],]
  data_2 = data.frame(data_1,module=paste("m",i,sep=""))
  Module_colour_Gene = rbind(Module_colour_Gene,data_2)
}
# > data.frame(data_1,module=paste("m",i,sep=""))
# Error in data.frame(data_1, module = paste("m", i, sep = "")) : 
#   参数值意味着不同的行数: 0, 1
# 这里因为unique(expredata_node$type)的颜色数量和Module的颜色数量不一致，要去前面的cyt改参数threshold = 0.01

Module_colour_Gene = Module_colour_Gene[-1,]
write.table(Module_colour_Gene,"Module_colour_Gene.txt",row.names = F,quote = F,sep="\t", col.names = T)

Module_Gene = Module_colour_Gene[,c(3,1)]
write.table(Module_Gene,"Module_Gene.txt",row.names = F,quote = F,sep="\t", col.names = T)


# ## 寻找每个模块的关键基因
# Module_colour = unique(Module_colour_Gene[,2:3])
# # 本次命名所有的mergedColors为mergedColors.my，另外，可以用以前的模块颜色moduleColors，因为新旧的module没有变化
# HubGenes <- data.frame(colour=rownames(data.frame(HubGenes=chooseTopHubInEachModule(dataExpr,moduleColors))),HubGenes=chooseTopHubInEachModule(dataExpr,moduleColors))
# HubGenes = merge(HubGenes, Module_colour, by.x = "colour",by.y = "type",all=FALSE)
# write.table(HubGenes,"HubGenes.txt",row.names = F,quote = F,sep="\t", col.names = T)


### 建立Crosstalk文本Module_crosstalk.txt
Module <- data.frame(unique(Module_Gene[,1]))#对模块列去重
dir.create("Crosstalk_text")
for(i in 1:nrow(Module)){
  a = Module_Gene[Module_Gene$Module == Module[i,1],]
  b = data.frame((a[,2]))
  d = t(b)
  write.table(d,paste("Crosstalk_text/",Module[i,1], '.txt', sep=""),row.names = F,quote = F,sep="\t")#生成txt文件    
}





########    关联表型     ########
# traitData = data.frame(Sample=substr(traitData[,1],1,16),traitData)
# traitData = traitData[,c(1,3:6)]
# write.table(traitData,"traitData.txt",row.names = F,quote = F,sep="\t", col.names = T)


### 加入关联表型数据
sample1 <- read.table('../GSE20681/sample1.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40

trait <- sample1[, c(1,3)]
# colnames(trait)[7] <- "T" 
# # rownames(dataExpr) TCGA.E9.A1R2.01A
# trait$case_id <- gsub("-", ".", trait$case_id)
# rownames(trait) <- trait[,1]
write.table(trait,"trait_1.txt",row.names = F,quote = F,sep="\t", col.names = T)

# 根据trait_1.txt改一下文件格式，更名为trait.txt
trait <- "trait.txt"
# GSM	Control	Case
# GSM1620819	1	0
# GSM1620828	1	0

# 读入表型数据，不是必须的
if(trait != "") {
  traitData <- read.table(file=trait, sep='\t', header=T, row.names=1,
                          check.names=FALSE, comment='',quote="")
  sampleName = rownames(dataExpr)
  traitData = traitData[match(sampleName, rownames(traitData)), ]
}


### 模块与表型数据关联
MEList = moduleEigengenes(dataExpr, # 旧模块颜色colors = moduleColors
                          colors = mergedColors.my)  
MEs = MEList$eigengenes
if (corType=="spearman") {  # corType=="pearsoon"
  modTraitCor = cor(MEs, traitData, method ="spearman",use = "all.obs")
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
} else {
  modTraitCorP = bicorAndPvalue(MEs, traitData,alternative ="two.sided", robustY=robustY)
  modTraitCor = modTraitCorP$bicor
  modTraitP   = modTraitCorP$p
}
## Warning in bicor(x, y, use = use, ...): bicor: zero MAD in variable 'y'.
## Pearson correlation was used for individual columns with zero (or missing)
## MAD.




### 绘制模块之间相关性图
MEList = moduleEigengenes(dataExpr, # 旧模块颜色colors = moduleColors
                          colors = mergedColors.my)  
MEs = MEList$eigengenes
MEs_col <- MEs
# 根据基因间表达量进行聚类所得到的各模块间的相关性图
# marDendro/marHeatmap 设置下、左、上、右的边距

nrow(MEs_col)
length(MEs_col)==3
rownames(MEs_col) <- rownames(traitData)
MEs_col <- MEs_col[1:3, ]
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)            
# plot Module Eigengene all

## 如果有表型数据(traitData)，也可以跟ME数据放一起，一起出图 不太好看 
x <- cbind(MEs, traitData)
x <- na.omit(x)
MEs_colpheno = orderMEs(x)
# colnames(MEs_colpheno)[1,7] <- c("MEpink", "")
plotEigengeneNetworks(MEs_colpheno, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(7,7,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)                    
# 保存图片 Eigengene adjacency heatmap



# signif表示保留几位小数
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData), 
               yLabels = colnames(MEs_col), 
               cex.lab = 1, 
               ySymbols = colnames(MEs_col), colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 0.6, zlim = c(-1,1),
               main = paste("Module-trait relationships"))    
# 保存图片 Module-trait relationships
# modTraitCor
# MEblue  MEturquoise MEyellow   0.89   6e-4
# 结论，3个模块，和疾病相关性最好

## 从上图可以看到MEmagenta与Insulin_ug_l相关
## 模块内基因与表型数据关联
# 性状跟模块虽然求出了相关性，可以挑选最相关的那些模块来分析，
# 但是模块本身仍然包含非常多的基因，还需进一步的寻找最重要的基因。
# 所有的模块都可以跟基因算出相关系数，所有的连续型性状也可以跟基因的表达
# 值算出相关系数。
# 如果跟性状显著相关基因也跟某个模块显著相关，那么这些基因可能就非常重要。

### 计算模块与基因的相关性矩阵
if (corType=="pearsoon") {
  geneModuleMembership = as.data.frame(cor(dataExpr, MEs, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
    as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}


# 计算性状与基因的相关性矩阵
## 只有连续型性状才能进行计算，如果是离散变量，在构建样品表时就转为0-1矩阵。
if (corType=="pearsoon") {
  geneTraitCor = as.data.frame(cor(dataExpr, traitData, use = "p"))
  geneTraitP = as.data.frame(corPvalueStudent(
    as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(dataExpr, traitData, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}
## Warning in bicor(x, y, use = use, ...): bicor: zero MAD in variable 'y'.
## Pearson correlation was used for individual columns with zero (or missing)
## MAD.


# 最后把两个相关性矩阵联合起来,指定感兴趣模块进行分析
table(mergedColors)

# module = "red"
# pheno = "60 days"
module = "turquoise"
pheno = "Case"
modNames = substring(colnames(MEs), 3) # 将MEbrown变成brown
# 获取关注的列
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(traitData))
# 获取模块内的基因
moduleGenes = moduleColors == module

sizeGrWindow(7, 7)
par(mfrow = c(1,1))
# 与性状高度相关的基因，也是与性状相关的模型的关键基因
# 方案中11、关键模块hub基因筛选内容 ####
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   # abs(geneTraitCor[moduleGenes, pheno_column]),
                   abs(geneModuleMembership[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)        
# 保存图片    Module_gene significance_brown        




sizeGrWindow(7, 7)
par(mfrow = c(1,1))
# 与性状高度相关的基因，也是与性状相关的模型的关键基因
# 方案中11、关键模块hub基因筛选内容 ####
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   # abs(geneTraitCor[moduleGenes, pheno_column]),
                   abs(geneModuleMembership[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)        
# 保存图片    Module membership vs. gene significance        




# 3 PPI model ####
setwd("/Users/emperor/Documents/DataR/HeartD/3PPI/")
Module_colour_Gene <- read.table("Module_colour_Gene.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE)
# 要blue yellow turquoise青色这三个模块基因
geneBlue <- Module_colour_Gene[grep("blue", Module_colour_Gene$type), ]
geneYellow <- Module_colour_Gene[grep("yellow", Module_colour_Gene$type), ]
geneTurquoise <- Module_colour_Gene[grep("turquoise", Module_colour_Gene$type), ]
genePPI <- rbind(geneBlue, geneYellow, geneTurquoise)
genePPI <- unique(genePPI) # 315


geneInteg=genePPI
# 出一张PPI的图，和一张hubgene的图，仅仅只有我的基因的那种PPI
# 读取ppi模块基因
String_hsa_PPIs_Symbol <- read.table("/Users/emperor/Documents/DataR/zzhong.Zhaoyuanyuan/PPI/String_human_PPIs_Symbol.txt",sep = "\t",comment.char = "#", stringsAsFactors = F,header = T, fill=TRUE)
String_hsa_PPIs_Symbol = String_hsa_PPIs_Symbol[String_hsa_PPIs_Symbol$combined_score > 100,]

#读入网络基因，取消科学计数法，并进行筛选
diff_Gene <- unique(data.frame(Gene = geneInteg$Symbol))  # 85
diff_Gene$Gene <- toupper(diff_Gene$Gene)
Gene1=merge(diff_Gene, String_hsa_PPIs_Symbol, by.x = "Gene",by.y = "protein1",all=FALSE)
Gene2=merge(diff_Gene, String_hsa_PPIs_Symbol, by.x = "Gene",by.y = "protein2",all=FALSE)
colnames(Gene2)=colnames(Gene1)
Gene=rbind(Gene1,Gene2)
table(Gene$Gene %in% toupper(geneInteg$Symbol))  # 640346  

#unique是去掉基因名字重复的基因，只保留一个
Gene=unique(data.frame(Gene=Gene[,2]))
Gene=rbind(Gene,diff_Gene) # 
Gene=unique(Gene) # 15890
#write.table(Gene, 'Gene.txt',sep = '\t', quote=F, row.names=F)
write.table(Gene, 'Gene.txt',sep = '\t', quote=F, row.names=F)


data1 = merge(Gene,String_hsa_PPIs_Symbol,by.x = "Gene",by.y = "protein1",all=FALSE)
data2 = merge(Gene,data1,by.x = "Gene",by.y = "protein2",all=FALSE)
colnames(data2)=c("Gene1","Gene2","score")
data2$Gene1 <-as.character(data2[,1])
data2$Gene2 <-as.character(data2[,2])

if (nrow(data2) > 0 ){
  table(data2$Gene1>data2$Gene2)
  z=data2[data2$Gene1>data2$Gene2,]
  zz=data2[data2$Gene1<data2$Gene2,]
  zzz=zz[,c(2,1,3)]
  colnames(zzz)=colnames(z)
  zzzz=unique(rbind(z,zzz))
  table(zzzz$Gene1>zzzz$Gene2)
  network=zzzz
}
#write.table(network,"ppis_network_PI.txt",row.names = F,quote = F,sep="\t")  # 900
write.table(network,"ppis_network.txt",row.names = F,quote = F,sep="\t")


# 找一下两个数据集中 都高表达或者都低表达的基因
DEG_sort_p1 <- read.table("../GSE20681/DEG_sort_p1.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE)
DEG_sort_p2 <- read.table("../GSE42148/DEG_sort_p2.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE)

DEG_sort_pp <- merge(DEG_sort_p1, DEG_sort_p2, by = "DEG", all= F)
DEG_sort_pp <- DEG_sort_pp[DEG_sort_pp$logFC.x>0 &DEG_sort_pp$logFC.y>0 | DEG_sort_pp$logFC.x<0 &DEG_sort_pp$logFC.y<0, ]
write.table(DEG_sort_pp,"DEG_sort_pp.txt",row.names = F,quote = F,sep="\t")





# 保证network的gene1和gene2均是差异基因
DEgene <- unique(data.frame(symbol = DEG_sort_pp$DEG))
netDEG1 <- merge(DEgene, network, by.x = "symbol", by.y = "Gene1", all = F)  # 12853
netDEG2 <- merge(DEgene, netDEG1, by.x = "symbol", by.y = "Gene2", all = F)  # 1060
netDEG2 <- unique(netDEG2)
netDEG2 <- netDEG2[order(netDEG2$score, decreasing = T), ]
netDEG2 <- unique(netDEG2)
write.table(netDEG2,"network_DEG.txt",row.names = F,quote = F,sep="\t")


# 为ppis_network_DEG中的所有的gene加上差异的logFC
logfc.gene <- rbind(data.frame(gen = netDEG2$symbol), data.frame(gen = netDEG2$symbol.y))
logfc.gene <- unique(logfc.gene)  # 74
table <- data.frame(gen=geneInteg[toupper(geneInteg$Symbol) %in% logfc.gene$gen, ])
table$gen <- toupper(table$gen)
# 加上var
DEG_sort <-read.table('DEG_sort_pp.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
table2 <- merge(table, DEG_sort, by.x = "gen.Symbol", by.y = "DEG", all = F)

write.table(table2,"network_table.txt",row.names = F,quote = F,sep="\t")
write.table(logfc.gene,"network.gene.txt",row.names = F,quote = F,sep="\t")
# 选用network_DEG 和 network_table 这两个文件直接导入cytoscape中作图
# 这个直接导入就行，不需要进行clusterONE分析什么的


network_DEG <- read.table("./network_DEG.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE)
finaHub <- read.table("network_DEG_degree.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE)
finaHub <- finaHub[order(finaHub$degree.layout,decreasing = T), ]
finaHub <- finaHub[1:20, ]

network_hub1 <- network_DEG[network_DEG$symbol %in% finaHub$name, ]
network_hub2 <- network_hub1[network_hub1$symbol.y %in% finaHub$name, ]
write.table(network_hub2,"hubnetwork.txt",row.names = F,quote = F,sep="\t")
# 这个直接导入就行，可以得到一个小的hub网络
# hubgene
geneInteg2 <- geneInteg
geneInteg2$symbol <- toupper(geneInteg2$symbol)
hubgene <- merge(finaHub, geneInteg2, by.x = "name", by.y = "symbol", all = F)
write.table(hubgene,"hubgene20.txt",row.names = F,quote = F,sep="\t")











# 4 Enrichment ####
setwd("/Users/emperor/Documents/DataR/HeartD/4Enrichment//")
# setwd("/Users/emperor/Documents/DataR/HeartD/4Enrichment.revised/")

Module_colour_Gene <- read.table("Module_colour_Gene.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE)

Blue_select <- read.table("Blue.KEGG_select.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE)
Yellow_select <- read.table("Yellow_KEGG_select.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE)
Turquoise_select <- read.table("Turquoise_KEGG_select.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE)

# 气泡图
library(ggplot2)
Inputdata <- Turquoise_select
Inputdata <- Inputdata[, c(1, 4, 7)]
colnames(Inputdata) = c("Term","Count","pvalue") 
Inputdata$Count <- as.numeric(Inputdata$Count)
Inputdata$pvalue <- as.numeric(Inputdata$pvalue)
Inputdata <- Inputdata[order(Inputdata$Count, decreasing = T), ]
pp=ggplot(Inputdata[c(1:15), ], aes(Count,Term))
pbubble=pp+geom_point(aes(size=Count,color=-log(pvalue,10)))
pbubble+scale_colour_gradient(low="#0099FF",high="#CC00CC")
pr = pbubble + scale_size_continuous(range = c(2,7))+ scale_colour_gradient(low = "#0099FF", high = "#CC00CC")+ labs(color = expression(-log[10](FDR)),size="Count", x = "Count",y = "")
pr + xlim(0,8) 
# +theme(axis.text.x = element_text(size = 11, family = "myFont", color = "black", vjust = 0.5, hjust = 0.5), axis.text.y = element_text(size = 11, family = "myFont", color = "black"))


# 本地分析
# 转id
library("org.Hs.eg.db")
keytypes(org.Hs.eg.db)
blue <- Module_colour_Gene[grep("blue", Module_colour_Gene$type), ]
yellow <- Module_colour_Gene[grep("yellow", Module_colour_Gene$type), ]

ensids <- yellow$Symbol  # esembl ENSG00000223972
ensids <- as.character(ensids)
cols <- c("ENSEMBL","ENTREZID")
resu <- select(org.Hs.eg.db, keys=ensids,columns=cols,keytype="SYMBOL")
gene <- unique(resu$ENTREZID)  # KIF21B, 23046
gene <- na.omit(gene)

library(clusterProfiler)
KEGG <- enrichKEGG(gene = gene, 
                   #organism = "mmu", 
                   organism = "hsa",
                   keyType = "kegg",
                   pvalueCutoff = 0.05, pAdjustMethod = "holm", minGSSize = 10,
                   qvalueCutoff = 0.2, use_internal_data = F)
KEGG_Pathway <- as.data.frame(KEGG @ result)
write.table(KEGG_Pathway, "KEGG_Pathway.txt",row.names = F,quote = F,sep="\t")

# # KOBAS网站上进行kegg分析  http://kobas.cbi.pku.edu.cn/kobas3/genelist/

# GO分析
library(clusterProfiler)
library(org.Hs.eg.db)
ego_CC <- enrichGO(gene = gene,
                   #OrgDb=org.Mm.eg.db,
                   OrgDb=org.Hs.eg.db,
                   ont = "CC",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
CC = as.data.frame(ego_CC @ result)
write.table(CC, "Blue_GO_CC.txt",row.names = F,quote = F,sep="\t")

ego_BP <- enrichGO(gene = gene,
                   #OrgDb=org.Mm.eg.db,
                   OrgDb=org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
BP = as.data.frame(ego_BP @ result)
write.table(BP, "Blue_GO_BP.txt" ,row.names = F,quote = F,sep="\t")

ego_MF <- enrichGO(gene = gene,
                   #OrgDb=org.Mm.eg.db,
                   OrgDb=org.Hs.eg.db,
                   ont = "MF", # 可以选择all
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
MF = as.data.frame(ego_MF @ result)
write.table(MF, "Blue_GO_MF.txt" ,row.names = F,quote = F,sep="\t")

# yellow
ego_all <- enrichGO(gene = gene,
                   #OrgDb=org.Mm.eg.db,
                   OrgDb=org.Hs.eg.db,
                   ont = "All", # 可以选择all
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = TRUE)
All = as.data.frame(ego_MF @ result)
write.table(All, "Blue_GO_MF.txt" ,row.names = F,quote = F,sep="\t")


# 用kobas的bioinfo这个网站做的go分析
Blue_GO_select <- read.table("1Blue.GO.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE)
Turquoise_GO_select <- read.table("3Turquoise.GO.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE)
# Yellow_GO_select <- read.table("Yellow_GO_select.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE)

GO <- Turquoise_GO_select
library(ggplot2) # GO柱状图 金明
library(forcats)
color = colorRampPalette(colors = c("blue","white","red"))(45)
ggplot(GO)+ 
  geom_bar(aes(Ratio,fct_reorder(GO.Term,Database),fill = Database),stat="identity")+
  theme_bw()+
  labs(x='Number of genes',y='')+
  scale_fill_manual(values=c(BlueGO = "#0B0BFF")) # darkblue
 # https://www.modb.pro/db/598862



# 4 GSVA ####
## try http:// if https:// URLs are not supported
# source("https://bioconductor.org/biocLite.R")
# options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
# biocLite("GSVA")
library(GSVA)
browseVignettes("GSVA")
browseVignettes("estimate")

# test
library(GSVA)
p <- 20000    ## number of genes
n <- 30       ## number of samples
nGS <- 100    ## number of gene sets
min.sz <- 10  ## minimum gene set size
max.sz <- 100 ## maximum gene set size
X <- matrix(rnorm(p*n), nrow=p, dimnames=list(1:p, 1:n))
dim(X)
gs <- as.list(sample(min.sz:max.sz, size=nGS, replace=TRUE)) ## sample gene set sizes
gs <- lapply(gs, function(n, p) sample(1:p, size=n, replace=FALSE), p) ## sample gene sets
es.max <- gsva(X, gs, mx.diff=FALSE, verbose=FALSE, parallel.sz=1)
es.dif <- gsva(X, gs, mx.diff=TRUE, verbose=FALSE, parallel.sz=1)
pheatmap::pheatmap(es.max)
pheatmap::pheatmap(es.dif)

# OUR data
# 下载基因集网站 https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp
library(GSVA) # BiocManager::install('GSVA')
library(GSEABase)
setwd("/Users/emperor/Documents/DataR/HeartD/4GSVA/")
Exp1 <- read.table("../GSE20681/Exp1.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE)
rownames(Exp1) <- Exp1[, 1]
Exp1 <- Exp1[, -1]
Exp1[1:3,1:4]
X <- as.matrix(Exp1)
c3gsc2 <- getGmt("/Users/emperor/Documents/DataR/HeartD/4GSVA/c3.all.v2023.1.Hs.symbols.gmt",
                 collectionType=BroadCollection(category="c3"),
                 geneIdType=SymbolIdentifier())
# @ 解决！如果出现这个问题，对矩阵进行as.matrix,也就是X <- as.matrix(Exp1)
# Error in (function (classes, fdef, mtable)  : 
#             unable to find an inherited method for function ‘gsva’ for signature ‘"data.frame", "GeneSetCollection"’
#           
es.dif <- gsva(X, c3gsc2, method="zscore", kcdf="Gaussian",min.sz = 1,
               tau=switch(method, gsva=1, ssgsea=0.25, NA),
               mx.diff=TRUE, verbose=FALSE, parallel.sz=1)

es.dif[1:3,1:4]
gseaplot(es.dif[[2]],'KEGG_CELL_CYCLE') 


# 4.3 GSEA ####
# 教程: https://zhuanlan.zhihu.com/p/373388304
# 手把手教你用R做GSEA分析
setwd("/Users/emperor/Documents/DataR/HeartD/4Enrichment//")
setwd("/Users/emperor/Documents/DataR/HeartD/4Enrichment.revised/")

Module_colour_Gene <- read.table("Module_colour_Gene.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE)

library("org.Hs.eg.db")
keytypes(org.Hs.eg.db)
blue <- Module_colour_Gene[grep("blue", Module_colour_Gene$type), ]
# yellow <- Module_colour_Gene[grep("yellow", Module_colour_Gene$type), ]
turquoise <- Module_colour_Gene[grep("turquoise", Module_colour_Gene$type), ]

ensids <- turquoise$Symbol  # esembl ENSG00000223972
ensids <- as.character(ensids)
cols <- c("ENSEMBL","ENTREZID")
# @解决问题，如果出现报错"select"没有适用于"c('OrgDb', 'AnnotationDb', 'envRefClass', '.environment', 'refClass', 'environment', 'refObject', 'AssayData')"目标对象的方法
# 是因为dplyr函数也有select这个函数，r不知道用哪个，所以要指明AnnotationDbi::
# 参考 https://blog.csdn.net/abc1026497385/article/details/108607681
resu <- AnnotationDbi::select(org.Hs.eg.db, keys=ensids,columns=cols,keytype="SYMBOL")
resu <- resu[, c(1,3)]
# 输入文件准备，两列，一列是entriz id，一列是logFC值
DEG_sort_p1 <- read.table("../GSE20681/DEG_sort_p1.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE)
DEG_sort_p1 <- DEG_sort_p1[, 1:2]
InputGSEA <- merge(resu,DEG_sort_p1, by.x = "SYMBOL", by.y = "DEG", all = F )
InputGSEA <- na.omit(InputGSEA)
InputGSEA$logFC<-sort(InputGSEA$logFC,decreasing = T)
InputGSEA <- InputGSEA[, -1] # entriz id, logFC
geneList = InputGSEA[,2]
names(geneList) = as.character(InputGSEA[,1])
geneList
# 383      84002        597       7852       7850      10008      63971      57162      92162 
# 1.0862611  1.0727203  0.9688986  0.8579301  0.7504825  0.7177506  0.7139636  0.6018112  0.5599930 
# 23528      59348       9658 
# 0.5558811 -0.6105315 -0.7608811 

#GSEA分析——Reactome
library(fgsea)
library(clusterProfiler)
library(enrichplot)
library(ReactomePA)
library(tidyverse)
library(org.Hs.eg.db)
library(biomaRt)
Go_Reactomeresult <- gsePathway(geneList, nPerm = 1000, minGSSize = 1, maxGSSize = 1000, pAdjustMethod = "BH",pvalueCutoff = 1)
# blue lipid+iron+insulin
write.table (Go_Reactomeresult, file ="1Blue_GSEA_Reactomeresult1.txt", sep="\t", row.names =TRUE)
gseaplot2(Go_Reactomeresult, c(1,14,18,20,22,72,135), pvalue_table = F)
# # yellow N 
# write.table (Go_Reactomeresult, file ="Yellow_GSEA_Reactomeresult1.txt", sep="\t", row.names =TRUE)
# gseaplot2(Go_Reactomeresult, c(2,4,13,36,60,62,90,148), pvalue_table = F)
# turquoise immune/炎症
write.table (Go_Reactomeresult, file ="3Turquoise_GSEA_Reactomeresult1.txt", sep="\t", row.names =TRUE)
gseaplot2(Go_Reactomeresult, c(2,3,8,14,51,90, 94,138,207,301,357), pvalue_table = F)


# gseaplot(Go_Reactomeresult,1,pvalue_table = TRUE)
# KEGG_gseresult <- gseKEGG(geneList, nPerm = 1000, minGSSize = 1, maxGSSize = 1000, pvalueCutoff=1)
# Go_gseresult <- gseGO(geneList, 'org.Hs.eg.db', keyType = "ENTREZID", ont="all", nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
# write.table (Go_gseresult, file ="Blue_GSEA_goresult.txt", sep="\t", row.names =TRUE)
# #波浪图
# ridgeplot(Go_Reactomeresult,10) #输出前十个结果








# 5 Random SVM ####

## 5.1 Randomforest ####
setwd("/Users/emperor/Documents/DataR/HeartD/5Randon_SVM/")
setwd("/Users/emperor/Documents/DataR/HeartD/5Randon_SVM.revised//")
exp1 <- read.table('../GSE20681/Exp1.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40

# validation
exp3 <- read.table('../9Validation/Exp3.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40

Module_colour_Gene <- read.table('./Module_colour_Gene.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
Blue <- Module_colour_Gene[grep("blue", Module_colour_Gene$type), ]
brown <- Module_colour_Gene[grep("brown", Module_colour_Gene$type), ]
Turquoise <- Module_colour_Gene[grep("turquoise", Module_colour_Gene$type), ]


need1 <- Blue
need1 <- brown
need1 <- Turquoise

ranforData2 <- exp1[exp1$Symbol %in% need1$Symbol, ]

# validatioon
ranforData2 <- exp3[exp3$Symbol %in% need1$Symbol, ]
rownames(ranforData2) <- ranforData2[,1]
ranforData2 <- ranforData2[, -1]
ranforData2 <- t(ranforData2)
ranforData2 <- data.frame(ranforData2)
ranforData2 <- data.frame(aa=rownames(ranforData2), ranforData2)


library(randomForest)
ranforData2$aa <- as.factor(ranforData2$aa) # @ 第一列要factor格式才可以，不然会出现“参数不是数值也不是逻辑值：回覆NA”报错
need.rf2 <- randomForest(aa ~ ., data=ranforData2, importance=TRUE, proximity=TRUE)
print(need.rf2)

need.mds2 <- cmdscale(1 - need.rf2$proximity, eig=TRUE)
plot(need.mds2$points, col = rep(c("red", "blue", "green"), each = 50))
need.rf2$importance
varImpPlot(need.rf2, main = "Top 15 - variable importance") # OK
resu2 <- data.frame(need.rf2$importance)
resu2 <- resu2[, 196:197]
resu2 <- data.frame(symbol = rownames(resu2), resu2)

write.table(resu2,"1Blue_RandomForest_Importance.txt",row.names = F,quote = F,sep="\t")
write.table(resu2,"2brown_RandomForest_Importance.txt",row.names = F,quote = F,sep="\t")
write.table(resu2,"3Turquoise_RandomForest_Importance.txt",row.names = F,quote = F,sep="\t")

# validation
write.table(resu2,"/Users/emperor/Documents/DataR/HeartD/9Validation/Validation_1Blue_RandomForest_Importance.txt",row.names = F,quote = F,sep="\t")


# 5.1.1 ROC-random forest ####
# 教程 https://www.bioinfo.online/articleList/20228952669.html
library(pROC)

sample1 <- read.table('../GSE20681/sample1.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
sample1 <- sample1[, c(1,3)]
ranforData3 <- merge(sample1, ranforData2, by.x = "GSM", by.y= "aa", all = F)
# ranforData3 是补充了sample的case和control分类，ranforData2才是构建随机森林的内容

sample3 <- read.table('../9Validation/sample3.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
sample3 <- sample3[, c(1,3)]
ranforData3 <- merge(sample3, ranforData2, by.x = "GSM", by.y= "aa", all = F)


library(RColorBrewer)
mycol<-brewer.pal(11, "Spectral")

plot.roc(ranforData3$Group, 
         print.auc=T,
         col="#00FFFF", # 颜色 mycol[9]
         predict(need.rf2, newdata = ranforData2,type="vote")[,2])

legend("bottomright", legend=c("TurquoiseRForest, AUC=0.563"),
       col="#00FFFF",lwd=2,cex = 0.8,
       x.intersp = 2, y.intersp=0.8, adj = c(0.1, 0.4))




# 5.1.2 cloudwind ####
RaF <- resu2
data2 <- RaF

# 环形柱状图 金明
# 教程 https://www.jianshu.com/p/865aa9023b27
df<-data.frame(individual=paste("feature",seq(1,24),sep=""),value=c(24:1))
df$id<-seq(1,nrow(df))
library(ggplot2)
p<-ggplot(df,aes(x=as.factor(id),y=value,fill=as.factor(id)))+
  geom_bar(stat="identity")
p+coord_polar()+theme_bw()+ theme(legend.position="none")+ylim(-2,24)

# our result 
RaF1 <- read.table('1Blue_RandomForest_Importance.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
RaF1 <- read.table('2brown_RandomForest_Importance.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
RaF1 <- read.table('3Turquoise_RandomForest_Importance.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40

RaF1 <- read.table('../9Validation/Validation_1Blue_RandomForest_Importance.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40


RaF1 <- RaF1[order(RaF1$MeanDecreaseGini, decreasing = T), ]
RaF1$MeanDecreaseGini <- RaF1$MeanDecreaseGini*10
# RaF1 <- as.factor(RaF1$symbol, levels=RaF1$symbol)
RaF1$id<-seq(1,nrow(RaF1))

p<-ggplot(RaF1[1:10, ], aes(x=as.factor(id), y=MeanDecreaseGini, fill=as.factor(symbol)))+
  geom_bar(stat="identity")
p+coord_polar()+theme_bw()+ theme(legend.position="none")
  # +ylim(-2,4)
write.table(RaF1,"1Blue_RandomForest_Importance.legend.txt",row.names = F,quote = F,sep="\t")
write.table(RaF1,"2brown_RandomForest_Importance.legend.txt",row.names = F,quote = F,sep="\t")
write.table(RaF1,"3Turquoise_RandomForest_Importance.legend.txt",row.names = F,quote = F,sep="\t")

write.table(RaF1,"/Users/emperor/Documents/DataR/HeartD/9Validation/Validation_1Blue_RandomForest_Importance.legend.txt",row.names = F,quote = F,sep="\t")


# 5.2 svm 三组 #### 
setwd("/Users/emperor/Documents/DataR/HeartD/5Randon_SVM/")
setwd("/Users/emperor/Documents/DataR/HeartD/5Randon_SVM.revised/")

exp1 <- read.table('../GSE20681/Exp1.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
exp1 <- read.table('../GSE20681/Exp1.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40

Module_colour_Gene <- read.table('./Module_colour_Gene.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
Blue <- Module_colour_Gene[grep("blue", Module_colour_Gene$type), ]
brown <- Module_colour_Gene[grep("brown", Module_colour_Gene$type), ]
Turquoise <- Module_colour_Gene[grep("turquoise", Module_colour_Gene$type), ]


need1 <- Blue
need2 <- brown
need3 <- Turquoise

svmData1 <- exp1[exp1$Symbol %in% need2$Symbol, ]
rownames(svmData1) <- svmData1[, 1]
svmData1 <- svmData1[, -1]

svmData1 <- t(svmData1)
#              ABHD5    ACOX1    ACSL4     ANXA3
# GSM518885 11.070190 8.445527 8.613062 10.373875
# GSM518886 10.012830 7.771118 7.860759  8.707333

# 分成训练集和测试集
smp.size <- floor(0.6*nrow(svmData1))
train.ind <- sample(seq_len(nrow(svmData1)), smp.size)
train <- svmData1[train.ind, ]
train1 <- data.frame(id=rownames(train), train) # 138 hang

test <- svmData1[-train.ind, ]
test1 <- data.frame(id=rownames(test), test) # 60 hang

#挤入临床信息 二分类
sample1 <- read.table('../GSE20681/sample1.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
sample1 <- sample1[, c(1,3)]
train2 <- merge(sample1, train1, by.x = "GSM", by.y = "id", all = F)
rownames(train2) <- train2[, 1]
train2 <- train2[, -1] # 去掉gsm
test2 <- merge(sample1, test1, by.x = "GSM", by.y = "id", all = F)
rownames(test2) <- test2[, 1]
test2 <- test2[, -1] # 去掉gsm


library(e1071)
train2$Group <- as.factor(train2$Group) # 如果报错Need numeric dependent variable for regression.，解决办法是加上as.factor
model <- svm(formula = Group ~ ., data = train2, kernel = "linear")
summary(model)



# 训练集的混淆矩阵
train.pred <- predict(model, train2)
table(real=train2$Group, predict=train.pred)
# blue      predict
# real      case control
# case      58       3
# control    0      57

# brown     predict
# real      case control
# case      46      11
# control   14      47

# Turquoise predict
# real      case control
# case      59       0
# control    0      59
confus.matrix <- table(real=train2$Group, predict=train.pred)
sum(diag(confus.matrix))/sum(confus.matrix) 
# blue 0.9745763  brown 0.7881356  turquoise 1.0 

# 测试集预测
# test2$Group <- as.factor(test2$Group)
test.pred <- predict(model, test2) # 一定要只有数据，不能有GSM

table(real=test2$Group, predict=test.pred)
confus.matrix2 <- table(real=test2$Group, predict=test.pred)
sum(diag(confus.matrix2))/sum(confus.matrix2) # 0.5125   0.475  0.5875


# 5.2.1 ROC ####
# 教程：https://www.jianshu.com/p/c926c191361e
# 提取模型预测值并进行格式处理
pred_1 <- as.factor(model$decision.values)
pred_1 <- as.ordered(pred_1)
modelroc_1 <- roc(train2$Group, pred_1)
modelroc_1
plot(modelroc_1, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2), grid.col=c("green", "red"), max.auc.polygon=TRUE, auc.polygon.col="skyblue", print.thres=TRUE)
# 可视化展示,使用add=TRUE将第二个模型添加到图形中
# plot.roc(modelroc_2, add=TRUE, col="green",print.thres=TRUE) 


# 5.2.2 SVM value ####
# 把基因和svm值组合起来
svmValue <- model[["x.scale"]]$`scaled:center`
svmValue <- data.frame(value=svmValue)
svmValue <- data.frame(symbol=rownames(svmValue), svmValue)

write.table(svmValue,"1Blue_SVM_Value.txt",row.names = F,quote = F,sep="\t")
write.table(svmValue,"2brown_SVM_Value.txt",row.names = F,quote = F,sep="\t")
write.table(svmValue,"3Turquoise_SVM_Value.txt",row.names = F,quote = F,sep="\t")


# 5.3 LASSO #### 
# 教程:https://www.jianshu.com/p/2995ae5f244e
# 教程更好： https://www.weinformatics.cn/44483852be/
setwd("/Users/emperor/Documents/DataR/HeartD/5Randon_SVM/")
setwd("/Users/emperor/Documents/DataR/HeartD/5Randon_SVM.revised/")

exp1 <- read.table('../GSE20681/Exp1.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
Module_colour_Gene <- read.table('./Module_colour_Gene.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
Blue <- Module_colour_Gene[grep("blue", Module_colour_Gene$type), ]
brown <- Module_colour_Gene[grep("brown", Module_colour_Gene$type), ]
Turquoise <- Module_colour_Gene[grep("turquoise", Module_colour_Gene$type), ]

need1 <- Blue
need3 <- Turquoise

LassoData1 <- exp1[exp1$Symbol %in% need1$Symbol, ]
rownames(LassoData1) <- LassoData1[, 1]
LassoData1 <- LassoData1[, -1]

sample1 <- read.table('../GSE20681/sample1.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
sample1 <- sample1[, c(1,3)]
meta <- sample1
meta$Group <- gsub("case", "1", meta$Group)
meta$Group <- gsub("control", "0", meta$Group)
meta$Group <- as.numeric(meta$Group)

# lasso
# cox1=c('RP4-616B8.5','RP11-389G6.3',"AP000696.2",'CTD-2377D24.6',"LINC01559" ,"LINC00629","AC005062.2","LINC01018")
# ##cox1为定义的任何向量，比如里面含有你所感兴趣的基因名称
# exprSet = exprSet[cox1,]#表达谱数据中提取你所感兴趣的基因表达
exprSet=as.data.frame(LassoData1)
dim(exprSet)
x=t(exprSet) # 样本在左rownames，基因在列，列名，colnames
y=meta$Group#提取病例的预后情况（0与1）列，meta是临床数据
library(glmnet)
model_lasso <- glmnet(x, y,nlambda=10, alpha=1)
print(model_lasso)
set.seed(13098)
cv_fit <- cv.glmnet(x=x, y=y, nlambda = 1000,alpha = 1)
plot(cv_fit)

# cv_fit$lambda.min blue
# [1] 0.008128992, log(lambda): -4.812318

# cv_fit$lambda.min yellow
# [1] 0.004071331

# cv_fit$lambda.min turquoise
# [1] 0.01013239 log(lambda): -4.592018

lasso.prob <- predict(cv_fit, newx=x , s=c(cv_fit$lambda.min,cv_fit$lambda.1se) )
head(lasso.prob)
re=cbind(y ,lasso.prob)
head(re)
re1=as.data.frame(re)
# write.csv(re,file='')
write.table(re1,"1Blue_LassoValue_sample.txt",row.names = F,quote = F,sep="\t")
# write.table(re1,"2LassoValue_Yellow_sample.txt",row.names = F,quote = F,sep="\t")
write.table(re1,"3Turquoise_LassoValue_sample.txt",row.names = F,quote = F,sep="\t")

# gene value ####
# 要基因对应的值  
# 两条虚线分别指示了两个特殊的λ值,一个是lambda.min,一个是lambda.1se,
# 这两个值之间的lambda都认为是合适的。lambda.1se构建的模型最简单，即使用的基因数量少，
# 而lambda.min则准确率更高一点，使用的基因数量更多一点。
# 用这两个λ值重新建模或者用min这个重新建模

model_lasso_min <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.min)
# model_lasso_1se <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.1se)
head(model_lasso_min$beta)
choose_gene_min=rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]
# choose_gene_1se=rownames(model_lasso_1se$beta)[as.numeric(model_lasso_1se$beta)!=0]

lassoValue <- model_lasso_min[["beta"]]@x
lassoValue <- data.frame(choose_gene_min,lassoValue)
write.table(lassoValue,"1Blue_LassoValue.txt",row.names = F,quote = F,sep="\t")
# write.table(lassoValue,"2LassoValue_Yellow.txt",row.names = F,quote = F,sep="\t")
write.table(lassoValue,"3Turquoise_LassoValue.txt",row.names = F,quote = F,sep="\t")



# 箱线图 ####
re=as.data.frame(re)
colnames(re)=c('event','prob_min','prob_1se')
re$event=as.factor(re$event)
library(ggpubr) 
p1 = ggboxplot(re, x = "event", y = "prob_min",
               color = "event", palette = "jco",
               add = "jitter")+ stat_compare_means()
p2 = ggboxplot(re, x = "event", y = "prob_1se",
               color = "event", palette = "jco",
               add = "jitter")+ stat_compare_means()
library(patchwork)
p1+p2
p1


# ROC-Lasso ####
library(ROCR)
library(caret)
# 自己预测自己
#min
pred_min <- prediction(re[,2], re[,1])
auc_min = performance(pred_min,"auc")@y.values[[1]]
perf_min <- performance(pred_min,"tpr","fpr")
plot(perf_min,colorize=FALSE, col="#00FFFF") 
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
text(0.8,0.2, labels = paste0("AUC = ",round(auc_min,3)))

#1se
pred_1se <- prediction(re[,3], re[,1])
auc_1se = performance(pred_1se,"auc")@y.values[[1]]
perf_1se <- performance(pred_1se,"tpr","fpr")
plot(perf_1se,colorize=FALSE, col="red") 
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
text(0.8,0.2, labels = paste0("AUC = ",round(auc_min,3)))







# 5.5 veen ####
setwd("/Users/emperor/Documents/DataR/HeartD/5Randon_SVM/Veen/")
setwd("/Users/emperor/Documents/DataR/HeartD/5Randon_SVM.revised/4veen/")

BlueRForest <- read.table('./1Blue_Random_median.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
BlueRForest.h <- BlueRForest[grep("highValue", BlueRForest$median.1.867151614), ]
colnames(BlueRForest.h)[1] <- "symbol"
BlueSVM <- read.table('./1Blue_SVM_median.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
BlueSVM.h <- BlueSVM[grep("highValue", BlueSVM$median.10.33748013), ]
colnames(BlueLasso)[1] <- "symbol"
BlueLasso <- read.table('./1Blue_LassoValue.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
colnames(BlueLasso)[1] <- "symbol"

Blue.h1 <- merge(BlueRForest.h, BlueSVM.h, by = "symbol", all = F)
Blue.h2 <- merge(Blue.h1, BlueLasso, by.x = "symbol", by.y = "symbol", all = F)
write.table(Blue.h2,"1Blue_mergeVeen10.txt",row.names = F,quote = F,sep="\t")


BlueRF.gene <-unique(BlueRForest.h$symbol)  # 51
BlueSVM.gene <-unique(BlueSVM.h$symbol)  # 51
BlueLasso.gene <-unique(BlueLasso$symbol)  # 48

library(VennDiagram)
venn.diagram(list(BlueRF=BlueRF.gene, BlueSVM=BlueSVM.gene, BlueLasso=BlueLasso.gene),
             resolution = 300, imagetype = "tiff", alpha=c(0.5,0.5,0.5),
             fill=c("darkred","darkgreen","darkblue"), 
             cat.cex=2,  # 标签字的大小
             # cat.fontfamily = 2, # 标签字字体
             cat.fontface =4,
             
             cex = 5, # 圈圈中数字的大小
             fontfamily="sans", # 数字的字体 "sans"
             units = "px", compression ="lzw", na = "stop", 
             main = NULL, sub = NULL, main.pos=c(0.5, 1.05), main.fontface = "plain",
             main.fontfamily = "serif", main.col = "black",
             main.cex = 1, main.just = c(0.5, 1),
             filename = "Venn.tiff")


library(grDevices)
venn <- venn.diagram(list(BlueRF=BlueRF.gene, BlueSVM=BlueSVM.gene, BlueLasso=BlueLasso.gene),
                     resolution = 300, 
                     # imagetype = "tiff", 
                     alpha=c(0.5,0.5,0.5),
                     fill=c("darkred","darkgreen","darkblue"), 
                     cat.cex=2,  # 标签字的大小
                     cat.fontfamily = 2, # 标签字字体
                     cat.fontface =4,
                     
                     cex = 5, # 圈圈中数字的大小
                     fontfamily=2, # 数字的字体 "sans"
                     units = "px", compression ="lzw", na = "stop", 
                     main = NULL, sub = NULL, main.pos=c(0.5, 1.05), main.fontface = "plain",
                     main.fontfamily = "serif", main.col = "black",
                     main.cex = 1, main.just = c(0.5, 1),
                     filename = NULL
)

pdf(file="venn_null.pdf")
grid.draw(venn)
dev.off()




TurquoiseRForest <- read.table('./3Turquoise_Random_median.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
TurquoiseRForest.h <- TurquoiseRForest[grep("highValue", TurquoiseRForest$median1.7210663), ]
colnames(TurquoiseRForest.h)[1] <- "symbol"

TurquoiseSVM <- read.table('./3Turquoise_SVM_median.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
TurquoiseSVM.h <- TurquoiseSVM[grep("highValue", TurquoiseSVM$median9.12110778), ]
colnames(TurquoiseSVM.h)[1] <- "symbol"

TurquoiseLasso <- read.table('./3Turquoise_Lasso_all.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
colnames(TurquoiseLasso)[1] <- "symbol"

Turquoise.h1 <- merge(TurquoiseRForest.h, TurquoiseSVM.h, by = "symbol", all = F)
Turquoise.h2 <- merge(Turquoise.h1, TurquoiseLasso, by = "symbol", all = F)
write.table(Turquoise.h2,"3Turquoise_mergeVeen8.txt",row.names = F,quote = F,sep="\t")

TurquoiseRF.gene <-unique(TurquoiseRForest.h$symbol)  # 57
TurquoiseSVM.gene <-unique(TurquoiseSVM.h$symbol)  # 57
TurquoiseLasso.gene <-unique(TurquoiseLasso$symbol)  # 42

library(VennDiagram)
venn.diagram(list(TurquoiseRF=TurquoiseRF.gene, TurquoiseSVM=TurquoiseSVM.gene, TurquoiseLasso=TurquoiseLasso.gene),
             resolution = 300, imagetype = "tiff", alpha=c(0.5,0.5,0.5),
             fill=c("darkred","darkgreen","darkblue"), 
             cat.cex=2,  # 标签字的大小
             # cat.fontfamily = 2, # 标签字字体
             cat.fontface =4,
             
             cex = 5, # 圈圈中数字的大小
             fontfamily="sans", # 数字的字体 "sans"
             units = "px", compression ="lzw", na = "stop", 
             main = NULL, sub = NULL, main.pos=c(0.5, 1.05), main.fontface = "plain",
             main.fontfamily = "serif", main.col = "black",
             main.cex = 1, main.just = c(0.5, 1),
             filename = "Venn.tiff")


library(grDevices)
venn <- venn.diagram(list(TurquoiseRF=TurquoiseRF.gene, TurquoiseSVM=TurquoiseSVM.gene, TurquoiseLasso=TurquoiseLasso.gene),
                     resolution = 300, 
                     # imagetype = "tiff", 
                     alpha=c(0.5,0.5,0.5),
                     fill=c("darkred","darkgreen","darkblue"), 
                     cat.cex=2,  # 标签字的大小
                     cat.fontfamily = 2, # 标签字字体
                     cat.fontface =4,
                     
                     cex = 5, # 圈圈中数字的大小
                     fontfamily=2, # 数字的字体 "sans"
                     units = "px", compression ="lzw", na = "stop", 
                     main = NULL, sub = NULL, main.pos=c(0.5, 1.05), main.fontface = "plain",
                     main.fontfamily = "serif", main.col = "black",
                     main.cex = 1, main.just = c(0.5, 1),
                     filename = NULL
)

pdf(file="venn_null.pdf")
grid.draw(venn)
dev.off()


# 
a <- read.table('./1Blue_mergeVeen10.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
a <- a[, c(1,2,4,6)]
colnames(a) <- c("gene", "MeanDecreaseGini", "SVMvalue", "lassoValue")
# b <- read.table('./Yellow_mergeVeen.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
# b <- b[, c(1,2,4,6)]
# colnames(b) <- c("gene", "MeanDecreaseGini", "SVMvalue", "lassoValue")

c <- read.table('./3Turquoise_mergeVeen8.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
c <- c[, c(1,2,4,6)]
colnames(c) <- c("gene", "MeanDecreaseGini", "SVMvalue", "lassoValue")

all <- rbind(a,  c)
write.table(all,"0all_mergeVeen.txt",row.names = F,quote = F,sep="\t")



YellowRForest <- read.table('./YellowRForest.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
YellowRForest.h <- YellowRForest[grep("highValue", YellowRForest$median.3.393555), ]
YellowS <- read.table('./YellowSVM.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
YellowS.h <- YellowS[grep("highValue", YellowS$median.9.122182), ]
YellowLasso <- read.table('./2LassoValue_Yellow.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
colnames(YellowLasso)[1] <- "symbol"

Yellow.h1 <- merge(YellowRForest.h, YellowS.h, by = "symbol", all = F)
Yellow.h2 <- merge(Yellow.h1, YellowLasso, by = "symbol", all = F)
write.table(Yellow.h2,"Yellow_mergeVeen.txt",row.names = F,quote = F,sep="\t")

YellowRF.gene <-unique(YellowRForest.h$symbol)  # 29
YellowSVM.gene <-unique(YellowS.h$symbol)  # 29
YellowLasso.gene <-unique(YellowLasso$symbol)  # 39

library(VennDiagram)
venn.diagram(list(YellowRF=YellowRF.gene, YellowSVM=YellowSVM.gene, YellowLasso=YellowLasso.gene),
             resolution = 300, imagetype = "tiff", alpha=c(0.5,0.5,0.5),
             fill=c("darkred","darkgreen","darkblue"), 
             cat.cex=2,  # 标签字的大小
             # cat.fontfamily = 2, # 标签字字体
             cat.fontface =4,
             
             cex = 5, # 圈圈中数字的大小
             fontfamily="sans", # 数字的字体 "sans"
             units = "px", compression ="lzw", na = "stop", 
             main = NULL, sub = NULL, main.pos=c(0.5, 1.05), main.fontface = "plain",
             main.fontfamily = "serif", main.col = "black",
             main.cex = 1, main.just = c(0.5, 1),
             filename = "Venn.tiff")


library(grDevices)
venn <- venn.diagram(list(YellowRF=YellowRF.gene, YellowSVM=YellowSVM.gene, YellowLasso=YellowLasso.gene),
                     resolution = 300, 
                     # imagetype = "tiff", 
                     alpha=c(0.5,0.5,0.5),
                     fill=c("darkred","darkgreen","darkblue"), 
                     cat.cex=2,  # 标签字的大小
                     cat.fontfamily = 2, # 标签字字体
                     cat.fontface =4,
                     
                     cex = 5, # 圈圈中数字的大小
                     fontfamily=2, # 数字的字体 "sans"
                     units = "px", compression ="lzw", na = "stop", 
                     main = NULL, sub = NULL, main.pos=c(0.5, 1.05), main.fontface = "plain",
                     main.fontfamily = "serif", main.col = "black",
                     main.cex = 1, main.just = c(0.5, 1),
                     filename = NULL
)

pdf(file="venn_null.pdf")
grid.draw(venn)
dev.off()





# 6 xcell ####
setwd("/Users/emperor/Documents/DataR/HeartD/ImmuneInfiltration.revised/")
# devtools::install_github('dviraran/xCell')
# survivalPDD 1581基因量不够
file.name<-'/Users/emperor/Documents/DataR/HeartD/GSE42148//Exp2.txt'
# file.name<-'/Users/emperor/Documents/DataR/UVM/1.consensis/EXP_UVM_PDD.txt'

expr = read.table(file.name, header=TRUE, as.is=TRUE, sep='\t')
rownames(expr) <- expr[, 1] # 用34949那个
expr <- expr[ , -1]
# expr <- expr[2000: 34000, ]

library(xCell)
# expr <- dplyr::distinct(expr, SYMBOL, .keep_all = T)

xCell.data$genes
xCell.data$signatures
scores01 = rawEnrichmentAnalysis(expr, xCell.data$signatures,
                                 xCell.data$genes, file.name = NULL, parallel.sz = 6, parallel.type = "Fork")
scores03 <- scores01/1000
write.table(scores03,"xcell_result_HeartCAD.revised.txt",row.names = F,quote = F,sep="\t")

# scores02 = xCellAnalysis(exp, rnaseq=F,cell.types.use = NULL)

scores04 <- t(scores03)
scores04 <- data.frame(scores04)
scores04 <- data.frame(id=rownames(scores04), scores04)
# scores04$id <- substr(scores04$id, 1, 15)

# 加入risk 
group <- read.table("../GSE42148/sample2.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE,check.names=F)


xCell.group <- merge(group, scores04, by.x = "GSM", by.y = "id", all = F)
write.table(xCell.group,"xCell.group_HeartCAD2.revised.txt",row.names = F,quote = F,sep="\t")





#correlation ####
setwd("/Users/emperor/Documents/DataR/HeartD/ImmuneInfiltration/xCell")
setwd("/Users/emperor/Documents/DataR/HeartD/ImmuneInfiltration.revised/")

Exp2 <- read.table("./Exp2.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE,check.names=F)

Blueg <- read.table("./1Blue_mergeVeen10-2LOC.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE,check.names=F)
Exp_blue <- Exp2[Exp2$Symbol %in% Blueg$symbol, ] # 11
rownames(Exp_blue) <- Exp_blue[, 1]
Exp_blue <- Exp_blue[, -1]
Exp_blue2 <- data.frame(t(Exp_blue))
Exp_blue3 <- data.frame(GSM = rownames(Exp_blue2), Exp_blue2)

Yellowg <- read.table("./Yellow_mergeVeen.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE,check.names=F)
Exp_Yellow <- Exp2[Exp2$Symbol %in% Yellowg$symbol, ] # 8
rownames(Exp_Yellow) <- Exp_Yellow[, 1]
Exp_Yellow <- Exp_Yellow[, -1]
Exp_Yellow2 <- data.frame(t(Exp_Yellow))
Exp_Yellow3 <- data.frame(GSM = rownames(Exp_Yellow2), Exp_Yellow2)

Turquoiseg <- read.table("./3Turquoise_mergeVeen8.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE,check.names=F)
Exp_Turquoise <- Exp2[Exp2$Symbol %in% Turquoiseg$symbol, ] # 9
rownames(Exp_Turquoise) <- Exp_Turquoise[, 1]
Exp_Turquoise <- Exp_Turquoise[, -1]
Exp_Turquoise2 <- data.frame(t(Exp_Turquoise))
Exp_Turquoise3 <- data.frame(GSM = rownames(Exp_Turquoise2), Exp_Turquoise2)

xCell.group <- read.table("./xCell_Select.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE,check.names=F)

corDATA.blue <- cbind(Exp_blue3, xCell.group)
corDATA.blue <- corDATA.blue[, -c(11,12)]


corDATA.yellow <- cbind(Exp_Yellow3, xCell.group)
corDATA.yellow <- corDATA.yellow[, c(1,10,11,12, 2:9, 13:76)]
corDATA.yellow <- corDATA.yellow[, -c(1,2,3)]


corDATA.turquoise <- cbind(Exp_Turquoise3, xCell.group)
corDATA.turquoise <- corDATA.turquoise[, -c(10,11)]




library(dplyr)
library(psych)
a <- corDATA.blue
# a <- corDATA.yellow
a <- corDATA.turquoise

# a <- data.frame( a[, c(1, 2:10)], a$CD8..Tem, a$Eosinophils, a$Macrophages.M2, a$NK.cells,
#            a$pDC, a$Tgd.cells, a$Th1.cells) # id riskscore cells
# rownames(a) <- a[,1]
a <- a[, -1]
for(i in 1:dim(a)[2]){
  a[,i]<-as.numeric(round(as.numeric(a[,i]),2))
}  # numeric化&保留小数点后2位
# corr.test:查找矩阵或data.frame元素之间的相关性，样本大小和概率值。
cor <- corr.test(a, use = "pairwise", method = "pearson", adjust = "fdr")
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat) 
  data.frame(row = rownames(cormat)[row(cormat)[ut]],column = rownames(cormat)[col(cormat)[ut]], cor =(cormat)[ut], p = pmat[ut] )
}
cor_text=flattenCorrMatrix(cor$r, cor$p)
write.table(cor_text,"3Turquoise_cor_text.txt",row.names = F,quote = F,sep="\t")


# 先分析一遍看看显著性
# 然后把显著性的基因和cell选出来，展示只展示这些显著的
# blue模块，
bbb <- corDATA.blue
gene.cor <- c("CTDSP2", "NLRP1", "RILPL2", "JUNB", "DHRS7")
cell.cor <- c("Macrophages.M2", "Smooth.muscle", "Th2.cells", "Tgd.cells","NK.cells")
bbb2 <- bbb[, colnames(bbb) %in%gene.cor]
bbb3 <- bbb[, colnames(bbb) %in%cell.cor]
bbb4 <- cbind(bbb2,bbb3)
b <- bbb4
for(i in 1:dim(b)[2]){
  b[,i]<-as.numeric(round(as.numeric(b[,i]),2))
}  # numeric化&保留小数点后2位
tdc <- cor (b, method="pearson")
corrplot(tdc, method = "ellipse", type = "upper",
         tl.col = "black", tl.cex = 0.6, tl.srt = 45,tl.pos = "lt")
corrplot(tdc, method = "number", type = "lower",
         tl.col = "n", tl.cex = 0.1, tl.pos = "n",number.cex = 0.6, # 字体
         add = T)

# 画图，按照对角线看对应关系
# install.packages("PerformanceAnalytics")
library(PerformanceAnalytics)
chart.Correlation(a, histogram=TRUE, pch=19)



library(corrplot)
tdc <- cor (a, method="pearson")
tdc1 <- tdc[1:9, ] # 取gene
tdc2 <- tdc1[, 10:17 ] # 取gene
tdc3 <- t(tdc2)
corrplot(tdc, method = "ellipse", type = "upper",
         tl.col = "black", tl.cex = 0.6, tl.srt = 45,tl.pos = "lt")
corrplot(tdc, method = "number", type = "lower",
         tl.col = "n", tl.cex = 0.1, tl.pos = "n",number.cex = 0.6, # 字体
         add = T)




#  gene 与 细胞相关性
# blue: PELI1	a.Th1.cells; CXCR4	a.pDC; CLEC4D	a.Eosinophils
# ANXA11	a.Macrophages.M2; ARHGAP25	a.Eosinophils
# turquoise JUNB	a.Tgd.cells; STK17B	a.pDC; MTMR3	a.Eosinophils; FAM113B	a.Eosinophils; TBC1D14	a.pDC; ACSL1	a.Th1.cells
corr_exp <- corDATA.blue
# corr_exp <- corDATA.yellow
corr_exp <- corDATA.turquoise
library(corrplot)
library(ggpubr)
colnames(corr_exp)
gene <- "MARCKS"
U1 <- ggscatter(corr_exp,x = gene, #x变量
                y = "Smooth.muscle" ,#y变量
                add = "reg.line",##拟合曲线
                #conf.int = TRUE,##置信区间阴影带
                cor.coef = TRUE, ##系数
                color = "black",
                add.params = list(color = "blue",fill = "lightgray"),
                #cor.method =  "spearman",#方法  符合正态分布采用，否则用建议使用非参数检验,spearman等
                cor.method = "pearson",
                cor.coef.size = 4,
                xlab = "MARCKS", ## x轴
                ylab = paste("Smooth.muscle"))##y轴 gsub("_","-",gene)
U1



gene <- "SLC40A1"
U2 <- ggscatter(corr_exp,x = gene, #x变量
                 y = "Smooth.muscle" ,#y变量
                 add = "reg.line",##拟合曲线
                 #conf.int = TRUE,##置信区间阴影带
                 cor.coef = TRUE, ##系数
                 color = "black",
                 add.params = list(color = "blue",fill = "lightgray"),
                 #cor.method =  "spearman",#方法  符合正态分布采用，否则用建议使用非参数检验,spearman等
                 cor.method = "pearson",
                 cor.coef.size = 4,
                 xlab = "SLC40A1", ## x轴
                 ylab = paste("Smooth.muscle"))##y轴 gsub("_","-",gene) # p=0.072
U2
gene <- "PELI1"
U3 <- ggscatter(corr_exp,x = gene, #x变量
                y = "Th1.cells" ,#y变量
                add = "reg.line",##拟合曲线
                #conf.int = TRUE,##置信区间阴影带
                cor.coef = TRUE, ##系数
                color = "black",
                add.params = list(color = "blue",fill = "lightgray"),
                #cor.method =  "spearman",#方法  符合正态分布采用，否则用建议使用非参数检验,spearman等
                cor.method = "pearson",
                cor.coef.size = 4,
                xlab = "PELI1", ## x轴
                ylab = paste("Th1.cells"))##y轴 gsub("_","-",gene) # p=0.072

U3


gene <- "STK17B"
U4 <- ggscatter(corr_exp,x = gene, #x变量
                 y = "Tgd.cells" ,#y变量
                 add = "reg.line",##拟合曲线
                 #conf.int = TRUE,##置信区间阴影带
                 cor.coef = TRUE, ##系数
                 color = "black",
                 add.params = list(color = "blue",fill = "lightgray"),
                 #cor.method =  "spearman",#方法  符合正态分布采用，否则用建议使用非参数检验,spearman等
                 cor.method = "pearson",
                 cor.coef.size = 4,
                 xlab = "STK17B", ## x轴
                 ylab = paste("Tgd.cells"))##y轴 gsub("_","-",gene) # p=0.072
U4


library(cowplot) # 合并图 金明
plot_grid(U1,U2,U3,U4, ncol = 2, axis = "1", align = "v") # ncol = 3,, labels = c("A", "B")

# plot_grid(U4,U5, ncol = 1, axis = "1", align = "v") # ncol = 3,, labels = c("A", "B")
# plot_grid(U6,U7,U8,U9,U10,U11, ncol = 3, axis = "1", align = "v") # ncol = 3,, labels = c("A", "B")


# 9 肿瘤纯度  estimate ####
setwd("/Users/emperor/Documents/DataR/HeartD/ImmuneInfiltration/Estimate/")
library(estimate)
Exp2 <- read.table("./Exp2.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE,check.names=F)
exprSet<-Exp2
rownames(exprSet) <- exprSet[, 1]
exprSet <- exprSet[, -1]

dat <- exprSet
dat <- t(dat)
dat <- data.frame(dat)
write.table(dat,"dat.txt",row.names = T,quote = F,sep="\t")

filterCommonGenes(input.f="./dat.txt",  # z这个只是路径到文件的名字
                  output.f="CAD.gct", 
                  id="GeneSymbol")
estimateScore(input.ds = "CAD.gct",
              output.ds="PR_estimate_score.gct", 
              platform="affymetrix")

plotPurity(scores="PR_estimate_score.gct", samples="s516", 
           platform="affymetrix")

scores=read.table("PR_estimate_score.gct",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=t(scores[,3:ncol(scores)])
scores
resustEstimate <- data.frame(id=rownames(scores), scores)
# resustEstimate$id <- substr(resustEstimate$id, 1, 15)
write.table(resustEstimate,"Estimate_score_CAD.txt",row.names = F,quote = F,sep="\t")


sample2 <- read.table("sample2.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE,check.names=F)
# risk$id <- gsub("-", ".", risk$id)
sample2 <- sample2[, c(1,3)]
estimate.sample <- merge(resustEstimate, sample2, by.x = "id", by.y = "GSM", all = F)
write.table(estimate.sample,"Estimate_score_CAD_sample.txt",row.names = F,quote = F,sep="\t")









# sankey 桑基图
# 桑基图
setwd("/Users/emperor/Documents/DataR/HeartD/ImmuneInfiltration/xCell/")
setwd("/Users/emperor/Documents/DataR/HeartD/ImmuneInfiltration.revised/")

sankey <- read.table("./0sankey.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE,check.names=F)

# devtools::install_github("fbreitwieser/sankeyD3")
library(sankeyD3)
sanger <- sankey
sanger$value <- "1"
library(dplyr)
nodes <- data.frame(name=c(as.character(sanger$node1), as.character(sanger$node2)) %>% unique())
sanger$node1 <- match(sanger$node1, nodes$name)-1 
sanger$node2 <- match(sanger$node2, nodes$name)-1
p<-sankeyNetwork(Links = sanger, Nodes = nodes, Source = "node1", Target = "node2",
                 Value = "value", NodeID = "name",nodeWidth =10,units = 'TWh',
                 height=300,width=300,colourScale=JS("d3.scaleOrdinal(d3.schemeCategory10);"),
                 numberFormat=".0f",fontSize = 8)  
p
saveNetwork(p,"0sankey.html")
library(webshot) 
# webshot::install_phantomjs()
webshot("file:///Users/emperor/Documents/DataR/HeartD/ImmuneInfiltration.revised/0sankey.html", 
        "0sankey.pdf")
# 链接：https://www.jianshu.com/p/d103b4b143f3






# 10 Drug ####
setwd("/Users/emperor/Documents/DataR/HeartD/6DrugGene_DGIdb/")

library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
GDSC2_Expr <- readRDS("/Users/emperor/Downloads/DataFiles/DataFiles/Training Data/GDSC2_Expr (RMA Normalized and Log Transformed).rds")
GDSC2_Res <- readRDS("/Users/emperor/Downloads/DataFiles/DataFiles/Training Data/GDSC2_Res.rds")
GDSC2_Res <- exp(GDSC2_Res) # 基因名样本名都在行名或者列名，方框只有数据

# our data
# 列名 也就是样本名不一致没关系，每一行都是基因，基因名在rownames，数据框全部都是数据
risk <- read.table("./risk.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE,check.names=F)
riskdata <- risk[, c(1,4:11)] # 只要样本名，基因
rownames(riskdata) <- riskdata[,1]
riskdata <- riskdata[, -1]
riskdata2 <- t(riskdata)
riskdata2 <- data.frame(riskdata2)
riskdata2 <- as.matrix(riskdata2) # 样本量很多，但是基因只有risk的5基因，运行比较快


calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = riskdata2,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 2, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )

library(data.table)
testPtype <- fread('./calcPhenotype_Output/DrugPredictions.csv', data.table = F)
testPtype[1:4, 1:4]
drug <- data.frame(colnames(testPtype))

# nccn 去查 uvm的化疗药物
Drug.nccn <- c("Temozolomide_1375","Paclitaxel_1080","Trametinib_1372","Gemcitabine_1190",
               "Cisplatin_1005","Vinblastine_1004", "Docetaxel_1007","Docetaxel_1819",
               "Selumetinib_1736",
               "Sorafenib_1085", "Dabrafenib_1373" )

testPtype2 <- testPtype[, colnames(testPtype) %in% Drug.nccn]
testPtype2 <- data.frame(testPtype[,1], testPtype2)
testPtype_select <- testPtype2
colnames(testPtype_select)[1] <- "id"

testPtype_select_risk <- merge(risk, testPtype_select, by = "id", all = F)
write.table(testPtype_select_risk,"testPtype_select_risk.txt",row.names = F,quote = F,sep="\t")



DGIdb <- read.table("./DrugGene_DGIdb.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE,check.names=F)
Gene.immune <- read.table("./Gene.immune.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE,check.names=F)
MER <-  merge(Gene.immune, DGIdb, by.x = "symbol", by.y = "gene_name", all = F)


# 10.1桑基图 sankey ####
setwd("/Users/emperor/Documents/DataR/HeartD/6DrugGene_DGIdb/")
sankey <- read.table("sankey.drug.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE,check.names=F)
# devtools::install_github("fbreitwieser/sankeyD3")
library(sankeyD3)
sanger <- sankey
sanger$value <- "1"
library(dplyr)
nodes <- data.frame(name=c(as.character(sanger$node1), as.character(sanger$node2)) %>% unique())
sanger$node1 <- match(sanger$node1, nodes$name)-1 
sanger$node2 <- match(sanger$node2, nodes$name)-1
p<-sankeyNetwork(Links = sanger, Nodes = nodes, Source = "node1", Target = "node2",
                 Value = "value", NodeID = "name",nodeWidth =10,units = 'TWh',
                 height=300,width=300,colourScale=JS("d3.scaleOrdinal(d3.schemeCategory10);"),
                 numberFormat=".0f",fontSize = 8)  
p
saveNetwork(p,"sankey.ImmuneCell.html")
library(webshot) 
# webshot::install_phantomjs()
webshot("file:///Users/emperor/Documents/DataR/HeartD/ImmuneInfiltration/xCell/sankey.ImmuneCell.html", 
        "sankey.ImmuneCell.pdf")
# 链接：https://www.jianshu.com/p/d103b4b143f3


sankey <- read.table("sankey.drug.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE,check.names=F)
head(sankey)
#     gene          drug score
# 1  CXCR4      CP-55940 43.01
# 2 STK17B     quercetin 99.99
# 3   JUNB Dexamethasone 76.00
# 4  CXCR4   BEVACIZUMAB 21.00

library(googleVis)
Sankeyy <- gvisSankey(sankey, from="gene", to="drug", weight="score",
                      options=list(
                        sankey="{link: {color: { fill: '#d799ae' } },
                    node: { color: { fill: '#a61d4c' },
                    label: { color: '#871b47' } }}"))
q <- plot(Sankeyy)





# 11 roc ####
setwd("/Users/emperor/Documents/DataR/HeartD/8roc/")
# setwd("/Users/emperor/Documents/DataR/HeartD/ImmuneInfiltration.revised/")

Exp2 <- read.table("../GSE42148/Exp2.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE,check.names=F)
# Gene.immune <- read.table("../6DrugGene_DGIdb/Gene.immune.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE,check.names=F)
Gene.immune <- read.table("./0sankey.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE,check.names=F)

Exp2.im <- Exp2[Exp2$Symbol %in% Gene.immune$node1, ]
rownames(Exp2.im) <- Exp2.im[, 1]  
Exp2.im <- Exp2.im[, -1]
Exp2.im <- t(Exp2.im)
Exp2.im <- data.frame(Exp2.im)
Exp2.im2 <- data.frame(GSM=rownames(Exp2.im), Exp2.im)
sample2 <- read.table("../GSE42148/sample2.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE,check.names=F)
sample2 <- sample2[, c(1,3)]
Exp2.im.proce <- merge(sample2, Exp2.im2, by.x = "GSM", by.y = "GSM", all = F)
Exp2.im.proce <- Exp2.im.proce[ , -1]

tmpt <- Exp2.im.proce

# 常规roc
# pdf("mul.ROC.pdf", width = 6, height = 6)
ra<-rainbow(20)

library(pROC)
rocobj<- roc(tmpt$Group =='case', as.numeric(tmpt$CTDSP2),  smooth=TRUE)#感兴趣的基因
plot(rocobj,col="#FF0000FF")
rocobj$auc 

rocobj<- roc(tmpt$Group =='case', as.numeric(tmpt$MARCKS), smooth=TRUE)#感兴趣的基因
plot(rocobj,col="#00FFB2",add=T)
rocobj$auc 

rocobj<- roc(tmpt$Group =='case', as.numeric(tmpt$PELI1), smooth=TRUE)#感兴趣的基因
plot(rocobj,col="#0066FF",add=T)
rocobj$auc 

legend("bottomright", legend=c("CTDSP2, AUC=0.8626", "STK17B, AUC=0.8112", "PELI1, AUC=0.8169"),
       col=c('#FF0000FF',"#00FFB2","#0066FF"),lwd=2,cex = 0.8,
       x.intersp = 2, y.intersp=0.8, adj = c(0.1, 0.4))



rocobj<- roc(tmpt$Group =='case', as.numeric(tmpt$NLRP1), smooth=TRUE)#感兴趣的基因
plot(rocobj,col="#7F00FF")
rocobj$auc 

rocobj<- roc(tmpt$Group =='case', as.numeric(tmpt$STK17B), smooth=TRUE)#感兴趣的基因
plot(rocobj,col="darkblue",add=T)
rocobj$auc 


rocobj<- roc(tmpt$Group =='case', as.numeric(tmpt$JUNB), smooth=TRUE)#感兴趣的基因
plot(rocobj,col="#CC00FF",add=T)
rocobj$auc 

rocobj<- roc(tmpt$Group =='case', as.numeric(tmpt$DHRS7), smooth=TRUE)#感兴趣的基因
plot(rocobj,col="#00B3FF",add=T)
rocobj$auc 


rocobj<- roc(tmpt$Group =='case', as.numeric(tmpt$RILPL2), smooth=TRUE)#感兴趣的基因
plot(rocobj,col="#001AFF",add=T)
rocobj$auc 

rocobj<- roc(tmpt$Group =='case', as.numeric(tmpt$SLC40A1), smooth=TRUE)#感兴趣的基因
plot(rocobj,col="#00FF66",add=T)
rocobj$auc 


legend("bottomright", legend=c("NLRP1, AUC=0.7964", "MARCKS, AUC=0.7632", "JUNB, AUC=0.7108",
                               "DHRS7, AUC=0.6652","RILPL2, AUC=0.6627","SLC40A1, AUC=0.6391"),
       col=c('#7F00FF',"darkblue","#CC00FF","#00B3FF","#001AFF","#00FF66"),lwd=2,cex = 0.8,
       x.intersp = 2, y.intersp=0.8, adj = c(0.1, 0.4))







# 12 validation ####
setwd("/Users/emperor/Documents/DataR/HeartD/9Validation/")
GSE20680_series_matrix <- read.table("GSE20680_series_matrix.txt.gz",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE)
GPL4133 <- read.table("GPL4133-12599.txt",sep = "\t",comment.char = "#", stringsAsFactors = F,header = T, fill=TRUE, quote="")
GPL4133 = data.frame(GPL4133[,c(1,10)])

trim <- function( x ) {
  gsub(" ", "", x)
}
GPL4133$ID = trim(GPL4133$ID)

Exp <- merge(GSE20680_series_matrix,GPL4133,by.x = "ID_REF",by.y = "ID",all=FALSE) # 45015
Exp = Exp[,c(197,2:196)]
colnames(Exp)[1] = "Symbol"
# Exp = data.frame(Exp[-grep("///",Exp$"Symbol"),]) #去除一个探针对应多个Symbol的行(去一对多???
Exp <- Exp[Exp$Symbol != "" ,]  # 32696
# Exp = na.omit(Exp)
meanfun <- function(x) {
  x1 <- data.frame(unique(x[,1]))
  colnames(x1) <- c("Symbol")
  for (i in 2:ncol(x)){
    x2 <- data.frame(tapply(x[,i],x[,1],mean))
    x2[,2] <- rownames(x2)
    colnames(x2) <- c(colnames(x)[i], "Symbol")
    x1 <- merge(x1,x2,by.x = "Symbol",by.y = "Symbol",all=FALSE)
  }
  return(x1)
}
Exp <- meanfun(Exp) # 19749
write.table(Exp,"Exp3.txt",row.names = F,quote = F,sep="\t")


## 12.1 Randomforest ####
# validation
exp3 <- read.table('../9Validation/Exp3.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40

Module_colour_Gene <- read.table('./Module_colour_Gene.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
Blue <- Module_colour_Gene[grep("blue", Module_colour_Gene$type), ]
brown <- Module_colour_Gene[grep("brown", Module_colour_Gene$type), ]
Turquoise <- Module_colour_Gene[grep("turquoise", Module_colour_Gene$type), ]


need1 <- Blue
need1 <- brown
need1 <- Turquoise

# validatioon
ranforData2 <- exp3[exp3$Symbol %in% need1$Symbol, ]
rownames(ranforData2) <- ranforData2[,1]
ranforData2 <- ranforData2[, -1]
ranforData2 <- t(ranforData2)
ranforData2 <- data.frame(ranforData2)
ranforData2 <- data.frame(aa=rownames(ranforData2), ranforData2)


library(randomForest)
ranforData2$aa <- as.factor(ranforData2$aa) # @ 第一列要factor格式才可以，不然会出现“参数不是数值也不是逻辑值：回覆NA”报错
need.rf2 <- randomForest(aa ~ ., data=ranforData2, importance=TRUE, proximity=TRUE)
print(need.rf2)

need.mds2 <- cmdscale(1 - need.rf2$proximity, eig=TRUE)
plot(need.mds2$points, col = rep(c("red", "blue", "green"), each = 50))
need.rf2$importance
varImpPlot(need.rf2, main = "Top 15 - variable importance") # OK
resu2 <- data.frame(need.rf2$importance)
resu2 <- resu2[, 196:197]
resu2 <- data.frame(symbol = rownames(resu2), resu2)



# validation
write.table(resu2,"/Users/emperor/Documents/DataR/HeartD/9Validation/Validation_1Blue_RandomForest_Importance.txt",row.names = F,quote = F,sep="\t")
write.table(resu2,"/Users/emperor/Documents/DataR/HeartD/9Validation/Validation_2brown_RandomForest_Importance.txt",row.names = F,quote = F,sep="\t")
write.table(resu2,"/Users/emperor/Documents/DataR/HeartD/9Validation/Validation_3turquoise_RandomForest_Importance.txt",row.names = F,quote = F,sep="\t")


# 12.1.1 ROC-random forest ####
# 教程 https://www.bioinfo.online/articleList/20228952669.html
library(pROC)

sample3 <- read.table('../9Validation/sample3.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
sample3 <- sample3[, c(1,3)]
ranforData3 <- merge(sample3, ranforData2, by.x = "GSM", by.y= "aa", all = F)


library(RColorBrewer)
mycol<-brewer.pal(11, "Spectral")

plot.roc(ranforData3$Group, 
         print.auc=T,
         col="#00FFFF", # 颜色 mycol[9]
         predict(need.rf2, newdata = ranforData2,type="vote")[,2])

legend("bottomright", legend=c("turquoiseRForest, AUC=0.563"),
       col="#00FFFF",lwd=2,cex = 0.8,
       x.intersp = 2, y.intersp=0.8, adj = c(0.1, 0.4))




# 12.1.2 cloudwind ####
RaF <- resu2
data2 <- RaF

# 环形柱状图 金明
# 教程 https://www.jianshu.com/p/865aa9023b27
df<-data.frame(individual=paste("feature",seq(1,24),sep=""),value=c(24:1))
df$id<-seq(1,nrow(df))
library(ggplot2)
p<-ggplot(df,aes(x=as.factor(id),y=value,fill=as.factor(id)))+
  geom_bar(stat="identity")
p+coord_polar()+theme_bw()+ theme(legend.position="none")+ylim(-2,24)

# our result 

RaF1 <- read.table('../9Validation/Validation_1Blue_RandomForest_Importance.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
RaF1 <- read.table('../9Validation/Validation_2brown_RandomForest_Importance.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
RaF1 <- read.table('../9Validation/Validation_3turquoise_RandomForest_Importance.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40


RaF1 <- RaF1[order(RaF1$MeanDecreaseGini, decreasing = T), ]
RaF1$MeanDecreaseGini <- RaF1$MeanDecreaseGini*10
# RaF1 <- as.factor(RaF1$symbol, levels=RaF1$symbol)
RaF1$id<-seq(1,nrow(RaF1))

p<-ggplot(RaF1[1:10, ], aes(x=as.factor(id), y=MeanDecreaseGini, fill=as.factor(symbol)))+
  geom_bar(stat="identity")
p+coord_polar()+theme_bw()+ theme(legend.position="none")
# +ylim(-2,4)

write.table(RaF1,"/Users/emperor/Documents/DataR/HeartD/9Validation/Validation_1Blue_RandomForest_Importance.legend.txt",row.names = F,quote = F,sep="\t")
write.table(RaF1,"/Users/emperor/Documents/DataR/HeartD/9Validation/Validation_2brown_RandomForest_Importance.legend.txt",row.names = F,quote = F,sep="\t")
write.table(RaF1,"/Users/emperor/Documents/DataR/HeartD/9Validation/Validation_3turquoise_RandomForest_Importance.legend.txt",row.names = F,quote = F,sep="\t")


# 12.2 svm 三组 #### 
exp3 <- read.table('./Exp3.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40

Module_colour_Gene <- read.table('./Module_colour_Gene.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
Blue <- Module_colour_Gene[grep("blue", Module_colour_Gene$type), ]
brown <- Module_colour_Gene[grep("brown", Module_colour_Gene$type), ]
Turquoise <- Module_colour_Gene[grep("turquoise", Module_colour_Gene$type), ]


need <- Blue
need <- brown
need <- Turquoise

svmData1 <- exp3[exp3$Symbol %in% need$Symbol, ]
rownames(svmData1) <- svmData1[, 1]
svmData1 <- svmData1[, -1]

svmData1 <- t(svmData1)
#              ABHD5    ACOX1    ACSL4     ANXA3
# GSM518885 11.070190 8.445527 8.613062 10.373875
# GSM518886 10.012830 7.771118 7.860759  8.707333

# 分成训练集和测试集
smp.size <- floor(0.6*nrow(svmData1))
train.ind <- sample(seq_len(nrow(svmData1)), smp.size)
train <- svmData1[train.ind, ]
train1 <- data.frame(id=rownames(train), train) # 138 hang

test <- svmData1[-train.ind, ]
test1 <- data.frame(id=rownames(test), test) # 60 hang

#挤入临床信息 二分类
sample3 <- read.table('../9Validation/sample3.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
sample3 <- sample3[, c(1,3)]
train2 <- merge(sample3, train1, by.x = "GSM", by.y = "id", all = F)
rownames(train2) <- train2[, 1]
train2 <- train2[, -1] # 去掉gsm
test2 <- merge(sample3, test1, by.x = "GSM", by.y = "id", all = F)
rownames(test2) <- test2[, 1]
test2 <- test2[, -1] # 去掉gsm


library(e1071)
train2$Group <- as.factor(train2$Group) # 如果报错Need numeric dependent variable for regression.，解决办法是加上as.factor
model <- svm(formula = Group ~ ., data = train2, kernel = "linear")
summary(model)



# 训练集的混淆矩阵
train.pred <- predict(model, train2)
table(real=train2$Group, predict=train.pred)
# blue      predict
# real      case control
# case      87       0
# control    1      29

# brown     predict
# real      case control
# case      82      3
# control   8      24

# Turquoise predict
# real      case control
# case      86       0
# control    0      31
confus.matrix <- table(real=train2$Group, predict=train.pred)
sum(diag(confus.matrix))/sum(confus.matrix) 
# blue 0.9745763  brown 0.7881356  turquoise 1.0 
# validation blue 0.991453 , BROWN 0.9059829  turquoise 1.0 

# 测试集预测
# test2$Group <- as.factor(test2$Group)
test.pred <- predict(model, test2) # 一定要只有数据，不能有GSM

table(real=test2$Group, predict=test.pred)
confus.matrix2 <- table(real=test2$Group, predict=test.pred)
sum(diag(confus.matrix2))/sum(confus.matrix2) # 0.5641026   0.6025641  0.7307692


# 12.2.1 ROC ####
# 教程：https://www.jianshu.com/p/c926c191361e
# 提取模型预测值并进行格式处理
pred_1 <- as.factor(model$decision.values)
pred_1 <- as.ordered(pred_1)
modelroc_1 <- roc(train2$Group, pred_1)
modelroc_1
plot(modelroc_1, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2), grid.col=c("green", "red"), max.auc.polygon=TRUE, auc.polygon.col="skyblue", print.thres=TRUE)
# 可视化展示,使用add=TRUE将第二个模型添加到图形中
# plot.roc(modelroc_2, add=TRUE, col="green",print.thres=TRUE) 


# 12.2.2 SVM value ####
# 把基因和svm值组合起来
svmValue <- model[["x.scale"]]$`scaled:center`
svmValue <- data.frame(value=svmValue)
svmValue <- data.frame(symbol=rownames(svmValue), svmValue)

write.table(svmValue,"Validation_1Blue_SVM_Value.txt",row.names = F,quote = F,sep="\t")
write.table(svmValue,"Validation_2brown_SVM_Value.txt",row.names = F,quote = F,sep="\t")
write.table(svmValue,"Validation_3Turquoise_SVM_Value.txt",row.names = F,quote = F,sep="\t")


# 12.3 LASSO #### 
# 教程:https://www.jianshu.com/p/2995ae5f244e
# 教程更好： https://www.weinformatics.cn/44483852be/

exp3 <- read.table('./Exp3.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40

Module_colour_Gene <- read.table('./Module_colour_Gene.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
Blue <- Module_colour_Gene[grep("blue", Module_colour_Gene$type), ]
brown <- Module_colour_Gene[grep("brown", Module_colour_Gene$type), ]
Turquoise <- Module_colour_Gene[grep("turquoise", Module_colour_Gene$type), ]

needLasso <- Blue
needLasso <- brown
needLasso <- Turquoise

LassoData1 <- exp3[exp3$Symbol %in% needLasso$Symbol, ]
rownames(LassoData1) <- LassoData1[, 1]
LassoData1 <- LassoData1[, -1]

sample3 <- read.table('../9Validation/sample3.txt',sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill = TRUE) #40
sample3 <- sample3[, c(1,3)]
meta <- sample3
meta$Group <- gsub("case", "1", meta$Group)
meta$Group <- gsub("control", "0", meta$Group)
meta$Group <- as.numeric(meta$Group)

# lasso
# cox1=c('RP4-616B8.5','RP11-389G6.3',"AP000696.2",'CTD-2377D24.6',"LINC01559" ,"LINC00629","AC005062.2","LINC01018")
# ##cox1为定义的任何向量，比如里面含有你所感兴趣的基因名称
# exprSet = exprSet[cox1,]#表达谱数据中提取你所感兴趣的基因表达
exprSet=as.data.frame(LassoData1)
dim(exprSet)
x=t(exprSet) # 样本在左rownames，基因在列，列名，colnames
y=meta$Group#提取病例的预后情况（0与1）列，meta是临床数据
library(glmnet)
model_lasso <- glmnet(x, y,nlambda=10, alpha=1)
print(model_lasso)
set.seed(13098)
cv_fit <- cv.glmnet(x=x, y=y, nlambda = 1000,alpha = 1)
plot(cv_fit)

# cv_fit$lambda.min blue
# [1] 0.008128992, log(lambda): -4.812318
# validation 0.03726663 log(lambda): -3.289657

# cv_fit$lambda.min brown
# [1] 0.004071331
# validation 0.05116375, log(lambda): -2.972724

# cv_fit$lambda.min turquoise
# [1] 0.01013239 log(lambda): -4.592018
# validation 0.02004464, log(lambda): -3.909793

lasso.prob <- predict(cv_fit, newx=x , s=c(cv_fit$lambda.min,cv_fit$lambda.1se) )
head(lasso.prob)
re=cbind(y ,lasso.prob)
head(re)
re1=as.data.frame(re)
# write.csv(re,file='')
# 结果有0，1，0表示control。1表示case
write.table(re1,"Validation_1Blue_LassoValue_sample.txt",row.names = F,quote = F,sep="\t")
write.table(re1,"Validation_2brown_LassoValue_sample.txt",row.names = F,quote = F,sep="\t")
write.table(re1,"Validation_3Turquoise_LassoValue_sample.txt",row.names = F,quote = F,sep="\t")

# 12.3.1 gene value ####
# 要基因对应的值  
# 两条虚线分别指示了两个特殊的λ值,一个是lambda.min,一个是lambda.1se,
# 这两个值之间的lambda都认为是合适的。lambda.1se构建的模型最简单，即使用的基因数量少，
# 而lambda.min则准确率更高一点，使用的基因数量更多一点。
# 用这两个λ值重新建模或者用min这个重新建模

model_lasso_min <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.min)
# model_lasso_1se <- glmnet(x=x, y=y, alpha = 1, lambda=cv_fit$lambda.1se)
head(model_lasso_min$beta)
choose_gene_min=rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]
# choose_gene_1se=rownames(model_lasso_1se$beta)[as.numeric(model_lasso_1se$beta)!=0]

lassoValue <- model_lasso_min[["beta"]]@x
lassoValue <- data.frame(choose_gene_min,lassoValue)
write.table(lassoValue,"Validation_1Blue_LassoValue.txt",row.names = F,quote = F,sep="\t")
write.table(lassoValue,"Validation_2brown_LassoValue.txt",row.names = F,quote = F,sep="\t")
write.table(lassoValue,"Validation_3Turquoise_LassoValue.txt",row.names = F,quote = F,sep="\t")



# 12.3.2 箱线图 ####
re=as.data.frame(re)
colnames(re)=c('event','prob_min','prob_1se')
re$event=as.factor(re$event)
library(ggpubr) 
p1 = ggboxplot(re, x = "event", y = "prob_min",
               color = "event", palette = "jco",
               add = "jitter")+ stat_compare_means()
p2 = ggboxplot(re, x = "event", y = "prob_1se",
               color = "event", palette = "jco",
               add = "jitter")+ stat_compare_means()
library(patchwork)
p1+p2
p1


# 12.3.3 ROC-Lasso ####
library(ROCR)
library(caret)
# 自己预测自己
#min
pred_min <- prediction(re[,2], re[,1])
auc_min = performance(pred_min,"auc")@y.values[[1]]
perf_min <- performance(pred_min,"tpr","fpr")
plot(perf_min,colorize=FALSE, col="#00FFFF") 
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
text(0.8,0.2, labels = paste0("AUC = ",round(auc_min,3)))

#1se
pred_1se <- prediction(re[,3], re[,1])
auc_1se = performance(pred_1se,"auc")@y.values[[1]]
perf_1se <- performance(pred_1se,"tpr","fpr")
plot(perf_1se,colorize=FALSE, col="red") 
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
text(0.8,0.2, labels = paste0("AUC = ",round(auc_min,3)))





