library(NOISeq)

mycounts<-read.csv("/Users/lee/Desktop/project/guo/mRNA_Expression.txt",sep="\t",head=T)

myFactors <- data.frame(Tissue=c("Control","Control","Control","Control","TALL",
	"TALL","TALL","TALL","TALL","TALL","TALL","TALL","TALL"))

anno <- read.csv("/Users/lee/Desktop/project/guo/anno_in_NOISeq.txt",sep="\t",head=F)

mylength <- as.numeric(anno[,5])
names(mylength) <- anno[,1]
mybiotypes <- as.character(anno[,6])
names(mybiotypes) <- anno[,1]
mychroms <- data.frame(row.names=anno[,1],Chr=anno[,2],GeneStart=anno[,3],GeneEnd=anno[,4])
mydata <- readData(data = mycounts,length = mylength,biotype = mybiotypes,chromosome = mychroms, factors = myFactors)

#生成报表，查看信息
QCreport(mydata, samples = NULL, factor = "Tissue", norm = FALSE)

#归一化
myRPKM = rpkm(mycounts, long = mylength, k = 0, lc = 1)

#低表达读段过滤，11447 features are to be kept
myfilt = filtered.data(myRPKM, factor = myFactors$Tissue, norm = FALSE, depth = NULL, method = 1, cv.cutoff = 100, cpm = 2, p.adj = "fdr")
mydata2 <- readData(data = myfilt,length = mylength,biotype = mybiotypes,chromosome = mychroms, factors = myFactors)

#再次生成报表
QCreport(mydata2, samples = NULL, factor = "Tissue", norm = TRUE)

#查看信息
head(assayData(mydata)$exprs)
head(pData(mydata))
head(featureData(mydata)@data)

set.seed(123)
mycounts2 = mycounts
mycounts2[, 1:4] = mycounts2[, 1:4] + runif(nrow(mycounts2) * 4, 3, 5)
myfactors = data.frame(myFactors, batch = c(rep(1, 4), rep(2, 9)))
mydata2 = readData(mycounts2, factors = myFactors)

 myPCA = dat(mydata2, type = "PCA")
 par(mfrow = c(1, 2))
 explo.plot(myPCA, factor = "Tissue")
 explo.plot(myPCA, factor = "batch")







