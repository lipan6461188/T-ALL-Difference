library(NOISeq)

mycounts<-read.csv("/Users/lee/Desktop/project/guo/mRNA_Expression2.txt",sep="\t",head=T)

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
myRPKM = rpkm(mycounts, long = mylength, k = 0.5, lc = 1)
myUQUA = uqua(mycounts, long = mylength, lc = 0.5, k = 0.5)
myTMM = tmm(mycounts, long = 1000, lc = 0)

#低表达读段过滤，11447 features are to be kept
myfilt <- filtered.data(myUQUA, factor = myFactors$Tissue, norm = TRUE, depth = NULL, method = 1, cv.cutoff = 100, cpm = 1, p.adj = "fdr")
mydata2 <- readData(data = myfilt,length = mylength,biotype = mybiotypes,chromosome = mychroms, factors = myFactors)

#再次生成报表
QCreport(mydata2, samples = NULL, factor = "Tissue", norm = TRUE)

#差异表达分析
mynoiseqbio = noiseqbio(input = mydata2, r = 20, plot = TRUE, a0per = 0.9, random.seed = 12345, norm = "n", filter = 0, factor="Tissue")

#选取差异表达的gene
mynoiseq.deg = degenes(mynoiseqbio, q = 0.95, M = NULL)
mynoiseq.deg1 = degenes(mynoiseqbio, q = 0.95, M = "up")
mynoiseq.deg2 = degenes(mynoiseqbio, q = 0.95, M = "down")

#把基因的名字和ID对应上去
index = match(table=anno[,1],x=rownames(mynoiseq.deg))
mynoiseq.deg <- cbind(mynoiseq.deg, GeneName=anno[index,7])

#存储文件
write.table(file="T-ALL mRNA差异表达.txt",x=mynoiseq.deg,sep="\t",col.names=T,row.names=T,quote=F)

#作图
DE.plot(mynoiseqbio, q = 0.95, graphic = "expr", log.scale = TRUE)
DE.plot(mynoiseqbio, q = 0.95, graphic = "MD")
DE.plot(mynoiseqbio, chromosomes = c(1, 2), log.scale = TRUE, join = FALSE, q = 0.95, graphic = "chrom")
DE.plot(mynoiseqbio, chromosomes = c(1:3), q = 0.8, graphic = "distr")

#两个文献中提到的gene
result <- mynoiseqbio@results[[1]]
Myb <- result[rownames(result)=="ENSG00000118513",]
Hbp1 <- result[rownames(result)=="ENSG00000105856",]
