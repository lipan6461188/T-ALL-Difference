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


#把reads count转成RPKM，RPKM表示 该基因的读段数/[ (样本总读段数/100,0000) * (基因的长度/1000) ]
myRPKM = mycounts
for(column in c(1:dim(mycounts)[2]))
{
	reads = mycounts[,column]
	total_reads_and_length = sum(mycounts[,column])/1000000  * ( mylength / 1000 )
	myRPKM[,column] = reads / total_reads_and_length
}

#归一化
myRPKM2 = rpkm(mycounts, long = mylength, k = 0, lc = 1)
myUQUA = uqua(mycounts, long = mylength, lc = 0.5, k = 0.5)
myTMM = tmm(mycounts, long = 1000, lc = 0)

#低表达读段过滤，11447 features are to be kept
myfilt <- filtered.data(myRPKM, factor = myFactors$Tissue, norm = TRUE, depth = NULL, method = 1, cv.cutoff = 100, cpm = 10, p.adj = "fdr")
mydata2 <- readData(data = myfilt,length = mylength,biotype = mybiotypes,chromosome = mychroms, factors = myFactors)

#再次生成报表
QCreport(mydata2, samples = NULL, factor = "Tissue", norm = TRUE)

#差异表达分析
mynoiseqbio = noiseqbio(input = mydata2, r = 20, plot = TRUE, a0per = 0.9, random.seed = 12345, norm = "n", filter = 0, factor="Tissue")

#选取差异表达的gene
mynoiseq.deg = degenes(mynoiseqbio, q = 0.95, M = NULL)
mynoiseq_up = degenes(mynoiseqbio, q = 0.95, M = "up")
mynoiseq_down = degenes(mynoiseqbio, q = 0.95, M = "down")

#把基因的名字和ID对应上去
index = match(table=anno[,1],x=rownames(mynoiseq_up))
mynoiseq_up <- cbind( GeneName=anno[index,7],mynoiseq_down)
mynoiseq_up <- mynoiseq_up[mynoiseq_up$log2FC>1,]
index = match(table=anno[,1],x=rownames(mynoiseq_down))
mynoiseq_down <- cbind( GeneName=anno[index,7],mynoiseq_down)
mynoiseq_down <- mynoiseq_down[mynoiseq_down$log2FC< -1,]

#存储文件
write.table(file="T-ALL-up-12-22.txt",x=mynoiseq_up,sep="\t",col.names=T,row.names=T,quote=F)
write.table(file="T-ALL-down-12-22.txt",x=mynoiseq_down,sep="\t",col.names=T,row.names=T,quote=F)

#作图
pdf("T-ALL差异表达12-22.pdf")
DE.plot(mynoiseqbio, q = 0.95, graphic = "expr", log.scale = TRUE)
DE.plot(mynoiseqbio, q = 0.95, graphic = "MD")
DE.plot(mynoiseqbio, chromosomes = c(1, 2), log.scale = TRUE, join = FALSE, q = 0.95, graphic = "chrom")
DE.plot(mynoiseqbio, chromosomes = c(1:3), q = 0.8, graphic = "distr")
dev.off()

#两个文献中提到的gene
result <- mynoiseqbio@results[[1]]
Myb <- result[rownames(result)=="ENSG00000118513",]
Hbp1 <- result[rownames(result)=="ENSG00000105856",]
