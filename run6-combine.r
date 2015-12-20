ROOT_PATH <- "/data/lipan/RNA_expresionn/annotation/"

map_file <- paste(ROOT_PATH,"map.txt",sep="")
map <- read.csv(map_file,sep="\t",head=F)
column <- as.data.frame(strsplit(as.character(map[,2])," "))
map[,2] <- t(column[1,])
map[,3] <- t(column[2,])

output_file_path <- paste(ROOT_PATH,"output2/",sep="")
all_files <- list.files(output_file_path)
all_files <- all_files[grepl(all_files,patter=".*ReadsPerGene.*",perl=T)]

genes <- ""
l <- list()
for(i in c(1:dim(map)[1]))
{
	file_path <- paste(output_file_path,as.character(map[i,1]),".ReadsPerGene.out.tab",sep="")
	myFile <- read.csv(file_path,sep="\t",head=F,skip=4)[,2] 
	if(i == 1){ genes <- as.character(read.csv(file_path,sep="\t",head=F,skip=4)[,1]) }
	sample <- as.character(map[i,3])
	if(sample %in% names(l))
	{
		l[[sample]] <- l[[sample]] + myFile
	}else{
		l[[sample]] <- myFile
	}
	print(i)
}

split_genes <- strsplit(split="\\.",x=genes)
genes <- c()
for(i in c(1:length(split_genes)))
{
	genes <- c(genes,split_genes[[i]][1])
}
frame <- data.frame(row.names=genes)
for(i in l)
{
	frame <- cbind(frame,i)
}
colnames(frame) <- names(l)

write.table(x=frame,file="mRNA_Expression.txt",sep="\t",quote=FALSE,col.names=TRUE,row.names=TRUE)



