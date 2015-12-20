from os import listdir,system
from re import match

path = '/data/lipan/RNA_expresionn/RNA-Seq'
fqDumpBaseCommand = "fastq-dump --split-files --split-3 --outdir " + path + " "
files = listdir(path)
for myfile in files:
	if match('SRR',myfile) == None:
		continue
	filepath = path + '/' +myfile
	system(fqDumpBaseCommand + filepath)
