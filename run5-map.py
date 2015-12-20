import os
import re
from multiprocessing import Pool

def match(x):
	mat = re.findall('.*fastq$',x)
	if len(mat) != 0:
		return True
	else:
		return False

def worker(eachFile):
	outFileNamePrefix = ROOT_PATH + '/annotation/output3/' + eachFile + '.'
	inputfile1 = ROOT_PATH+'/RNA-Seq/'+eachFile+'_1.fastq'
	inputfile2 = ROOT_PATH+'/RNA-Seq/'+eachFile+'_2.fastq'
	command = 'STAR --runThreadN 2 --runMode alignReads' + ' --genomeDir ' + index_file_path + ' --readFilesIn ' + inputfile1 + ' ' + inputfile2 + ' --quantMode GeneCounts --outFilterMultimapNmax 5 --outFilterMismatchNmax 2 --outFileNamePrefix ' + outFileNamePrefix
	os.system(command)

if __name__ == '__main__':
	ROOT_PATH = '/data/lipan/RNA_expresionn'
	all_files = os.listdir(ROOT_PATH+'/RNA-Seq')
	fastq_end = map(match,all_files)
	file_list = []
	for i in range(0,len(fastq_end)):
		if not fastq_end[i]:
			file_list.append(all_files[i])

	index_file_path = ROOT_PATH + '/annotation/index'

	i = 0
	while(i<len(file_list)):
		if(i+10>len(file_list)):
			p = Pool(len(file_list) - i)
			p.map(worker,file_list[i:])
			i += 10
		else:
			p = Pool(10)
			p.map(worker,file_list[i:i+10])
			i += 10
		p.join()
