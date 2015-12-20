import re

handle = open('gencode.v23.annotation.gtf','r')
output = open('anno_in_NOISeq.txt','w')

pattern = '^chr(\w+)\s+\w+\s+gene\s+(\w+)\s+(\w+).*gene_id "(\w+).*gene_type "(\S+)";'
reg = re.compile(pattern)
line = handle.readline()
while(line):
	if re.match('#',line) is not None:
		line = handle.readline()
		continue
	
	element = re.findall(reg,line)
	if len(element) == 0:
		line = handle.readline()
		continue
	element = element[0]
	mystr = '%s\t%s\t%s\t%s\t%d\t%s\n' % (element[3], element[0], element[1], element[2], int(element[2]) - int(element[1]) + 1,element[4])
	output.writelines(mystr)
	line = handle.readline()

handle.close()
output.close()
