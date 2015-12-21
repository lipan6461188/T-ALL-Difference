import re

h = open('./mRNA_Expression.txt','r')
line = h.readline().strip()
print line
line = h.readline().strip()
while(line):
	g = re.match('(\w+)\S+(\s+.*)',line)
	print '%s%s' % (g.group(1),g.group(2))
	line = h.readline().strip()

h.close()