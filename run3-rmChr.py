import re

I = open('GRCh38.primary_assembly.genome.fa','r')

line = I.readline().strip()

fg = False
while line:
        chrname = re.match('>(\w+)',line)
        if chrname != None:
                if chrname.group(1).startswith('chr'):
                        fg = True
                        print '>%s' % ( chrname.group(1) )
                else:
                        fg = False
                line = I.readline().strip()
                continue
        if re.match('#',line) != None:
                line = I.readline().strip()
                continue
        if fg:
                print line
        line = I.readline().strip()

I.close()
