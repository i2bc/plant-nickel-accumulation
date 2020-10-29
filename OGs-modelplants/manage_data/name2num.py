import sys
(script, numerotation, groupfile) = sys.argv

nfh = open(numerotation,'r')
num_sequences = nfh.read()
all_proteines = num_sequences.split("\n")

name2num = dict()
for protein in all_proteines:
	try:
		(numero, name) = protein.split (': ')
		name2num[name] = numero		
	except:
		pass 
nfh.close
print name2num["PPE_001G21990"]
print name2num["RC27699G00050"]


gfh = open(groupfile, 'r')
outfh = open(groupfile+'_renum', 'w')
all_groups_as_string = gfh.read()
all_groups_as_string = all_groups_as_string.rstrip()
all_groups_as_list = all_groups_as_string.split("\n")
for group in all_groups_as_list:
	proteins = group.split(';')
	line = ';'.join([name2num[p] for p in proteins]) + "\n"	
	outfh.write(line)	
