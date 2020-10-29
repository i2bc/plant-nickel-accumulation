
import sys
(script, numerotation, groupfile) = sys.argv

nfh = open(numerotation,'r')
num_sequences = nfh.read()
all_proteines = num_sequences.split("\n")

num2name = {}
for protein in all_proteines:
	try:
		(numero, name) = protein.split (': ')
		num2name[numero] = name
	except:
		pass 
nfh.close

print num2name["0_12"]
gfh = open(groupfile, 'r')
outfh = open(groupfile+'_renamed', 'w')
all_groups_as_string = gfh.read()
all_groups_as_string = all_groups_as_string.rstrip()
all_groups = all_groups_as_string.split("\n")
for group in all_groups:
	proteins = group.split(';')
	line = ';'.join([num2name[p] for p in proteins]) + "\n"	
	outfh.write(line)
