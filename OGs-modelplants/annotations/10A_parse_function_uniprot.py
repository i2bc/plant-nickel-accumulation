import untangle
swissprot= untangle.parse('annot_source/uniprot_sprot.xml')
fh = open ('annot_files/swissprot.fa', 'w')

for each_entry in  swissprot.uniprot.entry:
	fh.write (">"+each_entry.name.cdata+" ")
	try:		
		fh.write("swissprot:"+each_entry.protein.recommendedName.fullName.cdata+"\n")
	except IndexError:
		print('no name')
		
	fh.write (each_entry.sequence.cdata+"\n")
	

	


