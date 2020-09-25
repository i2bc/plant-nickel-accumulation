#!/usr/bin/python
# the script filters the main results of hmmscan (--tblout output) for evalues < cutoff and separates the results according to the number of domains (one or more than one)
def filter_by_evalue(min_evalue, infile, outfile1, outfile2):
	ifh = open(infile, "r")
	lignes = ifh.readlines()
	ofh1 = open(outfile1, "w")
	ofh2 = open(outfile2, "w")
	for l in lignes:
		if (not l.startswith('#')):
			hmmresult = l.split()
			fullevalue = hmmresult[4]
			bestevalue = hmmresult[7]
			nbdomain = hmmresult[15]
			if (float(fullevalue) < min_evalue and float(bestevalue) < min_evalue):
				if (int(nbdomain) == 1):
					ofh1.write(l)
				else:
					ofh2.write(l)
	ifh.close()       
	ofh1.close()
	ofh2.close()

###hmmscan --domtblout output fields
#0:target name    
#1:accession   
#2:tlen 
#3:query name           
#4:accession   
#5:qlen

##full sequence      
#6:E-value  
#7:score  
#8:bias
#9: nombre de domaine
#10: numero du domaine

##This domain    
#11:c-Evalue  
#12:i-Evalue  
#13:score  
#14:bias

#hmm coord  
#15:from    
#16:to

##ali coord  
#17:from
#18:to

##env coord
#19:from    
#20:to  
#21:acc 
#22:description of target
# the script filters the results of hmmscan (--domtblout output) for evalues < cutoff (see cutoff in the Snakefile)

def filter_dom_by_evalue (min_evalue, infile, outfile):
	ifh = open(infile, "r")
	lignes = ifh.readlines()
	ofh = open(outfile, "w")
	for l in lignes:
		if (not l.startswith('#')):
			hmmresult = l.split()
			fullevalue = float(hmmresult[6])
			domevalue = float(hmmresult[12])			
			if (fullevalue < min_evalue and domevalue < min_evalue):				
				ofh.write(l)
				
	ifh.close()       
	ofh.close()

# the script filters the results of hmmscan (--domtblout output) for query coverage > cutoff (see cutoff in the Snakefile)
def filter_dom_by_coverage (min_coverage, infile, outfile):
	import numpy
	ifh = open(infile, "r")
	lignes = ifh.readlines()
	ofh = open(outfile, "w")
	for l in lignes:
		if (not l.startswith('#')):
			hmmresult = l.split()
			qlen = int(hmmresult[5])
			ali_len = int(hmmresult[18]) - int(hmmresult[19])
			cov = ali_len/qlen*100
			if (cov > min_coverage):				
				ofh.write(l)
				
	ifh.close()       
	ofh.close()

# the script generate a tabular file (sequence length, target group, sequence ID) from the best hits of the sequences of all genomes
def add2group(min_coverage, infile, outfile):
	ifh = open(infile, "r")
	lignes = ifh.readlines()
	ofh = open(outfile, "a")
	lastseq = 'none'
	for l in lignes:
		if (not l.startswith('#')):
			hmmresult = l.split()
			if (lastseq != hmmresult[3]):
				target_group = hmmresult[0]
				query_seq = hmmresult[3]
				qlen = int(hmmresult[5])
				ali_len = int(hmmresult[18]) - int(hmmresult[19])
				cov = ali_len/qlen*100
				if (cov > min_coverage):					
					ofh.write(hmmresult[5]+'|'+target_group + ';' + query_seq + "\n")
				lastseq = hmmresult[3]
	ifh.close()       
	ofh.close()

def tab2gff (infile):
	ifh = open(infile, "r")
	lignes = ifh.readlines()

	for l in lignes:
		splitedl = l.split()
		contig = splitedl[0]
		lcontig = splitedl[1]
		strand = splitedl[2]
		if(strand == "plus"):
			strand = "+"
			start = splitedl[5]
			end = splitedl[6]
		else:
			strand = "-"
			start = int(lcontig) - int(splitedl[6]) +1
			end = int(lcontig) - int(splitedl[5]) +1
		lORF = splitedl[3]
		cov = splitedl[4]		 
		
		meth = splitedl[7]
		gff = contig+"\tsearchORF\tORF\t"+str(start)+"\t"+str(end)+"\t"+cov+"\t"+strand+"\t.\tmeth="+meth+"; length="+lORF
		print (gff)							
				
	ifh.close()       
