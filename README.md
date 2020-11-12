# Analyses for paper "Differential Expression on Publicly Available RNA-Seq Datasets Using Phylogenetic Comparative Methods to Understand Nickel Accumulation in Plants"

* Folder OGs-modelplants: Contains all scripts and  open source code from external softwares (see Dependencies file) needed to built and annotate orthologs groups .
	* input: fasta formatted proteomes, swissprot database, metacyc database
	* output: fasta formatted orthologs groups, flat files of groups numbers with associated annotations (putative functions, EC numbers and Gene Onthology entries)

* Folder OGs-nickelplants : contains the scripts needed to enriched the orthologous groups built from model plants (see OGs-modelplants) with putative proteins translated from de novo transcriptomes. 
	* inputs : transcriptome (fasta formatted contigs from de novo assemblies)
	* output : a list of peptide length, contig name, ortholog group number

* Folder differential-expression-analysis : contains everything needed to perform differential analysis 
	* inputs  : the expression table for each sample, metadata for each sample, list of contigs and associated OG number for all species, phylogenetic tree 
	* output : list of DE OGs, figures

