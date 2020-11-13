#!/usr/bin/Rscript
rm(list = ls())
#
G = c("Gpru", "Hbet","LhavRL", "Mper","Pcos", "Psem", "Grac", "Hkan", "Lmic","Nmon", "Pcle", "Pgab", "Pgra","Noccaea","AT")
maxcontig <- vector("numeric", length(G))
mediancontig <- vector("numeric", length(G))
maxorf <- vector("numeric", length(G))
medianorf <- vector("numeric", length(G))
nb <- vector("numeric", length(G))

for (i in 1:length(G))
{
	print(G[i])
	file = paste("/home/chris/COGs/",G[i],".tab",sep = "")
	print(file)
	A <- as.data.frame(read.table(file, sep = "\t"))
	# longueur des contigs
	A[,2] ->  l_contig
	l_contig = as.numeric(l_contig)
	# longueur de la plus grande ORF
	A[,3] ->  l_orf
	l_orf = as.numeric(l_orf)

	# diagramme des valeurs observées

	hist(l_contig, breaks = 100, main=paste("genome",G[i],sep=" of "),xlim=c(0,max(l_contig)),ylim=c(0,16000), col=rgb(0,1,0,1))  # first histogram
	par(new=TRUE)
	hist(l_orf, breaks = 100, main=paste("genome",G[i],sep=" of "),xlim=c(0,max(l_contig)),col=rgb(1,0,0,1),ylim=c(0,16000))
	
	maxcontig[i] =  max(l_contig)
	mediancontig[i] = median(l_contig)
	maxorf[i] =  max(l_orf)
	medianorf[i] = median(l_orf)
	nb[i] = length(l_contig)	
	#mean(l_orf)
	#sd(l_orf)	
}
A<-gl(15,15,15,labels=G)
data <-cbind(mediancontig,medianorf)
rownames(data) <- levels(A)
#barplot(mediancontig,names.arg=levels(A))
barplot(t(data),beside=T,legend.text=colnames(data),
col=c("grey50","grey80"),ylab="longueur(nt)")
warnings()
