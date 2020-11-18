#!/usr/bin/Rscript
rm(list = ls())
#
G = c("Scor" ,"Gpru", "Grac", "Hbet", "Hkan","Mper", "Ncfi", "Nmon", "Pcos", "Pgab", "Pgra", "Prev", "Phco", "Phlu", "Psem", "ath")
maxcontig <- vector("numeric", length(G))
mediancontig <- vector("numeric", length(G))
maxorf <- vector("numeric", length(G))
medianorf <- vector("numeric", length(G))
nb <- vector("numeric", length(G))

for (i in 1:length(G))
{
	print(G[i])
	file = paste("../../NickelPlants_data_article/COGs/orf_tab/",G[i],".tab",sep = "")
	print(file)
	A <- as.data.frame(read.table(file, sep = "\t"))
	# longueur des contigs
	A[,2] ->  l_contig
	l_contig = as.numeric(l_contig)
	# longueur de la plus grande ORF
	A[,4] ->  l_orf
	l_orf = as.numeric(l_orf)

	# diagramme des valeurs observées
	#hist(l_contig, breaks = 100, main=paste("genome",G[i],sep=" of "),xlim=c(0,max(l_contig)),ylim=c(0,16000), col=rgb(0,1,0,1))  # first histogram
	#par(new=TRUE)
	#hist(l_orf, breaks = 100, main=paste("genome",G[i],sep=" of "),xlim=c(0,max(l_contig)),col=rgb(1,0,0,1),ylim=c(0,16000))
	
	maxcontig[i] =  max(l_contig)
	print (c("contigs maximum length: ", maxcontig[i]))
	mediancontig[i] = median(l_contig)
	print (c("contigs median length: ",mediancontig[i]))
	maxorf[i] =  max(l_orf)
	print (c("ORF maximum length: ", maxorf[i]))
	medianorf[i] = median(l_orf)
	print (c("ORF median length: ",medianorf[i]))
	nb[i] = length(l_contig)	
	print (c("number of contigs: ", nb[i]))

}
A<-gl(15,15,15,labels=G)
mat <-cbind(mediancontig,medianorf)
rownames(mat) <- levels(A)

#save(t(mat),file="mat_speci.RData")
barplot(t(mat), beside = TRUE, main = "", horiz=F, col=c("grey50","grey80"),  xlab = "Species",cex.names=1, las=2, ylab="length (nt)")
legend("center", 
       legend = c("Contigs","ORFs"), 
       fill = c("grey50","grey80"), ncol = 2,
       cex = 0.75)


warnings()
