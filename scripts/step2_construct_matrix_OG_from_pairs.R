# Copyright (C) {2020} {MM, MG, SJ, CD}
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

## author : Mélina Gallopin from a initial version done by Marie Michel and Sarah Jelassi
## date : 21/01/2020
## aim : construc OG expression matrix
## input : count data, list of contigs associated to OGs
## output : OG expression count matrix and OG length matrix

rm(list=ls())
library(here)


library(DESeq2)
library(FactoMineR)
library(dplyr)
library(UpSetR)
library(tibble)
mypath <- here()

########################################################################
## Open files
########################################################################
dir.create(here("data"))
setwd(dir = here("raw_data", "length"))
myfiles_length <- list.files(pattern = ".RData", full.names=TRUE)
setwd(dir = here("raw_data", "sum_counts"))
myfiles_counts <- list.files(pattern = ".RData", full.names=TRUE)


########################################################################
## collect matrices for Phco
########################################################################  
name_spe="Phco"
setwd(mypath)
setwd("./raw_data")


# load counts 
setwd(dir = "./sum_counts")
load(myfiles_counts[grep(paste0("(",name_spe,")"), myfiles_counts)])
countData <- sum_counts_OGspartage
#countData[is.na(countData)] <-0
table(rowMeans(countData)==0)


# load length
setwd(dir = "../length")
load(myfiles_length[grep(paste0("(",name_spe,")"), myfiles_length)])
dim(sum_OGslength)
table(is.na(sum_OGslength))
matrice_length = sum_OGslength
# mymeans= rowMeans(sum_OGslength,na.rm=T)
# for(i in 1:nrow(sum_OGslength)){
#   matrice_length[i, is.na(sum_OGslength[i,])] <- mymeans[i]
# }
# table(is.na(matrice_length))

# reorder columns
countData <- countData[,order(colnames(countData))]
matrice_length <- matrice_length[,order(colnames(matrice_length))]
##sum_counts_OGspartage <- sum_counts_OGspartage[,order(colnames(sum_counts_OGspartage))]
#colnames(sum_counts_OGspartage) <- rownames(colData)

print("check")
print(colnames(matrice_length))
print(colnames(countData))



## write matrices 
mergedCounts <- countData
mergedCounts <- rownames_to_column(mergedCounts,var="OGid")
mergedLengths <- matrice_length
mergedLengths <- rownames_to_column(mergedLengths,var="OGid")
dim(mergedCounts)
dim(mergedLengths)

########################################################################
## collect matrices for Phco
########################################################################
#list_for_grep <- c("Hkan","Gpru", "Pgab", "Phco","Pcos","Nmon","Scor")
list_for_grep <- c("Hkan","Gpru","Pgab","Pcos","Nmon","Scor")
for( name_spe in list_for_grep){  
  
  setwd(mypath)
  setwd("./raw_data")
  
  
  # load counts 
  setwd(dir = "./sum_counts")
  load(myfiles_counts[grep(paste0("(",name_spe,")"), myfiles_counts)])
  countData <- sum_counts_OGspartage
  #countData[is.na(countData)] <-0
  table(rowMeans(countData)==0)
  
  
  # load length
  setwd(dir = "../length")
  load(myfiles_length[grep(paste0("(",name_spe,")"), myfiles_length)])
  dim(sum_OGslength)
  table(is.na(sum_OGslength))
  matrice_length = sum_OGslength
  # mymeans= rowMeans(sum_OGslength,na.rm=T)
  # for(i in 1:nrow(sum_OGslength)){
  #   matrice_length[i, is.na(sum_OGslength[i,])] <- mymeans[i]
  # }
  # table(is.na(matrice_length))
  
  
  
  countData<- countData[,order(colnames(countData))]
  matrice_length <- matrice_length[,order(colnames(matrice_length))]
  ##sum_counts_OGspartage <- sum_counts_OGspartage[,order(colnames(sum_counts_OGspartage))]
  #colnames(sum_counts_OGspartage) <- rownames(colData)
  
  print("check")
  print(colnames(matrice_length))
  print(colnames(countData))
  
  
  
  ########################################################################
  ## collect matrix count
  ########################################################################  
  print(dim(mergedCounts))
  print(dim(countData))
  head(countData)
  head(mergedCounts)
  countData <- rownames_to_column(countData,var="OGid")
  mergedCounts <- full_join(mergedCounts,countData,by="OGid")
  head(mergedCounts)
  print(dim(mergedCounts))
  matrice_length <- rownames_to_column(matrice_length,var="OGid")
  mergedLengths <- full_join(mergedLengths,matrice_length,by="OGid")
  
}
setwd(mypath)
setwd("./data")
print(getwd())
print(head(mergedCounts))
dim(mergedCounts)
write.csv(mergedCounts,file="counts.txt", row.names = FALSE)
counts <- column_to_rownames(mergedCounts, var = "OGid")
save(counts,file="counts.RData")

print(head(mergedLengths))
print(dim(mergedCounts))
print(dim(mergedLengths))
write.csv(mergedLengths,file="lengths.txt", row.names = FALSE)
lengths <- column_to_rownames(mergedLengths, var = "OGid")
save(lengths,file="lengths.RData")
