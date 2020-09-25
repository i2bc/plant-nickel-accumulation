# Copyright (C) {2020} {PB, MG}
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

## author : Mélina Gallopin and Paul Bastide
## date : 25/09/2020
## aim : phylogenetic differential expression analysis
## input : count data, length data (output of script `step2_construct_matrix_OG_from_pairs.R`) and condition table for each sample
## output : list of OGs differentially expressed

rm(list=ls())
library(here)
# devtools::install_github("pbastide/phylocompcodeR")
library(phylocompcodeR)

################################################################################
## File management
################################################################################
condName <- "Hyperaccu" # change here Hyperaccu or Tolerance
# condName <- "Tolerance"

## File
datestamp_day_real_nickel <- format(Sys.time(), "%Y-%m-%d")
results_directory <- paste0(datestamp_day_real_nickel, "_results_nickel_", condName)
dir.create(here(results_directory))

################################################################################
## read data
################################################################################
# counts contains the counts for each COG (sum of the counts assigned to each contigs associated with the given COG) for each sample
counts <- read.csv(here("data/counts.txt"), sep=",")
# info on samples
sample_annotation <- read.csv(here("raw_data/samples_annotation.txt"), sep="\t")
# length_matrix contains the "length" of each COG (sum of the lengths of the contigs associated with the given COG) for each sample
length_matrix <- read.csv(here("data/lengths.txt"), sep=",")

## Rename Psychotria_clementis -> Psychotria_revoluta in data
# colnames(comptages_COG) <- sub("Pcle", "Prev", colnames(comptages_COG))
# colnames(length_COG) <- sub("Pcle", "Prev", colnames(length_COG))
# rownames(sample_annotation) <- sub("Pcle", "Prev", rownames(sample_annotation))
# sample_annotation$ID.echantillon <- sub("Pcle", "Prev", sample_annotation$ID.echantillon)
# sample_annotation$ID.espece <- sub("Pcle", "Prev", sample_annotation$ID.espece)
# sample_annotation$Nom.espece <- sub("Psychotria clementis", "Psychotria revoluta", sample_annotation$Nom.espece)


# format input data
rownames(counts) <- counts[,1]
counts <- counts[,-1]
rownames(length_matrix) <- length_matrix[,1]
length_matrix <- length_matrix[,-1]
rownames(sample_annotation) <- sample_annotation$ID.echantillon

# remove NA
counts_noNA <- counts[complete.cases(counts),]
dim(counts_noNA)
length_noNA <- length_matrix[complete.cases(counts),colnames(counts_noNA)]
dim(length_noNA)
# check names consistency
colnames(counts_noNA)==colnames(length_noNA)

# select factors for DE analysis  
colnames(sample_annotation)
table(sample_annotation$Localisation,sample_annotation$Climat)
colData <- sample_annotation[colnames(counts_noNA),c("ID.espece","Famille", "Localisation", "Tolerance","Hyperaccu")]
# three populations are nickel tolerant ( AM ScorA, KP ScorC, GAL ScorB) one is not (PF scorD)
table(colData$Tolerance,colData$Hyperaccu)
# merge localisation and familly
colData$LocAndFam <- paste(colData$Famille,colData$Localisation)

################################################################################
## Tree
################################################################################
library(ape)
## Get the tree
tree <- read.tree(file = here("raw_data/TreeHyperMai2020.nwk"))
plot(tree)

## Species names
correspondances <- unique(sample_annotation[, c("Nom.espece", "ID.espece")])
# Format
correspondances[, "Nom.espece"] <- sub(" ", "_", correspondances[, "Nom.espece"])
# Typos
correspondances[, "Nom.espece"] <- sub("Senocio_coronatus", "Senecio_coronatus", correspondances[, "Nom.espece"])
correspondances[, "Nom.espece"] <- sub("Microthlaspi_perfoliatum", "Microthalspi_perfoliatum", correspondances[, "Nom.espece"])
correspondances[, "Nom.espece"] <- sub("Homalium_kanaliense", "Homalium_kanalense", correspondances[, "Nom.espece"])
correspondances[, "Nom.espece"] <- sub("Homalium_betulifoium", "Homalium_betulifolium", correspondances[, "Nom.espece"])
# Match
tree_data_cor <- match(tree$tip.label, correspondances[, "Nom.espece"])
data_tree_cor <- match(correspondances[, "Nom.espece"], tree$tip.label)
# Species in the tree NOT in data
tree$tip.label[is.na(tree_data_cor)]
# Species in data NOT in the tree
correspondances[is.na(data_tree_cor), "Nom.espece"]

## Format Tree
# Rename species with ID
tree$tip.label <- as.character(correspondances[match(tree$tip.label, correspondances[, "Nom.espece"]), "ID.espece"])
plot(tree)
# Add replicates
tree_rep <- tree
for (tip_label in tree$tip.label) {
  replis <- colnames(counts_noNA)[grepl(tip_label, colnames(counts_noNA))]
  for (rep in replis) {
    tree_rep <- phytools::bind.tip(tree_rep, tip.label = rep,
                                   where = which(tree_rep$tip.label == tip_label))
  }
}
# Remove original tips
tree_rep <- ape::drop.tip(tree_rep, tree$tip.label)
# Plot
plot(tree_rep)

################################################################################
# format compcodeR
################################################################################
## rename key column "condition" in colData
colnames(colData)[grep(condName, colnames(colData))] <- "condition"
info.parameters <- list(dataset = "nickel_cpd", uID = "1", tree = tree_rep)
cpd <- compData(count.matrix = counts_noNA, 
                sample.annotations = colData, 
                info.parameters = info.parameters,
                length.matrix = as.matrix(length_noNA))
## check obj conformity
check_compData(cpd)

## save data to rds file
dataset <- "nickel_cpd"
dataset_file <- here(file.path(results_directory, paste0(dataset, ".rds")))
saveRDS(cpd, file = dataset_file)

################################################################################
# Analyses
################################################################################

##################################################################
## phylolm - OU - RPKM - log2
##################################################################
method <- "phylolm"
method_name <- paste0(method, ".RPKM.log2.OU")

runDiffExp(data.file = dataset_file,
           result.extent = method_name, 
           Rmdfunction = paste0(method, ".createRmd"), 
           output.directory = here(results_directory), 
           norm.method = "TMM",
           model = "OUfixedRoot",
           measurement_error = TRUE,
           lengthNormalization = "RPKM",
           dataTransformation = "log2",
           extraDesignFactors = c("Localisation"))

generateCodeHTMLs(here(results_directory, paste0(dataset, "_", method_name, ".rds")), results_directory)
