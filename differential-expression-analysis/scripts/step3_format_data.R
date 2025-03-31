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
## aim : format data for phylogenetic differencial expression analysis
## input : - count data, length data (output of script `step2_construct_matrix_OG_from_pairs.R`)
##         - condition table for each sample
## output : a `phyloCompData` object with formatted expression data and tree.

rm(list=ls())
library(here)
library(compcodeR)

################################################################################
## File management
################################################################################
condName <- "Hyperaccu" 

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
stopifnot(isTRUE(all(colnames(counts_noNA)==colnames(length_noNA))))

# select factors for DE analysis  
colnames(sample_annotation)
table(sample_annotation$Localisation,sample_annotation$Climat)
colData <- sample_annotation[colnames(counts_noNA),c("ID.espece","Famille", "Localisation", "Tolerance","Hyperaccu")]
# three populations are nickel tolerant ( AM ScorA, KP ScorC, GAL ScorB) one is not (PF scorD)
table(colData$Tolerance,colData$Hyperaccu)
# merge localisation and familly
colData$LocAndFam <- paste(colData$Famille,colData$Localisation)
colData$id.species <- colData$ID.espece
colData$condition <- colData[[condName]]

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

## match tree and data
counts_noNA <- counts_noNA[, match(tree_rep$tip.label, colnames(counts_noNA))]
length_noNA <- length_noNA[, match(tree_rep$tip.label, colnames(length_noNA))]
colData <- colData[match(tree_rep$tip.label, rownames(colData)), ]

################################################################################
# format compcodeR
################################################################################
## rename key column "condition" in colData (done already)
## colnames(colData)[grep(condName, colnames(colData))] <- "condition"
info.parameters <- list(dataset = "nickel_cpd", uID = "1", tree = tree_rep)
cpd <- phyloCompData(count.matrix = counts_noNA, 
                     sample.annotations = colData, 
                     info.parameters = info.parameters,
                     length.matrix = as.matrix(length_noNA),
                     tree = tree_rep)
## check obj conformity
check_phyloCompData(cpd)

## save data to rds file
dataset <- "nickel_cpd"
## save in data and not in "results_directory"
#dataset_file <- here(file.path(results_directory, paste0(dataset, ".rds")))
dataset_file <- here(file.path("data", paste0(dataset, ".rds")))
saveRDS(cpd, file = dataset_file)
