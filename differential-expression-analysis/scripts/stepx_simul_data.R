# Copyright (C) {2021} {PB, MG}
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

# Program to simulate data using empirical dataset 
# Use input from script step3_format_data.R

# NOTE : This is a computation intensive script
# NOTE : The full result table from our own run (result of 04_simulations_data_analysis.R)
#        is provided in the repository for convenience

## Install correct version of compcodeR
# devtools::install_github("csoneson/compcodeR", ref = "phylocomp")
rm(list=ls())
library(here)
library(compcodeR)
library(knitr)
library(ggplot2)
library(ape)

################################################################################
## Load Data
################################################################################
## read data
data_nickel <- readRDS(here("data", "nickel_cpd.rds"))

counts_noNA <- data_nickel@count.matrix
length_noNA <- data_nickel@length.matrix
colData <- data_nickel@sample.annotations
# on a déjà la colonne Hyperaccu
#colnames(colData)[colnames(colData) == "condition"] <- "Hyperaccu"
tree <- data_nickel@tree

################################################################################
## Files management
################################################################################
# here_dir <- sub("save", "work", here())
here_dir <- here()

datestamp_day <- format(Sys.time(), "%Y-%m-%d")
simus_directory <- paste0(datestamp_day, "_simulations_nickel")
simus_directory <- file.path(here_dir, simus_directory)
dir.create(simus_directory)

################################################################################
## Analyse Data to get parameters
################################################################################
## run DESeq2 to get parameters estimation
library(DESeq2)
# Correct by families
dds <- DESeqDataSetFromMatrix(counts_noNA, colData, ~ condition )
# normalisation factors for library size
dds <-  estimateSizeFactors(dds)
size_fac <- sizeFactors(dds)
# compute normalization with length
mat_size_fac <- matrix(size_fac, ncol=length(size_fac), nrow=length(counts_noNA[,1]), byrow=T)
normFactors <- (mat_size_fac*length_noNA) / exp(rowMeans(log(mat_size_fac*length_noNA)))
normalizationFactors(dds) <- as.matrix(normFactors)
# analysis
dds <- DESeq2::DESeq(dds, fitType = 'parametric', test = 'Wald', betaPrior = TRUE)
dds <- DESeq2::replaceOutliersWithTrimmedMean(dds)
dds <- DESeq2::DESeq(dds, fitType = 'parametric', test = 'Wald', betaPrior = TRUE)
res <- DESeq2::results(dds, independentFiltering = TRUE, cooksCutoff = TRUE)

################################################################################
## Get empirical moments for simulation
################################################################################
dispersions_nickel <- DESeq2::dispersions(dds)
seqdepth_nickel <- mean(colSums(counts_noNA))
minfact_nickel <- min(colSums(counts_noNA)/median(colSums(counts_noNA)))
maxfact_nickel <- max(colSums(counts_noNA)/median(colSums(counts_noNA)))
relmeans_nickel <- rowMeans(counts_noNA[,colData$condition == "oui"])
n.vars <- dim(counts_noNA)[1]

plot(relmeans_nickel,dispersions_nickel)
plot(log(relmeans_nickel), log(dispersions_nickel))

plot(log(relmeans_nickel), log(dispersions_nickel),
     xlim=c(0,15),ylim=c(0,15))
abline(0,1)

hist(log(dispersions_nickel))

anyNA(dispersions_nickel)

################################################################################
## Format Tree
################################################################################
n <- length(unique(sub("_.*", "", tree$tip.label))) # number of species
N <- length(tree$tip.label) # total number of observations (with replicates)

tree$edge.length <- tree$edge.length / max(diag(vcv(tree))[1:N]) # normalize tree to unit height

id_species <- sub("_.*", "", tree$tip.label)
id_species <- sub("Scor[A-D]", "Scor", id_species) # ScorA, ScorB and ScorC are the same
id_species <- factor(id_species)
names(id_species) <- tree$tip.label

species_names <- sapply(1:length(colnames(length_noNA)), function(i) strsplit(colnames(length_noNA),"_")[[i]][1])
species_names <- sapply(1:length(species_names), function(i) strsplit(species_names,"[[:upper:]]$")[[i]][1])

################################################################################
## Format Lengths and compute empirical stats
################################################################################
unique_length <- length_noNA[, !duplicated(species_names)]

mean_lengths <- rowMeans(unique_length)
var_lengths <- apply(unique_length,1,var)
disp_lengths <- (var_lengths - mean_lengths) / mean_lengths^2 # size = mu^2/(var-mu), disp = 1 / size
disp_lengths[disp_lengths < 0] <- 1

plot(log(mean_lengths), log(var_lengths))
abline(0,1)

################################################################################
## Empirical fraction of data explained by the tree
################################################################################
# Normalisation TPM
nf <- edgeR::calcNormFactors(counts_noNA / length_noNA, method = 'TMM')
lib.size <- colSums(counts_noNA / length_noNA) * nf
data.norm <- sweep((counts_noNA + 0.5) / length_noNA, 2, lib.size + 1, '/')
data.norm <- data.norm * 1e6
# Transformation log2
data.trans <- log2(data.norm)
rownames(data.trans) <- rownames(counts_noNA)
# phylolm pagel lambda analysis
lambdas <- apply(data.trans, 1, function(x) phylolm::phylolm(x ~ colData$condition, phy = tree, model = "lambda")$optpar)
hist(lambdas)

# fit empirical distribution
estima_beta <- RobPer::betaCvMfit(lambdas[lambdas > 1e-7])
x <- seq(0, 1, 0.01)
hist(lambdas, breaks = 100, freq = FALSE)
lines(x, dbeta(x, estima_beta[1], estima_beta[2]))
set.seed(1289)
lambdas_sim <- rbeta(length(lambdas), estima_beta[1], estima_beta[2])
hist(lambdas_sim,  breaks = 100, freq = FALSE)

# OU analysis to get empirical selection strength
oufits <- apply(data.trans, 1, function(x) phylolm::phylolm(x ~ colData$condition, phy = tree,
                                                            model = "OUfixedRoot", measurement_error = TRUE,
                                                            lower.bound = list(sigma2_error = (.Machine$double.eps)^0.5)))

get_lambda_error <- function(sigma2, sigma2_error, h_tree) {
  return(sigma2 * h_tree / (sigma2_error + sigma2 * h_tree))
}
get_gamma <- function(phyfit) {
  return(phyfit$sigma2 / 2 / phyfit$optpar)
}
tree_height <- function(tree) {
  return(max(ape::node.depth.edgelength(tree)))
}
get_lambda_error_OU <- function(alpha, sigma2, sigma2_error) {
  tree_model <- phylolm::transf.branch.lengths(tree, "OUfixedRoot",
                                               parameters = list(alpha = alpha))$tree
  tilde_t <- tree_height(tree_model) / (2 * alpha)
  lambda_ou_error <- get_lambda_error(sigma2, sigma2_error, tilde_t)
  return(lambda_ou_error)
}
get_lambda_error_OU_lm <- function(phyfit) {
  get_lambda_error_OU(phyfit$optpar, phyfit$sigma2, phyfit$sigma2_error)
}
alphas <- sapply(oufits, function(x) x$optpar)
sig2_err <- sapply(oufits, function(x) x$sigma2_error)
gammas <- sapply(oufits, get_gamma)
lambdas_OU <- sapply(oufits, get_lambda_error_OU_lm)
hist(alphas, breaks = 100)
hist(log(2) / alphas, breaks = 100)
hist(gammas, breaks = 100)
hist(lambdas_OU, breaks = 100)

# fit empirical distribution lambda_OU
estim_beta_OU <- RobPer::betaCvMfit(lambdas_OU[lambdas_OU > 1e-7])
x <- seq(0, 1, 0.01)
hist(lambdas_OU, breaks = 100, freq = FALSE)
lines(x, dbeta(x, estim_beta_OU[1], estim_beta_OU[2]))
set.seed(1289)
lambdas_OU_sim <- rbeta(length(lambdas_OU), estim_beta_OU[1], estim_beta_OU[2])
hist(lambdas_OU_sim, breaks = 100, freq = FALSE)

# fit empirical distribution alpha
fitgamma <- robust::gammaRob(alphas[alphas < 50])
x <- seq(0, 50, 0.01)
hist(alphas, breaks = 100, freq = FALSE)
lines(x, dgamma(x, shape = fitgamma$estimate["shape"], scale = fitgamma$estimate["scale"]))
set.seed(1289)
alphas_OU_sim <- rgamma(length(alphas), shape = fitgamma$estimate["shape"], scale = fitgamma$estimate["scale"])
hist(alphas_OU_sim, breaks = 100, freq = FALSE)


################################################################################
## Parameters for the simulation
################################################################################
# samples.per.cond <- N / 2
n.diffexp <- 250

Nrep <- 5 # number of replicates
all_effect_size <- c(1.0, 2.0, 3.0)

selection.strength <- log(2) / 0.5 # half life is 50% of tree height
all_selection_strength <- list("emp", "sim", selection.strength)

#all_prop_var_tree <- list("emp", 0.6, 0.8, 0.9, 1.0)
all_prop_var_tree <- list("emp", "sim", 0.8)

all_fact_disp <- c(1)

#all_use_lengths <- c("with_lengths", "no_lengths")
all_use_lengths <- c("with_lengths")

################################################################################
## Conditions design
################################################################################
all_cond_types <- NULL
all_conds <- NULL

## hyperaccu 
id_cond <- colData$condition
id_cond[id_cond=="oui"] <- 1
id_cond[id_cond=="non"] <- 2
id_cond <- as.numeric(id_cond)
levels(id_cond) <- c(1, 2)

names(id_cond) <- rownames(colData)

# reorder
id_cond <- id_cond[match(tree$tip.label, names(id_cond))]

plot(tree)
tiplabels(pch = 21, col = id_cond, bg = id_cond)

all_cond_types <- c(all_cond_types, "condition")
all_conds[["condition"]] <- id_cond


################################################################################
## Simulations 
################################################################################
all_tree_types <- NULL

tree_type <- "real_tree"
all_tree_types <- c(all_tree_types, tree_type)

all_cond_types_tree <- NULL
all_cond_types_tree[[tree_type]] <- all_cond_types
all_model_process_tree <- NULL
all_model_process_tree[[tree_type]] <- c("BM", "OU")
# all_model_process_tree[[tree_type]] <- c("OU")

# cond_type <- all_cond_types_tree[[tree_type]][1]
# use_lengths <- all_use_lengths[1]
# model_process <- all_model_process_tree[[tree_type]][1]
# effect_size <- all_effect_size[1]
# prop_var_tree <- all_prop_var_tree[[1]]
# sel_strength <- all_selection_strength[1]
# fact_disp <- all_fact_disp[1]

for (cond_type in all_cond_types_tree[[tree_type]]) {
  for (use_lengths in all_use_lengths) {
    for (model_process in all_model_process_tree[[tree_type]]) {
      for (effect_size in all_effect_size) {
        for(prop_var_tree in all_prop_var_tree) {
          for(sel_strength in all_selection_strength) {
            for(fact_disp in all_fact_disp) {
         
          
              dataset <- paste(tree_type, cond_type, model_process, use_lengths, effect_size, prop_var_tree, sel_strength, fact_disp, sep = "_")
              
              pvt <- prop_var_tree
              sl <- sel_strength
              if ((pvt == "emp" || sl == "emp" ||  pvt == "sim" || sl == "sim") && pvt != sl) break
              if (pvt == "emp") {
                if (model_process == "BM") pvt <- lambdas
                if (model_process == "OU") pvt <- lambdas_OU
              } else if (pvt == "sim") {
                if (model_process == "BM") pvt <- lambdas_sim
                if (model_process == "OU") pvt <- lambdas_OU_sim
              }
              if (model_process == "OU") {
                if (sl == "emp") {
                  sl <- alphas
                } else if (sl == "sim") {
                  sl <- alphas_OU_sim
                }
              }
              if (model_process == "BM") sl <- 0.0
              
              set.seed(18570823)
              for (i in 1:Nrep) {
                simus_obj <- generateSyntheticData(dataset = dataset,
                                                   n.vars = n.vars,
                                                   # samples.per.cond = samples.per.cond,
                                                   n.diffexp = n.diffexp,
                                                   repl.id = i,
                                                   seqdepth = seqdepth_nickel,
                                                   minfact = minfact_nickel,
                                                   maxfact = maxfact_nickel,
                                                   relmeans = relmeans_nickel,
                                                   dispersions = dispersions_nickel * fact_disp,
                                                   fraction.upregulated = 0.5,
                                                   between.group.diffdisp = FALSE,
                                                   filter.threshold.total = 1,
                                                   filter.threshold.mediancpm = 0,
                                                   fraction.non.overdispersed = 0,
                                                   random.outlier.high.prob = 0,
                                                   random.outlier.low.prob = 0,
                                                   single.outlier.high.prob = 0,
                                                   single.outlier.low.prob = 0,
                                                   effect.size = effect_size,
                                                   output.file = file.path(simus_directory, paste0(dataset, "_", i, ".rds")),
                                                   tree = tree,
                                                   prop.var.tree = pvt,
                                                   id.condition = all_conds[[cond_type]],
                                                   id.species = id_species,
                                                   model.process = model_process,
                                                   selection.strength = sl,
                                                   lengths.relmeans = if (use_lengths == "with_lengths") mean_lengths else NULL,
                                                   lengths.dispersions = if (use_lengths == "with_lengths") disp_lengths else NULL
                )
                
                if (i == 1) { ## Show summary only for the first iteration
                  summarizeSyntheticDataSet(data.set = file.path(simus_directory, paste0(dataset, "_", i, ".rds")),
                                            output.filename = file.path(simus_directory, paste0(dataset, "_", i, "_check.html")))
                }
              }
            }
          }
        }
      }
    }
  }
}

save(Nrep,
     all_effect_size,
     all_use_lengths,
     all_tree_types,
     all_cond_types_tree,
     all_model_process_tree,
     all_conds,
     all_prop_var_tree,
     all_selection_strength,
     all_fact_disp,
     dispersions_nickel,
     relmeans_nickel,
     seqdepth_nickel,
     minfact_nickel,
     maxfact_nickel,
     mean_lengths,
     disp_lengths,
     lambdas,
     lambdas_OU,
     lambdas_OU_sim,
     lambdas_sim,
     alphas,
     alphas_OU_sim,
     n.diffexp,
     id_species,
     tree,
     file = file.path(simus_directory, "simulation_parameters.RData"))

