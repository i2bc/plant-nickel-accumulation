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

# Program to compute scores of the differential expression analysis methods.
# Use input from script 03_simulations_data_analysis.R

# NOTE: The full result table from our own run (result of 04_simulations_data_analysis_results.R)
#       is provided in the repository for convenience

library(here)
library(compcodeR)

################################################################################
# Files
################################################################################
condName <- "Hyperaccu"
datestamp_day_real_nickel <- "2025-03-17"
results_directory <- paste0(datestamp_day_real_nickel, "_results_nickel_", condName)
dataset <- "nickel_cpd"
dataset_file <- here(file.path(results_directory, paste0(dataset, ".rds")))

################################################################################
## Create result table
################################################################################

## Find all analysis files that deal with the simulation dataset
all_method_files <- list.files(path = results_directory,
                               pattern = paste0(dataset, "..*.rds"),
                               full.names = TRUE)

##################################################################
## Comparison scores
##################################################################
## File table with all analyses for a given dataset
file.table <- data.frame(input.files = all_method_files, stringsAsFactors = FALSE)

## Parameters and computed scores
th <- 1e-2
parameters <- list(incl.nbr.samples = NULL,
                   incl.replicates = NULL,
                   incl.dataset = dataset, 
                   incl.de.methods = NULL, 
                   fdr.threshold = th, tpr.threshold = th, typeI.threshold = th,
                   ma.threshold = th, fdc.maxvar = 1500, overlap.threshold = th,
                   fracsign.threshold = th, mcc.threshold = th, 
                   nbrtpfp.threshold = th, 
                   comparisons = c("auc",
                                   "mcc",
                                   "fdr", "tpr",
                                   "fdrvsexpr",
                                   "typeIerror", "fracsign", "nbrsign", "nbrtpfp",
                                   "maplot",
                                   "fdcurvesone",
                                   "rocone",
                                   "overlap", "sorensen",
                                   "scorevsexpr",
                                   "scorevssignal",
                                   "correlation"
                   ))

## Compute scores and save results in a table
runComparison(file.table = file.table,
              parameters = parameters,
              output.directory = here(results_directory, dataset),
              save.result.table = TRUE,
              knit.results = TRUE)

## get list of significative genes
signif_genes <- NULL
for (ff in all_method_files) {
  res <- readRDS(here(ff))
  signif_genes[[res@method.names$methodMod]] <- which(res@result.table$adjpvalue < th)
}
saveRDS(signif_genes, file = file.path(results_directory, paste0("list_all_signif_genes_.rds")))

library(UpSetR)
pdf(file = file.path(results_directory, paste0("upset_genes.pdf")), width = 11.7, height = 8.8)
upset(fromList(signif_genes), order.by = "freq",
      nsets = 5, text.scale = 2,
      sets.x.label = paste0("adjusted p-values < ", th))
dev.off()


