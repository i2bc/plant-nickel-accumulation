library(here)
library(compcodeR)
library(doParallel)

################################################################################
# Files
################################################################################
condName <- "Hyperaccu"
datestamp_day_real_nickel <- "2025-03-17"
results_directory <- paste0(datestamp_day_real_nickel, "_results_nickel_", condName)
dataset <- "nickel_cpd"
dataset_file <- here(file.path(results_directory, paste0(dataset, ".rds")))

################################################################################
# Analyses
################################################################################

## Parameters for all methods
# TPM normalisation
lnorm <- "TPM"
# log2 transformation
ltrans <- "log2"
# eBayes with trend
trend <- "with_trend"

####################################################################
## Limma cor
####################################################################
method <- "lengthNorm.limma"

method_name <- paste(method, lnorm, ltrans, trend, sep = ".")

method_info <- list(method = method,
                    lengthNormalization = lnorm,
                    transformation = ltrans,
                    methodProcess = NA,
                    eBayes = TRUE,
                    trend = trend,
                    cor = "cor",
                    regCor = TRUE,
                    nulldisp = NA,
                    ddf_method = NA)
method_info$methodMod <- paste0("limma_cor")

## Run analysis
el_time <- system.time(
  runDiffExp(data.file = dataset_file,
             result.extent = method_name,
             Rmdfunction = paste0(method, ".createRmd"),
             output.directory = here(results_directory),
             norm.method = "TMM",
             length.normalization = lnorm,
             data.transformation = ltrans,
             trend = (trend == "with_trend"),
             block.factor = "id.species",
             extra.design.covariates = c("Localisation"))
)
method_info$time <- el_time[3]

## Save parameters
res_file_name <- file.path(results_directory, paste0(dataset, "_", method_name, ".rds"))
res_tmp <- readRDS(res_file_name)
res_tmp@method.names <- c(res_tmp@method.names, method_info)
saveRDS(res_tmp, file = res_file_name)

## generate the analysis code
generateCodeHTMLs(res_file_name, results_directory)


##################################################################
## phylolm
##################################################################
method <- "phylolm"

foreach (sproc = c("BM", "OUfixedRoot")) %do% {
  
  ## Method id
  method_name <- paste(method, lnorm, ltrans, sproc, sep = ".")
  
  method_info <- list(method = method,
                      lengthNormalization = lnorm,
                      transformation = ltrans,
                      methodProcess = sproc,
                      eBayes = FALSE,
                      trend = NA,
                      cor = NA,
                      regCor = NA,
                      nulldisp = NA,
                      ddf_method = NA)
  method_info$methodMod <- paste0(method_info$method, "_",
                                  method_info$methodProcess)
  
  ## Run analysis
  el_time <- system.time(
    runDiffExp(data.file = dataset_file,
               result.extent = method_name,
               Rmdfunction = paste0(method, ".createRmd"),
               output.directory = here(results_directory),
               norm.method = "TMM",
               model = sproc,
               measurement_error = TRUE,
               extra.design.covariates = c("Localisation"),
               length.normalization = lnorm,
               data.transformation = ltrans,
               lower.bound = list(sigma2_error = (.Machine$double.eps)^0.5),
               upper.bound = list(lambda = 1 / (1 + (.Machine$double.eps)^0.5)))
  )
  method_info$time <- el_time[3]
  
  ## Save parameters
  res_file_name <- file.path(results_directory, paste0(dataset, "_", method_name, ".rds"))
  res_tmp <- readRDS(res_file_name)
  res_tmp@method.names <- c(res_tmp@method.names, method_info)
  saveRDS(res_tmp, file = res_file_name)
  
  ## generate the analysis code
  generateCodeHTMLs(res_file_name, results_directory)
}

##################################################################
## phylolimma - with lengths
##################################################################
method <- "phylolimma"
regCor <- TRUE
eBayes <- TRUE
ddf_method <- "Samples"

foreach (sproc = c("BM", "OUfixedRoot")) %do% {
  
  ## Method id
  method_name <- paste(method, lnorm, ltrans, sproc, trend, paste0("eBayes_", eBayes), paste0("regCor_", regCor), paste0("ddf_method_", ddf_method), sep = ".")
  
  method_info <- list(method = method,
                      lengthNormalization = lnorm,
                      transformation = ltrans,
                      methodProcess = sproc,
                      eBayes = eBayes,
                      trend = trend,
                      cor = NA,
                      regCor = ifelse(regCor, "reg", "no_reg"),
                      nulldisp = NA,
                      ddf_method = ddf_method)
  
  method_info$methodMod <- paste0("phylolimma_",
                                  method_info$methodProcess)
  
  ## Remove some combinations of parameters
  if(lnorm != "TPM" || ltrans != "log2") return(NULL)
  
  ## Run analysis
  el_time <- system.time(
    runDiffExp(data.file = dataset_file,
               result.extent = method_name,
               Rmdfunction = paste0(method, ".createRmd"),
               output.directory = results_directory,
               norm.method = "TMM",
               model = sproc,
               measurement_error = TRUE,
               extra.design.covariates = c("Localisation"),
               length.normalization = lnorm,
               data.transformation = ltrans,
               use.eBayes = eBayes,
               trend = (trend == "with_trend"),
               regularize.correlation = regCor,
               ddf.method = ddf_method)
  )
  method_info$time <- el_time[3]
  
  ## Save parameters
  res_file_name <- file.path(results_directory, paste0(dataset, "_", method_name, ".rds"))
  res_tmp <- readRDS(res_file_name)
  res_tmp@method.names <- c(res_tmp@method.names, method_info)
  saveRDS(res_tmp, file = res_file_name)
  
  ## generate the analysis code
  generateCodeHTMLs(file.path(results_directory, paste0(dataset, "_", method_name, ".rds")), results_directory)
}
