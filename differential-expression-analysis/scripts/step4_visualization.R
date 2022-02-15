rm(list=ls())
library(here)
# devtools::install_github("pbastide/phylocompcodeR")
library(phylocompcodeR)
library(stringr)

################################################################################
## File management
################################################################################
condName <- "Hyperaccu" # change here Hyperaccu or Tolerance
# condName <- "Tolerance"

## File
datestamp_day_real_nickel <- "2020-09-25"
results_directory <- paste0(datestamp_day_real_nickel, "_results_nickel_", condName)

dir.create(here("figures"))

################################################################################
## Data
################################################################################
dataset <- "nickel_cpd"
dataset_file <- here(file.path(results_directory, paste0(dataset, ".rds")))
cpd <- readRDS(file = dataset_file)
tree <- cpd@info.parameters$tree

################################################################################
## Analysis
################################################################################
method_name <- "phylolm.RPKM.log2.OU"

res <- readRDS(here(file.path(results_directory, paste0(dataset, "_", method_name, ".rds"))))
res_table <- res@result.table

## Select COGs
res_select <- rownames(res_table[res_table$adjpvalue <= 0.01 & abs(res_table$logFC) >= 1.5, ])

# RPKM
nf <- edgeR::calcNormFactors(cpd@count.matrix, method = 'TMM')
lib.size <- colSums(cpd@count.matrix) * nf
data.norm <- sweep((cpd@count.matrix + 0.5) / cpd@length.matrix, 2, lib.size + 1, '/')
data.norm <- data.norm * 1e9
rownames(data.norm) <- rownames(cpd@count.matrix)
## log2
data.trans <- log2(data.norm)
# data.trans <- data.norm

## Data for selected
data_select <- data.norm[res_select, ]
data_select <- as.data.frame(t(data_select))

data_select_trans <- data.trans[res_select, ]
data_select_trans <- as.data.frame(t(data_select_trans))

data_select_norm <- scale(data_select)
data_select_norm <- as.data.frame(data_select_norm)

data_select_trans_norm <- scale(data_select_trans)
data_select_trans_norm <- as.data.frame(data_select_trans_norm)

################################################################################
## Re-Do Analysis to get Coefficients
################################################################################
design_formula <- as.formula(" ~  Localisation + condition")
design_data <- res@sample.annotations[, c('Localisation', 'condition'), drop = FALSE]
design_data <- data.frame(apply(design_data, 2, as.factor))
design <- model.matrix(design_formula, design_data)
model <- "OUfixedRoot"

extract_residuals <- function(dat) {
  data_reg <- cbind(data.frame(expr = dat), design_data)
  levels(data_reg$condition) <- c(1, 2)
  res_phylolm <- phylolm::phylolm(paste("expr", paste(as.character(design_formula), collapse = "")), 
                                  data = data_reg, phy = tree, model = model, measurement_error = TRUE)
  preds <- names(res_phylolm$coefficients) != "condition2"
  resids <- res_phylolm$y - res_phylolm$X[, preds] %*% res_phylolm$coefficients[preds]
  return(resids)
}

data_select_trans_resids <- apply(data_select_trans, 2, extract_residuals)
rownames(data_select_trans_resids) <- rownames(data_select_trans)
data_select_trans_resids <- as.data.frame(data_select_trans_resids)

data_select_trans_resids_norm <- scale(data_select_trans_resids)
data_select_trans_resids_norm <- as.data.frame(data_select_trans_resids_norm)

data_select_resids <- 2^data_select_trans_resids
data_select_resids_norm <- scale(data_select_resids)
data_select_resids_norm <- as.data.frame(data_select_resids_norm)


################################################################################
## Plot - PCA
################################################################################
library(DESeq2)
library(ggplot2)
DESeq2.length.ds <- DESeq2::DESeqDataSetFromMatrix(countData = res@count.matrix,
                                                   colData = res@sample.annotations[, c("Famille", "condition")],
                                                   design = ~ Famille + condition)
## Size Factors
DESeq2.length.ds <-  estimateSizeFactors(DESeq2.length.ds)
size_fac <- sizeFactors(DESeq2.length.ds)
mat_size_fac <- matrix(size_fac, ncol = length(size_fac), nrow = nrow(res@count.matrix), byrow = T)
## Extra factors
extraNormFactor <- res@length.matrix
normFactors <- (mat_size_fac * extraNormFactor) / exp(rowMeans(log(mat_size_fac * extraNormFactor)))
normalizationFactors(DESeq2.length.ds) <- as.matrix(normFactors)

colData(DESeq2.length.ds)
rld <- rlog(DESeq2.length.ds)

##  rename oui non 
library(plyr)
rld@colData@listData$condition <- mapvalues(rld@colData@listData$condition, from = c("oui", "non"), to = c("hyper", "non-hyper"))

## plot PCA 
pcaData <- plotPCA(rld, intgroup=c("condition","Famille"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition, fill=condition, shape=Famille)) +
  geom_point(size=3, alpha = 0.8) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_color_viridis_d(end = 0.7, name = "Hyperaccu.",
                        breaks = c("hyper", "non-hyper"),
                        labels = c("Yes", "No")) +
  scale_fill_viridis_d(end = 0.7,
                       breaks = c("hyper", "non-hyper"),
                       labels = c("Yes", "No"),
                       guide = "none") + 
  #geom_text(aes(label=c(matrix(unlist(strsplit(rownames(colData),"_")),ncol=2,byrow = T)[,1])),hjust=.5, vjust=-.8, size=3) +
  #geom_density2d(alpha=.5) +
  # scale_shape(name = "Familly", solid = TRUE) +
  scale_shape_manual(name = "Familly", values = c(21:25, 8)) +
  guides(shape = guide_legend(override.aes = list(fill = "black"))) + 
  theme_bw() +
  coord_fixed()
# colorblindr::cvd_grid()
ggsave(here::here("figures/PCA.pdf"), width = 5, height = 4)


################################################################################
## Plot - Full Matrix
################################################################################
library(ggplot2)
library(ggtree)
library(tidytree)
library(tibble)
library(tidyr)
library(dplyr)
library(aplot)
library(cowplot)

conds <- data.frame(label = rownames(cpd@sample.annotations),
                    condition = as.factor(cpd@sample.annotations$condition),
                    localisation = as.factor(cpd@sample.annotations$Localisation))
gt <- ggtree(tree) + ylim(c(0, 52))
# gt <- gt + geom_tiplab(size = 1.9, linesize = .5, offset = 10)
gt <- gt %<+% conds + geom_tippoint(aes(shape = localisation, color = condition), size = 2, position = position_nudge(2))
gt <- revts(gt + theme_tree2() + labs(caption = "myr ago"))
gt <- gt + scale_color_viridis_d(name = "Hyperaccu.",
                                 end = 0.7,
                                 # values = c("red", "blue"),
                                 breaks = c("oui", "non"),
                                 labels = c("Yes", "No"),
                                 guide = guide_legend(title.position = "top",
                                                      title.hjust = 0.5,
                                                      order = 1))
gt <- gt + scale_shape_manual(name = "Location",
                              breaks = c("Afrique du sud", "Cuba", "France", "Nouvelle Caledonie"),
                              labels = c("South Africa", "Cuba", "France", "New Caledonia"),
                              values = 15:18)
gt <- gt + theme(text = element_text(size = 7))
# gt <- gt + geom_text(aes(label=node))
gt <- gt %>% rotate(52) %>% rotate(57) %>% rotate(73)
# gt

## Scaled
data_select_long <- data_select_norm %>% rownames_to_column("label") %>%
  pivot_longer(cols = !"label", names_to = "OG") 
data_select_long$OG <- factor(data_select_long$OG, levels = unique(sort(as.numeric(data_select_long$OG))))

p2 <- ggplot(data_select_long,
             aes(x = OG, y = label)) +
  geom_tile(aes(fill = value)) +
  scale_fill_viridis_c(option = 'magma', direction = -1, name = 'RPKM (scaled)',
                       guide = guide_colorbar(#title.position = "top",
                         title.hjust = 0.5,
                         order = 2)) +
  scale_x_discrete(position = "top") +
  # theme_tree2() +
  theme(text = element_text(size = 6),
        # legend.position = "bottom",
        axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5),
        axis.line.x.top = element_blank(),
        # axis.ticks.x.top = element_blank(),
        axis.title.y = element_blank())
y_lim_groups <- c(2, 5, 8, 11, 14, 17, 21, 23, 25, 30, 33, 36, 39, 42, 45)
coord_groups <- data.frame(x1 = 0.5, x2 = 31.5, y1 = y_lim_groups + 0.5, y2 = y_lim_groups + 0.5)
p2 <- p2 + geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2), data = coord_groups)
# p2

p3 <- p2 %>% insert_left(gt, width = 0.5)
p3
ggsave(here::here("figures/selected_OG.pdf"), plot = p3, width = 4.8, height = 4.8)

################################################################################
## Plot - Full Matrix - Residuals
################################################################################
gt <- ggtree(tree) + ylim(c(0, 52))
# gt <- gt + geom_tiplab(size = 1.9, linesize = .5, offset = 10)
gt <- gt %<+% conds + geom_tippoint(aes(color = condition), size = 2, position = position_nudge(2))
gt <- revts(gt + theme_tree2() + labs(caption = "myr ago"))
gt <- gt + scale_color_viridis_d(name = "Hyperaccu.",
                                 end = 0.7,
                                 # values = c("red", "blue"),
                                 breaks = c("oui", "non"),
                                 labels = c("Yes", "No"),
                                 guide = guide_legend(title.position = "top",
                                                      title.hjust = 0.5,
                                                      order = 1))
# gt <- gt + scale_shape_manual(name = "Location",
#                               breaks = c("Afrique du sud", "Cuba", "France", "Nouvelle Caledonie"),
#                               labels = c("South Africa", "Cuba", "France", "New Caledonia"),
#                               values = 15:18)
gt <- gt + theme(text = element_text(size = 7))
# gt <- gt + geom_text(aes(label=node))
gt <- gt %>% rotate(52) %>% rotate(57) %>% rotate(73)

data_select_long <- data_select_resids_norm %>% rownames_to_column("label") %>%
  pivot_longer(cols = !"label", names_to = "OG") 
data_select_long$OG <- factor(data_select_long$OG, levels = unique(sort(as.numeric(data_select_long$OG))))

p2 <- ggplot(data_select_long,
             aes(x = OG, y = label)) +
  geom_tile(aes(fill = value)) +
  scale_fill_viridis_c(option = 'magma', direction = -1, name = 'RPKM \n(scaled residuals)',
                       guide = guide_colorbar(#title.position = "top",
                         title.hjust = 0,
                         order = 2)) +
  scale_x_discrete(position = "top") +
  # theme_tree2() +
  theme(text = element_text(size = 6),
        # legend.position = "bottom",
        axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5),
        axis.line.x.top = element_blank(),
        # axis.ticks.x.top = element_blank(),
        axis.title.y = element_blank())
y_lim_groups <- c(2, 5, 8, 11, 14, 17, 21, 23, 25, 30, 33, 36, 39, 42, 45)
coord_groups <- data.frame(x1 = 0.5, x2 = 31.5, y1 = y_lim_groups + 0.5, y2 = y_lim_groups + 0.5)
p2 <- p2 + geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2), data = coord_groups)
# p2

p3 <- p2 %>% insert_left(gt, width = 0.5)
p3
ggsave(here::here("figures/selected_OG_residuals.pdf"), plot = p3, width = 4.8, height = 4.8)

################################################################################
## Plot - selected OGs
################################################################################
get_ticks <- function(llim) {
  xmin <- llim[1]
  xmax <- llim[2]
  # mmin <- floor(xmin)
  mmax <- floor((xmax - 10) / 10) * 10
  # return(c(0, (0 + mmax) / 2, mmax))
  return(c(0, mmax))
}

species <- as_tibble(tree)
species <- species[!is.na(species$label), c("label")]

## Data
conds <- data.frame(label = rownames(cpd@sample.annotations),
                    condition = as.factor(cpd@sample.annotations$condition),
                    localisation = as.factor(cpd@sample.annotations$Localisation),
                    familly = as.factor(cpd@sample.annotations$Famille))
fams <- data.frame(familly = unique(conds$familly),
                   mrca = sapply(unique(conds$familly),
                                 function(x)  ape::getMRCA(tree, as.character(conds[conds$familly == x, "label"]))))

tree <- ape::makeNodeLabel(tree)
for (i in 1:nrow(fams)) {
  tree$node.label <- sub(paste0("Node", fams[i, "mrca"] - length(tree$tip.label), "$"),
                         fams[i, "familly"], tree$node.label)
}

## Tree
gt <- ggtree(tree) + ylim(c(0, 48))
# gt <- gt + geom_tiplab(size = 1.9, linesize = .5, offset = 10)
gt <- gt %<+% conds + geom_tippoint(aes(color = condition), size = 1, position = position_nudge(2))
gt <- revts(gt + theme_tree2() + labs(caption = "myr ago"))
gt <- gt + scale_color_viridis_d(name = "Hyperaccu.",
                                 # values = c("red", "blue"),
                                 end = 0.7,
                                 breaks = c("oui", "non"),
                                 labels = c("Yes", "No"),
                                 guide = FALSE)
# gt <- gt + scale_shape_manual(name = "Location",
#                               breaks = c("Afrique du sud", "Cuba", "France", "Nouvelle Caledonie"),
#                               labels = c("South Africa", "Cuba", "France", "New Caledonia"),
#                               values = 15:18,
#                               guide = FALSE)
# gt <- gt + geom_text(aes(label=node))
gt <- gt %>% rotate(52) %>% rotate(57) %>% rotate(73)
gt <- gt %<+% fams + geom_text2(aes(label = label, subset = label %in% fams$familly),
                                color = 'steelblue',
                                hjust = 1.1, vjust = -1, size = 1.8)
# gt <- gt + geom_hilight(data = fams, mapping  = aes(node = mrca, fill = familly))
# gt <- gt + scale_fill_discrete(guide = FALSE)
# gt <- gt + geom_hilight(node = 85)
gt <- gt + theme(text = element_text(size = 10))
# gt

## Data
d <- filter(gt, isTip) %>% select(c(label, y))
data_traits <-  data_select
colnames(data_traits) <- paste0("OG_", colnames(data_traits))
data_traits <-  rownames_to_column(data_traits, "label")
data_traits <- left_join(data_traits, d, by='label')
data_traits <- left_join(data_traits, conds, by='label')
data_traits$label <- gsub("_", " ", data_traits$label)
## replace ech with s in label
data_traits$label <- str_replace(data_traits$label,"ech","s")

p1 <- ggplot(data_traits, aes(y, OG_13722)) +
  ggtitle("13722") +
  geom_point(aes(color = condition), size = 0.5) +
  coord_flip() + theme_tree2() +
  theme(legend.position = 'none',
        plot.title = element_text(size = 10),
        text = element_text(size = 10)) +
  scale_color_viridis_d(breaks = c("oui", "non"), labels = c("Yes", "No"), end = 0.7, guide = FALSE) +
  scale_y_continuous(breaks = get_ticks, expand = c(0, 0), limits = c(-10, max(data_traits$OG_13722) + 10)) +
  labs(caption = "RPKM") +
  ylim2(gt)
p2 <- ggplot(data_traits, aes(y, OG_4147)) +
  ggtitle("4147") +
  geom_point(aes(color = condition), size = 0.5) +
  coord_flip() + theme_tree2() +
  theme(legend.position = 'none',
        plot.title = element_text(size = 10),
        text = element_text(size = 10)) +
  scale_color_viridis_d(breaks = c("oui", "non"), labels = c("Yes", "No"), end = 0.7, guide = FALSE) +
  scale_y_continuous(breaks = get_ticks, expand = c(0, 0), limits = c(-10, max(data_traits$OG_4147) + 10)) +
  labs(caption = "RPKM") +
  ylim2(gt)
p3 <- ggplot(data_traits, aes(y, OG_7137)) +
  ggtitle("7137") +
  geom_point(aes(color = condition), size = 0.5) +
  coord_flip() + theme_tree2() +
  theme(legend.position = 'none',
        plot.title = element_text(size = 10),
        text = element_text(size = 10)) +
  scale_color_viridis_d(breaks = c("oui", "non"), labels = c("Yes", "No"), end = 0.7, guide = FALSE) +
  scale_y_continuous(breaks = get_ticks, expand = c(0, 0), limits = c(-10, max(data_traits$OG_7137) + 10)) +
  labs(caption = "RPKM") +
  ylim2(gt)

pLabel <- ggplot(data_traits, aes(x = y, y = 0, color = condition)) + geom_text(aes(label = label, hjust = 0), size = 2.1) +
  ylim(0, 0.1) +
  # ggtitle("R3") +
  coord_flip() + theme_tree2() +
  theme(legend.position='none') +
  scale_color_viridis_d(breaks = c("oui", "non"), labels = c("Yes", "No"), end = 0.7, guide = FALSE) + 
  ylim2(gt) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank())

p5 <- plot_grid(gt, p1, p2, p3, pLabel, ncol = 5, align = 'h', rel_widths = c(1, 1, 1, 1, 0.8))
p5
ggsave(here::here("figures/selected_OG_dots_3.pdf"), plot = p5, width = 4.8, height = 4.8)

p6 <- plot_grid(gt, p2, p3, pLabel, ncol = 4, align = 'h', rel_widths = c(1, 1.5, 1.5, 0.8))
p6
ggsave(here::here("figures/selected_OG_dots_2.pdf"), plot = p6, width = 4.8, height = 4.8)

