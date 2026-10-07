### ======================================================================= ###
### Title: Bumblebee Microbiome Analysis - Beta diversity

### Description: This script outlies the analysis of beta diversity used in 
### "Field-relevant stressors alter the bumblebee gut microbial community".
### Date of the last modification: Aug 25th, 2026
### ======================================================================= ###

### Housekeeping

rm(list=ls())

graphics.off()


### Installing packages and loading libraries

cran_packages <- c(
  "readxl",
  "writexl",
  "tidyverse",
  "ggpubr",
  "devtools",
  "phyloseq"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c("microbiome")

BiocManager::install(
  bioc_packages[!bioc_packages %in% rownames(installed.packages())]
)

library(tidyverse)
library(readxl)
library(writexl)
library(ggpubr)
library(microbiome)
library(vegan)
library(phyloseq)

### ======================================================================= ###
### Read in the data:

AllMergedRarefyRDS = readRDS(file.choose())

class(AllMergedRarefyRDS)
str(AllMergedRarefyRDS)

### Prepare the table:
metadata_cols = c("Row.names","Treatment", "Treatment2")
species_cols <- setdiff(colnames(AllMergedRarefyRDS), metadata_cols)

### Check for NAs or blanks before converting into numbers
sapply(AllMergedRarefyRDS[species_cols], function(x) any(x == "" | x == " ", na.rm = TRUE))

AllMergedRarefyRDS[species_cols] = lapply(AllMergedRarefyRDS[species_cols], as.numeric)

### Sanity check of the attributes
sapply(AllMergedRarefyRDS[species_cols], class)
sum(is.na(AllMergedRarefyRDS[species_cols]))
sum(AllMergedRarefyRDS[species_cols] == 0, na.rm = TRUE)

### Review the structure of the table
str(AllMergedRarefyRDS)

### Associate the metadata
metadata = AllMergedRarefyRDS[c("Row.names","Treatment", "Treatment2")]

### Convert RDS table into a numerical matrix
AllMergedRarefyRDS_matrix = as.matrix(AllMergedRarefyRDS[, species_cols])
str(AllMergedRarefyRDS_matrix)



### ======================================================================= ###
### 1. BETA DIVERSITY ####

### 1.1 Bray-Curtis distances ####
BrayDistance_metaMDS_T = metaMDSdist(AllMergedRarefyRDS_matrix, 
                                     distance = "bray", 
                                     autotransform = TRUE)
print (BrayDistance_metaMDS_T)

### Convert data into a matrix readable in Excel
BrayDistance_metaMDS_T_matrix <- as.matrix(BrayDistance_metaMDS_T)

### Save the table
write.csv(BrayDistance_metaMDS_T_matrix,"BrayDistance_metaMDS_T_matrix.csv", row.names = T) 


### 1.2 Bray - PERMANOVA ####
set.seed(111)
permanova_metaMDS_T <- adonis2(BrayDistance_metaMDS_T ~ Treatment, 
                            data = metadata, 
                            permutations = 999) # FDR correction is not possible in adonis2.

print (permanova_metaMDS_T)

### Save the table
capture.output (permanova_metaMDS_T, file = "permanova_metaMDS_T.txt")


### 1.2 Bray - Pairwise PERMANOVA ####

### Get all pairwise combinations of treatment groups
groups   <- metadata$Treatment
combos   <- combn(unique(groups), 2, simplify = FALSE)

### Loop over each pair and run adonis2
results_list <- lapply(combos, function(pair) {
  
  ### Subset to only the two groups in this pair
  keep     <- groups %in% pair
  sub_dist <- as.dist(as.matrix(BrayDistance_metaMDS_T)[keep, keep])
  sub_meta <- metadata[keep, , drop = FALSE]
  
  ### Run PERMANOVA on the subset
  set.seed(111)
  res <- adonis2(
    sub_dist ~ Treatment,
    data         = sub_meta,
    permutations = 999
  )
  
  ### Extract key statistics into a one-row data frame
  data.frame(
    group1  = pair[1],
    group2  = pair[2],
    R2      = res$R2[1],
    F_value = res$F[1],
    p_value = res$`Pr(>F)`[1]
  )
})

### Combine all rows

pairwise_df <- do.call(rbind, results_list)

### Apply BH correction across all pairwise p-values

pairwise_df$p_adj_BH <- p.adjust(pairwise_df$p_value, method = "BH")

print(pairwise_df)

### Save the table

write.csv(pairwise_df, "pairwise_permanova_BH.csv", row.names = FALSE)
capture.output (pairwise_df, file = "pairwise_permanova.txt")


### ======================================================================= ###
### 2. DISPERSION ####

### 2.1 Dispersion using the Bray-Curtis distance matrix ####
Betadispersion = betadisper(BrayDistance_metaMDS_T, sample_data(AllMergedRarefyRDS)$Treatment)
print (Betadispersion)


### ANOVA
anova_disp = anova(Betadispersion)

### Save the result
capture.output (anova_disp, file = "anova_dispersion.txt")

### Observe differences between groups
plot(Betadispersion)
boxplot(Betadispersion)

### 2.2 Pairwise of the dispersion ####
pairwise_disp = TukeyHSD(Betadispersion)
capture.output (pairwise_disp, file = "pairwise_dispersion.txt")



### ======================================================================= ###
### 3. BETA DIVERSITY - PLOTS ####


### Load libraries for the PCoA and UMAP

library(ape)
library(patchwork)
library(umap)
library(ggforce)  # for stat_ellipse or encircling groups
library(RColorBrewer)


### 3.1 PCoA ####

### Read the distance matrix
BrayDistance_metaMDS_T = read.csv("BrayDistance_metaMDS_T_matrix.csv", row.names = 1, check.names = FALSE)
dist_mat_PCoA = BrayDistance_metaMDS_T

### Run PCoA 
pcoa_result <- pcoa(dist_mat_PCoA)
print (pcoa_result)

### Extract coordinates
pcoa_coords <- as.data.frame(pcoa_result$vectors[, 1:2])
colnames(pcoa_coords) <- c("PCoA1", "PCoA2")

### Add treatment metadata to the PCoA
pcoa_coords$Treatment <- factor(AllMergedRarefyRDS$Treatment,
                                levels = c("Ctrl", "Par",
                                           "PFP", "NFP"),
                                labels = c("Control",
                                           "Parasite",
                                           "Pyrethroid+fungicide+\nparasite",
                                           "Neonicotinoid+fungicide+\nparasite"))

eig_percent <- round(pcoa_result$values$Relative_eig[1:2] * 100, 2)

### Define colours 
custom_colors <- c("Control" = "#7CAB7DFF", "Parasite" = "#0057B7FF",
                   "Pyrethroid+fungicide+\nparasite" = "#F9D662FF", "Neonicotinoid+fungicide+\nparasite" = "#FF7676FF")

### Plot each pair in a single plot

PCoA = ggplot(pcoa_coords, aes(x = PCoA1, y = PCoA2, color = Treatment, fill = Treatment, shape = Treatment)) +
  geom_point(size = 6, alpha = 0.8)+
  stat_ellipse(geom = "polygon", 
               alpha = 0.15,
               level = 0.95,
               fill = NA,
               linewidth = 0.8) +
  scale_color_manual(values = custom_colors,
                     labels = c("Control",
                                "Parasite",
                                "Pyrethroid+fungicide+\nparasite",
                                "Neonicotinoid+fungicide+\nparasite"))+
  scale_fill_manual(values = custom_colors,
                    labels = c("Control",
                               "Parasite",
                               "Pyrethroid+fungicide+\nparasite",
                               "Neonicotinoid+fungicide+\nparasite")) +
  scale_shape_manual(values = c("Control" = 21, 
                                "Parasite" = 22,
                                "Pyrethroid+fungicide+\nparasite" = 23,
                                "Neonicotinoid+fungicide+\nparasite" = 24),
                     labels = c("Control",
                                "Parasite",
                                "Pyrethroid+fungicide+\nparasite",
                                "Neonicotinoid+fungicide+\nparasite"))+
  labs(
    x = paste0("PCoA1 (", eig_percent[1], "%)"),
    y = paste0("PCoA2 (", eig_percent[2], "%)")
  ) +
  theme_classic() +
  theme(axis.title = element_text(size = 22),   # axis titles (PCoA1/PCoA2 labels)
        axis.text  = element_text(size = 24),   # axis tick numbers
        axis.text.x = element_text(margin = margin(t = 12)),  # space above x-axis text
        axis.text.y = element_text(margin = margin(r = 10)),   # space right of y-axis text
        legend.title = element_text(size = 20, face = "bold"), # "Treatment" legend title
        legend.text  = element_text(size = 20),
        legend.key.spacing.y = unit(0.3, "cm"),
        legend.spacing.y = unit(0.3, "cm"),
        legend.box.margin = margin(b = 10)) # legend entry labels


PCoA

### Save the plot
ggsave("PCoA_all.pdf", plot = PCoA, width = 14, height = 12, dpi = 300)



### ======================================================================= ###
### 3.2. UMAP (Supplementary material: figure S)####

### Read distance table
BrayDistance_metaMDS_T_matrix = read.csv("BrayDistance_metaMDS_T_matrix.csv", header = TRUE)

### Prepare the table
rownames(BrayDistance_metaMDS_T_matrix) <- BrayDistance_metaMDS_T_matrix[, 1]
BrayDistance_metaMDS_T_matrix <- BrayDistance_metaMDS_T_matrix[, -1]

### Convert into matrix
dist_mat_UMAP = as.matrix (BrayDistance_metaMDS_T_matrix)

### Define the pairwise comparisons
comparisons <- list(
  c("Ctrl", "Par"),
  c("Ctrl", "PFP"),
  c("Ctrl", "NFP"),
  c("Par",  "PFP"),
  c("Par",  "NFP"),
  c("PFP",  "NFP")
)

### Read metadata
sample_meta <- metadata
sample_meta <- sample_meta[match(rownames(dist_mat_UMAP), sample_meta$Row.names), ]

### Run UMAP for each pair and collect results
umap_list <- lapply(comparisons, function(pair) {
  
  keep <- sample_meta$Treatment %in% pair
  sub_dist <- as.matrix(dist_mat_UMAP)[keep, keep]
  sub_meta <- sample_meta[keep, ]
  
  umap_config <- umap.defaults
  umap_config$input <- "dist"
  
  umap_result <- umap(sub_dist, config = umap_config)
  
  data.frame(
    UMAP1 = umap_result$layout[, 1],
    UMAP2 = umap_result$layout[, 2],
    Treatment = sub_meta$Treatment,
    Comparison = paste(pair, collapse = " vs ")
  )
})

### Combine into one dataframe
umap_all <- bind_rows(umap_list)

### Set comparison order
umap_all$Comparison <- factor(umap_all$Comparison,
                              levels = c("Ctrl vs Par", "Ctrl vs PFP", "Ctrl vs NFP",
                                         "Par vs PFP", "Par vs NFP", "PFP vs NFP"))
### Reorganise treatment labels
umap_all$Treatment <- factor(umap_all$Treatment,
                             levels = c("Ctrl", "Par", "PFP", "NFP"),
                             labels = c("Control", 
                                        "Parasite", 
                                        "Pyrethroid+fungicide+\nparasite", 
                                        "Neonicotinoid+fungicide+\nparasite"))

### Plotting

UMAP <- ggplot(umap_all, aes(x = UMAP1, y = UMAP2, fill = Treatment, shape = Treatment)) +
  geom_point(size = 6, alpha = 0.8, color = "white") +
  stat_ellipse(aes (fill = NULL, color = Treatment),
               geom = "polygon", 
               alpha = 0.15,
               level = 0.95,
               fill = NA,
               linewidth = 0.8) +
  scale_fill_manual(values = c("Control" = "#7CAB7DFF",
                               "Parasite" = "#0057B7FF",
                               "Pyrethroid+fungicide+\nparasite" = "#F9D662FF",
                               "Neonicotinoid+fungicide+\nparasite" = "#FF7676FF")) +
  scale_color_manual(values = c("Control" = "#7CAB7DFF",
                                "Parasite" = "#0057B7FF",
                                "Pyrethroid+fungicide+\nparasite" = "#F9D662FF",
                                "Neonicotinoid+fungicide+\nparasite" = "#FF7676FF"),
                     guide = "none") +   # avoid a duplicate legend
  scale_shape_manual(values = c("Control" = 21, 
                                "Parasite" = 22,
                                "Pyrethroid+fungicide+\nparasite" = 23,
                                "Neonicotinoid+fungicide+\nparasite" = 24))+
  facet_wrap(~ Comparison, ncol = 3, scales = "free") +
  coord_cartesian(xlim = c(-4, 4), ylim = c(-4, 4)) +
  theme_classic() +
  theme(
    axis.text  = element_text(size = 14),
    axis.title = element_text(size = 20),
    strip.text = element_blank(), # no titles in every plot
    legend.title = element_text(face = "bold", size = 22),
    legend.text  = element_text(size = 20),
    axis.text.x = element_text(margin = margin(t = 12)),
  ) +
  labs(x = "UMAP1", y = "UMAP2")

UMAP

### Saving the plot
ggsave("UMAP_pairwise_all.pdf", plot = UMAP, width = 14, height = 12, dpi = 300)

