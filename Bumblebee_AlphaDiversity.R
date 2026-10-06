### ======================================================================= ###
### Title: Bumblebee Microbiome Analysis - Alpha diversity

### Description: This script outlies the analysis of alpha diversity used in 
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
  "devtools"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c("phyloseq", "microbiome")

BiocManager::install(
  bioc_packages[!bioc_packages %in% rownames(installed.packages())]
)

if (!requireNamespace("pairwiseAdonis", quietly = TRUE)) {
  devtools::install_github(
    "pmartinezarbizu/pairwiseAdonis/pairwiseAdonis"
  )
}

library(tidyverse)
library(readxl)
library(writexl)
library(ggpubr)
library(microbiome)
library(tidyr)
library(dplyr)
library(vegan)
library(ggplot2)


### ======================================================================= ###
### Read in the data:

AllMergedRarefyRDS = readRDS("AllMergedRarefy.rds")
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
### 1. ALPHA DIVERSITY ####

### 1.1 Measure RICHNESS and EVENNESS

diversity_measures <- data.frame(
  SampleID = rownames(AllMergedRarefyRDS_matrix),
  Observed = specnumber(AllMergedRarefyRDS_matrix),
  Shannon  = diversity(AllMergedRarefyRDS_matrix, index = "shannon"),
  Simpson  = diversity(AllMergedRarefyRDS_matrix, index = "simpson")   # this is Gini-Simpson (1 - D)
)

print (diversity_measures)
write.csv(diversity_measures, "alpha_diversity_measures.csv", row.names = TRUE)


### 1.2 Statistical analysis of the diversity

## 1.2.1. OBSERVED ####

### Assessing normality
diversity_measures$Treatment <- metadata$Treatment
model_Observed <- lm(Observed ~ Treatment, data = diversity_measures)
shapiro.test(residuals(model_Observed))

### ANOVA
anova_result_observed <- aov(Observed ~ Treatment, data = diversity_measures)
summary(anova_result_observed)


## 1.2.3. SHANNON ####

### Assessing normality
diversity_measures$Treatment <- metadata$Treatment
model_Shannon <- lm(Shannon ~ Treatment, data = diversity_measures)
shapiro.test(residuals(model_Shannon))

### ANOVA
anova_result_Shannon <- aov(Shannon ~ Treatment, data = diversity_measures)
summary(anova_result_Shannon)


## 1.2.4. SIMPSON ####

### Assessing normality
diversity_measures$Treatment <- metadata$Treatment
model_Simpson <- lm(Simpson ~ Treatment, data = diversity_measures)
shapiro.test(residuals(model_Simpson))

### ANOVA
anova_result_Simpson <- aov(Simpson ~ Treatment, data = diversity_measures)
summary(anova_result_Simpson)


### ======================================================================= ###
### 2. ALPHA DIVERSITY - PLOT ####

### Reorder treatments
diversity_measures$Treatment <- factor (diversity_measures$Treatment, 
                                        levels = c("Ctrl", "Par", "PFP", "NFP")
)

### Define colours
custom_colors <- c("Ctrl" = "#7CAB7DFF", "Par" = "#0057B7FF",
                   "PFP" = "#F9D662FF", "NFP" = "#FF7676FF")

### Plotting
alpha = ggplot (diversity_measures, aes(x = Treatment, y = Shannon, fill = Treatment)) +
  geom_boxplot (outlier.shape = NA, aes (fill = Treatment)) + # Remove the outliers from the plot but not from the data
  scale_fill_manual(values = custom_colors)+
  scale_color_manual(values = custom_colors) +
  geom_jitter(aes(fill = Treatment, shape = Treatment), 
              colour = "grey", size = 3, stroke = 1) +
  scale_shape_manual(values = c("Ctrl" = 21, "Par" = 22, "PFP" = 23, "NFP" = 24))+
  scale_x_discrete(labels = function(x) stringr::str_wrap(
    recode (x, "Ctrl"="Control", "Par"="Parasite","PFP"="Pyrethroid- fungicide-parasite",
            "NFP"="Neonicotinoid- fungicide-parasite"), width = 20))+
  labs(x = "Treatment",
       y = "Shannon Diversity") +
  ylim(0, 1.5)+
  theme_classic(base_size = 12) + # Size of the button label
  theme(legend.position = "none", # Remove legend on the right
        axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5), # X-axis tick labels (Treatment names)
        axis.text.y = element_text(size = 20),    # Y-axis tick labels
        axis.title.x = element_text(size = 22, margin = margin(t = 10)), # X-axis title
        axis.title.y = element_text(size = 22, margin = margin(r = 10)),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.box.margin = margin(b = 10))
alpha

### Saving the plot
ggsave("Alpha_diversity.pdf", plot = alpha, width = 12, height = 8, dpi = 300)

