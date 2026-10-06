### ======================================================================= ###
### Title: Bumblebee Microbiome Analysis - Relative abundance

### Description: This script outlies the analysis of alpha diversity used in 
### "Field-relevant stressors alter the bumblebee gut microbial community".
### Date of the last modification: Aug 26th, 2026
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
  "stringr"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(readxl)
library(writexl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)


### ======================================================================= ###
### Read in the data:

AllMergedRarefy = read.csv("AllMergedRarefy.csv", row.names = 1, check.names = FALSE)
raw <- AllMergedRarefy

data_t <- as.data.frame(t(raw), stringsAsFactors = FALSE)

### Restore sample IDs
data_t$Row.names <- rownames(data_t)
rownames(data_t) <- NULL

### Move metadata columns to the front
data <- data_t %>%
  select(Row.names, Treatment, Treatment2, everything())

### Convert all taxon columns from character to numeric
taxon_cols <- setdiff(colnames(data), c("Row.names", "Treatment", "Treatment2"))
data[taxon_cols] <- lapply(data[taxon_cols], as.numeric)

### Sanity check
str(data[, 1:6])

### Pivot table
df_long <- data %>%
  pivot_longer(
    cols = -c(Row.names, Treatment, Treatment2),
    names_to = "Taxonomy",
    values_to = "Count"
  )

### Split taxonomy into rank columns
df_long <- df_long %>%
  separate(Taxonomy,
           into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
           sep = ";",
           fill = "right",
           extra = "merge")


### ======================================================================= ###
### 1. Relative abundance ####
### This plot was done at Family level

### Rename unclassified taxa
df_long$Family <- ifelse(df_long$Family == "?" | is.na(df_long$Family),
                        paste0("Unknown_", df_long$Order),
                        df_long$Family)

### Aggregate counts to family level per sample
family_counts <- df_long %>%
  group_by(Row.names, Treatment, Family) %>%
  summarise(Count = sum(Count), .groups = "drop")

### Convert to relative abundance per sample
family_rel <- family_counts %>%
  group_by(Row.names) %>%
  mutate(RelAbund = Count / sum(Count)) %>%
  ungroup()


### ======================================================================= ###
### 1.2 Bar Plot ####

### Identify the top 4 most abundant genera
top_families <- family_rel %>%
  group_by(Family) %>%
  summarise(mean_abund = mean(RelAbund)) %>%
  arrange(desc(mean_abund)) %>%
  slice_head(n = 4) %>%
  pull(Family)

print (top_families)

### Cluster everything else into "Other"
family_rel$Family_grouped <- ifelse(family_rel$Family %in% top_families, family_rel$Family, "Other")

family_order <- c(top_families, "Other")
family_rel$Family_grouped <- factor(family_rel$Family_grouped, levels = family_order)

### Organise treatments
family_rel$Treatment <- factor(family_rel$Treatment,
                              levels = c("Ctrl", "Par", "PFP", "NFP"),
                              labels = c("Control",
                                         "Parasite",
                                         "Pyrethroid+\n fungicide+\nparasite",
                                         "Neonicotinoid+\n fungicide+\nparasite"))

### Summarise mean relative abundance per Treatment x Family group
df_summary <- family_rel %>%
  group_by(Treatment, Family_grouped) %>%
  summarise(RelAbund = sum(RelAbund), .groups = "drop") %>%
  group_by(Treatment) %>%
  mutate(RelAbund = RelAbund / sum(RelAbund))

### Define colours (4 genera + Other)
family_colors <- c(setNames(c("#C4CFBFFF","#81754EFF","#CBA660FF","#86551CFF"), top_families),
                  "Other" = "darkolivegreen")

### Plot stacked bar chart
Family_Rel_Ab = ggplot(df_summary, aes(x = Treatment, y = RelAbund, fill = Family_grouped)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = family_colors) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    x = "Treatment",
    y = "Relative Abundance",
    fill = "Family",
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    legend.title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 20),
    legend.key.spacing.y = unit(0.3, "cm"),
    legend.spacing.y = unit(0.3, "cm"),
  )

Family_Rel_Ab

### Save 
ggsave("Family_Rel_Abundance.pdf", plot = Family_Rel_Ab, width = 12, height = 8, dpi = 300)




### ======================================================================= ###
### 2. Counts table ####
### This table is a summary of the counts per taxa per group (core and non-core).

### Change the name of the data frame
df_long_1 <- df_long

### Verify number of unique taxa
unique(df_long_1$Genus[grepl("Schmidhempelia", df_long_1$Genus)])
unique(df_long_1$Species[grepl("Lactobacillus", df_long_1$Species)])

### Define and group core taxa
df_long_1 <- df_long_1 %>%
  mutate(
    Bacterial_group = case_when(
      Genus %in% c("Snodgrassella", "Candidatus Schmidhempelia", "Gilliamella", "Bombiscardovia") ~ "Core",
      Genus == "Lactobacillus" & Species %in% c("Lactobacillus apis (Firm-5)", 
                                                "Lactobacillus bombi", 
                                                "Lactobacillus bombicola") ~ "Core",
      TRUE ~ "Non_core"
    )
  )

### Summarise average and total counts per Treatment, per taxon
summary_table <- df_long_1 %>%
  group_by(Kingdom, Phylum, Class, Order, Family, Genus, Species, Bacterial_group, Treatment) %>%
  summarise(
    Average_counts_per_Treatment = mean(Count),
    Total_counts_per_Treatment = sum(Count),
    .groups = "drop"
  )

### Review table
head(summary_table)

### Save the table
write.csv(summary_table, "Core_NonCore_table_fromRarefy.csv", row.names = FALSE)

### Count unique taxa count per core and non-core

taxa_counts <- summary_table %>%
  distinct(Kingdom, Phylum, Class, Order, Family, Genus, Species, Bacterial_group) %>%
  count(Bacterial_group, name = "n_unique_taxa")

print(taxa_counts)
