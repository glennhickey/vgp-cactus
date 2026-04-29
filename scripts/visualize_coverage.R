#!/usr/bin/env Rscript

# Coverage Visualization Script for 577-way alignment
# Adapted from mammals/coverage/visualize_coverage.R
# Generates a total coverage plot colored by Lineage

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(RColorBrewer)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Args: <coverage_file> [annotations_file] [category] [subset_filter]
# category: "Total" (default), "Sex_Chroms", or "Autosomes"
# subset_filter: optional regex to filter the Orders column (e.g. "Chiroptera" for bats)
if (length(args) == 0) {
  input_file <- "vgp-577way-v1-chicken.maf.gz.cov.tsv"
  annotations_file <- "annotations.tsv"
  category <- "Total"
  cat("No input file specified. Using defaults:\n")
} else if (length(args) == 1) {
  input_file <- args[1]
  annotations_file <- "annotations.tsv"
  category <- "Total"
} else if (length(args) == 2) {
  input_file <- args[1]
  annotations_file <- args[2]
  category <- "Total"
} else {
  input_file <- args[1]
  annotations_file <- args[2]
  category <- args[3]
}

subset_filter <- if (length(args) >= 4) args[4] else NULL
subset_label <- if (!is.null(subset_filter)) tolower(str_replace_all(subset_filter, "[^A-Za-z0-9]", "_")) else NULL

cat("Coverage file:", input_file, "\n")
cat("Annotations file:", annotations_file, "\n")
cat("Category:", category, "\n")
if (!is.null(subset_filter)) cat("Subset filter (Orders):", subset_filter, "\n")

# Map category arg to contig filter value
contig_filter <- case_when(
  category == "Total" ~ "_Total_",
  category == "Sex_Chroms" ~ "_Sex_Chroms_",
  category == "Autosomes" ~ "_Autosomes_",
  TRUE ~ paste0("_", category, "_")
)
# Pretty name for titles/filenames
category_pretty <- case_when(
  category == "Total" ~ "Total",
  category == "Sex_Chroms" ~ "Sex Chromosomes",
  category == "Autosomes" ~ "Autosomes",
  TRUE ~ category
)
category_file <- tolower(str_replace_all(category_pretty, " ", "_"))
subset_file_suffix <- if (!is.null(subset_label)) paste0("_", subset_label) else ""
subset_title_suffix <- if (!is.null(subset_filter)) paste0(" [", subset_filter, "]") else ""

if (!file.exists(input_file)) stop("Error: File '", input_file, "' not found!")
if (!file.exists(annotations_file)) stop("Error: File '", annotations_file, "' not found!")

# Extract reference name from filename
# Pattern: vgp-577way-v1-REFERENCE.maf.gz.cov.tsv or .single.maf.gz.cov.tsv
filename <- basename(input_file)
ref_genome <- str_extract(filename, "v1-(.+?)\\.(single\\.)?maf", group = 1)
if (is.na(ref_genome)) {
  ref_genome <- "reference"
  cat("Warning: Could not extract reference genome from filename\n")
}
cat("Reference genome:", ref_genome, "\n")

# Read the coverage data
data <- read.table(input_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(data) <- c("contig", "max_gap", "len", "query", "aln", "ident",
                     "one_to_one_aln", "one_to_one_ident", "aln_bp", "ident_bp",
                     "one_to_one_aln_bp", "one_to_one_ident_bp")

# Filter to selected category
data <- data %>% filter(contig == contig_filter)
if (nrow(data) == 0) {
  cat("No data found for category:", category, "- skipping\n")
  quit(status = 0)
}
data$aln_gbp <- data$aln_bp / 1e9
data$ident_gbp <- data$ident_bp / 1e9
data$aln_mbp <- data$aln_bp / 1e6
data$ident_mbp <- data$ident_bp / 1e6
data$one_to_one_aln_mbp <- data$one_to_one_aln_bp / 1e6
data$one_to_one_ident_mbp <- data$one_to_one_ident_bp / 1e6

# Remove self-alignment (the query with highest coverage is the reference itself)
self_query <- data$query[which.max(data$aln_gbp)]
cat("Removing self-alignment:", self_query, "\n")
data <- data %>% filter(query != self_query)

# Read annotations file
header_line <- readLines(annotations_file, n = 1)
col_names <- strsplit(gsub("^#\\s*", "", header_line), "\t")[[1]]

annotations <- read.table(annotations_file, header = FALSE, sep = "\t",
                         stringsAsFactors = FALSE, comment.char = "#",
                         quote = "\"", fill = TRUE)
colnames(annotations) <- col_names

# Build accession -> Lineage and accession -> EnglishName lookups
annotations_subset <- annotations %>%
  select(accession, Lineage, Superorder, Orders, Family, EnglishName, ScientificName) %>%
  mutate(accession = str_trim(accession),
         Superorder = str_trim(Superorder),
         Orders = str_trim(Orders),
         Family = str_trim(Family),
         EnglishName = str_trim(EnglishName),
         ScientificName = str_trim(ScientificName)) %>%
  distinct(accession, .keep_all = TRUE)

# Split Fishes into Ray-finned, Cartilaginous, and Other Fishes
annotations_subset <- annotations_subset %>%
  mutate(Lineage = case_when(
    Lineage == "Fishes" & grepl("Actinopterygii", Superorder) ~ "Ray-finned Fishes",
    Lineage == "Fishes" & grepl("Chondrichthyes", Superorder) ~ "Cartilaginous Fishes",
    Lineage == "Fishes" ~ "Other Fishes",
    TRUE ~ Lineage
  ))

# Add manually added genomes not in annotations table
extra <- tribble(
  ~accession,            ~Lineage,    ~Superorder,   ~Orders,               ~Family,          ~EnglishName,           ~ScientificName,
  "GCA_000001405.15",    "Mammals",   NA,            "Primates",            "Hominidae",      "human (hg38)",         "Homo sapiens",
  "GCA_000001635.9",     "Mammals",   NA,            "Rodentia",            "Muridae",        "mouse (mm39)",         "Mus musculus",
  "GCF_037038585.1",     "Birds",     NA,            "Columbiformes",       "Columbidae",     "band-tailed pigeon",   "Patagioenas fasciata"
)
annotations_subset <- bind_rows(annotations_subset, extra)

# Join with coverage data by accession
data <- data %>%
  left_join(annotations_subset, by = c("query" = "accession"))

# Remove any other genomes of the same species as the reference
# (e.g. there are 2 humans and 2 mice in the alignment)
ref_species <- annotations_subset$ScientificName[annotations_subset$accession == self_query]
if (length(ref_species) == 1 && !is.na(ref_species)) {
  same_species <- annotations_subset$accession[annotations_subset$ScientificName == ref_species]
  same_species <- setdiff(same_species, self_query)  # already removed above
  if (length(same_species) > 0) {
    cat("Removing same-species genomes (", ref_species, "):", paste(same_species, collapse = ", "), "\n")
    data <- data %>% filter(!(query %in% same_species))
  }
}

# Check for missing annotations
missing_annot <- data %>% filter(is.na(Lineage)) %>% select(query) %>% distinct()
if (nrow(missing_annot) > 0) {
  cat("Warning: The following accessions have no Lineage annotation:\n")
  print(missing_annot$query)
  data$Lineage[is.na(data$Lineage)] <- "Unknown"
  data$EnglishName[is.na(data$EnglishName)] <- data$query[is.na(data$EnglishName)]
}

# Apply subset filter on Orders column if specified
if (!is.null(subset_filter)) {
  n_before <- nrow(data)
  data <- data %>% filter(grepl(subset_filter, Orders, ignore.case = TRUE))
  cat("Subset filter '", subset_filter, "': kept", nrow(data), "of", n_before, "species\n")
  if (nrow(data) == 0) stop("No species remaining after subset filter!")
}

# Use English name for display, fall back to accession
data$display_name <- ifelse(is.na(data$EnglishName) | data$EnglishName == "",
                            data$query, data$EnglishName)

# Choose color grouping: Family when subset filter is active, Lineage otherwise
if (!is.null(subset_filter)) {
  color_var <- "Family"
  data$color_group <- ifelse(is.na(data$Family) | data$Family == "", "Unknown", data$Family)
} else {
  color_var <- "Lineage"
  data$color_group <- data$Lineage
}

# Order groups by mean coverage
group_order <- data %>%
  group_by(color_group) %>%
  summarise(mean_cov = mean(aln_gbp)) %>%
  arrange(desc(mean_cov)) %>%
  pull(color_group)

cat("\n", color_var, "groups found (ordered by coverage):", paste(group_order, collapse = ", "), "\n")

# Fixed color palette for Lineage mode
all_lineage_colors <- c(
  "Mammals"              = "#66C2A5",
  "Birds"                = "#FC8D62",
  "Reptiles"             = "#8DA0CB",
  "Amphibians"           = "#E78AC3",
  "Ray-finned Fishes"    = "#A6D854",
  "Cartilaginous Fishes" = "#66BD63",
  "Other Fishes"         = "#D9EF8B",
  "Invertebrates"        = "#FFD92F",
  "Unknown"              = "#B3B3B3"
)

if (!is.null(subset_filter)) {
  # Generate a qualitative palette for Family groups
  n_groups <- length(group_order)
  if (n_groups <= 8) {
    pal <- brewer.pal(max(3, n_groups), "Set2")
  } else {
    pal <- colorRampPalette(brewer.pal(8, "Set2"))(n_groups)
  }
  group_colors <- setNames(pal[1:n_groups], group_order)
} else {
  group_colors <- all_lineage_colors[group_order]
}

cat("Color assignments:\n")
for (g in group_order) {
  cat("  ", g, ":", group_colors[g], "\n")
}

cat("\nData summary:\n")
cat("  Number of species:", nrow(data), "\n")
cat("  Coverage range:", round(min(data$aln_gbp), 2), "-", round(max(data$aln_gbp), 2), "%\n")

# Order species by color group then by coverage within each group (descending)
data <- data %>%
  mutate(color_group = factor(color_group, levels = group_order)) %>%
  arrange(color_group, desc(aln_gbp)) %>%
  mutate(display_name = factor(display_name, levels = display_name))

# Create color vector for x-axis labels
species_order <- levels(data$display_name)
label_colors <- sapply(species_order, function(sp) {
  g <- data$color_group[data$display_name == sp][1]
  group_colors[as.character(g)]
})

# Compute nice y-axis breaks that include a value near the maximum
y_max_mbp <- max(data$aln_mbp)
y_breaks <- pretty(c(0, y_max_mbp), n = 8)
# Ensure the top break is at or above the data max
if (max(y_breaks) < y_max_mbp) {
  y_breaks <- c(y_breaks, ceiling(y_max_mbp / 100) * 100)
}

# Adaptive plot dimensions and text size based on species count
n_species <- nrow(data)
plot_width <- max(8, min(20, n_species * 0.35))
x_text_size <- if (n_species <= 50) 7 else if (n_species <= 150) 5 else 3

# Create the plot
# Light bar = total aligned, darker overlay = identical matches
p <- ggplot(data, aes(x = display_name)) +
  geom_bar(aes(y = aln_mbp, fill = color_group), stat = "identity", alpha = 0.4) +
  geom_bar(aes(y = ident_mbp, fill = color_group), stat = "identity", alpha = 0.9) +
  scale_fill_manual(values = group_colors, breaks = group_order, name = color_var) +
  labs(
    title = paste0("Alignment Coverage: ", category_pretty, subset_title_suffix),
    subtitle = paste("Reference:", ref_genome),
    x = "Species",
    y = "Alignment Coverage (Mbp)\n(dark = identical, light = aligned)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = x_text_size,
                               color = label_colors),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 9)
  ) +
  scale_y_continuous(expand = c(0, 0), breaks = y_breaks,
                     limits = c(0, max(y_breaks) * 1.05))

# Save linear plot
output_file <- paste0(filename, "_coverage_", category_file, subset_file_suffix, ".png")
cat("Saving plot:", output_file, "\n")
ggsave(output_file, p, width = plot_width, height = 8, dpi = 300)

# Log-scaled version (in Mbp so values are positive on log scale)
p_log <- ggplot(data, aes(x = display_name)) +
  geom_bar(aes(y = aln_mbp, fill = color_group), stat = "identity", alpha = 0.4) +
  geom_bar(aes(y = ident_mbp, fill = color_group), stat = "identity", alpha = 0.9) +
  scale_fill_manual(values = group_colors, breaks = group_order, name = color_var) +
  labs(
    title = paste0("Alignment Coverage (log scale): ", category_pretty, subset_title_suffix),
    subtitle = paste("Reference:", ref_genome),
    x = "Species",
    y = "Alignment Coverage (Mbp, log scale)\n(dark = identical, light = aligned)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = x_text_size,
                               color = label_colors),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 9)
  ) +
  scale_y_log10(breaks = c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000),
                labels = scales::comma,
                limits = c(1, max(y_breaks) * 1.05),
                expand = c(0, 0))

output_file_log <- paste0(filename, "_coverage_", category_file, subset_file_suffix, "_log.png")
cat("Saving log plot:", output_file_log, "\n")
ggsave(output_file_log, p_log, width = plot_width, height = 8, dpi = 300)

# 1:1 version: dark bars = 1:1 aligned instead of identical
p_1to1 <- ggplot(data, aes(x = display_name)) +
  geom_bar(aes(y = aln_mbp, fill = color_group), stat = "identity", alpha = 0.4) +
  geom_bar(aes(y = one_to_one_aln_mbp, fill = color_group), stat = "identity", alpha = 0.9) +
  scale_fill_manual(values = group_colors, breaks = group_order, name = color_var) +
  labs(
    title = paste0("Alignment Coverage (1:1): ", category_pretty, subset_title_suffix),
    subtitle = paste("Reference:", ref_genome),
    x = "Species",
    y = "Alignment Coverage (Mbp)\n(dark = 1:1 aligned, light = all aligned)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = x_text_size,
                               color = label_colors),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 9)
  ) +
  scale_y_continuous(expand = c(0, 0), breaks = y_breaks,
                     limits = c(0, max(y_breaks) * 1.05))

output_file_1to1 <- paste0(filename, "_coverage_", category_file, subset_file_suffix, "_1to1.png")
cat("Saving 1:1 plot:", output_file_1to1, "\n")
ggsave(output_file_1to1, p_1to1, width = plot_width, height = 8, dpi = 300)

# 1:1 log version
p_1to1_log <- ggplot(data, aes(x = display_name)) +
  geom_bar(aes(y = aln_mbp, fill = color_group), stat = "identity", alpha = 0.4) +
  geom_bar(aes(y = one_to_one_aln_mbp, fill = color_group), stat = "identity", alpha = 0.9) +
  scale_fill_manual(values = group_colors, breaks = group_order, name = color_var) +
  labs(
    title = paste0("Alignment Coverage (1:1, log scale): ", category_pretty, subset_title_suffix),
    subtitle = paste("Reference:", ref_genome),
    x = "Species",
    y = "Alignment Coverage (Mbp, log scale)\n(dark = 1:1 aligned, light = all aligned)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = x_text_size,
                               color = label_colors),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 9)
  ) +
  scale_y_log10(breaks = c(1, 3, 10, 30, 100, 300, 1000, 3000, 10000),
                labels = scales::comma,
                limits = c(1, max(y_breaks) * 1.05),
                expand = c(0, 0))

output_file_1to1_log <- paste0(filename, "_coverage_", category_file, subset_file_suffix, "_1to1_log.png")
cat("Saving 1:1 log plot:", output_file_1to1_log, "\n")
ggsave(output_file_1to1_log, p_1to1_log, width = plot_width, height = 8, dpi = 300)

# Summary by color group
cat("\n=== Summary by", color_var, "===\n")
lineage_stats <- data %>%
  group_by(color_group) %>%
  summarise(
    n_species = n(),
    mean_coverage = mean(aln_gbp),
    median_coverage = median(aln_gbp),
    min_coverage = min(aln_gbp),
    max_coverage = max(aln_gbp)
  ) %>%
  mutate(across(where(is.numeric) & !n_species, ~ round(., 2))) %>%
  arrange(desc(mean_coverage))

print(lineage_stats)

cat("\n=== Done! ===\n")
