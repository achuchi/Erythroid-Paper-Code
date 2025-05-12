library(dplyr)
library(Seurat)
library(SeuratObject)
library(patchwork)
library(monocle3)
library(biomaRt)
library(SeuratWrappers)
library(ggplot2)

#------------------------------ DATA CLEANING


# Loading full expression matrix -- path is up to you

full_matrix <- read.table("C:/Users/aidan/Downloads/expressionmatrix_all.txt", 
                          header = TRUE, 
                          row.names = 1)


# Function to clean gene names of different symbols that don't play nice w/ Seurat
clean_gene_names <- function(names) {
  cleaned <- gsub("_.*|\\-.*", "", names)
  cleaned <- trimws(cleaned)
  cleaned[cleaned == ""] <- NA
  return(cleaned)
}

# Applying clean_gene_names function and removing NA values 
rownames(full_matrix) <- clean_gene_names(rownames(full_matrix))
valid_rows <- !is.na(rownames(full_matrix)) & rownames(full_matrix) != ""
full_matrix <- full_matrix[valid_rows, ]
full_matrix <- full_matrix[!duplicated(rownames(full_matrix)), ]

# Convert to sparse matrix for memory purposes 
count_matrix <- as(as.matrix(full_matrix), "dgCMatrix")

#--------------------------------- SEURAT

# Create Seurat object
full_seurat <- CreateSeuratObject(counts = count_matrix)

# Create a more sophisticated metadata assignment based on cell name patterns
cell_info <- colnames(full_seurat)

# Assign cell type
full_seurat$cell_type <- case_when(
  grepl("CFUE", cell_info) ~ "CFUE",
  grepl("BFUE", cell_info) ~ "BFUE",
  grepl("E25", cell_info) ~ "E25",
  grepl("E50", cell_info) ~ "E50",
  TRUE ~ "Other"
)

# Assign day information (combining all possible day labels into single label)
full_seurat$day <- case_when(
  grepl("D1|.D1|Day1", cell_info) ~ "Day1",
  grepl("D2|.D2|Day2", cell_info) ~ "Day2",
  grepl("D3|.D3|Day3", cell_info) ~ "Day3",
  TRUE ~ NA_character_
)

# Assign treatment information
full_seurat$treatment <- case_when(
  grepl("DEX|Dex", cell_info) ~ "DEX",
  grepl("NoDex", cell_info) ~ "NoDex",
  TRUE ~ "Control"
)

# Create a combined identity column
full_seurat$identity <- ifelse(is.na(full_seurat$day),
                               full_seurat$cell_type,
                               paste(full_seurat$treatment, full_seurat$cell_type, "Day", full_seurat$day))

# Check metadata assignment -- This is just a diagnostic check I do down the pipeline
print("Cell type distribution:")
print(table(full_seurat$cell_type))
print("Treatment distribution:")
print(table(full_seurat$treatment))
print("Day distribution:")
print(table(full_seurat$day, useNA = "ifany"))

# Fix empty or duplicate rownames
empty_names <- which(rownames(full_seurat) == "")
if(length(empty_names) > 0) {
  rownames(full_seurat)[empty_names] <- paste0("gene_", empty_names)
}

dup_names <- which(duplicated(rownames(full_seurat)))
if(length(dup_names) > 0) {
  rownames(full_seurat)[dup_names] <- paste0(rownames(full_seurat)[dup_names], "_dup", seq_along(dup_names))
}

# Standard SEURAT preprocessing
full_seurat <- NormalizeData(full_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
full_seurat <- FindVariableFeatures(full_seurat, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(full_seurat)
full_seurat <- ScaleData(full_seurat, features = all.genes)
full_seurat <- RunPCA(full_seurat, features = VariableFeatures(object = full_seurat))

#--------------------------BiomaRt

# Get mouse Ensembl database with simpler connection
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Get the gene IDs from your Seurat object
clean_ids <- rownames(full_seurat)

# Query gene name database with biomaRt
gene_mapping <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filters = "ensembl_gene_id",
  values = clean_ids,
  mart = mart
)

# Create a named vector for mapping
id_to_symbol <- setNames(gene_mapping$mgi_symbol, gene_mapping$ensembl_gene_id)

# Replace gene IDs with symbols where available
new_names <- clean_ids
names_to_change <- clean_ids %in% names(id_to_symbol)
new_names[names_to_change] <- id_to_symbol[clean_ids[names_to_change]]
rownames(full_seurat) <- new_names


# Standard cell cycle scoring with mouse genes
s.genes.mouse <- c("Mcm5", "Pcna", "Tyms", "Fen1", "Mcm2", "Mcm4", "Rrm1", "Ung", "Gins2", "Mcm6")
g2m.genes.mouse <- c("Cdk1", "Top2a", "Mki67", "Cenpf", "Aurkb", "Bub1", "Ccnb2", "Birc5")

# Check which genes exist in the dataset
s.genes.exist <- s.genes.mouse[s.genes.mouse %in% rownames(full_seurat)]
g2m.genes.exist <- g2m.genes.mouse[g2m.genes.mouse %in% rownames(full_seurat)]


# Run cell cycle scoring with known genes
full_seurat <- CellCycleScoring(full_seurat,
                                s.features = s.genes.exist,
                                g2m.features = g2m.genes.exist,
                                set.ident = TRUE)

# Run standard UMAP
full_seurat <- RunUMAP(full_seurat, dims = 1:30)


#------------------------------------HARMONY BATCH CORRECTION

library(harmony)

full_seurat$treatment_clean <- ifelse(is.na(full_seurat$treatment), "Unknown", full_seurat$treatment)
full_seurat$day_clean <- ifelse(is.na(full_seurat$day), "Unknown", full_seurat$day)


full_seurat <- RunHarmony(object = full_seurat, 
                          group.by.vars = "treatment_clean", 
                          reduction.use = "pca",
                          plot_convergence = TRUE,
                          max.iter.harmony = 20)


full_seurat <- RunUMAP(full_seurat, reduction = "harmony", dims = 1:30)

# Updated Seurat Trajectory w/ Batch Correction

cds <- as.cell_data_set(full_seurat, reduction = "umap") # cds = cell data set

cds <- cluster_cells(cds, reduction_method = 'UMAP', k = 20)
cds <- learn_graph(cds, 
                   use_partition = FALSE,
                   close_loop = FALSE,
                   learn_graph_control = list(
                     minimal_branch_len = 10,
                     geodesic_distance_ratio = 1/3
                   ))


p1_harmony <- DimPlot(full_seurat, reduction = "umap", group.by = "cell_type", pt.size = 1) + 
  ggtitle("UMAP with Harmony - Cell Type") +
  theme(legend.position = "right")

p2_harmony <- DimPlot(full_seurat, reduction = "umap", group.by = "treatment", pt.size = 1) + 
  ggtitle("UMAP with Harmony - Treatment") +
  theme(legend.position = "right")

p3_harmony <- DimPlot(full_seurat, reduction = "umap", group.by = "Phase", pt.size = 1) + 
  ggtitle("UMAP with Harmony - Cell Cycle") +
  theme(legend.position = "right")

p4_harmony <- DimPlot(full_seurat, reduction = "umap", group.by = "day", pt.size = 1) +
  ggtitle("UMAP with Harmony - Day") +
  theme(legend.position = "right")

combined_harmony_plot <- p1_harmony + p2_harmony + p3_harmony + p4_harmony
print(combined_harmony_plot)

print(p4_harmony)


#---------------------------------- Monocle3 w/ Harmony Adjustment

# Convert Harmony-corrected Seurat object to monocle3 cell_data_set

cds_harmony <- as.cell_data_set(full_seurat, reduction = 'umap')

# Cluster cells using the harmony-corrected UMAP
cds_harmony <- cluster_cells(cds_harmony, reduction_method = 'UMAP', k = 20)

# Learn the trajectory graph
cds_harmony <- learn_graph(cds_harmony, 
                           use_partition = FALSE,
                           close_loop = FALSE,
                           learn_graph_control = list(
                             minimal_branch_len = 15,        # Increased from 10
                             geodesic_distance_ratio = 1/4,  # Decreased from 1/3
                             prune_graph = TRUE              # Additional pruning
                           ))

# Monocle3 root with BFUE cells as they're developmentally earlier
bfue_g1_cells <- which(full_seurat$cell_type == "BFUE" & full_seurat$Phase == "G1")
if(length(bfue_g1_cells) > 5) {
  root_cells <- bfue_g1_cells
  print("Using BFUE G1 cells as root")
} else {
  # If not enough BFUE G1 cells, use all BFUE cells
  root_cells <- which(full_seurat$cell_type == "BFUE")
  print("Using all BFUE cells as root")
}

# Order cells
cds_harmony <- order_cells(cds_harmony, root_cells = colnames(cds_harmony)[root_cells])

# Basic visualization of the trajectory 
traj_plot1_harmony <- plot_cells(cds_harmony, 
                                 color_cells_by = "cell_type", 
                                 label_branch_points = TRUE,
                                 graph_label_size = 3,     
                                 cell_size = 0.8) + theme(legend.position = "right")

traj_plot2_harmony <- plot_cells(cds_harmony, 
                                 color_cells_by = "Phase",
                                 label_branch_points = TRUE,
                                 graph_label_size = 3,
                                 cell_size = 0.8) + theme(legend.position = "right")

traj_plot3_harmony <- plot_cells(cds_harmony, 
                                 color_cells_by = "pseudotime", 
                                 label_branch_points = TRUE,
                                 graph_label_size = 3,
                                 cell_size = 0.8) + theme(legend.position = "right")

print(traj_plot1_harmony)
print(traj_plot2_harmony)
print(traj_plot3_harmony)


# Get pseudotime values from the harmony-corrected trajectory
pt_values_harmony <- pseudotime(cds_harmony)
pseudotime_scaled_harmony <- (pt_values_harmony - min(pt_values_harmony)) / 
  (max(pt_values_harmony) - min(pt_values_harmony)) * 125

# Create plotting data 
plot_data_harmony <- data.frame(
  Pseudotime = pseudotime_scaled_harmony,
  Cell_Type = cds_harmony$cell_type,
  Treatment = cds_harmony$treatment,
  Day = cds_harmony$day,
  Phase = colData(cds_harmony)$Phase,
  Identity = cds_harmony$identity
)

# Add a simplified grouping for BFUE cells with days
plot_data_harmony$BFUE_Group <- case_when(
  plot_data_harmony$Cell_Type == "BFUE" & !is.na(plot_data_harmony$Day) & plot_data_harmony$Treatment == "DEX" ~ 
    paste("+ Dex Day", plot_data_harmony$Day),
  plot_data_harmony$Cell_Type == "BFUE" & !is.na(plot_data_harmony$Day) & plot_data_harmony$Treatment != "DEX" ~ 
    paste("no Dex Day", plot_data_harmony$Day),
  plot_data_harmony$Cell_Type == "BFUE" & is.na(plot_data_harmony$Day) ~ "BFUE (no day)",
  plot_data_harmony$Cell_Type == "CFUE" ~ "CFUE",
  TRUE ~ as.character(plot_data_harmony$Identity)
)

# For visualizations similar to your original code
# Setting factor levels for BFUE with day designations
bfue_with_day_harmony <- plot_data_harmony[!is.na(plot_data_harmony$Day) & plot_data_harmony$Cell_Type == "BFUE", ]
if(nrow(bfue_with_day_harmony) > 0) {
  bfue_groups_harmony <- unique(bfue_with_day_harmony$BFUE_Group)
  bfue_order_harmony <- c(
    "+ Dex Day Day3", "+ Dex Day Day2", "+ Dex Day Day1",
    "no Dex Day Day3", "no Dex Day Day2", "no Dex Day Day1"
  )
  bfue_order_harmony <- bfue_order_harmony[bfue_order_harmony %in% bfue_groups_harmony]
  
  # Create a factor for BFUE groups only
  plot_data_harmony$BFUE_Group_Factor <- factor(
    ifelse(plot_data_harmony$BFUE_Group %in% bfue_order_harmony, plot_data_harmony$BFUE_Group, NA),
    levels = bfue_order_harmony
  )
}

#-------------------------------"Harmony-Corrected Erythroid Developmental Trajectory"

# Plot 1: Standard pseudotime for BFUE cells with day designation
p1_harmony <- ggplot(subset(plot_data_harmony, !is.na(BFUE_Group_Factor)), 
                     aes(x = Pseudotime, y = BFUE_Group_Factor, color = Phase)) +
  geom_jitter(size = 1.5, width = 0, height = 0.2) +
  scale_color_manual(values = c(
    "G1" = "purple",
    "G2M" = "blue",
    "S" = "orange"
  )) +
  scale_x_continuous(breaks = seq(0, 125, 25)) +
  theme_classic() +
  labs(x = "Developmental progression (Harmony)", title = "BFUE Cells by Day & Treatment (Harmony)", y = "") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray90"),
    panel.grid.major.y = element_line(color = "gray90")
  )

# Create a comprehensive pseudotime plot for ALL cells with harmony-corrected data

# Group by cell type with customized point shapes and colors
all_cells_plot_harmony <- ggplot(plot_data_harmony, 
                                 aes(x = Pseudotime, y = Cell_Type, color = Cell_Type, shape = Phase)) +
  geom_jitter(size = 2, width = 0, height = 0.3, alpha = 0.8) +
  scale_color_brewer(palette = "Set1") +
  scale_shape_manual(values = c("G1" = 16, "S" = 17, "G2M" = 15)) +
  scale_x_continuous(breaks = seq(0, 125, 25)) +
  theme_classic() +
  labs(
    x = "Developmental Progression (Harmony Pseudotime)",
    y = "",
    title = "Harmony-Corrected Erythroid Developmental Trajectory",
    subtitle = "All cells arranged by pseudotime after batch correction"
  ) +
  theme(
    legend.position = "right",
    panel.grid.major.x = element_line(color = "gray90"),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12, face = "italic")
  )

# Add density curves on top to show cell distribution along pseudotime
density_plot_harmony <- ggplot(plot_data_harmony, aes(x = Pseudotime, fill = Cell_Type)) +
  geom_density(alpha = 0.5) +
  scale_fill_brewer(palette = "Set1") +
  scale_x_continuous(breaks = seq(0, 125, 25)) +
  theme_classic() +
  labs(x = "Harmony Pseudotime", y = "Density") +
  theme(legend.position = "none")

# Combine the plots with patchwork
combined_pseudotime_harmony <- density_plot_harmony / all_cells_plot_harmony +
  plot_layout(heights = c(1, 4))


#-----------------------------"Cell Type Distribution by Treatment after Harmony"

# Create comparison plots to allow direct comparison of original and harmony pseudotime
treatment_pseudotime_harmony <- ggplot(plot_data_harmony, 
                                       aes(x = Pseudotime, y = Cell_Type, color = Treatment)) +
  geom_jitter(size = 2, width = 0, height = 0.3, alpha = 0.8) +
  scale_x_continuous(breaks = seq(0, 125, 25)) +
  theme_classic() +
  labs(
    x = "Developmental Progression (Harmony Pseudotime)",
    y = "",
    title = "Cell Type Distribution by Treatment after Harmony",
    color = "Treatment"
  ) +
  theme(
    legend.position = "right",
    panel.grid.major.x = element_line(color = "gray90"),
    panel.grid.major.y = element_line(color = "gray90")
  )

# Display all plots
print(p1_harmony)
print(combined_pseudotime_harmony)
print(treatment_pseudotime_harmony)


#-------------------------------"Cell Type Distribution by Treatment and Day"


# Variation 1: Facet by day (shows cell types separated by day)
day_treatment_pseudotime_harmony <- ggplot(plot_data_harmony, 
                                           aes(x = Pseudotime, y = Cell_Type, color = Treatment)) +
  geom_jitter(size = 2, width = 0, height = 0.3, alpha = 0.8) +
  facet_wrap(~ Day, scales = "free_y", ncol = 1) +  # Facets by day
  scale_x_continuous(breaks = seq(0, 125, 25)) +
  theme_classic() +
  labs(
    x = "Developmental Progression (Harmony Pseudotime)",
    y = "",
    title = "Cell Type Distribution by Treatment and Day",
    color = "Treatment"
  ) +
  theme(
    legend.position = "right",
    panel.grid.major.x = element_line(color = "gray90"),
    panel.grid.major.y = element_line(color = "gray90"),
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(face = "bold")
  )

# Variation 2: Create a combined cell type and day identifier
plot_data_harmony$CellDay <- ifelse(is.na(plot_data_harmony$Day),
                                    paste(plot_data_harmony$Cell_Type, "(No Day)"),
                                    paste(plot_data_harmony$Cell_Type, "Day", plot_data_harmony$Day))

# Create ordered factor to control display order
cell_types <- unique(plot_data_harmony$Cell_Type)
days <- c("Day1", "Day2", "Day3")
combined_levels <- c()

# First add all cell types with specific days in order
for (cell in cell_types) {
  for (day in days) {
    combined_levels <- c(combined_levels, paste(cell, "Day", day))
  }
  # Add entry for cells without day information
  combined_levels <- c(combined_levels, paste(cell, "(No Day)"))
}

# Remove any levels that don't exist in your data
combined_levels <- combined_levels[combined_levels %in% unique(plot_data_harmony$CellDay)]

# Set factor levels
plot_data_harmony$CellDay <- factor(plot_data_harmony$CellDay, levels = combined_levels)

# Create plot with combined cell type and day on y-axis
cellday_treatment_pseudotime <- ggplot(plot_data_harmony, 
                                       aes(x = Pseudotime, y = CellDay, color = Treatment)) +
  geom_jitter(size = 2, width = 0, height = 0.3, alpha = 0.8) +
  scale_x_continuous(breaks = seq(0, 125, 25)) +
  theme_classic() +
  labs(
    x = "Developmental Progression (Harmony Pseudotime)",
    y = "",
    title = "Cell Type by Day and Treatment",
    color = "Treatment"
  ) +
  theme(
    legend.position = "right",
    panel.grid.major.x = element_line(color = "gray90"),
    panel.grid.major.y = element_line(color = "gray90"),
    axis.text.y = element_text(size = 9)
  )

# Print both options to see which you prefer
print(day_treatment_pseudotime_harmony)
print(cellday_treatment_pseudotime)



#----------------------------------PSEUDOTIME BINS SORTED BY CELL PHASE

# Create pseudotime bins
plot_data_harmony$PseudotimeBin <- cut(plot_data_harmony$Pseudotime, 
                                       breaks = seq(0, 125, 25),
                                       labels = c("0-25", "25-50", "50-75", "75-100", "100-125"),
                                       include.lowest = TRUE)

# Filter out cells with no assigned day
plot_data_harmony_days <- plot_data_harmony %>%
  filter(!is.na(Day))

# Summarize cell counts by cell type, bin, and phase
phase_summary <- plot_data_harmony_days %>%
  group_by(Cell_Type, PseudotimeBin, Phase) %>%
  summarise(Count = n(), .groups = 'drop')

# Get day counts for each cell type and bin
day_counts <- plot_data_harmony_days %>%
  group_by(Cell_Type, PseudotimeBin, Day) %>%
  summarise(Count = n(), .groups = 'drop')

# Get total cells per cell type and bin
celltype_counts <- phase_summary %>%
  group_by(Cell_Type, PseudotimeBin) %>%
  summarise(TotalCells = sum(Count), .groups = 'drop')

# Get total cells per bin
bin_totals <- celltype_counts %>%
  group_by(PseudotimeBin) %>%
  summarise(BinTotal = sum(TotalCells), .groups = 'drop')

# Set output directory
output_dir <- "C:/Users/aidan/OneDrive/Desktop/LI LAB Pseudotime Stuff"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Create pie charts for each bin
for(bin in levels(plot_data_harmony$PseudotimeBin)) {
  # Filter data for this bin
  bin_data <- phase_summary %>% filter(PseudotimeBin == bin)
  bin_total <- bin_totals %>% filter(PseudotimeBin == bin) %>% pull(BinTotal)
  
  # Skip if no data for this bin
  if(nrow(bin_data) == 0 || is.null(bin_total) || bin_total == 0) {
    next
  }
  
  # Create pie charts and day count text for each cell type
  pie_charts <- list()
  day_texts <- list()
  
  for(cell_type in unique(bin_data$Cell_Type)) {
    # Get data for this cell type
    ct_data <- bin_data %>% filter(Cell_Type == cell_type)
    ct_count <- celltype_counts %>% 
      filter(Cell_Type == cell_type, PseudotimeBin == bin) %>% 
      pull(TotalCells)
    
    # Get day counts for this cell type
    day_data <- day_counts %>% 
      filter(Cell_Type == cell_type, PseudotimeBin == bin) %>%
      arrange(Day)
    
    # Create day count text
    day_text <- ""
    if(nrow(day_data) > 0) {
      day_counts_str <- paste(day_data$Day, ": n=", day_data$Count, sep="", collapse="\n")
      day_text <- grid::textGrob(day_counts_str, gp = grid::gpar(fontsize = 10))
      day_texts[[cell_type]] <- day_text
    }
    
    if(nrow(ct_data) > 0) {
      # Create pie chart with count labels
      pie <- ggplot(ct_data, aes(x = "", y = Count, fill = Phase)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar("y", start = 0) +
        scale_fill_manual(values = c("G1" = "#FF9E9E", "G2M" = "#00C094", "S" = "#81B1FF")) +
        geom_text(aes(label = Count), position = position_stack(vjust = 0.5), color = "black", size = 4) +
        ggtitle(paste(cell_type, "- n =", ct_count)) +
        theme_void() +
        theme(plot.title = element_text(hjust = 0.5, size = 14))
      
      pie_charts[[cell_type]] <- pie
    }
  }
  
  # Create a separate image for this bin
  if(length(pie_charts) > 0) {
    # Add a title for this bin
    title <- paste0("Pseudotime Bin: ", bin, " (Total cells with days: ", bin_total, ")")
    
    # For each cell type, create a combined grob with pie chart and day count text
    combined_grobs <- list()
    for(cell_type in names(pie_charts)) {
      # Create a combined grob for this cell type
      combined_grob <- gridExtra::arrangeGrob(
        pie_charts[[cell_type]],
        day_texts[[cell_type]],
        ncol = 1,
        heights = c(4, 1)
      )
      combined_grobs[[cell_type]] <- combined_grob
    }
    
    # Create the plot with all cell types
    bin_plot <- gridExtra::grid.arrange(
      grobs = combined_grobs,
      ncol = length(combined_grobs),
      top = title
    )
    
    # Save the plot for this bin
    ggsave(
      filename = file.path(output_dir, paste0("CellCycle_Phase_Bin_", bin, "_days_only.png")),
      plot = bin_plot,
      width = 2.5 * length(pie_charts),
      height = 5,
      dpi = 300
    )
    
    # Print progress
    cat("Created plot for bin:", bin, "\n")
  }
}


# ------------------Bin sorted by DEX vs No DEX (ONLY DAYS)

# Create pseudotime bins
plot_data_harmony$PseudotimeBin <- cut(plot_data_harmony$Pseudotime, 
                                       breaks = seq(0, 125, 25),
                                       labels = c("0-25", "25-50", "50-75", "75-100", "100-125"),
                                       include.lowest = TRUE)

# Filter out cells with no assigned day
plot_data_harmony_days <- plot_data_harmony %>%
  filter(!is.na(Day))

# Summarize cell counts by cell type, bin, and treatment
treatment_summary <- plot_data_harmony_days %>%
  group_by(Cell_Type, PseudotimeBin, Treatment) %>%
  summarise(Count = n(), .groups = 'drop')

# Get total cells per cell type and bin
celltype_counts <- treatment_summary %>%
  group_by(Cell_Type, PseudotimeBin) %>%
  summarise(TotalCells = sum(Count), .groups = 'drop')

# Get total cells per bin
bin_totals <- celltype_counts %>%
  group_by(PseudotimeBin) %>%
  summarise(BinTotal = sum(TotalCells), .groups = 'drop')

# Set output directory
output_dir <- "C:/Users/aidan/OneDrive/Desktop/LI LAB Pseudotime Stuff"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Create a separate image for each bin
for(bin in levels(plot_data_harmony$PseudotimeBin)) {
  # Filter data for this bin
  bin_data <- treatment_summary %>% filter(PseudotimeBin == bin)
  bin_total <- bin_totals %>% filter(PseudotimeBin == bin) %>% pull(BinTotal)
  
  # Skip if no data for this bin
  if(nrow(bin_data) == 0 || is.null(bin_total) || bin_total == 0) {
    next
  }
  
  # Create pie charts for each cell type
  pie_charts <- list()
  
  for(cell_type in unique(bin_data$Cell_Type)) {
    # Get data for this cell type
    ct_data <- bin_data %>% filter(Cell_Type == cell_type)
    ct_count <- celltype_counts %>% 
      filter(Cell_Type == cell_type, PseudotimeBin == bin) %>% 
      pull(TotalCells)
    
    if(nrow(ct_data) > 0) {
      # Create pie chart with count labels
      pie <- ggplot(ct_data, aes(x = "", y = Count, fill = Treatment)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar("y", start = 0) +
        scale_fill_manual(values = c("Control" = "#FF9E9E", "DEX" = "#00BFC4")) +
        geom_text(aes(label = Count), position = position_stack(vjust = 0.5), color = "black", size = 4) +
        ggtitle(paste(cell_type, "- n =", ct_count)) +
        theme_void() +
        theme(plot.title = element_text(hjust = 0.5, size = 14))
      
      pie_charts[[cell_type]] <- pie
    }
  }
  
  # Create a separate image for this bin
  if(length(pie_charts) > 0) {
    # Add a title for this bin
    title <- paste0("Pseudotime Bin: ", bin, " (Total cells with days: ", bin_total, ")")
    
    # Create the plot
    bin_plot <- gridExtra::grid.arrange(
      grobs = pie_charts,
      ncol = length(pie_charts),
      top = title
    )
    
    # Save the plot for this bin
    ggsave(
      filename = file.path(output_dir, paste0("Treatment_Distribution_Bin_", bin, ".png")),
      plot = bin_plot,
      width = 2.5 * length(pie_charts),
      height = 4,
      dpi = 300
    )
    
    # Print progress
    cat("Created plot for bin:", bin, "\n")
  }
}


# Create pseudotime bins
plot_data_harmony$PseudotimeBin <- cut(plot_data_harmony$Pseudotime, 
                                       breaks = seq(0, 125, 25),
                                       labels = c("0-25", "25-50", "50-75", "75-100", "100-125"),
                                       include.lowest = TRUE)

# Filter out cells with no assigned day
plot_data_harmony_days <- plot_data_harmony %>%
  filter(!is.na(Day))

# Summarize cell counts by cell type, bin, and treatment
treatment_summary <- plot_data_harmony_days %>%
  group_by(Cell_Type, PseudotimeBin, Treatment) %>%
  summarise(Count = n(), .groups = 'drop')

# Get day counts for each cell type and bin
day_counts <- plot_data_harmony_days %>%
  group_by(Cell_Type, PseudotimeBin, Day) %>%
  summarise(Count = n(), .groups = 'drop')

# Get total cells per cell type and bin
celltype_counts <- treatment_summary %>%
  group_by(Cell_Type, PseudotimeBin) %>%
  summarise(TotalCells = sum(Count), .groups = 'drop')

# Get total cells per bin
bin_totals <- celltype_counts %>%
  group_by(PseudotimeBin) %>%
  summarise(BinTotal = sum(TotalCells), .groups = 'drop')

# Set output directory
output_dir <- "C:/Users/aidan/OneDrive/Desktop/LI LAB Pseudotime Stuff"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Create a separate image for each bin
for(bin in levels(plot_data_harmony$PseudotimeBin)) {
  # Filter data for this bin
  bin_data <- treatment_summary %>% filter(PseudotimeBin == bin)
  bin_total <- bin_totals %>% filter(PseudotimeBin == bin) %>% pull(BinTotal)
  
  # Skip if no data for this bin
  if(nrow(bin_data) == 0 || is.null(bin_total) || bin_total == 0) {
    next
  }
  
  # Create pie charts and day count text for each cell type
  pie_charts <- list()
  day_texts <- list()
  
  for(cell_type in unique(bin_data$Cell_Type)) {
    # Get data for this cell type
    ct_data <- bin_data %>% filter(Cell_Type == cell_type)
    ct_count <- celltype_counts %>% 
      filter(Cell_Type == cell_type, PseudotimeBin == bin) %>% 
      pull(TotalCells)
    
    # Get day counts for this cell type
    day_data <- day_counts %>% 
      filter(Cell_Type == cell_type, PseudotimeBin == bin) %>%
      arrange(Day)
    
    # Create day count text
    day_text <- ""
    if(nrow(day_data) > 0) {
      day_counts_str <- paste(day_data$Day, ": n=", day_data$Count, sep="", collapse="\n")
      day_text <- grid::textGrob(day_counts_str, gp = grid::gpar(fontsize = 10))
      day_texts[[cell_type]] <- day_text
    }
    
    if(nrow(ct_data) > 0) {
      # Create pie chart with count labels
      pie <- ggplot(ct_data, aes(x = "", y = Count, fill = Treatment)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar("y", start = 0) +
        scale_fill_manual(values = c("Control" = "#FF9E9E", "DEX" = "#00BFC4")) +
        geom_text(aes(label = Count), position = position_stack(vjust = 0.5), color = "black", size = 4) +
        ggtitle(paste(cell_type, "- n =", ct_count)) +
        theme_void() +
        theme(plot.title = element_text(hjust = 0.5, size = 14))
      
      pie_charts[[cell_type]] <- pie
    }
  }
  
  # Create a separate image for this bin
  if(length(pie_charts) > 0) {
    # Add a title for this bin
    title <- paste0("Pseudotime Bin: ", bin, " (Total cells with days: ", bin_total, ")")
    
    # For each cell type, create a combined grob with pie chart and day count text
    combined_grobs <- list()
    for(cell_type in names(pie_charts)) {
      # Create a combined grob for this cell type
      combined_grob <- gridExtra::arrangeGrob(
        pie_charts[[cell_type]],
        day_texts[[cell_type]],
        ncol = 1,
        heights = c(4, 1)
      )
      combined_grobs[[cell_type]] <- combined_grob
    }
    
    # Create the plot with all cell types
    bin_plot <- gridExtra::grid.arrange(
      grobs = combined_grobs,
      ncol = length(combined_grobs),
      top = title
    )
    
    # Save the plot for this bin
    ggsave(
      filename = file.path(output_dir, paste0("Treatment_Distribution_Bin_", bin, "_with_days.png")),
      plot = bin_plot,
      width = 2.5 * length(pie_charts),
      height = 5,
      dpi = 300
    )
    
    # Print progress
    cat("Created plot for bin:", bin, "\n")
  }
}
