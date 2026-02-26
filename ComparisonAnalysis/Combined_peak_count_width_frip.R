#########################################################################################################################################
#####  Combination Frips and peak summary results for HumanPBMC data with input 
#########################################################################################################################################

library(dplyr)
library(stringr)
library(dplyr)
library(readr)
library(stringr)

# Define file paths
peak_file <- "~/result_frip/HumanPBMC_peak_summary_with_input.csv"
frip_file <- "~/result_frip/HumanPBMC_frip_with_input.csv"
output_file <- "~/result_frip/HumanPBMC_combined_peak_frip_summary_with_input.csv"

# Read data
peak_df <- read_csv(peak_file)
frip_df <- read_csv(frip_file)

peak_df_updated <- peak_df %>%
  mutate(
    HM = case_when(
      `Histone/TF` %in% c("H3K27ac-b", "H3K27ac-s") ~ "H3K27ac",
      `Histone/TF` == "H3K27me3" ~ "H3K27me3",
      `Histone/TF` == "H3K4me1" ~ "H3K4me1",
      `Histone/TF` == "H3K4me2" ~ "H3K4me2",
      `Histone/TF` == "H3K4me3" ~ "H3K4me3",
      `Histone/TF` == "H3K9me3" ~ "H3K9me3",
      TRUE ~ `Histone/TF`
    ),
    .before = 1  # add new column before 1st one
  ) %>%
  rename(
    Histone = `Histone/TF`  # Reanming it
  )

# Addnew column and going to processing the FRiPs values 
add_new_columns <- function(frip_df) {
  frip_df %>%
    mutate(
      # New column: Histone mark + sample type (H3K27ac-b, H3K27ac-s, H3K9me3, )
      Histone_SampleType = case_when(
        str_detect(Sample, "H3K27ac-b") ~ "H3K27ac-b",
        str_detect(Sample, "H3K27ac-s") ~ "H3K27ac-s",
        str_detect(Sample, "H3K27me3") ~ "H3K27me3",
        str_detect(Sample, "H3K4me1") ~ "H3K4me1",
        str_detect(Sample, "H3K4me2") ~ "H3K4me2",
        str_detect(Sample, "H3K4me3") ~ "H3K4me3",
        str_detect(Sample, "H3K9me3") ~ "H3K9me3",
        TRUE ~ NA_character_
      ),
      
      # Secodly column:  cell type (B, CD4T, CD8T, DC, Mono, NK, other, otherT)
      CellType = case_when(
        str_detect(Sample, "_B$") ~ "B",
        str_detect(Sample, "_CD4T$") ~ "CD4T",
        str_detect(Sample, "_CD8T$") ~ "CD8T",
        str_detect(Sample, "_DC$") ~ "DC",
        str_detect(Sample, "_Mono$") ~ "Mono",
        str_detect(Sample, "_NK$") ~ "NK",
        str_detect(Sample, "_other$") ~ "other",
        str_detect(Sample, "_otherT$") ~ "otherT",
        TRUE ~ NA_character_
      )
    ) %>%
    # column reorganize 
    select(HM, Histone, Method, Sample, Histone_SampleType, CellType, everything())
}

# Function apply
frip_df_updated <- add_new_columns(frip_df)

# ensure the order of whole data
custom_order <- c("H3K27ac-b", "H3K27ac-s", "H3K27me3", "H3K4me1", "H3K4me2", "H3K4me3", "H3K9me3")
frip_df_updated <- frip_df_updated %>%
  mutate(
    Histone_SampleType = factor(Histone_SampleType, levels = custom_order)
  ) %>%
  arrange(Histone_SampleType)
# Operations make sure on single line
frip_df_updated <- frip_df_updated %>%
  select(HM, Histone = Histone_SampleType, Method, Sample, Total_Reads, Reads_in_Peaks, FRiP)


library(dplyr)
library(stringr)

allowed_celltypes <- c("otherT", "B", "CD4T", "CD8T", "DC", "Mono", "NK", "other")

peak_df <- peak_df_updated %>%
  mutate(celltype = str_extract(Sample, paste(allowed_celltypes, collapse = "|"))) %>%
  filter(celltype %in% allowed_celltypes) %>%
  distinct(HM, Histone, Method, celltype, .keep_all = TRUE)

frip_df <- frip_df_updated %>%
  mutate(celltype = str_extract(Sample, paste(allowed_celltypes, collapse = "|"))) %>%
  filter(celltype %in% allowed_celltypes) %>%
  distinct(HM, Histone, Method, celltype, .keep_all = TRUE)
merged_df <- inner_join(
  peak_df %>% select(HM, Histone, Method, celltype, Number_of_Peaks, Mean_Peak_Width),
  frip_df %>% select(HM, Histone, Method, celltype, Total_Reads, Reads_in_Peaks, FRiP),
  by = c("HM", "Histone", "Method", "celltype")
)

final_df <- merged_df %>%
  rename(Celltype = celltype) %>%
  select(HM, Histone, Method, Celltype, Number_of_Peaks, Mean_Peak_Width, Total_Reads, Reads_in_Peaks, FRiP)

final_df[is.na(final_df)] <- 0
# Save result
write_csv(final_df, output_file)


###===========================================================================================================###
###  Combination Frips and peak summary results for MouseBrain data with input 
###===========================================================================================================###
library(dplyr)
library(stringr)
library(dplyr)
library(readr)
library(stringr)

# Define file paths
peak_file <- "~/result_frip/MouseBrain_peak_summary_with_input.csv"
frip_file <- "~/result_frip/MouseBrain_frip_with_input.csv"
output_file <- "~/result_frip/MouseBrain_combined_peak_frip_summary_with_input.csv"

# Read data
peak_df <- read_csv(peak_file)
frip_df <- read_csv(frip_file)

peak_df_updated <- peak_df %>%
  mutate(
    HM = case_when(
      `Histone/TF` %in% c("H3K27ac-b", "H3K27ac-s") ~ "H3K27ac",
      `Histone/TF` == "H3K27me3" ~ "H3K27me3",
      `Histone/TF` == "H3K36me3" ~ "H3K36me3",
      `Histone/TF` == "H3K4me3" ~ "H3K4me3",
      `Histone/TF` == "Olig2" ~ "Olig2",
      `Histone/TF` == "Rad21" ~ "Rad21",
      TRUE ~ `Histone/TF`
    ),
    .before = 1  # add new column before 1st one
  ) %>%
  rename(
    Histone = `Histone/TF`  # Reanming it
  )

# Addnew column and going to processing the FRiPs values 
add_new_columns <- function(frip_df) {
  frip_df %>%
    mutate(
      # New column: Histone mark + sample type (H3K27ac-b, H3K27ac-s, H3K9me3, )
      Histone_SampleType = case_when(
        str_detect(Sample, "H3K27ac-b") ~ "H3K27ac-b",
        str_detect(Sample, "H3K27ac-s") ~ "H3K27ac-s",
        str_detect(Sample, "H3K27me3") ~ "H3K27me3",
        str_detect(Sample, "H3K36me3") ~ "H3K36me3",
        str_detect(Sample, "H3K4me3") ~ "H3K4me3",
        str_detect(Sample, "Olig2") ~ "Olig2",
        str_detect(Sample, "Rad21") ~ "Rad21",
        TRUE ~ NA_character_
      ),
      
      # Secodly column:  cell type (B, CD4T, CD8T, DC, Mono, NK, other, otherT)
      CellType = case_when(
        str_detect(Sample, "_Astrocytes$") ~ "Astrocytes",
        str_detect(Sample, "_Microglia$") ~ "Microglia",
        str_detect(Sample, "_mOL$") ~ "mOL",
        str_detect(Sample, "_Neurons1$") ~ "Neurons1",
        str_detect(Sample, "_Neurons2$") ~ "Neurons2",
        str_detect(Sample, "_Neurons3$") ~ "Neurons3",
        str_detect(Sample, "_OEC$") ~ "OEC",
        str_detect(Sample, "_OPC$") ~ "OPC",
        str_detect(Sample, "_VLMC$") ~ "VLMC",
        str_detect(Sample, "_Unknown$") ~ "Unknown",
        TRUE ~ NA_character_
      )
    ) %>%
    # column reorganize 
    select(HM, Histone, Method, Sample, Histone_SampleType, CellType, everything())
}

# Function apply
frip_df_updated <- add_new_columns(frip_df)

# ensure the order of whole data
custom_order <- c("H3K27ac-b", "H3K27ac-s", "H3K27me3", "H3K36me3", "H3K4me3", "Olig2", "Rad21")
frip_df_updated <- frip_df_updated %>%
  mutate(
    Histone_SampleType = factor(Histone_SampleType, levels = custom_order)
  ) %>%
  arrange(Histone_SampleType)
# Operations make sure on single line
frip_df_updated <- frip_df_updated %>%
  select(HM, Histone = Histone_SampleType, Method, Sample, Total_Reads, Reads_in_Peaks, FRiP)


library(dplyr)
library(stringr)

allowed_celltypes <- c("Astrocytes", "Microglia", "mOL", "Neurons1", "Neurons2", "Neurons3", "OEC", "OPC", "VLMC", "Unknown")

peak_df <- peak_df_updated %>%
  mutate(celltype = str_extract(Sample, paste(allowed_celltypes, collapse = "|"))) %>%
  filter(celltype %in% allowed_celltypes) %>%
  distinct(HM, Histone, Method, celltype, .keep_all = TRUE)

frip_df <- frip_df_updated %>%
  mutate(celltype = str_extract(Sample, paste(allowed_celltypes, collapse = "|"))) %>%
  filter(celltype %in% allowed_celltypes) %>%
  distinct(HM, Histone, Method, celltype, .keep_all = TRUE)
merged_df <- inner_join(
  peak_df %>% select(HM, Histone, Method, celltype, Number_of_Peaks, Mean_Peak_Width),
  frip_df %>% select(HM, Histone, Method, celltype, Total_Reads, Reads_in_Peaks, FRiP),
  by = c("HM", "Histone", "Method", "celltype")
)

final_df <- merged_df %>%
  rename(Celltype = celltype) %>%
  select(HM, Histone, Method, Celltype, Number_of_Peaks, Mean_Peak_Width, Total_Reads, Reads_in_Peaks, FRiP)

final_df[is.na(final_df)] <- 0
# Save result
write_csv(final_df, output_file)
