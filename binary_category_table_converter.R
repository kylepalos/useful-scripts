# Sometimes I need to investigate which genes are in which functional categories
# For example, a venn diagram or upset plot may suggest a certain category of genes are doing something interesting
# However, when there are more than a couple of categories, this can be tedious to investigate
# The following code generates a table that is easily parsed for categories of interest:


# Load required libraries
library(dplyr)
library(tidyr)


# Function to convert 2-column data frame to wide data format
convert_to_wide <- function(input_df) {
  # Group by gene_id and type, then summarize the counts
  summarized_df <- input_df %>%
    group_by(gene_id, type) %>%
    summarize(count = n()) %>%
    ungroup()
  
  # Pivot the data to wide format with 0s and 1s representing absence or presence in the category
  wide_df <- summarized_df %>%
    pivot_wider(names_from = type, values_from = count, values_fill = 0)
  
  return(wide_df)
}

# Don't need to write a function if you don't want to


# Generate a sample input data frame (for demonstration purposes)
sample_data <- data.frame(
  gene_id = c(rep("Gene_A", 2), rep("Gene_B", 3), rep("Gene_C", 2)),
  type = c("Category_1", "Category_2", "Category_1", "Category_3", "Category_2", "Category_3", "Category_1")
)

# Call the function to convert to wide data format
combined_wide <- convert_to_wide(sample_data)



# (Optional) Save the wide data format table to a CSV file for sharing on GitHub
# write.csv(combined_wide, file = "wide_data_format_table.csv", row.names = TRUE, quote = F)
