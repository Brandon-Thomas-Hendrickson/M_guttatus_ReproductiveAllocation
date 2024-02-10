# Load the necessary libraries
library(dplyr)

# Create the trait List
trait_list <- c("flower_number_1", "D2Flowering", "pollen_avg", "ff_height", "total_node_number", "prop_early_growth", "final_height", "veg_weight", "leaf_area")

# Load modelframe
data <- read.csv("/Users/brandonhendrickson/Documents/Github_Projects/M_guttatus_ReproductiveAllocation/Data/csv/raw_data.csv")

# Change Final Density to Factor
data$Final_Density <- as.factor(data$Final_Density)

# Modify the function to standardize a single trait
standardize <- function(data, trait) {
    data <- data %>%
        group_by(Final_Density) %>%
        mutate(!!paste0(trait, "_std") := ( .data[[trait]] - mean(.data[[trait]], na.rm = TRUE)) / sd(.data[[trait]], na.rm = TRUE))
    return(data)
}

# Apply the function to all elements in the trait list
for (trait in trait_list) {
    data <- standardize(data, trait)
}

# Save the data
write.csv(data, "/Users/brandonhendrickson/Documents/Github_Projects/M_guttatus_ReproductiveAllocation/Data/csv/standardized_data.csv", row.names = FALSE)

# Load modelframe
data <- read.csv("/Users/brandonhendrickson/Documents/Github_Projects/M_guttatus_ReproductiveAllocation/Data/csv/raw_data.csv")

# Create the trait List
trait_list <- c("flower_number_1", "D2Flowering", "pollen_avg", "ff_height", "total_node_number", "prop_early_growth", "final_height", "veg_weight","leaf_area")

# Change Final Density to Factor
data$Final_Density <- as.factor(data$Final_Density)

data <- data[data$F_LF == "F" & data$age != "1.E",]

# Modify the function to standardize a single trait
standardize <- function(data, trait) {
    data <- data %>%
        group_by(Final_Density) %>%
        mutate(!!paste0(trait, "_std") := ( .data[[trait]] - mean(.data[[trait]], na.rm = TRUE)) / sd(.data[[trait]], na.rm = TRUE))
    return(data)
}

# Apply the function to all elements in the trait list
for (trait in trait_list) {
    data <- standardize(data, trait)
}

# Save data_F
write.csv(data, "/Users/brandonhendrickson/Documents/Github_Projects/M_guttatus_ReproductiveAllocation/Data/csv/standardized_data_F.csv", row.names = FALSE)
