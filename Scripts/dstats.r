install.packages(c("nortest", "car", "readr","Rmisc","psych","fitdistrplus"))
library(nortest)
library(car)
library(readr)
library(Rmisc)
library(psych)
library(fitdistrplus)

#Load Dataframe ")
data <- read.csv("/Users/brandonhendrickson/Documents/Github_Projects/M_guttatus_ReproductiveAllocation/Data/csv/raw_data.csv")

#Make Trait List
trait_list <- list("yield","flower_number_1","D2Flowering","pollen_avg","ff_height","total_node_number","prop_early_growth","total_weight","leaf_area","flower_area")

# Describe mean, standard deviation, variance, range, and kurtosis
for (trait in trait_list) {
    summary <- describe(data[[trait]])
    print(paste("Summary for", trait, ":"))
    print(paste("Mean:", summary$mean))
    print(paste("Standard deviation:", summary$sd))
    print(paste("Kurtosis:", summary$kurtosis))
}

# Determine normality of variables using Anderson-Darling Test
for (trait in trait_list) {
    ad_test_result <- ad.test(data[[trait]])
    print(paste("Anderson-Darling test for", trait, ":"))
    print(paste("Trait Name:", trait))
    print(ad_test_result)
}

# Define the distributions and their starting parameters
distributions <- list(
    "norm" = NULL,
    "gamma" = list(shape = 1, rate = 1),
    "exp" = list(rate = 1),
    "logis" = list(location = 0, scale = 1),
    "geom" = list(prob = 0.5)
)

# Initialize an empty dataframe to store the results
summary_df <- data.frame(Trait = character(), Best_Model = character(), stringsAsFactors = FALSE)

for (trait in trait_list) {
    print(paste("Trait:", trait))
    
    # Remove NA values from the current trait
    trait_data <- na.omit(data[[trait]])
    
    # Check if the data is continuous or discrete
    is_discrete <- all(trait_data == round(trait_data))
    
    # Check if the data is greater than 0
    is_greater_than_zero <- all(trait_data > 0)
    
    # Fit the distributions to the data and store the results in a list
    fit_list <- lapply(names(distributions), function(dist) {
        if ((dist == "geom" && !is_discrete) || (dist %in% c("gamma", "exp") && !is_greater_than_zero)) {
            return(NULL)
        }
        fitdist(as.numeric(trait_data), dist, start = distributions[[dist]])
    })

    # Remove NULL elements from the list
    fit_list <- fit_list[!sapply(fit_list, is.null)]
    
    # Compare the fits
    compare_dist <- gofstat(fit_list)

    # Find the best model based on AIC
    best_model <- which.min(compare_dist$aic)
    
    # Add the result to the summary dataframe
    summary_df <- rbind(summary_df, data.frame(Trait = trait, Best_Model = names(best_model), stringsAsFactors = FALSE))
}

# Print the summary dataframe
print(summary_df)

# Equality of Varinces betweenn age classes
# Define a function to perform Levene's test and print the result
perform_levene_test <- function(trait, group, data) {
    formula <- as.formula(paste(trait, "~", group))
    result <- leveneTest(formula, data = data)
    print(result)
}

# Initialize an empty dataframe to store the results
all_tests_df <- data.frame(Trait = character(), Group = character(), P_Value = numeric(), stringsAsFactors = FALSE)

for (trait in trait_list) {
    # Perform the Levene's test for each group
    test_result_age <- perform_levene_test(trait, "age", data)
    test_result_F_LF <- perform_levene_test(trait, "F_LF", data)
    
    # Add the result to the dataframe
    all_tests_df <- rbind(all_tests_df, data.frame(Trait = trait, Group = "age", P_Value = test_result_age$"Pr(>F)"[1], stringsAsFactors = FALSE))
    all_tests_df <- rbind(all_tests_df, data.frame(Trait = trait, Group = "F_LF", P_Value = test_result_F_LF$"Pr(>F)"[1], stringsAsFactors = FALSE))
}

# Print the dataframe
print(all_tests_df)

test_model <- merge(all_tests_df, summary_df, by = "Trait")
# Remove the "2-mle-" string from the "Best_Model" column
test_model$Best_Model <- sub("\\d-mle-", "", test_model$Best_Model)


write.csv(test_model, "/Users/brandonhendrickson/Documents/Github_Projects/M_guttatus_ReproductiveAllocation/Data/csv/test_model.csv")