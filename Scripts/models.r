# Install the package
install.packages(c("geepack", "readr", "sjPlot", "MASS", "stats","emmeans"))

# Load the package
library(readr)
library(sjPlot)
library(MASS)
library(stats)
library(emmeans)

############################################################################################################
## General Response of Fitness, Phenology, and Reproductive Allocation to Growing Length and Competition ###
############################################################################################################

#Load modelframe
data <- read.csv("/Users/brandonhendrickson/Documents/Github_Projects/M_guttatus_ReproductiveAllocation/Data/csv/raw_data.csv")

#Identify numeric columns
numeric_cols <- sapply(data, is.numeric)

#Replace 0 with 0.000001
data[numeric_cols] <- lapply(data[numeric_cols], function(x) {
  ifelse(x == 0, 0.000001, x)
})
# Create the trait List
trait_list <- list("yield","flower_number_1","D2Flowering","pollen_avg")

# Initialize a list to store the models
glm_models <- list()

# Define the predictors
predictors <- c("age", "Final_Density","ID")

# Generate all combinations of predictors
combinations <- lapply(1:length(predictors), function(x) combn(predictors, x, simplify = FALSE))

# Flatten the list
combinations <- unlist(combinations, recursive = FALSE)

# Loop over each trait in the trait list
for (trait in trait_list) {
    # Remove NA values from the current trait
    data_clean <- data[!is.na(data[[trait]]), ]

    # Loop over each combination of predictors
    for (combination in combinations) {
        # Create a formula for the model with the current combination of predictors
        if (length(combination) == 3) {
            formula <- as.formula(paste(trait, "~ (", paste(combination, collapse = "+"), ")^3"))
        } else if (length(combination) == 2) {
            formula <- as.formula(paste(trait, "~", paste(combination, collapse = "*")))
        } else {
            formula <- as.formula(paste(trait, "~", paste(combination, collapse = "+")))
        }

        # Check if the trait is binary (0/1) to decide which family to use
        if (all(data_clean[[trait]] %in% c(0, 1))) {
            # Use binomial family for binary traits
            model <- glm(formula, data = data_clean, family = binomial(link = "logit"))
        } else {
            # Use Gamma family for continuous traits
            model <- glm(formula, data = data_clean, family = Gamma(link = "inverse"))
        }

        # Add the model to the list
        glm_models[[paste(trait, paste(combination, collapse = "_"), sep = ":")]] <- model
    }
}

# Initialize an empty dataframe to store the results
results_glm <- data.frame(Model = character(), AIC = numeric(), BIC = numeric(), stringsAsFactors = FALSE)

# Loop over each model in the glm_models list
for (model_name in names(glm_models)) {
    # Get the current model
    model <- glm_models[[model_name]]
    
    # Calculate AIC and BIC
    aic <- AIC(model)
    bic <- BIC(model)
    
    # Add the results to the dataframe
    results_glm <- rbind(results_glm, data.frame(Model = model_name, AIC = aic, BIC = bic, stringsAsFactors = FALSE))
}

# Split the "Model" column by ":"
split_models <- strsplit(results_glm$Model, ":")

# Get the "Response" and "Predictor" columns
results_glm$Response <- sapply(split_models, `[`, 1)
results_glm$Predictor <- sapply(split_models, `[`, 2)

# Print the results
print(results_glm)

# Find the minimum AIC for each "Response"
min_aic <- aggregate(AIC ~ Response, data = results_glm, FUN = min)

# Merge the minimum AIC values with the original dataframe
results_glm_min_aic <- merge(results_glm, min_aic, by = c("Response", "AIC"))

# Print the results
print(results_glm_min_aic)

tab_model(aov(glm_models[["D2Flowering:age_Final_Density_ID"]]),aov(glm_models[["flower_number_1:age_Final_Density_ID"]]),aov(glm_models[["pollen_avg:age_Final_Density"]]),aov(glm_models[["yield:age_Final_Density_ID"]]),show.aic = TRUE)
###############################################################
##Pollen Trade-Offs between Fittest and Least Fit Individuals##
###############################################################

# Create a list of predictors to include in the model
predictors <- c("age", "F_LF", "flower_number_1", "D2Flowering","Final_Density")

# Remove rows with NA in any of the specified columns
data_smalls <- data[complete.cases(data[, predictors]), ]

# Initialize a list to store the models
glm_models_pollen <- list()

# Remove NA values from "pollen_avg"
data_smalls <- data_smalls[!is.na(data_smalls$pollen_avg), ]

# Loop over the predictors
for (i in 1:length(predictors)) {
    # Create a formula for the model
    formula <- as.formula(paste("pollen_avg ~", predictors[i]))

    # Fit the GLM model
    model <- glm(formula, data = data_smalls, family = Gamma(link = "log"))

    # Add the model to the list
    glm_models_pollen[[paste(predictors[i], sep = "_")]] <- model
}

# Loop over the predictors
for (i in 1:(length(predictors) - 1)) {
    for (j in (i + 1):length(predictors)) {
        # Create a formula for the model with independent effects and interaction
        predictors_subset <- c(predictors[i], predictors[j])
        interactions_two_way <- combn(predictors_subset, 2)
        interaction_terms_two_way <- apply(interactions_two_way, 2, function(x) paste(x, collapse = ":"))
        formula <- as.formula(paste("pollen_avg ~", predictors[i], "+", predictors[j], "+", paste(interaction_terms_two_way, collapse = "+")))
        
        # Fit the GLM model
        model <- glm(formula, data = data_smalls, family = Gamma(link = "log"))
        
        # Add the model to the list
        glm_models_pollen[[paste(predictors[i], predictors[j], sep = ":")]] <- model
    }
}

# Loop over the predictors
for (i in 1:(length(predictors) - 2)) {
    for (j in (i + 1):(length(predictors) - 1)) {
        for (k in (j + 1):length(predictors)) {
            # Create a formula for the model with independent effects and interaction
            predictors_subset <- c(predictors[i], predictors[j], predictors[k])
            interactions_two_way <- combn(predictors_subset, 2)
            interactions_three_way <- combn(predictors_subset, 3)
            interaction_terms_two_way <- apply(interactions_two_way, 2, function(x) paste(x, collapse = ":"))
            interaction_terms_three_way <- apply(interactions_three_way, 2, function(x) paste(x, collapse = ":"))
            formula <- as.formula(paste("pollen_avg ~", predictors[i], "+", predictors[j], "+", predictors[k], "+", paste(interaction_terms_two_way, collapse = "+"), "+", paste(interaction_terms_three_way, collapse = "+")))
            
            # Fit the GLM model
            model <- glm(formula, data = data_smalls, family = Gamma(link = "log"))
            
            # Add the model to the list
            glm_models_pollen[[paste(predictors[i], predictors[j], predictors[k], sep = ":")]] <- model
        }
    }
}

# Loop over the predictors
for (i in 1:(length(predictors)-3)) {
    for (j in (i+1):(length(predictors)-2)) {
        for (k in (j+1):(length(predictors)-1)) {
            for (l in (k+1):length(predictors)) {
                # Create a formula for the model with independent effects and interaction
                predictors_subset <- c(predictors[i], predictors[j], predictors[k], predictors[l])
                interactions_two_way <- combn(predictors_subset, 2)
                interactions_three_way <- combn(predictors_subset, 3)
                interaction_terms_two_way <- apply(interactions_two_way, 2, function(x) paste(x, collapse = ":"))
                interaction_terms_three_way <- apply(interactions_three_way, 2, function(x) paste(x, collapse = ":"))
                formula <- as.formula(paste("pollen_avg ~", predictors[i], "+", predictors[j], "+", predictors[k], "+", predictors[l], "+", paste(interaction_terms_two_way, collapse = "+"), "+", paste(interaction_terms_three_way, collapse = "+")))
                
                # Fit the GLM model
                model <- glm(formula, data = data_smalls, family = Gamma(link = "log"))
                
                # Add the model to the list
                glm_models_pollen[[paste(predictors[i], predictors[j], predictors[k], predictors[l], sep = ":")]] <- model
            }
        }
    }
}

# Loop over the predictors
for (i in 1:(length(predictors)-4)) {
    for (j in (i+1):(length(predictors)-3)) {
        for (k in (j+1):(length(predictors)-2)) {
            for (l in (k+1):(length(predictors)-1)) {
                for (m in (l+1):length(predictors)) {
                    # Create a formula for the model with independent effects and interaction
                    predictors_subset <- c(predictors[i], predictors[j], predictors[k], predictors[l], predictors[m])
                    interactions_two_way <- combn(predictors_subset, 2)
                    interactions_three_way <- combn(predictors_subset, 3)
                    interaction_terms_two_way <- apply(interactions_two_way, 2, function(x) paste(x, collapse = ":"))
                    interaction_terms_three_way <- apply(interactions_three_way, 2, function(x) paste(x, collapse = ":"))
                    formula <- as.formula(paste("pollen_avg ~", predictors[i], "+", predictors[j], "+", predictors[k], "+", predictors[l], "+", predictors[m], "+", paste(interaction_terms_two_way, collapse = "+"), "+", paste(interaction_terms_three_way, collapse = "+")))
                    
                    # Fit the GLM model
                    model <- glm(formula, data = data_smalls, family = Gamma(link = "log"))
                    
                    # Add the model to the list
                    glm_models_pollen[[paste(predictors[i], predictors[j], predictors[k], predictors[l], predictors[m], sep = ":")]] <- model
                }
            }
        }
    }
}
# Initialize an empty dataframe
results_glm_pollen <- data.frame(Model = character(), AIC = numeric(), BIC = numeric(), stringsAsFactors = FALSE)

# Perform the tests on the models
for (i in names(glm_models_pollen)) {
# Calculate AIC and BIC
    aic <- AIC(glm_models_pollen[[i]])
    bic <- BIC(glm_models_pollen[[i]])
# Add the results to the data_smallsframe
    results_glm_pollen <- rbind(results_glm_pollen, data.frame(Model = i, AIC = aic, BIC = bic, stringsAsFactors = FALSE))
}

# Print the results
print(results_glm_pollen)

best<-subset(results_glm_pollen, AIC == min(AIC))[["Model"]]

tab_model(aov(glm_models_pollen[[best]]),show.aic = TRUE)
##########################################################################
## Flowering Production Trade-Offs between Fittestt and Least Fit Plants##
##########################################################################

# Create a list of predictors to include in the model
predictors <- c("age", "F_LF", "pollen_avg", "D2Flowering","Final_Density")

# Remove rows with NA in any of the specified columns
data_smalls <- data[complete.cases(data[, predictors]), ]

# Initialize a list to store the models
glm_models_flower <- list()

# Remove NA values from "flower_number_1"
data_smalls <- data[!is.na(data$flower_number_1), ]

# Loop over the predictors
for (i in 1:length(predictors)) {
    # Create a formula for the model
    formula <- as.formula(paste("flower_number_1 ~", predictors[i]))

    # Fit the GLM model
    model <- glm(formula, data = data_smalls, family = Gamma(link = "log"))

    # Add the model to the list
    glm_models_flower[[paste(predictors[i], sep = ":")]] <- model
}

# Loop over the predictors
for (i in 1:(length(predictors) - 1)) {
    for (j in (i + 1):length(predictors)) {
        # Create a formula for the model with independent effects and interaction
        predictors_subset <- c(predictors[i], predictors[j])
        interactions_two_way <- combn(predictors_subset, 2)
        interaction_terms_two_way <- apply(interactions_two_way, 2, function(x) paste(x, collapse = ":"))
        formula <- as.formula(paste("flower_number_1 ~", predictors[i], "+", predictors[j], "+", paste(interaction_terms_two_way, collapse = "+")))
        
        # Fit the GLM model
        model <- glm(formula, data = data_smalls, family = Gamma(link = "log"))
        
        # Add the model to the list
        glm_models_flower[[paste(predictors[i], predictors[j], sep = ":")]] <- model
    }
}

# Loop over the predictors
for (i in 1:(length(predictors) - 2)) {
    for (j in (i + 1):(length(predictors) - 1)) {
        for (k in (j + 1):length(predictors)) {
            # Create a formula for the model with independent effects and interaction
            predictors_subset <- c(predictors[i], predictors[j], predictors[k])
            interactions_two_way <- combn(predictors_subset, 2)
            interactions_three_way <- combn(predictors_subset, 3)
            interaction_terms_two_way <- apply(interactions_two_way, 2, function(x) paste(x, collapse = ":"))
            interaction_terms_three_way <- apply(interactions_three_way, 2, function(x) paste(x, collapse = ":"))
            formula <- as.formula(paste("flower_number_1 ~", predictors[i], "+", predictors[j], "+", predictors[k], "+", paste(interaction_terms_two_way, collapse = "+"), "+", paste(interaction_terms_three_way, collapse = "+")))
            
            # Fit the GLM model
            model <- glm(formula, data = data_smalls, family = Gamma(link = "log"))
            
            # Add the model to the list
            glm_models_flower[[paste(predictors[i], predictors[j], predictors[k], sep = ":")]] <- model
        }
    }
}

# Loop over the predictors
for (i in 1:(length(predictors)-3)) {
    for (j in (i+1):(length(predictors)-2)) {
        for (k in (j+1):(length(predictors)-1)) {
            for (l in (k+1):length(predictors)) {
                # Create a formula for the model with independent effects and interaction
                predictors_subset <- c(predictors[i], predictors[j], predictors[k], predictors[l])
                interactions_two_way <- combn(predictors_subset, 2)
                interactions_three_way <- combn(predictors_subset, 3)
                interaction_terms_two_way <- apply(interactions_two_way, 2, function(x) paste(x, collapse = ":"))
                interaction_terms_three_way <- apply(interactions_three_way, 2, function(x) paste(x, collapse = ":"))
                formula <- as.formula(paste("flower_number_1 ~", predictors[i], "+", predictors[j], "+", predictors[k], "+", predictors[l], "+", paste(interaction_terms_two_way, collapse = "+"), "+", paste(interaction_terms_three_way, collapse = "+")))
                
                # Fit the GLM model
                model <- glm(formula, data = data_smalls, family = Gamma(link = "log"))
                
                # Add the model to the list
                glm_models_flower[[paste(predictors[i], predictors[j], predictors[k], predictors[l], sep = ":")]] <- model
            }
        }
    }
}

# Loop over the predictors
for (i in 1:(length(predictors)-4)) {
    for (j in (i+1):(length(predictors)-3)) {
        for (k in (j+1):(length(predictors)-2)) {
            for (l in (k+1):(length(predictors)-1)) {
                for (m in (l+1):length(predictors)) {
                    # Create a formula for the model with independent effects and interaction
                    predictors_subset <- c(predictors[i], predictors[j], predictors[k], predictors[l], predictors[m])
                    interactions_two_way <- combn(predictors_subset, 2)
                    interactions_three_way <- combn(predictors_subset, 3)
                    interaction_terms_two_way <- apply(interactions_two_way, 2, function(x) paste(x, collapse = ":"))
                    interaction_terms_three_way <- apply(interactions_three_way, 2, function(x) paste(x, collapse = ":"))
                    formula <- as.formula(paste("flower_number_1 ~", predictors[i], "+", predictors[j], "+", predictors[k], "+", predictors[l], "+", predictors[m], "+", paste(interaction_terms_two_way, collapse = "+"), "+", paste(interaction_terms_three_way, collapse = "+")))
                    
                    # Fit the GLM model
                    model <- glm(formula, data = data_smalls, family = Gamma(link = "log"))
                    
                    # Add the model to the list
                    glm_models_flower[[paste(predictors[i], predictors[j], predictors[k], predictors[l], predictors[m], sep = ":")]] <- model
                }
            }
        }
    }
}
# Initialize an empty data_smallsframe
results_glm_flower <- data.frame(Model = integer(), AIC = numeric(), BIC = numeric(), stringsAsFactors = FALSE)

# Perform the tests on the models
for (i in seq_along(glm_models_flower)) {
# Calculate AIC and BIC
    aic <- AIC(glm_models_flower[[i]])
    bic <- BIC(glm_models_flower[[i]])
# Add the results to the data_smallsframe
    results_glm_flower <- rbind(results_glm_flower, data.frame(Model = i, AIC = aic, BIC = bic, stringsAsFactors = FALSE))
}

# Print the results
print(results_glm_flower)

best<-subset(results_glm_flower, AIC == min(AIC))[["Model"]]

tab_model(aov(glm_models_flower[[best]]), show.aic = TRUE)

######################################################################################
## Take the Least Square Means and Standard Error for Predictors from the Best Model##
######################################################################################

