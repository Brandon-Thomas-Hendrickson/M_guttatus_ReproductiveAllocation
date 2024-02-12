# Install the package
install.packages(c("ggplot2","geepack", "readr", "sjPlot", "MASS", "stats","emmeans","wesanderson","ggpubr","purrr"))

# Load the package
library(readr)
library(sjPlot)
library(MASS)
library(stats)
library(emmeans)
library(wesanderson)
library(ggpubr)
library(ggplot2)
library(purrr)

############################################################################################################
## General Response of Fitness, Phenology, and Reproductive Allocation to Growing Length and Competition ###
############################################################################################################

#Load modelframe
data <- read.csv("/Users/brandonhendrickson/Documents/Github_Projects/M_guttatus_ReproductiveAllocation/Data/csv/raw_data.csv")

# Convert Final_Density to a factor
data$Final_Density <- as.factor(data$Final_Density)

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
predictors <- c("age", "Final_Density")

# Generate all combinations of predictors
combinations <- lapply(1:length(predictors), function(x) combn(predictors, x, simplify = FALSE))

# Flatten the list
combinations <- unlist(combinations, recursive = FALSE)

# Loop over the trait list
for (trait in trait_list) {
    # Loop over the combinations of predictors
    for (predictor_combination in combinations) {
        # Create the formula with interaction
        formula <- as.formula(paste(trait, "~", paste(predictor_combination, collapse = " * ")))
        
        # Determine the family based on the trait
        if (trait %in% c("yield", "pollen_avg")) {
            family <- gaussian
            formula <- as.formula(paste("log(", trait, ") ~", paste(predictor_combination, collapse = " * ")))
        } else if (trait %in% c("D2Flowering", "flower_number_1")) {
            family <- Gamma(link = "log")
        }
        
        # Fit the model
        model <- glm(formula, data = data, family = family)
        
        # Store the model in the list
        glm_models[[paste(trait, paste(predictor_combination, collapse = ":"), sep = ":")]] <- model
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

for (i in 1:length(results_glm_min_aic)) {
    print(summary(aov(glm_models[[results_glm_min_aic$Model[i]]])))
}

# Initialize an empty dataframe to store the results
results_df <- data.frame(Predictor = character(), Trait = character(), Estimate = numeric(), Pr_F = numeric(), F_value = numeric(), Sum_Sq = numeric(), Mean_Sq = numeric(), stringsAsFactors = FALSE)

# Loop over each model in the models list
for(trait in results_glm_min_aic$Model) {
    # Get the model
    model <- glm_models[[trait]]
    
    # Get the ANOVA table
    anova_table <- summary(aov(model))
    
    # Loop over each coefficient and row in the ANOVA table
    for (i in 2:nrow(anova_table[[1]])) {
        # Get the estimate
        estimate <- coef(model)[i][[1]]
        
        # Get the "Pr(>F)", "F value", "Sum Sq", and "Mean Sq"
        pr_f <- anova_table[[1]][i-1, "Pr(>F)"]
        f_value <- anova_table[[1]][i-1, "F value"]
        sum_sq <- anova_table[[1]][i-1, "Sum Sq"]
        mean_sq <- anova_table[[1]][i - 1, "Mean Sq"]
        predictor <- rownames(anova_table[[1]])[i-1]
        
        # Add the result to the results dataframe
        results_df <- rbind(results_df, data.frame(Predictor = predictor, Trait = trait, Estimate = estimate, Pr_F = pr_f, F_value = f_value, Sum_Sq = sum_sq, Mean_Sq = mean_sq, stringsAsFactors = FALSE))
    }
}
###############################################################
##Pollen Trade-Offs between Fittest and Least Fit Individuals##
###############################################################
#Load modelframe
data <- read.csv("/Users/brandonhendrickson/Documents/Github_Projects/M_guttatus_ReproductiveAllocation/Data/csv/raw_data.csv")

# Convert Final_Density to a factor
data$Final_Density <- as.factor(data$Final_Density)

#Identify numeric columns
numeric_cols <- sapply(data, is.numeric)

#Replace 0 with 0.000001
data[numeric_cols] <- lapply(data[numeric_cols], function(x) {
  ifelse(x == 0, 0.000001, x)
})
# Create a list of predictors to include in the model
predictors <- c("age", "F_LF", "flower_number_1","Final_Density")

# Generate all combinations of predictors
combinations <- lapply(1:length(predictors), function(x) combn(predictors, x, simplify = FALSE))

# Flatten the list
combinations <- unlist(combinations, recursive = FALSE)

# Remove rows with NA in any of the specified columns
data_smalls <- data[complete.cases(data[, predictors]), ]

# Initialize a list to store the models
glm_models_pollen <- list()

# Remove NA values from "pollen_avg"
data_smalls <- data_smalls[!is.na(data_smalls$pollen_avg), ]

for (predictor_combination in combinations) {
    # Determine the family based on the trait
    family <- gaussian

    # Check if more than one predictor is used
    if (length(predictor_combination) > 1) {
        # Create the formula with all two-way and three-way interactions
        formula <- as.formula(paste("log(pollen_avg) ~ (", paste(predictor_combination, collapse = " + "), ")^3"))
    } else {
        # Create the formula without interactions
        formula <- as.formula(paste("log(pollen_avg) ~", paste(predictor_combination, collapse = " + ")))
    }

    # Fit the model
    model <- glm(formula, data = data_smalls, family = family)

    # Store the model in the list
    glm_models_pollen[[paste(predictor_combination, collapse = ":")]] <- model
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

best_pollen<-subset(results_glm_pollen, AIC == min(AIC))[["Model"]]

tab_model(aov(glm_models_pollen[[best_pollen]]),show.aic = TRUE)
##########################################################################
## Flowering Production Trade-Offs between Fittest and Least Fit Plants ##
##########################################################################
#Load modelframe
data <- read.csv("/Users/brandonhendrickson/Documents/Github_Projects/M_guttatus_ReproductiveAllocation/Data/csv/raw_data.csv")

# Convert Final_Density to a factor
data$Final_Density <- as.factor(data$Final_Density)

#Identify numeric columns
numeric_cols <- sapply(data, is.numeric)

#Replace 0 with 0.000001
data[numeric_cols] <- lapply(data[numeric_cols], function(x) {
  ifelse(x == 0, 0.000001, x)
})
# Create a list of predictors to include in the model
predictors <- c("age", "F_LF", "pollen_avg","Final_Density")

# Generate all combinations of predictors
combinations <- lapply(1:length(predictors), function(x) combn(predictors, x, simplify = FALSE))

# Flatten the list
combinations <- unlist(combinations, recursive = FALSE)

# Remove rows with NA in any of the specified columns
data_smalls <- data[complete.cases(data[, predictors]), ]

# Initialize a list to store the models
glm_models_flower <- list()

# Remove NA values from "pollen_avg"
data_smalls <- data_smalls[!is.na(data_smalls$flower_number_1), ]

for (predictor_combination in combinations) {
    # Determine the family based on the trait
    family <- gaussian

    # Check if more than one predictor is used
    if (length(predictor_combination) > 1) {
        # Create the formula with all two-way and three-way interactions
        formula <- as.formula(paste("flower_number_1 ~ (", paste(predictor_combination, collapse = " + "), ")^3"))
    } else {
        # Create the formula without interactions
        formula <- as.formula(paste("flower_number_1 ~", paste(predictor_combination, collapse = " + ")))
    }

    # Fit the model
    model <- glm(formula, data = data_smalls, family = Gamma)

    # Store the model in the list
    glm_models_flower[[paste(predictor_combination, collapse = ":")]] <- model
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

best_flower<-subset(results_glm_flower, AIC == min(AIC))[["Model"]]

tab_model(aov(glm_models_flower[[best_flower]]), show.aic = TRUE)

############################################################
## Day to First Flowering by Age, F_LF, and Final Density ##
############################################################

#Load modelframe
data <- read.csv("/Users/brandonhendrickson/Documents/Github_Projects/M_guttatus_ReproductiveAllocation/Data/csv/raw_data.csv")

# Convert Final_Density to a factor
data$Final_Density <- as.factor(data$Final_Density)

#Identify numeric columns
numeric_cols <- sapply(data, is.numeric)

#Replace 0 with 0.000001
data[numeric_cols] <- lapply(data[numeric_cols], function(x) {
  ifelse(x == 0, 0.000001, x)
})
# Create a list of predictors to include in the model
predictors <- c("age", "F_LF", "Final_Density")

# Generate all combinations of predictors
combinations <- lapply(1:length(predictors), function(x) combn(predictors, x, simplify = FALSE))

# Flatten the list
combinations <- unlist(combinations, recursive = FALSE)

# Remove rows with NA in any of the specified columns
data_smalls <- data[complete.cases(data[, predictors]), ]

# Initialize a list to store the models
glm_models_D2Flowering <- list()

# Remove NA values from "pollen_avg"
data_smalls <- data_smalls[!is.na(data_smalls$D2Flowering), ]

for (predictor_combination in combinations) {
    # Determine the family based on the trait
    family <- gaussian

    # Check if more than one predictor is used
    if (length(predictor_combination) > 1) {
        # Create the formula with all two-way and three-way interactions
        formula <- as.formula(paste("D2Flowering ~ (", paste(predictor_combination, collapse = " + "), ")^3"))
    } else {
        # Create the formula without interactions
        formula <- as.formula(paste("D2Flowering ~", paste(predictor_combination, collapse = " + ")))
    }

    # Fit the model
    model <- glm(formula, data = data_smalls, family = family)

    # Store the model in the list
    glm_models_D2Flowering[[paste(predictor_combination, collapse = ":")]] <- model
}
# Initialize an empty data_smallsframe
results_glm_D2Flowering <- data.frame(Model = integer(), AIC = numeric(), BIC = numeric(), stringsAsFactors = FALSE)

# Perform the tests on the models
for (i in seq_along(glm_models_D2Flowering)) {
# Calculate AIC and BIC
    aic <- AIC(glm_models_D2Flowering[[i]])
    bic <- BIC(glm_models_D2Flowering[[i]])
# Add the results to the data_smallsframe
    results_glm_D2Flowering <- rbind(results_glm_D2Flowering, data.frame(Model = i, AIC = aic, BIC = bic, stringsAsFactors = FALSE))
}

# Print the results
print(results_glm_D2Flowering)

best_D2Flowering<-subset(results_glm_D2Flowering, AIC == min(AIC))[["Model"]]

tab_model(aov(glm_models_D2Flowering[[best_D2Flowering]]), show.aic = TRUE)

######################################################################################
## Take the Least Square Means and Standard Error for Predictors from the Best Model##
######################################################################################
get_means <- function(model_name) {
    lsmeans <- emmeans(glm_models[[model_name]], ~ age + Final_Density)
    means <- as.data.frame(summary(lsmeans))
    return(means)
}

means_D2Flowering <- get_means("D2Flowering:age:Final_Density")
means_flower_number_1 <- get_means("flower_number_1:age:Final_Density")
means_yield <- get_means("yield:age:Final_Density")
means_pollen_avg <- get_means("pollen_avg:age:Final_Density")



p4 <- ggplot(means_D2Flowering, aes(x = Final_Density, y = emmean, color = age, group = age)) +
    geom_line(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.1, position = position_dodge(width = 0.5)) +
    geom_point(position = position_dodge(width = 0.5)) +
    scale_color_manual(values = wes_palette("FantasticFox1")) +
    labs(
        x = "Final Density",
        y = "Days to Flowering"
    ) +
    theme_light() +
    theme(legend.position = "none",
          text = element_text(family = "Times", size = 12),
          axis.text = element_text(size = 12))

p2 <- ggplot(means_flower_number_1, aes(x = Final_Density, y = emmean, color = age, group = age)) +
    geom_line(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.1, position = position_dodge(width = 0.5)) +
    geom_point(position = position_dodge(width = 0.5)) +
    scale_color_manual(values = wes_palette("FantasticFox1"),labels = c("1.E" = "0", "2.M" = "15", "3.L" = "30")) +
    labs(
        x = "Final Density",
        y = "Flower Number"
    ) +
    theme_light() +
    theme(text = element_text(family = "Times", size = 12),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 12),
          legend.position = c(0.4,0.8), 
          legend.justification = c(0,0)) +
    guides(color = guide_legend(title = NULL, direction = "horizontal"))

p3 <- ggplot(means_yield, aes(x = Final_Density, y = emmean, color = age, group = age)) +
    geom_line(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.1, position = position_dodge(width = 0.5)) +
    geom_point(position = position_dodge(width = 0.5)) +
    scale_color_manual(values = wes_palette("FantasticFox1")) +
    labs(
        x = "Final Density",
        y = "Seed Mass"
    ) +
    theme_light() +
    theme(legend.position = "none",
          text = element_text(family = "Times", size = 12),
          axis.text = element_text(size = 12))

p1 <- ggplot(means_pollen_avg, aes(x = Final_Density, y = emmean, color = age, group = age)) +
    geom_line(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.1, position = position_dodge(width = 0.5)) +
    geom_point(position = position_dodge(width = 0.5)) +
    scale_color_manual(values = wes_palette("FantasticFox1")) +
    labs(
        x = "Final Density",
        y = "Pollen Number"
    ) +
    theme_light() +
    theme(legend.position = "none",
          text = element_text(family = "Times", size = 12),
          axis.text = element_text(size = 12))

# Arrange the plots
plots <- ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)

# Save the plots as a PDF
ggsave("traits_finalDensity_age.pdf", plots, width = 8, height = 6, units = "in")


    lsmeans <- emmeans(glm_models_pollen[["age:F_LF:flower_number_1:Final_Density"]], ~ age + Final_Density + F_LF + flower_number_1)
    means <- as.data.frame(summary(lsmeans))
    
    # Calculate the average emmean and SE for each combination of age and Final_Density
    average_means <- aggregate(cbind(emmean, SE) ~ age + Final_Density + F_LF + flower_number_1, data = means, FUN = mean)
p1 <- ggplot(average_means, aes(x = Final_Density, y = emmean, color = F_LF, group = F_LF)) +
    geom_line(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.1, position = position_dodge(width = 0.5)) +
    geom_point(position = position_dodge(width = 0.5)) +
    scale_color_manual(values = wes_palette("FantasticFox1"),labels = c("F" = "Fittest", "LF" = "Least Fit","Other"="Other")) +
    labs(
        x = "Final Density",
        y = "Pollen Number"
    ) +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5),
          strip.text = element_blank(),
          legend.position = c(0.98, 0.65),
          legend.justification = c(1, 0),
          text = element_text(family = "Times", size = 12),
          axis.text = element_text(size = 12),
          legend.text = element_text(size = 12)) +
    guides(color = guide_legend(title = NULL, direction = "horizontal")) + 
    facet_wrap(~age)

    lsmeans <- emmeans(glm_models_flower[[15]], ~ age + Final_Density + F_LF + pollen_avg)
    means <- as.data.frame(summary(lsmeans))
    
    # Calculate the average emmean and SE for each combination of age and Final_Density
    average_means <- aggregate(cbind(emmean, SE) ~ age + Final_Density + F_LF + pollen_avg, data = means, FUN = mean)

p2 <- ggplot(average_means, aes(x = Final_Density, y = emmean, color = F_LF, group = F_LF)) +
    geom_line(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.1, position = position_dodge(width = 0.5)) +
    geom_point(position = position_dodge(width = 0.5)) +
    scale_color_manual(values = wes_palette("FantasticFox1"), labels = c("F" = "Fittest", "LF" = "Least Fit", "Other" = "Other")) +
    labs(
        x = "Final Density",
        y = "Flower Number"
    ) +
    theme_light() +
    theme(
        plot.title = element_text(hjust = 0.5),
        strip.text = element_blank(),
        legend.position = "none",
        text = element_text(family = "Times", size = 12),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12)
    ) +
    guides(color = guide_legend(title = NULL)) +
    facet_wrap(~age)

    lsmeans <- emmeans(glm_models_D2Flowering[[7]], ~ age + Final_Density + F_LF)
    means <- as.data.frame(summary(lsmeans))
    
    # Calculate the average emmean and SE for each combination of age and Final_Density
    average_means <- aggregate(cbind(emmean, SE) ~ age + Final_Density + F_LF, data = means, FUN = mean)

p3 <- ggplot(average_means, aes(x = Final_Density, y = emmean, color = F_LF, group = F_LF)) +
    geom_line(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.1, position = position_dodge(width = 0.5)) +
    geom_point(position = position_dodge(width = 0.5)) +
    scale_color_manual(values = wes_palette("FantasticFox1"), labels = c("F" = "Fittest", "LF" = "Least Fit", "Other" = "Other")) +
    labs(
        x = "Final Density",
        y = "Days to Flowering"
    ) +
    theme_light() +
    theme(
        plot.title = element_text(hjust = 0.5),
        strip.text = element_blank(),
        legend.position = "none",
        text = element_text(family = "Times", size = 12),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12)
    ) +
    guides(color = guide_legend(title = NULL)) +
    facet_wrap(~age)

plots <- ggarrange(p1, p2, p3, nrow = 3)

# Save the plots as a PDF
ggsave("traits_finalDensity_age_Fitness.pdf", plots, width = 10, height = 6, units = "in")