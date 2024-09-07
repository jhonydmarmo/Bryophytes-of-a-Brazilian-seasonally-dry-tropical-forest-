################################################## BOUNDING BOX #########################################################

# Install and load the necessary package
if (!requireNamespace("sf", quietly = TRUE)) install.packages("sf")
library(sf) # For spatial data manipulation

# Define the path to the shapefile
shapefile_path <- "C:\\path_to_directory\\caatinga.shp"

# Read the Caatinga shapefile
caatinga <- st_read(shapefile_path)

# Calculate the bounding box of the shapefile
bbox <- st_bbox(caatinga)

# Print the bounding box
print(bbox)

################################################## COORDINATES FILTERING #########################################################

# Install necessary packages if not already installed
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("sf", quietly = TRUE)) install.packages("sf")
if (!requireNamespace("writexl", quietly = TRUE)) install.packages("writexl")

# Load necessary libraries
library(readxl)    # For reading Excel files
library(sf)        # For spatial data manipulation
library(writexl)   # For saving Excel files

# Define file paths
mosses_path <- "C:\\path_to_directory\\bryophyta.xlsx"
liverworts_path <- "C:\\path_to_directory\\marchantiophyta.xlsx"
hornworts_path <- "C:\\path_to_directory\\anthocerotophyta.xlsx"

# Read the Excel files with the coordinates
mosses <- read_excel(mosses_path)
liverworts <- read_excel(liverworts_path)
hornworts <- read_excel(hornworts_path)

# Convert the coordinates to an sf object, keeping the CRS of the shapefile
coordinates_mosses <- st_as_sf(mosses, coords = c("decimalLongitude", "decimalLatitude"), crs = st_crs(caatinga))
coordinates_liverworts <- st_as_sf(liverworts, coords = c("decimalLongitude", "decimalLatitude"), crs = st_crs(caatinga))
coordinates_hornworts <- st_as_sf(hornworts, coords = c("decimalLongitude", "decimalLatitude"), crs = st_crs(caatinga))

# Perform the filtering: check which coordinates are within the Caatinga shapefile for each group
indices_inside_mosses <- st_within(coordinates_mosses, caatinga, sparse = FALSE)
indices_inside_liverworts <- st_within(coordinates_liverworts, caatinga, sparse = FALSE)
indices_inside_hornworts <- st_within(coordinates_hornworts, caatinga, sparse = FALSE)

# Filter the rows of the original dataframes based on the indices
coordinates_inside_mosses <- mosses[unlist(indices_inside_mosses), ]
coordinates_inside_liverworts <- liverworts[unlist(indices_inside_liverworts), ]
coordinates_inside_hornworts <- hornworts[unlist(indices_inside_hornworts), ]

# Check the results
print(coordinates_inside_mosses)
print(coordinates_inside_liverworts)
print(coordinates_inside_hornworts)

# Save the filtered coordinates to Excel files
write_xlsx(coordinates_inside_mosses, "filtered_coordinates_mosses.xlsx")
write_xlsx(coordinates_inside_liverworts, "filtered_coordinates_liverworts.xlsx")
write_xlsx(coordinates_inside_hornworts, "filtered_coordinates_hornworts.xlsx")



########################################## TAXONOMIC NAME UPDATE ##############################################

# Install the U.Taxonstand package (if not already installed)
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")  # Package for installing R packages from GitHub
}
devtools::install_github("ecoinfor/U.Taxonstand")  # Install the U.Taxonstand package from GitHub

# Load necessary packages
library(U.Taxonstand)  # For standardizing species names and taxonomic matching
library(openxlsx)      # For reading from and writing to Excel files

# Define file paths
database_path <- "path_to_your_database/databaseb.xlsx"
species_path <- "path_to_your_species_list/species.xlsx"

# Load your taxonomy database
database <- read.xlsx(database_path)  # Read the taxonomy database from an Excel file

# Load your species list
species <- read.xlsx(species_path)    # Read the species list from an Excel file

# Run the main function to standardize species names
result <- nameMatch(spList = species, spSource = database, 
                    author = TRUE, genusPairs = NULL, Append = FALSE)  # Standardize and match species names

# Save the results to an Excel file
write.xlsx(result, "Standardized_Taxonomy_Results.xlsx", overwrite = TRUE)  # Save the results to an Excel file


######################################### CHAO1 RICHNESS ESTIMATION FOR TOTAL SPECIES NUMBER ###########################

######################################### CHAO1 RICHNESS ESTIMATION FOR TOTAL SPECIES NUMBER ###########################

# Load necessary libraries
library(vegan)   # For community ecology analyses, including richness estimators like Chao1
library(ggplot2) # For data visualization and plotting

# Read the community matrix from a CSV file
community_matrix <- read.csv("C:\\path_to_directory\\community_matrix.csv", row.names = 1)
community_matrix <- community_matrix[, -1]  # Remove the first column if it is not needed for analysis

# Estimate Chao1 richness for the entire community matrix
# `estaccumR` estimates species richness and performs accumulation curves
est_chao1 <- estaccumR(community_matrix, permutations = 100)  # Perform 100 permutations for estimation

# Visualize the results
options(max.print = 10000)  # Set the maximum print option to view large outputs
summary(est_chao1, display = "chao")  # Display summary statistics of Chao1 estimates
str(est_chao1)  # Display the structure of the Chao1 estimate object

# Prepare data for plotting
# Extract Chao1 and observed richness estimates from the summary
results <- summary(est_chao1, display = c("S", "chao"))
res_chao <- cbind(results$chao[, 1:4], results$S[, 2:4])
res_chao <- as.data.frame(res_chao)  # Convert the results to a data frame
colnames(res_chao) <- c("Samples", "Chao", "C_lower", "C_upper", "Richness", "R_lower", "R_upper")

# Plot Chao1 richness estimates and observed richness
ggplot(res_chao, aes(y = Richness, x = Samples)) +
  geom_point(aes(y = Chao, x = Samples + 0.1), size = 4, 
             color = "darkorange", alpha = 0.7) +  # Plot Chao1 estimates
  geom_point(aes(y = Richness, x = Samples), size = 4, 
             color = "cyan4", alpha = 0.7) +  # Plot observed richness
  geom_label(y = 550, x = 255, label = "Estimated richness") +  # Label for estimated richness
  geom_label(y = 420, x = 255, label = "Observed richness") +   # Label for observed richness
  geom_line(aes(y = Chao, x = Samples), color = "darkorange") +  # Line for Chao1 estimates
  geom_line(aes(y = Richness, x = Samples), color = "cyan4") +   # Line for observed richness
  geom_linerange(aes(ymin = C_lower, ymax = C_upper, x = Samples + 0.1), color = "darkorange") +  # Error bars for Chao1
  geom_linerange(aes(ymin = R_lower, ymax = R_upper, x = Samples), color = "cyan4") +  # Error bars for observed richness
  scale_x_continuous(breaks = seq(0, max(res_chao$Samples, na.rm = TRUE), by = 70)) +  # Customize x-axis breaks
  labs(x = "Number of samples", y = "Estimated richness - Chao 1")  # Label axes

######################################### OBSERVED AND ESTIMATED RICHNESS BY LOCALITY ##########################################

# Install necessary package
install.packages("vegan")  # For community ecology analyses, including richness estimators

# Load necessary library
library(vegan) # Used for community ecology analyses, including richness estimators

# Estimate Chao1 richness for each locality
chao1_estimates <- estimateR(community_matrix, method = "chao")
print(chao1_estimates)  # Print the estimated richness
print(str(chao1_estimates))  # Print the structure of the estimates object
write.csv(chao1_estimates, file = "chao1_richness.csv", row.names = TRUE)  # Save the estimates to a CSV file

# Perform the Shapiro-Wilk test for observed richness
shapiro_result_observed <- shapiro.test(chao1_estimates["S.obs", ])
print(shapiro_result_observed)  # Print the result of the Shapiro-Wilk test for observed richness

# Perform the Shapiro-Wilk test for Chao1 estimated richness
shapiro_result_estimated_chao1 <- shapiro.test(chao1_estimates["S.chao1", ])
print(shapiro_result_estimated_chao1)  # Print the result of the Shapiro-Wilk test for Chao1 estimated richness

# Perform the Wilcoxon signed-rank test for paired samples
wilcoxon_result <- wilcox.test(x = chao1_estimates["S.obs", ], 
                               y = chao1_estimates["S.chao1", ],
                               paired = TRUE)
print(wilcoxon_result)  # Print the result of the Wilcoxon signed-rank test

##################################################### PCA ########################################################

# Loading required libraries
library(FactoMineR)  # For performing Principal Component Analysis (PCA)
library(factoextra)  # For visualization of PCA results

#### PCA Analysis ####

# Performing PCA on the environmental data
res.pca <- PCA(environmental_data, graph = FALSE)

# Retrieve and print eigenvalues to understand the variance explained by each principal component
eig <- get_eig(res.pca)
print(eig)

# Visualize variables with colors based on their contributions to the principal components
fviz_pca_var(res.pca, col.var = "contrib", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE)

# Create a biplot of individuals and variables with the same color specifications
fviz_pca_biplot(res.pca, col.var = "contrib", 
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
                repel = TRUE)

# Extract variable loadings on the principal components
variable_loadings <- get_pca_var(res.pca)$var$coord

# Calculate Spearman correlation between the original variables and the principal components
correlation_matrix <- cor(environmental_data, variable_loadings, method = "spearman")

# Display the correlation matrix to see how the original variables correlate with the principal components
print(correlation_matrix)

# Save the Spearman correlation matrix as a CSV file
write.csv(correlation_matrix, file = "Spearman_Correlation_Var.csv", row.names = TRUE)

######################################################### SAMPBIAS ########################################################
library(sampbias)    # For sampling bias analysis
library(sf)          # For reading and manipulating spatial data in shapefile format

# Reading occurrence data
species <- read.csv("./data/occurrence_data.csv", sep = ";", header = TRUE)

# Reading the Caatinga shapefile using sf
caatinga_sf <- st_read("G:\\.shortcut-targets-by-id\\1-eUDhwgIhXRRCbbUt_qa3tilbg97UlsY\\Jhonyd\\Mestrado\\Overview - artigo 1\\SHAPES\\Caatinga (shapefile) - TerrabrasilisINPE\\biome_border.shp")

# Preparing the sampling bias matrix
# Note: Ensure that 'sampbias_matrix' is correctly defined. If you need to create it, ensure it matches the required format for `calculate_bias`.

# Setting parameters for sampling bias calculation
bryophytes_samp <- calculate_bias(
  x = species,          # Assuming 'species' contains the occurrence data
  res = 0.0287,         # Resolution for the bias calculation
  buffer = 1,           # Buffer size
  verbose = TRUE,       # Print detailed output
  restrict_sample = caatinga_sf  # Restrict to the area defined by the shapefile
)

# Generating a statistical summary and visualizing the data
summary(bryophytes_samp)
plot(bryophytes_samp)

# Projecting and visualizing the effects of bias in space
proj <- project_bias(bryophytes_samp)
map_bias(proj, gaz = NULL, sealine = FALSE, type = "sampling_rate")

# Saving the graph of estimated bias weights
png("estimated_bias_weights.png", width = 2000, height = 2000, res = 500)
plot(bryophytes_samp)
dev.off()

# Saving visualization of bias effects in space
png("bias_effects_in_space.png", width = 2000, height = 2000, res = 500)
map_bias(proj, gaz = NULL, sealine = FALSE, type = "sampling_rate")
dev.off()


########################################### GENERALIZED LINEAR MIXED MODELS (GLMM) ######################################
# Install and load necessary packages
if (!requireNamespace("glmmTMB", quietly = TRUE)) {
  install.packages("glmmTMB")
}
if (!requireNamespace("MuMIn", quietly = TRUE)) {
  install.packages("MuMIn")
}
if (!requireNamespace("doParallel", quietly = TRUE)) {
  install.packages("doParallel")
}
if (!requireNamespace("foreach", quietly = TRUE)) {
  install.packages("foreach")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("lmtest", quietly = TRUE)) {
  install.packages("lmtest")
}
if (!requireNamespace("spdep", quietly = TRUE)) {
  install.packages("spdep")
}

library(glmmTMB)    # For fitting GLMMs
library(MuMIn)      # For model selection and comparison
library(doParallel) # For parallel processing
library(foreach)    # For parallel processing
library(dplyr)      # For data manipulation
library(lmtest)     # For diagnostic tests
library(spdep)      # For spatial dependence tests

# Load the dataset
data <- read.csv("path_to_data.csv")

# Scale the environmental variables
scaled_data <- data %>%
  mutate(across(c(bio2, bio4, bio15, bio18, elev, ai), scale))  # Scaling selected environmental variables

# Ensure that 'locality_id' is a factor
scaled_data$locality_id <- as.factor(scaled_data$locality_id)

# Define the explanatory variables
variables <- c("bio2", "bio4", "bio15", "bio18", "elev", "ai")

# Create all possible combinations of the explanatory variables
combinations <- unlist(lapply(1:length(variables), function(x) combn(variables, x, simplify = FALSE)), recursive = FALSE)

# Set up parallel processing
numCores <- detectCores() - 1
cl <- makeCluster(numCores)
registerDoParallel(cl)

# Function to fit a glmmTMB model with Negative Binomial family and BFGS optimizer
fit_model <- function(vars) {
  formula_current <- as.formula(paste("S.obs ~", paste(vars, collapse = " + "), "+ (1 | locality_id)"))
  model <- glmmTMB(
    formula_current, 
    family = nbinom2(), 
    data = scaled_data,
    verbose = TRUE,
    control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS", maxit = 1000, reltol = 1e-8))
  )
  return(model)
}

# Fit all possible models in parallel using foreach
all_models <- foreach(vars = combinations, .packages = c("glmmTMB")) %dopar% {
  fit_model(vars)
}

# Name the fitted models
names(all_models) <- paste("model", seq_along(all_models), sep = "_")

# Stop the cluster
stopCluster(cl)

# Model comparison using AIC with MuMIn and ordering by AICc
model_comparison <- model.sel(all_models)
ordered_model_comparison <- model_comparison[order(model_comparison$AICc), ]

# Display the best model based on AICc
best_model <- ordered_model_comparison[1, ]
print(best_model)

# Display the ordered model comparison table
print(ordered_model_comparison)

# Save model comparison results to a text file
output_file <- "best_models.txt"
output_text <- capture.output(print(ordered_model_comparison))
writeLines(output_text, con = output_file)
cat("Model comparison results have been saved to", output_file, "\n")

################################################## Best Model Analysis ################################################

# Extract the name of the best model based on AICc
best_model_name <- rownames(model_comparison)[1]
best_model <- all_models[[best_model_name]]

# Display summary of the best model
summary(best_model)

# Extract residuals from the best model
residuals_glmm <- residuals(best_model, type = "pearson")

# Shapiro-Wilk test for normality of residuals
shapiro_result <- shapiro.test(residuals_glmm)
print(shapiro_result)

# Breusch-Pagan test for heteroscedasticity
bptest_result <- bptest(best_model)
print(bptest_result)

# Create a spatial weights matrix based on geographic coordinates
coords <- data.frame(x = scaled_data$decimalLongitude, y = scaled_data$decimalLatitude)
nb <- knn2nb(knearneigh(coords, k = 4))  # 4 nearest neighbors
lw <- nb2listw(nb, style = "W")  # Spatial weights matrix

# Calculate Moran's I for model residuals
moran_result <- moran.test(residuals_glmm, lw)
print(moran_result)

# Calculate Geary's C for model residuals
geary_result <- geary.test(residuals_glmm, lw)
print(geary_result)

# Check for multicollinearity using VIF
vif_results <- vif(best_model)
print(vif_results)

########################################### PERMUTATIONAL MULTIVARIATE ANALYSIS OF VARIANCE (PERMANOVA) ##################

# Load necessary libraries
library(vegan) # Used for ecological data analysis, including PERMANOVA

# Scale environmental data variables to ensure comparability
environmental_data$bio2 <- scale(environmental_data$bio2)
environmental_data$bio4 <- scale(environmental_data$bio4)
environmental_data$bio15 <- scale(environmental_data$bio15)
environmental_data$elev <- scale(environmental_data$elev)
environmental_data$ia <- scale(environmental_data$ia)
environmental_data$bio18 <- scale(environmental_data$bio18)

# Calculate the distance matrix from the community matrix
distance_matrix <- vegdist(community_matrix)

# Perform PERMANOVA to test the effect of environmental variables on community composition
permanova_result <- adonis2(distance_matrix ~ bio2 + bio4 + bio15 + bio18 + bio19 + elev + ia, 
                            data = environmental_data, 
                            permutations = 1000)

# Print the PERMANOVA results
print(permanova_result)

# Save the PERMANOVA results to a CSV file for further examination
write.csv(as.data.frame(permanova_result), "permanova_results.csv", row.names = TRUE)
