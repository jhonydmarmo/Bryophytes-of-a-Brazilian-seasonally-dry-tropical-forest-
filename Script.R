################################################## BOUNDING BOX #########################################################
# Load necessary libraries
library(sf) # Used for spatial data manipulation

# Read the Caatinga shapefile
caatinga <- st_read("C:\\path_to_directory\\caatinga.shp")
bbox <- st_bbox(caatinga)
print(bbox)

################################# FILTERING COORDINATES FROM THE RAW SPREADSHEET #######################################
# Load necessary libraries
library(readxl) # Used for reading Excel files
library(sf) # Used for spatial data manipulation
library(writexl) # Used for saving Excel files

# Read the Excel files with the coordinates
mosses <- read_excel("C:\\path_to_directory\\bryophyta.xlsx")
liverworts <- read_excel("C:\\path_to_directory\\marchantiophyta.xlsx")
hornworts <- read_excel("C:\\path_to_directory\\anthocerotophyta.xlsx")

# Convert the coordinates to an sf object, keeping the CRS of the shapefile
coordinates_mosses <- st_as_sf(mosses, coords = c("decimalLongitude", "decimalLatitude"), crs = st_crs(caatinga))
coordinates_liverworts <- st_as_sf(liverworts, coords = c("decimalLongitude", "decimalLatitude"), crs = st_crs(caatinga))
coordinates_hornworts <- st_as_sf(hornworts, coords = c("decimalLongitude", "decimalLatitude"), crs = st_crs(caatinga))

# Perform the filtering: check which coordinates are within the Caatinga shapefile for each group
indices_inside_mosses <- st_within(coordinates_mosses, caatinga, sparse = FALSE)
indices_inside_liverworts <- st_within(coordinates_liverworts, caatinga, sparse = FALSE)
indices_inside_hornworts <- st_within(coordinates_hornworts, caatinga, sparse = FALSE)

# Filter the rows of the original dataframes based on the indices
coordinates_inside_mosses <- mosses[indices_inside_mosses, ]
coordinates_inside_liverworts <- liverworts[indices_inside_liverworts, ]
coordinates_inside_hornworts <- hornworts[indices_inside_hornworts, ]

# Check the results
print(coordinates_inside_mosses)
print(coordinates_inside_liverworts)
print(coordinates_inside_hornworts)

# Save the filtered coordinates to Excel files
write_xlsx(coordinates_inside_mosses, "filtered_coordinates_mosses.xlsx")
write_xlsx(coordinates_inside_liverworts, "filtered_coordinates_liverworts.xlsx")
write_xlsx(coordinates_inside_hornworts, "filtered_coordinates_hornworts.xlsx")

######################################### CHAO1 RICHNESS ESTIMATION FOR TOTAL SPECIES NUMBER ###########################
# Load necessary libraries
library(vegan) # Used for community ecology analyses, including richness estimators
library(ggplot2) # Used for data visualization

# Read the community matrix from a CSV file
community_matrix <- read.csv("C:\\path_to_directory\\community_matrix.csv", row.names = 1)

# Estimate Chao1 richness for the entire community matrix
est_chao1 <- estaccumR(community_matrix, permutations = 100)

# Visualize the results
options(max.print = 10000)
summary(est_chao1, display = "chao")
str(est_chao1)

# Prepare data for plotting
results <- summary(est_chao1, display = c("S", "chao"))
res_chao <- cbind(results$chao[, 1:4], results$S[, 2:4])
res_chao <- as.data.frame(res_chao)
colnames(res_chao) <- c("Samples", "Chao", "C_lower", "C_upper", "Richness", "R_lower", "R_upper")

# Chao 1 plot
ggplot(res_chao, aes(y = Richness, x = Samples)) +
  geom_point(aes(y = Chao, x = Samples + 0.1), size = 4, 
             color = "darkorange", alpha = 0.7) +
  geom_point(aes(y = Richness, x = Samples), size = 4, 
             color = "cyan4", alpha = 0.7) +
  geom_label(y = 550, x = 255, label = "Estimated richness") +
  geom_label(y = 420, x = 255, label = "Observed richness") + 
  geom_line(aes(y = Chao, x = Samples), color = "darkorange") +
  geom_line(aes(y = Richness, x = Samples), color = "cyan4") +
  geom_linerange(aes(ymin = C_lower, ymax = C_upper, x = Samples + 0.1), color = "darkorange") +
  geom_linerange(aes(ymin = R_lower, ymax = R_upper, x = Samples), color = "cyan4") +
  scale_x_continuous(breaks = seq(0, max(res_chao$Samples, na.rm = TRUE), by = 70)) +
  labs(x = "Number of samples", y = "Estimated richness - Chao 1")


######################################### ESTIMATOR CHAO1 BY LOCALITY ##########################################
# Load necessary libraries
install.packages("vegan")
install.packages("stats")
library(vegan) # Used for community ecology analyses, including richness estimators
library(stats) # Used for statistical tests

# Estimate Chao1 richness for each locality
chao1_estimates <- estimateR(community_matrix, method = "chao")
print(chao1_estimates)
print(str(chao1_estimates))
write.csv(chao1_estimates, file = "chao1_richness.csv", row.names = TRUE)

# Performing the Shapiro-Wilk test for observed richness
shapiro_result_observed <- shapiro.test(chao1_estimates["S.obs", ])
print(shapiro_result_observed)

# Performing the Shapiro-Wilk test for Chao1 estimated richness
shapiro_result_estimated_chao1 <- shapiro.test(chao1_estimates["S.chao1", ])
print(shapiro_result_estimated_chao1)

# Performing the Wilcoxon test for paired samples
wilcoxon_result <- wilcox.test(x = chao1_estimates["S.obs", ], 
                               y = chao1_estimates["S.chao1", ],
                               paired = TRUE)

# Printing the Wilcoxon test results
print(wilcoxon_result)

############################################SAMPBIAS###############################################
Installing necessary packages
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("raster", quietly = TRUE)) install.packages('raster', repos='https://rspatial.r-universe.dev')

# Loading packages
library(devtools)
devtools::install_github("azizka/sampbias") # Sampbias for sampling bias analysis
library(sampbias)
library(maptools)
library(sp)
library(rgdal)
library(rgeos)
library(raster)
library(sf)
library(tidyr)
library(ggplot2)
library(terra)
library(ggspatial)

# Loading data
bryophytes <- read_excel("C:\path_to_directory\bryophyta.xlsx")
caatinga_sf <- st_read("C:\path_to_directory\caatinga.shp")

# Performing sampling bias analysis with specific adjustments for the Caatinga
bryophytes_samp <- calculate_bias(
  x = bryophytes,
  res = 0.0287,
  buffer = 1,
  verbose = TRUE,
  restrict_sample = caatinga_sf
)

# Generating a statistical summary of the sampbias object
summary(bryophytes_samp)

# Creating a plot of the estimated bias weights
plot(bryophytes_samp)

# Projecting the estimated bias effects in space
proj <- project_bias(bryophytes_samp)

# Visualizing the bias effects in space
map_bias(proj, gaz = NULL, sealine = FALSE, type = "sampling_rate")

########################################### GENERALIZED LINEAR MIXED MODELS (GLMM) ######################################

# GLMM WITH ALL VARIABLES SELECTED BY PCA

# Load necessary libraries
library(lme4) # Used for fitting linear mixed-effects models
library(scales) # Used for scaling and rescaling data
library(DHARMa) # Used for diagnostic tools for hierarchical (mixed) regression models
library(car) # Used for companion to applied regression
library(lmtest) # Used for testing linear regression models
library(MuMIn) # Used for model selection and model averaging
library(ggplot2) # Used for data visualization

# Read the richness data from a CSV file
glmm_richness <- read.csv("C:\\path_to_directory\\glmm_richness.csv", row.names = 1)

# Rescale variables correctly
glmm_richness$bio2 <- scale(glmm_richness$bio2)
glmm_richness$bio4 <- scale(glmm_richness$bio4)
glmm_richness$bio15 <- scale(glmm_richness$bio15)
glmm_richness$elev <- scale(glmm_richness$elev)
glmm_richness$ia <- scale(glmm_richness$ia)
glmm_richness$bio18 <- scale(glmm_richness$bio18)
glmm_richness$bio19 <- scale(glmm_richness$bio19)

# Model with just one variable
glmm_model <- glmer.nb(S.obs ~ bio2 + (1 | locality_id), 
                                  data = glmm_richness, 
                                  control = glmerControl(optimizer = "bobyqa"))
summary(glmm_model)

# Adding variables step by step and checking for convergence
glmm_model <- update(glmm_model, . ~ . + bio4)
summary(glmm_model)

glmm_model <- update(glmm_model, . ~ . + bio15)
summary(glmm_model)

glmm_model <- update(glmm_model, . ~ . + elev)
summary(glmm_model)

glmm_model <- update(glmm_model, . ~ . + ia)
summary(glmm_model)

glmm_model <- update(glmm_model, . ~ . + bio18)
summary(glmm_model)

glmm_model <- update(glmm_model, . ~ . + bio19)
summary(glmm_model)

# Perform model selection using dredge
options(na.action = "na.fail")
model_selection <- dredge(glmm_model)

# Display the best models
top_models <- get.models(model_selection, subset = delta < 2)
print(top_models)

# Relative importance of variables
variable_importance <- sw(model_selection)
print(variable_importance)

# Create a dataframe with variable importance
importance_df <- data.frame(
  Variable = c("elev", "ia", "bio4", "bio2", "bio15", "bio18", "bio19"),
  Importance = c(0.93, 0.92, 0.80, 0.77, 0.31, 0.30, 0.30)
)

# Order the dataframe by importance
importance_df <- importance_df[order(importance_df$Importance, decreasing = TRUE), ]

# Create the visualization
ggplot(importance_df, aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Relative Importance of Variables",
       x = "Variables",
       y = "Relative Importance (Sum of Akaike Weights)")

# Calculate marginal and conditional R² using MuMIn
r2_glmm <- r.squaredGLMM(glmm_model)
print(r2_glmm)

# Simulate residuals using DHARMa for diagnostics
simulated_residuals_glmm <- simulateResiduals(fittedModel = glmm_model)

# Plot simulated residuals
plot(simulated_residuals_glmm)

# Test for homoscedasticity of residuals
dispersion_test_glmm <- testDispersion(simulated_residuals_glmm)
print(dispersion_test_glmm)

# Test for normality of standardized residuals
shapiro_test_glmm <- shapiro.test(simulated_residuals_glmm$scaledResiduals)
print(shapiro_test_glmm)

# Test for collinearity among predictor variables
vif_vals_glmm <- vif(glmm_model)
print(vif_vals_glmm)

# Durbin-Watson test for autocorrelation
data_for_test <- data.frame(
  residuals = simulated_residuals_glmm$scaledResiduals,
  fitted_values = fitted(glmm_model)
)
dw_test <- dwtest(residuals ~ fitted_values, data = data_for_test, alternative = "two.sided")
print(dw_test)

# GLM WITH MOST IMPORTANT VARIABLES

# Rescale variables correctly
glmm_richness$bio2 <- scale(glmm_richness$bio2)
glmm_richness$bio4 <- scale(glmm_richness$bio4)
glmm_richness$elev <- scale(glmm_richness$elev)
glmm_richness$ia <- scale(glmm_richness$ia)

# Fit a new model with important variables (importance above 0.6)
# Model with just one variable
simplified_glmm_model <- glmer.nb(S.obs ~ bio2 + (1 | locality_id), 
                                  data = glmm_richness, 
                                  control = glmerControl(optimizer = "bobyqa"))
summary(simplified_glmm_model)
# Add one variable at a time to the simplified model
simplified_glmm_model <- update(simplified_glmm_model, . ~ . + bio4)
summary(simplified_glmm_model)

simplified_glmm_model <- update(simplified_glmm_model, . ~ . + elev)
summary(simplified_glmm_model)

simplified_glmm_model <- update(simplified_glmm_model, . ~ . + ia)
summary(simplified_glmm_model)

# Calculate marginal and conditional R² using MuMIn
r2_glmm_important <- r.squaredGLMM(simplified_glmm_model)
print(r2_glmm_important)

# Simulate residuals using DHARMa for diagnostics
simulated_residuals_important <- simulateResiduals(fittedModel = simplified_glmm_model)

# Plot simulated residuals
plot(simulated_residuals_important, main = "Residuals Important Model")

# Test for homoscedasticity of residuals
dispersion_test_important <- testDispersion(simulated_residuals_important)
print(dispersion_test_important)

# Test for normality of standardized residuals
shapiro_test_important <- shapiro.test(simulated_residuals_important$scaledResiduals)
print(shapiro_test_important)

# Test for collinearity among predictor variables
vif_vals_important <- vif(simplified_glmm_model)
print(vif_vals_important)

# Durbin-Watson test for autocorrelation
data_for_test_important <- data.frame(
  residuals = simulated_residuals_important$scaledResiduals,
  fitted_values = fitted(simplified_glmm_model)
)
dw_test_important <- dwtest(residuals ~ fitted_values, data = data_for_test_important, alternative = "two.sided")
print(dw_test_important)
########################################### PERMUTATIONAL MULTIVARIATE ANALYSIS OF VARIANCE (PERMANOVA) ##################

# Load necessary libraries
library(vegan) # Used for ecological data analysis
library(scales) # Used for scaling and rescaling data

# Scale environmental data variables
environmental_data$bio2 <- scale(environmental_data$bio2)
environmental_data$bio4 <- scale(environmental_data$bio4)
environmental_data$bio15 <- scale(environmental_data$bio15)
environmental_data$elev <- scale(environmental_data$elev)
environmental_data$ia <- scale(environmental_data$ia)
environmental_data$bio18 <- scale(environmental_data$bio18)
environmental_data$bio19 <- scale(environmental_data$bio19)

# Calculate the distance matrix
distance_matrix <- vegdist(community_matrix)

# Perform PERMANOVA
permanova_result <- adonis2(distance_matrix ~ bio2 + bio4 + bio15 + bio18 + bio19 + elev + ia, data = environmental_data, permutations = 1000)
print(permanova_result)

# Save the PERMANOVA results to a CSV file
write.csv(permanova_result, "permanova_results.csv")