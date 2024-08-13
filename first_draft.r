library(dplyr)
library(duckdb)
library(ggplot2)
library(ebirdst)
library(devtools)
library(ratlas)
# library(ggeffects)
# library(lme4)
library(httr) # To download the data
library(sf) # To work with spatial data


#### AVINET TRAITS ####
# Sptraits
birdlife <- read.csv("/home/claire/BDQC-GEOBON/data/PCI/ELEData/ELEData/TraitData/AVONET1_BirdLife.csv")
# colnames(birdlife)
spTraits <- birdlife[, c("Species1", "Family1", "Order1", "Habitat", "Migration", "Trophic.Level", "Trophic.Niche", "Primary.Lifestyle")]
# unique(spTraits$Species1)

#### SAPTIAL SCALE ####
# Define the path to the shapefile
shp_path <- "/home/claire/BDQC-GEOBON/data/QUEBEC_regions/CERQ_SHP/CR_NIV_02_S.shp"
# Read the shapefile
map_shp <- st_read(shp_path)
# Extract the centroids of the geometries
centroids <- st_centroid(map_shp)
# Extract the latitude values from the centroids
latitudes <- st_coordinates(centroids)[, 2]
# Create a new variable with the integer part of the latitude values
map_shp_to_plot <- map_shp %>%
    mutate(latitude_int = floor(latitudes))
map_shp_to_plot <- map_shp_to_plot %>%
    mutate(latitude_cat = cut(latitude_int, breaks = 10, labels = FALSE))
map_shp <- map_shp %>%
    select(geometry, FID02)
# Plot the shapefile
ggplot() +
    geom_sf(data = map_shp)


#### EXTRACT ATLAS OCCURRENCES ####
atlas_local <- function(parquet_file,
                        tblname = "atlas") {
    requireNamespace("duckdbfs")
    atlas <- duckdbfs::open_dataset(parquet_file, tblname = tblname)
    atlas
}
atlas <- atlas_local("/home/claire/BDQC-GEOBON/data/atlas/atlas_2024-07-16.parquet") |> filter(within_quebec == "t")
spinatlas <- atlas |>
    filter(within_quebec == "t") |>
    group_by(valid_scientific_name) |>
    summarize(cnt = n()) |>
    filter(cnt > 5000) |>
    select(valid_scientific_name) |>
    collect()
mb <- atlas |>
    group_by(valid_scientific_name) |>
    summarize(cnt = count()) |>
    collect()
res <- ebirdst_runs |> filter(is_resident == FALSE)
res <- res$scientific_name
species <- res[res %in% spinatlas$valid_scientific_name]
selected_spTraits <- spTraits[spTraits$Species1 %in% species, ] # Data filterd
julian_dates <- read.csv("spqt.csv")
julian_dates <- julian_dates |>
    filter(jquant >= 31 & jquant <= 181)
julian_dates <- left_join(julian_dates, selected_spTraits, by = c("species" = "Species1"))
julian_dates_mean_trophic_niche <- julian_dates |>
    group_by(FID02, Trophic.Niche, year_obs) |>
    summarize(jj_date = mean(jquant))
julian_dates_mean_trophic_niche_sf <- left_join(map_shp_to_plot, julian_dates_mean_trophic_niche, by = "FID02")
julian_dates_sf <- left_join(map_shp_to_plot, julian_dates, by = "FID02")
# Group the data by Trophic.Level and fit linear regression for each group
lm_by_trophic <- julian_dates %>%
    group_by(Trophic.Niche) %>%
    do(model = lm(jquant ~ year_obs, data = .))
# Extract coefficients for each model
coefs <- lm_by_trophic %>%
    summarize(intercept = coef(model)[1], slope = coef(model)[2])
# Print coefficients
print(coefs)
# Plot the linear regression lines for each Trophic.Level
julian_dates <- julian_dates %>%
    filter(!is.na(Trophic.Niche))
ggplot(data = julian_dates, aes(x = year_obs, y = jquant, color = Trophic.Niche)) +
    geom_smooth(method = "lm", se = TRUE) +
    labs(title = "Linear Regression by Trophic Level", x = "Year", y = "Julian Date") +
    theme_bw()
### Filterd data to avoid lm in empty FID
filtered_julian_dates <- julian_dates %>%
    group_by(FID02) %>%
    filter(n_distinct(year_obs) >= 3) %>%
    ungroup()
# Model
lm_by_trophic_fid <- filtered_julian_dates %>%
    group_by(Trophic.Niche, FID02) %>%
    do(model = lm(jquant ~ year_obs, data = .))
# Predict data
# Create a sequence of years from 1990 to 2023
years <- 1990:2023
# Create a data frame to store the predictions
predictions <- data.frame()
# Iterate over each model in lm_by_trophic_fid and predict values for each year
for (i in 1:nrow(lm_by_trophic_fid)) {
    model <- lm_by_trophic_fid$model[[i]]
    trophic_niche <- lm_by_trophic_fid$Trophic.Niche[i]
    fid02 <- lm_by_trophic_fid$FID02[i]
    # Create a data frame with the years for prediction
    predict_data <- data.frame(year_obs = years)
    # Predict the values
    predict_data$predicted_jquant <- predict(model, newdata = predict_data)
    # Add Trophic.Niche and FID02 for identification
    predict_data$Trophic.Niche <- trophic_niche
    predict_data$FID02 <- fid02
    # Combine predictions with existing data
    predictions <- rbind(predictions, predict_data)
}
predictions_arrival <- left_join(map_shp_to_plot, predictions, by = "FID02")
