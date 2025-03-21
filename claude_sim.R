# Fish Reflectance Spectra Simulation
# This script generates simulated reflectance spectra for 51 fish species
# with 3 individuals per species and varying color patches

# Set seed for reproducibility
set.seed(42)

# Define the parameters
num_species <- 51
individuals_per_species <- 3
wavelength_range <- seq(400, 800, by = 4) # 400-800nm in 4nm steps

# Create a list of species names
species_names <- c(
  "Acanthurus blochii", "Acanthurus dussumieri", "Acanthurus nigrofuscus", 
  "Acanthurus olivaceus", "Acanthurus thompsoni", "Acanthurus triostegus",
  "Chaetodon auriga", "Chaetodon lunula", "Chaetodon multicinctus",
  "Chaetodon ornatissimus", "Chaetodon quadrimaculatus", "Chaetodon reticulatus",
  "Chaetodon trifascialis", "Ctenochaetus strigosus", "Forcipiger flavissimus",
  "Gomphosus varius", "Halichoeres ornatissimus", "Labroides phthirophagus",
  "Macropharyngodon geoffroy", "Melichthys niger", "Melichthys vidua",
  "Monotaxis grandoculis", "Mulloidichthys flavolineatus", "Mulloidichthys vanicolensis",
  "Naso hexacanthus", "Naso lituratus", "Paracirrhites arcatus",
  "Paracirrhites forsteri", "Parupeneus multifasciatus", "Parupeneus pleurostigma",
  "Pseudocheilinus tetrataenia", "Rhinecanthus rectangulus", "Scarus psittacus",
  "Scarus rubroviolaceus", "Stethojulis balteata", "Sufflamen bursa",
  "Thalassoma ballieui", "Thalassoma duperrey", "Thalassoma trilobatum",
  "Zebrasoma flavescens", "Zebrasoma veliferum", "Chromis hanui",
  "Chromis ovalis", "Chromis vanderbilti", "Chromis verater",
  "Abudefduf abdominalis", "Abudefduf sordidus", "Centropyge potteri",
  "Dascyllus albisella", "Plectroglyphidodon johnstonianus", "Zanclus cornutus"
)

# Define common color patches for reef fish
color_patches <- c("Blue", "Yellow", "Green", "Red", "Black", "White", "Orange", "Purple", "Brown", "Pink")

# Function to create base reflectance curves for different colors
create_base_curve <- function(color, wavelength) {
  # Different reflectance patterns based on color
  # Each color has its characteristic spectral signature
  
  reflectance <- numeric(length(wavelength))
  
  if (color == "Blue") {
    # Blue: high reflectance in 450-490nm, lower elsewhere
    reflectance <- 10 + 90 * exp(-0.5 * ((wavelength - 470) / 50)^2)
  } 
  else if (color == "Yellow") {
    # Yellow: high reflectance above 520nm, low below
    reflectance <- 10 + 90 / (1 + exp(-0.1 * (wavelength - 520)))
  } 
  else if (color == "Green") {
    # Green: peak around 520-550nm
    reflectance <- 10 + 90 * exp(-0.5 * ((wavelength - 530) / 60)^2)
  } 
  else if (color == "Red") {
    # Red: high reflectance above 600nm
    reflectance <- 10 + 90 / (1 + exp(-0.1 * (wavelength - 600)))
  } 
  else if (color == "Black") {
    # Black: low reflectance across spectrum
    reflectance <- 5 + 10 * sin(wavelength/100)
  } 
  else if (color == "White") {
    # White: high reflectance across spectrum
    reflectance <- 90 + 10 * sin(wavelength/200)
  } 
  else if (color == "Orange") {
    # Orange: between red and yellow
    reflectance <- 10 + 90 / (1 + exp(-0.1 * (wavelength - 580)))
  } 
  else if (color == "Purple") {
    # Purple: peaks at both blue and red ends
    blue_component <- 40 * exp(-0.5 * ((wavelength - 450) / 40)^2)
    red_component <- 30 / (1 + exp(-0.1 * (wavelength - 650)))
    reflectance <- 10 + blue_component + red_component
  } 
  else if (color == "Brown") {
    # Brown: gradual increase above 550nm, but not as high as red/yellow
    reflectance <- 20 + 40 / (1 + exp(-0.05 * (wavelength - 600)))
  } 
  else if (color == "Pink") {
    # Pink: like red but with some blue component
    blue_component <- 30 * exp(-0.5 * ((wavelength - 450) / 60)^2)
    red_component <- 60 / (1 + exp(-0.1 * (wavelength - 600)))
    reflectance <- 20 + blue_component + red_component
  } 
  else {
    # Default curve if color not recognized
    reflectance <- 50 * rep(1, length(wavelength))
  }
  
  # Ensure all values are within 0-100 range
  reflectance <- pmax(0, pmin(100, reflectance))
  
  return(reflectance)
}

# Function to add species-specific variations to base curves
modify_curve_for_species <- function(base_curve, species_index) {
  # Add species-specific modifications to make curves unique per species
  # Using the species index to create consistent but varied patterns
  
  # Scale factor varies by species (0.8 to 1.2)
  scale_factor <- 0.8 + (species_index / num_species) * 0.4
  
  # Shift factor varies by species (-50 to +50 nm)
  shift <- -50 + (species_index / num_species) * (100)
  
  # Apply modifications
  n <- length(base_curve)
  shifted_indices <- round(seq(1, n, length.out = n) + (shift/4))
  shifted_indices <- pmin(pmax(1, shifted_indices), n)
  
  # Create shifted and scaled curve
  modified_curve <- base_curve[shifted_indices] * scale_factor
  
  # Ensure all values are within 0-100 range
  modified_curve <- pmax(0, pmin(100, modified_curve))
  
  return(modified_curve)
}

# Function to add individual variation within species
add_individual_variation <- function(species_curve, individual_index) {
  # Add noise to simulate individual variation within a species
  
  # Noise amplitude varies by individual (1-5%)
  noise_amplitude <- 1 + individual_index * 2
  
  # Generate noise
  noise <- rnorm(length(species_curve), mean = 0, sd = noise_amplitude)
  
  # Add noise to curve
  varied_curve <- species_curve + noise
  
  # Ensure all values are within 0-100 range
  varied_curve <- pmax(0, pmin(100, varied_curve))
  
  return(varied_curve)
}

# Assign colors to species (2-4 colors per species)
set.seed(42) # Reset seed for reproducibility
species_colors <- list()
for (i in 1:num_species) {
  # Each species has 2-4 color patches
  num_colors <- sample(2:4, 1)
  species_colors[[i]] <- sample(color_patches, num_colors)
}

# Create data frame to store all the data
reflectance_data <- data.frame()

# Generate data for each species
for (s in 1:num_species) {
  species_name <- species_names[s]
  colors <- species_colors[[s]]
  
  for (c in colors) {
    # Get base reflectance curve for this color
    base_curve <- create_base_curve(c, wavelength_range)
    
    # Modify curve for this specific species
    species_curve <- modify_curve_for_species(base_curve, s)
    
    # Generate data for each individual
    for (i in 1:individuals_per_species) {
      # Add individual variation
      individual_curve <- add_individual_variation(species_curve, i)
      
      # Add to data frame
      individual_name <- paste0("Ind", i)
      
      for (w in 1:length(wavelength_range)) {
        reflectance_data <- rbind(reflectance_data, data.frame(
          Species = species_name,
          Individual = individual_name,
          Patch = c,
          Wavelength = wavelength_range[w],
          Reflectance = individual_curve[w]
        ))
      }
    }
  }
  
  # Print progress
  if (s %% 10 == 0) {
    cat("Processed", s, "of", num_species, "species\n")
  }
}

# Write to CSV
write.csv(reflectance_data, "Simulated_Fish_Reflectance.csv", row.names = FALSE)

cat("Simulation complete. Data saved to 'Simulated_Fish_Reflectance.csv'\n")
cat("Generated data for", num_species, "species with", individuals_per_species, "individuals each\n")
cat("Total records:", nrow(reflectance_data), "\n")

# Print a summary of the data structure
cat("\nData Summary:\n")
cat("Number of species:", length(unique(reflectance_data$Species)), "\n")
cat("Number of color patches:", length(unique(reflectance_data$Patch)), "\n")
cat("Wavelength range:", min(reflectance_data$Wavelength), "to", max(reflectance_data$Wavelength), "nm\n")
cat("Reflectance range:", min(reflectance_data$Reflectance), "to", max(reflectance_data$Reflectance), "\n")