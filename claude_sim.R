# Simulate reflectance data for multiple fish species
# This script generates realistic spectral reflectance data for various fish species
# with different color patches and natural variation

# Load required packages
library(tidyverse)

# Set random seed for reproducibility
set.seed(42)

# Define species data with family and color patches
species_data <- list(
  "Acanthuridae" = list(
    "Acanthurus blochii" = c("Blue", "Yellow"),
    "Acanthurus dussumieri" = c("Blue", "Yellow"),
    "Acanthurus nigrofuscus" = c("Brown/Olive", "Blue"),
    "Acanthurus olivaceus" = c("Blue", "Yellow", "White"),
    "Acanthurus thompsoni" = c("Blue", "Black"),
    "Acanthurus triostegus" = c("White", "Black"),
    "Ctenochaetus strigosus" = c("Brown/Olive", "Blue"),
    "Naso hexacanthus" = c("Blue", "White"),
    "Naso lituratus" = c("Blue", "Yellow", "Black"),
    "Zebrasoma flavescens" = c("Yellow"),
    "Zebrasoma veliferum" = c("Yellow", "Black", "White")
  ),
  "Chaetodontidae" = list(
    "Chaetodon auriga" = c("White", "Yellow", "Black"),
    "Chaetodon lunula" = c("Yellow", "White", "Black"),
    "Chaetodon multicinctus" = c("White", "Black", "Yellow"),
    "Chaetodon ornatissimus" = c("White", "Orange", "Black"),
    "Chaetodon quadrimaculatus" = c("Yellow", "Black", "White"),
    "Chaetodon reticulatus" = c("White", "Black", "Yellow"),
    "Chaetodon trifascialis" = c("White", "Black"),
    "Forcipiger flavissimus" = c("Yellow", "Black")
  ),
  "Labridae" = list(
    "Gomphosus varius" = c("Green", "Blue"),
    "Halichoeres ornatissimus" = c("Green", "Blue", "Orange"),
    "Labroides phthirophagus" = c("Blue", "Black"),
    "Macropharyngodon geoffroy" = c("Green", "Blue", "Orange"),
    "Thalassoma ballieui" = c("Green", "Blue", "Red"),
    "Thalassoma duperrey" = c("Green", "Blue", "Red"),
    "Thalassoma trilobatum" = c("Green", "Blue", "Red")
  ),
  "Balistidae" = list(
    "Melichthys niger" = c("Black", "Blue"),
    "Melichthys vidua" = c("Blue", "White"),
    "Rhinecanthus rectangulus" = c("White", "Black", "Yellow"),
    "Sufflamen bursa" = c("Yellow", "Black", "White")
  ),
  "Mullidae" = list(
    "Mulloidichthys flavolineatus" = c("Yellow", "White"),
    "Mulloidichthys vanicolensis" = c("Red", "Yellow", "White"),
    "Parupeneus multifasciatus" = c("Yellow", "Red", "White", "Violet"),
    "Parupeneus pleurostigma" = c("White", "Red", "Black")
  ),
  "Cirrhitidae" = list(
    "Paracirrhites arcatus" = c("Red", "White"),
    "Paracirrhites forsteri" = c("Red", "White")
  ),
  "Scaridae" = list(
    "Scarus psittacus" = c("Green", "Blue"),
    "Scarus rubroviolaceus" = c("Green", "Blue", "Red")
  ),
  "Pomacentridae" = list(
    "Chromis hanui" = c("Brown/Olive", "Blue"),
    "Chromis ovalis" = c("Brown/Olive", "Blue"),
    "Chromis vanderbilti" = c("Blue", "Yellow"),
    "Chromis verater" = c("Blue", "Brown/Olive"),
    "Abudefduf abdominalis" = c("White", "Black", "Yellow"),
    "Abudefduf sordidus" = c("White", "Black"),
    "Dascyllus albisella" = c("Black", "White"),
    "Plectroglyphidodon johnstonianus" = c("White", "Black", "Yellow")
  ),
  "Zanclidae" = list(
    "Zanclus cornutus" = c("White", "Black", "Yellow")
  ),
  "Lethrinidae" = list(
    "Monotaxis grandoculis" = c("White", "Yellow")
  ),
  "Pomacanthidae" = list(
    "Centropyge potteri" = c("Orange", "Blue")
  )
)

# Extract all species and calculate total number
all_species <- unlist(lapply(species_data, names))
n_species <- length(all_species)
n_ind <- 3  # Number of individuals per species
wavelengths <- seq(400, 800, by = 5)  # Wavelength range in nm
n_wavelengths <- length(wavelengths)

# Define color patches and their spectral characteristics
patches <- list(
  "UV" = list(peak = 350, width = 30, height = 0.8),
  "Violet" = list(peak = 420, width = 40, height = 0.85),
  "Blue" = list(peak = 470, width = 50, height = 0.9),
  "Green" = list(peak = 550, width = 60, height = 0.85),
  "Yellow" = list(peak = 580, width = 45, height = 0.9),
  "Orange" = list(peak = 620, width = 40, height = 0.85),
  "Red" = list(peak = 650, width = 35, height = 0.8),
  "Brown/Olive" = list(peak = 600, width = 70, height = 0.7),
  "Black" = list(peak = 500, width = 200, height = 0.2),
  "White" = list(peak = 500, width = 200, height = 0.9)
)

# Create empty CSV file with headers
write.csv(data.frame(
  Family = character(),
  Species = character(),
  Individual = integer(),
  Patch = character(),
  Wavelength = numeric(),
  Reflectance = numeric(),
  stringsAsFactors = FALSE
), "Simulated_Fish_Reflectance.csv", row.names = FALSE)

# Function to generate reflectance spectrum for a patch
generate_spectrum <- function(wavelengths, peak, width, height, noise = 0.02) {
  # Generate Gaussian curve
  spectrum <- height * exp(-(wavelengths - peak)^2 / (2 * width^2))
  
  # Add noise
  spectrum <- spectrum * (1 + rnorm(length(wavelengths), 0, noise))
  
  # Ensure values are between 0 and 1
  spectrum <- pmax(0, pmin(1, spectrum))
  
  return(spectrum)
}

# Process species in batches for memory efficiency
batch_size <- 5
for (family in names(species_data)) {
  family_species <- species_data[[family]]
  
  for (species in names(family_species)) {
    # Get color patches for this species
    species_patches <- family_species[[species]]
    
    # Generate data for each individual and patch
    for (ind in 1:n_ind) {
      # Add random brightness variation per individual (±10%)
      brightness_factor <- 1 + rnorm(1, 0, 0.1)
      
      for (patch in species_patches) {
        # Get patch parameters
        params <- patches[[patch]]
        
        # Add random wavelength shift (±20 nm)
        shifted_peak <- params$peak + rnorm(1, 0, 20)
        
        # Generate spectrum
        reflectance <- generate_spectrum(
          wavelengths = wavelengths,
          peak = shifted_peak,
          width = params$width,
          height = params$height * brightness_factor
        )
        
        # Add to batch data
        batch_data <- data.frame(
          Family = family,
          Species = species,
          Individual = ind,
          Patch = patch,
          Wavelength = wavelengths,
          Reflectance = reflectance * 100  # Convert to percentage
        )
        
        # Append to CSV file
        write.table(batch_data, "Simulated_Fish_Reflectance.csv", 
                   append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE)
      }
    }
  }
}

# Print summary
cat("Simulation complete. Data saved to 'Simulated_Fish_Reflectance.csv'\n")
cat(sprintf("Generated data for %d species with %d individuals each\n", n_species, n_ind))

# Print data summary
data <- read.csv("Simulated_Fish_Reflectance.csv")
cat("\nData Summary:\n")
cat(sprintf("Number of families: %d\n", length(unique(data$Family))))
cat(sprintf("Number of species: %d\n", length(unique(data$Species))))
cat(sprintf("Number of color patches: %d\n", length(unique(data$Patch))))
cat(sprintf("Wavelength range: %d to %d nm\n", min(data$Wavelength), max(data$Wavelength)))
cat(sprintf("Reflectance range: %.1f to %.1f\n", min(data$Reflectance), max(data$Reflectance)))