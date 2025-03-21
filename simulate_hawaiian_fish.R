set.seed(123)
library(tidyverse)  # For data manipulation
library(pavo)       # For color validation

# Wavelength range: 350-800 nm, 4 nm steps (113 bands)
wavelengths <- seq(350, 800, by = 4)

# Species and colors (unchanged)
species_colors <- list(
  "Acanthurus blochii" = c("Blue", "Yellow"),
  "Acanthurus dussumieri" = c("Yellow", "Blue"),
  "Acanthurus nigrofuscus" = c("Brown/Olive"),
  "Acanthurus olivaceus" = c("Orange", "Blue"),
  "Acanthurus thompsoni" = c("Brown/Olive"),
  "Acanthurus triostegus" = c("Black", "White"),
  "Chaetodon auriga" = c("Yellow", "Black", "White"),
  "Chaetodon lunula" = c("Yellow", "Black", "White"),
  "Chaetodon multicinctus" = c("White", "Yellow", "Black"),
  "Chaetodon ornatissimus" = c("Yellow", "White", "Black"),
  "Chaetodon quadrimaculatus" = c("Yellow", "Black", "White"),
  "Chaetodon trifascialis" = c("White", "Yellow", "Black"),
  "Chaetodon unimaculatus" = c("Yellow", "Black", "White"),
  "Chromis agilis" = c("Yellow", "Black"),
  "Chromis hanui" = c("Yellow"),
  "Chromis ovalis" = c("Yellow/UV"),
  "Chromis vanderbilti" = c("Black", "Yellow"),
  "Cirrhitus pinnulatus" = c("Red", "White"),
  "Coris flavovittata" = c("Green", "Blue/red"),
  "Coris gaimard" = c("Blue/red", "Orange-UV", "Green"),
  "Ctenochaetus hawaiiensis" = c("Black", "Orange"),
  "Ctenochaetus strigosus" = c("Yellow"),
  "Dascyllus albisella" = c("Black", "White"),
  "Fistularia commersonii" = c("Green"),
  "Forcipiger flavissimus" = c("Yellow", "Black"),
  "Forcipiger longirostris" = c("Yellow", "Black"),
  "Gomphosus varius" = c("Green", "Blue/red"),
  "Halichoeres ornatissimus" = c("Green", "Orange"),
  "Hemitaurichthys polylepis" = c("Yellow", "Black"),
  "Heniochus diphreutes" = c("White", "Black", "Yellow"),
  "Myripristis kuntee" = c("Red"),
  "Naso lituratus" = c("Blue/UV", "Orange"),
  "Paracirrhites arcatus" = c("Red", "Yellow"),
  "Paracirrhites forsteri" = c("Yellow", "Black"),
  "Parupeneus cyclostomus" = c("Yellow", "Blue"),
  "Parupeneus multifasciatus" = c("Red/UV", "Violet"),
  "Parupeneus porphyreus" = c("Violet"),
  "Pervagor spilosoma" = c("Yellow", "Black"),
  "Plectorhinchus vittatus" = c("Yellow", "Black"),
  "Plectroglyphidodon imparipennis" = c("Brown/Olive"),
  "Plectroglyphidodon johnstonianus" = c("Yellow"),
  "Priacanthus meeki" = c("Red"),
  "Pygoplites diacanthus" = c("Blue", "Yellow", "Black"),
  "Rhinecanthus aculeatus" = c("White", "Black", "Orange", "Blue"),
  "Sargocentron diadema" = c("Red", "White"),
  "Sargocentron xantherythrum" = c("Red", "White"),
  "Scarus dubius" = c("Green", "Blue"),
  "Scarus psittacus" = c("Green"),
  "Sufflamen bursa" = c("Brown/Olive", "Yellow"),
  "Thalassoma duperrey" = c("Green", "Blue/red"),
  "Zebrasoma flavescens" = c("Yellow")
)

# Color definitions with Gaussian peaks
colors <- list(
  "UV" = list(peak = 360, sigma = 40, amp = 50, base = 10),
  "Violet" = list(peak = 400, sigma = 50, amp = 65, base = 15),  # Figure 2B: 400 nm
  "Blue" = list(peak = 470, sigma = 50, amp = 80, base = 15),    # Figure 2B: ~450-505 nm
  "Green" = list(peak = 525, sigma = 60, amp = 70, base = 20),   # Figure 2B: 505 nm
  "Yellow" = list(peak = 575, sigma = 60, amp = 85, base = 15),  # Approx. from Yellow step R50
  "Orange" = list(peak = 620, sigma = 60, amp = 85, base = 15),  # Approx. from Orange step R50
  "Red" = list(peak = 675, sigma = 60, amp = 75, base = 15),     # Approx. from Red step R50
  "Brown/Olive" = list(peak = 600, sigma = 60, amp = 55, base = 20),  # Approx. from Brown step
  "Black" = list(peak = NA, sigma = NA, amp = 15, base = 15),    # Flat
  "White" = list(peak = NA, sigma = NA, amp = 85, base = 85),    # Flat
  "Blue/UV" = list(peaks = c(470, 360), sigmas = c(50, 40), amps = c(80, 50), base = 15),
  "Yellow/UV" = list(peaks = c(575, 360), sigmas = c(60, 40), amps = c(85, 50), base = 15),
  "Orange-UV" = list(peaks = c(620, 360), sigmas = c(60, 40), amps = c(85, 50), base = 15),
  "Red/UV" = list(peaks = c(675, 360), sigmas = c(60, 40), amps = c(75, 50), base = 15),
  "Blue/red" = list(peaks = c(470, 675), sigmas = c(50, 60), amps = c(80, 60), base = 15),
  "Labriform-green" = list(peaks = c(525, 450), sigmas = c(60, 50), amps = c(70, 55), base = 20),
  "Labriform-purple" = list(peaks = c(400, 650), sigmas = c(50, 60), amps = c(65, 60), base = 15)
)

# Curve generation function (Gaussian only)
peak_curve <- function(wl, peak, sigma, amp, base = 0) {
  if (is.na(peak)) rep(amp, length(wl))  # Flat for Black/White
  else base + amp * exp(-((wl - peak)^2) / (2 * sigma^2))
}

# Add noise and light smoothing
add_noise_and_smooth <- function(curve, wl) {
  noisy_curve <- curve + rnorm(length(curve), 0, 2)  # Reduced to 2% for smoothness
  smoothed_curve <- numeric(length(wl))
  for (i in 1:length(wl)) {
    window <- pmax(1, i-1):pmin(length(wl), i+1)  # 3-point moving average (~5 nm)
    smoothed_curve[i] <- mean(noisy_curve[window], na.rm = TRUE)
  }
  pmax(0, pmin(100, smoothed_curve))
}

# Simulate data for one species
simulate_species <- function(species, patches, wl) {
  data_list <- list()
  for (i in 1:3) {  # 3 individuals
    ind_patches <- list()
    for (p in patches) {
      col <- colors[[p]]
      if (is.null(col$peaks)) {  # Single peak or flat
        peak <- col$peak + rnorm(1, 0, 15)
        amp <- col$amp * (1 + rnorm(1, 0, 0.15))
        curve <- peak_curve(wl, peak, col$sigma, amp, col$base)
      } else {  # Complex (dual peaks)
        curve <- peak_curve(wl, col$peaks[1] + rnorm(1, 0, 15), col$sigmas[1], col$amps[1] * (1 + rnorm(1, 0, 0.15)), col$base) +
          peak_curve(wl, col$peaks[2] + rnorm(1, 0, 15), col$sigmas[2], col$amps[2] * (1 + rnorm(1, 0, 0.15)), 0)
      }
      curve <- add_noise_and_smooth(curve, wl)
      ind_patches[[p]] <- curve
    }
    data_list[[paste0("Ind", i)]] <- ind_patches
  }
  data <- expand.grid(Species = species, Individual = paste0("Ind", 1:3), Patch = patches, Wavelength = wl)
  data$Reflectance <- numeric(nrow(data))
  for (i in 1:3) {
    for (p in patches) {
      idx <- data$Species == species & data$Individual == paste0("Ind", i) & data$Patch == p
      data$Reflectance[idx] <- data_list[[i]][[p]]
    }
  }
  data
}

# Generate full dataset
all_data <- do.call(rbind, lapply(names(species_colors), function(s) simulate_species(s, species_colors[[s]], wavelengths)))
write.csv(all_data, "Hawaiian_Reef_Fish_Reflectance.csv", row.names = FALSE)

# Color validation with pavo
validate_colors <- function(data) {
  mean_specs <- data %>%
    group_by(Patch, Wavelength) %>%
    summarise(Reflectance = mean(Reflectance), .groups = "drop") %>%
    pivot_wider(names_from = Wavelength, values_from = Reflectance)
  
  rspec_data <- as.rspec(mean_specs, lim = c(350, 800))
  rgb_colors <- spec2rgb(rspec_data, visual = "cie2", percent = TRUE)
  
  cat("\nSimulated Colors (RGB) vs. Expected:\n")
  for (i in 1:nrow(rgb_colors)) {
    patch <- mean_specs$Patch[i]
    rgb <- rgb_colors[i, ]
    cat(sprintf("%s: RGB(%d, %d, %d) - Expected: %s\n", patch, rgb$r, rgb$g, rgb$b, patch))
  }
  
  plot(rspec_data, main = "Simulated Mean Spectra per Patch (350-800 nm)")
  invisible(rgb_colors)
}

if (requireNamespace("pavo", quietly = TRUE)) {
  validate_colors(all_data)
} else {
  cat("\nNote: Install 'pavo' for color validation.\n")
}