set.seed(123)
library(tidyverse)  # For data manipulation

wavelengths <- seq(400, 800, by = 4)  # 101 bands, ULTRIS-like

# Species and colors from Table 1 (unchanged)
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

# Color definitions (unchanged from last revision)
colors <- list(
  "UV" = list(type = "peak", peak = 360, sigma = 30, amp = 50),
  "Violet" = list(type = "peak", peak = 425, sigma = 35, amp = 65),
  "Blue" = list(type = "peak", peak = 470, sigma = 35, amp = 80),
  "Green" = list(type = "peak", peak = 525, sigma = 40, amp = 70),
  "Yellow" = list(type = "step", R50 = 515, k = 0.15, amp = 85),
  "Orange" = list(type = "step", R50 = 575, k = 0.15, amp = 85),
  "Red" = list(type = "step", R50 = 650, k = 0.15, amp = 75),
  "Brown/Olive" = list(type = "step", R50 = 550, k = 0.1, amp = 55),
  "Black" = list(type = "flat", amp = 15),
  "White" = list(type = "flat", amp = 85),
  "Blue/UV" = list(type = "complex", peaks = c(470, 360), sigmas = c(35, 30), amps = c(80, 50)),
  "Yellow/UV" = list(type = "complex", R50 = 515, k = 0.15, amp = 85, uv_peak = 360, uv_amp = 50),
  "Orange-UV" = list(type = "complex", R50 = 575, k = 0.15, amp = 85, uv_peak = 360, uv_amp = 50),
  "Red/UV" = list(type = "complex", R50 = 650, k = 0.15, amp = 75, uv_peak = 360, uv_amp = 50),
  "Blue/red" = list(type = "complex", peaks = c(470, 650), sigmas = c(35, 40), amps = c(80, 60)),
  "Labriform-green" = list(type = "complex", peaks = c(525, 450), sigmas = c(40, 35), amps = c(70, 55)),
  "Labriform-purple" = list(type = "complex", peaks = c(425, 600), sigmas = c(35, 40), amps = c(65, 60))
)

# Curve generation functions
peak_curve <- function(wl, peak, sigma, amp) amp * exp(-((wl - peak)^2) / (2 * sigma^2))
step_curve <- function(wl, R50, k, amp) amp / (1 + exp(-k * (wl - R50)))
flat_curve <- function(wl, amp) rep(amp, length(wl))

# Reduced smoothing (5 nm FWHM)
smooth_curve <- function(curve, wl, fwhm = 5) {
  sigma <- fwhm / 2.355
  weights <- dnorm(wl, mean = wl, sd = sigma)
  conv <- convolve(curve, rev(weights), type = "open")
  conv[(length(wl)+1):(2*length(wl))]
}

# Add noise function
add_noise <- function(curve) curve + rnorm(length(curve), 0, 3)  # SD = 3% noise

# Simulate data for one species (no fading)
simulate_species <- function(species, patches, wl) {
  data_list <- list()
  for (i in 1:3) {  # 3 individuals, no fading
    ind_patches <- list()
    for (p in patches) {
      col <- colors[[p]]
      if (col$type == "peak") {
        peak <- col$peak + rnorm(1, 0, 15)
        amp <- col$amp * (1 + rnorm(1, 0, 0.15))  # No fade multiplier
        curve <- peak_curve(wl, peak, col$sigma, amp)
      } else if (col$type == "step") {
        R50 <- col$R50 + rnorm(1, 0, 15)
        amp <- col$amp * (1 + rnorm(1, 0, 0.15))
        curve <- step_curve(wl, R50, col$k, amp)
      } else if (col$type == "flat") {
        amp <- col$amp * (1 + rnorm(1, 0, 0.15))
        curve <- flat_curve(wl, amp)
      } else {  # complex
        if (!is.null(col$peaks)) {
          curve <- peak_curve(wl, col$peaks[1] + rnorm(1, 0, 15), col$sigmas[1], col$amps[1] * (1 + rnorm(1, 0, 0.15))) +
            peak_curve(wl, col$peaks[2] + rnorm(1, 0, 15), col$sigmas[2], col$amps[2] * (1 + rnorm(1, 0, 0.15)))
        } else {
          curve <- step_curve(wl, col$R50 + rnorm(1, 0, 15), col$k, col$amp * (1 + rnorm(1, 0, 0.15))) +
            peak_curve(wl, col$uv_peak, 30, col$uv_amp * (1 + rnorm(1, 0, 0.15)))
        }
      }
      curve <- pmax(0, pmin(100, smooth_curve(curve, wl)))  # Apply reduced smoothing
      curve <- add_noise(curve)  # Add noise
      curve <- pmax(0, pmin(100, curve))  # Re-clamp after noise
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

# Save the data
write.csv(all_data, "Hawaiian_Reef_Fish_Reflectance.csv", row.names = FALSE)

# Color Validation Notes:
# Several approaches were attempted to validate colors using the pavo package:
#
# 1. Direct spec2rgb conversion:
#    - Tried using pavo::spec2rgb(rspec_data) directly
#    - Failed due to wavelength range mismatch with built-in visual systems
#
# 2. CIE2 color space:
#    - Attempted using visual = "cie2" in vismodel
#    - Required wavelength range 360-830nm
#    - Data only covers 400-800nm, causing mismatch errors
#
# 3. Standard human vision:
#    - Tried visual = "segment" in vismodel
#    - Still encountered wavelength range issues
#
# 4. Attempted solutions:
#    - Interpolation to different wavelength ranges (380-700nm, 360-830nm)
#    - Different step sizes (1nm vs 5nm)
#    - Various data formats (matrix vs data frame)
#    - None resolved the fundamental wavelength mismatch
#
# Future improvements could include:
# 1. Extending wavelength range in simulation to match CIE2 (360-830nm)
# 2. Using alternative color space conversion methods
# 3. Implementing custom color conversion using specific cone sensitivities
# 4. Using colorspace package directly instead of pavo for RGB conversion

# # Optional color validation (requires pavo package) - commented out for now
# if (requireNamespace("pavo", quietly = TRUE)) {
#   validate_colors <- function(data, colors) {
#     # Calculate mean spectra per patch
#     mean_specs <- data %>%
#       group_by(Patch, Wavelength) %>%
#       summarise(Reflectance = mean(Reflectance), .groups = "drop")
#     
#     # Print wavelength range
#     cat("\nWavelength range in data:", 
#         sprintf("%.1f to %.1f nm\n", 
#                 min(mean_specs$Wavelength), 
#                 max(mean_specs$Wavelength)))
#     
#     # Filter to standard visible range (400-700nm)
#     mean_specs <- mean_specs %>%
#       filter(Wavelength >= 400 & Wavelength <= 700) %>%
#       # Interpolate to 5nm steps
#       group_by(Patch) %>%
#       complete(Wavelength = seq(400, 700, by = 5)) %>%
#       mutate(Reflectance = approx(Wavelength[!is.na(Reflectance)], 
#                                  Reflectance[!is.na(Reflectance)], 
#                                  xout = Wavelength)$y) %>%
#       ungroup()
#     
#     # Create data frame in pavo format
#     pavo_data <- data.frame(
#       wl = sort(unique(mean_specs$Wavelength))
#     )
#     
#     # Add each patch as a column
#     for(patch in unique(mean_specs$Patch)) {
#       patch_data <- mean_specs %>%
#         filter(Patch == patch) %>%
#         arrange(Wavelength) %>%
#         pull(Reflectance)
#       pavo_data[[patch]] <- patch_data
#     }
#     
#     # Convert to rspec object and calculate RGB values
#     rspec_data <- pavo::as.rspec(pavo_data)
#     rgb_colors <- pavo::spec2rgb(rspec_data)
#     
#     # Print results
#     cat("\nSimulated Colors (RGB) vs. Expected:\n")
#     patch_names <- unique(data$Patch)
#     for (i in 1:length(patch_names)) {
#       patch <- patch_names[i]
#       rgb <- rgb_colors[i,]
#       cat(sprintf("%s: RGB(%d, %d, %d) - Expected: %s\n", 
#                   patch, rgb[1]*255, rgb[2]*255, rgb[3]*255, patch))
#     }
#     
#     # Plot spectra with custom colors
#     par(mfrow = c(2, 1))
#     plot(rspec_data, main = "Simulated Mean Spectra per Patch")
#     
#     # Plot RGB representation using the calculated colors
#     barplot(t(rgb_colors), beside = TRUE, 
#             names.arg = patch_names,
#             main = "RGB Values per Patch",
#             col = c("red", "green", "blue"),
#             legend.text = c("R", "G", "B"),
#             args.legend = list(x = "topright"))
#     par(mfrow = c(1, 1))
#     
#     # Return the RGB colors for potential further use
#     return(rgb_colors)
#   }
#   
#   # Run validation if pavo is available
#   rgb_values <- validate_colors(all_data, colors)
# } else {
#   cat("\nNote: Color validation skipped (pavo package not available)\n")
#   cat("To enable color validation, install pavo package and XQuartz\n")
# }