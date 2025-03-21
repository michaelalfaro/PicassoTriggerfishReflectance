set.seed(123)
wavelengths <- seq(400, 800, by = 4)  # 101 bands, ULTRIS-like

# Full list of 51 species from Table 1 with colors
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

# Color definitions
colors <- list(
  "UV" = list(type = "peak", peak = 360, sigma = 30, amp = 40),
  "Violet" = list(type = "peak", peak = 425, sigma = 50, amp = 60),
  "Blue" = list(type = "peak", peak = 470, sigma = 50, amp = 85),
  "Green" = list(type = "peak", peak = 525, sigma = 50, amp = 70),
  "Yellow" = list(type = "step", R50 = 515, k = 0.1, amp = 90),
  "Orange" = list(type = "step", R50 = 575, k = 0.1, amp = 90),
  "Red" = list(type = "step", R50 = 650, k = 0.1, amp = 80),
  "Brown/Olive" = list(type = "step", R50 = 550, k = 0.08, amp = 50),
  "Black" = list(type = "flat", amp = 10),
  "White" = list(type = "flat", amp = 80),
  "Blue/UV" = list(type = "complex", peaks = c(470, 360), sigmas = c(50, 30), amps = c(85, 40)),
  "Yellow/UV" = list(type = "complex", R50 = 515, k = 0.1, amp = 90, uv_peak = 360, uv_amp = 40),
  "Orange-UV" = list(type = "complex", R50 = 575, k = 0.1, amp = 90, uv_peak = 360, uv_amp = 40),
  "Red/UV" = list(type = "complex", R50 = 650, k = 0.1, amp = 80, uv_peak = 360, uv_amp = 40),
  "Blue/red" = list(type = "complex", peaks = c(470, 650), sigmas = c(50, 50), amps = c(85, 60)),
  "Labriform-green" = list(type = "complex", peaks = c(525, 450), sigmas = c(50, 40), amps = c(70, 50)),
  "Labriform-purple" = list(type = "complex", peaks = c(425, 600), sigmas = c(40, 50), amps = c(60, 50))
)

# Curve generation functions
peak_curve <- function(wl, peak, sigma, amp) amp * exp(-((wl - peak)^2) / (2 * sigma^2))
step_curve <- function(wl, R50, k, amp) amp / (1 + exp(-k * (wl - R50)))
flat_curve <- function(wl, amp) rep(amp, length(wl))

# Smoothing function (10 nm FWHM)
smooth_curve <- function(curve, wl, fwhm = 10) {
  sigma <- fwhm / 2.355
  weights <- dnorm(wl, mean = wl, sd = sigma)
  conv <- convolve(curve, rev(weights), type = "open")
  conv[(length(wl)+1):(2*length(wl))]
}

# Simulate data for one species
simulate_species <- function(species, patches, wl) {
  data_list <- list()
  for (i in 1:3) {  # 3 individuals
    fade <- (i == 3)  # Last individual faded
    ind_patches <- list()
    for (p in patches) {
      col <- colors[[p]]
      if (col$type == "peak") {
        peak <- col$peak + rnorm(1, 0, 20)
        amp <- col$amp * (1 + rnorm(1, 0, 0.1)) * (if (fade) 0.7 else 1)
        curve <- peak_curve(wl, peak, col$sigma, amp)
      } else if (col$type == "step") {
        R50 <- col$R50 + rnorm(1, 0, 20)
        amp <- col$amp * (1 + rnorm(1, 0, 0.1)) * (if (fade) 0.7 else 1)
        curve <- step_curve(wl, R50, col$k, amp)
      } else if (col$type == "flat") {
        amp <- col$amp * (1 + rnorm(1, 0, 0.1)) * (if (fade) 0.7 else 1)
        curve <- flat_curve(wl, amp)
      } else {  # complex
        if (!is.null(col$peaks)) {
          curve <- peak_curve(wl, col$peaks[1] + rnorm(1, 0, 20), col$sigmas[1], col$amps[1] * (1 + rnorm(1, 0, 0.1)) * (if (fade) 0.7 else 1)) +
            peak_curve(wl, col$peaks[2] + rnorm(1, 0, 20), col$sigmas[2], col$amps[2] * (1 + rnorm(1, 0, 0.1)) * (if (fade) 0.7 else 1))
        } else {
          curve <- step_curve(wl, col$R50 + rnorm(1, 0, 20), col$k, col$amp * (1 + rnorm(1, 0, 0.1)) * (if (fade) 0.7 else 1)) +
            peak_curve(wl, col$uv_peak, 30, col$uv_amp * (1 + rnorm(1, 0, 0.1)) * (if (fade) 0.7 else 1))
        }
      }
      curve <- pmax(0, pmin(100, smooth_curve(curve, wl)))
      ind_patches[[p]] <- curve
    }
    data_list[[paste0("Ind", i)]] <- ind_patches
  }
  data <- expand.grid(Species = species, Individual = paste0("Ind", 1:3), 
                      Patch = patches, Wavelength = wl)
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