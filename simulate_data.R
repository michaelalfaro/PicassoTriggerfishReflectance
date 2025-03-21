set.seed(123)  # For reproducibility
wavelengths <- seq(400, 800, by = 5)  # 81 wavelengths
n_ind <- 5
n_rep <- 3
n_patch <- 4
patches <- c("White", "Black", "Orange", "Blue")

# Base reflectance curves (0-100% scale)
white_base <- 80 * exp(-((wavelengths - 500)^2) / (2 * 100^2)) + 70
black_base <- rep(10, length(wavelengths))
orange_base <- 90 * exp(-((wavelengths - 620)^2) / (2 * 80^2))
blue_base <- 85 * exp(-((wavelengths - 470)^2) / (2 * 60^2))

simulate_reflectance <- function(base_curve, ind_var, fade = FALSE) {
  ind_curve <- base_curve * (1 + rnorm(1, 0, 0.1))
  if (fade) ind_curve <- ind_curve * 0.7
  ind_curve <- ind_curve + rnorm(1, 0, 20) * (wavelengths - 500) / 1000
  ind_curve <- pmax(0, pmin(100, ind_curve))
  reps <- replicate(n_rep, ind_curve + rnorm(length(wavelengths), 0, 2))
  return(reps)
}

data_list <- list()
fade_ind <- sample(1:n_ind, 1)
for (i in 1:n_ind) {
  fade <- (i == fade_ind)
  data_list[[i]] <- list(
    White = simulate_reflectance(white_base, i, fade),
    Black = simulate_reflectance(black_base, i, fade),
    Orange = simulate_reflectance(orange_base, i, fade),
    Blue = simulate_reflectance(blue_base, i, fade)
  )
}

data <- expand.grid(
  Individual = paste0("Ind", 1:n_ind),
  Patch = patches,
  Replicate = 1:n_rep,
  Wavelength = wavelengths
)
data$Reflectance <- numeric(nrow(data))
for (i in 1:n_ind) {
  for (p in patches) {
    for (r in 1:n_rep) {
      idx <- data$Individual == paste0("Ind", i) & data$Patch == p & data$Replicate == r
      data$Reflectance[idx] <- data_list[[i]][[p]][, r]
    }
  }
}

write.csv(data, "Picasso_Triggerfish_Reflectance.csv", row.names = FALSE)