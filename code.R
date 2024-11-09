#                               FICK's PRINCIPLE 
################################################################################
# References: 
# [1]
# Fick, Adolf. 1870. Ueber die Messung des Blutquantums in den Herzventrikeln. 
# Sitzungsberichte der Physikalisch-Medizinischen Gesellschaft zu Würzburg: 16 

# Function to calculate oxygen consumption (VO2) using the Fick principle
calculate_VO2 <- function(Q = 5000, CaO2 = 20, CvO2 = 15) {
  # Convert oxygen content from per 100 mL blood to per mL blood
  CaO2_mL_per_mL <- CaO2 / 100
  CvO2_mL_per_mL <- CvO2 / 100
  
  # Calculate VO2 in mL O2 per minute & Return the result
  VO2 <- Q * (CaO2_mL_per_mL - CvO2_mL_per_mL)
  return(VO2) 
}

# Call the function 
VO2 <- calculate_VO2(Q = 5000, CaO2 = 20, CvO2 = 10) 

# Display the result
cat("Oxygen consumption (VO2) is", VO2, "mL O2 per minute\n")


                            # Severinghaus Equations 
################################################################################
# References:
# [1]
# Severinghaus JW. Simple, accurate equations for human blood O2 dissociation computations. 
# J Appl Physiol 46: 599–602, 1979. doi: 10.1152/jappl.1979.46.3.599.

# [2]
# Ellis RK. Determination of PO2 from saturation. J Appl Physiol 67: 902–902, 1989. 
# doi: 10.1152/jappl.1989.67.2.902.

# Calculating fractional oxyhemoglobin saturation (S) from partial pressure of O2 (pO2 in mmHg)
calc_O2_sat <- function(pO2 = 100) {
  saturation <- (((pO2 ^ 3 + 150 * pO2) ^ -1 * 23400) + 1) ^ -1
  return(saturation)
}

calc_O2_sat()

# Above equation can be inverted to calculate pO2 from fractional O2 saturation 
calc_pO2 <- function(S = 0.36) {
  A <- 11700 * (S^-1 - 1)^-1
  B <- (50^3 + A^2)^0.5
  PO2st <- (B + A)^(1/3) - (B - A)^(1/3)
  return(PO2st)
}

calc_pO2() 

# Calculating O2 content based on pO2, and Hb where the saturation equation ist nested
calc_O2_cont <- function(pO2 = 100, Hb = 15) {
  # Calculate O2 content in mL O2 per dL blood
  O2_content <- (1.34 * Hb * ((((pO2 ^ 3 + 150 * pO2) ^ -1 * 23400) + 1) ^ -1)) + (0.0031 * pO2)
  return(O2_content)
}

calc_O2_cont()



#             Severinghaus Implementation by Mairbäurl (unfinished)
################################################################################
# [1]
# Severinghaus, J.W., 1979. Simple, accurate equations for human blood O2 
# dissociation computations. J. Appl. Physiol. 46, 599–602.
# https://doi.org/10.1152/jappl.1979.46.3.599

# [2]
# Ellis, R.K., 1989. Determination of PO2 from saturation. 
# J. Appl. Physiol. 67, 902–902. https://doi.org/10.1152/jappl.1989.67.2.902

# [3]
# Okada, Y., Tyuma, I., Sugimoto, T., 1977. Evaluation of Severinghaus’ Equation 
# and Its Modification for 2, 3-Dpg. Jpn. J. Physiol. 27, 135–144. 
# https://doi.org/10.2170/jjphysiol.27.135

# [4]
# Mairbäurl, H., Weber, R.E., 2012. Oxygen transport by hemoglobin. 
# Compr. Physiol. 2, 1463–1489. https://doi.org/10.1002/cphy.c080113

# Required libraries to run the code 
library(ggplot2)

S <- seq(0, 1, 0.01)

# Reference [1]
Severinghaus_S <- function(x){ (((x^3 + 150*x)^-1 * 23400) + 1)^-1 }

# Note that S must be calculated as fractional saturation not percentage.
# Reference [2] used this modification as it avoids problems with log function in R 
Severinghaus_PO2st <- function(S) {
  A <- 11700 * (S^-1 - 1)^-1
  B <- (50^3 + A^2)^0.5
  PO2st <- (B + A)^(1/3) - (B - A)^(1/3)
  return(PO2st)
}

PO2st <- Severinghaus_PO2st(S)

# Conditions that can be changed
P50  <- 27.1 # Okkada P50 of standard blood = 27.1 but Severingahus = 26.86
pH   <- 7.4  # Standard 7.4
Temp <- 37.0 # Standard 37
BE   <- 0    # Standard 0 ?
DPG  <- 1.02 # accoreding to table 1 Okada 1.02 is norm. human blood with P50=27.1

# Reference [3]
logP50 <- log(P50) + 0.48*(7.4-pH) + 0.024*(Temp-37) + 0.0013*BE + (0.135*DPG)-0.116

# Reference [4] revert log with exp() and store observed P50 value as P50obs
P50obs <- exp(logP50)
PO2act <- PO2st * P50obs / 26.86

curves <- data.frame(S, PO2act, PO2st)[complete.cases(S, PO2act, PO2st), ]

# Plotting both PO2act and PO2st on the x-axis and S on the y-axis
ggplot(curves, aes(x = PO2act, y = S)) +
  geom_line(aes(color = "PO2act")) +
  geom_line(aes(x = PO2st, color = "PO2st")) +
  ylab("Fractional Saturation (S)") +
  xlab("Partial Pressure of Oxygen (PO2)") +
  labs(color = "PO2 Type") +
  scale_x_continuous(
    name = expression("PO"[2] * " (mmHg)"),
    breaks = seq(0, 100, 10),
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  geom_vline(xintercept = 21, linetype = "dashed", color = "black") +
  theme_bw() +
  theme(text = element_text(size = 9))




#                         Henderson Hasselblach Equation 
################################################################################
# Function to calculate pH using the Henderson-Hasselbalch equation
calculate_pH <- function(pKa, A_minus_conc, HA_conc) {
  # Check that concentrations are positive
  if (A_minus_conc <= 0 || HA_conc <= 0) {
    stop("Concentrations must be positive values.")
  }
  
  # Calculate the ratio [A-]/[HA]
  ratio <- A_minus_conc / HA_conc
  
  # Calculate pH
  pH <- pKa + log10(ratio)
  
  return(pH)
}

# Given values
pKa <- 4.76
A_minus_conc <- 0.1  # [CH3COOm]
HA_conc <- 0.1       # [CH3COOH]

# Calculate pH
pH_value <- calculate_pH(pKa, A_minus_conc, HA_conc)

# Display the result
cat("The pH of the buffer solution is:", pH_value, "\n")




