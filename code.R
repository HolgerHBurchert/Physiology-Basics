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




