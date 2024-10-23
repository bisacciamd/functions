require(tidyverse)
library(REDCapTidieR)

# Function to prompt user for offline or online rendering
redcap_mode <- function() {
  mode <- readline(prompt = "Would you like to use offline (0) or online (1) rendering? (0/1): ")
  if (mode == "0") {
    return("offline")
  } else if (mode == "1") {
    return("online")
  } else {
    cat("Invalid input. Please enter 'O' for offline or 'N' for online.\n")
    return(redcap_mode())
  }
}

# Prompt user for mode
mode <- redcap_mode()

if (mode == "offline") {
  # Load data from RDS file
  d <- readRDS("data/d.rds")
} else if (mode == "online") {
  # Load data from REDCap using URI and token
  d <- read_redcap(
    redcap_uri = readLines("data/uri.secret", warn = F),
    token = readLines("data/token.secret", warn = F)
  )
  # Saves to an RDS file for offline use
  saveRDS(d, "data/d.rds")
}

# Process the data
d %>%
  make_labelled() %>%
  bind_tibbles()
