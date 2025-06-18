#!/usr/bin/env Rscript

# ==============================================================================
# Automatic Navigation Generator for _quarto.yml
# ==============================================================================
# This script reads the runs_metadata.yml file and generates the navigation
# section of _quarto.yml automatically.
#
# Usage: Rscript update_navigation.R
# ==============================================================================

library(yaml)
library(glue)

# Configuration
METADATA_FILE <- "data/runs_metadata.yml"
QUARTO_CONFIG_FILE <- "_quarto.yml"
BACKUP_SUFFIX <- ".backup"

# Read metadata
cat("📚 Reading metadata from", METADATA_FILE, "\n")
if (!file.exists(METADATA_FILE)) {
  stop("❌ Metadata file not found: ", METADATA_FILE)
}

metadata <- read_yaml(METADATA_FILE)

# Read existing _quarto.yml
cat("📖 Reading existing Quarto config from", QUARTO_CONFIG_FILE, "\n")
if (!file.exists(QUARTO_CONFIG_FILE)) {
  stop("❌ Quarto config file not found: ", QUARTO_CONFIG_FILE)
}

quarto_config <- read_yaml(QUARTO_CONFIG_FILE)

# Create backup
backup_file <- paste0(QUARTO_CONFIG_FILE, BACKUP_SUFFIX)
file.copy(QUARTO_CONFIG_FILE, backup_file, overwrite = TRUE)
cat("💾 Backup created:", backup_file, "\n")

# Function to create navigation menu for a category
create_category_menu <- function(category_name, runs, display_name) {
  menu_items <- list()
  
  # Add overview item
  menu_items[[1]] <- list(
    text = "Overview",
    href = glue("reports/{category_name}/index.qmd")
  )
  
  # Add separator
  menu_items[[2]] <- list(text = "---")
  
  # Add individual run items
  for (i in seq_along(runs)) {
    run <- runs[[i]]
    menu_items[[i + 2]] <- list(
      text = run$display_name,
      href = glue("reports/{category_name}/{run$id}.qmd")
    )
  }
  
  return(list(
    text = display_name,
    menu = menu_items
  ))
}

# Generate new navigation
new_navbar_left <- list()

# Add main page
new_navbar_left[[1]] <- list(
  text = "Main",
  href = "index.qmd"
)

# Add category menus
category_info <- list(
  metagenomics = "Metagenomics",
  bacteria = "Genomic Bacteria", 
  virus = "Genomic Virus"
)

navbar_index <- 2
for (category in names(metadata$runs)) {
  if (category %in% names(category_info)) {
    new_navbar_left[[navbar_index]] <- create_category_menu(
      category, 
      metadata$runs[[category]], 
      category_info[[category]]
    )
    navbar_index <- navbar_index + 1
  }
}

# Add about page
new_navbar_left[[navbar_index]] <- list(
  text = "About",
  href = "about.qmd"
)

# Update the configuration
quarto_config$website$navbar$left <- new_navbar_left

# Write updated configuration
cat("✏️  Writing updated navigation to", QUARTO_CONFIG_FILE, "\n")
write_yaml(quarto_config, QUARTO_CONFIG_FILE)

# Print summary
total_runs <- sum(sapply(metadata$runs, length))
cat("🎉 Navigation updated successfully!\n")
cat("📊 Summary:\n")
for (category in names(metadata$runs)) {
  cat("  -", category_info[[category]], ":", length(metadata$runs[[category]]), "runs\n")
}
cat("📈 Total runs:", total_runs, "\n")
cat("💾 Backup saved as:", backup_file, "\n")
cat("🔄 Reload your Quarto preview to see the updated navigation.\n")
