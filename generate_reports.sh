#!/bin/bash

# ==============================================================================
# Complete FOI Reports Setup Script
# ==============================================================================
# This script generates all individual reports and updates navigation
#
# Usage: bash setup_all.sh
# ==============================================================================

echo "🚀 FOI Reports - Complete Setup"
echo "================================"
echo ""

# Check if required files exist
if [ ! -f "data/runs_metadata.yml" ]; then
    echo "❌ Error: data/runs_metadata.yml not found"
    echo "   Please create the metadata file first"
    exit 1
fi

if [ ! -f "_quarto.yml" ]; then
    echo "❌ Error: _quarto.yml not found"
    echo "   Please make sure you're in the Quarto project directory"
    exit 1
fi

# Create necessary directories
echo "📁 Creating directory structure..."
mkdir -p data
mkdir -p reports/metagenomics
mkdir -p reports/bacteria
mkdir -p reports/virus
mkdir -p assets/images

echo "📄 Step 1: Generating individual reports..."
echo "==========================================="
Rscript generate_reports.R

if [ $? -ne 0 ]; then
    echo "❌ Error generating reports"
    exit 1
fi

echo ""
echo "🧭 Step 2: Updating navigation..."
echo "================================="
Rscript update_navigation.R

if [ $? -ne 0 ]; then
    echo "❌ Error updating navigation"
    exit 1
fi

echo ""
echo "✨ Step 3: Creating category overview pages..."
echo "=============================================="

# Create overview pages for each category
echo "📝 Creating metagenomics overview..."
cat > reports/metagenomics/index.qmd << 'EOF'
---
title: "Metagenomics Reports"
description: "Overview of all metagenomics sequencing runs and analysis reports"
---

## Overview

This section contains sequencing reports and analyses for metagenomics studies conducted at FOI. Each report represents a complete sequencing run with associated metadata and analysis results.

## Available Reports

```{r}
#| echo: false
#| warning: false

library(yaml)
library(knitr)
library(dplyr)

# Read metadata
metadata <- read_yaml("../../data/runs_metadata.yml")
metagenomics_runs <- metadata$runs$metagenomics

# Create table of runs
runs_table <- data.frame(
  `Run ID` = sapply(metagenomics_runs, function(x) x$display_name),
  `Date` = sapply(metagenomics_runs, function(x) x$date),
  `Type` = sapply(metagenomics_runs, function(x) x$type),
  `Samples` = sapply(metagenomics_runs, function(x) x$sample_number),
  `Description` = sapply(metagenomics_runs, function(x) substr(x$description, 1, 50)),
  `Report` = sapply(metagenomics_runs, function(x) 
    paste0("[View Report](", x$id, ".qmd)")),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

kable(runs_table, format = "markdown")
```

## Recent Activity

The most recent metagenomics analysis was completed on `r max(sapply(metagenomics_runs, function(x) x$date))`.

## Analysis Pipeline

Our metagenomics pipeline includes:

- **Quality Control**: FastQC, MultiQC
- **Taxonomic Classification**: Kraken2, Bracken  
- **Functional Analysis**: HUMAnN3
- **Diversity Analysis**: QIIME2
- **Assembly**: MEGAHIT, MetaSPAdes
- **Binning**: MetaBAT2, CONCOCT

## Contact

For questions about metagenomics analyses, contact the FOI Bioinformatics Team.
EOF

echo "📝 Creating bacteria overview..."
cat > reports/bacteria/index.qmd << 'EOF'
---
title: "Bacterial Genomics Reports"
description: "Overview of all bacterial genomics sequencing runs and analysis reports"
---

## Overview

This section contains sequencing reports and analyses for bacterial genomics studies conducted at FOI. Each report represents whole genome sequencing of bacterial isolates with comprehensive genomic characterization.

## Available Reports

```{r}
#| echo: false
#| warning: false

library(yaml)
library(knitr)
library(dplyr)

# Read metadata
metadata <- read_yaml("../../data/runs_metadata.yml")
bacteria_runs <- metadata$runs$bacteria

# Create table of runs
runs_table <- data.frame(
  `Run ID` = sapply(bacteria_runs, function(x) x$display_name),
  `Date` = sapply(bacteria_runs, function(x) x$date),
  `Type` = sapply(bacteria_runs, function(x) x$type),
  `Samples` = sapply(bacteria_runs, function(x) x$sample_number),
  `Description` = sapply(bacteria_runs, function(x) substr(x$description, 1, 50)),
  `Report` = sapply(bacteria_runs, function(x) 
    paste0("[View Report](", x$id, ".qmd)")),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

kable(runs_table, format = "markdown")
```

## Recent Activity

The most recent bacterial genomics analysis was completed on `r max(sapply(bacteria_runs, function(x) x$date))`.

## Analysis Pipeline

Our bacterial genomics pipeline includes:

- **Quality Control**: FastQC, Trimmomatic
- **Assembly**: SPAdes
- **Annotation**: Prokka
- **Species Identification**: ANI calculation
- **AMR Detection**: ABRicate, ResFinder
- **Phylogenetics**: Core genome analysis
- **Typing**: MLST, cgMLST

## Contact

For questions about bacterial genomics analyses, contact the FOI Bioinformatics Team.
EOF

echo "📝 Creating virus overview..."
cat > reports/virus/index.qmd << 'EOF'
---
title: "Viral Genomics Reports"
description: "Overview of all viral genomics sequencing runs and analysis reports"
---

## Overview

This section contains sequencing reports and analyses for viral genomics studies conducted at FOI. Each report represents viral genome sequencing with consensus assembly and variant analysis.

## Available Reports

```{r}
#| echo: false
#| warning: false

library(yaml)
library(knitr)
library(dplyr)

# Read metadata
metadata <- read_yaml("../../data/runs_metadata.yml")
virus_runs <- metadata$runs$virus

# Create table of runs
runs_table <- data.frame(
  `Run ID` = sapply(virus_runs, function(x) x$display_name),
  `Date` = sapply(virus_runs, function(x) x$date),
  `Type` = sapply(virus_runs, function(x) x$type),
  `Samples` = sapply(virus_runs, function(x) x$sample_number),
  `Description` = sapply(virus_runs, function(x) substr(x$description, 1, 50)),
  `Report` = sapply(virus_runs, function(x) 
    paste0("[View Report](", x$id, ".qmd)")),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

kable(runs_table, format = "markdown")
```

## Recent Activity

The most recent viral genomics analysis was completed on `r max(sapply(virus_runs, function(x) x$date))`.

## Analysis Pipeline

Our viral genomics pipeline includes:

- **Quality Control**: FastQC, Trimmomatic
- **Host Removal**: BWA-MEM alignment
- **Viral Detection**: BLAST, custom databases
- **Consensus Assembly**: iVar
- **Variant Calling**: LoFreq, iVar
- **Phylogenetics**: IQ-TREE, NextStrain

## Contact

For questions about viral genomics analyses, contact the FOI Bioinformatics Team.
EOF

echo ""
echo "🎉 SETUP COMPLETE!"
echo "=================="
echo ""
echo "✅ Generated individual reports for all runs"
echo "✅ Updated navigation in _quarto.yml"  
echo "✅ Created category overview pages"
echo ""
echo "📁 Directory structure:"
echo "   ├── data/runs_metadata.yml"
echo "   ├── reports/"
echo "   │   ├── metagenomics/ ($(ls reports/metagenomics/*.qmd 2>/dev/null | wc -l) files)"
echo "   │   ├── bacteria/ ($(ls reports/bacteria/*.qmd 2>/dev/null | wc -l) files)"
echo "   │   └── virus/ ($(ls reports/virus/*.qmd 2>/dev/null | wc -l) files)"
echo ""
echo "🚀 Next steps:"
echo "   1. Run: quarto preview"
echo "   2. Check your navigation menus"
echo "   3. Customize individual reports as needed"
echo ""
echo "💡 To add new runs:"
echo "   1. Update data/runs_metadata.yml"
echo "   2. Run: bash setup_all.sh"
echo ""