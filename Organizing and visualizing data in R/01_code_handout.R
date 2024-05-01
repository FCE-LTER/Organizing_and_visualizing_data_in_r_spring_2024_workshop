# Code handout for FCE LTER Organizing and Visualizing Data in R
# Workshop

# Author: Gabriel E. Kamener
# Last Update: 2024-03-14

# Organization: Florida Coastal Everglades LTER Program
# Information Manager: Gabriel Kamener (gkamener@fiu.edu)
# Website: https://fcelter.fiu.edu
# Site Manager: Michael Rugge (ruggem@fiu.edu)
# Lead Principal Investigator: John Kominoski (jkominos@fiu.edu)
# Website: https://fcelter.fiu.edu
# GitHub site: https://github.com/FCE-LTER


#Install needed package in next line if not already installed
#install.packages("tidyverse")

# Load packages -----------------------------------------------------------
library(tidyverse)


# Load data ---------------------------------------------------------------


# Indexing and subsetting data frames -------------------------------------

# First row, first column


# First row, all columns


# All rows, first column


# Rows 1-3, columns 5-6


# Return a vector


# First row, WATER_DEPTH_CM column


# Subset first ten rows of peri_df in a new data frame


# Repeat with last ten rows


# Factors -----------------------------------------------------------------

# Create factor levels for FLOATING_SP1 column




# Formatting dates --------------------------------------------------------

# Create date vector from character vector with ymd()


# Try with paste method and "-" as separators


# Now apply to bad_dates

# Use that as an argument for ymd()

# Create DATE column from those values

# Review structure of modified data frame

# Get summary

# Subset records with missing dates



# Data Manipulation -------------------------------------------------------

## Select ####

# Note: dpylr's select function is often masked when loading other packages
# after loading dplyr. We can preface the function with "dplyr::" to ensure
# we are using the select function from that package.

# Select columns of interest



## Filter ####

# Filter to include only September - December sampling


# Filter another way to exclude records where peri TP
# values are missing



## Mutate ####


## Mutate new column containing percent organic values


## Pivoting ####


# Pivot diatom counts to long


# Pivot back to wide


## Group and summarize ####

# summarize mean water depth by wetland basin and year



# Pivot wide water depths summary


# Filter and summarize the total diatoms counted for each sample


# Calculate diatom relative percent abundance


## Count ####

# Count number of times each floating plant species was
# number one for periphyton substrate


# Count presence of diatom taxa across samples


# Exporting data ####



# Activity: create your own new variable ----------------------------------------------------



# Visualizing data with ggplot ----------------------------------------------

## The basic ggplot template ####
#ggplot(data = <DATA>, mapping = aes(<MAPPINGS>)) +  <GEOM_FUNCTION>()




## Scatter plots ####


## Boxplots ####


## Trend lines ####


## Stacked bar plot ####


# Filter the relative abundance data for 2014 samples
# in the W3B wetland basin

# Create stacked bar plot


## Customized plots ####


# Create a stacked bar plot with custom labels for axes and legend 



# Group any taxa that make up less than 2% relative abundance
# into an "OTHER" group


# Create a stacked bar plot with custom labels for axes and legend
# using cutoff data



## Save plot to file ####


# Save with default size


# Save with custom size


# Activity: create your own plot ----------------------------------------------------




