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
peri_df <- read_csv("data/raw/FCE1210_CERP_Periphyton.csv",
                    na = "-9999")

diatoms_df <- read_csv("data/raw/diatom_counts.csv",
                       na = "NA")

bad_dates_df <- read_csv("data/raw/bad_dates.csv",
                         na = "-9999")


# Indexing and subsetting data frames -------------------------------------

# First row, first column
peri_df[1, 1]

# First row, all columns
peri_df[1, ]

# All rows, first column
peri_df[, 1]

# Rows 1-3, columns 5-6
peri_df[c(1, 2, 3), c(5, 6)] 

peri_df[1:3, 5:6] 

# Return a vector
peri_df[[1, 1]]

# First row, WATER_DEPTH_CM column
peri_df[1, "WATER_DEPTH_CM"]

# Subset first ten rows of peri_df in a new data frame
first_ten_rows <- head(peri_df, 10)

first_ten_rows

# Repeat with last ten rows
last_ten_rows <- tail(peri_df, 10)

view(last_ten_rows)

# Factors -----------------------------------------------------------------

# Create factor levels for FLOATING_SP1 column
peri_df$FLOATING_SP1 <- factor(peri_df$FLOATING_SP1)

nlevels(peri_df$FLOATING_SP1)

levels(peri_df$FLOATING_SP1)


# Formatting dates --------------------------------------------------------

# Create date vector from character vector with ymd()
my_date <- ymd("2015-01-01")

str(my_date)

# Try with paste method and "-" as separators
my_date <- ymd(paste("2015", "1", "1", sep = "-")) 

str(my_date)

# Now apply to bad_dates
paste(bad_dates_df$YEAR, bad_dates_df$MONTH, bad_dates_df$DAY, sep = "-")

# Use that as an argument for ymd()
ymd(paste(bad_dates_df$YEAR, bad_dates_df$MONTH, bad_dates_df$DAY, sep = "-"))

# Create DATE column from those values
bad_dates_df$DATE <- ymd(paste(bad_dates_df$YEAR, bad_dates_df$MONTH, bad_dates_df$DAY, sep = "-"))

# Review structure of modified data frame
str(bad_dates_df)

# Get summary
summary(bad_dates_df$DATE)

# Subset records with missing dates
missing_dates <- bad_dates_df[is.na(bad_dates_df$DATE), c("YEAR", "MONTH", "DAY")]

head(missing_dates)


# Data Manipulation -------------------------------------------------------

## Select ####

# Note: dpylr's select function is often masked when loading other packages
# after loading dplyr. We can preface the function with "dplyr::" to ensure
# we are using the select function from that package.

# Select columns of interest
selected_peri <- peri_df %>%
  dplyr::select(TAG_ID,
                OBS_DATE,
                PRIMARY_SAMPLING_UNIT,
                OBS_DATE,
                FIELD_REPLICATE,
                EPISODE,
                WETLAND_BASIN,
                WATER_DEPTH_CM,
                FLOATING_SP1,
                PERI_AFDM_G_PER_M2,
                PERI_TP_UG_PER_G_DRY_MASS,
                PERI_PROP_ORGANIC
  )


## Filter ####

# Filter to include only September - December sampling
filtered_peri <- selected_peri %>%
  filter(month(OBS_DATE) > 8)

# Filter another way to exclude records where peri TP
# values are missing
filter_missing_peri_tp <- selected_peri %>%
  filter(!is.na(PERI_TP_UG_PER_G_DRY_MASS))


## Mutate ####

## Mutate new column containing percent organic values
mutated_peri <- filtered_peri %>%
  mutate(peri_percent_organic = PERI_PROP_ORGANIC*100)

mutated_peri %>%
  select(PERI_PROP_ORGANIC,
         peri_percent_organic)

## Pivoting ####

# Pivot diatom counts to long
diatoms_counts_long <- diatoms_df %>%
  pivot_longer(-c(TAG_ID:LSU_NAME),
               values_to = "SPECIMENS_COUNTED",
               names_to = "TAXON_CODE")

# Pivot back to wide
diatom_counts_wide <- diatoms_counts_long %>%
  pivot_wider(names_from = "TAXON_CODE",
              values_from = "SPECIMENS_COUNTED",
              names_sort = TRUE,
              values_fill = 0)


## Group and summarize ####

# summarize mean water depth by wetland basin and year
summarized_water_depths <- filtered_peri %>%
  group_by(WETLAND_BASIN,
           YEAR = year(OBS_DATE)) %>%
  summarize(mean_water_depth_cm = mean(WATER_DEPTH_CM))


# Pivot wide water depths summary
summarized_water_depths_wide <- summarized_water_depths %>%
  pivot_wider(names_from = WETLAND_BASIN,
              values_from = mean_water_depth_cm)

summarized_water_depths_wide

# Filter and summarize the total diatoms counted for each sample
diatoms_counts_filtered <- diatoms_counts_long %>%
  filter(!is.na(SPECIMENS_COUNTED)
         & SPECIMENS_COUNTED != 0)

total_diatoms <- diatoms_counts_long %>%
  group_by(TAG_ID) %>%
  summarize(TOTAL_COUNT = sum(SPECIMENS_COUNTED))

total_diatoms

# Calculate diatom relative percent abundance
rel_abund <- diatoms_counts_filtered %>%
  # Using dplyr's left_join function to join with the total_diatoms data frame
  left_join(., total_diatoms, by = "TAG_ID") %>%
  mutate(RELATIVE_PCT_ABUND = SPECIMENS_COUNTED/TOTAL_COUNT*100) %>%
  select(TAG_ID,
         OBS_DATE,
         PRIMARY_SAMPLING_UNIT,
         WETLAND_BASIN,
         TAXON_CODE,
         RELATIVE_PCT_ABUND
  )

rel_abund

## Count ####

# Count number of times each floating plant species was
# number one for periphyton substrate
filtered_peri %>%  count(FLOATING_SP1, sort = TRUE)


# Count presence of diatom taxa across samples
diatoms_counts_filtered %>%
  count(TAXON_CODE, sort = TRUE)


# Exporting data ####

write_csv(filtered_peri, "data/intermediate/filtered_peri.csv")

write_csv(summarized_water_depths_wide, "data/final/summarized_water_depths_wide.csv")


# Activity: create your own new variable ----------------------------------------------------



# Visualizing data with ggplot ----------------------------------------------

## The basic ggplot template ####
#ggplot(data = <DATA>, mapping = aes(<MAPPINGS>)) +  <GEOM_FUNCTION>()

ggplot(data = filtered_peri)

ggplot(data = filtered_peri, mapping = aes(x = WATER_DEPTH_CM, y = PERI_AFDM_G_PER_M2))


## Scatter plots ####
ggplot(data = filtered_peri, mapping = aes(x = WATER_DEPTH_CM, y = PERI_AFDM_G_PER_M2)) +
  geom_point()


## Boxplots ####
ggplot(data = filtered_peri, mapping = aes(x = WETLAND_BASIN, y = WATER_DEPTH_CM)) +
  geom_boxplot()


## Trend lines ####
ggplot(data = filtered_peri, mapping = aes(x = year(OBS_DATE), y = WATER_DEPTH_CM)) +
  geom_line()

ggplot(data = summarized_water_depths, mapping = aes(x = YEAR, y = mean_water_depth_cm, color = WETLAND_BASIN)) +
  geom_line() +
  scale_x_continuous(breaks = c(2005:2014))

ggplot(data = summarized_water_depths, mapping = aes(x = YEAR, y = mean_water_depth_cm, color = WETLAND_BASIN)) +
  geom_line() +
  scale_x_continuous(breaks = c(2005:2014))


## Stacked bar plot ####


# Filter the relative abundance data for 2014 samples
# in the W3B wetland basin
rel_abund_2014_W3B <- rel_abund %>%
  filter(year(OBS_DATE) == 2014
         & WETLAND_BASIN == "W3B")

# Create stacked bar plot
ggplot(rel_abund_2014_W3B,
       aes(fill = TAXON_CODE,
           y = RELATIVE_PCT_ABUND,
           x = PRIMARY_SAMPLING_UNIT)) + 
  geom_bar(position="stack", stat="identity")


## Customized plots ####


# Create a stacked bar plot with custom labels for axes and legend 
ggplot(rel_abund_2014_W3B,
       aes(fill = TAXON_CODE,
           y = RELATIVE_PCT_ABUND,
           x = PRIMARY_SAMPLING_UNIT)) + 
  geom_bar(position="stack", stat="identity") +
  labs(title = "Relative % abundance at W3B sites in 2014",
       x = "Primary Sampling Unit",
       y = "Relative Percent Abundance",
       fill = "Diatom taxon code") +
  theme(axis.title.x = element_text(size = 11),
        axis.title.y = element_text(size = 15))


# Group any taxa that make up less than 2% relative abundance
# into an "OTHER" group
rel_abund_2014_w3b_cut <- rel_abund %>%
  filter(year(OBS_DATE) == 2014
         & WETLAND_BASIN == "W3B") %>%  
  mutate(TAXON_CODE = if_else(RELATIVE_PCT_ABUND < 2,
                              "OTHER",
                              TAXON_CODE)) %>%
  group_by(PRIMARY_SAMPLING_UNIT,
           TAXON_CODE) %>%
  summarize(RELATIVE_PCT_ABUND = sum(RELATIVE_PCT_ABUND))

# Create a stacked bar plot with custom labels for axes and legend
# using cutoff data
plot_2014_w3b_cut <- ggplot(rel_abund_2014_w3b_cut,
                            aes(fill = TAXON_CODE,
                                y = RELATIVE_PCT_ABUND,
                                x = PRIMARY_SAMPLING_UNIT)) +
  geom_bar(position="stack", stat="identity") +
  labs(title = "Relative % abundance at W3B sites in 2014",
       x = "Primary Sampling Unit",
       y = "Relative Percent Abundance",
       fill = "Diatom taxon code") +
  theme(plot.title = element_text(size = 15), 
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))


plot_2014_w3b_cut <- ggplot(rel_abund_2014_w3b_cut,
                            aes(fill = TAXON_CODE,
                                y = RELATIVE_PCT_ABUND,
                                x = PRIMARY_SAMPLING_UNIT)) +
  geom_bar(position="stack", stat="identity") +
  labs(title = "Relative % abundance at W3B sites in 2014",
       x = "Primary Sampling Unit",
       y = "Relative Percent Abundance",
       fill = "Diatom taxon code") +
  theme(plot.title = element_text(size = 15),
        axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12))


plot_2014_w3b_cut

## Save plot to file ####

# Save with default size
ggsave("fig/W3B_2014_relative_abundance.png",
       plot = plot_2014_w3b_cut)

# Save with custom size
ggsave("fig/W3B_2014_relative_abundance_custom_size.png",
       plot = plot_2014_w3b_cut,
       width = 46,
       height = 26,
       units = "cm",
       limitsize = FALSE)


# Activity: create your own plot ----------------------------------------------------




