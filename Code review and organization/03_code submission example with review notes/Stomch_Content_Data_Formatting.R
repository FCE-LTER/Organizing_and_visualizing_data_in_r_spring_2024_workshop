library(tidyverse)
library(openxlsx)

#----------Load data---------------------------------
#set wd for data input
setwd("~/FIU/Research/Current Food Web Study/Data/Stomach Contents")

#read in data
stomachs.raw <- read.xlsx("dpe_stomach_contents_raw_data.xlsx")
stomachs.raw <- stomachs.raw[-c(4),] #remove MC28 duplicate empty stomach 
#output working directory
setwd("~/FIU/Research/Current Food Web Study/Analyses/Stomach Contents")

#look at prey groups a combine if needed to
prey.groups <- levels(as.factor(stomachs.raw$Diet.Item)) %>% as.data.frame()
#View(prey.groups)
#write.xlsx(prey.groups, "DPE_prey_groups.xlsx")


#------Condense Prey Groups--------------------------------------------------------------------
#Amphipoda
stomachs.raw$Diet.Item[stomachs.raw$Diet.Item %in% c("Amphipoda ", "Amphipoda  ")] <- "Amphipoda"

#Belo/Naucor
stomachs.raw$Diet.Item[stomachs.raw$Diet.Item %in% c("Belostoma", "Pelocoris")] <- "Belo/Naucor"

#Centrarchid
stomachs.raw$Diet.Item[stomachs.raw$Diet.Item %in% c("Elasomma evergladei")] <- "Centrarchid"


#"Chironomid", "Chironomid ", and "Chironomid larvae" need to all be "Chironomid"
stomachs.raw$Diet.Item[stomachs.raw$Diet.Item %in% c("Chironomid ", "Chironomid larvae")] <- "Chironomid"

#Coleoptera
stomachs.raw$Diet.Item[stomachs.raw$Diet.Item %in% c("Hydrilla", "Cybister")] <- "Coleoptera"

#Copepoda
stomachs.raw$Diet.Item[stomachs.raw$Diet.Item %in% c("Copepod")] <- "Copepoda"

#Cyprinodontid
stomachs.raw$Diet.Item[stomachs.raw$Diet.Item %in% 
                         c("Gambusia holbrooki", "Gambusia hoolbroki", 
                           "Heterandria formosa", "Lucania goodei")] <- "Cyprinodontid"

#Detritus
stomachs.raw$Diet.Item[stomachs.raw$Diet.Item %in% c(" Detritus", "Vasc Detritus", "Vascular Detritus")] <- "Detritus"

#Diptera
stomachs.raw$Diet.Item[stomachs.raw$Diet.Item %in% c("Culicidae")] <- "Diptera"

#Hemiptera
stomachs.raw$Diet.Item[stomachs.raw$Diet.Item %in% c("Corixidae", "Gerridae", "Hebridae", "Mesovelidae", 
                                                     "Nepidae", "Notonectidae", "Lubridae")] <- "Hemiptera"

#Hydrachnidia - Aquatic Mites
stomachs.raw$Diet.Item[stomachs.raw$Diet.Item %in% 
                         c("Blue Mite", "Brown Mite", "Mite", "Orange Mite", "Red Mite")] <- "Hydrachnidia"

#Hymenoptera
stomachs.raw$Diet.Item[stomachs.raw$Diet.Item %in% c("Formicidae", "Halictidae")] <- "Hymenoptera"


#Misc Invert
stomachs.raw$Diet.Item[stomachs.raw$Diet.Item %in% c("Misc Insecta", "Flatworm", "Plecoptera")] <- "Misc Invert"

#Mollusca
stomachs.raw$Diet.Item[stomachs.raw$Diet.Item %in% c("Bivalvia", "Haita spp", "Melanoides", "Mollusca",
                                                     "Physid Snail", "Planorbella duryi", "PSECOL",
                                                     "Snail  ", "Snail Physid")] <- "Mollusca"

#Odonata
stomachs.raw$Diet.Item[stomachs.raw$Diet.Item %in% c("Anisoptera", "Coenagridae", "Zygoptera")] <- "Odonata"

#Oligochaeta
stomachs.raw$Diet.Item[stomachs.raw$Diet.Item %in% c("Annelida")] <- "Oligochaeta"

#Vascular Plants
stomachs.raw$Diet.Item[stomachs.raw$Diet.Item %in% c("Utricularia gibba")] <- "Vasc Plant"

#---------reformat data and attach covariates---------------------------------------------------
#prey individual sum like prey items
#some individuals now have duplicate entries per diet item after condensing
stomachs.counts.condensed <- stomachs.raw %>% 
  select(ID, Diet.Item, Count) %>% 
  group_by(ID, Diet.Item) %>% summarise(Count = sum(Count))

#pivot wider for count data
stomachs.counts <- stomachs.counts.condensed %>% filter(Diet.Item != "-") %>% 
  pivot_wider(id_cols = ID, names_from = Diet.Item, values_from = Count, values_fill = 0) %>% 
  select(order(colnames(.))) %>% 
  select(ID, everything())

#do the same for volume data
#sum like diet items per individual
stomachs.volume.condensed <- stomachs.raw %>% 
  mutate_at(vars(Volume), funs(as.numeric)) %>% 
  select(ID, Diet.Item, Volume) %>% 
  group_by(ID, Diet.Item) %>% summarise(Volume = sum(Volume))

#pivot wider for volume data
stomachs.volume <- stomachs.volume.condensed %>% filter(Diet.Item != "-") %>% 
  pivot_wider(id_cols = ID, names_from = Diet.Item, values_from = Volume) %>% 
  mutate_at(vars(-c("ID")), funs(as.numeric))


#---------------------------------attach covariates------------------------------------------------------------------
#change wd to read in processed samples dataframe
setwd("~/FIU/Research/Current Food Web Study/Data/Processed Samples")

#read in processed samples data
processed.samples <- read.csv(file = "processed_samples_qaqc.csv")

#isolate sample that we dissected stomachs for
processed.stomachs <- processed.samples %>% filter(ID %in% stomachs.counts$ID) %>% distinct(ID, .keep_all = T) 

#we dissected two stomachs that don't appear in the processed samples datasheet
missing.stomachs <- stomachs.counts %>% filter(!ID %in% processed.stomachs$ID)
#they also do not appear in the physical datasheets
#omitting them from analyses

#merge stomach data with processed sample data to add covariates
stomachs.merge <- merge(processed.stomachs, stomachs.counts, by = "ID")

#----------------Export to Excel-------------------------------------------------------------------------------------
#output working directory
setwd("~/FIU/Research/Current Food Web Study/Analyses/Stomach Contents")

#trim off some columns write data out to Excel
stomachs.merge %>% select(-c(Sample_Type, Month, Day, Method, Throw, Stomach_Kept, Comments, Freezer, Tube_Weight, Total_Mass, Sorted,
                             neg_80_freezer, Freeze_Dried, Dried, Crushed, Isotope, Decarb, Entered_by, Checked_by, Count)) %>% 
                  select(ID, size_class, Species, everything()) %>% 
                write.xlsx("stomachs_count_merge.xlsx")
