library(tidyverse)
library(vegan)
library(openxlsx)
library(ggrepel)
library(RVAideMemoire)
library(pander)
library(ggpubr)
library(ggbreak)
library(vegan)

############load post-invasion data and set output working directory#####################################################################
#set wd to folder with merged stomach content data
#setwd("G:/Other Computers/My Laptop/Documents/FIU/Research/Current Food Web Study/Analyses/Stomach Contents")
setwd("~/FIU/Research/Current Food Web Study/Analyses/Stomach Contents")

#read in data
stomach.counts <- read.xlsx("stomachs_count_merge.xlsx")

#set output wd
#setwd("G:/Other Computers/My Laptop/Documents/FIU/Research/Current Food Web Study/Analyses/Stomach Contents/Diet")
setwd("~/FIU/Research/Current Food Web Study/Analyses/Stomach Contents/Diet")

#-----sample size tables----------------------------------------------------------------------------------------------------------------------------------------
sample.sizes.wide <- stomach.counts %>% group_by(size_class, Slough, Habitat, Season) %>% tally() %>% 
  pivot_wider(names_from = c(Slough, Habitat, Season), values_from = n) %>% 
  select(size_class, SRS_Marsh_Wet, SRS_NP_Wet, SRS_Pond_Wet,
                     SRS_Marsh_Dry, SRS_NP_Dry, SRS_Pond_Dry,
                     TSL_Marsh_Wet, TSL_NP_Wet, TSL_Pond_Wet,
                     TSL_Marsh_Dry, TSL_Pond_Dry)

write.xlsx(sample.sizes, "stomach_content_sample_sizes.xlsx")

sample.sizes.long <- stomach.counts %>% group_by(size_class, Slough, Habitat, Season) %>% tally() %>% 
  arrange(Slough, Season, Habitat, size_class)

paired.loftus.sample.sizes <- stomach.counts %>% 
  filter(Site %in% c ("BJB", "LJB", "AHCFPB", "6A", "6B", "6C")) %>% 
  group_by(size_class) %>% 
  tally()
  
srs.sample.sizes <- stomach.counts %>% 
  filter(Slough == "SRS") %>% 
  group_by(size_class) %>% 
  tally()


#-----format data for multivariate analyses----------------------------------------------------------------------------------------------------------------------
#remove covariates besides size class, habitat, season, and slough and remove empty stomachs
no.empty.stomachs <- stomach.counts %>% 
  filter(Empty == 0) %>% 
  #select(-c(Species, Year, Site, Length, Wet_Weight, Sex, Empty)) %>% 
  remove_rownames() %>% 
  column_to_rownames(var = "ID")
#this removes any empty rows

#isolate prey cols
prey.com <- no.empty.stomachs %>% select(Amphipoda:Vasc.Plant)

#look at column sums to assess rare prey groups
sort(apply(prey.com, 2, sum), decreasing = F)

#remove rare prey groups
prey.no.rare <- prey.com %>% select(-c(Argulus, Centrarchid, Vasc.Plant, Oligochaeta))

#drop zero rows that are created by removing rare prey items
#get row sums
row.sums <- apply(prey.no.rare, 1, sum)
#create T/F vector for where rows > 0
no.zero.rows <- row.sums > 0
#select for non-zero rows
prey.no.rare.no.zero <- prey.no.rare[no.zero.rows,]

post.covariates <- no.empty.stomachs %>% rownames_to_column(var = "ID") %>% select(ID:Habitat)
post.covariates.no.zero <- post.covariates[no.zero.rows,]
#log transform the data
#prey.log.com <- log(prey.no.rare.no.zero+1)
#------Run NMDS and calculate estimated variance explained-----------------------------------------------------------------------------------
set.seed(999)
post.nmds <- metaMDS(prey.no.rare.no.zero, k=2, distance = "horn", noshare = 0.1, autotransform = F, plot = T)
post.nmds

#quick visulization
ordiplot(post.nmds, type = "text", main = sprintf("k = 2; stress = %.3f", post.nmds$stress))

#calculate stress values by hand
#reconstruct community distance as in the nmds, accounting for zeros
comm.dist <- metaMDSredist(post.nmds, distance = "horn")
#calculate ordination distance
ord.dist <- vegdist(post.nmds$points[,1:2], method = "euclidean")
#extract residuals from a monotonic regression
dhat <- isoreg(comm.dist, ord.dist)
#plot monotonic regression
plot(dhat)
#calculate stress w/ residuals of this regression and ordination distances
reg.stress <- sqrt(residuals(dhat)^2/sum(ord.dist^2))
stress.label <- sprintf("Stress = %.2f", reg.stress*100)
stressplot(post.nmds)
text(1.5, 0.5, stress.label)

#----Simulation Stress Test-------------------------------------------------------------------------------------------
reps <- 1000
stressTest <- oecosimu(comm = round(prey.no.rare.no.zero), method = "quasiswap_count", 
                       nestfun = metaMDS, autotransform = FALSE, k = 2, noshare = 0.1,
                       distance = "horn", nsimul = reps, parallel = 2, statistic = "stress", 
                       alternative = "less", trace = TRUE, maxit = 1000,
                       trymax = 20, sratmax = 0.9999999)


#Take a look at the results
stressTest
stressTest$oecosimu$statistic
stressTest$oecosimu$simulated
stressTest$oecosimu$means
stressTest$oecosimu$z
stressTest$oecosimu$pval

#Plot the results using base graphics
hist(as.vector(stressTest$oecosimu$simulated), xlim = c(0,max(stressTest$oecosimu$simulated)+.05), xlab = "Ecological null model stress value",ylab = "Frequency of stress value", main = "", breaks = 7)
abline(v = stressTest$oecosimu$statistic, col = "red", lty = 2) 

#This result indicates that we can reject null hypothesis that observed stress
#is less than the simulated stress values

#Replot the results using more elegant graphics in ggplot2
simVector<-as.vector(stressTest$oecosimu$simulated)
simVector<-as.data.frame(simVector)
simVector<-cbind(simVector,seq(1:1000))

A<-ggplot(simVector,aes(x=simVector))
A<-A + geom_histogram(col="black",size=.25)
A<-A + theme_minimal()
A<-A + xlim(c(0.1, .32)) + ylim(c(0,500))
A<-A + geom_segment(aes(x = stressTest$oecosimu$statistic, y = 0, xend = stressTest$oecosimu$statistic), yend = 450,linetype=2,color = "red",size=.5)
A<-A + labs(x="Stress value")
A<-A + labs(y="Frequency observed")
A<-A + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor= element_blank(), axis.line = element_line(colour = "black",size=.25))
A<-A + scale_y_continuous(expand = c(0, 0),limits=c(0,100))
A<-A + scale_x_continuous(limits=c(.2,0.30))
A<-A + theme(axis.ticks = element_line(colour = "black", size = .25))
A<-A + theme(axis.text=element_text(colour="black"))
A<-A + theme(text=element_text(size=12))
A

(stress.sim.plot <- ggplot(simVector,aes(x=simVector))+
  geom_histogram(col="black",size=.25)+
  geom_segment(aes(x = stressTest$oecosimu$statistic, y = 0, xend = stressTest$oecosimu$statistic), 
                yend = 450,linetype=2,color = "red",size=.5)+
  annotate(geom = "text", label = "p < 0.001", x = 0.2325, y = 150)+
  annotate(geom = "text", label = "Stress = 0.227", x = 0.234, y = 160)+
  labs(x="Stress value", y="Frequency observed")+
    theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor= element_blank(), 
                       axis.line = element_line(colour = "black",size=.25))+
    theme(axis.ticks = element_line(colour = "black", size = .25), axis.text=element_text(colour="black"),
          text=element_text(size=12))

)

ggsave(filename = "stress_sim_plot.tiff", plot = stress.sim.plot, dpi  = 600, units = "mm", height = 150, width = 200)

#------------NMDS output format and plot------------------------------------------------------------------------------
post.nmds.consumers <- as.data.frame(post.nmds$points) %>% rownames_to_column(var = "ID") %>% 
  merge(., post.covariates.no.zero, by = "ID")
#write this out to Excel for future use in SIBER
write.xlsx(post.nmds.consumers, "post_nmds_consumers.xlsx")

#centroids per size class
post.nmds.spp.centroid <- post.nmds.consumers %>% group_by(size_class) %>% 
  summarise(MDS1 = mean(MDS1), MDS2 = mean(MDS2))

post.nmds.prey <- as.data.frame(post.nmds$species) %>% rownames_to_column(var = "Prey Group")

(post.nmds.plot <- ggplot()+
  #stat_ellipse(data = post.nmds.consumers, aes(x= MDS1, y = MDS2, color = Slough))+
  geom_text_repel(data = post.nmds.spp.centroid, aes(x = MDS1, y = MDS2, label = size_class),
                  max.overlaps = 17, size = 3)+
  geom_text_repel(data = post.nmds.prey, aes(x = MDS1, y = MDS2, label = `Prey Group`, color = "red"), size = 3)+
    coord_fixed()+
    theme_bw()+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(color = 'black'),
          axis.title = element_text(size = 10), axis.text = element_text(size = 8),
          #legend.title = element_text(size = 20), legend.text = element_text(size = 20),
          legend.position = "none")
)

#come back to this later - finalize plot before saving
ggsave(filename = "nmds_consumers.tiff", plot = post.nmds.plot, dpi = 600, width = 140, height = 140, units = "mm")

#------Post-Invasion SIMPER--------------------------------------------------------------------------------------------- 
habitat.season.slough <- post.covariates.no.zero %>%  unite("habdate_slough", c(Habitat, Season, Slough), sep = "_") %>% 
  select(habdate_slough)
habitat.season.slough <- as.factor(habitat.season.slough$habdate_slough)

#run SIMPER
post.simper <- simper(prey.no.rare.no.zero, habitat.season.slough)  
#results summary
post.sim.sum <- summary(post.simper) %>% lapply(., as.data.frame)
#export to Excel
write.xlsx(post.sim.sum, file = "post_simper.xlsx", row.names = T)

#determine prey groups that contribute to top 75% of diets in each habitat-season-slough
#isolate results per habitat-season-slough level
#SRS
srs.pond.wet <- post.simper$Pond_Wet_SRS_Marsh_Wet_SRS$ava
srs.pond.dry <- post.simper$Pond_Wet_SRS_Pond_Dry_SRS$avb
srs.marsh.wet <- post.simper$Marsh_Wet_SRS_Marsh_Dry_SRS$ava
srs.marsh.dry <- post.simper$Marsh_Dry_SRS_Marsh_Dry_TSL$ava
srs.np.wet <- post.simper$NP_Wet_SRS_NP_Dry_SRS$ava
srs.np.dry <- post.simper$NP_Dry_SRS_Pond_Dry_TSL$ava
#TSL
tsl.pond.wet <- post.simper$Pond_Wet_TSL_Marsh_Dry_SRS$ava
tsl.pond.dry <- post.simper$Pond_Wet_SRS_Pond_Dry_TSL$avb
tsl.marsh.wet <- post.simper$Marsh_Wet_TSL_NP_Dry_SRS$ava
tsl.marsh.dry <- post.simper$Marsh_Dry_TSL_Marsh_Wet_TSL$ava
tsl.np.wet <- post.simper$NP_Wet_TSL_Marsh_Wet_TSL$ava
#no tsl.np.dry because the landscape had dried and near pond habitat in TSL did not have water

#--------SRS reformat SIMPER results---------------------------------------------------------------------
#srs-pond-wet
#calculate percent contribution per prey group and reorder based on Per_Contri
srs.pond.wet.abun <- data.frame(srs.pond.wet) %>% mutate("Per_Contri" = srs.pond.wet/sum(srs.pond.wet)) %>% 
  rename("Av_Abun" = srs.pond.wet) %>% arrange(desc(Per_Contri)) %>% mutate("C_sum" = cumsum(Per_Contri)) %>% 
  mutate("Habitat_Season_Slough" = "srs.pond.wet") %>% rownames_to_column("Prey Group")
#filter for those that make up 75% 
srs.pond.wet.abun75 <- srs.pond.wet.abun %>% filter(C_sum < 0.76)
#srs-np-wet-
#calculate percent contribution per prey group and reorder based on Per_Contri
srs.np.wet.abun <- data.frame(srs.np.wet) %>% mutate("Per_Contri" = srs.np.wet/sum(srs.np.wet)) %>% 
  rename("Av_Abun" = srs.np.wet) %>% arrange(desc(Per_Contri)) %>% mutate("C_sum" = cumsum(Per_Contri)) %>% 
  mutate("Habitat_Season_Slough" = "srs.np.wet") %>% rownames_to_column("Prey Group")
#filter for those that make up 75% 
srs.np.wet.abun75 <- srs.np.wet.abun %>% filter(C_sum < 0.76)
#srs-marsh-wet
#calculate percent contribution per prey group and reorder based on Per_Contri
srs.marsh.wet.abun <- data.frame(srs.marsh.wet) %>% mutate("Per_Contri" = srs.marsh.wet/sum(srs.marsh.wet)) %>% 
  rename("Av_Abun" = srs.marsh.wet) %>% arrange(desc(Per_Contri)) %>% mutate("C_sum" = cumsum(Per_Contri)) %>% 
  mutate("Habitat_Season_Slough" = "srs.marsh.wet") %>% rownames_to_column("Prey Group")
#filter for those that make up 75% 
srs.marsh.wet.abun75 <- srs.marsh.wet.abun %>% filter(C_sum < 0.76)
#srs-pond-dry
#calculate percent contribution per prey group and reorder based on Per_Contri
srs.pond.dry.abun <- data.frame(srs.pond.dry) %>% mutate("Per_Contri" = srs.pond.dry/sum(srs.pond.dry)) %>% 
  rename("Av_Abun" = srs.pond.dry) %>% arrange(desc(Per_Contri)) %>% mutate("C_sum" = cumsum(Per_Contri)) %>% 
  mutate("Habitat_Season_Slough" = "srs.pond.dry") %>% rownames_to_column("Prey Group")
#filter for those that make up 75% 
srs.pond.dry.abun75 <- srs.pond.dry.abun %>% filter(C_sum < 0.76)
#srs-np-dry-
#calculate percent contribution per prey group and reorder based on Per_Contri
srs.np.dry.abun <- data.frame(srs.np.dry) %>% mutate("Per_Contri" = srs.np.dry/sum(srs.np.dry)) %>% 
  rename("Av_Abun" = srs.np.dry) %>% arrange(desc(Per_Contri)) %>% mutate("C_sum" = cumsum(Per_Contri)) %>% 
  mutate("Habitat_Season_Slough" = "srs.np.dry") %>% rownames_to_column("Prey Group")
#filter for those that make up 75% 
srs.np.dry.abun75 <- srs.np.dry.abun %>% filter(C_sum < 0.76)
#srs-marsh-dry
#calculate percent contribution per prey group and reorder based on Per_Contri
srs.marsh.dry.abun <- data.frame(srs.marsh.dry) %>% mutate("Per_Contri" = srs.marsh.dry/sum(srs.marsh.dry)) %>% 
  rename("Av_Abun" = srs.marsh.dry) %>% arrange(desc(Per_Contri)) %>% mutate("C_sum" = cumsum(Per_Contri)) %>% 
  mutate("Habitat_Season_Slough" = "srs.marsh.dry") %>% rownames_to_column("Prey Group")
#filter for those that make up 75% 
srs.marsh.dry.abun75 <- srs.marsh.dry.abun %>% filter(C_sum < 0.76)
#-------TSL reformat SIMPER results----------------------------------------------------------------------------------------------------------
#tsl-pond-wet
#calculate percent contribution per prey group and reorder based on Per_Contri
tsl.pond.wet.abun <- data.frame(tsl.pond.wet) %>% mutate("Per_Contri" = tsl.pond.wet/sum(tsl.pond.wet)) %>% 
  rename("Av_Abun" = tsl.pond.wet) %>% arrange(desc(Per_Contri)) %>% mutate("C_sum" = cumsum(Per_Contri)) %>% 
  mutate("Habitat_Season_Slough" = "tsl.pond.wet") %>% rownames_to_column("Prey Group")
#filter for those that make up 75% 
tsl.pond.wet.abun75 <- tsl.pond.wet.abun %>% filter(C_sum < 0.76)
#tsl-np-wet-
#calculate percent contribution per prey group and reorder based on Per_Contri
tsl.np.wet.abun <- data.frame(tsl.np.wet) %>% mutate("Per_Contri" = tsl.np.wet/sum(tsl.np.wet)) %>% 
  rename("Av_Abun" = tsl.np.wet) %>% arrange(desc(Per_Contri)) %>% mutate("C_sum" = cumsum(Per_Contri)) %>% 
  mutate("Habitat_Season_Slough" = "tsl.np.wet") %>% rownames_to_column("Prey Group")
#filter for those that make up 75% 
tsl.np.wet.abun75 <- tsl.np.wet.abun %>% filter(C_sum < 0.76)
#tsl-marsh-wet
#calculate percent contribution per prey group and reorder based on Per_Contri
tsl.marsh.wet.abun <- data.frame(tsl.marsh.wet) %>% mutate("Per_Contri" = tsl.marsh.wet/sum(tsl.marsh.wet)) %>% 
  rename("Av_Abun" = tsl.marsh.wet) %>% arrange(desc(Per_Contri)) %>% mutate("C_sum" = cumsum(Per_Contri)) %>% 
  mutate("Habitat_Season_Slough" = "tsl.marsh.wet") %>% rownames_to_column("Prey Group")
#filter for those that make up 75% 
tsl.marsh.wet.abun75 <- tsl.marsh.wet.abun %>% filter(C_sum < 0.76)

#tsl-pond-dry
#calculate percent contribution per prey group and reorder based on Per_Contri
tsl.pond.dry.abun <- data.frame(tsl.pond.dry) %>% mutate("Per_Contri" = tsl.pond.dry/sum(tsl.pond.dry)) %>% 
  rename("Av_Abun" = tsl.pond.dry) %>% arrange(desc(Per_Contri)) %>% mutate("C_sum" = cumsum(Per_Contri)) %>% 
  mutate("Habitat_Season_Slough" = "tsl.pond.dry") %>% rownames_to_column("Prey Group")
#filter for those that make up 75% 
tsl.pond.dry.abun75 <- tsl.pond.dry.abun %>% filter(C_sum < 0.76)
#no tsl-np-dry
#tsl-marsh-dry
#calculate percent contribution per prey group and reorder based on Per_Contri
tsl.marsh.dry.abun <- data.frame(tsl.marsh.dry) %>% mutate("Per_Contri" = tsl.marsh.dry/sum(tsl.marsh.dry)) %>% 
  rename("Av_Abun" = tsl.marsh.dry) %>% arrange(desc(Per_Contri)) %>% mutate("C_sum" = cumsum(Per_Contri)) %>% 
  mutate("Habitat_Season_Slough" = "tsl.marsh.dry") %>% rownames_to_column("Prey Group")
#filter for those that make up 75% 
tsl.marsh.dry.abun75 <- tsl.marsh.dry.abun %>% filter(C_sum < 0.76)
#---SIMPER plot------------------------------------------------------------------------------------------------------------------------------
#recombine all simper results
post.simper.plot.data <- bind_rows(srs.pond.wet.abun, srs.np.wet.abun, srs.marsh.wet.abun, srs.pond.dry.abun, srs.np.dry.abun, 
                              srs.marsh.dry.abun,
                              tsl.pond.wet.abun, tsl.np.wet.abun, tsl.marsh.wet.abun, tsl.pond.dry.abun, tsl.marsh.dry.abun) %>% 
  separate(Habitat_Season_Slough, into  = c("Slough", "Habitat", "Season"), sep = "[.//]") %>% 
  mutate(Slough = recode_factor(Slough, srs = "SRS", tsl = "TSL"),
         Habitat = recode_factor(Habitat, marsh = "Marsh", np = "Near Pond", pond = "Pond"),
         Season = recode_factor(Season, wet = "Wet", dry = "Dry"))

(post.simper.plot <- ggplot(data = post.simper.plot.data, aes(x = `Prey Group`, y  = Av_Abun))+
    geom_bar(stat = "identity", alpha = 0.8)+
    labs(y = "Av. Abun", x = "Prey Group")+
    #facet_wrap(~ Habitat + Season + Slough) +
    facet_grid(Habitat ~ Slough + Season)+
    scale_y_continuous(expand = c(0,0), limits = c(0, 6.1))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          strip.text.x = element_text(size = 5),
          strip.text.y = element_text(size = 5),
          #panel.border = element_blank(), panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(), axis.line = element_line(color = 'black'),
          axis.title = element_text(size = 5), axis.text = element_text(size = 4))
          #legend.title = element_text(size = 4), legend.text = element_text(size = 4))
  )

#export as .tiff file - 
ggsave("post_simper_plot.tiff", plot = post.simper.plot, width = 160, height = 90, units = "mm", dpi = 600)

(post.simper.meiofauna.plot <- ggplot(data = post.simper.plot.data %>% 
                                        filter(`Prey Group` %in% c("Amphipoda", "Chironomid", "Cladocera",
                                                                   "Copepoda", "Hydrachnidia", "Ostracoda")),
                                      aes(x = `Prey Group`, y  = Av_Abun, fill = Season, color = Season))+
    geom_bar(stat = "identity", alpha = 0.8)+
    labs(y = "Av. Abun", x = "Prey Group")+
    #facet_wrap(~ Habitat + Season + Slough) +
    facet_grid(Habitat ~ Slough)+
    scale_y_continuous(expand = c(0,0), limits = c(0, 6.1))+
    scale_fill_viridis_d(begin = 0.2, end = 0.88)+
    scale_color_viridis_d(begin = 0.2, end = 0.88)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          strip.text.x = element_text(size = 5),
          strip.text.y = element_text(size = 5),
          #panel.border = element_blank(), panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(), axis.line = element_line(color = 'black'),
          axis.title = element_text(size = 5), axis.text = element_text(size = 4),
          axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1))
  #legend.title = element_text(size = 4), legend.text = element_text(size = 4)))
)
#export as .tiff file - 
ggsave("post_simper_meiofauna_plot.tiff", plot = post.simper.meiofauna.plot, width = 160, height = 90, units = "mm", dpi = 600)
#---------PERMANOVA--------------------------------------------------------------------------------------------------------------------------
post.permanova <- adonis(prey.no.rare.no.zero ~ habitat.season.slough, method = "horn", permutations = 1000)
#look at results
post.permanova$aov.tab

#let's try a different model structure - this should be a better approach
post.permonva2 <- adonis( #wow my spelling was spot on here... don't want to wait to re-run this to fix it
  prey.no.rare.no.zero ~ post.covariates.no.zero$size_class*post.covariates.no.zero$Season*post.covariates.no.zero$Slough*post.covariates.no.zero$Habitat,
  method = "horn", permutations = 1000)
#look at results
post.permonva2$aov.tab
write.xlsx(post.permonva2$aov.tab, "post_permanova.xlsx", rowNames = T)

#-----Pairwise PERMANOVAs---------------------------------------------------------------------------------------------------------------
#create combined factor
#add size_class to habitat.season.invasion
slough.habitat.season.size.class <- post.covariates.no.zero %>% select(size_class, Slough, Habitat, Season, ) %>% 
  unite(size_class:Season, col = "ont.slough.habitat.season", sep = "_")
slough.habitat.season.size.class <- as.factor(slough.habitat.season.size.class$ont.slough.habitat.season)

#pairwise permanova - took about 15 hours on my old laptop
slough.habitat.season.pairwise.permanova <- pairwise.perm.manova(vegdist(prey.no.rare.no.zero, method = "horn"), 
                                                                 slough.habitat.season.size.class,
                                                    progress = T, F = T, R2 = T)

post.pairwise.permanova <- slough.habitat.season.pairwise.permanova
#format output
#have to add one row and one column to prepare to add r2 to upper triangle
post.pairwise.p <- post.pairwise.permanova[["p.value"]] %>% t() %>% as.data.frame() %>%
  mutate(`A.cal2_SRS_Marsh_Wet` = NA) %>% 
  select(`A.cal2_SRS_Marsh_Wet`, everything()) %>% 
  t() %>% as.data.frame() %>% 
  mutate(`T.mar1_TSL_Pond_Wet` = NA) %>% 
  round(3)


post.pairwise.r2 <- post.pairwise.permanova[["R2.value"]] %>%  t() %>% as.data.frame() %>%
  mutate(`A.cal2_Marsh_Wet_Post-Invasion` = NA) %>% 
  select(`A.cal2_Marsh_Wet_Post-Invasion`, everything()) %>% 
  t() %>% as.data.frame() %>% 
  mutate(`T.mar1_Pond_Wet_Pre-Invasion` = NA) %>% 
  round(3) %>% t() %>% as.data.frame()

#lower triangle as p-values and upper triangle as R2
post.pairwise.p.r2 <- post.pairwise.p
post.pairwise.p.r2[upper.tri(post.pairwise.p.r2)] <- post.pairwise.r2[upper.tri(post.pairwise.r2)]
#bold p-values below 0.05
emphasize.strong.cells(which(lower.tri(post.pairwise.p < 0.05)))

#subset by certain comparisons (e.g., slough, habitat, season)

#slough shifts
post.pairwise.slough.wide <- post.pairwise.p %>% t() %>% as.data.frame() %>% 
  rownames_to_column("slough_habitat_season") %>% 
  filter(grepl("SRS", slough_habitat_season)) %>% 
  separate("slough_habitat_season", into = c("size_class", "slough", "habitat", "season"), sep = "_") %>% 
  select(-slough) %>% 
  unite("habitat_season", c(size_class, habitat, season), sep = "_") %>% 
  column_to_rownames("habitat_season") %>% 
  select(contains("TSL")) %>% 
  select_all(~gsub("TSL_", "", .)) %>% 
  select(matches(row.names(.))) %>% 
  rownames_to_column("SRS") %>% 
  filter(SRS %in% colnames(.)) %>% 
  column_to_rownames("SRS") 

post.pairwise.slough.long <- post.pairwise.slough.wide %>% as.matrix() %>% diag() %>% as.data.frame() %>% 
  mutate(Comparison = "SRS-TS") %>% 
  rownames_to_column("spp_habitat_season") %>% 
  separate(spp_habitat_season, into = c("Size Class", "Habitat", "Season"), sep = "_") %>% 
  rename(p = ".") %>% 
  select(`Size Class`, Season, Habitat, Comparison, p) %>% 
  arrange(`Size Class`, Season, Habitat, Comparison, p) %>% 
  filter(!is.na(p))

post.pairwise.slough.diet.shifts <- post.pairwise.slough.long %>% 
  filter(p<=0.05)

#proportion of slough comparisons that shifted
nrow(post.pairwise.slough.diet.shifts)/nrow(post.pairwise.slough.long)
#8.3%


#habitat shifts
#marsh-np shifts
marsh.np.pairwise.wide <- post.pairwise.p %>% t() %>% as.data.frame() %>% 
  rownames_to_column("slough_habitat_season") %>% 
  filter(grepl("Marsh", slough_habitat_season)) %>% 
  separate("slough_habitat_season", into = c("size_class", "slough", "habitat", "season"), sep = "_") %>% 
  select(-habitat) %>% 
  unite("slough_season", c(size_class, slough, season), sep = "_") %>% 
  column_to_rownames("slough_season") %>% 
  select(contains("NP")) %>% 
  select_all(~gsub("NP_", "", .)) %>% 
  select(matches(row.names(.))) %>% 
  rownames_to_column("Marsh") %>% 
  filter(Marsh %in% colnames(.)) %>% 
  column_to_rownames("Marsh")

marsh.np.pairwise.long <- marsh.np.pairwise.wide %>% as.matrix() %>% diag() %>% as.data.frame() %>% 
  mutate(Comparison = "Marsh-NP") %>% 
  rownames_to_column("spp_slough_season") %>% 
  separate(spp_slough_season, into = c("Size Class", "Slough", "Season"), sep = "_") %>% 
  rename(p = ".") %>% 
  select(`Size Class`, Slough, Season, Comparison, p) %>% 
  arrange(`Size Class`, Slough, Season, Comparison, p) %>% 
  filter(!is.na(p))

marsh.np.pairwise.diet.shifts <- marsh.np.pairwise.long %>% 
  filter(p<=0.05)

#proportion of marsh-np comparisons that shifted
nrow(marsh.np.pairwise.diet.shifts)/nrow(marsh.np.pairwise.long)
#0%

#marsh-pond shifts
marsh.pond.pairwise.wide <- post.pairwise.p %>% t() %>% as.data.frame() %>% 
  rownames_to_column("slough_habitat_season") %>% 
  filter(grepl("Marsh", slough_habitat_season)) %>% 
  separate("slough_habitat_season", into = c("size_class", "slough", "habitat", "season"), sep = "_") %>% 
  select(-habitat) %>% 
  unite("slough_season", c(size_class, slough, season), sep = "_") %>% 
  column_to_rownames("slough_season") %>% 
  select(contains("Pond")) %>% 
  select_all(~gsub("Pond_", "", .)) %>% 
  select(matches(row.names(.))) %>% 
  rownames_to_column("Marsh") %>% 
  filter(Marsh %in% colnames(.)) %>% 
  column_to_rownames("Marsh")

marsh.pond.pairwise.long <- marsh.pond.pairwise.wide %>% as.matrix() %>% diag() %>% as.data.frame() %>% 
  mutate(Comparison = "Marsh-Pond") %>% 
  rownames_to_column("spp_slough_season") %>% 
  separate(spp_slough_season, into = c("Size Class", "Slough", "Season"), sep = "_") %>% 
  rename(p = ".") %>% 
  select(`Size Class`, Slough, Season, Comparison, p) %>% 
  arrange(`Size Class`, Slough, Season, Comparison, p) %>% 
  filter(!is.na(p))

marsh.pond.pairwise.diet.shifts <- marsh.pond.pairwise.long %>% 
  filter(p<=0.05)

#proportion of marsh-pond comparisons that shifted
nrow(marsh.pond.pairwise.diet.shifts)/nrow(marsh.pond.pairwise.long)
#12.7%

#np-pond shifts
np.pond.pairwise.wide <- post.pairwise.p %>% t() %>% as.data.frame() %>% 
  rownames_to_column("slough_habitat_season") %>% 
  filter(grepl("NP", slough_habitat_season)) %>% 
  separate("slough_habitat_season", into = c("size_class", "slough", "habitat", "season"), sep = "_") %>% 
  select(-habitat) %>% 
  unite("slough_season", c(size_class, slough, season), sep = "_") %>% 
  column_to_rownames("slough_season") %>% 
  select(contains("Pond")) %>% 
  select_all(~gsub("Pond_", "", .)) %>% 
  select(matches(row.names(.))) %>% 
  rownames_to_column("NP") %>% 
  filter(NP %in% colnames(.)) %>% 
  column_to_rownames("NP")

np.pond.pairwise.long <- np.pond.pairwise.wide %>% as.matrix() %>% diag() %>% as.data.frame() %>% 
  mutate(Comparison = "NP-Pond") %>% 
  rownames_to_column("spp_slough_season") %>% 
  separate(spp_slough_season, into = c("Size Class", "Slough", "Season"), sep = "_") %>% 
  rename(p = ".") %>% 
  select(`Size Class`, Slough, Season, Comparison, p) %>% 
  arrange(`Size Class`, Slough, Season, Comparison, p) %>% 
  filter(!is.na(p))

np.pond.pairwise.diet.shifts <- np.pond.pairwise.long %>% 
  filter(p<=0.05)

nrow(np.pond.pairwise.diet.shifts)/nrow(np.pond.pairwise.long)
#19.0%

#overall proportion of habitat shifts
(nrow(marsh.pond.pairwise.diet.shifts)+nrow(marsh.np.pairwise.diet.shifts)+nrow(np.pond.pairwise.diet.shifts))/
  (nrow(marsh.pond.pairwise.long)+nrow(marsh.np.pairwise.long)+nrow(np.pond.pairwise.long))
#11%

#seasonal shifts in diet
season.pairwise.wide <- post.pairwise.p %>% #t() %>% as.data.frame() %>% 
  rownames_to_column("slough_habitat_season") %>% 
  filter(grepl("Wet", slough_habitat_season)) %>% 
  separate("slough_habitat_season", into = c("size_class", "slough", "habitat", "season"), sep = "_") %>% 
  select(-season) %>% 
  unite("slough_habitat", c(size_class, slough, habitat), sep = "_") %>% 
  column_to_rownames("slough_habitat") %>% 
  select(contains("Dry")) %>% 
  select_all(~gsub("_Dry", "", .)) %>% 
  select(matches(row.names(.))) %>% 
  rownames_to_column("Wet") %>% 
  filter(Wet %in% colnames(.)) %>% 
  column_to_rownames("Wet") 

season.pairwise.long <- season.pairwise.wide %>% as.matrix() %>% diag() %>% as.data.frame() %>% 
  mutate(Comparison = "Wet-Dry") %>% 
  rownames_to_column("spp_slough_habitat") %>% 
  separate(spp_slough_habitat, into = c("Size Class", "Slough", "Habitat"), sep = "_") %>% 
  rename(p = ".") %>% 
  select(`Size Class`, Slough, Habitat, Comparison, p) %>% 
  arrange(`Size Class`, Slough, Habitat, Comparison, p) %>% 
  filter(!is.na(p))

season.pairwise.diet.shifts <- season.pairwise.long %>% 
  filter(p<=0.05)

#proportion of season comparisons that shifted
nrow(season.pairwise.diet.shifts)/nrow(season.pairwise.long)
#20.6%

write.xlsx(list(
  "Slough" = post.pairwise.slough.diet.shifts,
  "Season" = season.pairwise.diet.shifts,
  "Habitat" = bind_rows(
    marsh.np.pairwise.diet.shifts, marsh.pond.pairwise.diet.shifts, np.pond.pairwise.diet.shifts
  )), file = "post_spatiotemporal_diet_pairwise_permanovas.xlsx"
)
#------Prey Proportions Mean and SD----------------------------------------------------------------------------------------------------------
#calculate mean proportion and sd of proportion of each prey group
#create matrix to write over in loop  
post.prey.pro <- prey.com

#convert counts to proportions and create new matrix
for (i in 1:nrow(post.prey.pro)){
  for (j in 1:ncol(post.prey.pro)){
    post.prey.pro[i,j] <- post.prey.pro[i,j]/row.sums[i]
  }
} 

#this creates some NaN values because of attempts to divide 0 by 0
#change these to 0s
#also some Inf values (not sure why 1/1 is giving this) but change these to 1
post.prey.pro <- post.prey.pro %>% 
  mutate(across(.cols = everything(), ~ ifelse(is.infinite(.x), 1, .x))) %>% 
  mutate(across(.cols = everything(), ~ ifelse(is.na(.x), 0, .x)))

#double-check by summing first row - should equal 1
sum(post.prey.pro[1,])

post.prey.mean.sd <- post.prey.pro %>% rownames_to_column(var = "ID") %>% 
  merge(., post.covariates, by = "ID") %>% select(-c(ID)) %>% 
  group_by(Slough, Season, Habitat, size_class) %>% 
  summarise(across(everything(), list(mean = mean, sd = sd)))

#export to Excel
write.xlsx(post.prey.mean.sd, "post_prey_mean_sd.xlsx")

post.prey.mean.species.slough <- post.prey.pro %>% rownames_to_column(var = "ID") %>% 
  merge(., post.covariates, by = "ID") %>% select(-c(ID, Year, size_class, Season, Habitat, Site)) %>% 
  group_by(Slough, Species) %>% 
  summarise(across(everything(), list(mean = mean))) 

post.prey.mean.species.slough.codes <- post.prey.mean.species.slough[,3:32]

#create codes for food-web diagram
post.prey.mean.species.slough.codes[post.prey.mean.species.slough.codes > 0 & post.prey.mean.species.slough.codes < 0.05] <- 1
post.prey.mean.species.slough.codes[post.prey.mean.species.slough.codes >= 0.05 & post.prey.mean.species.slough.codes < 0.3] <- 2
post.prey.mean.species.slough.codes[post.prey.mean.species.slough.codes >= 0.3 & post.prey.mean.species.slough.codes < 0.5] <- 3
post.prey.mean.species.slough.codes[post.prey.mean.species.slough.codes >= 0.5] <- 4
post.prey.mean.species.slough.codes <- bind_cols(post.prey.mean.species.slough[,1:2], post.prey.mean.species.slough.codes)

#export to Excel
write.xlsx(post.prey.mean.species.slough.codes, "post_prey_mean_species_slough_codes.xlsx")

#Group at higher levels
#Functional groups - detritus, algae, plants, zooplankton, macroinverts, fishes
prey.pro.functional <- post.prey.pro %>% rownames_to_column(var = "ID") %>% 
  merge(., post.covariates, by = "ID") %>% select(-c(ID)) %>% 
  mutate(Producers = Vasc.Plant+Green.Algae,
         Detritus = Detritus+Periphyton,
                `Herb. Inverts` = Collembola+Ephemeroptera+Mollusca,
                `Omni. Inverts` = Amphipoda+Trichoptera+Chironomid+Diptera+Coleoptera+Hymenoptera+Oligochaeta+Cladocera+Copepoda+Ostracoda+Misc+Misc.Invert,
                `Carn. Inverts` = Hemiptera+Odonata+Araneae+`Belo/Naucor`+Hydrachnidia+Argulus,
                Decapods = Procambarus+Palaemonetes,
                `Sm. Fishes` = Cyprinodontid,
                `Lg. Fishes` = Centrarchid+Misc.Fish
  ) %>% 
    #Zooplankton = Amphipoda + Cladocera + Collembola + Copepoda + Hydrachnidia + Ostracoda,
    #Macroinvertebrates = Araneae + Argulus + `Belo/Naucor` + Chironomid + Coleoptera + Collembola + Diptera +
    #                     Ephemeroptera + Hemiptera + Hymenoptera + Misc + Misc.Invert + Mollusca + Odonata + Oligochaeta +
    #                     Palaemonetes + Procambarus + Trichoptera,
    #Fishes = Centrarchid + Cyprinodontid + Misc.Fish) %>% 
  #select(Detritus, Green.Algae, Vasc.Plant, Periphyton, Zooplankton, Macroinvertebrates, Fishes) %>% 
  #rename("Algae" = Green.Algae, "Plants" = Vasc.Plant) %>% 
  #rownames_to_column(var = "ID") %>% 
  #merge(., post.covariates, by = "ID") %>% select(-c(ID)) %>% 
  select(size_class, Slough, Season, Habitat, Producers, Detritus, 
         `Herb. Inverts`, `Omni. Inverts`, `Carn. Inverts`, Decapods,
         `Sm. Fishes`, `Lg. Fishes`) %>% 
  group_by(Slough, Season, Habitat, size_class) %>% 
  summarise(across(where(is.numeric), list(mean = mean))) %>% 
  select_all(~str_replace(., "_mean", "")) %>% 
  mutate(across(where(is.numeric), round, digits = 2)) %>% 
  select(size_class, Slough, Season, Habitat, everything()) %>% 
  rename("Size Class" = size_class)

#export to Excel
write.xlsx(prey.pro.functional, "post_prey_pro_functional.xlsx")

#############################################################################################################################################
#########Pre vs Post Invasion Comparisons####################################################################################################
#############################################################################################################################################
#########load and format pre-invasion data###################################################################################################
#load historic, pre-invasion data
setwd("G:/Other Computers/My Laptop/Documents/FIU/Research/BL_food_web/Stomach Contents")

pre.stomachs.raw <- read.csv(file = "All_Specimen_by_Number_Correct.csv")

#change wd back to desired output directory
setwd("G:/Other Computers/My Laptop/Documents/FIU/Research/Current Food Web Study/Analyses/Stomach Contents/Diet")

#------------------create size_class column for these data------------------------------------------------------------------------------
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==1 & pre.stomachs.raw$Size==1]<-"G_hol1"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==1 & pre.stomachs.raw$Size==2]<-"G_hol2"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==2 & pre.stomachs.raw$Size==1]<-"P_lat1"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==2 & pre.stomachs.raw$Size==2]<-"P_lat2"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==3 & pre.stomachs.raw$Size==1]<-"C_var1"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==4 & pre.stomachs.raw$Size==1]<-"F_con1"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==4 & pre.stomachs.raw$Size==2]<-"F_con2"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==5 & pre.stomachs.raw$Size==1]<-"J_flo1"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==5 & pre.stomachs.raw$Size==2]<-"J_flo2"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==6 & pre.stomachs.raw$Size==1]<-"H_for1"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==6 & pre.stomachs.raw$Size==2]<-"H_for2"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==7 & pre.stomachs.raw$Size==1]<-"L_goo1"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==7 & pre.stomachs.raw$Size==2]<-"L_goo2"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==8 & pre.stomachs.raw$Size==1]<-"F_chr1"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==8 & pre.stomachs.raw$Size==2]<-"F_chr2"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==9 & pre.stomachs.raw$Size==1]<-"N_cry1"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==10 & pre.stomachs.raw$Size==1]<-"L_sic1"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==11 & pre.stomachs.raw$Size==1]<-"A_nat1"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==11 & pre.stomachs.raw$Size==2]<-"A_nat2"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==11 & pre.stomachs.raw$Size==3]<-"A_nat3"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==11 & pre.stomachs.raw$Size==4]<-"A_nat4"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==15 & pre.stomachs.raw$Size==1]<-"A_xen1"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==17 & pre.stomachs.raw$Size==1]<-"L_pla1"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==17 & pre.stomachs.raw$Size==2]<-"L_pla2"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==20 & pre.stomachs.raw$Size==1]<-"M_sal1"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==30 & pre.stomachs.raw$Size==1]<-"N_pet1"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==32 & pre.stomachs.raw$Size==1]<-"N_mac1"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==33 & pre.stomachs.raw$Size==1]<-"E_suc1"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==33 & pre.stomachs.raw$Size==2]<-"E_suc2"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==33 & pre.stomachs.raw$Size==3]<-"E_suc3"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==37 & pre.stomachs.raw$Size==1]<-"N_gyr1"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==38 & pre.stomachs.raw$Size==1]<-"L_gul1"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==38 & pre.stomachs.raw$Size==2]<-"L_gul2"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==38 & pre.stomachs.raw$Size==3]<-"L_gul3"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==38 & pre.stomachs.raw$Size==4]<-"L_gul4"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==39 & pre.stomachs.raw$Size==1]<-"L_mac1"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==40 & pre.stomachs.raw$Size==1]<-"L_mic1"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==40 & pre.stomachs.raw$Size==2]<-"L_mic2"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==41 & pre.stomachs.raw$Size==1]<-"L_pun1"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==41 & pre.stomachs.raw$Size==2]<-"L_pun2"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==41 & pre.stomachs.raw$Size==3]<-"L_pun3"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==42 & pre.stomachs.raw$Size==1]<-"L_mar1"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==46 & pre.stomachs.raw$Size==1]<-"E_eve1"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==47 & pre.stomachs.raw$Size==1]<-"E_glo1"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==52 & pre.stomachs.raw$Size==1]<-"C_bim1"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==55 & pre.stomachs.raw$Size==1]<-"A_cal2"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==56 & pre.stomachs.raw$Size==1]<-"E_fus"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==60 & pre.stomachs.raw$Size==1]<-"C_bat1"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==77 & pre.stomachs.raw$Size==1]<-"T_mar1"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==84 & pre.stomachs.raw$Size==1]<-"C_uro1"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==84 & pre.stomachs.raw$Size==2]<-"C_uro2"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==86 & pre.stomachs.raw$Size==1]<-"B_bel1"
pre.stomachs.raw$size_class[pre.stomachs.raw$Species==86 & pre.stomachs.raw$Size==2]<-"B_bel2"

pre.stomachs.no.size.class <- pre.stomachs.raw %>% filter(is.na(size_class)) 
#there is 1, species number 85 is missing from historic metadata and n = 1
#will exclude from analyses
pre.stomachs.size.classes <- pre.stomachs.raw %>% filter(Species != 85)

#-------------add habdate variable to historic data------------------------------------------------------------------------------------
pre.stomachs.size.classes$HabDate[pre.stomachs.size.classes$Date==1 & pre.stomachs.size.classes$Habitat...codes==1]<-"Pond_Wet"
pre.stomachs.size.classes$HabDate[pre.stomachs.size.classes$Date==2 & pre.stomachs.size.classes$Habitat...codes==1]<-"Pond_Dry"
pre.stomachs.size.classes$HabDate[pre.stomachs.size.classes$Date==1 & pre.stomachs.size.classes$Habitat...codes==2]<-"Spikerush_Wet"
pre.stomachs.size.classes$HabDate[pre.stomachs.size.classes$Date==2 & pre.stomachs.size.classes$Habitat...codes==2]<-"Spikerush_Dry"
pre.stomachs.size.classes$HabDate[pre.stomachs.size.classes$Date==1 & pre.stomachs.size.classes$Habitat...codes==3]<-"Sawgrass_Wet"
pre.stomachs.size.classes$HabDate[pre.stomachs.size.classes$Date==2 & pre.stomachs.size.classes$Habitat...codes==3]<-"Sawgrass_Dry"
pre.stomachs.size.classes$HabDate[pre.stomachs.size.classes$Date==1 & pre.stomachs.size.classes$Habitat...codes == 5 | 
                                    pre.stomachs.size.classes$Habitat...codes == 6 |
                                    pre.stomachs.size.classes$Habitat...codes == 7]<-"Pond_Wet"
pre.stomachs.size.classes$HabDate[pre.stomachs.size.classes$Date==2 & pre.stomachs.size.classes$Habitat...codes == 4]<-"Pond_Dry"

filter(pre.stomachs.size.classes, is.na(HabDate)) #should be 0 rows

#------Reformat to match post-invasion habitats and prey cols--------------------------------------------------------------------------------------------------------
#remove sawgrass data
pre.stomachs.no.sawgrass <- pre.stomachs.size.classes %>% filter(!grepl("Sawgrass", HabDate))

#combine prey columns to match the post-invasion prey columns 
pre.stomachs.post.prey <- pre.stomachs.no.sawgrass %>% 
  filter(Empty == 0) %>% #remove empty stomachs
  rename("Amphipoda" = H..azteca, "Araneae" = Aranae, "Argulus" = Argulus.sp, "Belo/Naucor" = Belo.Naucor, "Centrarchid" = Centrachidae, 
         "Cyprinodontid" = Cyprinodont, "Ephemeroptera" = Ephemoptera, "Hydrachnidia" = Mites, "Green.Algae" = Green.algae, "Misc" = Misc.,
          "Oligochaeta" = Oligocheata, "Vasc.Plant" = Vasc..Plants, #rename prey groups to match post-invasion
         "Length" = SL, "Wet_Weight" = Weight) %>% #rename groups and variables to match post-invasion
  mutate(Detritus = Det.Alg.mix+Vasc..Detritus, Hemiptera = Hemiptera+Homoptera, Misc.Fish = Misc..Fish+A..natalis, 
         Misc.Invert = Mysids+Misc..Insecta) %>% #combine groups where needed
  select(-c(Det.Alg.mix, Vasc..Detritus, Homoptera, Misc..Fish, A..natalis, Mysids, Misc..Insecta, #drop vestigial prey groups after mutate 
            Species, Size, Date, Habitat...codes, code, Habitat...Orig, TL, Specimen, Empty)) %>% 
  select(order(colnames(.))) %>% #alphabetize columns
  mutate(Invasion_Status = "Pre-Invasion") %>%
  select(size_class, Invasion_Status, HabDate, Length, Wet_Weight, everything()) %>%  #pull covariates out in front
  separate(HabDate, into = c("Habitat", "Season"), sep = "_")  #split HabDate into Habitat and Season
  
  
#reformat post-invasion data and bind w/ pre-invasion data
pre.post.stomachs <- stomach.counts %>% filter(Empty == 0, Slough == "SRS", Habitat != "NP") %>% #Remove empty stomachs, TSL, and NP data
  select(-c(Species, Year, Slough, Site, Sex, Empty)) %>% #drop variables we don't want to retain 
  column_to_rownames(var = "ID") %>% #move IDs to rownames, probably won't need these but preserving them here just in case
  mutate(Invasion_Status = "Post-Invasion") %>% 
  select(size_class, Invasion_Status, Habitat, Season, Length, Wet_Weight, everything()) %>% #reorder columns
  bind_rows(pre.stomachs.post.prey) #bind w/ pre-invasion data


write.csv(pre.post.stomachs, file = "pre_post_diet_matrix.csv")
#------pre/post NMDS format and run--------------------------------------------------------------------------------------
#isolate covariates
pre.post.covariates <- pre.post.stomachs %>% select(size_class:Wet_Weight)

#isolate prey groups
pre.post.prey.com <- pre.post.stomachs %>% select(Amphipoda:Vasc.Plant)

#historic data has some counts < 1 for trace amounts, but the post-invasion data does not do this
#changing 0 < values < 1 to == 1
pre.post.prey.com[pre.post.prey.com < 1 & pre.post.prey.com > 0] <- 1

#look at column sums to assess rare prey groups
sort(apply(pre.post.prey.com, 2, sum), decreasing = F) #some of these still aren't integers... 

#remove rare prey groups - not removing any to start

#drop zero rows 
#get row sums
pre.post.row.sums <- apply(pre.post.prey.com, 1, sum)
#create T/F vector for where rows > 0
pre.post.no.zero.rows <- pre.post.row.sums > 0
#select for non-zero rows
pre.post.prey.no.zero <- pre.post.prey.com[pre.post.no.zero.rows,]

#make covariates the match
pre.post.covariates.no.zero <- pre.post.covariates[pre.post.no.zero.rows,] %>% rownames_to_column(var = "ID")
#log transform the data
pre.post.prey.log.com <- log(pre.post.prey.no.zero+1)

#Run NMDS
set.seed(999)
pre.post.nmds <- metaMDS(pre.post.prey.log.com, k=2, distance = "horn", noshare = 0.1, autotransform = F, plot = T)
pre.post.nmds

#quick visulization
ordiplot(pre.post.nmds, type = "text", main = sprintf("k = 2; stress = %.3f", pre.post.nmds$stress))

#calculate stress values by hand
#reconstruct community distance as in the nmds, accounting for zeros
pre.post.comm.dist <- metaMDSredist(pre.post.nmds, distance = "horn") 
#calculate ordination distance
pre.post.ord.dist <- vegdist(pre.post.nmds$points[,1:2], method = "euclidean")
#extract residuals from a monotonic regression
pre.post.dhat <- isoreg(pre.post.comm.dist, pre.post.ord.dist)
#plot monotonic regression
plot(pre.post.dhat)
#calculate stress w/ residuals of this regression and ordination distances
pre.post.reg.stress <- sqrt(residuals(pre.post.dhat)^2/sum(pre.post.ord.dist^2))
pre.post.stress.label <- sprintf("Stress = %.2f", pre.post.reg.stress*100)
stressplot(pre.post.nmds)

#------pre/post-invasion NMDS format and plot------------------------------------------------------------------------------------------------------
pre.post.nmds.consumers <- as.data.frame(pre.post.nmds$points) %>% rownames_to_column(var = "ID") %>% 
  merge(., pre.post.covariates.no.zero, by = "ID") %>% 
  mutate(Habitat = case_when(Habitat == "Spikerush" ~ "Marsh", 
                             Habitat == "Pond" ~ "Pond",
                             T ~ Habitat)) %>% 
  arrange(desc(ID))
#write this out to Excel for future use in SIBER
write.xlsx(pre.post.nmds.consumers, "pre_post_nmds_consumers.xlsx")

#centroids per size class
pre.post.nmds.spp.centroid <- pre.post.nmds.consumers %>% group_by(size_class, Invasion_Status) %>% 
  summarise(MDS1 = mean(MDS1), MDS2 = mean(MDS2)) 

pre.post.nmds.prey <- as.data.frame(pre.post.nmds$species) %>% rownames_to_column(var = "Prey Group")

(pre.post.nmds.plot <- ggplot()+
    stat_ellipse(data = pre.post.nmds.consumers, aes(x=MDS1, y = MDS2, color = Invasion_Status))+
    stat_conf_ellipse(data = pre.post.nmds.consumers, aes(x= MDS1, y = MDS2, color = Invasion_Status))+
    geom_text_repel(data = pre.post.nmds.spp.centroid, aes(x = MDS1, y = MDS2, label = size_class, color = "blue"),
                    max.overlaps = 30, size = 2)+
    geom_text_repel(data = pre.post.nmds.prey, aes(x = MDS1, y = MDS2, label = `Prey Group`, color = "red"), size = 2, max.overlaps = 15)+
    #scale_fill_viridis_d(option = "D", direction = 1, labels = c("Consumers", "Post-Invasion", "Pre-Invasion", "Prey"))+
    scale_color_viridis_d(option = "D", direction = 1, labels = c("Consumers", "Post-Invasion", "Pre-Invasion", "Prey"))+
    coord_fixed()+
    theme_bw()+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(color = 'black'),
          axis.title = element_text(size = 8), axis.text = element_text(size = 5),
          legend.title = element_text(size = 5), legend.text = element_text(size = 5),
          #legend.position = "none"
          )
)

#save plot
ggsave(filename = "pre_post_nmds_consumers.tiff", plot = pre.post.nmds.plot, dpi = 600, width = 140, height = 140, units = "mm")

#----Pre/Post-Invasion SIMPER--------------------------------------------------------------------------------------------------------------
pre.post.habitat.season.invasion <- pre.post.covariates.no.zero %>%  mutate(Invasion_Status = gsub("-Invasion", "", Invasion_Status), 
                                                                            Habitat = case_when(Habitat %in% c("Spikerush", "Marsh") ~ "Marsh",
                                                                                                Habitat == "Pond" ~ "Pond")) %>% 
  unite("habdate_invasion", c(Habitat, Season, Invasion_Status), sep = "_") %>% 
  select(habdate_invasion)

pre.post.habitat.season.invasion <- as.factor(pre.post.habitat.season.invasion$habdate_invasion)
levels(pre.post.habitat.season.invasion)

#run SIMPER
pre.post.simper <- simper(pre.post.prey.no.zero, pre.post.habitat.season.invasion)  
#results summary
pre.post.sim.sum <- summary(pre.post.simper) %>% lapply(., as.data.frame)
#export to Excel
write.xlsx(pre.post.sim.sum, file = "pre_post_simper.xlsx", row.names = T)

#isolate results per habitat-season-invasion level
#Pre-Invasion
pond.wet.pre <- pre.post.simper$Pond_Wet_Pre_Marsh_Dry_Pre$ava
pond.dry.pre <- pre.post.simper$Pond_Dry_Pre_Marsh_Dry_Pre$ava
marsh.wet.pre <- pre.post.simper$Marsh_Wet_Pre_Marsh_Dry_Pre$ava
marsh.dry.pre <- pre.post.simper$Marsh_Dry_Post_Marsh_Dry_Pre$avb

#Post-Invasion
pond.wet.post <- pre.post.simper$Pond_Wet_Post_Marsh_Dry_Post$ava
pond.dry.post <- pre.post.simper$Pond_Dry_Post_Marsh_Dry_Post$ava
marsh.wet.post <- pre.post.simper$Marsh_Wet_Post_Marsh_Dry_Post$ava
marsh.dry.post <- pre.post.simper$Marsh_Dry_Post_Pond_Wet_Pre$ava

#--------Reformat Pre-Invasion SIMPER results----------------------------------------------------------------------------------------------
#pond-wet-pre
#calculate percent contribution per prey group and reorder based on Per_Contri
pond.wet.pre.abun <- data.frame(pond.wet.pre) %>% mutate("Per_Contri" = pond.wet.pre/sum(pond.wet.pre)) %>% 
  rename("Av_Abun" = pond.wet.pre) %>% arrange(desc(Per_Contri)) %>% mutate("C_sum" = cumsum(Per_Contri)) %>% 
  mutate("Habitat_Season_Invasion" = "pond.wet.pre") %>% rownames_to_column("Prey Group")
#filter for those that make up 75% 
pond.wet.pre.abun75 <- pond.wet.pre.abun %>% filter(C_sum < 0.76)

#pond-dry-pre
#calculate percent contribution per prey group and reorder based on Per_Contri
pond.dry.pre.abun <- data.frame(pond.dry.pre) %>% mutate("Per_Contri" = pond.dry.pre/sum(pond.dry.pre)) %>% 
  rename("Av_Abun" = pond.dry.pre) %>% arrange(desc(Per_Contri)) %>% mutate("C_sum" = cumsum(Per_Contri)) %>% 
  mutate("Habitat_Season_Invasion" = "pond.dry.pre") %>% rownames_to_column("Prey Group")
#filter for those that make up 75% 
pond.dry.pre.abun75 <- pond.dry.pre.abun %>% filter(C_sum < 0.76)

#marsh-wet-pre
#calculate percent contribution per prey group and reorder based on Per_Contri
marsh.wet.pre.abun <- data.frame(marsh.wet.pre) %>% mutate("Per_Contri" = marsh.wet.pre/sum(marsh.wet.pre)) %>% 
  rename("Av_Abun" = marsh.wet.pre) %>% arrange(desc(Per_Contri)) %>% mutate("C_sum" = cumsum(Per_Contri)) %>% 
  mutate("Habitat_Season_Invasion" = "marsh.wet.pre") %>% rownames_to_column("Prey Group")
#filter for those that make up 75% 
marsh.wet.pre.abun75 <- marsh.wet.pre.abun %>% filter(C_sum < 0.76)

#marsh-dry-pre
#calculate percent contribution per prey group and reorder based on Per_Contri
marsh.dry.pre.abun <- data.frame(marsh.dry.pre) %>% mutate("Per_Contri" = marsh.dry.pre/sum(marsh.dry.pre)) %>% 
  rename("Av_Abun" = marsh.dry.pre) %>% arrange(desc(Per_Contri)) %>% mutate("C_sum" = cumsum(Per_Contri)) %>% 
  mutate("Habitat_Season_Invasion" = "marsh.dry.pre") %>% rownames_to_column("Prey Group")
#filter for those that make up 75% 
marsh.dry.pre.abun75 <- marsh.dry.pre.abun %>% filter(C_sum < 0.76)

#--------Reformat Post-Invasion SIMPER results----------------------------------------------------------------------------------------------
#pond-wet-post
#calculate percent contribution per Prey group and reorder based on Per_Contri
pond.wet.post.abun <- data.frame(pond.wet.post) %>% mutate("Per_Contri" = pond.wet.post/sum(pond.wet.post)) %>% 
  rename("Av_Abun" = pond.wet.post) %>% arrange(desc(Per_Contri)) %>% mutate("C_sum" = cumsum(Per_Contri)) %>% 
  mutate("Habitat_Season_Invasion" = "pond.wet.post") %>% rownames_to_column("Prey Group")
#filter for those that make up 75% 
pond.wet.post.abun75 <- pond.wet.post.abun %>% filter(C_sum < 0.76)

#pond-dry-post
#calculate percent contribution per Prey group and reorder based on Per_Contri
pond.dry.post.abun <- data.frame(pond.dry.post) %>% mutate("Per_Contri" = pond.dry.post/sum(pond.dry.post)) %>% 
  rename("Av_Abun" = pond.dry.post) %>% arrange(desc(Per_Contri)) %>% mutate("C_sum" = cumsum(Per_Contri)) %>% 
  mutate("Habitat_Season_Invasion" = "pond.dry.post") %>% rownames_to_column("Prey Group")
#filter for those that make up 75% 
pond.dry.post.abun75 <- pond.dry.post.abun %>% filter(C_sum < 0.76)

#marsh-wet-post
#calculate percent contribution per Prey group and reorder based on Per_Contri
marsh.wet.post.abun <- data.frame(marsh.wet.post) %>% mutate("Per_Contri" = marsh.wet.post/sum(marsh.wet.post)) %>% 
  rename("Av_Abun" = marsh.wet.post) %>% arrange(desc(Per_Contri)) %>% mutate("C_sum" = cumsum(Per_Contri)) %>% 
  mutate("Habitat_Season_Invasion" = "marsh.wet.post") %>% rownames_to_column("Prey Group")
#filter for those that make up 75% 
marsh.wet.post.abun75 <- marsh.wet.post.abun %>% filter(C_sum < 0.76)

#marsh-dry-post
#calculate percent contribution per Prey group and reorder based on Per_Contri
marsh.dry.post.abun <- data.frame(marsh.dry.post) %>% mutate("Per_Contri" = marsh.dry.post/sum(marsh.dry.post)) %>% 
  rename("Av_Abun" = marsh.dry.post) %>% arrange(desc(Per_Contri)) %>% mutate("C_sum" = cumsum(Per_Contri)) %>% 
  mutate("Habitat_Season_Invasion" = "marsh.dry.post") %>% rownames_to_column("Prey Group")
#filter for those that make up 75% 
marsh.dry.post.abun75 <- marsh.dry.post.abun %>% filter(C_sum < 0.76)

#----Pre/Post-Invasion SIMPER plot------------------------------------------------------------------------------------------------------------------------------
#recombine all simper results
pre.post.simper.plot.data <- bind_rows(pond.wet.pre.abun, marsh.wet.pre.abun, pond.dry.pre.abun, marsh.dry.pre.abun,
                                   pond.wet.post.abun, marsh.wet.post.abun, pond.dry.post.abun, marsh.dry.post.abun) %>% 
  separate(Habitat_Season_Invasion, into  = c("Habitat", "Season", "Invasion Status"), sep = "[.//]") %>% 
  mutate(`Invasion Status` = recode_factor(`Invasion Status`, pre = "Pre-Invasion", post = "Post-Invasion"),
         Habitat = recode_factor(Habitat, marsh = "Marsh", pond = "Pond"),
         Season = recode_factor(Season, wet = "Wet", dry = "Dry"))

(pre.post.simper.plot <- ggplot(data = pre.post.simper.plot.data, aes(x = `Prey Group`, y  = Av_Abun))+
    geom_bar(stat = "identity", alpha = 0.8)+
    labs(y = "Av. Abun", x = "Prey Group")+
    #facet_wrap(~ `Habitat` + Season + Slough) +
    facet_grid(`Invasion Status` ~ Habitat + Season)+
    scale_y_continuous(expand = c(0,0), limits = c(0, 6.1))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          strip.text.x = element_text(size = 5),
          strip.text.y = element_text(size = 5),
          #panel.border = element_blank(), panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(), axis.line = element_line(color = 'black'),
          axis.title = element_text(size = 5), axis.text = element_text(size = 4))
  #legend.title = element_text(size = 4), legend.text = element_text(size = 4))
)

#export as .tiff file - 
ggsave("pre_post_simper_plot.tiff", plot = pre.post.simper.plot, width = 160, height = 90, units = "mm", dpi = 600)

#-------Pre/Post SIMPER Plot by Trophic Guilds------------------------------------------------------------------------------------------
pre.post.simper.plot.functional.data <- pre.post.simper.plot.data %>% 
  mutate(Guild = case_when(
    `Prey Group` %in% c("Vasc.Plant", "Green.Algae", "Periphyton") ~ "Producers",
    `Prey Group` == "Detritus" ~ "Detritus",
    `Prey Group` %in% c("Collembola", "Ephemeroptera", "Mollusca") ~ "Herb. Inverts",
    `Prey Group` %in% c("Amphipoda", "Chironomid", "Cladocera", "Coleoptera", "Copepoda", "Diptera", "Hymenoptera",
                       "Misc.Invert", "Misc", "Oligochaeta", "Ostracoda", "Trichoptera") ~ "Omni. Inverts",
    `Prey Group` %in% c("Araneae", "Argulus", "Belo/Naucor", "Hemiptera", "Hydrachnidia", "Odonata") ~ "Carn. Inverts",
    `Prey Group` %in% c("Palaemonetes", "Procambarus") ~ "Decapods",
    `Prey Group` %in% c("Cyprinodontid", "Misc.Fish") ~ "Small Fishes",
    `Prey Group` == "Centrarchid" ~ "Large Fishes"
  ),
  Order = case_when(
    Guild == "Detritus" ~ 1,
    Guild == "Producers" ~ 2,
    Guild == "Herb. Inverts" ~ 3,
    Guild == "Omni. Inverts" ~ 4,
    Guild == "Carn. Inverts" ~ 5,
    Guild == "Decapods" ~ 6,
    Guild == "Small Fishes" ~ 7,
    Guild == "Large Fishes" ~ 8
  ),
  Guild = fct_reorder(Guild, Order),
  `Prey Group` = fct_reorder(`Prey Group`, desc(Av_Abun)))


windowsFonts(Times = windowsFont("Times New Roman"))


(pre.post.simper.functional.plot <- ggplot(data = pre.post.simper.plot.functional.data, aes(x = Guild, y  = Av_Abun, fill = Guild))+
  geom_bar(stat = "identity", alpha = 0.8)+
  labs(y = "Av. Abun.", x = "Prey Guilds")+
  #facet_wrap(~ `Habitat` + Season + Slough) +
  facet_grid(Habitat + Season ~ `Invasion Status`, scales = "free")+
  scale_y_continuous(expand = c(0,0))+
  #scale_y_break(c(10,14))+
  scale_x_discrete(limits = c("Detritus", "Producers", "Herb. Inverts",
                              "Omni. Inverts", "Carn. Inverts", "Decapods",
                              "Small Fishes", "Large Fishes"))+
  scale_fill_viridis_d()+
  theme_bw()+
  theme(text = element_text(family = "Times", size = 8),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        strip.text.x = element_text(size = 8),
        strip.text.y = element_text(size = 8),
        #panel.border = element_blank(), panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(), axis.line = element_line(color = 'black'),
        axis.title = element_text(size = 8), axis.text = element_text(size = 8),
        axis.text.x.bottom = element_blank(), axis.ticks.x = element_blank())
#legend.title = element_text(size = 4), legend.text = element_text(size = 4))
)


#export as .tiff file - 
ggsave("pre_post_simper_functional_plot.tiff", plot = pre.post.simper.functional.plot, 
       width = 160, height = 90, units = "mm", dpi = 600)

#pre/post simper plot by species filled by guild
(pre.post.simper.spp.functional.plot <- ggplot(data = pre.post.simper.plot.functional.data, 
                                               aes(x = `Prey Group`, y  = Av_Abun, fill = Guild))+
    geom_bar(stat = "identity", alpha = 0.8)+
    labs(y = "Av. Abun.", x = NULL)+
    #facet_wrap(~ `Habitat` + Season + Slough) +
    facet_grid(Habitat + Season ~ `Invasion Status`, scales = "free")+
    scale_y_continuous(expand = c(0,0))+
    #scale_y_break(c(10,14))+
    scale_fill_viridis_d()+
    theme_bw()+
    theme(text = element_text(family = "Times", size = 8),
          axis.text.x = element_text(angle = 60, vjust = 1.05, hjust = 1),
          strip.text.x = element_text(size = 8),
          strip.text.y = element_text(size = 8),
          #panel.border = element_blank(), panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(), axis.line = element_line(color = 'black'),
          axis.title = element_text(size = 8), axis.text = element_text(size = 5),
          #axis.text.x.bottom = element_blank(), axis.ticks.x = element_blank())
  #legend.title = element_text(size = 4), legend.text = element_text(size = 4))
)
)
#export as .tiff file
ggsave("pre_post_simper_spp_functional_plot.tiff", plot = pre.post.simper.spp.functional.plot, 
       width = 160, height = 90, units = "mm", dpi = 600)
#---------Pre/Post PERMANOVA--------------------------------------------------------------------------------------------------------------------------
#this doesn't result in the full model running
#I don't think it likes that all size classes aren't present in each level of the other factors
pre.post.permanova2 <- adonis(pre.post.prey.no.zero ~ Invasion_Status*Season*Habitat*size_class, 
                              data = pre.post.covariates.no.zero, method = "horn", permutations = 1000)

pre.post.permanova2$aov.tab
write.xlsx(pre.post.permanova2$aov.tab, "pre_post_permanova.xlsx", rowNames = T)


#add size_class to habitat.season.invasion
pre.post.hab.season.invasion.size.class <- pre.post.covariates.no.zero %>% select(size_class, Habitat, Season, Invasion_Status) %>% 
  unite(size_class:Invasion_Status, col = "ont.hab.season.invasion", sep = "_")
pre.post.hab.season.invasion.size.class <- as.factor(pre.post.hab.season.invasion.size.class$ont.hab.season.invasion)
#pairwise permanova - took about 15 hours on my old laptop
#pre.post.pairwise.permanova <- pairwise.perm.manova(vegdist(pre.post.prey.no.zero, method = "horn"), pre.post.hab.season.invasion.size.class,
#                                                    progress = T, F = T, R2 = T)

#have to add one row and one column to prepare to add r2 to upper triangle
pre.post.pairwise.p <- pre.post.pairwise.permanova[["p.value"]] %>% t() %>% as.data.frame() %>%
  mutate(`A.cal2_Marsh_Wet_Post-Invasion` = NA) %>% 
  select(`A.cal2_Marsh_Wet_Post-Invasion`, everything()) %>% 
  t() %>% as.data.frame() %>% 
  mutate(`T.mar1_Pond_Wet_Pre-Invasion` = NA) %>% 
  round(3)

  
pre.post.pairwise.r2 <- pre.post.pairwise.permanova[["R2.value"]] %>%  t() %>% as.data.frame() %>%
  mutate(`A.cal2_Marsh_Wet_Post-Invasion` = NA) %>% 
  select(`A.cal2_Marsh_Wet_Post-Invasion`, everything()) %>% 
  t() %>% as.data.frame() %>% 
  mutate(`T.mar1_Pond_Wet_Pre-Invasion` = NA) %>% 
  round(3) %>% t() %>% as.data.frame()

#lower triangle as p-values and upper triangle as R2
pre.post.pairwise.p.r2 <- pre.post.pairwise.p
pre.post.pairwise.p.r2[upper.tri(pre.post.pairwise.p.r2)] <- pre.post.pairwise.r2[upper.tri(pre.post.pairwise.r2)]
#bold p-values below 0.05
emphasize.strong.cells(which(lower.tri(pre.post.pairwise.p < 0.05)))

#make a similar object per habitat per season where rows are pre and columns are post 
#only size classes that are in both pre- and post-invasion data
#post-invasion marsh-wet size class
#pre-invasion size classes are Genus_species, not Genus.species
#Need to change, for example, G_hol to G.hol
#pull out rownames to edit them
pre.post.pairwise.sc.format <- pre.post.pairwise.p %>% as.data.frame() %>% rownames_to_column("onthabdate_invasion") %>% 
  mutate(onthabdate_invasion = gsub("Spikerush", "Marsh", onthabdate_invasion))
#change the second character to "."
substring(pre.post.pairwise.sc.format$onthabdate_invasion, 2, 2) <- "."
#check to make sure it worked
levels(as.factor(pre.post.pairwise.sc.format$onthabdate_invasion))
#put that column back to rownames
pre.post.pairwise.sc.format <- pre.post.pairwise.sc.format %>% column_to_rownames(var = "onthabdate_invasion") 
#set colnames to match the rownames 
colnames(pre.post.pairwise.sc.format) <- rownames(pre.post.pairwise.sc.format)

#marsh-wet
#post-invasion marsh-wet size classes
post.marsh.wet.size.class <- pre.post.pairwise.sc.format %>% as.data.frame() %>% rownames_to_column("onthabdate_invasion") %>% 
  separate(onthabdate_invasion, into = c("size_class", "Habitat", "Season", "Invasion_Status"), sep = "_") %>% 
  filter(Habitat == "Marsh", Season == "Wet", Invasion_Status == "Post-Invasion") %>% 
  select(size_class)

#pre-invasion marsh-wet size classes
pre.marsh.wet.size.class <- pre.post.pairwise.sc.format %>% as.data.frame() %>% rownames_to_column("onthabdate_invasion") %>% 
  separate(onthabdate_invasion, into = c("size_class", "Habitat", "Season", "Invasion_Status"), sep = "_") %>% 
  filter(Habitat == "Marsh", Season == "Wet", Invasion_Status == "Pre-Invasion") %>% 
  select(size_class)

#share marsh-wet size classes in pre- and post-invasion
pre.post.marsh.wet.size.class <- pre.marsh.wet.size.class %>% filter(size_class %in% c(post.marsh.wet.size.class$size_class))

#create symmetric matrix of p.values and r2 where rows are pre-invasion and columns are post-invasion for the same size class
pairwise.marsh.wet.p <- pre.post.pairwise.sc.format %>% as.data.frame() %>% rownames_to_column("onthabdate_invasion") %>% 
  separate(onthabdate_invasion, into = c("size_class", "Habitat", "Season", "Invasion_Status"), sep = "_") %>% 
  filter(size_class %in% c(pre.post.marsh.wet.size.class$size_class), Habitat %in% c("Marsh", "Spikerush"), 
                           Season == "Wet", Invasion_Status == "Pre-Invasion") %>% 
  unite("onthabdate_invasion", size_class:Invasion_Status, sep = "_") %>% column_to_rownames(var = "onthabdate_invasion") %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column("onthabdate_invasion") %>% 
  separate(onthabdate_invasion, into = c("size_class", "Habitat", "Season", "Invasion_Status"), sep = "_") %>% 
  filter(size_class %in% c(pre.post.marsh.wet.size.class$size_class), Habitat %in% c("Marsh", "Spikerush"), 
                           Season == "Wet", Invasion_Status == "Post-Invasion") %>% 
  unite("onthabdate_invasion", size_class:Invasion_Status, sep = "_") %>% column_to_rownames(var = "onthabdate_invasion") %>% 
  t() %>% as.data.frame()
  
#marsh-Dry
#post-invasion marsh-Dry size classes
post.marsh.Dry.size.class <- pre.post.pairwise.sc.format %>% as.data.frame() %>% rownames_to_column("onthabdate_invasion") %>% 
  separate(onthabdate_invasion, into = c("size_class", "Habitat", "Season", "Invasion_Status"), sep = "_") %>% 
  filter(Habitat == "Marsh", Season == "Dry", Invasion_Status == "Post-Invasion") %>% 
  select(size_class)

#pre-invasion marsh-Dry size classes
pre.marsh.Dry.size.class <- pre.post.pairwise.sc.format %>% as.data.frame() %>% rownames_to_column("onthabdate_invasion") %>% 
  separate(onthabdate_invasion, into = c("size_class", "Habitat", "Season", "Invasion_Status"), sep = "_") %>% 
  filter(Habitat == "Marsh", Season == "Dry", Invasion_Status == "Pre-Invasion") %>% 
  select(size_class)

#share marsh-Dry size classes in pre- and post-invasion
pre.post.marsh.Dry.size.class <- pre.marsh.Dry.size.class %>% filter(size_class %in% c(post.marsh.Dry.size.class$size_class))

#create symmetric matrix of p.values and r2 where rows are pre-invasion and columns are post-invasion for the same size class
pairwise.marsh.Dry.p <- pre.post.pairwise.sc.format %>% as.data.frame() %>% rownames_to_column("onthabdate_invasion") %>% 
  separate(onthabdate_invasion, into = c("size_class", "Habitat", "Season", "Invasion_Status"), sep = "_") %>% 
  filter(size_class %in% c(pre.post.marsh.Dry.size.class$size_class), Habitat %in% c("Marsh", "Spikerush"), 
         Season == "Dry", Invasion_Status == "Pre-Invasion") %>% 
  unite("onthabdate_invasion", size_class:Invasion_Status, sep = "_") %>% column_to_rownames(var = "onthabdate_invasion") %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column("onthabdate_invasion") %>% 
  separate(onthabdate_invasion, into = c("size_class", "Habitat", "Season", "Invasion_Status"), sep = "_") %>% 
  filter(size_class %in% c(pre.post.marsh.Dry.size.class$size_class), Habitat %in% c("Marsh", "Spikerush"), 
         Season == "Dry", Invasion_Status == "Post-Invasion") %>% 
  unite("onthabdate_invasion", size_class:Invasion_Status, sep = "_") %>% column_to_rownames(var = "onthabdate_invasion") %>% 
  t() %>% as.data.frame()
  
#Pond-wet
#post-invasion Pond-wet size classes
post.Pond.wet.size.class <- pre.post.pairwise.sc.format %>% as.data.frame() %>% rownames_to_column("onthabdate_invasion") %>% 
  separate(onthabdate_invasion, into = c("size_class", "Habitat", "Season", "Invasion_Status"), sep = "_") %>% 
  filter(Habitat == "Pond", Season == "Wet", Invasion_Status == "Post-Invasion") %>% 
  select(size_class)

#pre-invasion Pond-wet size classes
pre.Pond.wet.size.class <- pre.post.pairwise.sc.format %>% as.data.frame() %>% rownames_to_column("onthabdate_invasion") %>% 
  separate(onthabdate_invasion, into = c("size_class", "Habitat", "Season", "Invasion_Status"), sep = "_") %>% 
  filter(Habitat == "Pond", Season == "Wet", Invasion_Status == "Pre-Invasion") %>% 
  select(size_class)

#share Pond-wet size classes in pre- and post-invasion
pre.post.Pond.wet.size.class <- pre.Pond.wet.size.class %>% filter(size_class %in% c(post.Pond.wet.size.class$size_class))

#create symmetric matrix of p.values and r2 where rows are pre-invasion and columns are post-invasion for the same size class
pairwise.Pond.wet.p <- pre.post.pairwise.sc.format %>% as.data.frame() %>% rownames_to_column("onthabdate_invasion") %>% 
  separate(onthabdate_invasion, into = c("size_class", "Habitat", "Season", "Invasion_Status"), sep = "_") %>% 
  filter(size_class %in% c(pre.post.Pond.wet.size.class$size_class), Habitat %in% c("Pond", "Spikerush"), 
         Season == "Wet", Invasion_Status == "Pre-Invasion") %>% 
  unite("onthabdate_invasion", size_class:Invasion_Status, sep = "_") %>% column_to_rownames(var = "onthabdate_invasion") %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column("onthabdate_invasion") %>% 
  separate(onthabdate_invasion, into = c("size_class", "Habitat", "Season", "Invasion_Status"), sep = "_") %>% 
  filter(size_class %in% c(pre.post.Pond.wet.size.class$size_class), Habitat %in% c("Pond", "Spikerush"), 
         Season == "Wet", Invasion_Status == "Post-Invasion") %>% 
  unite("onthabdate_invasion", size_class:Invasion_Status, sep = "_") %>% column_to_rownames(var = "onthabdate_invasion") %>% 
  t() %>% as.data.frame()

#Pond-Dry
#post-invasion Pond-Dry size classes
post.Pond.Dry.size.class <- pre.post.pairwise.sc.format %>% as.data.frame() %>% rownames_to_column("onthabdate_invasion") %>% 
  separate(onthabdate_invasion, into = c("size_class", "Habitat", "Season", "Invasion_Status"), sep = "_") %>% 
  filter(Habitat == "Pond", Season == "Dry", Invasion_Status == "Post-Invasion") %>% 
  select(size_class)

#pre-invasion Pond-Dry size classes
pre.Pond.Dry.size.class <- pre.post.pairwise.sc.format %>% as.data.frame() %>% rownames_to_column("onthabdate_invasion") %>% 
  separate(onthabdate_invasion, into = c("size_class", "Habitat", "Season", "Invasion_Status"), sep = "_") %>% 
  filter(Habitat == "Pond", Season == "Dry", Invasion_Status == "Pre-Invasion") %>% 
  select(size_class)

#share Pond-Dry size classes in pre- and post-invasion
pre.post.Pond.Dry.size.class <- pre.Pond.Dry.size.class %>% filter(size_class %in% c(post.Pond.Dry.size.class$size_class))

#create symmetric matrix of p.values and r2 where rows are pre-invasion and columns are post-invasion for the same size class
pairwise.Pond.Dry.p <- pre.post.pairwise.sc.format %>% as.data.frame() %>% rownames_to_column("onthabdate_invasion") %>% 
  separate(onthabdate_invasion, into = c("size_class", "Habitat", "Season", "Invasion_Status"), sep = "_") %>% 
  filter(size_class %in% c(pre.post.Pond.Dry.size.class$size_class), Habitat %in% c("Pond", "Spikerush"), 
         Season == "Dry", Invasion_Status == "Pre-Invasion") %>% 
  unite("onthabdate_invasion", size_class:Invasion_Status, sep = "_") %>% column_to_rownames(var = "onthabdate_invasion") %>% 
  t() %>% as.data.frame() %>% 
  rownames_to_column("onthabdate_invasion") %>% 
  separate(onthabdate_invasion, into = c("size_class", "Habitat", "Season", "Invasion_Status"), sep = "_") %>% 
  filter(size_class %in% c(pre.post.Pond.Dry.size.class$size_class), Habitat %in% c("Pond", "Spikerush"), 
         Season == "Dry", Invasion_Status == "Post-Invasion") %>% 
  unite("onthabdate_invasion", size_class:Invasion_Status, sep = "_") %>% column_to_rownames(var = "onthabdate_invasion") %>% 
  t() %>% as.data.frame()

write.xlsx(
  list("Marsh_Wet" = pairwise.marsh.wet.p, "Marsh_Dry" = pairwise.marsh.Dry.p, 
       "Pond_Wet" = pairwise.Pond.Dry.p, "Pond_Dry" = pairwise.Pond.Dry.p,
       "All" = pre.post.pairwise.p.r2
  ), "pre_post_pairwise_permanova.xlsx", rowNames = T
)

#within size class pre/post comparisons only
#marsh wet
pairwise.marsh.wet.long <- pairwise.marsh.wet.p %>% as.matrix() %>% diag() %>% as.data.frame() %>% 
  mutate(Species_Habitat_Season_Invasion = row.names(pairwise.marsh.wet.p)) %>% 
  separate(Species_Habitat_Season_Invasion, into = c("Size Class", "Habitat", "Season", "Invasion"), sep = "_") %>% 
  rename("p" = ".") %>% select(`Size Class`, Habitat, Season, p, -Invasion) %>% 
  filter(!is.na(p))

#marsh dry
pairwise.marsh.dry.long <- pairwise.marsh.Dry.p %>% as.matrix() %>% diag() %>% as.data.frame() %>% 
  mutate(Species_Habitat_Season_Invasion = row.names(pairwise.marsh.Dry.p)) %>% 
  separate(Species_Habitat_Season_Invasion, into = c("Size Class", "Habitat", "Season", "Invasion"), sep = "_") %>% 
  rename("p" = ".") %>% select(`Size Class`, Habitat, Season, p, -Invasion) %>% 
  filter(!is.na(p))

#pond wet
pairwise.pond.wet.long <- pairwise.Pond.wet.p %>% as.matrix() %>% diag() %>% as.data.frame() %>% 
  mutate(Species_Habitat_Season_Invasion = row.names(pairwise.Pond.wet.p)) %>% 
  separate(Species_Habitat_Season_Invasion, into = c("Size Class", "Habitat", "Season", "Invasion"), sep = "_") %>% 
  rename("p" = ".") %>% select(`Size Class`, Habitat, Season, p, -Invasion) %>% 
  filter(!is.na(p))

#pond dry
pairwise.pond.dry.long <- pairwise.Pond.Dry.p %>% as.matrix() %>% diag() %>% as.data.frame() %>% 
  mutate(Species_Habitat_Season_Invasion = row.names(pairwise.Pond.Dry.p)) %>% 
  separate(Species_Habitat_Season_Invasion, into = c("Size Class", "Habitat", "Season", "Invasion"), sep = "_") %>% 
  rename("p" = ".") %>% select(`Size Class`, Habitat, Season, p, -Invasion) %>% 
  filter(!is.na(p))

#combine
pre.post.size.class.p <- bind_rows(pairwise.marsh.wet.long, pairwise.marsh.dry.long,
                                   pairwise.pond.dry.long, pairwise.pond.wet.long)

write.xlsx(pre.post.size.class.p, "pre_post_size_class_permanova_comparisons.xlsx")

#proportion of size classes with statistically different diets pre- vs post-invasion
sum(pre.post.size.class.p$p<=0.05)/nrow(pre.post.size.class.p)

#------Pre/Post-Invasion Prey Proportions Mean and SD----------------------------------------------------------------------------------------------------------
#calculate mean proportion and sd of proportion of each prey group
#create matrix to write over in loop  
pre.post.prey.pro <- pre.post.prey.com

#convert counts to proportions and create new matrix
for (i in 1:nrow(pre.post.prey.pro)){
  for (j in 1:ncol(pre.post.prey.pro)){
    pre.post.prey.pro[i,j] <- pre.post.prey.pro[i,j]/pre.post.row.sums[i]
  }
}

#this creates some NaN values because of attempts to divide 0 by 0
#change these to 0s
pre.post.prey.pro[is.na(pre.post.prey.pro)] <- 0
pre.post.prey.pro[is.infinite(pre.post.prey.pro)] <- 0
#double-check by summing first row - should equal 1
sum(pre.post.prey.pro[1,])

pre.post.covariates.id <- pre.post.covariates %>% rownames_to_column(var = "ID")

pre.post.prey.mean.sd <- pre.post.prey.pro %>% rownames_to_column(var = "ID") %>% 
  merge(., pre.post.covariates.id, by = "ID") %>% select(-c(ID)) %>% 
  mutate(Habitat = gsub("Spikerush", "Marsh", Habitat)) %>% 
  mutate(size_class = gsub("_", ".", size_class)) %>% 
  group_by(Invasion_Status, Season, Habitat, size_class) %>% 
  summarise(across(everything(), list(mean = mean, sd = sd))) %>% 
  arrange(size_class, Habitat, Season, Invasion_Status) %>% 
  select(size_class, Invasion_Status, Habitat, Season, contains("Length"), contains("Weight"), everything()) %>% 
  mutate(across(is.numeric, round, digits = 2)) %>% 
  unite("Amphipoda", c(contains("Amphipoda")), sep = " \u00B1 ") %>% 
  unite("Araneae", c(contains("Araneae")), sep = " \u00B1 ") %>% 
  unite("Argulus", c(contains("Argulus")), sep = " \u00B1 ") %>% 
  unite("Belo/Naucor", c(contains("Belo/Naucor")), sep = " \u00B1 ") %>% 
  unite("Centrarchid", c(contains("Centrarchid")), sep = " \u00B1 ") %>% 
  unite("Chironomid", c(contains("Chironomid")), sep = " \u00B1 ") %>% 
  unite("Cladocera", c(contains("Cladocera")), sep = " \u00B1 ") %>% 
  unite("Coleoptera", c(contains("Coleoptera")), sep = " \u00B1 ") %>% 
  unite("Collembola", c(contains("Collembola")), sep = " \u00B1 ") %>% 
  unite("Copepoda", c(contains("Copepoda")), sep = " \u00B1 ") %>% 
  unite("Cyprinodontid", c(contains("Cyprinodontid")), sep = " \u00B1 ") %>% 
  unite("Detritus", c(contains("Detritus")), sep = " \u00B1 ") %>% 
  unite("Diptera", c(contains("Diptera")), sep = " \u00B1 ") %>% 
  unite("Ephemeroptera", c(contains("Ephemeroptera")), sep = " \u00B1 ") %>% 
  unite("Green.Algae", c(contains("Green.Algae")), sep = " \u00B1 ") %>% 
  unite("Hemiptera", c(contains("Hemiptera")), sep = " \u00B1 ") %>% 
  unite("Hydrachnidia", c(contains("Hydrachnidia")), sep = " \u00B1 ") %>% 
  unite("Hymenoptera", c(contains("Hymenoptera")), sep = " \u00B1 ") %>% 
  unite("Misc", c(contains("Misc_")), sep = " \u00B1 ") %>% 
  unite("Misc.Fish", c(contains("Misc.Fish")), sep = " \u00B1 ") %>% 
  unite("Misc.Invert", c(contains("Misc.Invert")), sep = " \u00B1 ") %>% 
  unite("Mollusca", c(contains("Mollusca")), sep = " \u00B1 ") %>% 
  unite("Odonata", c(contains("Odonata")), sep = " \u00B1 ") %>% 
  unite("Oligochaeta", c(contains("Oligochaeta")), sep = " \u00B1 ") %>% 
  unite("Ostracoda", c(contains("Ostracoda")), sep = " \u00B1 ") %>% 
  unite("Palaemonetes", c(contains("Palaemonetes")), sep = " \u00B1 ") %>% 
  unite("Periphyton", c(contains("Periphyton")), sep = " \u00B1 ") %>% 
  unite("Procambarus", c(contains("Procambarus")), sep = " \u00B1 ") %>% 
  unite("Trichoptera", c(contains("Trichoptera")), sep = " \u00B1 ") %>% 
  unite("Vasc.Plant", c(contains("Vasc.Plant")), sep = " \u00B1 ") %>% 
  unite("Length", c(contains("Length")), sep = " \u00B1 ") %>% 
  unite("Wet_Weight", c(contains("Weight")), sep = " \u00B1 ")  
  
pre.post.prey.n <- pre.post.prey.pro %>% rownames_to_column(var = "ID") %>% 
  merge(., pre.post.covariates.id, by = "ID") %>% select(-c(ID)) %>% 
  mutate(Habitat = gsub("Spikerush", "Marsh", Habitat)) %>% 
  mutate(size_class = gsub("_", ".", size_class)) %>% 
  group_by(Invasion_Status, Season, Habitat, size_class) %>% 
  tally() %>% 
  arrange(size_class, Habitat, Season, Invasion_Status) 

pre.post.prey.mean.sd.n <- bind_cols(pre.post.prey.mean.sd, pre.post.prey.n$n) %>% 
  rename(N = `...37`) %>% select(size_class, Invasion_Status, Habitat, Season, N, everything())

#export to Excel
write.xlsx(pre.post.prey.mean.sd.n, "pre_post_prey_mean_sd.xlsx")

#----------------group by functional guilds---------------------------------------------------------------------------------------------------
pre.post.guild.prey.pro <- pre.post.prey.pro %>% rownames_to_column(var = "ID") %>% 
  merge(., pre.post.covariates.id, by = "ID") %>% select(-c(ID)) %>% 
  mutate(Habitat = gsub("Spikerush", "Marsh", Habitat),
         size_class = gsub("_", ".", size_class),
         Producers = Vasc.Plant+Green.Algae,
         Detritus = Detritus+Periphyton,
         `Herb. Inverts` = Collembola+Ephemeroptera+Mollusca,
         `Omni. Inverts` = Amphipoda+Trichoptera+Chironomid+Diptera+Coleoptera+Hymenoptera+Oligochaeta+Cladocera+Copepoda+Ostracoda+Misc+Misc.Invert,
         `Carn. Inverts` = Hemiptera+Odonata+Araneae+`Belo/Naucor`+Hydrachnidia+Argulus,
         Decapods = Procambarus+Palaemonetes,
         `Sm. Fishes` = Cyprinodontid,
         `Lg. Fishes` = Centrarchid+Misc.Fish
         ) %>% 
  select(size_class, Invasion_Status, Habitat, Season, Producers, Detritus, `Herb. Inverts`, `Omni. Inverts`, 
         `Carn. Inverts`, Decapods, `Sm. Fishes`, `Lg. Fishes`) %>% 
  group_by(size_class, Invasion_Status, Habitat, Season) %>% 
  summarise(across(everything(), list(mean = mean))) %>% 
  set_names(~str_replace_all(., "_mean", ""))

write.xlsx(pre.post.guild.prey.pro, "pre_post_guild_prey_pro.xlsx")

pre.post.guild.long <- pre.post.guild.prey.pro %>% 
  pivot_longer(cols = c(Producers, Detritus, `Herb. Inverts`, `Omni. Inverts`, `Carn. Inverts`, Decapods, 
                        `Sm. Fishes`, `Lg. Fishes`), names_to = "Prey Type", values_to = "Proportion") %>% 
  unite(`Invasion Status - Habitat - Season`, c(Invasion_Status, Habitat, Season), sep = " - ")

pre.post.guild.pairwise <- pre.post.guild.prey.pro %>% 
  group_by(size_class) %>% filter(Invasion_Status)

pre.post.guild.long1 <- pre.post.guild.long[1:1120, ]
pre.post.guild.long2 <- pre.post.guild.long[1121:2240,]

#let's try to plot this now as facet_wrap or facet_grid pie charts
(pre.post.guild.pie.chart <- ggplot(pre.post.guild.long, aes(x = factor(1), y = Proportion, fill = `Prey Type`), width = 100)+
    geom_bar(stat = "identity")+
    geom_text(aes(label = round(Proportion, 2)*100), position = position_stack(vjust = 0.5), size = 0.1)+
    coord_polar(theta = "y")+
    facet_grid(cols = vars(size_class), rows = vars(`Invasion Status - Habitat - Season`), labeller = "label_parsed")+
    labs(x = "Diet Proportions", y = NULL)+
    theme_void()+
    theme(axis.title = element_text(size = 5), axis.title.x = element_text(angle = 90), axis.title.y = element_text(angle = 90),
          axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
          legend.key.size = unit(0.05, "mm"), legend.text = element_text(size = 5), legend.title = element_text(size = 5),
          strip.text.x = element_text(size = 5, angle = 90), strip.text.y = element_text(size = 5),
          legend.position = "bottom")
  
)

ggsave("pre_post_guild_pie_chart.tiff", plot = pre.post.guild.pie.chart, height = 174, width = 234, units = "mm", dpi = 600)


(pre.post.guild.pie.first.half <- ggplot(pre.post.guild.long1, aes(x = factor(1), y = Proportion, fill = `Prey Type`), width = 100)+
    geom_bar(stat = "identity")+
    geom_text(aes(label = round(Proportion, 2)*100), position = position_stack(vjust = 0.5), size = 0.1)+
    scale_y_continuous(expand = c(0,0))+
    scale_x_discrete(expand = c(0,0))+
    coord_polar(theta = "y")+
    facet_grid(cols = vars(size_class), rows = vars(`Invasion Status - Habitat - Season`), labeller = "label_parsed")+
    labs(x = "Diet Proportions", y = NULL)+
    theme_void()+
    theme(axis.title = element_text(size = 5), axis.title.x = element_text(angle = 90), axis.title.y = element_text(angle = 90),
          axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
          legend.key.size = unit(0.05, "mm"), legend.text = element_text(size = 5), legend.title = element_text(size = 5),
          strip.text.x = element_text(size = 5, angle = 90), strip.text.y = element_text(size = 5),
          legend.position = "bottom")
  
)

(pre.post.guild.pie.second.half <- ggplot(pre.post.guild.long2, aes(x = factor(1), y = Proportion, fill = `Prey Type`), width = 100)+
    geom_bar(stat = "identity")+
    geom_text(aes(label = round(Proportion, 2)*100), position = position_stack(vjust = 0.5), size = 0.1)+
    scale_y_continuous(expand = c(0,0))+
    scale_x_discrete(expand = c(0,0))+
    coord_polar(theta = "y")+
    facet_grid(cols = vars(size_class), rows = vars(`Invasion Status - Habitat - Season`), labeller = "label_parsed")+
    labs(x = "Diet Proportions", y = NULL)+
    theme_void()+
    theme(axis.title = element_text(size = 5), axis.title.x = element_text(angle = 90), axis.title.y = element_text(angle = 90),
          axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
          legend.key.size = unit(0.05, "mm"), legend.text = element_text(size = 5), legend.title = element_text(size = 5),
          strip.text.x = element_text(size = 5, angle = 90), strip.text.y = element_text(size = 5),
          legend.position = "bottom")
  
)

pre.post.guild.pie.stacked <- ggarrange(pre.post.guild.pie.first.half, pre.post.guild.pie.second.half, ncol = 1, nrow = 2, legend = "bottom",
                                        common.legend = T)

ggsave("pre_post_guild_pie_stacked.tiff", plot = pre.post.guild.pie.stacked, width = 174, height = 234, units = "mm", dpi = 600)


pre.post.guild.long.q1 <- pre.post.guild.long[1:(nrow(pre.post.guild.long)/4), ]
pre.post.guild.long.q2 <- pre.post.guild.long[((nrow(pre.post.guild.long)/4)+1):(nrow(pre.post.guild.long)/2), ]
pre.post.guild.long.q3 <- pre.post.guild.long[1121:1680, ]
pre.post.guild.long.q4 <- pre.post.guild.long[1681:2240, ]

(pre.post.guild.pie.q1 <- ggplot(pre.post.guild.long.q1, aes(x = factor(1), y = Proportion, fill = `Prey Type`), width = 100)+
    geom_bar(stat = "identity")+
    geom_text(aes(label = round(Proportion, 2)*100), position = position_stack(vjust = 0.5), size = 0.1)+
    scale_y_continuous(expand = c(0,0))+
    scale_x_discrete(expand = c(0,0))+
    coord_polar(theta = "y")+
    facet_grid(cols = vars(size_class), rows = vars(`Invasion Status - Habitat - Season`), labeller = "label_parsed")+
    labs(x = "Diet Proportions", y = NULL)+
    theme_void()+
    theme(axis.title = element_text(size = 5), axis.title.x = element_text(angle = 90), axis.title.y = element_text(angle = 90),
          axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
          legend.key.size = unit(0.05, "mm"), legend.text = element_text(size = 5), legend.title = element_text(size = 5),
          strip.text.x = element_text(size = 5, angle = 90), strip.text.y = element_text(size = 5),
          legend.position = "bottom")
  
)

(pre.post.guild.pie.q2 <- ggplot(pre.post.guild.long.q2, aes(x = factor(1), y = Proportion, fill = `Prey Type`), width = 100)+
    geom_bar(stat = "identity")+
    geom_text(aes(label = round(Proportion, 2)*100), position = position_stack(vjust = 0.5), size = 0.1)+
    scale_y_continuous(expand = c(0,0))+
    scale_x_discrete(expand = c(0,0))+
    coord_polar(theta = "y")+
    facet_grid(cols = vars(size_class), rows = vars(`Invasion Status - Habitat - Season`), labeller = "label_parsed")+
    labs(x = "Diet Proportions", y = NULL)+
    theme_void()+
    theme(axis.title = element_text(size = 5), axis.title.x = element_text(angle = 90), axis.title.y = element_text(angle = 90),
          axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
          legend.key.size = unit(0.05, "mm"), legend.text = element_text(size = 5), legend.title = element_text(size = 5),
          strip.text.x = element_text(size = 5, angle = 90), strip.text.y = element_text(size = 5),
          legend.position = "bottom")
  
)

(pre.post.guild.pie.q3 <- ggplot(pre.post.guild.long.q3, aes(x = factor(1), y = Proportion, fill = `Prey Type`), width = 100)+
    geom_bar(stat = "identity")+
    geom_text(aes(label = round(Proportion, 2)*100), position = position_stack(vjust = 0.5), size = 0.1)+
    scale_y_continuous(expand = c(0,0))+
    scale_x_discrete(expand = c(0,0))+
    coord_polar(theta = "y")+
    facet_grid(cols = vars(size_class), rows = vars(`Invasion Status - Habitat - Season`), labeller = "label_parsed")+
    labs(x = "Diet Proportions", y = NULL)+
    theme_void()+
    theme(axis.title = element_text(size = 5), axis.title.x = element_text(angle = 90), axis.title.y = element_text(angle = 90),
          axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
          legend.key.size = unit(0.05, "mm"), legend.text = element_text(size = 5), legend.title = element_text(size = 5),
          strip.text.x = element_text(size = 5, angle = 90), strip.text.y = element_text(size = 5),
          legend.position = "bottom")
  
)

(pre.post.guild.pie.q4 <- ggplot(pre.post.guild.long.q4, aes(x = factor(1), y = Proportion, fill = `Prey Type`), width = 100)+
    geom_bar(stat = "identity")+
    geom_text(aes(label = round(Proportion, 2)*100), position = position_stack(vjust = 0.5), size = 0.1)+
    scale_y_continuous(expand = c(0,0))+
    scale_x_discrete(expand = c(0,0))+
    coord_polar(theta = "y")+
    facet_grid(cols = vars(size_class), rows = vars(`Invasion Status - Habitat - Season`), labeller = "label_parsed")+
    labs(x = "Diet Proportions", y = NULL)+
    theme_void()+
    theme(axis.title = element_text(size = 5), axis.title.x = element_text(angle = 90), axis.title.y = element_text(angle = 90),
          axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
          legend.key.size = unit(0.05, "mm"), legend.text = element_text(size = 5), legend.title = element_text(size = 5),
          strip.text.x = element_text(size = 5, angle = 90), strip.text.y = element_text(size = 5),
          legend.position = "bottom")
  
)

pre.post.guild.pie.q.stacked <- ggarrange(pre.post.guild.pie.q1, pre.post.guild.pie.q2,
                                          pre.post.guild.pie.q3, pre.post.guild.pie.q4, 
                                          ncol = 1, nrow = 4, legend = "bottom",
                                          common.legend = T)

ggsave("pre_post_guild_pie_q4.tiff", plot = pre.post.guild.pie.q.stacked, width = 174, height = 234, units = "mm", dpi = 600)


#library(ggstatsplot)
#might explore this option in the future

#------Table for Food-Web Review Paper----------------------------------------------------------------------
species.prey.slough <- stomach.counts %>% 
  filter(Slough == "SRS", Habitat == "Marsh") %>% 
  group_by(Species) %>% 
  summarise(across(c(Amphipoda:Vasc.Plant), sum)) %>% 
  column_to_rownames(var = "Species")

#create matrix to write over
spp.prey.pro <- species.prey.slough

#convert counts to proportions and create new matrix
for (i in 1:nrow(spp.prey.pro)){
  for (j in 1:ncol(spp.prey.pro)){
    spp.prey.pro[i,j] <- spp.prey.pro[i,j]/rowSums(spp.prey.pro)[i]
  }
} 

spp.prey.cats <- spp.prey.pro 
spp.prey.cats[spp.prey.cats > 0 & spp.prey.cats < 0.05 ] <- 1
spp.prey.cats[spp.prey.cats >= 0.05 & spp.prey.cats < 0.3] <- 2  
spp.prey.cats[spp.prey.cats >= 0.3 & spp.prey.cats < 0.5] <- 3
spp.prey.cats[spp.prey.cats >= 0.5 & spp.prey.cats < 1] <- 4

setwd("C:/Users/peter/Documents/FIU/Research/Everglades Food-Web Diagram")
write.xlsx(spp.prey.cats, file = "SRS_slough_species_prey_categories.xlsx", rowNames = T)  
