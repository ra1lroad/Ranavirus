#Ranavirus In Juvenile Metamorphs Analysis 
#Last Updated: Dec 25, 2023

#-----------------------------------------------------------------------------------------------------------------

#Loading and setting up the data, renaming some columns, adding a few columns, data wrangling etc. 

#Loading Packages
library(tidyverse)
library(ggplot2)
library(readr)
library(stringr)
library(lubridate)
library(RColorBrewer)
library(rstatix)
library(dplyr)
library(ggpubr)
library(gridExtra)

#set working directory
setwd("https://github.com/ra1lroad/Ranavirus.git")

#loading in the data
MetaTransmit <- read_csv("MetaTransmit.csv")

#---------------Data wrangling, cleaning up names etc

#getting rid of a couple extra lines that shouldn't be there
MetaTransmit <- MetaTransmit[-c(439:442),]

#renaming
names(MetaTransmit)[1] <- "Simple_ID"
names(MetaTransmit)[2] <- "Exposure_ID"
names(MetaTransmit)[3] <- "Death_timing"
names(MetaTransmit)[4] <- "Pre_exp_ID"
names(MetaTransmit)[5] <- "Exposure_type"
names(MetaTransmit)[6] <- "Exposure_start"
names(MetaTransmit)[7] <- "Exposure_duration"
names(MetaTransmit)[8] <- "Animal_type"
names(MetaTransmit)[9] <- "Death_date"
names(MetaTransmit)[12] <- "Erin_ID"

#making a few alterations
MetaTransmit <- MetaTransmit %>% 
  #Converting to Numerical, pulling out the "hrs"
  mutate(Exposure_duration_numerical = parse_number(Exposure_duration)) %>% 
  #making the start date a date 
  mutate(Exposure_start = as.Date(mdy(Exposure_start))) %>%
  #setting mean copy final as a numeric value
  mutate(MeanCopyFinal = as.numeric(MeanCopyFinal)) %>% 
  #Log of MeanCopyFinal
  #there's an issue with the zeroes.... they're -inf. To solve this, I add + 1 to everything.
  #It's not a great solution but it's what I have right now...
  mutate(logMeanCopyFinalPlus1 = log10(MeanCopyFinal + 1)) %>% 
  #MeanCopyFinal + 1 (for the glm linear model)
  mutate(MeanCopyFinalPlus1 = MeanCopyFinal + 1) %>% 
  #create the column Trio ID (for pair matching) and modify it from Exposure ID
  mutate(Trio_ID = str_replace(Exposure_ID, pattern = ".S[HQ].", replacement = "."))

#Larval Treatment IDs -- adding a new column
MetaTransmit$Larval_condition <- ifelse(grepl("ALS", MetaTransmit$Exposure_ID), "ALS", NA)
MetaTransmit$Larval_condition[grepl("AHS", MetaTransmit$Exposure_ID)] = "AHS"
MetaTransmit$Larval_condition[grepl("EHS", MetaTransmit$Exposure_ID)] = "EHS"
MetaTransmit$Larval_condition[grepl("ELS", MetaTransmit$Exposure_ID)] = "ELS"

#Creating a single column for Status using Status and Status_1 (from the redos)
MetaTransmit$Final_status[grepl("Neg", MetaTransmit$Status)] = "Neg"
MetaTransmit$Final_status[grepl("Pos", MetaTransmit$Status)] = "Pos"
MetaTransmit$Final_status[grepl("Neg", MetaTransmit$Status_1)] = "Neg"
MetaTransmit$Final_status[grepl("Pos", MetaTransmit$Status_1)] = "Pos"

#creating a binomial of infection status
MetaTransmit$Final_status_binomial[grepl("Pos", MetaTransmit$Final_status)] = "1"
MetaTransmit$Final_status_binomial[grepl("Neg", MetaTransmit$Final_status)] = "0"

#doing some factoring: 
MetaTransmit <- MetaTransmit %>% 
  #factor exposure duration in order
  mutate(Exposure_duration = factor(Exposure_duration, 
                                    levels = c("1hr", "6hr", "12hr", "24hr", "48hr"))) %>% 
  #factor larval condition in order
  mutate(Larval_condition = factor(Larval_condition, 
                                   levels = c("AHS", "ALS", "EHS", "ELS"))) %>% 
  
  #factor animal type in order
  mutate(Animal_type = factor(Animal_type, 
                              levels = c("Naive", "Infected","Control (pseudo-naive)", 
                                         "Control (pseudo-infected)")))

#-----------------------------------------------------------------------------------------------------------------
#Subseting out a few smaller data sets for future analysis. 
#Specifically: 
#  - Naive
#- Naive Sequential
#- Zeroes removed for viral load
#- Naive Shared
#- Zeroes removed for viral load
#- Infected
#- Control animals
#- Control pseudo-naive
#- Control pseudo-infected 

#----------Subset naive data
MetaTransmit_naive <- MetaTransmit %>% 
  subset(Animal_type== "Naive",)

#Subset Sequential Naives
MetaTransmitNaiveSQ <- MetaTransmit_naive %>% 
  subset(Exposure_type == "Sequential",)
#subset Shared Naives
MetaTransmitNaiveSH <- MetaTransmit_naive %>% 
  subset(Exposure_type == "Shared",)

# --------Getting rid of 0s for the shared and sequential data

#Naives Overall

#SHARED
#getting rid of the 0s as recommended
NoZeroesNaiveSH <- MetaTransmitNaiveSH %>% 
  subset(!(MeanCopyFinal == 0)) %>% 
  #Adding a variable called logmeancopyfinal without the zeroes
  mutate(logMeanCopyFinal = log(MeanCopyFinal)) 

#SEQUENTIAL
#getting rid of the 0s as Jesse Recommended
NoZeroesNaiveSQ <- MetaTransmitNaiveSQ %>% 
  subset(!(MeanCopyFinal == 0)) %>% 
  #Adding a variable called logmeancopyfinal without the zeroes
  mutate(logMeanCopyFinal = log(MeanCopyFinal)) 

#-----------Subseting out virus animals
MetaTransmit_virus <- MetaTransmit %>% 
  subset(MetaTransmit$Animal_type=="Infected",)

#------------Subsetting  control animals 
MetaTransmit_controlPN <- MetaTransmit %>% 
  subset(Animal_type=="Control (pseudo-naive)",)
MetaTransmit_controlPV <- MetaTransmit %>% 
  subset(MetaTransmit$Animal_type=="Control (pseudo-infected)",)

#-----------------------------------------------------------------------------------------------------------------

#PLOTS: 

#Plotting the start dates
ggplot(MetaTransmit, aes(x=Exposure_start, fill = Exposure_type)) +
  geom_bar()
#is that one start date in the middle an issue? Seems like it might have been recorded wrong. 

#-----------------------------------------------------------------------------------------------------------------

#Suspected Ranavirus Mortalities (Focal Naive)

#facet labels
# New facet label names for exposure type variable
exposure.labs <- c("Indirect Contact", "Direct Contact")
names(exposure.labs) <- c("Sequential", "Shared")

#plot naive animal mortality
ggplot(MetaTransmit_naive, aes(x= Exposure_duration , y= lengths(Exposure_duration), fill = Death_timing)) + 
  geom_bar(stat= "identity") +
  guides(fill=guide_legend(title=NULL)) + 
  facet_wrap(~Exposure_type, 
             labeller = labeller(Exposure_type = exposure.labs)
  ) +
  labs(x = "Time Between Exposures         Exposure Period", y = "Susceptible Animal Counts") + 
  theme(axis.title.y = element_text(size = 18, hjust = 0.7),
        axis.title.x = element_text(size = 18, hjust = 0.3),
        strip.text.x = element_text(size = 20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.text=element_text(size=15),
        axis.title = element_text(size=20, color = "black"),
        legend.position = "bottom") + 
  scale_fill_brewer(palette="Set2")

ggsave(
  'suspectedMortalities.jpg',
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = NA,
  height = NA,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL
)

#-----------------------------------------------------------------------------------------------------------------

#Infection Status (Focal Naive)
#plot naive animal mortality
ggplot(MetaTransmit_naive, aes(x= Exposure_duration , y= lengths(Exposure_duration), fill = Final_status)) + 
  geom_bar(stat= "identity") +
  guides(fill=guide_legend(title=NULL)) + 
  facet_wrap(~Exposure_type, 
             labeller = labeller(Exposure_type = exposure.labs)
  ) +
  labs(x = "Time Between Exposures         Exposure Period", y = "Susceptible Animal Counts") + 
  theme(axis.title.y = element_text(size = 18, hjust = 0.7),
        axis.title.x = element_text(size = 18, hjust = 0.3),
        strip.text.x = element_text(size = 20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.text=element_text(size=15),
        axis.title = element_text(size=20, color = "black"),
        legend.position = "bottom") + 
  #scale_fill_brewer(palette="Set2")
  scale_fill_manual(breaks = c("Neg", "Pos", NA), 
                    values=c("lightblue", "darkred", "grey"))

ggsave(
  'virusStatusNaive.jpg',
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = NA,
  height = NA,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL
)

summary(MetaTransmitNaiveSH$Final_status["Neg"])
summary(MetaTransmitNaiveSH$Final_status["Pos"])

#summary(MetaTransmitNaiveSH$Final_status_binomial)
length(MetaTransmitNaiveSH$Final_status_binomial)
summary(length(MetaTransmitNaiveSH$Final_status_binomial))

summary(MetaTransmitNaiveSH$Exposure_type)
summary(MetaTransmitNaiveSQ$Exposure_type)

#-----------------------------------------------------------------------------------------------------------------

#Virus Animals
#plot mortality
ggplot(MetaTransmit_virus, aes(x= Exposure_duration , y= lengths(Exposure_duration), fill = Final_status)) + 
  geom_bar(stat= "identity") +
  guides(fill=guide_legend(title=NULL)) + 
  labs(x = "Exposure Period", y = "Infective Animal Counts") + 
  theme(axis.title.y = element_text(size = 18, hjust = 0.7),
        axis.title.x = element_text(size = 18, hjust = 0.5),
        strip.text.x = element_text(size = 20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.text=element_text(size=15),
        axis.title = element_text(size=20, color = "black"),
        legend.position = "bottom") + 
  #scale_fill_brewer(palette="Set2")
  scale_fill_manual(breaks = c("Neg", "Pos", NA), 
                    values=c("lightblue", "darkred", "grey"))

ggsave(
  'virusStatusVirus.jpg',
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = NA,
  height = NA,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL
)

#-----------------------------------------------------------------------------------------------------------------

#Control Naive Animals
#plot naive animal mortality
ggplot(MetaTransmit_controlPN, aes(x= Exposure_duration , y= lengths(Exposure_duration), fill = Final_status)) + 
  geom_bar(stat= "identity") +
  guides(fill=guide_legend(title=NULL)) + 
  facet_wrap(~Exposure_type, 
             labeller = labeller(Exposure_type = exposure.labs)
  ) +
  labs(x = "Time Between Exposures         Exposure Period", y = "Control Susceptible Counts") + 
  theme(axis.title.y = element_text(size = 18, hjust = 0.7),
        axis.title.x = element_text(size = 18, hjust = 0.3),
        strip.text.x = element_text(size = 20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.text=element_text(size=15),
        axis.title = element_text(size=20, color = "black"),
        legend.position = "bottom") + 
  #scale_fill_brewer(palette="Set2")
  scale_fill_manual(breaks = c("Neg", "Pos", NA), 
                    values=c("lightblue", "darkred", "grey"))

ggsave(
  'virusStatusControlPseudonaive.jpg',
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = NA,
  height = NA,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL
)

#-----------------------------------------------------------------------------------------------------------------

#Control Pseudovirus Animals: 
#plot naive animal mortality
ggplot(MetaTransmit_controlPV, aes(x= Exposure_duration , y= lengths(Exposure_duration), fill = Final_status)) + 
  geom_bar(stat= "identity") +
  guides(fill=guide_legend(title=NULL)) + 
  labs(x = "Exposure Period", y = "Control Pseudo-infectives Counts") + 
  theme(axis.title.y = element_text(size = 15, hjust = 0.7),
        axis.title.x = element_text(size = 18, hjust = 0.5),
        strip.text.x = element_text(size = 20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        legend.text=element_text(size=15),
        axis.title = element_text(size=20, color = "black"),
        legend.position = "bottom") + 
  #scale_fill_brewer(palette="Set2")
  scale_fill_manual(breaks = c("Neg", "Pos", NA), 
                    values=c("lightblue", "darkred", "grey"))

ggsave(
  'virusStatusControlPseudovirus.jpg',
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = NA,
  height = NA,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL
)

#-----------------------------------------------------------------------------------------------------------------

#Performing a few data visualizations on viral load and exposure duration: 
#-------------Boxplots

#Violin Box Plot Shared
ggplot(NoZeroesNaiveSH, aes(x=Exposure_duration, y= log(MeanCopyFinal)))+
  geom_violin()+
  geom_boxplot(alpha=0.3) + 
  labs(x = "Exposure Period", y = "Log of Viral Load") 
ggsave("NaiveSHViralLoadViolinNoZeroes.jpg", plot = last_plot())

#Violin Box Plot Sequential
ggplot(NoZeroesNaiveSQ, aes(x=Exposure_duration, y= log(MeanCopyFinal)))+
  geom_violin()+
  geom_boxplot(alpha=0.3) + 
  labs(x = "Time Between Exposures", y = "Log of Viral Load")
ggsave("NaiveSQViralLoadViolinNoZeroes.jpg", plot = last_plot())

# ------------------
#5 number summaries - a bit clunky but here it is. 

#SEQUENTIAL
simpleSQ <- NoZeroesNaiveSQ %>% 
  select(Exposure_duration, logMeanCopyFinal) %>% 
  mutate(Exposure_duration = as.factor(Exposure_duration))
#1hr
simpleSQ1hr <- simpleSQ %>% 
  subset(Exposure_duration == "1hr") 
summary(simpleSQ1hr$logMeanCopyFinal)
#6hr
simpleSQ6hr <- simpleSQ %>% 
  subset(Exposure_duration == "6hr") 
summary(simpleSQ6hr$logMeanCopyFinal)
#12hr
simpleSQ12hr <- simpleSQ %>% 
  subset(Exposure_duration == "12hr") 
summary(simpleSQ12hr$logMeanCopyFinal)
#24hr
simpleSQ24hr <- simpleSQ %>% 
  subset(Exposure_duration == "24hr") 
summary(simpleSQ24hr$logMeanCopyFinal)
#48hr
simpleSQ48hr <- simpleSQ %>% 
  subset(Exposure_duration == "48hr") 
summary(simpleSQ48hr$logMeanCopyFinal)

#SHARED
#5 number summaries - a bit clunky but here it is. 
simpleSH <- NoZeroesNaiveSH %>% 
  select(Exposure_duration, logMeanCopyFinal) %>% 
  mutate(Exposure_duration = as.factor(Exposure_duration))
#1hr
simpleSH1hr <- simpleSH %>% 
  subset(Exposure_duration == "1hr") 
summary(simpleSH1hr$logMeanCopyFinal)
#6hr
simpleSH6hr <- simpleSH %>% 
  subset(Exposure_duration == "6hr") 
summary(simpleSH6hr$logMeanCopyFinal)
#12hr
simpleSH12hr <- simpleSH %>% 
  subset(Exposure_duration == "12hr") 
summary(simpleSH12hr$logMeanCopyFinal)
#24hr
simpleSH24hr <- simpleSH %>% 
  subset(Exposure_duration == "24hr") 
summary(simpleSH24hr$logMeanCopyFinal)
#48hr
simpleSH48hr <- simpleSH %>% 
  subset(Exposure_duration == "48hr") 
summary(simpleSH48hr$logMeanCopyFinal)

#-----------------------------------------------------------------------------------------------------------------
#violin boxplots of shared and sequential facet_wrapped
ggplot(NoZeroesNaive, aes(x=Exposure_duration, y= log(MeanCopyFinal), 
                          fill = Exposure_duration)) +
  geom_violin()+
  geom_boxplot(alpha=0.3) + 
  facet_wrap(~Exposure_type, 
             labeller = labeller(Exposure_type = exposure.labs)
  ) +
  labs(x = "Time Between Exposures              Exposure Period", y = "Log of Viral Load", size = 20) +
  scale_fill_brewer(palette = "YlOrRd") +
  theme(legend.position = "none", 
        axis.text = element_text(size = 15), 
        axis.title = element_text(size = 17), 
        strip.text.x = element_text(size = 20))
ggsave("ViralLoadViolinNoZeroes.jpg", plot = last_plot())

#-----------------------------------------------------------------------------------------------------------------
#Linear Models for Viral Load and Exposure Durations by Experiment: 
#assuming no intercept
#V ≈ b*Duration
#log(V) ≈ log(b) + log(Duration)

#the linear model shared
log_transformed_lm_SH = glm(log(MeanCopyFinal) ~ log(Exposure_duration_numerical), data = NoZeroesNaiveSH)
summary(log_transformed_lm_SH)
plot(log_transformed_lm_SH)

#the linear model sequential
log_transformed_lm_SQ = glm(log(MeanCopyFinal) ~ log(Exposure_duration_numerical), data = NoZeroesNaiveSQ)
summary(log_transformed_lm_SQ)
plot(log_transformed_lm_SQ)

#-----------------------------------------------------------------------------------------------------------------
#Paired Pearson's Linear Model: 
#The goal here is to compare the viral load of the infected animal with the paired naive animal 
#to see if an animal gets a more severe infection from an infected animal with a higher viral load. 

#Ok so this is pairing the infected animal with the naive animal by ID for the Shared housing experiment. 
Shared_grouped <- MetaTransmit %>% 
  select(Exposure_ID, Exposure_type, logMeanCopyFinalPlus1, Animal_type, Exposure_duration) %>% 
  filter(Exposure_type == "Shared") %>% 
  pivot_wider(names_from = Animal_type, values_from = logMeanCopyFinalPlus1)

shared_pair_plot <- Shared_grouped %>% 
  ggplot(aes(x = Infected, y = Naive)) +
  geom_point(width = 5, height = 4, dpi = 300, units = "in", device='png') + 
  geom_smooth(method = lm, se = FALSE) + 
  stat_cor(method = "pearson", show.legend = FALSE) + 
  xlab("Infective Frog Log Viral Load") +
  ylab("") +
  theme_minimal() +
  ggtitle("Direct Contact Pairs") +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 15), 
        plot.title = element_text(size = 20, hjust = .5))
ggsave("SharedPair.jpg", plot = last_plot())

#Sequential Housing
Sequential_grouped <- MetaTransmit %>% 
  select(Trio_ID, Exposure_type, logMeanCopyFinalPlus1, Animal_type, Exposure_duration) %>% 
  filter(Exposure_type == "Shared" & Animal_type == "Infected" | 
           Exposure_type == "Sequential" & Animal_type == "Naive") %>% 
  select(Trio_ID, Animal_type, logMeanCopyFinalPlus1, Exposure_duration) %>% 
  pivot_wider(names_from = Animal_type, values_from = logMeanCopyFinalPlus1)

sequential_paired_plot <- Sequential_grouped %>% 
  ggplot(aes(x = Infected, y = Naive)) +
  geom_point(width = 5, height = 4, dpi = 300, units = "in", device='png') + 
  geom_smooth(method = lm, se = FALSE) + 
  stat_cor(method = "pearson", show.legend = FALSE) +
  xlab("Infective Frog Log Viral Load") +
  ylab("Susceptible Frog Log Viral Load") +
  theme_minimal() +
  ggtitle("Indirect Contact Pairs") +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 15), 
        plot.title = element_text(size = 20, hjust = .5))
ggsave("SequentialPair.jpg", plot = last_plot())

#combining the plots
#putting all the plots together
require(gridExtra)
Combined_pair_plots <- grid.arrange(sequential_paired_plot, shared_pair_plot, ncol =2, nrow =1)

ggsave("Sequential_and_Shared_Pair_BothPlots.jpg", plot = Combined_pair_plots)

#-----------------------------------------------------------------------------------------------------------------
#Testing for an effect of salinity and heat treatments on viral load: 
#Because all of the animals came out of Nicole's salinity and heat treatments as tadpoles.
#--------------Visualizing it first
#SHARED
ggplot(data = NoZeroesNaiveSH, aes(x=Larval_condition, y = log(MeanCopyFinal), fill = Larval_condition)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "BrBG") +
  theme(legend.position = "none")
#SEQUENTIAL
ggplot(data = NoZeroesNaiveSQ, aes(x=Larval_condition, y = log(MeanCopyFinal), fill = Larval_condition)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "BrBG") +
  theme(legend.position = "none")

#--------------1-way ANOVA
set.seed(1234)

#SHARED
# Compute the analysis of variance
res_aov_SH <- aov(log(MeanCopyFinal) ~ Larval_condition, data = NoZeroesNaiveSH)
# Summary of the analysis
summary(res_aov_SH)

#SEQUENTIAL
# Compute the analysis of variance
res_aov_SQ <- aov(log(MeanCopyFinal) ~ Larval_condition, data = NoZeroesNaiveSQ)
# Summary of the analysis
summary(res_aov_SQ)



