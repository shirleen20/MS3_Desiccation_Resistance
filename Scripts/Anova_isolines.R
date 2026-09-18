#_____________________________________01/08/2025_________________________________________________

# Perform ANOVA analysis on isolines.

library("tidyverse")
library("ggpubr")
library("Rmisc")
library("car")
library("multcomp")
library("multcompView")
library("lme4") # uses lmer() function for linear mixed models
library("emmeans")
library("broom")

rm(list=ls())


#_______________________________________Set 1 isolines ________________________________


# Import data for 1st batch isolines


DF <- read_csv("Processed_data/Batch_corrected_dataset.csv") 


# Prepare data frames to be used for linear modelling
# Note: 1st batch isoline analysed in two tranches DF1 tranch 1 and DF2 tranch 2
# DF_common has lines that were common in both tranches

DF1 <- DF %>% 
  dplyr::filter(Line %in% c("BR5","BR34","BR39","BR50","MB9","MB46","MB50","MB64","NAB2","NAB11","NAB28","UT1","UT3")) %>% 
  dplyr::filter(Generation %in% c("7", "8", "9", "10", "11", "12"))

DF2 <- DF %>% 
  dplyr::filter(Line %in% c("BR5","BR34","BR39","BR50","MB9","MB46","MB50","MB64","NAB2","NAB11","NAB28","UT1","UT3")) %>% 
  dplyr::filter(Generation %in% c("13", "15", "16", "17", "18")) #%>% 
#dplyr::filter(Line != "UT3") # removed UT3 from DF3 because it only has scores for G13 


DF3 <- DF %>% 
  dplyr::filter(Line %in% c("BR50","MB9","MB50","NAB11","NAB28","UT1","UT3")) 


#_________________Set 1 Tranch 1 analysis______________


# Do two way anova for 1st set isolines tranch 1 as.numeric(Generation)* Line

M1 <- lm(rtime ~ as.numeric(Generation) * Line, data = DF1, contrasts = list(Line = contr.sum)) %>% 
  car::Anova(type = 2) %>% broom::tidy() %>% 
  dplyr::mutate(Stars = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*",ifelse(p.value < 0.1, ".","")))))  


M1.Emmeans <- lm(rtime ~ as.numeric(Generation) * Line, data = DF1, contrasts = list(Line = contr.sum)) %>% 
  emmeans(~Line) %>% cld(Letter = "abcdefgh") %>% as_tibble() %>% bind_rows(.id = "id")



#______________Set 1 seven lines present in both tranches of set 1 __________________________

# "BR50","MB9","MB50","NAB11","NAB28","UT1","UT3"

M2 <- lm(rtime ~ as.numeric(Generation) * Line, data = DF3, contrasts = list(Line = contr.sum)) %>% 
  car::Anova(type = 2) %>% broom::tidy() %>% 
  dplyr::mutate(Stars = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*",ifelse(p.value < 0.1, ".","")))))  


M2.Emmeans <- lm(rtime ~ as.numeric(Generation) * Line, data = DF3, contrasts = list(Line = contr.sum)) %>% 
  emmeans(~Line) %>% cld(Letter = "abcdefgh") %>% as_tibble() %>% bind_rows(.id = "id")



#_______________________________Isofemale line anova - 2nd batch_______________________________________

#Prepare data frames to be used for linear modelling


DF4 <- DF %>% 
  dplyr::filter(Line %in% c("AS9","AS19","AS36","CT7","CT30","SYD6", "SYD13","SYD18","SYD53"))


M4 <- lm(rtime ~ as.numeric(Generation) * Line, data = DF4, contrasts = list(Line = contr.sum)) %>% 
  car::Anova(type = 2) %>% broom::tidy() %>% 
  dplyr::mutate(Stars = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*",ifelse(p.value < 0.1, ".","")))))  


#_______________END______________


