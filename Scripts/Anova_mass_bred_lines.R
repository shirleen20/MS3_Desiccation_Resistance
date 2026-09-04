#_____________________________________30/07/2025_________________________________________________

# Perform ANOVA analysis on mass bred lines 

library("tidyverse")
library("ggpubr")
library("Rmisc")
library("car")
library("multcompView")
library("ggpubr")
library("emmeans")
library("broom")


rm(list=ls())


#_____________________________Import data __________________________________________


DF <- read_csv("Processed_data/Correctedcompletedataset.csv") 


#_____________Prepare data frames for mass bred lines for linear modelling__________________

DF1 <- DF %>% dplyr::filter(Line %in% c("CNnew", "SDnew", "CTnew"))      # all new mass bred lines

DF2 <- DF %>% dplyr::filter(Line %in% c("CNnew", "CNold", "SDnew", "SDold", "SDolder", "CTnew")) # all mass bred lines

DF3 <- DF %>% dplyr::filter(Line %in% c("CNnew", "CNold", "SDnew", "SDold", "CTnew"))  # all mass bred lines except S06


#______________________________Test the main effects of rtime ~ Line _____________________________

# One-way anova on corrected values (excluding S06)

M1 <- lm(rtime ~ Line, data = DF3, contrasts = list(Line = contr.sum)) %>% 
  car::Anova(type = 2) %>% broom::tidy() %>% 
  dplyr::mutate(Stars = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*",ifelse(p.value < 0.1, ".","")))))  



#___________________________Test the main effects of as.numeric(Generation) * Locality______________________________


# Two-way anova as.numeric(Generation) * Locality on corrected values (excluding S06)

M2 <- lm(rtime ~ as.numeric(Generation) * Locality, data = DF3, contrasts = list(Locality = contr.sum)) %>% 
  car::Anova(type = 2) %>% broom::tidy() %>% 
  dplyr::mutate(Stars = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*",ifelse(p.value < 0.1, ".","")))))  


M2.Emmeans <- lm(rtime ~ as.numeric(Generation) * Locality, data = DF3, contrasts = list(Locality = contr.sum)) %>% 
  emmeans(~interaction(Locality)) %>% 
  cld(Letter = "abcdefgh") %>% as_tibble() %>% bind_rows(.id = "id") 



#________________________Test the main effects of three new mass bred Line * Locality _____________________________

# One-way anova on corrected values (all new lines DF1)

M3 <- lm(rtime ~ Locality, data = DF1, contrasts = list(Locality = contr.sum)) %>% 
  car::Anova(type = 2) %>% broom::tidy() %>% 
  dplyr::mutate(Stars = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*",ifelse(p.value < 0.1, ".","")))))  


M3.Emmeans <- lm(rtime ~ Locality, data = DF1, contrasts = list(Locality = contr.sum)) %>% 
  emmeans(~Locality) %>% 
  cld(Letter = "abcdefgh") %>% as_tibble() %>% bind_rows(.id = "id") 


#________________________Test the main effects of all mass bred (no S06) Line * Locality _____________________________

# One-way anova on corrected values (all new lines DF3)

M4 <- lm(rtime ~ Locality, data = DF3, contrasts = list(Locality = contr.sum)) %>% 
  car::Anova(type = 2) %>% broom::tidy() %>% 
  dplyr::mutate(Stars = ifelse(p.value < 0.001, "***", ifelse(p.value < 0.01, "**", ifelse(p.value < 0.05, "*",ifelse(p.value < 0.1, ".","")))))  


M4.Emmeans <- lm(rtime ~ Locality, data = DF3, contrasts = list(Locality = contr.sum)) %>% 
  emmeans(~Locality) %>% 
  cld(Letter = "abcdefgh") %>% as_tibble() %>% bind_rows(.id = "id") 


#_______________END______________

