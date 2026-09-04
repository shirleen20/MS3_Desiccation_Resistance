#_____________________________________________30/07/2025_____________________________________________

# We analyse values for median knock down times recorded for desiccation resistance bioassay
# Correct values for batch effects using S06 (use combined desiccation resistance data for mass breds, old & new isofemale lines)

library("tidyverse")
library("dplyr")

rm(list=ls())

#___________________Import data___________________________________

DF <- read_csv("Data/Combined_data_raw.csv") 

# Notes:
# In DF, hours refers to the median knock down time and Mean refers to the mean knockdown time for DR
# Add to DF the actual generation numbers for S06 i.e:
# Bioassay <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)
# Generation <- c(120, 122, 123, 124, 125, 126, 128, 129, 130, 131, 132, 133, 134)


#____________________________________Correct for Batch effects_____________________________________

# filter S06 data to work on the S06 subset:

DF_S06 <- DF %>%
  dplyr::filter(Locality == "S06") %>%  
  dplyr::filter(Hours!= "NA") # removes all the missing values in the Hours column (there aren't any NA for S06). 


#_____________________Calculate global mean and experimental mean values_____________________________

# Use the filtered S06 data to find global mean and experimental mean: 
# note in Df_S06 below column "mean" refers to global mean of S06 and 
# "emean" refers to experimental mean forS06

Df_S06 <- DF_S06 %>% 
  dplyr::mutate(mean= round(mean(Hours, na.rm = TRUE),2)) %>% #calculate global mean knockdown time 
  dplyr::group_by(Bioassay) %>%  #groups the observations (Hours) by the Bioassay (NB: "Bioassay" is not the same as "Generation") 
  dplyr::mutate(emean = round(mean(Hours, na.rm = TRUE),2)) %>% #calculate experimental mean for S06
  dplyr::ungroup() %>%  # Ungroup the dataframe, because we want to carry the two calculated values (emean and mean)
  dplyr::group_by(Bioassay, emean, mean) %>% 
  dplyr::summarise()


#___________________Calculate the corrected time for batch effects____________________________________

# Create main dataframe by joining Df_S06 that has experimental and global mean values of 
# S06 using Bioassay as an identifier:

DF_main <- DF %>% 
  left_join(Df_S06, by = "Bioassay") %>% # Combine DF with the S06's global mean and experimental mean knock time times to calculate the adjusted time response. 
  dplyr::group_by(Bioassay) %>%          # Combine by experiments to calculate the adjusted time for each experiment
  dplyr::mutate(rtime = round(Hours/(emean/mean), 2)) # use round() to reduce decimal places. 

Correctedcompletedataset <- DF_main %>% 
  dplyr::mutate(Line = recode(Line,  CanberraNew = "CNnew",CanberraOld  = "CNold", CapeTribulationNew = "CTnew",
                              SydneyNew = "SDnew",SydneyOld  = "SDold", S06 = "SDolder")) %>% 
  dplyr::mutate(Locality = recode(Locality,  S06 = "Sydney")) #%>% 
  #write.csv("Processed_data/Correctedcompletedataset.csv")    # Save dataframe to be used for analysis:


# Save corrected data

Corrected_data <- read_csv("Processed_data/Correctedcompletedataset.csv")


#_______________________________END______________________________


