#_____________________________________30/07/2025_________________________________________________

# We make plots for mass bred and isolines in  Figure1 for MS3
# we also calculate the mean of median survival time for SDolder compared to other mass breds 
# we than calculate the mean of median survival time for SDnew and CTnew compared to all other isolines

library("tidyverse")
library("ggpubr")
library("Rmisc")
library("car")
library("ggpubr")
library("patchwork")
library("cowplot")

rm(list=ls())

#__________________________Import data corrected for batch effects____________________________

DF <- read_csv("Processed_data/Batch_corrected_dataset.csv")

# Here we want to show all mass breds except S06 using corrected values but S06 using uncorrected values
# 1st create dataframe with corrected values for all mass breds except S06:

DF_mass <- DF %>% 
  dplyr::filter(Line %in% c("CNnew", "CNold", "SDnew", "SDold", "CTnew")) %>% 
  dplyr::select(Bioassay, Generation, Locality, Line, Replicate, rtime)

# 2nd create dataframe with uncorrected values for S06:

DF_massS06 <- DF %>% 
  dplyr::filter(Line %in% c("SDolder")) %>% 
  dplyr::select(Bioassay, Generation, Locality, Line, Replicate, Hours) %>% 
  dplyr::rename(rtime = Hours) # this df contains uncorrected values for So6 but we rename the column Hour to rtime so we can merge data with DF_mass for plotting

# Merge into one dataframe containing corrected values for all other mass breds and uncorrected values for S06:

DF.mass <- rbind(DF_mass, DF_massS06) %>% 
  summarySE(groupvars = c("Line", "Generation"), measurevar = "rtime") 

# Calculate overall mean rtime for all other mass bred lines corrected values and SDolder uncorrected values

Overall.mean <- DF.mass %>%
  summarySE(groupvars = c("Line"), measurevar = "rtime") %>% 
  dplyr::select(Line, rtime, se)


# Calculate the mean of median survival time for SDolder compared to other mass breds

#Mass <- c("CNnew", "CNold", "SDnew", "SDold", "CTnew")
#Older <- c("SDolder")

#SDoldervsothers <- DF.mass %>%
#  dplyr::mutate(Type = ifelse(Line %in% Mass, "Mass",  "Older"))%>% 
#  summarySE(groupvars = c("Type"), measurevar = "rtime") # calculate mean of median knockdown time uncorrected for batch

#_____________make plot_____________________________

DF.massp1 <- DF.mass %>% 
  dplyr::mutate(Line = factor(as.factor(Line), levels = c("SDnew","SDold", "SDolder", "CNnew","CNold", "CTnew"))) %>% 
  dplyr::mutate(Generation = factor(as.factor(Generation), levels = c("0","1","2","3","4","5","6","7","8","9","10",
                                                                      "21","22","23","24","28","29","30","31", 
                                                                      "120", "122", "123", "124", "125", "126", "128", "129", "130", "131", "132", "133", "134"))) %>% 
  mutate(Line2 = dplyr::recode(Line, "CTnew" = "CTnew (37.72±1.88)", "SDnew" = "SDnew (43.38±2.59)","SDold" = "SDold (38.63±2.01)","SDolder"= "SDolder (32.31±0.65)*", "CNnew" = "CNnew (42.05±2.78)","CNold" = "CNold (35.51±2.45)"))

# Prepare numeric Generation column
DF.massp1 <- DF.massp1 %>%
  mutate(GenNum = as.numeric(as.character(Generation)))

# Compute global y-range (include errorbars) 
# include errorbars so nothing gets clipped
all_y_min <- min(DF.massp1$rtime - DF.massp1$se, na.rm = TRUE)
all_y_max <- max(DF.massp1$rtime + DF.massp1$se, na.rm = TRUE)

# ensure 30 and 50 are inside the limits
y_low  <- min(all_y_min, 30)
y_high <- max(all_y_max, 50)

# optionally add a tiny margin so points don't sit on the axis line
margin <- 0.02 * (y_high - y_low)
y_low  <- y_low - margin
y_high <- y_high + margin

# common y scale that'll be applied to all plots
common_y_scale <- scale_y_continuous(limits = c(y_low, y_high),breaks = c(30, 40, 50),expand = expansion(mult = c(0,0)))

# helper: small theme chunk 
common_theme <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 8),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 8, color = "black"),
    axis.ticks.y = element_line(),
    axis.text.x = element_text(angle = 0, hjust = 1),
    panel.border = element_rect(colour = "black"),
    plot.margin = margin(5,5,5,5))

#  plotting functions using common_y_scale
make_first_col_plot <- function(df, label) {
  df %>% 
    filter(Line2 == label) %>%
    ggplot(aes(x = GenNum, y = rtime, group = 1)) +
    geom_point(colour = "steelblue", size = 0.5) +
    geom_line(colour = "steelblue") +
    geom_errorbar(aes(ymin = rtime - se, ymax = rtime + se),width = 0.1, size = 0.5, alpha = 0.5) +
    scale_x_continuous(breaks = 0:10, limits = c(0,10)) +
    common_y_scale +
    labs(title = label, x = NULL, y = NULL) +
    common_theme
}

make_other_plot <- function(df, label) {
  df %>% 
    filter(Line2 == label) %>%
    ggplot(aes(x = GenNum, y = rtime, group = 1)) +
    geom_point(colour = "steelblue", size = 0.5) +
    geom_line(colour = "steelblue") +
    geom_errorbar(aes(ymin = rtime - se, ymax = rtime + se),width = 0.1, size = 0.5, alpha = 0.5) +
    scale_x_continuous(breaks = 21:31,limits = c(21, 31)) +
    common_y_scale +
    labs(title = label, x = NULL, y = NULL) +
    common_theme
}

make_SDolder_plot <- function(df, label) {
  df %>% 
    filter(Line2 == label) %>%
    ggplot(aes(x = GenNum, y = rtime, group = 1)) +
    geom_point(colour = "steelblue", size = 0.5) +
    geom_line(colour = "steelblue") +
    geom_errorbar(aes(ymin = rtime - se, ymax = rtime + se),width = 0.1, size = 0.5, alpha = 0.5) +
    scale_x_continuous(breaks = seq(120, 135, 5), limits = c(120, 135)) +
    common_y_scale +
    labs(title = label, x = NULL, y = NULL) +
    common_theme
}

# create plots (as before)
p_SDnew   <- make_first_col_plot(DF.massp1, "SDnew (43.38±2.59)")
p_CNnew   <- make_first_col_plot(DF.massp1, "CNnew (42.05±2.78)")
p_CTnew   <- make_first_col_plot(DF.massp1, "CTnew (37.72±1.88)")

p_SDold   <- make_other_plot(DF.massp1, "SDold (38.63±2.01)")
p_CNold   <- make_other_plot(DF.massp1, "CNold (35.51±2.45)")

p_SDolder <- make_SDolder_plot(DF.massp1, "SDolder (32.31±0.65)*")

# Stack columns (unchanged)
col1 <- p_SDnew / p_CNnew / p_CTnew + plot_layout(heights = rep(1,3))
col2 <- p_SDold / p_CNold / plot_spacer() + plot_layout(heights = rep(1,3))
col3 <- p_SDolder / plot_spacer() / plot_spacer() + plot_layout(heights = rep(1,3))


#create a clean empty placeholder plot for empty cells
empty_plot <- ggplot() + theme_void()

# layout (3 rows x 3 cols):# Row 1 (top): SDnew, SDold, SDolder, # Row 2 (mid): CNnew, CNold, empty, # Row 3 (bot): CTnew, empty,  empty
grid_plots <- list( p_SDnew,p_SDold,p_SDolder,p_CNnew,p_CNold,empty_plot, p_CTnew,empty_plot, empty_plot)

# combine into a regular grid with cowplot::plot_grid
# ncol = 3 makes a 3x3 grid and ensures each cell has the same width/height

main_grid <- plot_grid(plotlist = grid_plots, ncol = 3, align = "v", axis = "tblr")# align vertically, draw axes where present(top,bottom,left,right as appropriate)

# add a centered title using cowplot
title_row <- ggdraw() + draw_label("Mass bred lines", fontface = "bold", size = 8, x = 0.5, hjust = 0.5)

# stack title and grid
p1 <- plot_grid(title_row, main_grid, ncol = 1, rel_heights = c(0.06, 1))


#_______________________________________________

#2. Line graph showing effect of domestication in set 1 isofemale lines

DF_set1 <- DF %>% 
  dplyr::filter(Locality %in% c("Brisbane", "Mareeba", "Narrabri", "UtcheeCreek")) %>% 
  summarySE(groupvars = c("Line", "Generation"), measurevar = "rtime") 

# Calculate overall mean rtime for all set1 isolines corrected values 

Overall.mean <- DF_set1 %>%
  summarySE(groupvars = c("Line"), measurevar = "rtime") %>% 
  dplyr::select(Line, rtime, se)


p2 <- DF_set1 %>% 
  dplyr::mutate(Line = factor(as.factor(Line), levels = c("BR5","BR34","BR39","BR50","MB9","MB46","MB50","MB64",
                                                          "NAB2","NAB11","NAB28","UT1","UT3")),
  #dplyr::mutate(Generation = factor(as.factor(Generation), levels = c("7","9","10","11","12","13","15","16","17","18"))) %>%
  Generation = as.numeric(as.character(Generation))) %>%  # convert to numeric
   mutate(Line = dplyr::recode(Line,
                              "BR5" = "BR5 (36.55±1.22)","BR34" = "BR34 (34.36±2.97)","BR39" = "BR39 (34.66±1.81)","BR50" = "BR50 (32.38±1.08)",
                              "MB9" = "MB9 (36.65±1.15)", "MB46" = "MB46 (35.48±2.81)", "MB50" = "MB50 (36.31±2.15)", "MB64" = "MB64 (36.27±2.30)",
                              "NAB2" = "NB2 (32.24±2.47)", "NAB11" = "NB11 (31.86±1.08)","NAB28" = "NB28 (39.72±1.16)",
                              "UT1" = "UT1 (32.59±1.74)", "UT3" = "UT3 (37.51±1.75)"))  %>% 
  ggplot(aes(x = Generation , y = rtime, group=1)) +
  geom_point(colour = "steelblue", size = 0.5) + 
  geom_line(colour="steelblue") +
  geom_errorbar(aes(ymin = rtime - se, ymax = rtime + se), width = 0.1, size = 0.5, alpha = 0.5) +
  labs(title="Set 1 isofemale lines", x= "Isofemale lines (Set 1)", y="Mortality (hours)")+
  facet_wrap(~ Line, ncol = 4, labeller = label_wrap_gen(multi_line = TRUE)) +
  scale_x_continuous(breaks = 7:18, limits = c(7, 18)) +
  scale_y_continuous(limits = c(y_low, y_high),breaks = c(30, 40, 50),expand = expansion(mult = c(0,0)))+
  theme_bw() +
  theme(strip.background = element_blank(),
    strip.text = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_blank(),
    panel.border = element_rect(colour = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 8)) 
  


#2. Line graph showing effect of domestication in set 2 isofemale lines

DF_set2 <- DF %>% 
  dplyr::filter(Line %in% c("AS9","AS19","AS36","CT7","CT30","SYD6", "SYD13","SYD18","SYD53")) %>% 
  summarySE(groupvars = c("Line", "Generation"), measurevar = "rtime") 

# Calculate overall mean rtime for all set2isolines corrected values 

Overall.mean <- DF_set2 %>%
  summarySE(groupvars = c("Line"), measurevar = "rtime") %>% 
  dplyr::select(Line, rtime, se)

# Compute global y-range (include errorbars) 
# include errorbars so nothing gets clipped
all_y_min <- min(DF_set2$rtime - DF_set2$se, na.rm = TRUE)
all_y_max <- max(DF_set2$rtime + DF_set2$se, na.rm = TRUE)

# ensure 30 and 50 are inside the limits
y_low  <- min(all_y_min, 30)
y_high <- max(all_y_max, 50)

# optionally add a tiny margin so points don't sit on the axis line
margin <- 0.02 * (y_high - y_low)
y_low  <- y_low - margin
y_high <- y_high + margin


p3 <- DF_set2 %>% 
  dplyr::mutate(Line = factor(as.factor(Line), levels = c("AS9","AS19","AS36","CT7","CT30","SYD6", "SYD13","SYD18","SYD53")),
  Generation = as.numeric(as.character(Generation))) %>% # convert to numeric
  mutate(Line = dplyr::recode(Line,
                              "AS9" = "AS9 (38.91±3.05)","AS19" = "AS19 (36.60±2.05)","AS36" = "AS36 (35.25±1.36)",
                              "CT7" = "CT7 (37.37±2.47)", "CT30" = "CT30 (39.38±2.73)", 
                              "SYD6" = "SD6 (31.32±2.33)", "SYD13" = "SD13 (40.53±2.70)", 
                              "SYD18" = "SD18 (35.54±3.26)", "SYD53" = "SD53 (36.36±1.58)"))  %>%
  ggplot(aes(x = Generation , y = rtime, group=1)) +
  geom_point(colour = "steelblue", size = 0.5) + 
  geom_line(colour="steelblue") +
  geom_errorbar(aes(ymin = rtime - se, ymax = rtime + se), width = 0.1, size = 0.5, alpha = 0.5) +
  labs(title="Set 2 isofemale lines", x= "Isofemale lines (Set 2)", y="Mortality (hours)")+
  facet_wrap(~ Line, ncol = 4, labeller = label_wrap_gen(multi_line = TRUE)) +
  scale_x_continuous(breaks = 1:9, limits = c(2, 9)) +
  scale_y_continuous(limits = c(y_low, y_high),breaks = c(30, 40, 50),expand = expansion(mult = c(0,0)))+
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        panel.border = element_rect(colour = "black"),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 8)) 



# Arrange plots 

Figure1 <- ggarrange(p1,p2,p3, nrow = 3, heights =c(4.0,4.5,3.4))


# Annotate and save the figure by adding a common labels

Figure1_annotated <- annotate_figure(Figure1,bottom = text_grob("Generation", color = "black",hjust = 0.5, face = "bold", size = 8),
                                     left = text_grob("Median survival time (hours)", color = "black",rot = 90,face = "bold", size = 8))

ggsave(plot = Figure1_annotated, width = 7.0, height = 9.0, units = "in", dpi = 300,filename = "Figures/Figure 1 showing desiccation resistance in mass bred vs isofemale lines.jpg")

#___________________________________________________________________________________________________________________________

#_________________________________________________________________________________________________________________________




#__________Figure S1: mass bred lines values not corrected for batch effects__________________________


rm(list=ls())


#__________________________Import data corrected for batch effects____________________________

DF <- read_csv("Processed_data/Batch_corrected_dataset.csv")

# Here we want to show all mass breds using uncorrected values

DF_mass <- DF %>% 
  dplyr::filter(Line %in% c("CNnew", "CNold", "SDnew", "SDold", "CTnew", "SDolder")) %>% 
  dplyr::select(Bioassay, Generation, Locality, Line, Replicate, Hours)
  

DF.mass <- DF_mass %>% 
  summarySE(groupvars = c("Line", "Generation"), measurevar = "Hours") 

# Calculate overall mean "Hours for all other mass bred lines uncorrected values

Overall.mean <- DF.mass %>%
  summarySE(groupvars = c("Line"), measurevar = "Hours")%>% 
  dplyr::select(Line, Hours, se)


#_____________make plot_____________________________

DF.massp1 <- DF.mass %>% 
  dplyr::mutate(Line = factor(as.factor(Line), levels = c("SDnew","SDold", "SDolder", "CNnew","CNold", "CTnew"))) %>% 
  dplyr::mutate(Generation = factor(as.factor(Generation), levels = c("0","1","2","3","4","5","6","7","8","9","10",
                                                                      "21","22","23","24","28","29","30","31", 
                                                                      "120", "122", "123", "124", "125", "126", "128", "129", "130", "131", "132", "133", "134"))) %>% 
  mutate(Line2 = dplyr::recode(Line, "CTnew" = "CTnew (37.72±1.74)", "SDnew" = "SDnew (42.88±1.76)","SDold" = "SDold (41.25±2.46)","SDolder"= "SDolder (32.31±0.65)*", "CNnew" = "CNnew (45.00±4.09)","CNold" = "CNold (38.00±3.18)"))

# Prepare numeric Generation column
DF.massp1 <- DF.massp1 %>%
  mutate(GenNum = as.numeric(as.character(Generation)))

# Compute global y-range (include errorbars) 
# include errorbars so nothing gets clipped
all_y_min <- min(DF.massp1$Hours - DF.massp1$se, na.rm = TRUE)
all_y_max <- max(DF.massp1$Hours + DF.massp1$se, na.rm = TRUE)

# ensure 30 and 50 are inside the limits
y_low  <- min(all_y_min, 30)
y_high <- max(all_y_max, 53)

# optionally add a tiny margin so points don't sit on the axis line
margin <- 0.02 * (y_high - y_low)
y_low  <- y_low - margin
y_high <- y_high + margin

# common y scale that'll be applied to all plots
common_y_scale <- scale_y_continuous(limits = c(y_low, y_high),breaks = c(30, 40, 50),expand = expansion(mult = c(0,0)))

# helper: small theme chunk 
common_theme <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 8),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.ticks.y = element_line(),
        axis.text.x = element_text(angle = 0, hjust = 1),
        panel.border = element_rect(colour = "black"),
        plot.margin = margin(5,5,5,5))

#  plotting functions using common_y_scale
make_first_col_plot <- function(df, label) {
  df %>% 
    filter(Line2 == label) %>%
    ggplot(aes(x = GenNum, y = Hours, group = 1)) +
    geom_point(colour = "steelblue", size = 0.5) +
    geom_line(colour = "steelblue") +
    geom_errorbar(aes(ymin = Hours - se, ymax = Hours + se),width = 0.1, size = 0.5, alpha = 0.5) +
    scale_x_continuous(breaks = 0:10, limits = c(0,10)) +
    common_y_scale +
    labs(title = label, x = NULL, y = NULL) +
    common_theme
}

make_other_plot <- function(df, label) {
  df %>% 
    filter(Line2 == label) %>%
    ggplot(aes(x = GenNum, y = Hours, group = 1)) +
    geom_point(colour = "steelblue", size = 0.5) +
    geom_line(colour = "steelblue") +
    geom_errorbar(aes(ymin = Hours - se, ymax = Hours + se),width = 0.1, size = 0.5, alpha = 0.5) +
    scale_x_continuous(breaks = 21:31,limits = c(21, 31)) +
    common_y_scale +
    labs(title = label, x = NULL, y = NULL) +
    common_theme
}

make_SDolder_plot <- function(df, label) {
  df %>% 
    filter(Line2 == label) %>%
    ggplot(aes(x = GenNum, y = Hours, group = 1)) +
    geom_point(colour = "steelblue", size = 0.5) +
    geom_line(colour = "steelblue") +
    geom_errorbar(aes(ymin = Hours - se, ymax = Hours + se),width = 0.1, size = 0.5, alpha = 0.5) +
    scale_x_continuous(breaks = seq(120, 135, 5), limits = c(120, 135)) +
    common_y_scale +
    labs(title = label, x = NULL, y = NULL) +
    common_theme
}

# create plots (as before)
p_SDnew   <- make_first_col_plot(DF.massp1, "SDnew (42.88±1.76)")
p_CNnew   <- make_first_col_plot(DF.massp1, "CNnew (45.00±4.09)")
p_CTnew   <- make_first_col_plot(DF.massp1, "CTnew (37.72±1.74)")

p_SDold   <- make_other_plot(DF.massp1, "SDold (41.25±2.46)")
p_CNold   <- make_other_plot(DF.massp1, "CNold (38.00±3.18)")

p_SDolder <- make_SDolder_plot(DF.massp1, "SDolder (32.31±0.65)*")

# Stack columns (unchanged)
col1 <- p_SDnew / p_CNnew / p_CTnew + plot_layout(heights = rep(1,3))
col2 <- p_SDold / p_CNold / plot_spacer() + plot_layout(heights = rep(1,3))
col3 <- p_SDolder / plot_spacer() / plot_spacer() + plot_layout(heights = rep(1,3))


#create a clean empty placeholder plot for empty cells
empty_plot <- ggplot() + theme_void()

# layout (3 rows x 3 cols):# Row 1 (top): SDnew, SDold, SDolder, # Row 2 (mid): CNnew, CNold, empty, # Row 3 (bot): CTnew, empty,  empty
grid_plots <- list( p_SDnew,p_SDold,p_SDolder,p_CNnew,p_CNold,empty_plot, p_CTnew,empty_plot, empty_plot)

# combine into a regular grid with cowplot::plot_grid
# ncol = 3 makes a 3x3 grid and ensures each cell has the same width/height

main_grid <- plot_grid(plotlist = grid_plots, ncol = 3, align = "v", axis = "tblr")# align vertically, draw axes where present(top,bottom,left,right as appropriate)

# add a centered title using cowplot
title_row <- ggdraw() + draw_label("Mass bred lines", fontface = "bold", size = 8, x = 0.5, hjust = 0.5)

# stack title and grid
Figure1 <- plot_grid(title_row, main_grid, ncol = 1, rel_heights = c(0.06, 1))


# Annotate and save the figure by adding a common labels

FigureS1_annotated <- annotate_figure(Figure1,bottom = text_grob("Generation", color = "black",hjust = 0.5, face = "bold", size = 8),
                                     left = text_grob("Median survival time (hours)", color = "black",rot = 90,face = "bold", size = 8))

ggsave(plot = FigureS1_annotated, width = 7.0, height = 4.5, units = "in", dpi = 300,filename = "Figures/Figure S1 showing desiccation resistance in mass bred uncorrected for batch effects.jpg")


#_______________________________________________________________________________________________________



################################################################


#______________________________________________________________________________________________________

# Calculate the mean of median survival time for SDnew and CTnew compared vs all isolines as one group till G12

MassSDnew <- c("SDnew")
MassCTnew <- c("CTnew")
Isoline <- c("BR5","BR34","BR39","BR50","MB9","MB46","MB50","MB64", "NAB2","NAB11","NAB28","UT1","UT3", 
             "AS9","AS19","AS36","CT7","CT30","SYD6", "SYD13","SYD18","SYD53")

Massvsisolines <- DF %>%
  dplyr::filter(Line %in% c("SDnew", "CTnew", "BR5","BR34","BR39","BR50","MB9","MB46","MB50","MB64", 
                            "NAB2","NAB11","NAB28","UT1","UT3", "AS9","AS19","AS36","CT7","CT30","SYD6","SYD13","SYD18","SYD53")) %>% 
  dplyr::mutate(Type = ifelse(Line %in% MassSDnew, "MassSDnew", 
                              ifelse(Line %in% MassCTnew, "MassCTnew", "Isoline"))) %>% 
  dplyr::filter(Generation %in% c("0","1","2","3","4","5","6","7","8","9","10","11","12")) %>%
  summarySE(groupvars = c("Type"), measurevar = "rtime") # calculate mean of median knockdown time uncorrected for batch

# Calculate the mean of median survival time for SDnew and CTnew as one group compared vs all isolines as one group till G12

Mass <- c("SDnew","CTnew")

Massvsisolines01 <- DF %>%
  dplyr::filter(Line %in% c("SDnew", "CTnew", "BR5","BR34","BR39","BR50","MB9","MB46","MB50","MB64", 
                            "NAB2","NAB11","NAB28","UT1","UT3", "AS9","AS19","AS36","CT7","CT30","SYD6","SYD13","SYD18","SYD53")) %>% 
  dplyr::mutate(Type = ifelse(Line %in% Mass, "Mass",  "Isoline"))%>% 
  dplyr::filter(Generation %in% c("0","1","2","3","4","5","6","7","8","9","10","11","12")) %>%
  summarySE(groupvars = c("Type"), measurevar = "rtime") # calculate mean of median knockdown time uncorrected for batch



#_______________________________END_____________________

