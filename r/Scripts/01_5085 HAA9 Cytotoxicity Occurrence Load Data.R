## Scripts to support "Evaluating regulatory scenarios to limit U.S. nationwide exposure to cytotoxic HAAs"
## Peterson, Raseman, Stanford, Bruce, Klintworth, Reckhow, 2023, AWWA Water Science


# Load libraries
library(tidyverse)
# library(rpart)
# library(rpart.plot)
library(patchwork)
# library(ggbeeswarm)
library(scales)
# library(ggrepel)
library(RColorBrewer)

# Update ggplot theme
theme_set(
  theme_classic(base_size=8) %+replace% 
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.363636*2), # Add full border to plots, size is approximate of tick size
      axis.line = element_blank(), # Match axis line and panel border line
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black", size = rel(1.2)),
      plot.title = element_text(color = "black", size = rel(1.3), vjust = 1.5),
      axis.ticks = element_line(color = "black"),
      strip.background = element_blank(), # No facet label background
      legend.title = element_text(color = "black", size = rel(1)),
      legend.text = element_text(color = "black", size = rel(1)),
      strip.text = element_text(color = "black", size = rel(1), vjust = 1),
      plot.subtitle = element_text(color = "black", size = rel(1.2), vjust = 1)
    )
)


# source("C:/Users/epeterson/OneDrive - Hazen and Sawyer/Documents/R Projects/Theme_Set.R")
# theme_set(theme_EP_pub())

# Import DBP properties
index <- read_csv("Supporting Files/DBP Index.csv") %>%
  mutate(DBP = case_when(
    DBP == "CAA" ~ "MCAA",
    DBP == "BAA" ~ "MBAA",
    T ~ DBP
    )
  ) %>%
  mutate(subclass = fct_reorder(subclass, halogens)) %>%
  mutate(DBP = fct_reorder(DBP, order))

# Import DBP in vitro toxicity values
index.tox <- read_csv("Supporting Files/Toxicity Index Values.csv")

# Import PWSID treatment and water quality meta data
df.PWSID <- read_csv("Data/PWSID characteristics.csv", na = c("", "NA", 0))

# Import DBP occurrence data (LAA species concentrations (ug/L) at location of max LAA for each PWSID under HAA5, HAA6Br, and HAA9 scenarios)
# Note: each PWSID has 3 entries because the location of the max LAA can vary between HAA5, HAA6Br, and HAA9 regulatory scenarios.
df.DBP <- read_csv("Data/HAA Species Concentrations at Max LAA.csv", na = c("", "NA", 0)) %>%
  select(!c(FacilityID, SamplePointID))

# Merge PWSID characteristics and DBP occurrence data for analysis 
df <- full_join(df.PWSID, df.DBP, multiple = "all") %>%
  filter(!STG2_MCL == "Violation") %>%
  pivot_longer(c(BCAA:TCAA), names_to = "DBP", values_to = "c.mass") %>%
  left_join(., index, by = "DBP") %>%
  mutate(c.mol = c.mass/MW, .after = c.mass) %>%
  mutate(CAT = c.mol * CTI, .after = c.mol) %>%
  mutate(subclass = fct_reorder(subclass, halogens))
  

# Summarize data
df.sum <- df %>%
# df.sum_no.mHAA <- df %>%
# filter(!DBP %in% c("CAA", "BAA")) %>% #Exlude mHAAs option
  group_by(PWSID, STG2_MCL, Adv.Treat, Chloramines, Bromide_ppb, TOC_ppm, Impaired.Source, HAA.Group) %>%
  summarize(c.mass_HAA5 = sum(c.mass[group == "HAA5"], na.rm = T), #Compute groups
            c.mass_HAA6Br = sum(c.mass[!n.br == 0], na.rm = T),
            c.mass_HAA9 = sum(c.mass, na.rm = T),
            c.mol_HAA5 = sum(c.mol[group == "HAA5"], na.rm = T),
            c.mol_HAA6Br = sum(c.mol[!n.br == 0], na.rm = T),
            c.mol_HAA9 = sum(c.mol, na.rm = T),
            CAT_HAA5 = sum(CAT[group == "HAA5"], na.rm = T),
            CAT_HAA6Br = sum(CAT[!n.br == 0], na.rm = T),
            CAT_HAA9 = sum(CAT, na.rm = T),
            c.mass_MBAA = c.mass[DBP == "MBAA"],
            CAT_MBAA = CAT[DBP == "MBAA"]
            ) %>%
  ungroup() %>%
  mutate(CAT_high = ifelse(CAT_HAA9 >= quantile(CAT_HAA9, 0.95, na.rm = T), T, F)) %>% # Categorize CAT by percentile
  mutate(HAA6Br_MCL = ifelse(c.mass_HAA6Br >= 20, T, F)) %>% #Binary for hypothetical MCL
  mutate(HAA6Br_high = ifelse(c.mass_HAA6Br >= quantile(c.mass_HAA6Br, 0.95, na.rm = T), T, F)) %>% # Categorize HAA6Br by percentile
  mutate(HAA5_MCL = ifelse(c.mass_HAA5 >= 60, T, F)) %>% #Binary for existing MCL
  mutate(HAA9_MCL = ifelse(c.mass_HAA9 >= 60, T, F)) %>% #Binary for hypothetical MCL
  mutate(HAA9_high = ifelse(c.mass_HAA9 >= quantile(c.mass_HAA9, 0.95, na.rm = T), T, F)) %>% # Categorize HAA9 by percentile
  mutate(TOC.cat = case_when( # Define source water TOC bins
    is.na(TOC_ppm) ~ "NR",
    TOC_ppm < 1 ~ "< 1",
    TOC_ppm < 2 ~ "1-2",
    TOC_ppm < 3 ~ "2-3",
    TOC_ppm < 4 ~ "3-4",
    T ~ "> 4"
  )) %>%
  mutate(Br.cat = case_when( # Define source water bromide bins
    is.na(Bromide_ppb) ~ "NR",
    Bromide_ppb < 20 ~ "< 20",
    Bromide_ppb < 50 ~ "20-50",
    Bromide_ppb < 80 ~ "50-80",
    Bromide_ppb < 120 ~ "80-120",
    T ~ "> 120"
  )) %>%
  mutate(TOC.cat = factor(TOC.cat, levels = c("NR", "< 1", "1-2", "2-3", "3-4", "> 4"))) %>% # Arrange TOC bin factor levels
  mutate(Br.cat = factor(Br.cat, levels = c("NR", "< 20", "20-50", "50-80", "80-120", "> 120"))) # Arrange bromide bin factor levels
  

# Summarize data excluding mHAA concentrations
df.sum_no.mHAA <- df %>%
  filter(!DBP %in% c("MCAA", "MBAA")) %>% 
  group_by(PWSID, STG2_MCL, Adv.Treat, Chloramines, Bromide_ppb, TOC_ppm, Impaired.Source, HAA.Group) %>%
  summarize(c.mass_HAA5 = sum(c.mass[group == "HAA5"], na.rm = T), #Compute groups
            c.mass_HAA6Br = sum(c.mass[!n.br == 0], na.rm = T),
            c.mass_HAA9 = sum(c.mass, na.rm = T),
            c.mol_HAA5 = sum(c.mol[group == "HAA5"], na.rm = T),
            c.mol_HAA6Br = sum(c.mol[!n.br == 0], na.rm = T),
            c.mol_HAA9 = sum(c.mol, na.rm = T),
            CAT_HAA5 = sum(CAT[group == "HAA5"], na.rm = T),
            CAT_HAA6Br = sum(CAT[!n.br == 0], na.rm = T),
            CAT_HAA9 = sum(CAT, na.rm = T)
  ) %>%
  ungroup() %>%
  mutate(CAT_high = ifelse(CAT_HAA9 >= quantile(CAT_HAA9, 0.95, na.rm = T), T, F)) %>% # Categorize CAT by percentile
  mutate(HAA6Br_MCL = ifelse(c.mass_HAA6Br >= 20, T, F)) %>% #Binary for hypothetical MCL
  mutate(HAA6Br_high = ifelse(c.mass_HAA6Br >= quantile(c.mass_HAA6Br, 0.95, na.rm = T), T, F)) %>% # Categorize HAA6Br by percentile
  mutate(HAA5_MCL = ifelse(c.mass_HAA5 >= 60, T, F)) %>% #Binary for existing MCL
  mutate(HAA9_MCL = ifelse(c.mass_HAA9 >= 60, T, F)) %>% #Binary for hypothetical MCL
  mutate(HAA9_high = ifelse(c.mass_HAA9 >= quantile(c.mass_HAA9, 0.95, na.rm = T), T, F)) %>% # Categorize HAA6Br by percentile
  mutate(TOC.cat = case_when( # Define source water TOC bins
    is.na(TOC_ppm) ~ "NR",
    TOC_ppm < 1 ~ "< 1",
    TOC_ppm < 2 ~ "1-2",
    TOC_ppm < 3 ~ "2-3",
    TOC_ppm < 4 ~ "3-4",
    T ~ "> 4"
  )) %>%
  mutate(Br.cat = case_when( # Define source water bromide bins
    is.na(Bromide_ppb) ~ "NR",
    Bromide_ppb < 20 ~ "< 20",
    Bromide_ppb < 50 ~ "20-50",
    Bromide_ppb < 80 ~ "50-80",
    Bromide_ppb < 120 ~ "80-120",
    T ~ "> 120"
  )) %>%
  mutate(TOC.cat = factor(TOC.cat, levels = c("NR", "< 1", "1-2", "2-3", "3-4", "> 4"))) %>% # Arrange TOC bin factor levels
  mutate(Br.cat = factor(Br.cat, levels = c("NR", "< 20", "20-50", "50-80", "80-120", "> 120"))) # Arrange bromide bin factor levels
