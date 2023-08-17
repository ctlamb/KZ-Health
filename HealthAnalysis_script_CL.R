# Analyzing Klinse-Za health samples # 

## load packages ##
library(here)
library(rstatix)
library(hrbrthemes)
library(gghalves)
library(lubridate)
library(lme4)
library(readxl)
library(ggeffects)
library(broom.mixed)
library(MuMIn)
library(sjPlot)
library(ggpubr)
library(FSA)
library(corrplot)
library(scales)
library(tidyverse)
library(tidylog)

options(scipen = 999)
# Assumptions
### GLM
# Linearity: The relationship between X and the mean of Y is linear.
# Homoscedasticity: The variance of residual is the same for any value of X.
# Independence: Observations are independent of each other.
# Normality: For any fixed value of X, Y is normally distributed..errors are normal

### t test
# The data are continuous.
# The sample data have been randomly sampled from a population.
# There is homogeneity of variance (i.e., the variability of the data in each group is similar).
# The distribution is approximately normal.

## Kruskal Wallis
# observations in the data set are independent of each other.
# distribution of the population should not be necessarily normal and the variances should not be necessarily equal.
# observations must be drawn from the population by the process of random sampling.

#*********************** #
### LOAD DATA #####
#*********************** #

## read in entire annual capture results sheet
AnnualCaptData <- read_excel("data/KZONCP_Health_Year3_results_221214.xlsx", sheet = "annual_capture_results", na = c("unk", "Unk", "wild", "NA", "na", "NRH", "NotSentYet", "pending", "no result", "No sample"))

## cap locs
cap.locs <- read_excel("data/Copy of KZCaptureByZone_SM_220311.xlsx", sheet = "4ClaytonNoNulls")

# fecal cort
FecalCort <- read_excel("data/KZ_FecalCort_PenWild.xlsx", na = c("unk", "Unk", "wild", "NA", "na", "NRH", "NotSentYet", "pending"))

# fecal cort
FecalCort.id <- read_excel("data/genetic_samples_identified.xlsx", sheet = "success_fecal_samples")

# fecal nit
fecal.nit <- read_excel("data/fecal_nitrogen_data.xlsx", na = c("unk", "Unk", "wild", "NA", "na", "NRH", "NotSentYet", "pending")) %>%
  filter(study_area %in% "Klinse-Za")

# other herds cort and nutrients
cort.other <- read_excel("data/ONCP_HCC_WithCollectionAttributes_220704.xlsx", sheet = "Combined", na = c("unk", "Unk", "wild", "NA", "na", "NRH", "NotSentYet", "pending"))

nutr.other <- read_excel("data/ONCP_health_RawLabData_220315_CL.xlsx", sheet = "Trace_nutrients")

bor <- read_excel("data/Caribou_data_trace_minerals.xlsx", na = c("NA", "No Sample"))

# other herds Pathogens
pathogens.other <- read_excel("data/ONCP_health_RawLabData_220315_CL.xlsx", sheet = "Erysipelothrix") %>%
  mutate(name = "Erysipelothrix") %>%
  select(Wii_ID = `Animal ID`, Herd = `Study Area`, Sex, name, value = Interpretation) %>%
  rbind(
    read_excel("data/ONCP_health_RawLabData_220315_CL.xlsx", sheet = "Neospora") %>%
      mutate(name = "Neospora") %>%
      select(Wii_ID = `WII Animal ID`, Herd = `Study Area`, Sex, name, value = Result)
  ) %>%
  rbind(
    read_excel("data/ONCP_health_RawLabData_220315_CL.xlsx", sheet = "Herpes_IBR") %>%
      mutate(name = "Herpes_IBR") %>%
      select(Wii_ID = `WII Animal ID`, Herd = `Study Area`, Sex, name, value = Result)
  )


#*********************** #
### PREP DATA ####
#*********************** #

## clean up columns and derive some new variables
individual_data <- AnnualCaptData %>%
  dplyr::select(
    WLH_ID = `WLH ID`,
    Wii_ID = `WII Animal ID`,
    Date,
    Year,
    herd = herd_area,
    Calf_result,
    stay_in_pen,
    weight = `weight (kg)`,
    condition_clean,
    body_fat_perc,
    age_clean = est_age_clean,
    Toxoplasma:`Haptoglobin g/L`,
    -ERYSIP,
    hair_cort
  ) %>%
  mutate(calf = case_when(
    Calf_result %in% "Alive" ~ 1,
    is.na(Calf_result) ~ NA_real_,
    TRUE ~ 0
  )) %>%
  mutate(
    calf_succ = case_when(
      calf == 1 ~ "Yes",
      calf == 0 ~ "No"
    ),
    calf_3 = case_when(
      Calf_result %in% "np" ~ "not preg.",
      Calf_result %in% "Alive" ~ "alive",
      Calf_result %in% "Neonatal death" ~ "alive",
      Calf_result %in% "Aborted/Stillborn" ~ "abort",
      TRUE ~ Calf_result
    ),
    preg = case_when(
      calf_3 %in% "not preg." ~ "not_preg",
      calf_3 %in% c("alive", "abort") ~ "preg",
      is.na(calf_3) ~ NA_character_
    ),
    loc = case_when(
      is.na(stay_in_pen) ~ "Free-ranging",
      TRUE ~ "Pen"
    )
  ) %>%
  mutate(
    calf_3 = fct_relevel(calf_3, "not preg.", "abort", "alive"),
    age_clean = fct_relevel(age_clean, "young", "mature", "old"),
    condition_clean = fct_relevel(condition_clean, "poor", "fair", "good")
  )

## add in capture location
individual_data <- individual_data %>%
  left_join(cap.locs %>%
    select(
      Wii_ID = `Animal_Id`,
      Year,
      cap.loc = SURVEYUNIT
    ) %>%
    distinct())

# Mo has some <MDL values, min detection limit, which is 0.9. Used midpoint between mdl and 0, so 0.45
### need to find out MDL
individual_data <- individual_data %>%
  mutate(`Mo (ng/mL)` = case_when(
    `Mo (ng/mL)` %in% "<MDL" ~ 0.45,
    TRUE ~ as.numeric(`Mo (ng/mL)`)
  ))

nutr.other <- nutr.other %>%
  mutate(
    `Mo (ng/mL)` = case_when(
      `Mo (ng/mL)` %in% "<MDL" ~ 0.45,
      TRUE ~ as.numeric(`Mo (ng/mL)`)
    ),
    `Co (ng/mL)` = case_when(
      `Co (ng/mL)` %in% "<MDL" ~ 0.2,
      TRUE ~ as.numeric(`Co (ng/mL)`)
    )
  )

bor <- bor %>%
  mutate(
    `Molybden ppm` = case_when(
      `Molybden ppm` %in% "BDL" ~ (0.45 / 1000),
      TRUE ~ as.numeric(`Molybden ppm`)
    ),
    `HCC (pg/mg)` = case_when(
      `HCC (pg/mg)` %in% "BDL" ~ (0.2),
      TRUE ~ as.numeric(`HCC (pg/mg)`)
    )
  )

## some 2019 animals caught twice, only 1 set of health data taken
individual_data <- individual_data %>%
  filter(month(Date) %in% 3:4) %>%
  group_by(Wii_ID, Year) %>%
  add_count(name = "records.yr") %>%
  ungroup()

## summarise sample sizes
individual_data %>%
  # drop_na(`Se (ug/mL)`)%>%
  distinct(Wii_ID, Year, .keep_all = TRUE) %>%
  group_by(loc) %>%
  summarise(
    individuals = n_distinct(Wii_ID),
    an_yrs = n()
  )

## denote location the year before
individual_data <- individual_data %>%
  left_join(
    individual_data %>%
      select(Wii_ID, Year, loc.lastyr = loc) %>%
      mutate(Year = Year + 1),
    by = c("Wii_ID", "Year")
  ) %>%
  mutate(loc.lastyr = replace_na(loc.lastyr, "Free-ranging"))

## unique individuals
unique(individual_data$Wii_ID) %>% length()

## unique individuals that spent time in the pen
individual_data %>%
  drop_na(stay_in_pen) %>%
  group_by(Wii_ID) %>%
  filter(stay_in_pen == max(stay_in_pen)) %>%
  ungroup() %>%
  summarise(
    n = n_distinct(Wii_ID),
    n.stays = sum(stay_in_pen),
    max.stays = max(stay_in_pen)
  )

## assess repeat detections of same animal for RE comment during R1
individual_data%>%
  distinct(Wii_ID,Year)%>%
  group_by(Wii_ID)%>%
  summarise(n=n())%>%
  group_by(n)%>%
  summarise(count=n())

individual_data%>%
  distinct(Wii_ID,Year)%>%
  group_by(Wii_ID)%>%
  summarise(n=n())%>%
  ggplot(aes(x=n))+
  geom_histogram()


### Correlation amongst all health metrics
M <- cor(
  individual_data %>%
    filter(`Haptoglobin g/L` < 1, hair_cort < 200) %>%
    select(Neospora:hair_cort, -Negative_all3) %>%
    rename(Aphaherpesvirus = AHC_IBR, Eryysipelothrix = Erysip_simplified, `Hair cortisol` = hair_cort) %>%
    mutate(across(where(is.character),
      .fns = ~ case_when(
        . %in% c("Positive", "positive") ~ 1,
        . %in% c("Negative", "negative") ~ 0,
        TRUE ~ NA_real_
      )
    )),
  use = "complete.obs"
)

order.df <- tibble(
  name = c(
    "Neospora", "Aphaherpesvirus", "Eryysipelothrix", "Mn (ng/mL)", "Fe (ug/mL)", "Co (ng/mL)", "Cu (ug/mL)", "Zn (ug/mL)",
    "Se (ug/mL)", "Mo (ng/mL)", "Haptoglobin g/L", "Hair cortisol"
  ),
  color = c(
    rep("orchid", times = 3),
    rep("black", times = 7),
    "forestgreen",
    "tan2"
  )
)

pdf(file = here::here("output", "plots", "correlations.pdf"))
corrplot(M, type = "upper", p.mat = 1 - abs(M), sig.level = 0.9, col = COL2("BrBG"), order = "original", tl.col = order.df$color, insig = "blank", addCoef.col = "black")
dev.off()

## mean and SD of times in pen
individual_data %>%
  drop_na(stay_in_pen) %>%
  group_by(Wii_ID) %>%
  filter(stay_in_pen == max(stay_in_pen)) %>%
  ungroup() %>%
  summarize(
    mean = mean(stay_in_pen),
    sd = sd(stay_in_pen),
    se = sd / sqrt(n())
  )

#*********************** #
### TRACE MINERALS #####
#*********************** #

## pivot to longer format
CalvingNutri.long <- individual_data %>%
  select(Wii_ID, Year, age_clean, calf_succ, calf_3, preg, age_clean, stay_in_pen, loc, loc.lastyr, `Mn (ng/mL)`:`Mo (ng/mL)`, cap.loc) %>%
  drop_na(calf_succ) %>%
  pivot_longer(`Mn (ng/mL)`:`Mo (ng/mL)`)


## compare to other herds
nutr.plot <- nutr.other %>%
  select(Wii_ID = `WII Animal ID`, Sex, Herd = `Study Area`, `Mn (ng/mL)`:`Mo (ng/mL)`) %>%
  pivot_longer(`Mn (ng/mL)`:`Mo (ng/mL)`) %>%
  rbind(CalvingNutri.long %>%
    mutate(
      Sex = "F",
      Herd = case_when(
        loc == "Pen" ~ "Klinse-Za",
        loc == "Free-ranging" ~ "Klinse-Za"
      )
    ) %>%
    select(Wii_ID, Sex, Herd, name, value)) %>%
  mutate(Herd = relevel(factor(Herd), ref = "Klinse-Za")) %>%
  rbind(bor %>%
    mutate(
      `Mn (ng/mL)` = `Manganese ppm` %>% as.numeric() * 1000,
      `Fe (ug/mL)` = `Iron ppm` %>% as.numeric(),
      `Co (ng/mL)` = `Cobalt ppb` %>% as.numeric(),
      `Cu (ug/mL)` = `Copper ppm` %>% as.numeric(),
      `Zn (ug/mL)` = `Zinc ppm` %>% as.numeric(),
      `Se (ug/mL)` = `Selenium ppm` %>% as.numeric(),
      `Mo (ng/mL)` = `Molybden ppm` * 1000,
      Herd = "BC Boreal"
    ) %>%
    select(Wii_ID = `Caribou ID`, Sex, Herd, `Mn (ng/mL)`:`Mo (ng/mL)`) %>%
    pivot_longer(`Mn (ng/mL)`:`Mo (ng/mL)`)) %>%
  mutate(Herd = relevel(factor(Herd), ref = "Klinse-Za")) %>%
  filter(Sex == "F") ## 8 records

## clean up data
nutr.plot2 <- nutr.plot %>%
  filter(!Herd %in% "Scott") %>% ## only 3 animals
  mutate(cull = case_when(
    name %in% c("Fe (ug/mL)", "Mn (ng/mL)") & value > 100 ~ 1,
    name %in% c("Se (ug/mL)") & value > 0.3 ~ 1,
    name %in% c("Mo (ng/mL)") & value > 30 ~ 1,
    TRUE ~ 0
  )) %>%
  filter(!cull == 1) ## 9 rows removed of 3,645



nutr.compare <- ggplot(
  data = nutr.plot2,
  aes(Herd, value, col = Herd == "Klinse-Za")
) +
  geom_boxplot(position = position_dodge(0.3), width = .6, outlier.shape = NA, alpha = 0.4, show_guide = FALSE) + ## plot first to get order right.
  geom_rect(
    data = tibble(
      name = unique(nutr.plot %>% pull(name)),
      xmin = "Klinse-Za", xmax = "BC Boreal",
      Herd = "Wolverine",
      value = 0,
      ymin = c(NA, NA, NA, 0.7, 1.1, 0.05, NA), ymax = c(NA, NA, NA, 1.8, 2.5, 0.14, NA)
    ),
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    color = NA,
    alpha = .2
  ) +
  geom_boxplot(position = position_dodge(0.3), width = .6, outlier.shape = NA, alpha = 0.4, show_guide = FALSE) +
  geom_jitter(width = 0.1, alpha = .1, show_guide = FALSE) +
  theme_ipsum() +
  labs(x = "", y = "Mineral level") +
  theme(
    axis.title.x = element_text(size = 18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = element_text(size = 11, angle = 35, hjust = 1),
    axis.text.y = element_text(size = 13),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    strip.text.x.top = element_text(size = 13),
    legend.position = "none"
  ) +
  scale_y_continuous(breaks=breaks_pretty(n=3))+
  guides(fill = guide_legend(title = "Location"), color = guide_legend(title = "Location")) +
  facet_wrap(vars(name), scales = "free_y")

nutr.comp.betas <- nutr.plot %>%
  filter(Herd != "Scott", value < 1000) %>%
  group_by(name) %>%
  dunn_test(value ~ Herd) %>%
  filter(group1 %in% "Klinse-Za") %>%
  mutate(p.val = case_when(
    p.adj < 0.05 ~ "*",
    TRUE ~ NA_character_
  )) %>%
  ggplot(aes(x = statistic, y = name, group = group2, color = group2, label = p.val)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(position = position_dodge(width = -0.5), size = 2) +
  geom_text(position = position_dodge(width = -0.5), vjust = 0.1) +
  labs(x = "H Statistic (relative to Klinse-Za)", y = "Term", color = "") +
  theme_ipsum() +
  theme(
    axis.title.x = element_text(size = 18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE))

## compare between all KZ vs first cap only KZ
nutr.comparepooled <- nutr.plot %>%
  filter(value < 1000, Herd=="Klinse-Za") %>%
  group_by(name)%>%
  summarise(median=median(value),
            sd=sd(value),
            lcl=case_when(median-(1.96*sd)>0~median-(1.96*sd),
                           TRUE~0),
            ucl=median+(1.96*sd))%>%
  mutate(pooled=paste0(median%>%round(2),
                       " (",
                       lcl%>%round(2),
                       "-",
                       ucl%>%round(2),
                       ")"))%>%
  select(type=name, pooled)%>%
  left_join(nutr.plot %>%
          filter(value < 1000, Herd=="Klinse-Za") %>%
          distinct(Wii_ID, name, .keep_all = TRUE)%>%
          group_by(name)%>%
          summarise(median=median(value),
                    sd=sd(value),
                    lcl=case_when(median-(1.96*sd)>0~median-(1.96*sd),
                                  TRUE~0),
                    ucl=median+(1.96*sd))%>%
            mutate(single=paste0(median%>%round(2),
                                 " (",
                                 lcl%>%round(2),
                                 "-",
                                 ucl%>%round(2),
                                 ")"))%>%
            select(type=name, single),
          by="type")


## stats if KZ nutrients are different than other herds

## compare KZ to ONCP as a whole.
nutr.plot %>%
  filter(!Herd %in% c("Scott", "BC Boreal"), value < 1000) %>%
  mutate(Herd = case_when(
    Herd == "Klinse-Za" ~ "KZ",
    TRUE ~ "ONCP"
  )) %>%
  group_by(name) %>%
  kruskal_test(value ~ Herd) %>%
  arrange(p) %>%
  filter(p < 0.05)




## plot together
ggarrange(nutr.compare, nutr.comp.betas,
  ncol = 2, nrow = 1,
  labels = "AUTO",
  widths = c(1.5, 1)
) %>%
  annotate_figure(top = text_grob("Trace Mineral Comparisons", color = "black", face = "bold", size = 14))

ggsave(here::here("output", "plots", "mineral_compare_plate.png"), width = 12, height = 8, bg = "white")



# Reproductive outcome all classes
nutr.repro <- ggplot(CalvingNutri.long %>%
  drop_na(stay_in_pen, value, calf_3), aes(x = calf_3, y = value, fill = name, color = name)) +
  # ggdist::stat_halfeye( width =1, .width = 0, alpha=0.6,justification = -.2, point_colour = NA)+
  geom_boxplot(position = position_dodge(0.3), width = .3, alpha = 0.4, show_guide = FALSE) +
  # geom_point(position=position_dodge(0.3),alpha = .2, show_guide=FALSE)+
  labs(y = "Mineral level", x = "Reproductive outcome") +
  theme_ipsum() +
  theme(
    axis.title.x = element_text(size = 18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = element_text(size = 13, angle = 35, hjust = 1),
    axis.text.y = element_text(size = 13),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.position = "none"
  ) +
  scale_y_continuous(breaks = breaks_pretty(n = 3)) +
  facet_wrap(vars(name), scales = "free_y", ncol = 4)

# ggsave(here::here("output","plots","Repro.nutrient.png"),width=10,height=5, bg="white")


## statistics for pen only

# Calf success outcome binary
CalvingNutri.long %>%
  drop_na(stay_in_pen, value, calf_succ) %>%
  group_by(name) %>%
  kruskal_test(value ~ calf_succ) %>%
  arrange(p)

# Pregnancy
CalvingNutri.long %>%
  drop_na(stay_in_pen, value, preg) %>%
  group_by(name) %>%
  kruskal_test(value ~ preg) %>%
  arrange(p)

# Reproductive outcome all classes
CalvingNutri.long %>%
  drop_na(stay_in_pen, value, calf_3) %>%
  group_by(name) %>%
  dunn_test(value ~ calf_3) %>%
  arrange(p)



### plot nutrient levels as a function of number of stints in the pen
ggplot(CalvingNutri.long %>% drop_na(stay_in_pen), aes(x = stay_in_pen, y = value)) +
  geom_point(aes(x = stay_in_pen, y = value, fill = Wii_ID, color = Wii_ID, group = Wii_ID),
    alpha = 0.5
  ) +
  geom_path(aes(x = stay_in_pen, y = value, fill = Wii_ID, color = Wii_ID, group = Wii_ID),
    alpha = 0.5
  ) +
  geom_smooth(aes(x = stay_in_pen, y = value),
    method = "lm",
    color = "black",
    se = FALSE,
    linetype = "dashed"
  ) +
  labs(x = "Stays in pen (n)", y = "Nutrient level") +
  theme_ipsum() +
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.position = "none"
  ) +
  facet_wrap(vars(name), scales = "free_y")

# ggsave(here::here("output","plots","pentime.nutrient.png"),width=7,height=7, bg="white")

## check annual variaton
ggplot(CalvingNutri.long, aes(x = Year, y = value)) +
  geom_point(aes(x = Year, y = value, fill = Wii_ID, color = Wii_ID, group = Wii_ID),
    alpha = 0.5
  ) +
  geom_path(aes(x = Year, y = value, fill = Wii_ID, color = Wii_ID, group = Wii_ID),
    alpha = 0.5
  ) +
  geom_smooth(aes(x = Year, y = value),
    color = "black",
    se = FALSE,
    linetype = "dashed"
  ) +
  labs(x = "Year", y = "Nutrient level") +
  theme_ipsum() +
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.position = "none"
  ) +
  facet_wrap(vars(name), scales = "free_y")

## annual trend, penned vs unpenned
nutr.time <- ggplot(CalvingNutri.long, aes(x = Year, y = value, fill = loc.lastyr, color = loc.lastyr)) +
  geom_point(alpha = 0.5) +
  geom_path(aes(group = Wii_ID), alpha = 0.5) +
  geom_smooth(aes(x = Year, y = value),
    method = "lm",
    color = "grey40",
    se = FALSE,
    linetype = "dashed"
  ) +
  labs(x = "Year", y = "Mineral level") +
  theme_ipsum() +
  theme(
    axis.title.x = element_text(size = 18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    strip.text.y.right = element_text(size = 10),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.position = "none"
  ) +
  scale_x_continuous(breaks = breaks_pretty(n = 2)) +
  scale_y_continuous(breaks = breaks_pretty(n = 2)) +
  facet_grid(name ~ loc.lastyr, scales = "free_y")

# ggsave(here::here("output","plots","annual.nutrient.png"),width=5,height=10, bg="white")




## statistics

## need RE for year?
m1 <- lmer(value ~ stay_in_pen + (1 | Year), data = CalvingNutri.long %>% drop_na(stay_in_pen, value, age_clean) %>% filter(name == "Zn (ug/mL)"))
m1a <- lm(value ~ stay_in_pen, data = CalvingNutri.long %>% drop_na(stay_in_pen, value, age_clean) %>% filter(name == "Zn (ug/mL)"))
anova(m1, m1a)

m1 <- lmer(value ~ stay_in_pen + (1 | Year), data = CalvingNutri.long %>% drop_na(stay_in_pen, value, age_clean) %>% filter(name == "Se (ug/mL)"))
m1a <- lm(value ~ stay_in_pen, data = CalvingNutri.long %>% drop_na(stay_in_pen, value, age_clean) %>% filter(name == "Se (ug/mL)"))
anova(m1, m1a)

m1 <- lmer(value ~ stay_in_pen + (1 | Year), data = CalvingNutri.long %>% drop_na(stay_in_pen, value, age_clean) %>% filter(name == "Mn (ng/mL)"))
m1a <- lm(value ~ stay_in_pen, data = CalvingNutri.long %>% drop_na(stay_in_pen, value, age_clean) %>% filter(name == "Mn (ng/mL)"))
anova(m1, m1a)

## Evidence for RE for year. Use RE

## need RE for cap.loc?
m1 <- lmer(value ~ stay_in_pen + (1 | cap.loc), data = CalvingNutri.long %>% drop_na(stay_in_pen, value, age_clean, cap.loc) %>% filter(name == "Zn (ug/mL)"))
m1a <- lm(value ~ stay_in_pen, data = CalvingNutri.long %>% drop_na(stay_in_pen, value, age_clean, cap.loc) %>% filter(name == "Zn (ug/mL)"))
anova(m1, m1a)

m1 <- lmer(value ~ stay_in_pen + (1 | cap.loc), data = CalvingNutri.long %>% drop_na(stay_in_pen, value, age_clean, cap.loc) %>% filter(name == "Se (ug/mL)"))
m1a <- lm(value ~ stay_in_pen, data = CalvingNutri.long %>% drop_na(stay_in_pen, value, age_clean, cap.loc) %>% filter(name == "Se (ug/mL)"))
anova(m1, m1a)

m1 <- lmer(value ~ stay_in_pen + (1 | cap.loc), data = CalvingNutri.long %>% drop_na(stay_in_pen, value, age_clean, cap.loc) %>% filter(name == "Mn (ng/mL)"))
m1a <- lm(value ~ stay_in_pen, data = CalvingNutri.long %>% drop_na(stay_in_pen, value, age_clean, cap.loc) %>% filter(name == "Mn (ng/mL)"))
anova(m1, m1a)

## No Evidence for RE for cap loc

## need RE for individual?
m1 <- lmer(value ~ stay_in_pen + (1 | Year), data = CalvingNutri.long %>% drop_na(stay_in_pen, value, age_clean) %>% filter(name == "Zn (ug/mL)"))
m1a <- lmer(value ~ stay_in_pen + (1 | Year) + (1 | Wii_ID), data = CalvingNutri.long %>% drop_na(stay_in_pen, value, age_clean) %>% filter(name == "Zn (ug/mL)"))
anova(m1, m1a)

m1 <- lmer(value ~ stay_in_pen + (1 | Year), data = CalvingNutri.long %>% drop_na(stay_in_pen, value, age_clean) %>% filter(name == "Se (ug/mL)"))
m1a <- lmer(value ~ stay_in_pen + (1 | Year) + (1 | Wii_ID), data = CalvingNutri.long %>% drop_na(stay_in_pen, value, age_clean) %>% filter(name == "Se (ug/mL)"))
anova(m1, m1a)

m1 <- lmer(value ~ stay_in_pen + (1 | Year), data = CalvingNutri.long %>% drop_na(stay_in_pen, value, age_clean) %>% filter(name == "Mn (ng/mL)"))
m1a <- lmer(value ~ stay_in_pen + (1 | Year) + (1 | Wii_ID), data = CalvingNutri.long %>% drop_na(stay_in_pen, value, age_clean) %>% filter(name == "Mn (ng/mL)"))
anova(m1, m1a)

## No Evidence for RE for individual

## check stay in pen stats
## controlling for age
CalvingNutri.long %>%
  drop_na(stay_in_pen, age_clean, value) %>%
  group_by(name) %>%
  do(tidy(lmer(value ~ stay_in_pen + age_clean + calf_succ + (1 | Year), data = .))) %>%
  filter(term %in% "stay_in_pen") %>%
  arrange(-abs(statistic)) %>%
  mutate(sig = case_when(abs(statistic) >= 1.96 ~ "*", TRUE ~ "")) ## Se higher


## change through time
CalvingNutri.long %>%
  drop_na(age_clean, value) %>%
  group_by(name) %>%
  do(tidy(lm(value ~ Year + loc.lastyr + age_clean + calf_succ, data = .))) %>%
  filter(term != "(Intercept)") %>%
  arrange(-abs(statistic)) %>%
  mutate(sig = case_when(abs(statistic) >= 1.96 ~ "*", TRUE ~ "")) ## Co and Zn declining through time, Fe increasing through time, Co lower in old animals, Se higher in animals in pen last yr


## correlation between stays in pen and year

# across individuals
CalvingNutri.long %>%
  drop_na(stay_in_pen) %>%
  select(Year, stay_in_pen) %>%
  cor() ## 0.46, not bad

## by individual
CalvingNutri.long %>%
  drop_na(stay_in_pen) %>%
  select(Wii_ID, Year, stay_in_pen) %>%
  group_by(Wii_ID) %>%
  summarise(cor = cor(Year, stay_in_pen)) %>%
  drop_na() %>%
  ungroup() %>%
  summarise(mean = mean(cor)) ## higly correlated

## change through time accounting for stays in pen, only for penned animals
CalvingNutri.long %>%
  drop_na(age_clean, value, stay_in_pen) %>%
  group_by(name) %>%
  do(tidy(lm(value ~ Year + age_clean + calf_succ + stay_in_pen, data = .))) %>%
  filter(term != "(Intercept)") %>%
  arrange(-abs(statistic)) %>%
  mutate(sig = case_when(abs(statistic) >= 1.96 ~ "*", TRUE ~ "")) ## Only thing changing with stays in pen is higher selenium


## plot change through time betas with standardized betas
nutr.beta <- CalvingNutri.long %>%
  drop_na(age_clean, value) %>%
  group_by(name) %>%
  mutate(across(c("Year", "value"), scale)) %>%
  do(tidy(lm(value ~ Year + loc.lastyr + age_clean + calf_succ, data = .))) %>%
  filter(term != "(Intercept)") %>%
  mutate(term = case_when(
    term == "loc.lastyrPen" ~ "Location prev. yr (pen:unpenned)",
    term == "calf_succYes" ~ "Calf (yes:no)",
    term == "age_cleanold" ~ "Ageclass (old:young)",
    term == "age_cleanmature" ~ "Ageclass (mature:young)",
    TRUE ~ term
  )) %>%
  arrange(-abs(statistic)) %>%
  mutate(sig = case_when(abs(statistic) >= 1.96 ~ "*", TRUE ~ "")) %>%
  ggplot(aes(x = estimate, xmin = estimate - (std.error * 1.96), xmax = estimate + (std.error * 1.96), y = term, group = name, color = name)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_linerange(position = position_dodge(width = -0.5)) +
  geom_point(position = position_dodge(width = -0.5)) +
  labs(x = "Estimate (standardized)", y = "Term", color = "Mineral") +
  theme_ipsum() +
  theme(
    axis.title.x = element_text(size = 18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.position = "right"
  )
# ggsave(here::here("output","plots","nutrient.betas.png"),width=8,height=6, bg="white")


## check interaction
CalvingNutri.long %>%
  drop_na(age_clean, value) %>%
  mutate(Year = Year / 2000) %>%
  group_by(name) %>%
  do(tidy(lm(value ~ Year * loc.lastyr + age_clean + calf_succ, data = .))) %>%
  filter(term != "(Intercept)") %>%
  arrange(-abs(statistic)) %>%
  mutate(sig = case_when(abs(statistic) >= 1.96 ~ "*", TRUE ~ "")) ## no evidence for that changes through time are different between penned and unpenned

## check interaction again using individual models
m1 <- lm(value ~ Year + loc.lastyr + age_clean + calf_succ, data = CalvingNutri.long %>% filter(name == "Zn (ug/mL)"))
m2 <- lm(value ~ Year * loc.lastyr + age_clean + calf_succ, data = CalvingNutri.long %>% filter(name == "Zn (ug/mL)"))
anova(m1, m2) ## no effect of interaction

plot(ggpredict(m1, terms = c("Year", "loc.lastyr"))) ## Zn declining in both groups, Zn lower in pen, but binned age may be an effect
plot_model(m1, "diag")
summary(m1)


m1 <- lm(value ~ Year + loc.lastyr + age_clean + calf_succ, data = CalvingNutri.long %>% filter(name == "Co (ng/mL)"))
m2 <- lm(value ~ Year * loc.lastyr + age_clean + calf_succ, data = CalvingNutri.long %>% filter(name == "Co (ng/mL)"))
anova(m1, m2) ## no effect of interaction

plot(ggpredict(m1, terms = c("Year", "loc.lastyr"))) ## Co declining in both groups, Co lower in pen, but binned age may be an effect
plot_model(m1, "diag")
summary(m1)




## does capture location matter?
CalvingNutri.long %>%
  drop_na(age_clean, value) %>%
  mutate(Year = Year / 2000) %>%
  group_by(name) %>%
  do(tidy(lm(value ~ Year + loc.lastyr + age_clean + calf_succ + cap.loc, data = .))) %>%
  filter(str_detect(term, "cap.loc")) %>%
  arrange(-abs(statistic)) %>%
  mutate(sig = case_when(abs(statistic) >= 1.96 ~ "*", TRUE ~ "")) ## most RE's not significant and some high, some low



## plot together
ggarrange(ggarrange(nutr.repro, nutr.beta, ncol = 1, nrow = 2, labels = c("A", "C")), nutr.time,
  ncol = 2, nrow = 1,
  labels = c(NA, "B"),
  widths = c(1.5, 1)
) %>%
  annotate_figure(top = text_grob("Trace Minerals", color = "black", face = "bold", size = 14))

ggsave(here::here("output", "plots", "mineral_plate.png"), width = 11, height = 9, bg = "white")


#### Take home####
## Co and Zn declining through time, Fe, Cu, and Zn lower in pen, Fe, Co, and Zn lower in old animals.
## Se increasing with increased penning. All other temporal effects are happening irrespective of capture location or penning.
## looked at   + cap.loc and (1|cap.loc) but didn't find much. slightly higher Co in Rochfort and lower Cu (**)


#*********************** #
### HAPTAGLOBIN #####
#*********************** #

haptoglobin <- individual_data %>%
  select(Wii_ID, calf_succ, calf_3, preg, age_clean, stay_in_pen, `Haptoglobin g/L`, Year, loc, loc.lastyr) %>%
  drop_na(calf_succ)

## mean

hapto.comparepooled <- haptoglobin %>%
  drop_na(`Haptoglobin g/L`)%>%
  summarise(median=median(`Haptoglobin g/L`),
            sd=sd(`Haptoglobin g/L`),
            lcl=case_when(median-(1.96*sd)>0~median-(1.96*sd),
                          TRUE~0),
            ucl=median+(1.96*sd))%>%
  mutate(type="Haptoglobin g/L",
         pooled=paste0(median%>%round(2),
                       " (",
                       lcl%>%round(2),
                       "-",
                       ucl%>%round(2),
                       ")"))%>%
  select(type, pooled)%>%
  left_join(haptoglobin %>%
              drop_na(`Haptoglobin g/L`)%>%
              distinct(Wii_ID,.keep_all = TRUE)%>%
              summarise(median=median(`Haptoglobin g/L`),
                        sd=sd(`Haptoglobin g/L`),
                        lcl=case_when(median-(1.96*sd)>0~median-(1.96*sd),
                                      TRUE~0),
                        ucl=median+(1.96*sd))%>%
              mutate(type="Haptoglobin g/L",
                     single=paste0(median%>%round(2),
                                   " (",
                                   lcl%>%round(2),
                                   "-",
                                   ucl%>%round(2),
                                   ")"))%>%
              select(type, single),
              
            by="type")

haptoglobin %>%
  drop_na(`Haptoglobin g/L`) %>%
  summarise(
    mean = mean(`Haptoglobin g/L`),
    min = min(`Haptoglobin g/L`),
    max = max(`Haptoglobin g/L`),
    n = n()
  )

## plot
ggplot(haptoglobin, aes(y = calf_succ, x = `Haptoglobin g/L`)) +
  ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0, alpha = 0.6, justification = -.2, point_colour = NA) +
  geom_boxplot(position = position_dodge(0.3), width = .15, outlier.shape = NA, alpha = 0.4, show_guide = FALSE) +
  geom_point(position = position_dodge(0.3), alpha = .2, show_guide = FALSE) +
  labs(x = "Haptoglobin (g/L)", y = "Calving success") +
  theme_ipsum() +
  theme(
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.position = "none"
  )

hapto.repro <- ggplot(haptoglobin %>% filter(loc == "Pen"), aes(x = calf_3, y = `Haptoglobin g/L`)) +
  geom_boxplot(width = .5, outlier.shape = NA, alpha = 0.4, show_guide = FALSE) +
  geom_jitter(height = 0.1, width = 0.1, alpha = .2, show_guide = FALSE) +
  labs(y = "Haptoglobin (g/L)", x = "Reproductive outcome") +
  theme_ipsum() +
  theme(
    axis.title.x = element_text(size = 18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = element_text(size = 13, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 13),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.position = "none"
  )


# Calf success outcome binary
haptoglobin %>%
  filter(loc == "Pen") %>%
  kruskal_test(`Haptoglobin g/L` ~ calf_succ)

# Pregnancy
haptoglobin %>%
  filter(loc == "Pen") %>%
  kruskal_test(`Haptoglobin g/L` ~ preg)

# Reproductive outcome all classes
haptoglobin %>%
  filter(loc == "Pen") %>%
  dunn_test(`Haptoglobin g/L` ~ calf_3) %>%
  arrange(p)



hapto.time <- ggplot(haptoglobin) +
  geom_point(aes(x = Year, y = `Haptoglobin g/L`, fill = loc.lastyr, color = loc.lastyr, group = Wii_ID),
    alpha = 0.5
  ) +
  geom_path(aes(x = Year, y = `Haptoglobin g/L`, fill = loc.lastyr, color = loc.lastyr, group = Wii_ID),
    alpha = 0.5
  ) +
  geom_smooth(aes(x = Year, y = `Haptoglobin g/L`),
    method = "lm",
    se = FALSE,
    linetype = "dashed",
    color = "grey50"
  ) +
  labs(x = "Year", y = "Haptoglobin (g/L)") +
  theme_ipsum() +
  theme(
    axis.title.x = element_text(size = 18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.position = "none"
  ) +
  facet_wrap(vars(loc.lastyr))




## effect through time
haptoglobin %>%
  do(tidy(lm(`Haptoglobin g/L` ~ Year * loc.lastyr, .))) %>%
  filter(term != "(Intercept)") %>%
  mutate(sig = case_when(abs(statistic) >= 1.96 ~ "*", TRUE ~ ""))

m5 <- lm(`Haptoglobin g/L` ~ loc.lastyr * Year, data = haptoglobin)
plot(ggpredict(m5, terms = c("Year", "loc.lastyr")))
plot_model(m5, "diag") ## VIF is wonky
summary(m5)


## effect through repeated penning
haptoglobin %>%
  filter(loc == "Pen") %>%
  do(tidy(lm(`Haptoglobin g/L` ~ Year + stay_in_pen, .))) %>%
  filter(term != "(Intercept)") %>%
  mutate(sig = case_when(abs(statistic) >= 1.96 ~ "*", TRUE ~ "")) ## no effect

#### Take home####
## calf effect not significant
## seems like haptoglobin is decreasing through successive stays in pen fits similar to year for penned animals


## plot together
ggarrange(hapto.repro, hapto.time,
  ncol = 2, nrow = 1,
  labels = "AUTO",
  widths = c(1, 1.5)
) %>%
  annotate_figure(top = text_grob("Haptoglobin", color = "black", face = "bold", size = 14))

ggsave(here::here("output", "plots", "haptoglobin_plate.png"), width = 10, height = 5, bg = "white")



#*********************** #
### HAIR CORTISOL #####
#*********************** #

HairCort <- individual_data %>%
  drop_na(hair_cort) %>%
  mutate(weight = as.numeric(weight))

## mean
cort.comparepooled <- individual_data %>%
  drop_na(hair_cort)%>%
  summarise(median=median(hair_cort),
            sd=sd(hair_cort),
            lcl=case_when(median-(1.96*sd)>0~median-(1.96*sd),
                          TRUE~0),
            ucl=median+(1.96*sd))%>%
  mutate(type="Hair cortisol (pg/mg)",
         pooled=paste0(median%>%round(2),
                       " (",
                       lcl%>%round(2),
                       "-",
                       ucl%>%round(2),
                       ")"))%>%
  select(type, pooled)%>%
  left_join(individual_data %>%
              drop_na(hair_cort)%>%
              distinct(Wii_ID,.keep_all = TRUE)%>%
              summarise(median=median(hair_cort),
                        sd=sd(hair_cort),
                        lcl=case_when(median-(1.96*sd)>0~median-(1.96*sd),
                                      TRUE~0),
                        ucl=median+(1.96*sd))%>%
              mutate(type="Hair cortisol (pg/mg)",
                     single=paste0(median%>%round(2),
                                   " (",
                                   lcl%>%round(2),
                                   "-",
                                   ucl%>%round(2),
                                   ")"))%>%
              select(type, single),
            
            by="type")

# remove outlier
HairCort <- HairCort %>%
  filter(hair_cort < 200)


HairCort %>%
  group_by(loc) %>%
  summarise(
    mean = median(hair_cort),
    min = min(hair_cort),
    max = max(hair_cort),
    n = n()
  )

HairCort %>%
  group_by(calf_3) %>%
  drop_na(calf_3) %>%
  summarise(
    mean.cort = mean(hair_cort),
    sd.cort = sd(hair_cort) / sqrt(n())
  ) %>%
  arrange(mean.cort)

HairCort %>%
  group_by(age_clean) %>%
  drop_na(age_clean) %>%
  summarise(
    mean.cort = mean(hair_cort),
    sd.cort = sd(hair_cort) / sqrt(n())
  ) %>%
  arrange(age_clean)

HairCort %>%
  group_by(condition_clean) %>%
  drop_na(condition_clean) %>%
  summarise(
    mean.cort = mean(hair_cort),
    sd.cort = sd(hair_cort) / sqrt(n())
  ) %>%
  arrange(condition_clean)



# correlation with body condition
cort.condition <- ggplot(data = HairCort %>% drop_na(condition_clean), aes(condition_clean, hair_cort)) +
  geom_boxplot(width = .5, outlier.shape = NA, alpha = 0.4, show_guide = FALSE) +
  geom_jitter(height = 0.1, width = 0.1, alpha = .2, show_guide = FALSE) +
  theme_ipsum() +
  labs(x = "Body condition", y = "") +
  theme(
    axis.title.x = element_text(size = 18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 1, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = element_text(size = 13, angle = 25, hjust = 1),
    axis.text.y = element_text(size = 13),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(title = "Location"), color = guide_legend(title = "Location"))

#
# ## test the realtionship between hair cort and body codition via dunn_test
HairCort %>%
  drop_na(condition_clean) %>%
  dunn_test(hair_cort ~ condition_clean)




# # Kruskall-Wallace but revise body condition to be just 'good' and 'not good' ad remove 'Unk'
HairCort %>%
  drop_na(condition_clean) %>%
  mutate(condition_clean = case_when(condition_clean == "good" ~ "good", TRUE ~ "poor/fair")) %>%
  kruskal_test(hair_cort ~ condition_clean)


# Correlation with weight
cort.mass <- ggplot(data = HairCort, aes(as.numeric(weight), hair_cort)) +
  geom_point(alpha = .8, show_guide = FALSE) +
  theme_ipsum() +
  labs(x = "Mass (kg)", y = "") +
  theme(
    axis.title.x = element_text(size = 18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 1, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.position = "bottom"
  )


ggplot(data = HairCort, aes(as.numeric(weight), hair_cort)) +
  geom_point(alpha = .8, show_guide = FALSE) +
  theme_ipsum() +
  labs(x = "Mass (kg)", y = "") +
  theme(
    axis.title.x = element_text(size = 15, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.position = "bottom"
  ) +
  geom_smooth(formula = y ~ x, method = "lm", linetype = "dashed", se = FALSE, color = "black")


HairCort %>%
  # filter(hair_cort<20)%>% ##same result even with those two upper values removed
  do(tidy(lm(hair_cort ~ weight, .))) %>%
  filter(term != "(Intercept)")


# correlation with body fat percentage
ggplot(data = HairCort, aes(body_fat_perc, hair_cort)) +
  geom_point(alpha = .8, show_guide = FALSE) +
  theme_ipsum() +
  labs(x = "Fat (%)", y = "Hair cortisol") +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.position = "bottom"
  ) +
  geom_smooth(formula = y ~ x, method = "lm", linetype = "dashed", se = FALSE, color = "black")


HairCort %>%
  # filter(hair_cort<20)%>% ##same result even with those two upper values removed
  do(tidy(lm(hair_cort ~ body_fat_perc, .))) %>%
  filter(term != "(Intercept)")

### calf outcomes
## pen only
cort.calf <- ggplot(data = HairCort %>% drop_na(calf_3, stay_in_pen), aes(x = calf_3, y = hair_cort)) +
  geom_boxplot(width = .5, outlier.shape = NA, alpha = 0.4, show_guide = FALSE) +
  geom_jitter(height = 0.1, width = 0.1, alpha = .2, show_guide = FALSE) +
  theme_ipsum() +
  labs(y = "Hair cortisol (pg/mg)", x = "Reproductive outcome") +
  theme(
    axis.title.x = element_text(size = 18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = element_text(size = 13, angle = 25, hjust = 1),
    axis.text.y = element_text(size = 13),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.position = "bottom"
  )


# Calf success outcome binary
HairCort %>%
  drop_na(calf_succ) %>%
  kruskal_test(hair_cort ~ calf_succ)

# Pregnancy
HairCort %>%
  drop_na(preg) %>%
  kruskal_test(hair_cort ~ preg)

# Reproductive outcome all classes
HairCort %>%
  drop_na(calf_3) %>%
  dunn_test(hair_cort ~ calf_3) %>%
  arrange(p) ## weak effect of Alive vs stillborn



## not sig for preg/not or calf result

####### Compare with other herds ######

## add in latest KZ pen data
cort.plot <- cort.other %>%
  mutate(Year = year(collection_date)) %>%
  select(WLH_ID = wlh_id, Year, Herd = study_area, hair_cort = hcc_clean) %>%
  rbind(HairCort %>%
    select(WLH_ID, Year, Herd = herd, hair_cort) %>%
    mutate(Herd = "Klinse-Za")) %>%
  rbind(bor %>%
    mutate(
      Herd = "BC Boreal",
      Year = year(`Capture (dd/mm/yy)`),
      WLH_ID = str_sub(`Caribou ID`, 1, 5)
    ) %>%
    select(WLH_ID,
      Year,
      Herd,
      hair_cort = `HCC (pg/mg)`
    )) %>%
  mutate(herd.stats = relevel(factor(Herd), ref = "Klinse-Za"))


herd.order <- cort.plot %>%
  group_by(Herd) %>%
  summarise(cort.med = median(`hair_cort`)) %>%
  arrange(cort.med)

cort.compare <- ggplot(
  data = cort.plot %>% filter(!Herd %in% "Scott"), ## only one record
  aes(fct_relevel(Herd, herd.order$Herd), hair_cort, color = Herd == "Klinse-Za")
) +
  geom_boxplot(position = position_dodge(0.3), width = .4, outlier.shape = NA, alpha = 0.4, show_guide = FALSE) +
  geom_jitter(height = 0.1, width = 0.1, alpha = .2, show_guide = FALSE) +
  coord_cartesian(ylim = c(0, 10)) +
  theme_ipsum() +
  labs(x = "Caribou herd", y = "Hair cortisol (pg/mg)") +
  theme(
    axis.title.x = element_text(size = 18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = element_text(size = 13, angle = 35, hjust = 1),
    axis.text.y = element_text(size = 13),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.position = "none"
  ) +
  annotate(geom = "text", x = 3.5, y = 9, label = "6 ONCP, 22 boreal, & 13 KZ\nrecords not shown\nthat ranged from 11-48", color = "black", size = 2.5) +
  guides(fill = guide_legend(title = "Location"), color = guide_legend(title = "Location"))

## what data is missing in plot?
cort.plot %>%
  filter(!Herd %in% "Scott") %>%
  filter(hair_cort > 10) %>%
  group_by(Herd) %>%
  summarise(
    n = n(),
    min = min(hair_cort),
    max = max(hair_cort)
  )


## stats on how herds are different from one another
# glm(hair_cort~herd.stats,data=cort.plot%>%filter(!Herd%in%"Scott"))%>%summary()
# glm(hair_cort~herd.stats=="Klinse-Za",data=cort.plot)%>%summary()

cort.plot %>%
  mutate(herd = case_when(
    herd.stats %in% "Klinse-Za" ~ "KZ",
    herd.stats %in% "BC Boreal" ~ "BCbor",
    TRUE ~ "ONCP"
  )) %>%
  ggplot(aes(x = herd, y = hair_cort)) +
  geom_boxplot()

cort.plot %>%
  mutate(herd = case_when(
    herd.stats %in% "Klinse-Za" ~ "KZ",
    herd.stats %in% "BC Boreal" ~ "BCbor",
    TRUE ~ "ONCP"
  )) %>%
  dunn_test(hair_cort ~ herd)

cort.plot %>%
  dunn_test(hair_cort ~ herd.stats) %>%
  filter(group1 %in% "Klinse-Za")



####### Cort change during successive penning ######
cort.penstay <- ggplot(individual_data %>% arrange(Wii_ID, stay_in_pen) %>% drop_na(hair_cort) %>% filter(hair_cort < 200), aes(x = stay_in_pen, y = hair_cort)) +
  geom_point(aes(group = Wii_ID),
    alpha = 0.5, color = "#00BFC4"
  ) +
  geom_path(aes(group = Wii_ID),
    alpha = 0.5, color = "#00BFC4"
  ) +
  geom_smooth(
    data = individual_data %>% arrange(Wii_ID, stay_in_pen) %>% drop_na(hair_cort) %>% filter(hair_cort < 200),
    method = "lm",
    se = FALSE,
    linetype = "dashed",
    color = "grey50"
  ) +
  labs(x = "Stays in pen (n)", y = "") +
  theme_ipsum() +
  theme(
    axis.title.x = element_text(size = 18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 1, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.position = "none"
  )


cort.time <- ggplot(individual_data %>% arrange(Wii_ID, stay_in_pen) %>% drop_na(hair_cort) %>% filter(hair_cort < 200), aes(x = Year, y = hair_cort)) +
  geom_point(aes(group = Wii_ID),
    alpha = 0.5, color = "#00BFC4"
  ) +
  geom_path(aes(group = Wii_ID),
    alpha = 0.5, color = "#00BFC4"
  ) +
  geom_smooth(
    data = individual_data %>% arrange(Wii_ID, stay_in_pen) %>% drop_na(hair_cort) %>% filter(hair_cort < 200),
    method = "lm",
    formula = y ~ poly(x, 2),
    se = FALSE,
    linetype = "dashed",
    color = "grey50"
  ) +
  labs(x = "Year", y = "") +
  theme_ipsum() +
  theme(
    axis.title.x = element_text(size = 18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 1, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.position = "none"
  )


HairCort %>%
  filter(loc == "Pen") %>%
  do(broom.mixed::tidy(lm(hair_cort ~ stay_in_pen, .))) %>%
  filter(term %in% "stay_in_pen") %>%
  arrange(-abs(statistic)) %>%
  mutate(sig = case_when(abs(statistic) >= 1.96 ~ "*", TRUE ~ ""))


## change through time?
HairCort %>%
  do(broom.mixed::tidy(lm(hair_cort ~ Year, .))) %>%
  arrange(-abs(statistic)) %>%
  mutate(sig = case_when(abs(statistic) >= 1.96 ~ "*", TRUE ~ "")) ## no


## stays in pen vs time?
year <- HairCort %>%
  drop_na(stay_in_pen, hair_cort) %>%
  filter(loc == "Pen") %>%
  lm(hair_cort ~ Year, .)

year2 <- HairCort %>%
  drop_na(stay_in_pen, hair_cort) %>%
  filter(loc == "Pen") %>%
  lm(hair_cort ~ Year + I(Year^2), .)

stays <- HairCort %>%
  drop_na(stay_in_pen, hair_cort) %>%
  filter(loc == "Pen") %>%
  lm(hair_cort ~ stay_in_pen, .)


anova(stays, year)
AIC(year, stays, year2)


## plot together
ggarrange(cort.compare, cort.mass, cort.condition, cort.calf, cort.penstay, cort.time,
  labels = "AUTO",
  ncol = 3,
  nrow = 2
) %>%
  annotate_figure(top = text_grob("Hair Cortisol", color = "black", face = "bold", size = 14))

ggsave(here::here("output", "plots", "hair.cort_plate.png"), width = 10, height = 8, bg = "white")


m1 <- lm(hair_cort ~ Year * loc.lastyr, HairCort)
plot(ggpredict(m1, terms = c("Year", "loc.lastyr")))

ggplot(HairCort, aes(x = Year, y = hair_cort, color = loc.lastyr, shape = loc.lastyr)) +
  geom_jitter(width = 0.1) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2))

cort.time2 <- ggplot(HairCort, aes(x = factor(Year), y = hair_cort, fill = factor(loc.lastyr), shape = factor(loc.lastyr))) +
  geom_boxplot(position = position_dodge(0.3), width = .4, outlier.shape = NA, alpha = 0.4) +
  geom_point(position = position_dodge(0.3), height = 0, width = 0, alpha = .6, show_guide = FALSE) +
  theme_ipsum() +
  labs(x = "", y = "Hair cortisol (pg/mg)", fill = "Location", shape = "Location") +
  theme(
    axis.title.x = element_text(size = 15, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15)
  )


cort.compare.time <- ggplot(cort.plot, aes(x = Year, y = hair_cort, color = Herd)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_ipsum() +
  labs(x = "Year", y = "Hair cortisol (pg/mg)") +
  theme(
    axis.title.x = element_text(size = 15, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15)
  )+
  labs(color="Subpopulation")

cort.plot %>%
  group_by(Herd) %>%
  do(tidy(glm(hair_cort ~ Year, data = .))) %>%
  filter(term != "(Intercept)")

cort.plot %>%
  do(broom.mixed::tidy(lm(hair_cort ~ Year, .))) %>%
  arrange(-abs(statistic)) %>%
  mutate(sig = case_when(abs(statistic) >= 1.96 ~ "*", TRUE ~ "")) ## no


ggarrange(cort.time2, cort.compare.time,
  labels = "AUTO",
  ncol = 1,
  nrow = 2
) %>%
  annotate_figure(top = text_grob("Hair Cortisol Trends", color = "black", face = "bold", size = 14))

ggsave(here::here("output", "plots", "hair.cort_plate_trends.png"), width = 7, height = 8, bg = "white")





#### Take home####

## not enough data to compare penned an unpenned, just use penned.

## caribou in good body condition generally had lower cort, and those with heavier weights had lower cort, both weakly sig (p=0.09, but p<0.05 with kruskal wallis)

## Caribou F in pen who were pregnant generally had lower hair cort, but the effect wasn't significant
## removing an outlier (cort>7, confirmed outlier via qqplot) suggests p<0.03 and no violation of normality
## BUT, a possible other outliers, cor=5.6 and cor=1.91, when removed makes p=0.22. Fragile result

## KZ pen stress levels similar to other herds, of 6 herds, ranks 2nd lowest

## not enough data for RE on individual. But raw regression shows hair cor increasing through successive penning, but this is driven by a single point.
## 315S in 2018. Remove this point and regression no longer significant.
## overall dont think there is a penning effect on cort. No effect through time either.



#*********************** #
### PATHOGENS #####
#*********************** #

pathogens <- individual_data %>%
  select(Wii_ID, Year, calf_succ, preg, calf_3, stay_in_pen, Toxoplasma:Erysip_simplified, loc) %>%
  drop_na(calf_succ) %>%
  pivot_longer(Toxoplasma:Erysip_simplified) %>%
  mutate(
    value = case_when(
      value %in% c("Positive", "positive") ~ "Positive",
      value %in% c("Negative", "negative") ~ "Negative",
      TRUE ~ value
    ),
    name = case_when(
      name == "Erysip_simplified" ~ "Erysipelothrix",
      name == "AHC_IBR" ~ "Alphaherpesvirus",
      TRUE ~ name
    )
  ) %>%
  drop_na(value)



##prevelance
patho.prevalence <- pathogens%>%
  group_by(name)%>%
  #distinct(Wii_ID,.keep_all=TRUE)%>%
  distinct(Wii_ID,Year, .keep_all=TRUE)%>%
  summarize(pos=sum(value=="Positive"),
            neg=sum(value=="Negative"),
            total=pos+neg,
            prev=pos/total,
            sd=sqrt((prev*(1-prev))/total),
            lcl=case_when(prev-(1.96*sd)>0~prev-(1.96*sd),
                          TRUE~0),
            ucl=prev+(1.96*sd))


patho.comparepooled  <- patho.prevalence%>%
  mutate(pooled=paste0(prev%>%round(2),
                       " (",
                       lcl%>%round(2),
                       "-",
                       ucl%>%round(2),
                       ")"))%>%
  select(name,pooled)%>%
  left_join(pathogens%>%
              group_by(name)%>%
              distinct(Wii_ID,.keep_all=TRUE)%>%
              summarize(pos=sum(value=="Positive"),
                        neg=sum(value=="Negative"),
                        total=pos+neg,
                        prev=pos/total,
                        sd=sqrt((prev*(1-prev))/total),
                        lcl=case_when(prev-(1.96*sd)>0~prev-(1.96*sd),
                                      TRUE~0),
                        ucl=prev+(1.96*sd))%>%
              mutate(single=paste0(prev%>%round(2),
                                   " (",
                                   lcl%>%round(2),
                                   "-",
                                   ucl%>%round(2),
                                   ")"))%>%
              select(name,single),
            by="name")%>%
  rename(type=name)

##print all compares
nutr.comparepooled%>%
  rbind(patho.comparepooled)%>%
  rbind(cort.comparepooled)%>%
  rbind(hapto.comparepooled)%>%
  write_csv(here::here("output/tables/compare_pooled.csv"))
  


## plot
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall = k))


pathogens.repro <- pathogens %>%
  group_by(name, calf_succ) %>%
  count(value) %>%
  mutate(prop = (n / sum(n))) %>%
  ggplot(aes(y = calf_succ, x = value, color = prop, size = prop)) +
  geom_point(alpha = .5, show_guide = FALSE) +
  geom_label(aes(label = n)) +
  labs(y = "Calving success", x = "Test result") +
  theme_ipsum() +
  theme(
    axis.title.x = element_text(size = 18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.position = "none"
  ) +
  facet_wrap(vars(name), scales = "free_x")

# ggsave(here::here("output","plots","pathogens.png"),width=6,height=6, bg="white")

## statistics

library(yardstick)
pathogens %>%
  filter(loc == "Pen") %>%
  group_by(name) %>%
  drop_na(value, calf_succ) %>%
  mutate(
    value = case_when(
      value %in% "Positive" ~ "1",
      value %in% "Negative" ~ "0"
    ) %>%
      as.factor(),
    calf_succ = case_when(
      calf_succ %in% "Yes" ~ "1",
      calf_succ %in% "No" ~ "0"
    ) %>%
      as.factor()
  ) %>%
  metrics(value, calf_succ)


pathogens %>%
  filter(loc == "Pen") %>%
  group_by(name) %>%
  drop_na(value, calf_succ) %>%
  mutate(
    value = case_when(
      value %in% "Positive" ~ "1",
      value %in% "Negative" ~ "0"
    ) %>%
      as.factor(),
    preg = case_when(
      preg %in% "preg" ~ 1,
      preg %in% "not_preg" ~ 0,
      TRUE ~ NA_real_
    ) %>%
      as.factor()
  ) %>%
  metrics(value, preg)




## compare to other herds
pathogens.plot <- pathogens.other %>%
  mutate(name = case_when(
    name == "Herpes_IBR" ~ "Alphaherpesvirus",
    TRUE ~ name
  )) %>%
  rbind(pathogens %>%
    mutate(
      Sex = "F",
      Herd = case_when(
        loc == "Pen" ~ "Klinse-Za",
        loc == "Free-ranging" ~ "Klinse-Za"
      )
    ) %>%
    select(Wii_ID, Sex, Herd, name, value)) %>%
  mutate(
    Herd = relevel(factor(Herd), ref = "Klinse-Za"),
    value = case_when(
      value == "borderline negative" ~ "Negative",
      value == "borderline positive" ~ "Positive",
      value == "positive" ~ "Positive",
      value == "negative" ~ "Negative",
      TRUE ~ value
    )
  )

##prevelance
pathogens.plot%>%
  mutate(group=case_when(!Herd%in%c("Klinse-Za", "Boreal")~"ONCP",
                         TRUE~Herd))%>%
  group_by(name, group)%>%
  #distinct(Wii_ID,.keep_all=TRUE)%>%
  summarize(pos=sum(value=="Positive"),
            neg=sum(value=="Negative"),
            total=pos+neg,
            prev=pos/total,
            sd=sqrt((prev*(1-prev))/total),
            lcl=prev-(1.96*sd),
            ucl=prev+(1.96*sd))

pathogens.compare <- pathogens.plot %>%
  filter(
    !Herd %in% "Scott" & ## only one record
      Sex == "F",
    name != "Toxoplasma"
  ) %>% ## not in ONCP data
  group_by(Herd, name) %>%
  count(value) %>%
  pivot_wider(names_from = value, values_from = n) %>%
  replace(is.na(.), 0) %>%
  rbind(tibble(
    Herd = "BC Boreal",
    name = unique(pathogens.plot %>% filter(name != "Toxoplasma") %>% pull(name)),
    Negative = c(
      (149 + 40 + 20 + 17) - (21 + 17 + 10 + 7),
      243 - 5,
      (155 + 39 + 25 + 1) - (98 + 24 + 15 + 1)
    ),
    Positive = c(
      (21 + 17 + 10 + 7),
      5,
      (98 + 24 + 15 + 1)
    )
  )) %>%
  mutate(Pos = round((Positive / sum(Positive + Negative)) * 100, 0)) %>%
  ungroup() %>%
  mutate(Herd = relevel(as.factor(Herd), ref = "Klinse-Za")) %>%
  ggplot(aes(Herd, Pos, fill = Herd == "Klinse-Za")) +
  geom_col() +
  theme_ipsum() +
  labs(x = "", y = "Pathogen prevelance (%)") +
  theme(
    axis.title.x = element_text(size = 18, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = element_text(size = 13, angle = 40, hjust = 1),
    axis.text.y = element_text(size = 13),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.position = "none"
  ) +
  guides(fill = guide_legend(title = "Location"), color = guide_legend(title = "Location")) +
  facet_wrap(vars(name), scales = "free_y") +
  scale_y_continuous(expand = c(0.1, 1)) +
  geom_text(aes(label = paste0("n=", Positive + Negative)), vjust = -0.5, size = 3)




pathogens.comp.betas <- pathogens.plot %>%
  filter(
    !Herd %in% "Scott" & ## only one record
      Sex == "F",
    name != "Toxoplasma"
  ) %>% ## not in ONCP data
  group_by(name) %>%
  mutate(
    value = case_when(
      value %in% "Positive" ~ 1,
      value %in% "Negative" ~ 0,
      TRUE ~ NA_real_
    ),
    cull = sum(value)
  ) %>%
  ungroup() %>%
  filter(cull > 1) %>%
  group_by(name) %>%
  do(tidy(glm(value ~ Herd != "Klinse-Za", data = ., family = "binomial"))) %>%
  filter(term != "(Intercept)") %>%
  mutate(term = "ONCP") %>%
  ggplot(aes(x = estimate, xmin = estimate - (std.error * 1.96), xmax = estimate + (std.error * 1.96), y = name)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_linerange(position = position_dodge(width = -0.5)) +
  geom_point(position = position_dodge(width = -0.5)) +
  labs(x = "Estimate (relative to Klinse-Za, standardized)", y = "Term", color = "") +
  theme_ipsum() +
  theme(
    axis.title.x = element_text(size = 15, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 18, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(nrow = 3, byrow = TRUE))

## plot together
ggarrange(pathogens.comp.betas, pathogens.repro,
  ncol = 2, nrow = 1,
  labels = c("B", "C"),
  widths = c(0.7, 1)
) %>%
  ggarrange(pathogens.compare, .,
    ncol = 1, nrow = 2,
    labels = c("A", NA),
    heights = c(1, 1)
  ) %>%
  annotate_figure(top = text_grob("Pathogens", color = "black", face = "bold", size = 14))

ggsave(here::here("output", "plots", "pathogen_compare_plate.png"), width = 10, height = 8, bg = "white")


#### Take home####
## Disease prevalence noted in the population, but no notable effects on reproduction, or effects from penning
## accuracy and kappa both quite poor




##### FECAL ANALYSES WHERE WE DON'T KNOW INDIVIDUAL


#*********************** #
### FECAL CORTISOL ####
#*********************** #

# remove missing cortisol values and clean up month and survey type column
FecalCort <- FecalCort %>%
  drop_na(cortisol) %>%
  mutate(
    month = month(date, label = TRUE),
    doy = yday(date),
    year = year(date)
  )

## get doy 0, i.e., min of days surveyed
FecalCort <- FecalCort %>%
  mutate(doy.zero = doy - min(doy))

## add in Sex and ID from genetics
FecalCort <- FecalCort %>%
  left_join(FecalCort.id %>%
    rename(wii_sample_id = sample_ID))


# summary stats
## by pen status
FecalCort %>%
  group_by(survey_type) %>%
  summarise(
    cort = mean(cortisol),
    n = n()
  )

## by pen status and season
FecalCort %>%
  group_by(month, survey_type) %>%
  summarise(cort = median(cortisol))

FecalCort %>%
  drop_na(Sex) %>%
  group_by(survey_type, Sex) %>%
  summarise(
    cort = median(cortisol),
    n = n()
  )

## check if males and females are statisticcally different?
FecalCort %>%
  drop_na(Sex) %>%
  do(broom.mixed::tidy(lmer(cortisol ~ doy.zero * survey_type + Sex + (1 | year) + (1 | Year), .))) %>%
  filter(effect == "fixed" & term != "(Intercept)") %>%
  arrange(-abs(statistic)) %>%
  mutate(sig = case_when(abs(statistic) >= 1.96 ~ "*", TRUE ~ ""))


FecalCort %>%
  drop_na(Sex) %>%
  ggplot(aes(y = cortisol, x = doy, color = Sex)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(vars(survey_type)) +
  theme_ipsum() +
  labs(x = "Day of year", y = "Cortisol") +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.position = "bottom"
  ) +
  ylim(0, 1000)

ggplot(FecalCort %>% drop_na(Sex) %>% filter(survey_type == "Free-ranging"), aes(y = Sex, x = cortisol, fill = survey_type, color = survey_type)) +
  ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0, alpha = 0.6, justification = -.2, point_colour = NA) +
  geom_boxplot(position = position_dodge(0.3), width = .15, outlier.shape = NA, alpha = 0.4, show_guide = FALSE) +
  geom_point(position = position_dodge(0.3), alpha = .2, show_guide = FALSE) +
  theme_ipsum() +
  labs(x = "Fecal glucocorticoid metabolites (ng/g)", y = "Sex") +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(title = "Location"), color = guide_legend(title = "Location"))



## no evidence that males have different cortisol, keep together in analysis

## plot
ggplot(FecalCort, aes(y = month, x = cortisol, fill = survey_type, color = survey_type)) +
  ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0, alpha = 0.6, justification = -.2, point_colour = NA) +
  geom_boxplot(position = position_dodge(0.3), width = .15, outlier.shape = NA, alpha = 0.4, show_guide = FALSE) +
  geom_point(position = position_dodge(0.3), alpha = .2, show_guide = FALSE) +
  theme_ipsum() +
  labs(x = "Fecal glucocorticoid metabolites (ng/g)", y = "Sampling month") +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(title = "Location"), color = guide_legend(title = "Location"))


ggplot(FecalCort, aes(y = cortisol, x = doy %>% as.Date(origin = "2013-12-31"), fill = survey_type, color = survey_type)) +
  geom_point(alpha = .2, show_guide = FALSE) +
  theme_ipsum() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, alpha = 0.5) +
  labs(y = "Fecal glucocorticoid metabolites (ng/g)", x = "Time") +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(title = "Location"), color = guide_legend(title = "Location"))


## statistical analysis

## need year RE?
m1 <- lmer(cortisol ~ doy.zero * survey_type + (1 | year), data = FecalCort)
m1a <- lm(cortisol ~ doy.zero * survey_type, data = FecalCort)
anova(m1, m1a)

## yes use RE
pen.cort <- lmer(cortisol ~ doy.zero * survey_type + (1 | year), data = FecalCort)

pen.cort %>%
  tidy() %>%
  filter(effect == "fixed" & term != "(Intercept)") %>%
  arrange(-abs(statistic)) %>%
  mutate(sig = case_when(abs(statistic) >= 1.96 ~ "*", TRUE ~ ""))

ranef(pen.cort)

plot(ggpredict(pen.cort, terms = c("doy.zero", "survey_type")))


#### Take home####
## no effect of sex, keep sexes pooled.
## pen and unpenned start off at same faecal cort levels but through time diverge such that pen is lower

#*********************** #
### FECAL NITROGEN ####
#*********************** #
# join fecal nitrogen and cortisol data -- this represents fecal nitrogen data from the pen/wild summer sampling (not capture). It also is a subset of the fecal
# cort data we have so do not analyze feca cortisol here

FecalNit <- FecalCort %>%
  left_join(fecal.nit %>% select(wii_sample_id, fecal_nitro_perc, survey_type.nit = survey_type), by = "wii_sample_id")

## differences in % N in or out of pen?
FecalNit %>%
  drop_na(fecal_nitro_perc) %>%
  kruskal_test(fecal_nitro_perc ~ survey_type.nit)

FecalNit %>%
  drop_na(fecal_nitro_perc) %>%
  group_by(survey_type.nit) %>%
  summarise(
    mean = mean(fecal_nitro_perc),
    sd = sd(fecal_nitro_perc) / sqrt(n())
  )

## check if males and females are stastically different?
FecalNit %>%
  drop_na(fecal_nitro_perc) %>%
  lmer(fecal_nitro_perc ~ doy.zero * survey_type + Sex + (1 | year), data = .) %>%
  tidy() %>%
  filter(effect == "fixed" & term != "(Intercept)") %>%
  arrange(-abs(statistic)) %>%
  mutate(sig = case_when(abs(statistic) >= 1.96 ~ "*", TRUE ~ ""))


FecalNit %>%
  drop_na(Sex) %>%
  ggplot(aes(y = fecal_nitro_perc, x = doy, color = Sex)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(vars(survey_type))

ggplot(FecalNit %>%
  drop_na(Sex) %>% filter(survey_type == "Free-ranging"), aes(y = Sex, x = fecal_nitro_perc, fill = survey_type, color = survey_type)) +
  ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0, alpha = 0.6, justification = -.2, point_colour = NA) +
  geom_boxplot(position = position_dodge(0.3), width = .15, outlier.shape = NA, alpha = 0.4, show_guide = FALSE) +
  geom_point(position = position_dodge(0.3), alpha = .2, show_guide = FALSE) +
  theme_ipsum() +
  labs(x = "Fecal nitrogen (%)", y = "Sex") +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(title = "Location"), color = guide_legend(title = "Location"))


## seems similar, pool.


### plot
ggplot(FecalNit, aes(y = month, x = fecal_nitro_perc, fill = survey_type, color = survey_type)) +
  ggdist::stat_halfeye(adjust = .5, width = .6, .width = 0, alpha = 0.6, justification = -.2, point_colour = NA) +
  geom_boxplot(position = position_dodge(0.3), width = .15, outlier.shape = NA, alpha = 0.4, show_guide = FALSE) +
  geom_point(position = position_dodge(0.3), alpha = .2, show_guide = FALSE) +
  theme_ipsum() +
  labs(x = "Fecal nitrogen (%)", y = "Sampling month") +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(title = "Location"), color = guide_legend(title = "Location"))

ggplot(FecalNit, aes(y = fecal_nitro_perc, x = doy %>% as.Date(origin = "2013-12-31"), fill = survey_type, color = survey_type)) +
  geom_point(alpha = .2, show_guide = FALSE) +
  theme_ipsum() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, alpha = 0.5) +
  labs(y = "Fecal nitrogen (%)", x = "Day of year") +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(title = "Location"), color = guide_legend(title = "Location"))

## statistical analysis
## need year RE?
m1 <- lmer(fecal_nitro_perc ~ doy.zero * survey_type + (1 | year), data = FecalNit)
m1a <- lm(fecal_nitro_perc ~ doy.zero * survey_type, data = FecalNit)
anova(m1, m1a)

## yes need RE
pen.nit <- lmer(fecal_nitro_perc ~ doy.zero * survey_type + (1 | year), data = FecalNit)


pen.nit %>%
  tidy() %>%
  filter(effect == "fixed", term != "(Intercept)") %>%
  arrange(-abs(statistic)) %>%
  mutate(sig = case_when(abs(statistic) >= 1.96 ~ "*", TRUE ~ ""))

ranef(pen.nit)

plot(ggpredict(pen.nit, terms = c("doy.zero", "survey_type")))

## noticed trend in RE. Include as fixed effect
pen.nit <- lm(fecal_nitro_perc ~ doy.zero * survey_type + year * survey_type, data = FecalNit %>% mutate(year = year - 2014))

pen.nit %>%
  tidy() %>%
  filter(term != "(Intercept)") %>%
  arrange(-abs(statistic)) %>%
  mutate(sig = case_when(abs(statistic) >= 1.96 ~ "*", TRUE ~ ""))

plot(ggpredict(pen.nit, terms = c("survey_type", "year")))
summary(pen.nit)


ggplot(FecalNit, aes(y = fecal_nitro_perc, x = as.factor(year), fill = survey_type, color = survey_type)) +
  geom_boxplot(position = position_dodge(0.8), alpha = .2) +
  geom_jitter(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.1), alpha = .2, show_guide = FALSE) +
  theme_ipsum() +
  labs(y = "Fecal nitrogen (%)", x = "Year") +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(title = "Location"), color = guide_legend(title = "Location"))


## plot together


ggplot(FecalNit %>% mutate(value = fecal_nitro_perc, name = "Fecal nitrogen (%)") %>% select(name, doy, value, survey_type) %>%
  rbind(FecalCort %>% mutate(value = cortisol, name = "Fecal glucocorticoid metabolites (ng/g)") %>% select(name, doy, value, survey_type)), aes(y = value, x = doy %>% as.Date(origin = "2013-12-31"), fill = survey_type, color = survey_type)) +
  geom_point(alpha = .2, show_guide = FALSE) +
  theme_ipsum() +
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, alpha = 0.5) +
  labs(y = "Value", x = "Day of year", title = "Fecal Health Metrics") +
  theme(
    axis.title.x = element_text(size = 15, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(size = 15, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.position = "bottom"
  ) +
  guides(fill = guide_legend(title = "Location"), color = guide_legend(title = "Location")) +
  facet_wrap(vars(name), nrow = 2, scales = "free_y")


ggsave(here::here("output", "plots", "fecal.cort-nit.season.png"), width = 4, height = 7, bg = "white")


#### Take home####
## no effect of sex, keep sexes pooled.
## pen slightly higher fecal nitrogen
