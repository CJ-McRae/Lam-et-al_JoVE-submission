####JoVe manuscript data analysis####

#Lam et al.
#Effective techniques for feeding and ex situ culture of a brooding scleractinian coral, Pocillopora acuta

#February 15 2023 (original submission)
#April 24 2023 (revision 1)

#clear workspace
rm(list = ls()) 

#set language to English
Sys.setenv(LANG = "en")

####_______________####

####<Reproductive Output>####

####1) Libraries####

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(car)

####2) Data####

data_RO <- read_csv("data_RO.csv")

#class check

summary(data_RO)

data_RO$Gregorian_month <-as.factor(data_RO$Gregorian_month)
data_RO$lunar_month_year <-as.factor(data_RO$lunar_month_year)
data_RO$colony_ID <-as.factor(data_RO$colony_ID)
data_RO$treatment <- as.factor(data_RO$treatment)
data_RO$feeding_treatment <-as.factor(data_RO$feeding_treatment)
data_RO$temp_treatment <-as.factor(data_RO$temp_treatment)
data_RO$tank_ID <- as.factor(data_RO$tank_ID)
data_RO$RO <- as.integer(data_RO$RO)

summary(data_RO)
head(data_RO)


#colony size comparison
TLE.anova <- aov(TLE~treatment, data=data_RO)
summary(TLE.anova)

mean(data_RO$TLE)
#21.34958
sd(data_RO$TLE)
#2.757824



####3) Plot####

#set x axis order
data_RO$treatment <- factor(data_RO$treatment, 
                            levels=c("FH", "NH", "FL", "NL"))

#RO box plot
RO_box <- 
  ggplot(data=data_RO, aes(x=treatment, y=RO, fill=treatment)) + 
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values = c("brown3", "lightpink", "royalblue3", "skyblue"))+
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), 
             size=1.3) +
  labs(y = "# of planulae released", x="Treatment", title ="")+
  scale_y_continuous(limits=c(0, 600), breaks=c(0, 100, 200, 300, 400, 500, 600))+
  theme_classic()

RO_box + theme(text=element_text(size=14,  family="serif")) +
  theme(legend.position = "none")


#remove outlier
data_RO_no.outlier <- 
  data_RO %>%  
  filter(!row_number() %in% c(20))

RO_box_no.outlier <- 
  ggplot(data=data_RO_no.outlier, aes(x=treatment, y=RO, fill=treatment)) + 
  geom_boxplot(outlier.shape = NA)+
  scale_fill_manual(values = c("brown3", "lightpink", "royalblue3", "skyblue"))+
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7), 
             size=1.3) +
  labs(y = "# of planulae released", x="Treatment", title ="")+
  scale_y_continuous(limits=c(0, 400), breaks=c(0, 100, 200, 300, 400))+
  theme_classic()

RO_box_no.outlier + theme(text=element_text(size=14,  family="serif")) +
  theme(legend.position = "none")

#mean RO and sd for each treatment
RO_avg_no.outlier <- data_RO_no.outlier%>% filter(!is.na(RO))

RO_avg_no.outlier_SUMMARY <-
  RO_avg_no.outlier%>%
  group_by(treatment)%>%
  summarise(mean.treat = mean(RO), sd.treat =sd(RO))

RO_avg_no.outlier_SUMMARY


####4) Analysis####

#RO analysis without outlier colony #L1N2

#take a look at the data
plot(data_RO_no.outlier$RO)
hist(data_RO_no.outlier$RO)

#glmer approach (feed*temp)
RO_no.outlier.mm <- glmer(RO~feeding_treatment*temp_treatment + (1|tank_ID), data=data_RO_no.outlier, family = "poisson")
summary(RO_no.outlier.mm)

#check assumptions 
plot(fitted(RO_no.outlier.mm),residuals(RO_no.outlier.mm))
hist(residuals(RO_no.outlier.mm))
qqnorm(residuals(RO_no.outlier.mm))
vif(RO_no.outlier.mm)

RO_no.outlier.mm_posthoc <- emmeans(RO_no.outlier.mm, pairwise~feeding_treatment*temp_treatment)
RO_no.outlier.mm_posthoc 


#just for reference: all RO data (i.e., no outlier removed)

#take a look at the data
plot(data_RO$RO)
hist(data_RO$RO)

#glmer approach (feed*temp)
RO.mm <- glmer(RO~feeding_treatment*temp_treatment + (1|tank_ID), data=data_RO, family = "poisson")
summary(RO.mm)

#check assumptions 
plot(fitted(RO.mm),residuals(RO.mm))
hist(residuals(RO.mm))
qqnorm(residuals(RO.mm))
vif(RO.mm)

RO.mm_posthoc <- emmeans(RO.mm, pairwise~feeding_treatment*temp_treatment)
RO.mm_posthoc 


####________________####

####<Reproductive Timing>####


####1) Libraries####

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(car)
library(Hmisc)


####2) Data####

data_RT <- read_csv("data_RT.csv")

#class check
summary(data_RT)

data_RT$lunar_month <-as.factor(data_RT$lunar_month)
data_RT$lunar_day <-as.numeric(data_RT$lunar_day)
data_RT$treatment <- as.factor(data_RT$treatment)
data_RT$feeding_treatment <-as.factor(data_RT$feeding_treatment)
data_RT$temp_treatment <-as.factor(data_RT$temp_treatment)
data_RT$tank_id <- as.factor(data_RT$tank_id)
data_RT$colony_id <-as.factor(data_RT$colony_id)
data_RT$RO <- as.integer(data_RT$RO)

summary(data_RT)
head(data_RT)

HF.subset <-
  data_RT%>%
  filter(treatment =="HF")

HN.subset <-
  data_RT%>%
  filter(treatment =="HN")

LF.subset <-
  data_RT%>%
  filter(treatment =="LF")

LN.subset <-
  data_RT%>%
  filter(treatment =="LN")


####3) Plots####

###first calculate weighted mean and sd

#remove NA rows
data_RT_no.NA <- 
  data_RT %>% 
  filter(!is.na(RO))

#add in monthly sum to assess RT weight based on RO
data_RT_no.NA_sum <- 
  data_RT_no.NA%>%
  group_by(colony_id) %>%
  mutate(month.sum = sum(RO , na.rm= T))

#calculate RT weight
data_RT_no.NA_sum <- 
  data_RT_no.NA_sum%>%
  mutate(RT_weight = RO/month.sum)

#calculate weighted mean of lunar day
WT.MLD <- 
  data_RT_no.NA_sum%>%
  group_by(colony_id, treatment, feeding_treatment, temp_treatment, tank_id)%>%
  summarise(W.MLD = weighted.mean(lunar_day,RT_weight))

#weighted mean and sd per season for treatment
data_RT_no.NA_sum_noNA <- data_RT_no.NA_sum %>% filter(!is.na(RT_weight))

WT.MLD.SD <- 
  data_RT_no.NA_sum_noNA%>%
  group_by(treatment)%>%
  summarise(W.MLD = weighted.mean(lunar_day,RT_weight), wt.var = wtd.var(lunar_day, RT_weight), wt.sd = sqrt(wt.var))

WT.MLD.SD 

# A tibble: 4 × 4
#treatment W.MLD wt.var wt.sd
#<fct>        <dbl>  <dbl> <dbl>
#1 HF         7.94   3.76   1.94
#2 HN         6.52   6.08   2.47
#3 LF         11.1   6.51   2.55
#4 LN         8.75   18.6   4.32


HF_barplot <- 
  ggplot(HF.subset, aes(fill=colony_id, y=RO, x=lunar_day)) + 
  geom_bar(position="stack", stat="identity") +
  scale_y_continuous(limits=c(0, 400), breaks=c(0, 100, 200, 300, 400))+
  scale_x_continuous(limits=c(1, 29), breaks=c(1, 5, 10, 15, 20, 25, 29))+
  labs(y = "# of planulae released", x="Lunar day", title ="Fed 28⁰C")+
  theme_bw()+
  scale_fill_manual(values=c("darksalmon", "light salmon","coral","red", "firebrick","darkred"))+
  theme(legend.position="none")+
  theme(text=element_text(size=14,  family="serif"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.2))

HF_barplot

HF_barplot +  geom_vline(xintercept=7.94, linetype="dashed", 
                         color = "black", size=0.5) 


HN_barplot <- 
  ggplot(HN.subset, aes(fill=colony_id, y=RO, x=lunar_day)) + 
  geom_bar(position="stack", stat="identity") +
  scale_y_continuous(limits=c(0, 400), breaks=c(0, 100, 200, 300, 400))+
  scale_x_continuous(limits=c(1, 29), breaks=c(1, 5, 10, 15, 20, 25, 29))+
  labs(y = "# of planulae released", x="Lunar day", title ="Unfed 28⁰C")+
  theme_bw()+
  scale_fill_manual(values=c( "pink", "violet","plum", "hotpink","magenta","magenta4"))+
  theme(legend.position="none")+
  theme(text=element_text(size=14,  family="serif"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.2))

HN_barplot

HN_barplot +  geom_vline(xintercept=6.52, linetype="dashed", 
                         color = "black", size=0.5) 


LF_barplot <- 
  ggplot(LF.subset, aes(fill=colony_id, y=RO, x=lunar_day)) + 
  geom_bar(position="stack", stat="identity") +
  scale_y_continuous(limits=c(0, 400), breaks=c(0, 100, 200, 300, 400))+
  scale_x_continuous(limits=c(1, 29), breaks=c(1, 5, 10, 15, 20, 25, 29))+
  labs(y = "# of planulae released", x="Lunar day", title ="Fed 24⁰C")+
  theme_bw()+
  scale_fill_manual(values=c("dodgerblue4","royalblue2","royalblue4", "darkblue","blue", "navy"))+
  theme(legend.position="none")+
  theme(text=element_text(size=14,  family="serif"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.2))

LF_barplot

LF_barplot +  geom_vline(xintercept=11.1, linetype="dashed", 
                         color = "black", size=0.5) 


LN_barplot <- 
  ggplot(LN.subset, aes(fill=colony_id, y=RO, x=lunar_day)) + 
  geom_bar(position="stack", stat="identity") +
  scale_y_continuous(limits=c(0, 400), breaks=c(0, 100, 200, 300, 400))+
  scale_x_continuous(limits=c(1, 29), breaks=c(1, 5, 10, 15, 20, 25, 29))+
  labs(y = "# of planulae released", x="Lunar day", title ="Unfed 24⁰C")+
  theme_bw()+
  scale_fill_manual(values=c("steelblue1","skyblue", "steelblue2","steelblue3", "steelblue4","skyblue2"))+
  theme(legend.position="none")+
  theme(text=element_text(size=14,  family="serif"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.2))

LN_barplot

LN_barplot +  geom_vline(xintercept=8.75, linetype="dashed", 
                         color = "black", size=0.5) 


####4) Analysis####

WT.MLD.mm <- lmer(W.MLD~feeding_treatment*temp_treatment + (1|tank_id), data=WT.MLD)
summary(WT.MLD.mm)

#check assumptions
plot(fitted(WT.MLD.mm),residuals(WT.MLD.mm))
hist(residuals(WT.MLD.mm))
qqnorm(residuals(WT.MLD.mm))
vif(WT.MLD.mm)

WT.MLD_posthoc1 <- emmeans(WT.MLD.mm, pairwise~feeding_treatment*temp_treatment)
WT.MLD_posthoc1



####_______________####

####<Artemia Density>####

####1) Libraries####

library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(car)

####2) Data####

art_data <- read_csv("art_data.csv")


#class check

summary(art_data)

art_data$date <- as.Date(art_data$date)
art_data$trial<- as.factor(art_data$trial)
art_data$time <- as.factor(art_data$time)
art_data$temp <- as.factor(art_data$temp)
art_data$replicate <- as.factor(art_data$replicate)

summary(art_data)
head(art_data)

#subsets

#trial1
art1 <- 
  art_data%>%
  filter(trial == "trial_1")

#trial2
art2 <- 
  art_data%>%
  filter(trial == "trial_2")

#trial3
art3 <- 
  art_data%>%
  filter(trial == "trial_3")


####3) Plot####

#set order for x axis
art_data$time <- factor(art_data$time, levels=c("Before", "After"))

art_boxplot <- 
  ggplot(data=art_data, aes(x=temp, y=density, fill=time)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
             size=1) +
  scale_fill_manual(values = c("cadetblue", "grey"))+
  scale_y_continuous(limits=c(0, 80), breaks=c(0, 20, 40, 60, 80))+
  labs(y = "Density (indiv/mL)", x="Temperature (°C)", title="")+
  theme_bw()+
  theme(text=element_text(size=14,  family="serif"))+
  facet_wrap(~date)

art_boxplot


####4)Analysis####

#trial1
art1_anova <- aov(density~temp*time, data=art1)
summary(art1_anova)
plot(art1_anova)
TukeyHSD(art1_anova)

#trial2
art2_anova <- aov(density~temp*time, data=art2)
summary(art2_anova)
plot(art2_anova)
TukeyHSD(art2_anova)

#trial3
art3_anova <- aov(sqrt(density)~temp*time, data=art3)
summary(art3_anova)
plot(art3_anova)
TukeyHSD(art3_anova)


#<_____END_____>#
