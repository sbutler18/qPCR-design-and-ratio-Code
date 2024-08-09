library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpmisc)
library(tidyverse)
library(purrr)

setwd("/Volumes/RSMASFILES2")
getwd()

CTs <- read.csv("09222023_PROBEEFF_t1t2_data.csv")

#groups the samples by start letter or original sample before dilution

head(CTs)


CTs2<- CTs %>% 
  mutate(start_letter = substr(Sample.Name, 1, 3)) 

head(CTs2)

#split versions for CT "TLC"


CTs2 %>%
  filter(Cт!="Undetermined") %>%
  filter(Target.Name=="TLC_8") %>%
  filter(Sample.Name!="NA") %>%
  filter(Cт.SD!="NA") %>%
  ggplot( 
    aes(Quantity, y = as.numeric(Cт))
  ) + 
  geom_point(aes(color=Target.Name)) + ggtitle("TLC")+
  facet_grid(start_letter~.)+
  scale_x_continuous(trans='log10') +
  geom_smooth(method='lm', color="black") +
  
  stat_poly_eq(formula = y ~ x,
               eq.with.lhs = "italic(hat(y))~`=`~",
               aes(label = paste( ..rr.label.., ..p.value.label.., sep = "~~~")), 
               parse = TRUE, label.y = "bottom", label.x = "left", color="black",
               rr.digits = 3, , size = 2) +
  
  stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
               label.y = 0.15,
               eq.with.lhs = "italic(hat(y))~`=`~",
               eq.x.rhs = "~italic(x)",
               formula = y ~ x, parse = TRUE, size = 2,
               label.x = "left", color="black")

#all samples in one graph

CTs2 %>%
  filter(Cт!="Undetermined") %>%
  filter(Target.Name=="TLC_8") %>%
  filter(Sample.Name!="NA") %>%
  filter(Cт.SD!="NA") %>%
  ggplot( 
    aes(Quantity, y = as.numeric(Cт))
  ) + 
  geom_point(aes(color=Target.Name)) + ggtitle("TLC")+
  scale_x_continuous(trans='log10') +
  geom_smooth(method='lm', color="black") +
  
  stat_poly_eq(formula = y ~ x,
               eq.with.lhs = "italic(hat(y))~`=`~",
               aes(label = paste( ..rr.label.., ..p.value.label.., sep = "~~~")), 
               parse = TRUE, label.y = "bottom", label.x = "left", color="black",
               rr.digits = 3, , size = 2) +
  
  stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
               label.y = 0.15,
               eq.with.lhs = "italic(hat(y))~`=`~",
               eq.x.rhs = "~italic(x)",
               formula = y ~ x, parse = TRUE, size = 2,
               label.x = "left", color="black")



#Split version for CAM

CTs2 %>%
  filter(Cт!="Undetermined") %>%
  filter(Target.Name=="CAM_4") %>%
  filter(Sample.Name!="NA") %>%
  filter(Cт.SD!="NA") %>%
  ggplot( 
    aes(Quantity, y = as.numeric(Cт))
  ) + 
  geom_point(aes(color=Target.Name)) + ggtitle("TLC")+
  facet_grid(start_letter~.)+
  scale_x_continuous(trans='log10') +
  geom_smooth(method='lm', color="black") +
  
  stat_poly_eq(formula = y ~ x,
               eq.with.lhs = "italic(hat(y))~`=`~",
               aes(label = paste( ..rr.label.., ..p.value.label.., sep = "~~~")), 
               parse = TRUE, label.y = "bottom", label.x = "left", color="black",
               rr.digits = 3, , size = 2) +
  
  stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
               label.y = 0.15,
               eq.with.lhs = "italic(hat(y))~`=`~",
               eq.x.rhs = "~italic(x)",
               formula = y ~ x, parse = TRUE, size = 2,
               label.x = "left", color="black")

CTs2 %>%
  filter(Cт!="Undetermined") %>%
  filter(Target.Name=="CAM_4") %>%
  filter(Sample.Name!="NA") %>%
  filter(Cт.SD!="NA") %>%
  ggplot( 
    aes(Quantity, y = as.numeric(Cт))
  ) + 
  geom_point(aes(color=Target.Name)) + ggtitle("TLC")+
  scale_x_continuous(trans='log10') +
  geom_smooth(method='lm', color="black") +
  
  stat_poly_eq(formula = y ~ x,
               eq.with.lhs = "italic(hat(y))~`=`~",
               aes(label = paste( ..rr.label.., ..p.value.label.., sep = "~~~")), 
               parse = TRUE, label.y = "bottom", label.x = "left", color="black",
               rr.digits = 3, , size = 2) +
  
  stat_poly_eq(aes(label = paste(..eq.label.., sep = "~~~")), 
               label.y = 0.15,
               eq.with.lhs = "italic(hat(y))~`=`~",
               eq.x.rhs = "~italic(x)",
               formula = y ~ x, parse = TRUE, size = 2,
               label.x = "left", color="black")

#slopes? I could not get these to work, just get NA for slope


CTs %>%
filter(Ct!="Undetermined") %>%
filter(Target=="CAM") %>%
filter(Sample_Names!="NA") %>%
filter(CT_SD!="NA") %>%
  group_by(Sample_Names) %>% 
  mutate(rownum = row_number()) %>% 
  summarise(slope = lm(as.numeric(Ct) ~ log10(Dilution))$coefficients[2]) 



CTs %>%
filter(Ct!="Undetermined") %>%
filter(Target=="TLC") %>%
filter(Sample_Names!="NA") %>%
filter(CT_SD!="NA") %>%
  group_by(Sample_Names) %>% 
  mutate(rownum = row_number()) %>% 
 summarise(slope = lm(as.numeric(Ct) ~ log10(Dilution))$coefficients[2]) 
