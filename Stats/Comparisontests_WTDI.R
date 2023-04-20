# want to load source data, plot, and perform stats
# 1) Compare distributions of 1:0 and 0:10 WT:DI ratios - Kruskal-Wallis test?
# 2) All pairwise combinations between WT:DI ratios - comparison between 1:100 and other treatment groups - Student two-sample t-test?
  

# clear console screen
cat("\014");

# only have to run these once
install.packages("tidyverse");  # install package
install.packages("zoo");
install.packages("ggsci");
install.packages("egg");

# load a bunch of packages first!
library(zoo); # moving averages
library(tidyverse);
# library(lubridate);
library(ggplot2);
library(dplyr);
library(ggsci);
library(egg);

# set working directory
setwd("/Users/jharris387/Dropbox (Personal)/Cycling_ms_share/Stats/")


# 'Source Data 3': Figure 3B of ms
my_file <- "./../Source_Data/Figure3B_Data_fixed.csv"

data <- read_csv(my_file,head=TRUE)  # read in the csv data file

# pivot into longer format
long_data <- pivot_longer(data,cols = c('replicate 1', 'replicate 2',  'replicate 3'), names_to = 'replicate', values_to = 'total viral yield')

colnames(long_data) <- c('treatment','replicate','total_viral_yield')


# kruskal.test(total_viral_yield ~ treatment, data = long_data) 
wc_test <- wilcox.test(total_viral_yield~treatment, data=long_data, na.rm=TRUE, paired=FALSE, exact=FALSE,conf.int = TRUE)
print(wc_test)

#To subset data (pull out data points of interest)
#data_subset <- data %>% subset(data$treatment == "treatment 1" | "treatment 2")




# now reduce to relevant variables
# group by the month data
df_groupby_week <- us_data %>% select(MMWR_week,end_week,age_group,covid19_deaths) %>%
  filter(age_group =="All Ages") %>%
  mutate(date=mdy(end_week)) %>%
  arrange(date) %>%
  group_by(date);

df_all_ages <- df_groupby_week %>% summarise(total_deaths_per_week = sum(covid19_deaths,na.rm = TRUE)/2);

# Deaths due to COVID-19 by week from 24 May 2020 to 19 Sept 2020: data are shifted by 21 days to estimate the shape of the incident cases from the beginning of May to the end of August.

# let's start at week 22 = 2020-05-30 up to week 36 = 2020-09-05
g1 <- ggplot(df_all_ages, aes(date, total_deaths_per_week, group = 1)) +
  geom_bar(fill="#FFDDAF", width=6, stat="identity",alpha=0.75)+
  geom_line(position = "identity",colour = "#1261A0", size = 0.75,alpha=0.9)+
  labs(x = "Month", y = "COVID-19 deaths (per week)")+
  # labs() +
  # theme()
  scale_x_date(limits = as.Date(c("2020-05-24","2020-09-19")),breaks = as.Date(c("2020-06-13","2020-07-12","2020-08-13","2020-09-13")), labels =c("May","June","July","Aug"),expand = c(0, 0))+
  scale_y_continuous(limits = c(0,10000),expand = c(0, 0))+
  theme_bw()+
  theme(plot.tag = element_text(),
        plot.margin = unit(c(1,0.5,0.5,0.5), "lines"),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12))


# now load the median age data
# Median age of COVID-19 infection (positive RT-PCR) beginning of week from 1 May 2020 to 28 Aug 2020
# extracted from Beohmer et al.: https://www.cdc.gov/mmwr/volumes/69/wr/mm6939e1.htm
my_file <- "./../Data/MMWR_Boehmeretal_medianages_4USregions_MayAug2020.csv";

us_medianage_data <- read_csv(my_file);  # read in the csv data file

g2 <- ggplot(us_medianage_data, aes(x=end_week, y=median_positive, color=Region)) +
  geom_line(size = 0.75,alpha=0.9)+
  labs(x = "Month", y = "Median Age")+
  scale_x_date(limits = as.Date(c("2020-05-01","2020-08-28")),breaks = as.Date(c("2020-05-22","2020-06-22","2020-07-22","2020-08-22")), labels =c("May","June","July","Aug"),expand = c(0, 0))+
  scale_y_continuous(limits = c(30,55),expand = c(0, 0))+
  scale_color_manual(values=c('#000000','#56B3E9','#E69F00','#009E73','#D55E00'),limits = c("US","Midwest","Northeast","South","West"))+
  theme_bw()+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "lines"),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12));

g3 <- ggarrange(g1, g2, nrow=2)
# ggsave("FigS1_medianage_changes_US.pdf", g3,height=11, width=8)
# export as .png in plot window so that text size is saved
# added panel labels A and B later in adobe illustrator and saved as .eps
