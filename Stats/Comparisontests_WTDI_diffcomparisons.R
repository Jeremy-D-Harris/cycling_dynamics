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

# data <- read_csv(my_file,head=TRUE)  # read in the csv data file
data <- read_csv(my_file)  # read in the csv data file
colnames(data) <- c('treatment','replicate 1', 'replicate 2',  'replicate 3')
# view(data)

# pivot into longer format
long_data <- pivot_longer(data,cols = c('replicate 1', 'replicate 2',  'replicate 3'), names_to = 'replicate', values_to = 'total viral yield')
colnames(long_data) <- c('treatment','replicate','total_viral_yield')


# kruskal.test(total_viral_yield ~ treatment, data = long_data) 
wilcox_test_Figure3B_data <- wilcox.test(total_viral_yield~treatment, data=long_data, na.rm=TRUE, paired=FALSE, exact=FALSE,conf.int = TRUE)
print(wilcox_test_Figure3B_data)

#
graph_Figure3B_data <- ggplot(data = long_data, aes(x=treatment, y=total_viral_yield,fill=replicate)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_line(size = 0.75,alpha=0.9)+
  labs(x = "Ratio WT to DI virus inoculum", y = "Total viral yield (GE/mL)")+
  theme_bw()+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "lines"),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12));
  # geom_bar(fill="#FFDDAF", width=0.5, stat="identity",alpha=0.75)
  # geom_line(position = "identity",colour = "#1261A0", size = 0.75,alpha=0.9)
  
ggarrange(graph_Figure3B_data, nrow=1)

ggsave("graph_Figure3B_data.pdf", graph_Figure3B_data,height=4.5, width=6)


###############
# now compare WT:DI = 1:100 with:
# WT:DI = 1:1
# WT:DI = 1:1
# WT:DI = 1:0.1
# WT:DI = 1:0

# 'Source Data 3': viral yield comparisons -> 'yield_comparisons'
my_file_comparisons <- "./../Source_Data/yield_comparisons.csv"

data_comparisons <- read_csv(my_file_comparisons,show_col_types = FALSE)  # read in the csv data file
colnames(data_comparisons) <- c('treatment','replicate 1', 'replicate 2',  'replicate 3')
# show_col_types = FALSE - this was due to an error, since cell 1,1 is empty -> gets filled with '...1'
# 1. WT:DI = 1:0
# 2. WT:DI = 1:0.1
# 3. WT:DI = 1:1
# 4. WT:DI = 1:10
# 5. WT:DI = 1:100


# compare treatments with WT
long_data_comparisons <- pivot_longer(data_comparisons,cols = c('replicate 1', 'replicate 2',  'replicate 3'), names_to = 'replicate', values_to = 'total viral yield')
colnames(long_data_comparisons) <- c('treatment','replicate','total_viral_yield')


kuskal_wallis_multiple_means <- kruskal.test(total_viral_yield~treatment, data = long_data_comparisons)
# kuskal_wallis_multiple_means <-  long_data %>% kruskal.test(total_viral_yield~treatment)
print(kuskal_wallis_multiple_means)


graph_compare_withWT <- ggplot(data = long_data_comparisons, aes(x=treatment, y=total_viral_yield,fill=replicate)) +
  geom_bar(stat="identity", position=position_dodge())+
  # geom_line(size = 0.75,alpha=0.9)+
  labs(x = "Ratio WT to DI virus inoculum", y = "Total viral yield (GE/mL)")+
  theme_bw()+
  scale_fill_brewer(palette="Blues")+
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "lines"),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12));
# geom_bar(fill="#FFDDAF", width=0.5, stat="identity",alpha=0.75)
# geom_line(position = "identity",colour = "#1261A0", size = 0.75,alpha=0.9)

ggarrange(graph_compare_withWT, nrow=1)

ggsave("graph_compare_withWT.pdf", graph_compare_withWT,height=4.5, width=8)

# data_combine
# 1. WT:DI = 1:0
# 2. WT:DI = 0:10
# 3. WT:DI = 1:0
# 4. WT:DI = 1:0.1
# 5. WT:DI = 1:1
# 6. WT:DI = 1:10
# 7. WT:DI = 1:100



# combine treatments with and without WT virus
data_combine <-merge(x=data,y=data_comparisons,all= T)
print(data_combine)
class(data_combine)


# pivot into longer format
long_data_combine <- pivot_longer(data_combine,cols = c('replicate 1', 'replicate 2',  'replicate 3'), names_to = 'replicate', values_to = 'total viral yield')
colnames(long_data_combine) <- c('treatment','replicate','total_viral_yield')




kuskal_wallis_multiple_means_combine <- kruskal.test(total_viral_yield~treatment, data = long_data_combine)
# kuskal_wallis_multiple_means <-  long_data %>% kruskal.test(total_viral_yield~treatment)
print(kuskal_wallis_multiple_means_combine)


graph_compare_combine <- ggplot(data = long_data_combine, aes(x=treatment, y=total_viral_yield,group=replicate,fill=treatment)) +
  geom_bar(stat="identity", colour="black",position=position_dodge(),,show.legend = FALSE)+
  scale_y_log10()+
  # geom_line(size = 0.75,alpha=0.9)+
  labs(x = "Ratio WT to DI virus inoculum", y = "Total viral yield (GE/mL)")+
  theme_bw()+
  # scale_fill_brewer(palette="Blues")+
  scale_fill_manual(values=c("#63ABBD","#0F2080","#BDB8AD","#EE442F","#A95AA1","#5F1A4A","#63ABBD"))
  theme(plot.margin = unit(c(1,0.5,0.5,0.5), "lines"),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 12));
# geom_bar(fill="#FFDDAF", width=0.5, stat="identity",alpha=0.75)
# geom_line(position = "identity",colour = "#1261A0", size = 0.75,alpha=0.9)

ggarrange(graph_compare_combine, nrow=1)

ggsave("graph_compare_combine.pdf", graph_compare_combine,height=4.5, width=11)

# ----- for reference -----
# data_combine
# 1. WT:DI = 1:0
# 2. WT:DI = 0:10
# 3. WT:DI = 1:0
# 4. WT:DI = 1:0.1
# 5. WT:DI = 1:1
# 6. WT:DI = 1:10
# 7. WT:DI = 1:100


## combine subsets of treatments for pair-wise comparisons
# compare_WT_only_treatments <- data_combine[],subset = rownames(data_comparisons) == c(" WT 1 DI 0"))

# combine treatments with and without WT virus
compare_WT_only_treatments <-merge(x=data_combine[1,],y=data_combine[7,],all= T)
print(compare_WT_only_treatments)

long_compare_WT_only_treatments <- pivot_longer(compare_WT_only_treatments,cols = c('replicate 1', 'replicate 2',  'replicate 3'), names_to = 'replicate', values_to = 'total viral yield')
colnames(long_compare_WT_only_treatments) <- c('treatment','replicate','total_viral_yield')

wilcox_test_compare_WT_only_treatments <- wilcox.test(total_viral_yield~treatment, data=long_compare_WT_only_treatments, na.rm=TRUE, paired=FALSE, exact=FALSE,conf.int = TRUE)
print(wilcox_test_compare_WT_only_treatments)


# compare DI = 10 (alone) to WT only in the first group
# data_combine[1,]
# data_combine[6,]
compare_DI_only_WT_only_treatments <-merge(x=data_combine[1,],y=data_combine[6,],all= T)
print(compare_DI_only_WT_only_treatments)

long_compare_DI_only_WT_only_treatments <- pivot_longer(compare_DI_only_WT_only_treatments,cols = c('replicate 1', 'replicate 2',  'replicate 3'), names_to = 'replicate', values_to = 'total viral yield')
colnames(long_compare_DI_only_WT_only_treatments) <- c('treatment','replicate','total_viral_yield')

wilcox_test_compare_DI_only_WT_only_treatments <- wilcox.test(total_viral_yield~treatment, data=long_compare_DI_only_WT_only_treatments, na.rm=TRUE, paired=FALSE, exact=FALSE,conf.int = TRUE)
print(wilcox_test_compare_DI_only_WT_only_treatments)















long_data_two_vs_one <- pivot_longer(long_data_combine,cols = c('replicate 1', 'replicate 2',  'replicate 3'), names_to = 'replicate', values_to = 'total viral yield')
colnames(long_data_two_vs_one) <- c('treatment','replicate','total_viral_yield')
# view(long_data_five_vs_one)

# kruskal.test(total_viral_yield ~ treatment, data = long_data) 
wilcox_test_two_vs_one <- wilcox.test(total_viral_yield~treatment, data=long_data_two_vs_one, na.rm=TRUE, paired=FALSE, exact=FALSE,conf.int = TRUE)
print(wilcox_test_two_vs_one)

# # here's student's t-test for comparison!
# t_test_five_vs_one <- t.test(total_viral_yield~treatment,data=long_data_five_vs_one,
#        alternative = c("two.sided", "less", "greater"),
#        mu = 0, paired = FALSE, var.equal = FALSE,
#        conf.level = 0.95)
# print(t_test_five_vs_one)


## 5 vs. 2
treatment_five_vs_two <- subset(data_comparisons,subset = rownames(data_comparisons) %in% c('2','5'))

long_data_five_vs_two <- pivot_longer(treatment_five_vs_two,cols = c('replicate 1', 'replicate 2',  'replicate 3'), names_to = 'replicate', values_to = 'total viral yield')
colnames(long_data_five_vs_two) <- c('treatment','replicate','total_viral_yield')

wilcox_test_five_vs_two <- wilcox.test(total_viral_yield~treatment, data=long_data_five_vs_two, na.rm=TRUE, paired=FALSE, exact=FALSE,conf.int = TRUE)
print(wilcox_test_five_vs_two)

t_test_five_vs_two <- t.test(total_viral_yield~treatment,data=long_data_five_vs_two,
                             alternative = c("two.sided", "less", "greater"),
                             mu = 0, paired = FALSE, var.equal = FALSE,
                             conf.level = 0.95)
print(t_test_five_vs_two)


## 5 vs. 3
treatment_five_vs_three <- subset(data_comparisons,subset = rownames(data_comparisons) %in% c('3','5'))

long_data_five_vs_three <- pivot_longer(treatment_five_vs_three,cols = c('replicate 1', 'replicate 2',  'replicate 3'), names_to = 'replicate', values_to = 'total viral yield')
colnames(long_data_five_vs_three) <- c('treatment','replicate','total_viral_yield')

wilcox_test_five_vs_three <- wilcox.test(total_viral_yield~treatment, data=long_data_five_vs_three, na.rm=TRUE, paired=FALSE, exact=FALSE,conf.int = TRUE)
print(wilcox_test_five_vs_three)

## 5 vs. 4
treatment_five_vs_four <- subset(data_comparisons,subset = rownames(data_comparisons) %in% c('4','5'))

long_data_five_vs_four <- pivot_longer(treatment_five_vs_four,cols = c('replicate 1', 'replicate 2',  'replicate 3'), names_to = 'replicate', values_to = 'total viral yield')
colnames(long_data_five_vs_four) <- c('treatment','replicate','total_viral_yield')

wilcox_test_five_vs_four <- wilcox.test(total_viral_yield~treatment, data=long_data_five_vs_four, na.rm=TRUE, paired=FALSE, exact=FALSE,conf.int = TRUE)
print(wilcox_test_five_vs_four)

## 5 vs. 5
treatment_five_vs_five <- subset(data_comparisons,subset = rownames(data_comparisons) %in% c('5','5')) # for sanity
treatment_five_vs_five <- treatment_five_vs_five[rep(1, 2),]
treatment_five_vs_five[2,1] <- "Pretend Different" 

long_data_five_vs_five <- pivot_longer(treatment_five_vs_five,cols = c('replicate 1', 'replicate 2',  'replicate 3'), names_to = 'replicate', values_to = 'total viral yield')
colnames(long_data_five_vs_five) <- c('treatment','replicate','total_viral_yield')

wilcox_test_five_vs_five <- wilcox.test(total_viral_yield~treatment, data=long_data_five_vs_five, na.rm=TRUE, paired=FALSE, exact=FALSE,conf.int = TRUE)
print(long_data_five_vs_five)

# pivot each into longer format


# g2 <- ggplot(us_medianage_data, aes(x=end_week, y=median_positive, color=Region)) +
#   geom_line(size = 0.75,alpha=0.9)+
#   labs(x = "Month", y = "Median Age")+
#   scale_x_date(limits = as.Date(c("2020-05-01","2020-08-28")),breaks = as.Date(c("2020-05-22","2020-06-22","2020-07-22","2020-08-22")), labels =c("May","June","July","Aug"),expand = c(0, 0))+
#   scale_y_continuous(limits = c(30,55),expand = c(0, 0))+
#   scale_color_manual(values=c('#000000','#56B3E9','#E69F00','#009E73','#D55E00'),limits = c("US","Midwest","Northeast","South","West"))+
#   theme_bw()+
#   theme(plot.margin = unit(c(1,0.5,0.5,0.5), "lines"),
#         axis.text = element_text(size = 11),
#         axis.title = element_text(size = 12));

# g3 <- ggarrange(g1, g2, nrow=2)
# ggsave("FigS1_medianage_changes_US.pdf", g3,height=11, width=8)
# export as .png in plot window so that text size is saved
# added panel labels A and B later in adobe illustrator and saved as .eps
