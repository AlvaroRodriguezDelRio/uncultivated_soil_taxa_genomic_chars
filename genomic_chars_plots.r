
library(dplyr)
library(data.table)
library(stringr)
library(dplyr)
library(tidyverse)
library(ggridges)
library(ggpubr)
library(patchwork)
library(ggplot2)
library(ggnewscale)
library(purrr)
library(ggfortify)
library(randomForest)
library(caTools)
library(corrplot)


data = read.csv('supplementary_data.csv',
                header = T,check.names = FALSE)

######
# correlations of items
######

d_cor = data %>% 
  dplyr::select(!c(Uncultivated,proportion_RS_genomes))


c = cor(d_cor,use = "pairwise.complete.obs")
pval <- psych::corr.test(d_cor, adjust="none")$p

corrplot(c, method = "circle", type = "upper",p.mat=pval,
         sig.level = 0.05,  insig = "blank",diag = F,
         order = 'hclust',tl.col = 'black')


####
# see difference between completely uncultivated and the rest
###

# remove outgroups 
data_all_boxplot = data %>% 
  dplyr::select(!proportion_RS_genomes) %>% 
  mutate(`genome size` = `genome size` / 1000000) %>% 
  gather(key = 'item',value = 'value',-c(Uncultivated)) %>% 
  filter(!(item == "doubling time" & value > 40)) %>% 
  filter(!(item == "biosynthetic clusters" & value > 20)) %>% 
  filter(!(item == "salinity optimum" & value > 5)) 


p_box = ggplot(data_all_boxplot,
               aes(x = Uncultivated,y = value,
                   fill = Uncultivated))+
  geom_boxplot(outliers = F,alpha = 0.7)+
  facet_wrap(~item,scales = 'free_x',ncol = 2)+
  theme_classic()+
  stat_compare_means(aes(group=Uncultivated),label = "p.signif", method = "wilcox.test",
                     ref.group = "Cultivated representative",hide.ns = TRUE, tip.length = 0,paired = F)+
  coord_flip()+
  scale_fill_manual(
    values = c("Cultivated representative" = "orange",
               "Uncultivated" = "lightblue"))+
  scale_color_manual(
    values = c("Cultivated representative" = "orange",
               "Uncultivated" = "lightblue"))+
  geom_jitter(aes(x = Uncultivated,y = value, color = Uncultivated),
              position = position_jitter(0.2),alpha = 0.02)+
  theme(legend.position = 'none')+
  theme(
    axis.text = element_text(size = 14),      # x & y tick labels
    axis.title = element_text(size = 16),       # axis titles (optional)
    strip.text = element_text(size = 16)
  )+
  xlab('')


p_box

# Wilcoxon tests results 
data_all_boxplot = data %>% 
  dplyr::select(!proportion_RS_genomes) %>% 
  gather(key = 'item',value = 'value',-c(Uncultivated)) 

pvals = data_all_boxplot %>% 
  filter(!is.na(value)) %>% 
  group_by(item) %>%
  reframe(
    p_value = wilcox.test(value ~ Uncultivated)$p.value,
  ) %>%
  print(n = 100)


mean_cult = data_all_boxplot %>% 
  filter(Uncultivated != "Uncultivated") %>% 
  group_by(item) %>%
  summarise(mean_cult = mean(value)) 


mean_uncult = data_all_boxplot %>% 
  filter(Uncultivated == "Uncultivated") %>% 
  group_by(item) %>%
  summarise(mean_uncult = mean(value)) 


pvals %>% 
  left_join(mean_cult,by = 'item') %>% 
  left_join(mean_uncult,by = 'item') %>% 
  print(n = 100)


# 10% threshold for uncultivated 
data_all_boxplot_10 = data %>% 
  mutate(Uncultivated = case_when(proportion_RS_genomes < 0.1 ~'Uncultivated',
                                  .default = 'Cultivated representative')) %>% 
  dplyr::select(!proportion_RS_genomes) %>% 
  mutate(`genome size` = `genome size` / 1000000) %>% 
  gather(key = 'item',value = 'value',-c(Uncultivated)) %>% 
  filter(!(item == "doubling time" & value > 40)) %>% 
  filter(!(item == "biosynthetic clusters" & value > 20)) %>% 
  filter(!(item == "salinity optimum" & value > 5)) 


p_box_10 = ggplot(data_all_boxplot_10,
               aes(x = Uncultivated,y = value,
                   fill = Uncultivated))+
  geom_boxplot(outliers = F,alpha = 0.7)+
  facet_wrap(~item,scales = 'free_x',ncol = 2)+
  theme_classic()+
  stat_compare_means(aes(group=Uncultivated),label = "p.signif", method = "wilcox.test",
                     ref.group = "Cultivated representative",hide.ns = TRUE, tip.length = 0,paired = F)+
  coord_flip()+
  scale_fill_manual(
    values = c("Cultivated representative" = "orange",
               "Uncultivated" = "lightblue"))+
  scale_color_manual(
    values = c("Cultivated representative" = "orange",
               "Uncultivated" = "lightblue"))+
  geom_jitter(aes(x = Uncultivated,y = value, color = Uncultivated),
              position = position_jitter(0.2),alpha = 0.02)+
  theme(legend.position = 'none')+
  theme(
    axis.text = element_text(size = 14),      # x & y tick labels
    axis.title = element_text(size = 16),       # axis titles (optional)
    strip.text = element_text(size = 16)
  )+
  xlab('')


p_box_10


# violin plot 
ggplot(data_all_boxplot,
       aes(x = Uncultivated,y = value,
           fill = Uncultivated))+
  geom_violin(outliers = F,alpha = 0.7)+
  facet_wrap(~item,scales = 'free_x',ncol = 2)+
  theme_classic()+
  stat_compare_means(aes(group=Uncultivated),label = "p.signif", method = "wilcox.test",
                     ref.group = "Cultivated representative",hide.ns = TRUE, tip.length = 0,paired = F)+
  coord_flip()+
  scale_fill_manual(
    values = c("Cultivated representative" = "orange",
               "Uncultivated" = "lightblue"))+
  scale_color_manual(
    values = c("Cultivated representative" = "orange",
               "Uncultivated" = "lightblue"))+
  geom_jitter(aes(x = Uncultivated,y = value, color = Uncultivated),
              position = position_jitter(0.2),alpha = 0.02)+
  theme(legend.position = 'none')+
  theme(
    axis.text = element_text(size = 14),      # x & y tick labels
    axis.title = element_text(size = 16),       # axis titles (optional)
    strip.text = element_text(size = 16)
  )+
  xlab('')


###
# distribution proportion uncultivated per genus
###

ggplot(data)+
  geom_histogram(aes(x = proportion_RS_genomes))+
  theme_classic()+
  xlab('Proportion of reference genomes per genus')


####
# correlation doubling time with the rest
###

data %>% 
  gather(key = 'item',value = 'value',-c(Uncultivated,`doubling time`)) %>% 
  group_by(item) %>% 
  summarise(c = cor.test(value,`doubling time`,method = 'spearman')$estimate,
            p = cor.test(value,`doubling time`,method = 'spearman')$p.value) %>% 
  arrange((p)) %>% 
  print(n = 100)


#########
# ML, most important variables for prediction 
#########

data_ml = data
tips = data_ml$tip

proportion_RS_genomes = data_ml$proportion_RS_genomes
data_ml$proportion_RS_genomes = NULL


# model with all 
set.seed(1)
data_ml$Uncultivated <- as.factor(data_ml$Uncultivated)
original_names <- names(data_ml)
names(data_ml) <- make.names(original_names)

rf_model <- randomForest(Uncultivated ~ ., data = data_ml, ntree = 500, 
                         mtry = 3, importance = TRUE)

print(rf_model)
imp = as.data.frame(importance(rf_model))

imp$item = rownames(imp)

p_ml = ggplot(imp)+
  geom_histogram(aes(y = `MeanDecreaseAccuracy`, 
                     x = reorder(item, `MeanDecreaseAccuracy`)),
                 stat = 'identity',alpha = 0.5)+
  coord_flip()+
  theme_classic()+
  xlab('')+
  ylab('Mean decrease\naccuracy')+
  theme(
    axis.text = element_text(size = 14),      # x & y tick labels
    axis.title = element_text(size = 16),       # axis titles (optional)
    strip.text = element_text(size = 16)
  )

p_ml

#########
# ML, most important variables for prediction, 10% threshold
#########

data_ml = data
tips = data_ml$tip


data_ml = data_ml %>% 
  mutate(Uncultivated = case_when(proportion_RS_genomes < 0.1 ~'Uncultivated',
                                  .default = 'Cultivated representative')) 

proportion_RS_genomes = data_ml$proportion_RS_genomes
data_ml$proportion_RS_genomes = NULL


# model with all 
set.seed(1)
data_ml$Uncultivated <- as.factor(data_ml$Uncultivated)


original_names <- names(data_ml)
names(data_ml) <- make.names(original_names)

rf_model <- randomForest(Uncultivated ~ ., data = data_ml, ntree = 500, 
                         mtry = 3, importance = TRUE)

print(rf_model)
imp = as.data.frame(importance(rf_model))

imp$item = rownames(imp)
head(imp,n = 100)

p_ml_10 = ggplot(imp)+
  geom_histogram(aes(y = `MeanDecreaseAccuracy`, 
                     x = reorder(item, `MeanDecreaseAccuracy`)),
                 stat = 'identity',alpha = 0.5)+
  coord_flip()+
  theme_classic()+
  xlab('')+
  ylab('Mean decrease\naccuracy')+
  theme(
    axis.text = element_text(size = 14),      # x & y tick labels
    axis.title = element_text(size = 16),       # axis titles (optional)
    strip.text = element_text(size = 16)
  )

p_ml_10

###
#combine
###


layout <- "
AAAAB
"

p_box + p_ml + plot_layout(design = layout)+
  plot_annotation(tag_levels = 'A')

layout <- "
AAAAB
"

p_box_10 + p_ml_10 + plot_layout(design = layout)+
  plot_annotation(tag_levels = 'A')

