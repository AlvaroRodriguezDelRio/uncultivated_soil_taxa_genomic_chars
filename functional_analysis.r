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
library(ggrepel)


kos = read.csv('data_functional.csv')


###
# Plot coverages
###



dp = kos %>% 
  dplyr::select(general,kp_desc,mean_cultivated,mean_uncultivated) %>% 
  unique() %>% 
  group_by(general,kp_desc) %>% 
  summarise(mu = mean(mean_uncultivated),
            mc = mean(mean_cultivated)) %>% 
  mutate(enriched = case_when(mu > mc~"Uncultivated",
                              .default = 'Cultivated')) %>% 
  filter(!is.na(general)) 


pcorr = ggplot(dp)+
  geom_point(aes(x = mu * 100, y = mc * 100, color = enriched),size = 3)+
  geom_abline()+ 
  geom_text_repel(data = dp %>% 
                    filter(mu -mc>0.04 |
                             mu -mc < -0.1),
                  aes(x = mu * 100,y = mc * 100, label = kp_desc),
                  box.padding = 1,
                  nudge_y = 1,
                  #segment.curvature = -0.1,
                  segment.ncp = 2,
                  #segment.angle = 50,
                  arrow = arrow(length = unit(0.015, "npc")),
                  seed = 100
  )+
  scale_colour_manual(
    values = c("Cultivated" = "orange",
               "Uncultivated" = "lightblue"))+
  ylab("% genomes (cultivated)")+
  xlab("% genomes (uncultivated)")+
  theme_classic()+
  theme(legend.position = 'none')


pcorr


#########
####
# tsne plot 
####
##########

ref_genomes = read.table("lin2n_refseq.tab",sep = '\t')
names(ref_genomes) = c("tip","genome_orig","n_genomes_gtdb")
head(ref_genomes)

ref_genomes = ref_genomes %>% 
  group_by(tip) %>% 
  mutate(tot = sum(n_genomes_gtdb)) %>% 
  mutate(proportion_RS_genomes = n_genomes_gtdb / tot) %>% 
  filter(genome_orig == 'RS') %>% 
  dplyr::select(tip,proportion_RS_genomes)


size = read.table('len_per_lin.tab',sep = '\t')
names(size) = c('tip','mean_size','median_size','n_genomes')

######
# ko tsne
######

tsne = read.csv('tsne_results.kp.genus.csv')

data = tsne %>% 
  dplyr::left_join(ref_genomes, by = 'tip') %>% 
  dplyr::left_join(size, by = 'tip') 


# quit unknown species
data = data %>% 
  filter(!grepl("_$", tip)) %>% 
  filter(!is.na(tSNE1))


head(data)
data$taxonomic_level <- sapply(strsplit(data$tip, ";"), function(x) substr(tail(x, 1), 1, 1))


data = data%>% 
  mutate(uncultivated = case_when(proportion_RS_genomes == 0~'Uncultivated',
                                  .default = 'Cultivated representative'))

###
# plot
###

# color by phylum 
data = data %>% 
  mutate(phylum = tstrsplit(tip, "\\;")[[2]])


# phylum colors only if certain number   
npp = data %>% 
  group_by(phylum) %>% 
  mutate(n = n()) %>% 
  mutate(name = case_when(n > 50 ~ phylum,
                          .default = 'Other'))



# coloring per phylum
npp$name = as.factor(npp$name)
default_colors <- scale_color_hue()$palette(length(unique(npp$name)))
species_levels <- levels(npp$name)
label_colors <- setNames(default_colors, species_levels)
label_colors["Other"] <- "grey"


ggplot(npp)+
  geom_point(aes(x = tSNE1     , y = tSNE2, color = name,
                 size = mean_size /  1000000),alpha = 0.6)+
  scale_shape_manual(values=c(15, 19))+
  scale_color_manual(values = label_colors)+
  labs(color = "Phylum", size = "Genome size (Mbp)")+
  guides(
    color = guide_legend(override.aes = list(size = 5)))+
  theme_classic()


# coloring by cultivated / uncultivated
p1 = ggplot(data)+
  geom_point(aes(x = tSNE1     , y = tSNE2, color = uncultivated,
                 size = mean_size /1000000),alpha = 0.6)+
  scale_colour_manual(
    values = c("Cultivated representative" = "orange",
               "Uncultivated" = "lightblue"))+
  labs(color = "Cultivated", size = "Genome size (Mbp)")+
  guides(
    color = guide_legend(override.aes = list(size = 5)))+
  theme_classic()


p1 
p2 = ggplot(data,aes(x = tSNE1, color = uncultivated,fill = uncultivated))+
  geom_density(alpha = 0.5)+
  scale_colour_manual(
    values = c("Cultivated representative" = "orange",
               "Uncultivated" = "lightblue"))+
  scale_fill_manual(
    values = c("Cultivated representative" = "orange",
               "Uncultivated" = "lightblue"))+
  theme_classic()+
  theme(legend.position = 'none')



p3 = ggplot(data,aes(x = tSNE2, color = uncultivated,fill = uncultivated))+
  geom_density(alpha = 0.5)+
  scale_colour_manual(
    values = c("Cultivated representative" = "orange",
               "Uncultivated" = "lightblue"))+
  scale_fill_manual(
    values = c("Cultivated representative" = "orange",
               "Uncultivated" = "lightblue"))+
  coord_flip()+
  theme_classic()+
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p3

layout <- "
AAAAAA#
BBBBBBC
BBBBBBC
BBBBBBC
BBBBBBC
BBBBBBC
BBBBBBC
"

ppca = p2 + p1 + p3 + plot_layout(design = layout,guides = "collect")
ppca



#########
# combine figures
##########

layout <- "
1111X44
2222344
2222355
2222355
"


ppca+ pcorr + guide_area() +
  plot_layout(design = layout,guides = "collect")+
  #plot_annotation(tag_levels = 'A')+
  theme(legend.position = "right")


