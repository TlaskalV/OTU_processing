####Knihovny####
#test fce git
# library(swirl)
library(Rmisc)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(readr)
library(xlsx)

setwd("D:/Vojta/Disk/Laborka/Data/Manualy/Skripty_R_Python/OTU_processing_test")

####Nacteni souboru####

otus <- read_delim("D:/Vojta/Disk/Laborka/Data/Manualy/Skripty_R_Python/OTU_processing_test/otus.csv", 
                   ";", escape_double = FALSE, trim_ws = TRUE)
View(otus)

###

tax <- read_delim("D:/Vojta/Disk/Laborka/Data/Manualy/Skripty_R_Python/OTU_processing_test/tax.csv", 
                  ";", escape_double = FALSE, trim_ws = TRUE)
View(tax)

###

samples <- read_delim("D:/Vojta/Disk/Laborka/Data/Manualy/Skripty_R_Python/OTU_processing_test/samples.csv", 
                      ";", escape_double = FALSE, trim_ws = TRUE)
View(samples)

###

div <- read_delim("D:/Vojta/Disk/Laborka/Data/Manualy/Skripty_R_Python/OTU_processing_test/div.csv", 
                  ";", escape_double = FALSE, trim_ws = TRUE)
View(div)

####Uprava vstupnich tabulek####

colnames(otus)[which(names(otus) == "Group name:")] <- "cluster_name"
colnames(tax)[which(names(tax) == "SEQ TITLE")] <- "cluster_name"
colnames(tax)[which(names(tax) == "Coverage[%]")] <- "Coverage"
colnames(tax)[which(names(tax) == "Similarity[%]")] <- "Similarity"
colnames(div)[which(names(div) == "X1")] <- "sample"
samples$year <- as.factor(samples$year)


otus_w <- tbl_df(otus)
tax_w <- tbl_df(tax)
samples_w <- tbl_df(samples)
div_w <- tbl_df(div)
head(otus_w)
head(tax_w)
head(samples_w)
head(div_w)


titles_tax <- select(tax_w, `cluster_name`) %>% print
titles_otus <- select(otus_w, `cluster_name`) %>% print

tax_w <- mutate(tax_w, 
                     `cluster_name` = sapply(strsplit(tax_w$`cluster_name`, 
                                                               split = '|', fixed = TRUE), function(x) (x[2]))) %>%
  print

otus_w <- mutate(otus_w, 
                `cluster_name` = sapply(strsplit(otus_w$`cluster_name`, 
                                              split = '|', fixed = TRUE), function(x) (x[2]))) %>%
  print

table_all <- left_join(otus_w, tax_w, by = "cluster_name") %>%   # kombinace tabulek, v pripade vice radku v taxonomii oproti poctu clusteru v otus vezme pouze prislusne radky z taxonomie
  filter(superkingdom == "Bacteria") %>%  # filtrace pouze bakterii, vyradi NO HIT a Archea a Eukaryota 
  gather(sample, count, contains("VOJ"), na.rm = FALSE) %>%  # preformatovani vzorku
  full_join(samples_w, by = "sample") %>%
  select(cluster_name, Accession, Similarity, Coverage, phylum, class, order, genus, sample, count, tree, year) %>%
  replace_na(list(phylum = "unclass", class = "unclass", order = "unclass", genus = "unclass")) %>% # nahrazeni NAs
  print

titles_bact_only <- distinct(table_all, `cluster_name`) %>% print
write.table(titles_bact_only, file = "titles_bact_only.csv", row.names = FALSE, col.names = TRUE, dec = ".", sep = ";")

sequence_counts <- group_by(table_all, sample) %>%   # pocty sekvenci ve vzorcich
  summarize(sum(count)) %>%
  print
write.table(sequence_counts, file = "sequence_counts.csv", row.names = FALSE, col.names = TRUE, dec = ".", sep = ";")

###RODY####

table_genus <- group_by(table_all, genus, sample) %>% # tabulka procenta, serazeno
  summarize(seq_number = sum(count)) %>%
  ungroup %>%
  group_by(sample) %>%
  mutate(per = (seq_number / sum(seq_number) * 100)) %>%
  ungroup %>%
  arrange(desc(per), genus) %>%
  print

treshold_genus <- group_by(table_genus, genus) %>%   # pocty sekvenci ve vzorcich
  filter(per > 3) %>%                                   # TRESHOLD
  summarise(treshold_count = n()) %>%
  print

###RODY TREE####

table_genus_tree <- left_join(table_genus, treshold_genus, by = "genus") %>%
  replace_na(list(treshold_count = 0)) %>% # nahrazeni NA nulou
  arrange(desc(treshold_count), genus) %>%
  full_join(samples_w, by = "sample") %>%
  select(-id) %>%
  print %>%
  group_by(genus, tree) %>% # podle ceho grupovat
  mutate(avg = mean(per)) %>% # prumery
  print 

stat_genus_tree <- summarySE(table_genus_tree, measurevar = "per", groupvars = c("genus", "tree"), conf.interval = 0.95) # sumarni statistika

graf_genus_abundant_tree <- filter(table_genus_tree, treshold_count > 1) %>% # vybrat kolikrat mely nad treshold
  group_by(genus, tree) %>%
  summarise(avg = mean(per)) %>%
  print

graf_genus_rare_tree <- filter(table_genus_tree, treshold_count <= 1) %>% # vybrat kolikrat mely pod treshold
  group_by(genus, tree) %>%
  summarise(avg = mean(per)) %>%
  print %>%
  ungroup %>%
  group_by(tree) %>%
  mutate(sum_per_tree = sum(avg)) %>% # soucet v ramci stromu pro other
  print %>%
  summarise(avg = sum(avg)) %>% # soucet v ramci stromu pro other
  print

graf_genera_tree_final <- bind_rows(graf_genus_abundant_tree, graf_genus_rare_tree) %>% # spojeni abundantnich a other
  replace_na(list(genus = " other")) %>% # pojmenovani other
  filter(genus == " other" | genus == "unclass") %>% # filtr other a unclass
  group_by(tree) %>%
  summarise("avg" = sum(avg)) %>% # soucet pro other a unclass
  ungroup %>%
  bind_rows(filter(graf_genus_abundant_tree, genus != "unclass")) %>% # spojeni souctu s tabulkou abundatnich bez unclass
  replace_na(list(genus = " other")) %>% # opet nutne pojmenovat, mezera pro lepsi trideni v ggplotu
  left_join(stat_genus_tree, by = c("genus", "tree")) %>% # spojeni se statistikou
  print
write.table(graf_genera_tree_final, file = "graf_genera_tree_final.csv", row.names = FALSE, col.names = TRUE, dec = ".", sep = ";")
# write.xlsx(graf_genera_year_final, "graf_genera_year_final.xlsx", sheetName="genera")

graf_genera_tree_gg <- ggplot(graf_genera_tree_final, aes(tree)) +
  geom_bar(aes(fill = genus, weight = avg), position = "fill", width = .5) +
  # scale_y_discrete (limits = "Acidobacteria", "Bacteroidetes") +
  scale_fill_brewer(palette = "Set3", type = "div") +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  labs(x = "Tree", y = "Relative abundance") 
graf_genera_tree_gg
ggsave(plot = graf_genera_tree_gg, filename = "graf_genera_tree_gg.pdf", dpi = 300, height = 210, width = 297, units = "mm")

graf_genera_tree_gg <- ggplot(graf_genera_tree_final, aes(genus, fill = tree)) +  # horizontalni graf po rodech
  geom_bar(aes(weight = avg), position = "dodge", width = .8) +
  scale_fill_brewer(palette = "Set1", type = "div", name = "Tree species", labels = c("Quercus sp.", "Abies sp.", "Picea sp.")) +
  # scale_y_continuous(labels = scales::percent) +
  geom_errorbar(aes(ymin = avg - se, ymax = avg + se),
                size = .7,    
                width = .3,
                position = position_dodge(width = .9)) +
  theme_bw() +
  theme(axis.text.y = y_label, axis.text.x = x_label, axis.title = axis_titles, plot.title = plot_title, legend.title = legend_title, legend.text = legend_text) + 
  xlab("Bacterial genera") + 
  ylab("Relative abundance (%)") + 
  ggtitle("Genera per tree") +
  coord_flip()
graf_genera_tree_gg
ggsave(plot = graf_genera_tree_gg, filename = "graf_genera_tree_gg.pdf", dpi = 300, height = 210, width = 297, units = "mm")

#### Uprava popisku ggplot ####

y_label <- element_text(face = "bold.italic", color = "black", size = 12)
x_label <- element_text(face = "bold", color = "black", size = 12)
axis_titles <- element_text(size = 14, face = "plain")
plot_title <- element_text(size = 16, face = "bold")
legend_title <- element_text(size = 14, face = "bold")
legend_text <- element_text(size = 12, face = "italic")

###RODY YEAR####

table_genus_year <- left_join(table_genus, treshold_genus, by = "genus") %>%
  replace_na(list(treshold_count = 0)) %>% # nahrazeni NA nulou
  arrange(desc(treshold_count), genus) %>%
  full_join(samples_w, by = "sample") %>%
  select(-id) %>%
  print %>%
  group_by(genus, year) %>% # podle ceho grupovat
  mutate(avg = mean(per)) %>% # prumery
  print 

stat_genus_year <- summarySE(table_genus_year, measurevar = "per", groupvars = c("genus", "year"), conf.interval = 0.95)

graf_genus_abundant_year <- filter(table_genus_year, treshold_count > 1) %>% # vybrat kolikrat mely nad treshold
  group_by(genus, year) %>%
  summarise(avg = mean(per)) %>%
  print

graf_genus_rare_year <- filter(table_genus_year, treshold_count <= 1) %>% # vybrat kolikrat mely pod treshold
  group_by(genus, year) %>%
  summarise(avg = mean(per)) %>%
  print %>%
  ungroup %>%
  group_by(year) %>%
  mutate(sum_per_year = sum(avg)) %>% # soucet v ramci dtromu pro other
  print %>%
  summarise(avg = sum(avg)) %>% # soucet v ramci dtromu pro other
  print

graf_genera_year_final <- bind_rows(graf_genus_abundant_year, graf_genus_rare_year) %>% # spojeni abundantnich a other
  replace_na(list(genus = " other")) %>% # pojmenovani other
  filter(genus == " other" | genus == "unclass") %>% # filtr other a unclass
  group_by(year) %>%
  summarise("avg" = sum(avg)) %>% # soucet pro other a unclass
  ungroup %>%
  bind_rows(filter(graf_genus_abundant_year, genus != "unclass")) %>% # spojeni souctu s tabulkou abundatnich bez unclass
  replace_na(list(genus = " other")) %>% # opet nutne pojmenovat, mezera pro lepsi trideni v ggplotu
  left_join(stat_genus_year, by = c("genus", "year")) %>%
  print
write.table(graf_genera_year_final, file = "graf_genera_year_final.csv", row.names = FALSE, col.names = TRUE, dec = ".", sep = ";")
# write.xlsx(graf_genera_year_final, "graf_genera_year_final.xlsx", sheetName="genera")

graf_genera_year_gg <- ggplot(graf_genera_year_final, aes(year)) +
  geom_bar(aes(fill = genus, weight = avg), position = "fill", width = .5) +
  # scale_y_discrete (limits = "Acidobacteria", "Bacteroidetes") +
  scale_fill_brewer(palette = "Set3", type = "div") +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  labs(x = "Year", y = "Relative abundance") 
graf_genera_year_gg
ggsave(plot = graf_genera_year_gg, filename = "graf_genera_year_gg.pdf", dpi = 300, height = 210, width = 297, units = "mm")

graf_genera_year_gg <- ggplot(graf_genera_year_final, aes(genus, fill = year)) +  # horizontalni graf po rodech
  geom_bar(aes(weight = avg), position = "dodge", width = .8) +
  scale_fill_brewer(palette = "Set1", type = "div", name = "Year") +
  # scale_y_continuous(labels = scales::percent) +
  geom_errorbar(aes(ymin = avg - se, ymax = avg + se),
                size = .7,    
                width = .3,
                position = position_dodge(width = .9)) +
  theme_bw() +
  theme(axis.text.y = y_label, axis.text.x = x_label, axis.title = axis_titles, plot.title = plot_title, legend.title = legend_title, legend.text = legend_text) + 
  xlab("Bacterial genera") + 
  ylab("Relative abundance (%)") + 
  ggtitle("Genera per year") +
  coord_flip()
graf_genera_year_gg
ggsave(plot = graf_genera_year_gg, filename = "graf_genera_year_gg.pdf", dpi = 300, height = 210, width = 297, units = "mm")

####KMENY####

proteo_class <- group_by(table_all, class, sample) %>% # procenta pro class u Proteobacteria
  summarise(seq_number = sum(count)) %>%
  ungroup %>%
  group_by(sample) %>%
  mutate(per = (seq_number / sum(seq_number) * 100)) %>%
  ungroup %>%
  filter(class == "Alphaproteobacteria" | 
           class == "Betaproteobacteria" | 
           class == "Gammaproteobacteria" | 
           class == "Deltaproteobacteria" | 
           class ==  "Epsilonproteobacteria" |
           class == "Acidithiobacillia" | 
           class == "Oligoflexia") %>% # vyfiltrovat pouze Proteobacteria 
  rename(phylum = class) %>% # pro shodu hlavicek v dalsim bloku
  arrange(desc(per), phylum) %>% # serazeni
  print
# write.table(proteo_class, file = "proteo_class.csv", row.names = FALSE, col.names = TRUE, dec = ".", sep = ";")
  
table_phylum <- group_by(table_all, phylum, sample) %>% # tabulka procenta, serazeno
  summarize(seq_number = sum(count)) %>%
  ungroup %>%
  group_by(sample) %>%
  mutate(per = (seq_number / sum(seq_number) * 100)) %>%
  ungroup %>%
  filter(phylum != "Proteobacteria") %>%
  bind_rows(proteo_class) %>% # pripojeni class od Proteobacteria
  arrange(desc(per), phylum) %>%
  print
# write.table(table_phylum, file = "table_phylum.csv", row.names = FALSE, col.names = TRUE, dec = ".", sep = ";")

treshold_phylum <- group_by(table_phylum, phylum) %>%   # pocty sekvenci ve vzorcich
  filter(per > 1) %>%                                   # TRESHOLD
  summarise(treshold_count = n()) %>%
  print

###KMENY TREE####

table_phylum_tree <- left_join(table_phylum, treshold_phylum, by = "phylum") %>%
  replace_na(list(treshold_count = 0)) %>% # nahrazeni NA nulou
  arrange(desc(treshold_count), phylum) %>%
  full_join(samples_w, by = "sample") %>%
  select(-id) %>%
  print %>%
  group_by(phylum, tree) %>% # podle ceho grupovat
  mutate(avg = mean(per)) %>% # prumery
  print 

stat_phylum_tree <- summarySE(table_phylum_tree, measurevar = "per", groupvars = c("phylum", "tree"), conf.interval = 0.95)

graf_phylum_abundant_tree <- filter(table_phylum_tree, treshold_count > 1) %>% # vybrat kolikrat mely nad treshold
  group_by(phylum, tree) %>%
  summarise(avg = mean(per)) %>%
  print

graf_phylum_rare_tree <- filter(table_phylum_tree, treshold_count <= 1) %>% # vybrat kolikrat mely pod treshold
  group_by(phylum, tree) %>%
  summarise(avg = mean(per)) %>%
  print %>%
  ungroup %>%
  group_by(tree) %>%
  mutate(sum_per_tree = sum(avg)) %>% # soucet v ramci dtromu pro other
  print %>%
  summarise(avg = sum(avg)) %>% # soucet v ramci dtromu pro other
  print

graf_phyla_tree_final <- bind_rows(graf_phylum_abundant_tree, graf_phylum_rare_tree) %>% # spojeni abundantnich a other
  replace_na(list(phylum = " other")) %>%
  left_join(stat_phylum_tree, by = c("phylum", "tree")) %>%
  print
write.table(graf_phyla_tree_final, file = "graf_phyla_tree_final.csv", row.names = FALSE, col.names = TRUE, dec = ".", sep = ";")
# write.xlsx(graf_phyla_tree_final, "graf_phyla_tree_final.xlsx", sheetName="phyla")

graf_phyla_tree_gg <- ggplot(graf_phyla_tree_final, aes(tree)) +
  geom_bar(aes(fill = phylum, weight = avg), position = "fill", width = .5) +
  # scale_y_discrete (limits = "Acidobacteria", "Bacteroidetes") +
  scale_fill_brewer(palette = "Set3", type = "div") +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  labs(x = "Tree", y = "Relative abundance") 
graf_phyla_tree_gg

graf_phyla_tree_gg <- ggplot(graf_phyla_tree_final, aes(phylum, fill = tree)) +  # horizontalni graf po rodech
  geom_bar(aes(weight = avg), position = "dodge", width = .8) +
  scale_fill_brewer(palette = "Set1", type = "div", name = "Tree species", labels = c("Quercus sp.", "Abies sp.", "Picea sp.")) +
  #  scale_y_continuous(labels = scales::percent) +
  geom_errorbar(aes(ymin = avg - se, ymax = avg + se),
                size = .7,    
                width = .3,
                position = position_dodge(width = .9)) +
  theme_bw() +
  theme(axis.text.y = y_label, axis.text.x = x_label, axis.title = axis_titles, plot.title = plot_title, legend.title = legend_title, legend.text = legend_text) + 
  xlab("Bacterial phyla") + 
  ylab("Relative abundance (%)") + 
  ggtitle("Phyla per tree") +
  coord_flip()
graf_phyla_tree_gg
ggsave(plot = graf_phyla_tree_gg, filename = "graf_phyla_tree_gg.pdf", dpi = 300, height = 210, width = 297, units = "mm")

###KMENY YEAR####

table_phylum_year <- left_join(table_phylum, treshold_phylum, by = "phylum") %>%
  replace_na(list(treshold_count = 0)) %>% # nahrazeni NA nulou
  arrange(desc(treshold_count), phylum) %>%
  full_join(samples_w, by = "sample") %>%
  select(-id) %>%
  print %>%
  group_by(phylum, year) %>% # podle ceho grupovat
  mutate(avg = mean(per)) %>% # prumery
  print 

stat_phylum_year <- summarySE(table_phylum_year, measurevar = "per", groupvars = c("phylum", "year"), conf.interval = 0.95)

graf_phylum_abundant_year <- filter(table_phylum_year, treshold_count > 1) %>% # vybrat kolikrat mely nad treshold
  group_by(phylum, year) %>%
  summarise(avg = mean(per)) %>%
  print

graf_phylum_rare_year <- filter(table_phylum_tree, treshold_count <= 1) %>% # vybrat kolikrat mely pod treshold
  group_by(phylum, year) %>%
  summarise(avg = mean(per)) %>%
  print %>%
  ungroup %>%
  group_by(year) %>%
  mutate(sum_per_year = sum(avg)) %>% # soucet v ramci dtromu pro other
  print %>%
  summarise(avg = sum(avg)) %>% # soucet v ramci stromu pro other
  print

graf_phyla_year_final <- bind_rows(graf_phylum_abundant_year, graf_phylum_rare_year) %>% # spojeni abundantnich a other
  replace_na(list(phylum = " other")) %>%
  left_join(stat_phylum_year, by = c("phylum", "year")) %>%
  print
write.table(graf_phyla_year_final, file = "graf_phyla_year_final.csv", row.names = FALSE, col.names = TRUE, dec = ".", sep = ";")
# write.xlsx(graf_phyla_year_final, "graf_phyla_year_final.xlsx", sheetName="Sheet1")

graf_phyla_year_gg <- ggplot(graf_phyla_year_final, aes(year)) +
  geom_bar(aes(fill = phylum, weight = avg), position = "fill", width = .5) +
  # scale_y_discrete (limits = "Acidobacteria", "Bacteroidetes") +
  scale_fill_brewer(palette = "Set3", type = "div") +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  labs(x = "Year", y = "Relative abundance") 
graf_phyla_year_gg

graf_phyla_year_gg <- ggplot(graf_phyla_year_final, aes(phylum, fill = year)) +  # horizontalni graf po rodech
  geom_bar(aes(weight = avg), position = "dodge", width = .8) +
  scale_fill_brewer(palette = "Set1", type = "div", name = "Year") +
  #  scale_y_continuous(labels = scales::percent) +
  geom_errorbar(aes(ymin = avg - se, ymax = avg + se),
                size = .7,    
                width = .3,
                position = position_dodge(width = .9)) +
  theme_bw() +
  theme(axis.text.y = y_label, axis.text.x = x_label, axis.title = axis_titles, plot.title = plot_title, legend.title = legend_title, legend.text = legend_text) + 
  xlab("Bacterial phyla") + 
  ylab("Relative abundance (%)") + 
  ggtitle("Phyla per year") +
  coord_flip()
graf_phyla_year_gg
ggsave(plot = graf_phyla_year_gg, filename = "graf_phyla_year_gg.pdf", dpi = 300, height = 210, width = 297, units = "mm")

###DIVERZITA YEAR####

table_div <- left_join(div_w, samples_w, by = "sample") # vytvoreni klicove tabulky pro diverzitu

stat_div_year_chao <- summarySE(table_div, measurevar = "Chao-1", groupvars = "year", conf.interval = 0.95) # statistika pro jednotliva meritka
stat_div_year_chao <- mutate(stat_div_year_chao, indice = "Chao-1") %>% 
  rename(avg_val = `Chao-1`)

stat_div_year_sw <- summarySE(table_div, measurevar = "Shannon-Wiener Diversity Index", groupvars = "year", conf.interval = 0.95)
stat_div_year_sw <- mutate(stat_div_year_sw, indice = "Shannon-Wiener Diversity Index")  %>%
  rename(avg_val = `Shannon-Wiener Diversity Index`) 

stat_div_year_even <- summarySE(table_div, measurevar = "Evenness", groupvars = "year", conf.interval = 0.95)
stat_div_year_even <- mutate(stat_div_year_even, indice = "Evennness") %>%
  rename(avg_val = `Evenness`) 

stat_div_year_80 <- summarySE(table_div, measurevar = "Species Richness - 80% diversity", groupvars = "year", conf.interval = 0.95)
stat_div_year_80 <- mutate(stat_div_year_80, indice = "80% diversity") %>%
  rename(avg_val = `Species Richness - 80% diversity`) 

indices_all_year <- bind_rows(stat_div_year_chao, stat_div_year_sw) %>% # spojeni statistiky dohromady
  bind_rows(stat_div_year_even) %>%
  bind_rows(stat_div_year_80)
write.table(indices_all_year, file = "indices_all_year.csv", row.names = FALSE, col.names = TRUE, dec = ".", sep = ";")

indices_year <- ggplot(indices_all_year, aes(year, avg_val, fill = year)) + # facet vrap graf pro indices vyse
  geom_bar(stat = "identity", width = 0.5) +
  geom_errorbar(aes(ymin = avg_val - se, ymax = avg_val + se), colour = "black", width = .2) +
  facet_wrap(~ indice, scales = "free") +
  scale_fill_brewer(palette = "Set1", name = "Year") +
  xlab("Year") + 
  ylab("Diversity values") +
  theme(legend.position = "none", axis.text.y = y_label, axis.text.x = x_label, axis.title = axis_titles, plot.title = plot_title, strip.text = element_text(size = 16))
indices_year
ggsave(plot = indices_year, filename = "indices_year.pdf", dpi = 300, height = 210, width = 297, units = "mm")

ggplot(stat_div_year_chao, aes(year, `Chao-1`, colours = )) + # funguje pro jednotlive indices
  geom_bar(stat = "identity", width = 0.5) +
  geom_errorbar(aes(ymin = `Chao-1` - se, ymax = `Chao-1` + se), colour = "black", width = .2)

###DIVERZITA TREE####

stat_div_tree_chao <- summarySE(table_div, measurevar = "Chao-1", groupvars = "tree", conf.interval = 0.95) # statistika pro jednotliva meritka
stat_div_tree_chao <- mutate(stat_div_tree_chao, indice = "Chao-1") %>% 
  rename(avg_val = `Chao-1`)

stat_div_tree_sw <- summarySE(table_div, measurevar = "Shannon-Wiener Diversity Index", groupvars = "tree", conf.interval = 0.95)
stat_div_tree_sw <- mutate(stat_div_tree_sw, indice = "Shannon-Wiener Diversity Index")  %>%
  rename(avg_val = `Shannon-Wiener Diversity Index`) 

stat_div_tree_even <- summarySE(table_div, measurevar = "Evenness", groupvars = "tree", conf.interval = 0.95)
stat_div_tree_even <- mutate(stat_div_tree_even, indice = "Evennness") %>%
  rename(avg_val = `Evenness`) 

stat_div_tree_80 <- summarySE(table_div, measurevar = "Species Richness - 80% diversity", groupvars = "tree", conf.interval = 0.95)
stat_div_tree_80 <- mutate(stat_div_tree_80, indice = "80% diversity") %>%
  rename(avg_val = `Species Richness - 80% diversity`) 

indices_all_tree <- bind_rows(stat_div_tree_chao, stat_div_tree_sw) %>% # spojeni statistiky dohromady
  bind_rows(stat_div_tree_even) %>%
  bind_rows(stat_div_tree_80)
write.table(indices_all_tree, file = "indices_all_tree.csv", row.names = FALSE, col.names = TRUE, dec = ".", sep = ";")

indices_tree <- ggplot(indices_all_tree, aes(tree, avg_val, fill = tree)) + # facet vrap graf pro indices vyse
  geom_bar(stat = "identity", width = 0.5) +
  geom_errorbar(aes(ymin = avg_val - se, ymax = avg_val + se), colour = "black", width = .2) +
  facet_wrap(~ indice, scales = "free") +
  scale_fill_brewer(palette = "Set1", name = "Tree") +
  xlab("Tree") + 
  ylab("Diversity values") +
  theme(legend.position = "none", axis.text.y = y_label, axis.text.x = x_label, axis.title = axis_titles, plot.title = plot_title, strip.text = element_text(size = 16))
indices_tree
ggsave(plot = indices_tree, filename = "indices_tree.pdf", dpi = 300, height = 210, width = 297, units = "mm")

ggplot(stat_div_tree_chao, aes(tree, `Chao-1`, colours = )) + # funguje pro jednotlive indices
  geom_bar(stat = "identity", width = 0.5) +
  geom_errorbar(aes(ymin = `Chao-1` - se, ymax = `Chao-1` + se), colour = "black", width = .2)
