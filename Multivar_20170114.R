
library(vegan)
library(vegan3d)
library(BiodiversityR) # for netsed.npmanova
library(RColorBrewer)
library(ggplot2)
library(readr)
library(xlsx)
library(Rmisc)
library(dplyr)
library(tidyr)
library(agricolae) # anova a post hoc test
library(gridExtra)
library(grid)
library(spdep)
library(fields)

#### Uprava tabulek ####

setwd("D:/Vojta/Disk/Laborka/Data/Projekty/Pralesy_Zofin/Zofin_2013_working/New_database_analysis_R_script")

otus <- read_delim("D:/Vojta/Disk/Laborka/Data/Projekty/Pralesy_Zofin/Zofin_2013_working/New_database_analysis_R_script/otus.csv", 
                   ";", escape_double = FALSE, trim_ws = TRUE)
View(otus)
colnames(otus)[which(names(otus) == "Group name:")] <- "cluster_name"

env <- read.table(file = "D:/Vojta/Disk/Laborka/Data/Projekty/Pralesy_Zofin/Zofin_2013_working/New_database_analysis_R_script/env.csv", 
                  header = TRUE, sep = ";", dec = ".")
View(env)
dim(env)

colnames(env)[colnames(env) == "Ergosterol..ppm."] <- "Ergosterol"
env$Year <- as.factor(env$Year) 
env$Tree <- as.factor(env$Tree) 
str(env)

otus_w <- tbl_df(otus)
env_w <- tbl_df(env)

env_NA <- group_by(env_w, Tree, Year) %>% # tato fce nahradi NAs prumery pro shluky vzorku Tree+Year
  mutate(lig = ifelse(is.na(lig), mean(lig, na.rm = TRUE), lig)) %>%
  mutate(erg = ifelse(is.na(erg), mean(erg, na.rm = TRUE), erg)) %>%
  mutate(water = ifelse(is.na(water), mean(water, na.rm = TRUE), water)) %>%
  mutate(fb_ratio = ifelse(is.na(fb_ratio), mean(fb_ratio, na.rm = TRUE), fb_ratio)) %>%
  mutate(bac_copy = ifelse(is.na(bac_copy), mean(bac_copy, na.rm = TRUE), bac_copy)) %>%
  mutate(fun_copy = ifelse(is.na(fun_copy), mean(fun_copy, na.rm = TRUE), fun_copy)) %>%
  ungroup() %>%
  select(-dens)


attach(env_NA)
detach(env_NA)


otus_percent <- gather(otus_w, sample, count, contains("VOJ")) %>%  # procenta pro OTUs
  group_by(sample) %>%
  mutate(seq_number = sum(count)) %>%
  ungroup %>%
  mutate(per = (count / seq_number * 100))

treshold_otus <- group_by(otus_percent, cluster_name) %>%  
  filter(per > 0.5) %>%                                   # TRESHOLD pro procenta
  summarise(treshold_count = n()) %>%
  print

otus_multivar_titles <- left_join(otus_percent, treshold_otus, by = "cluster_name") %>%
  replace_na(list(treshold_count = 0)) %>% # nahrazeni NA nulou
  filter(treshold_count > 10) %>%  # kolikrat je OTU nad treshold
  distinct(cluster_name) %>%
  print
# write.table(otus_multivar, file = "otus_multivar.csv", row.names = FALSE, col.names = TRUE, dec = ".", sep = ";")

otus_multivar <- select(otus_percent, c(cluster_name, sample, per)) %>%
  spread(sample, per) %>%
  select(-VOJ112) %>% #### pro tento vzorek nejsou env data
  right_join(otus_multivar_titles, by = "cluster_name") # filtrovani pouze OTUs nad treshold
otus_multivar <- as.data.frame(otus_multivar) # zpet na data frame

cluster <- otus_multivar[,1]
rownames(otus_multivar) = cluster # pojmenovani radku podle prvniho sloupce
otus_multivar <- t(otus_multivar[ , 2:ncol(otus_multivar)]) # transpozice a vypusteni prvniho sloupce

spechell <- decostand(otus_multivar, "hellinger") #hellinger transformation, multivar OTUs
spechell_all <- decostand(otus_permut, "hellinger") #hellinger transformation, all OTUs

#### Abundant OTUs ####

otus_abund <- select(otus_percent, c(-seq_number, -count, -No.)) %>% 
  spread(sample, per) %>%
  left_join(tax_w, by = "cluster_name") %>%  # spojeni procent s tax
  gather(sample, per, contains("VOJ"), na.rm = FALSE) %>%  # preformatovani vzorku
  full_join(samples_w, by = "sample") %>%  # spojeni s charakteristikou vzorku
  select(cluster_name, Accession, Similarity, Coverage, phylum, class, genus, sample, per, tree, year) %>%
  replace_na(list(phylum = "unclass", class = "unclass", genus = "unclass")) %>% # nahrazeni NAs
  semi_join(otus_multivar_titles, by = "cluster_name") %>%  # filtr abundant OTUs z multivar
  summarySE(measurevar = "per", groupvars = c("year", "tree", "cluster_name"), conf.interval = 0.95)

otus_abund_avg <- select(otus_abund, year, tree, cluster_name, per) %>% # tabulka prumernych hodnot pro skupiny tree_year
  unite(sample_group, year, tree) %>% # spojeni skupin
  spread(sample_group, per)
write.table(otus_abund_avg, file = "otus_abund_avg.csv", row.names = FALSE, col.names = TRUE, dec = ".", sep = ";")

otus_abund_sd <- select(otus_abund, year, tree, cluster_name, se) %>% # tabulka SE hodnot pro skupiny tree_year
  unite(sample_group, year, tree) %>%
  spread(sample_group, se)
write.table(otus_abund_sd, file = "otus_abund_sd.csv", row.names = FALSE, col.names = TRUE, dec = ".", sep = ";")


#### NMDS ####

#vychazime z matice vzdalenosti
vegdist(otus_multivar, method = "bray")   # bez transformace
vegdist(decostand(otus_multivar, "hellinger"), method = "bray")   #hellinger transf, bray curtis diss
specdist = vegdist(otus_multivar, method = "bray", diag = T, upper = T)   #zobrazeni cele matice vcetne diagonaly, kde jsou nuly 
image(as.matrix(specdist)) #zobrazeni

set.seed(321)
mdsord = metaMDS(comm = otus_multivar, distance = "bray", trace = FALSE, k = 2, trymax = 200) #bez transformace
mdsord = metaMDS(comm = spechell, distance = "bray", trace = FALSE) #transf, nejmensi stres znamena nejvetsi R2
# matici vzdalenosti si spocte samo, bray je default, neni treba pouzivat vegdist

mdsord 
plot(mdsord, type = "t")
plot(mdsord, display = "sites", type = "t") # pouze druhy
stressplot(mdsord)

# Primy pristup ke skorum

scores(mdsord) 
scores(mdsord, display = "sites") 

mdsord$points
mdsord$species

#### ZOBRAZENI #####

##START fungujiciho skriptu##
pdf("NMDS_NAG_izo_transformed.pdf", paper = 'a4r')
plot(mdsord, disp = "sites", type = "n")
ordisurf(mdsord, pH, main = "", col = "forestgreen", display = "sites", type="n") #prida izocary pro promenne
#ordiellipse(mdsord, main="",col="black", cex=0.05, draw = "lines", Year)
pchlist = c(3:13)
points(mdsord, disp = "sites", pch = 16, col = heat.colors(3))
text(mdsord, labels=env$Year, display="sites",cex=0.4, pos=1)
graf = envfit(mdsord, cbind(CBH, bG, EC, F.B, NAG, PLFAB)) #funguje lepe nez radky vyse
plot(graf)
dev.off()

pdf("NMDS_final.pdf", paper = 'a4r')
tiff("NMDS_tree_year2.tiff", height = 210, width = 297, units = "mm", compression = "lzw", res = 300)
# par(mfrow = c(1, 2)) #dva grafy vedle sebe
ordisurf(mdsord, water, main = "", col = "forestgreen", display = "sites", type = "n") #prida izocary pro promenne
# pchlist = c("#E495A5", "#D89F7F", "#BDAB66" ,"#96B56C" ,"#65BC8C", "#39BEB1", "#55B8D0", "#91ACE1", "#C29DDE" ,"#DE94C8")
plot(mdsord, disp = "sites", type = "n") #, xlim = c(-1.5,1.5), ylim = c(-1,1))
# mycolors <- (terrain.colors(4))
mycolors3 <- brewer.pal(3, "Set2")
mycolors4 <- brewer.pal(4, "Set1")
pchlist = c(15, 16, 17)
points(mdsord, disp = "sites", pch = pchlist[env$Tree], cex = 1.4, col = mycolors4[env$Year])
points(mdsord, disp = "sites", cex = 1.4, col = mycolors3[env$Tree], pch = 16)
for (i in 1:nlevels(env$Year))
  ordihull(mdsord, env$Year, show.groups=levels(env$Year)[i],col=brewer.pal(4, "Spectral")[i])
for (i in 1:nlevels(env$Tree))
  ordihull(mdsord, env$Tree, show.groups=levels(env$Tree)[i], label = FALSE, col = brewer.pal(3, "Set2")[i])
for (i in 1:nlevels(env$Year))
  ordiellipse(mdsord, env$Year, show.groups=levels(env$Year)[i], kind = c("sd"), conf=0.95, label = FALSE, col = brewer.pal(4, "Spectral")[i])
for (i in 1:nlevels(env$Tree))
  ordiellipse(mdsord, env$Tree, show.groups=levels(env$Tree)[i], kind = c("sd"), conf=0.95, label = FALSE, col = brewer.pal(3, "Spectral")[i])
text(mdsord,labels = env$Tree, display = "sites",cex = 0.8, pos = 2)
text(mdsord,labels = env$ID, display = "sites",cex = 0.8, pos = 1)
legend("bottomright", pch = pchlist, legend = c("Abies","Fagus","Picea"))
legend("bottomleft",col = mycolors4, pch = 19, legend = c("1975","1997","2008","2013"))
legend("bottomleft",col = mycolors3, pch = 19, legend=c("Fagus","Abies","Picea"))
ordihull(mdsord, grp, lty = 2, col = "darkgreen") #je treba spocitat clustery v kodu nize
set.seed(321)
Year <- envfit(mdsord, env_NA$Year, permutations = 9999)
graf <- envfit(mdsord, permutations = 9999, cbind(pH, N, C, Tree, fb_ratio, lig, erg, water)) #funguje lepe nez radky vyse
plot(graf)
plot(Year)
dev.off()

tiff("NMDS_species.tiff", height = 210, width = 297, units = "mm", compression = "lzw", res = 300)
plot(mdsord, disp="species", type="n")
points(mdsord, display="sp", pch=17, cex=1)
text(mdsord, display="sp", cex=0.8, pos=1, select = mycol, col="brown") #mykolyticke
text(mdsord, display="sp", cex=0.8, pos=1, select = cel, col="darkgreen") #celulolyticke
ordihull(mdsord, grp, lty = 2, col = "darkgreen")
Year = envfit(mdsord, env$Year)
graf = envfit(mdsord, cbind(CEL, bG, NAG, pH, CBH, Tree)) #funguje lepe nez radky vyse
plot(graf)
plot(Year)
dev.off()

orditorp(mdsord,display="sites",cex=0.5,air=0.01)  #nahradi body textem
ordicluster(mdsord,main="",col="red", lwd=1, hclust(vegdist(spechell))) #prida clustery
ordihull(mdsord, Tree)
for (i in 1:nlevels(env$Year))
  ordihull(mdsord, env$Year, show.groups=levels(env$Year)[i],col=terrain.colors(length(table(env$Year)))[i]) #obarvi obalky dle zadani

#### ggplot2 ####

#data for plotting NMDS points
zof.NMDS.data <- select(env_NA, ID, Sample, Tree, Year) # there are other ways of doing this. But this is the way I do it for ease of plotting
zof.NMDS.data$NMDS1<-mdsord$points[ ,1] # this puts the NMDS scores for the plots into a new dataframe. you could put them into an existing one if you preferred.
zof.NMDS.data$NMDS2<-mdsord$points[ ,2]

# data for the envfit arrows
env.scores.zof <- as.data.frame(scores(graf, display = "vectors")) # extracts relevant scores from envfit
env.scores.zof <- rbind(env.scores.zof, scores(Year, display = "vectors"))
env.scores.zof <- cbind.data.frame(env.scores.zof, env.variables = rownames(env.scores.zof), stringsAsFactors = FALSE) #and then gives them their names
class(env.scores.zof[7,3])
env.scores.zof[4,3] <- "tree"

# Hulls - grouping
grp.1975 <- zof.NMDS.data[zof.NMDS.data$Year == "1975", ][chull(zof.NMDS.data[zof.NMDS.data$Year == 
                                                                                "1975", c("NMDS1", "NMDS2")]), ]  # hull values for grp A
grp.1997 <- zof.NMDS.data[zof.NMDS.data$Year == "1997", ][chull(zof.NMDS.data[zof.NMDS.data$Year == 
                                                                                "1997", c("NMDS1", "NMDS2")]), ] # hull values for grp B
grp.2008 <- zof.NMDS.data[zof.NMDS.data$Year == "2008", ][chull(zof.NMDS.data[zof.NMDS.data$Year == 
                                                                                "2008", c("NMDS1", "NMDS2")]), ]
grp.2013 <- zof.NMDS.data[zof.NMDS.data$Year == "2013", ][chull(zof.NMDS.data[zof.NMDS.data$Year == 
                                                                                "2013", c("NMDS1", "NMDS2")]), ]
hull.data <- rbind(grp.1975, grp.1997, grp.2008, grp.2013)  #combine grp.a and grp.b
hull.data

# finally plotting
(zof.nmds.gg1 <- ggplot(data = zof.NMDS.data, aes(y = NMDS2, x = NMDS1)) + 
    geom_point(aes(shape = Tree, colour = Year), show.legend = TRUE, size = 3) +
    # stat_ellipse(aes(color = Year), type = "norm", size = 1, level = 0.95, geom = "path") +
    geom_segment(data = env.scores.zof,
                 aes(x = 0, xend = 2.5*NMDS1, y = 0, yend = 2.5*NMDS2),
                 arrow = arrow(length = unit(0.25, "cm")), colour = "#d68f0c") +
    geom_text(data = env.scores.zof,
              aes(x = 2.5*NMDS1, y = 2.5*NMDS2, label = env.variables),
              size = 5,
              hjust = -0.3) +
    geom_polygon(data = hull.data, aes(x = NMDS1, y = NMDS2, fill = Year, group = Year), alpha = 0.05, show.legend = FALSE) + 
    guides(fill = "none") +
    scale_shape_discrete(name = "Tree species", labels = c("F. sylvatica", "A. alba", "P. abies")) +
    scale_colour_manual(values = mycolors4 , name = "Age class", labels = c(">38", "38-16", "15-5", "<5")) +
    scale_x_continuous(breaks = seq(-1, 1.5, by = 1)) +
    scale_y_continuous(breaks = seq(-1, 1.5, by = 1)) +
    theme_bw() +
    theme(axis.title.x = element_text(face = "plain", colour = "black", size = 16), axis.text.x  = element_text(angle = 0, vjust = 0.5, size = 16)) +
    theme(axis.title.y = element_text(face="plain", colour="black", size = 16), axis.text.y  = element_text(angle=0, vjust=0.5, size=16)) +
    theme(legend.text = element_text(colour="black", size = 16, face = "italic")) +
    theme(legend.title = element_text(colour="black", size = 16, face="bold")) +
    ggtitle("Figure 2"))

ggsave(plot = zof.nmds.gg1, filename = "zof.nmds.gg1.pdf", dpi = 300, height = 210, width = 297, units = "mm")

#### PcoA ####

pcoa <- cmdscale(s.dis, list. = TRUE, k = 3)
ordiplot(pcoa)
ordisurf(pcoa, env$N, add = TRUE)

##### DCA ####

ccout1 = decorana(spechell)
ccout1 = decorana(otus_multivar)
ccout1
#podil os na celkove variabilite

cceig = ccout1$evals / cca(spec1)$tot.chi #celkova variabilita je treba ziskat pomoci CA
barplot(rbind(cceig,cumsum(cceig)),beside=T)
legend("topleft",fill=c("grey25","grey75"),legend=c("each axis","cumulative"))


plot(ccout1)
plot(ccout1,display="sp",type="n") #vzorky, clustery
text(ccout1,display="sp",cex=0.4)

pdf("DCA_transformed_samples.pdf", paper='a4r')
tiff("DCA_tree_year.tiff", height = 210, width = 297, units = "mm", compression = "lzw", res = 300)
ordiplot (ccout1, display = 'sites', type = 'n')  #BEST, obarvi dle mesice
mycolors3 <- brewer.pal(3, "Set2")
mycolors4 <- brewer.pal(4, "Spectral")
pchlist=c(15, 16, 17)
points(ccout1, disp="sites",pch=pchlist[env$Tree], cex=1.4, col=mycolors4[env$Year])
points(ccout1, disp="sites", cex=1.4, col=mycolors3[env$Tree], pch=16)
for (i in 1:nlevels(env$Year))
  ordihull(ccout1, env$Year, show.groups=levels(env$Year)[i],col=brewer.pal(4, "Spectral")[i])
for (i in 1:nlevels(env$Tree))
  ordihull(ccout1, env$Tree, show.groups=levels(env$Tree)[i], label = FALSE, col = brewer.pal(3, "Set2")[i])
text(ccout1,labels=env$Tree, display="sites",cex=0.8, pos=1)
text(ccout1,labels=env$Year, display="sites",cex=0.8, pos=2)
legend("bottomright", pch=pchlist, legend=c("Quercus","Abies","Picea"))
legend("bottomleft",col=mycolors4, pch=19, legend=c("1975","1997","2008","2013"))
legend("bottomleft",col=mycolors3, pch=19, legend=c("Quercus","Abies","Picea"))
ordihull(ccout1, grp, lty = 2, col = "darkgreen") #je treba spocitat clustery v kodu nize


points (ccout1, col = env$Year, pch=16, cex=1)
text(ccout1,labels=env$Year, display="sites", cex=0.5, pos=1)
legend("bottomleft",pch=16,cex=0.8, pt.cex=c(0.035,1.15,1.58)/1.3,title="bact biomass (ug g-1)", legend=c("1","33","45"))

Year = envfit(mdsord, env$Year)
graf = envfit(mdsord, cbind(CEL, bG, NAG, pH, CBH, Tree)) #funguje lepe nez radky vyse
plot(graf)
plot(Year)
dev.off()

pdf("DCA_transformed_clusters.pdf", paper='a4r')
tiff("DCA_clusters.tiff", height = 210, width = 297, units = "mm", compression = "lzw", res = 300)
ordiplot (ccout1, display = 'species', type = 'n')
points(ccout1, display="sp", pch=17, cex=1)
text(ccout1, display="sp", cex=0.4, pos=1, select = mycol, col="brown") #mykolyticke
text(ccout1, display="sp", cex=0.4, pos=1, select = cel, col="darkgreen") #celulolyticke

Year = envfit(mdsord, env$Year)
graf = envfit(mdsord, cbind(CEL, bG, NAG, pH, CBH, Tree)) #funguje lepe nez radky vyse
plot(graf)
plot(Year)
dev.off()

ordirgl(ccout1) #3D DCA

ordisurf(ccout1,NAG,main="",col="forestgreen", display = "sites", type="n") #izocary
ordisurf(ccout1,F.B,main="",col="forestgreen", display = "sites", type="n") #izocary

ordiplot (ccout1, display = 'sp', type = 'n')
orditorp (ccout1, display = 'sp') #pouze clustery


ordihull(ccout1, env$Year) #obalky
for (i in 1:nlevels(env$Year))
  ordihull(ccout1 ,env$Year,show.groups=levels(env$Year)[i],col=rainbow(length(table(env$Year)))[i])

#####KLASIFIKACE####

dis <- vegdist(spec1, method="bray") #hellinger transf, bray curtis diss
clus <- hclust(dis, "single") #aglomerativni clusterovani
plot(clus) # hrozi zde zretezeni, vzdy je nejaky bod blizko
plot(clus,hang=-1) #vsechny dosahuji na nulu

#vysledek dramaticky zavisi nejen na volbe vzdalenosti, ale i na shlukovacim kriteriu 
clus <- hclust(dis, "complete") #bez zretezeni, abych se spojil, musim byt dostatecne blizko k nejvzdalenejsimu
plot(clus)

clus <- hclust(dis, "average")
plot(clus)

#prehled shlukovacich kriterii na 
help(hclust) #polozka method

#robustni a doporucitelna Wardova metoda 
pdf("Clustered_ward_transformed.pdf", paper='a4r')
tiff("Clustered_ward_transformed.tiff", height = 210, width = 297, units = "mm", compression = "lzw", res = 300)
clus <- hclust(dis, "ward.D") #asi nejpouzitelnejsi
plot(clus)

#dalsi prace s dendrogramem
#plot(clus)
rect.hclust(clus, 3) #pocet klastru k vyznaceni v dendrogramu
dev.off()

cor(cophenetic(hclust(dis)), dis) #korelace dendrogramu s dissimilarity matrix
heatmap(as.matrix(dis))

#ulozit klasifikacni vektor
grp <- cutree(clus, 3) # uteti dendrogramu na urovni nami zvolene
grp

#lze s nim dal pracovat jako s libovolnou jinou promennou
pdf("bG_classification.pdf", paper='a4r', family="sans")
tiff("bG_classification_text.tiff", height = 210, width = 297, units = "mm", compression = "lzw", res = 300)
boxplot(bG ~ grp, notch=F, main="bG in classification", xlab="clusters", ylab="bG values (standardized)") #nasazeni env dat na compositional data, nutna spojita data
dev.off()

table(grp, bG)
#lze barplot

#propojit klasifikaci s ordinaci, pozor, castecne dukaz kruhem

mdsord = metaMDS(comm = specdata, distance = "bray", trace = FALSE) 
plot(mdsord, display="sites")
ordicluster(mdsord, clus, col="salmon") #bez useknuti vyssi urovne hierarchie
ordihull(mdsord, grp, lty = 2, col = "darkgreen")

plot(ccout1, display="sites")
ordicluster(ccout1, clus, col="blue", prune=4) #useknout vyssi urovne hierarchie


#minimum spanning tree, spojeni objektu s nejmensi vzdalenosti

mst <- spantree(dis, toolong = 1)
plot(mst, ord=mdsord)


####K-means clustering####
# lisi se od Wardovy metody, rovnou deli objekt na 5 skupin

xx=kmeans(spec1, 4) #zadam kolik clusteru chci mit
xx

xx$centers
xx$cluster

#opet lze s vysledky klasifikace dal pracovat
grp = xx$clusters
boxplot(bG ~ grp, notch=T) 


####PCA#######

#data jsou pokryvnosti druhu, cili ve srovnatelnych jednotkach, neni nezbytne standardizovat

rdout=rda(spec1) #vlastni vypocet
rdout=rda(spechell) #vlastni vypocet


summary(rdout)
eigenvals(rdout) / sum(eigenvals(rdout))
barplot(as.numeric(eigenvals(rdout) / sum(eigenvals(rdout))))
barplot(cumsum(eigenvals(rdout) / sum(eigenvals(rdout))))
barplot(cumsum(eigenvals(rdout)[1:10] / sum(eigenvals(rdout))))

biplot(rdout)
text(rdout,display="sp")

#####INDIC SPEC####

library(labdsv)

iva <- indval(spec1, grp)  #indicator species pro clustery, library labdsv
group = iva$maxcls[iva$pval<=0.05]
indval=iva$indcls[iva$pval<=0.05]
pvalue=iva$pval[iva$pval<=0.05]
freq=apply(spec1>0,2,sum)[iva$pval<=0.05]
indicator.species=data.frame(group=group,indval= indval,pvalue= pvalue,freq= freq)
indicator.species
write.table(indicator.species,"indicator.species.txt", sep="\t")

library(indicspecies)

otus_multivar_indic <- as.data.frame(otus_permut)
# otus_multivar_indic <- otus_multivar_indic[,c(ncol(otus_multivar_indic),1:(ncol(otus_multivar_indic)-1))]
# otus_multivar_indic <- sapply(otus_multivar_indic[,1], as.factor) # nefunguje, je nutne do excel a potom importovat

#run the Indicator Species Analysis using the function "multipatt". The parameter "IndVal.g" calls the indicator species function using a verison that corrects for unequal group size, whereas "duleg" stes whether to use combintaions of a priori classes or not
set.seed(321)
indval <- multipatt(otus_multivar_indic, env$Year, func = "IndVal.g", duleg = TRUE, control = how(nperm = 9999))
#output - "A" is the specifity - "B" is the fidelity
indval
#short summary only showing significant species
summary(indval, indvalcomp = TRUE)
coverage(otus_multivar_indic, indval)
#longer summary showing all species
summary(indval, alpha = 1)
class(indval)

indic_sig <- read.csv(header = FALSE, "D:/Vojta/Disk/Laborka/Data/Projekty/Pralesy_Zofin/Zofin_2013_working/New_database_analysis_R_script/indic_sig.csv")
indicatory_tax <- dplyr::rename(indic_sig, `cluster_name` = `V1`) %>%
  left_join(tax_w, by = "cluster_name")
write.table(indicatory_tax, file = "indicatory_tax.csv", row.names = FALSE, col.names = TRUE, dec = ".", sep = ";")

#https://gist.github.com/eclarke/d9e51f7b6c797cd84256d45504542a95
# random forest model

##### PERMANOVA ####

otus_permut <- select(otus_percent, c(cluster_name, sample, per)) %>%
  spread(sample, per) %>%
  select(-VOJ112) #### pro tento vzorek nejsou env data
otus_permut <- as.data.frame(otus_permut)
# cluster <- otus_permut[,1]
rownames(otus_permut) = cluster # pojmenovani radku podle prvniho sloupce
otus_permut <- t(otus_permut[ , 2:ncol(otus_permut)]) # transpozice a vypusteni prvniho sloupce

s.dis <- vegdist(decostand(otus_permut, "hell"), "bray") # all OTUs, hellinger for degradation of extreme values, Bray Curtis dissimilar matrix
s.dis <- vegdist(spechell, "bray") # multivar OTUs

set.seed(321)
adonis(s.dis ~ env[,1], permutations = 9999) # PERMANOVA for testing influence of 1rd column on the community

set.seed(321)
adonis(s.dis ~ env$Year, permutations = 9999) # strata default - NULL
adonis(s.dis ~ env$Tree, permutations = 9999)
adonis(s.dis ~ env$space, permutations = 9999)
adonis(s.dis ~ env$pH, permutations = 9999)
adonis(s.dis ~ env$N, permutations = 9999)
adonis(s.dis ~ env$C, permutations = 9999)
adonis(s.dis ~ env_NA$lig, permutations = 9999)
adonis(s.dis ~ env_NA$water, permutations = 9999)
adonis(s.dis ~ env_NA$erg, permutations = 9999)

set.seed(321)
adonis(s.dis ~ env$Year, strata = env$Tree, permutations = 999) #strata - permutation only within Shrub categories, correct approach, but rises p value - less significant, p-value corrected
adonis(s.dis ~ env$Tree, strata = env$Year, permutations = 999) 

#nested.npmanova(response ~ groups+plots, data = data, method = "bray", permutations = 999) #corrected F-ratio

set.seed(321)
adonis(s.dis ~ env$Year * Tree, permutations = 9999, strata = env$Year)
adonis(s.dis ~ env$Year + Tree, permutations = 9999, strata = env$Year)
adonis(s.dis ~ env_NA$Tree + env_NA$Year * env_NA$pH , permutations = 9999, strata = NULL) 
adonis(s.dis ~ env_NA$Year + env_NA$pH + env_NA$lig + env_NA$water + env_NA$Tree + env_NA$N + env_NA$C, permutations = 9999, strata = NULL)
adonis(s.dis ~ env_NA$Year * env_NA$pH + env_NA$water + env_NA$lig + env_NA$Tree + env_NA$N + env_NA$C + env_NA$fb_ratio + env_NA$erg, permutations = 9999, strata = NULL) # asi nejlepsi model, odkryva factory postupne a jejich vliv odcita, neni dobre pouzivat factor z modelu jako stratu a obracene (Jari)
adonis(s.dis ~ env$Tree / Year, permutations = 999, strata = env$Year) #second nested with first
adonis(s.dis ~ env$NAG * PLFAT, permutations = 999) #test for interaction 
adonis(s.dis ~ env$NAG + PLFAT, permutations = 999) #test for additive effect
adonis(s.dis ~ env$Year / pH, permutations = 999, strata = env$Year)

RsquareAdj(0.01404, 120, 1)

##### Mantel test/ geografie #####

set.seed(321)
mantel(s.dis, vegdist(env$bG, method = "euc"), method = "pearson", permutations = 9999, strata = NULL, na.rm = FALSE)

coord <- cbind(env$x_pata, env$y_pata) #koordinaty
colnames(coord) <- c("y_pata", "x_pata")
rownames(coord) <- env$Id
coord_dist <- vegdist(coord, method = "euclidean", na.rm = TRUE) # distance pro koordinaty
mantel(s.dis, coord_dist, method = "pearson", permutations = 9999, na.rm = FALSE) # mantel pro koordinaty a OTUs
ade4::mantel.randtest(m1 = s.dis, m2 = coord_dist, nrepet = 999) # Mantel z jineho balicku

zof.resid <- resid(lm(as.matrix(s.dis) ~ ., data = as.data.frame(coord))) # detrending
zof.dis <- dist(zof.resid) # Compute the detrended species distance matrix

zof_mantel_cor <- mantel.correlog(D.eco = s.dis, D.geo = NULL, XY = coord, n.class = 0, break.pts = NULL, cutoff = FALSE, r.type = "pearson", nperm = 1000, mult = "holm", progressive = TRUE)
print(zof_mantel_cor)
plot(zof_mantel_cor)

# Moranovo I
zof.connectivity <- adegenet::chooseCN(xy = coord, type = 5, d1 = 0, d2 = "dmin", plot.nb = TRUE, result.type = "listw", edit.nb = FALSE)
zof.connectivity # konektivita

set.seed(321)
moran.test(x = mdsord$points[,1], listw = zof.connectivity, alternative = "greater", randomisation = TRUE) # testovani prvni osy z NMDS
zof.pcoa1.mctest <- moran.mc(x = mdsord$points[,1], listw = zof.connectivity, alternative = "greater", nsim = 1000)
zof.pcoa1.mctest

moran.test(x = mdsord$points[,2], listw = zof.connectivity, alternative = "greater", randomisation = TRUE) # testovani druhe osy z NMDS
zof.pcoa2.mctest <- moran.mc(x = mdsord$points[,2], listw = zof.connectivity, alternative = "greater", nsim = 1000)
zof.pcoa2.mctest
plot(zof.pcoa2.mctest) # Plot of density of permutations
moran.plot(x = mdsord$points[,2], listw = zof.connectivity)

# PCNM vektory z koordinat
distmat <- rdist(cbind(env$y_pata, env$x_pata))    # dist matrix of GPS coordinates
spatvectors <- pcnm(distmat) # PCNM vektory
spatvectors
ordiR2step(rdout, spatvectors, direction = "forward", permutations = how(nperm = 999), trace = TRUE, R2permutations = 1000) # nefunguje
sigvect <- forward.sel(spechell_all, spatvectors$vectors, nperm = 9999, alpha = 0.05) # forward selection signifikatnich vektoru, nefunguje
print(sigvect)         		# Shows significant vectors

rdout = rda(spechell)
plot(rdout)

space <- read.table("clipboard",header=T)   #Significant spatial vectors
space <- as.matrix(space)


#ANOSIM test
set.seed(321)
anos <- anosim(s.dis, env$Year, permutations = 999, distance = "bray", strata = NULL, parallel = getOption("mc.cores"))
summary(anos)
plot(anos)
#je lepsi volit adonis2 (podle helpu)
set.seed(321)
adonis2(s.dis ~ env$Year, permutations = 999, method = "bray",
        sqrt.dist = FALSE, add = FALSE, by = "terms",
        parallel = getOption("mc.cores"))

##### ANOVA + post hoc ####

#One-way ANOVA
aov1 <- filter(table_phylum_year, phylum == "Verrucomicrobia") # zadat ktere taxa
aov1 <- filter(table_genus_year, genus == "Steroidobacter") # zadat ktere taxa
aov1 <- table_div # zadat ktere indices
aov1 <- aov(aov1$per ~ aov1$year) # procenta vs roky/stromy
aov1 <- aov(aov1$`Chao-1` ~ aov1$year) # procenta vs roky/stromy
aov1 <- aov(env$fb_ratio ~ env$Year)
summary(aov1)

#TukeysHSD post hoc test
outhsd <- HSD.test(aov1, "Month", group = TRUE, console = TRUE, main = "nazev")
outhsd
#graf
bar.group(outhsd$groups, ylim = c(0,10), density = 10, border = "blue") #here we go!!!
means <- outhsd$means

#Fischer LSD post hoc test
outlsd <- LSD.test(aov1, "env$Tree", p.adj = "bonferroni")
outlsd
bar.group(outlsd$groups, ylim = c(0,10), density = 10, border = "blue") #here we go!!!

#### Uprava popisku ggplot ####

y_label <- element_text(face = "bold.italic", color = "black", size = 10)
x_label <- element_text(face = "bold", color = "black", size = 9)
axis_titles <- element_text(size = 12, face = "plain")
plot_title <- element_text(size = 14, face = "bold")
legend_title <- element_text(size = 12, face = "bold")
legend_text <- element_text(size = 11, face = "italic")

title <- textGrob("Supplemetary Figure 1", gp = gpar(fontface = "bold", cex = 1.5))

##### barplot pH, C, N, qPCR ####

stat_ph_year <- summarySE(env, measurevar = "pH", groupvars = "Year", conf.interval = 0.95)
stat_ph_year <- mutate(stat_ph_year, indice = "pH")  %>%
  dplyr::rename(avg_val = `pH`) 

stat_c_year <- summarySE(env, measurevar = "C", groupvars = "Year", conf.interval = 0.95)
stat_c_year <- mutate(stat_c_year, indice = "C content")  %>%
  dplyr::rename(avg_val = `C`)

stat_n_year <- dplyr::rename(env, n_cont = `N`) %>%
  summarySE(measurevar = "n_cont", groupvars = "Year", conf.interval = 0.95)
stat_n_year <- mutate(stat_n_year, indice = "N content")  %>%
  dplyr::rename(avg_val = `n_cont`)

#stat_ratio_year <- summarySE(env, measurevar = "fb_ratio", groupvars = "Year", conf.interval = 0.95, na.rm = TRUE)
stat_ratio_year <- group_by(env_w, Year) %>% #median, lepsi pro qPCR kvuli odlehlzm hodnotam
  summarise_each(funs(median(., na.rm = TRUE), max(., na.rm = TRUE), min(., na.rm = TRUE), mean(., na.rm = TRUE), sd(., na.rm = TRUE)), fb_ratio)
stat_ratio_year <- mutate(stat_ratio_year, indice = "F/B ratio (median)")  %>%
  dplyr::rename(avg_val = `median`) # nutne pro yobrayeni v ggplot, jedna se ale o median

stat_bac_year <- summarySE(env, measurevar = "bac_copy", groupvars = "Year", conf.interval = 0.95, na.rm = TRUE)
stat_bac_year <- mutate(stat_bac_year, indice = "Bac copy")  %>%
  dplyr::rename(avg_val = `bac_copy`) 

stat_fun_year <- summarySE(env, measurevar = "fun_copy", groupvars = "Year", conf.interval = 0.95, na.rm = TRUE)
stat_fun_year <- mutate(stat_fun_year, indice = "Fun copy")  %>%
  dplyr::rename(avg_val = `fun_copy`)


stat_ph_tree <- summarySE(env, measurevar = "pH", groupvars = "Tree", conf.interval = 0.95)
stat_ph_tree <- mutate(stat_ph_tree, indice = "pH")  %>%
  dplyr::rename(avg_val = `pH`) 

stat_c_tree <- summarySE(env, measurevar = "C", groupvars = "Tree", conf.interval = 0.95)
stat_c_tree <- mutate(stat_c_tree, indice = "C content")  %>%
  dplyr::rename(avg_val = `C`)

stat_n_tree <- dplyr::rename(env, n_cont = `N`) %>%
  summarySE(measurevar = "n_cont", groupvars = "Tree", conf.interval = 0.95)
stat_n_tree <- mutate(stat_n_tree, indice = "N content")  %>%
  dplyr::rename(avg_val = `n_cont`) 

#stat_ratio_tree <- summarySE(env, measurevar = "fb_ratio", groupvars = "Tree", conf.interval = 0.95, na.rm = TRUE)
stat_ratio_tree <- group_by(env_w, Tree) %>% #median, lepsi pro qPCR kvuli odlehlzm hodnotam
  summarise_each(funs(median(., na.rm = TRUE), max(., na.rm = TRUE), min(., na.rm = TRUE), mean(., na.rm = TRUE), sd(., na.rm = TRUE)), fb_ratio)
stat_ratio_tree <- mutate(stat_ratio_tree, indice = "F/B ratio (median)")  %>%
  dplyr::rename(avg_val = `median`) # nutne pro yobrayeni v ggplot, jedna se ale o median

chemics_all_year <- bind_rows(stat_ph_year, stat_c_year) %>% # spojeni statistiky dohromady
  bind_rows(stat_n_year) %>%
  bind_rows(stat_ratio_year)
chemics_all_tree <- bind_rows(stat_ph_tree, stat_c_tree) %>%
  bind_rows(stat_n_tree) %>%
  bind_rows(stat_ratio_tree)

chemics_graf_year <- ggplot(chemics_all_year, aes(Year, avg_val, fill = Year)) + # facet vrap graf pro chemii vyse
  geom_bar(stat = "identity", width = 0.5) +
  geom_errorbar(aes(ymin = avg_val - se, ymax = avg_val + se), colour = "black", width = .2) +
  facet_wrap(~ indice, scales = "free") +
  scale_fill_brewer(palette = "Set1", name = "Year") +
  scale_x_discrete(labels = c(">38", "38-16", "15-5", "<5")) +
  xlab("Age class") + 
  ylab("Values") +
  theme(legend.position = "none", axis.text.y = y_label, axis.text.x = x_label, axis.title = axis_titles, plot.title = plot_title, strip.text = element_text(size = 16), strip.background = element_blank())
chemics_graf_year
ggsave(plot = chemics_graf_year, filename = "chemics_graf_year.pdf", dpi = 300, height = 210, width = 297, units = "mm")

chemics_graf_tree <- ggplot(chemics_all_tree, aes(Tree, avg_val, fill = Tree)) + # facet vrap graf pro indices vyse
  geom_bar(stat = "identity", width = 0.5) +
  geom_errorbar(aes(ymin = avg_val - se, ymax = avg_val + se), colour = "black", width = .2) +
  facet_wrap(~ indice, scales = "free") +
  scale_fill_brewer(palette = "Set1", name = "Year") +
  scale_x_discrete(labels = c("FS", "AA", "PA")) +
  xlab("Tree") + 
  ylab("Values") +
  theme(legend.position = "none", axis.text.y = y_label, axis.text.x = x_label, axis.title = axis_titles, plot.title = plot_title, strip.text = element_text(size = 16), strip.background = element_blank())
chemics_graf_tree
ggsave(plot = chemics_graf_tree, filename = "chemics_graf_tree.pdf", dpi = 300, height = 210, width = 297, units = "mm")

##### scatterplot ####
scat_rat_n_fb <- ggplot(env_w, aes(x = N, y = fb_ratio, color = Year)) + 
  geom_point(shape = 19) + 
  stat_ellipse(na.rm = TRUE) + # 95% CI ellipse
  scale_fill_brewer(palette = "Set1", type = "div", name = "Age class", labels = c(">38", "38-16", "15-5", "<5")) +
  theme_bw() +
  theme(legend.position = "none", legend.text = legend_text, legend.title = legend_title, 
        axis.text.y = y_label, axis.text.x = x_label, axis.title = axis_titles, 
        plot.title = plot_title, strip.background = element_blank(), strip.text = element_blank(), panel.border = element_blank())
scat_rat_n_fb

scat_rat_c_n <- ggplot(env_w, aes(x = N, y = C, color = Year)) + 
  geom_point(shape = 19) + 
  stat_ellipse(na.rm = TRUE) + # 95% CI ellipse
  scale_fill_brewer(palette = "Set1", type = "div", name = "Age class", labels = c(">38", "38-16", "15-5", "<5")) +
  theme_bw() +
  theme(legend.position = "none", legend.text = legend_text, legend.title = legend_title, 
        axis.text.y = y_label, axis.text.x = x_label, axis.title = axis_titles, 
        plot.title = plot_title, strip.background = element_blank(), strip.text = element_blank(), panel.border = element_blank())
scat_rat_c_n

scat_rat_erg_fb <- ggplot(env_w, aes(x = erg, y = fb_ratio, color = Year)) + 
  geom_point(shape = 19) + 
  stat_ellipse(na.rm = TRUE) + # 95% CI ellipse
  scale_fill_brewer(palette = "Set1", type = "div", name = "Age class", labels = c(">38", "38-16", "15-5", "<5")) +
  theme_bw() +
  theme(legend.position = "none", legend.text = legend_text, legend.title = legend_title, 
        axis.text.y = y_label, axis.text.x = x_label, axis.title = axis_titles, 
        plot.title = plot_title, strip.background = element_blank(), strip.text = element_blank(), panel.border = element_blank())
scat_rat_erg_fb

scat_rat_erg_N <- ggplot(env_w, aes(x = erg, y = N, color = Year)) + 
  geom_point(shape = 19) + 
  stat_ellipse(na.rm = TRUE) + # 95% CI ellipse
  scale_fill_brewer(palette = "Set1", type = "div", name = "Age class", labels = c(">38", "38-16", "15-5", "<5")) +
  theme_bw() +
  theme(legend.position = "none", legend.text = legend_text, legend.title = legend_title, 
        axis.text.y = y_label, axis.text.x = x_label, axis.title = axis_titles, 
        plot.title = plot_title, strip.background = element_blank(), strip.text = element_blank(), panel.border = element_blank())
scat_rat_erg_N

scat_rat_lig_N <- ggplot(env_w, aes(x = lig, y = N, color = Year)) + 
  geom_point(shape = 19) + 
  stat_ellipse(na.rm = TRUE) + # 95% CI ellipse
  scale_fill_brewer(palette = "Set1", type = "div", name = "Age class", labels = c(">38", "38-16", "15-5", "<5")) +
  theme_bw() +
  theme(legend.position = "none", legend.text = legend_text, legend.title = legend_title, 
        axis.text.y = y_label, axis.text.x = x_label, axis.title = axis_titles, 
        plot.title = plot_title, strip.background = element_blank(), strip.text = element_blank(), panel.border = element_blank())
scat_rat_lig_N

scat_rat_lig_erg <- ggplot(env_w, aes(x = lig, y = erg, color = Year)) + 
  geom_point(shape = 19) + 
  stat_ellipse(na.rm = TRUE) + # 95% CI ellipse
  scale_fill_brewer(palette = "Set1", type = "div", name = "Age class", labels = c(">38", "38-16", "15-5", "<5")) +
  scale_colour_brewer(palette = "Set1" , name = "Age class", labels = c(">38", "38-16", "15-5", "<5")) +
  theme_bw() +
  theme(legend.position = "right", legend.text = legend_text, legend.title = legend_title, 
        axis.text.y = y_label, axis.text.x = x_label, axis.title = axis_titles, 
        plot.title = plot_title, strip.background = element_blank(), strip.text = element_blank(), panel.border = element_blank())
scat_rat_lig_erg

scat_rat_lig_water <- ggplot(env_w, aes(x = lig, y = water, color = Year)) + 
  geom_point(shape = 19) + 
  stat_ellipse(na.rm = TRUE) + # 95% CI ellipse
  scale_fill_brewer(palette = "Set1", type = "div", name = "Age class", labels = c(">38", "38-16", "15-5", "<5")) +
  theme_bw() +
  theme(legend.position = "none", legend.text = legend_text, legend.title = legend_title, 
        axis.text.y = y_label, axis.text.x = x_label, axis.title = axis_titles, 
        plot.title = plot_title, strip.background = element_blank(), strip.text = element_blank(), panel.border = element_blank())
scat_rat_lig_water

scat_rat_lig_dens <- ggplot(env_w, aes(x = lig, y = dens, color = Year)) + 
  geom_point(shape = 19) + 
  stat_ellipse(na.rm = TRUE) + # 95% CI ellipse
  scale_fill_brewer(palette = "Set1", type = "div", name = "Age class", labels = c(">38", "38-16", "15-5", "<5")) +
  theme_bw() +
  theme(legend.position = "none", legend.text = legend_text, legend.title = legend_title, 
        axis.text.y = y_label, axis.text.x = x_label, axis.title = axis_titles, 
        plot.title = plot_title, strip.background = element_blank(), strip.text = element_blank(), panel.border = element_blank())
scat_rat_lig_dens

scat_rat_erg_dens <- ggplot(env_w, aes(x = erg, y = dens, color = Year)) + 
  geom_point(shape = 19) + 
  stat_ellipse(na.rm = TRUE) + # 95% CI ellipse
  scale_fill_brewer(palette = "Set1", type = "div", name = "Age class", labels = c(">38", "38-16", "15-5", "<5")) +
  theme_bw() +
  theme(legend.position = "none", legend.text = legend_text, legend.title = legend_title, 
        axis.text.y = y_label, axis.text.x = x_label, axis.title = axis_titles, 
        plot.title = plot_title, strip.background = element_blank(), strip.text = element_blank(), panel.border = element_blank())
scat_rat_erg_dens

scat_rat_fb_water <- ggplot(env_w, aes(x = fb_ratio, y = water, color = Year)) + 
  geom_point(shape = 19) + 
  stat_ellipse(na.rm = TRUE) + # 95% CI ellipse
  scale_fill_brewer(palette = "Set1", type = "div", name = "Age class", labels = c(">38", "38-16", "15-5", "<5")) +
  theme_bw() +
  theme(legend.position = "none", legend.text = legend_text, legend.title = legend_title, 
        axis.text.y = y_label, axis.text.x = x_label, axis.title = axis_titles, 
        plot.title = plot_title, strip.background = element_blank(), strip.text = element_blank(), panel.border = element_blank())
scat_rat_fb_water

scatter_chemics <- grid.arrange(scat_rat_n_fb, scat_rat_c_n, scat_rat_erg_fb, scat_rat_erg_N, scat_rat_lig_N, scat_rat_lig_erg, scat_rat_lig_water, scat_rat_lig_dens, scat_rat_erg_dens)
ggsave(plot = scatter_chemics, filename = "scatter_chemics.pdf", dpi = 300, height = 210, width = 297, units = "mm", limitsize = FALSE)

##### Boxplots ####
env_w_spread <- gather(env_w, indice, value, 
                       c(`C`, `N`, `pH`, `water`, `lig`))

boxplot_chemics_c <- ggplot(env_w, aes(Year, C, fill = Tree)) +
  geom_boxplot() + 
  facet_wrap(~ Tree) +
  scale_fill_brewer(palette = "Set1", type = "div", name = "Tree species", labels = c("Fagus sp.", "Abies sp.", "Picea sp.")) +
  scale_x_discrete(labels = c(">38", "38-16", "15-5", "<5")) +
  labs(title = "Carbon content") +
  xlab("Age class") + 
  ylab("content [%]") +
  theme_bw() +
  # geom_jitter(width = 0.25, aes(colour = Tree)) +
  theme(legend.position = "none", legend.text = legend_text, legend.title = legend_title, 
        axis.text.y = y_label, axis.text.x = x_label, axis.title = axis_titles, 
        plot.title = plot_title, strip.background = element_blank(), strip.text = element_blank(), panel.border = element_blank())
boxplot_chemics_c

boxplot_chemics_n <- ggplot(env_w, aes(Year, N, fill = Tree)) +
  geom_boxplot() + 
  facet_wrap(~ Tree) +
  scale_fill_brewer(palette = "Set1", type = "div", name = "Tree species", labels = c("Fagus sp.", "Abies sp.", "Picea sp.")) +
  scale_x_discrete(labels = c(">38", "38-16", "15-5", "<5")) +
  labs(title = "Nitrogen content") +
  xlab("Age class") + 
  ylab("content [%]") +
  theme_bw() +
  # geom_jitter(width = 0.25, aes(colour = Tree)) +
  theme(legend.position = "right", legend.text = legend_text, legend.title = legend_title, 
        axis.text.y = y_label, axis.text.x = x_label, axis.title = axis_titles, 
        plot.title = plot_title, strip.background = element_blank(), strip.text = element_blank(), panel.border = element_blank())
boxplot_chemics_n

boxplot_chemics_ph <- ggplot(env_w, aes(Year, pH, fill = Tree)) +
  geom_boxplot() + 
  facet_wrap(~ Tree) +
  scale_fill_brewer(palette = "Set1", type = "div", name = "Tree species", labels = c("Fagus sp.", "Abies sp.", "Picea sp.")) +
  scale_x_discrete(labels = c(">38", "38-16", "15-5", "<5")) +
  labs(title = "pH") +
  xlab("Age class") + 
  ylab("value") +
  theme_bw() +
  # geom_jitter(width = 0.25, aes(colour = Tree)) +
  theme(legend.position = "none", legend.text = legend_text, legend.title = legend_title, 
        axis.text.y = y_label, axis.text.x = x_label, axis.title = axis_titles, 
        plot.title = plot_title, strip.background = element_blank(), strip.text = element_blank(), panel.border = element_blank())
boxplot_chemics_ph

boxplot_chemics_water <- ggplot(env_w, aes(Year, water, fill = Tree)) +
  geom_boxplot() + 
  facet_wrap(~ Tree) +
  scale_fill_brewer(palette = "Set1", type = "div", name = "Tree species", labels = c("Fagus sp.", "Abies sp.", "Picea sp.")) +
  scale_x_discrete(labels = c(">38", "38-16", "15-5", "<5")) +
  labs(title = "Water content") +
  xlab("Age class") + 
  ylab("content [%]") +
  theme_bw() +
  # geom_jitter(width = 0.25, aes(colour = Tree)) +
  theme(legend.position = "right", legend.text = legend_text, legend.title = legend_title, 
        axis.text.y = y_label, axis.text.x = x_label, axis.title = axis_titles, 
        plot.title = plot_title, strip.background = element_blank(), strip.text = element_blank(), panel.border = element_blank())
boxplot_chemics_water

boxplot_chemics_lig <- ggplot(env_w, aes(Year, lig, fill = Tree)) +
  geom_boxplot() + 
  facet_wrap(~ Tree) +
  scale_fill_brewer(palette = "Set1", type = "div", name = "Tree species", labels = c("Fagus sp.", "Abies sp.", "Picea sp.")) +
  scale_x_discrete(labels = c(">38", "38-16", "15-5", "<5")) +
  labs(title = "Lignin content") +
  xlab("Age class") + 
  ylab("content [%]") +
  theme_bw() +
  # geom_jitter(width = 0.25, aes(colour = Tree)) +
  theme(legend.position = "none", legend.text = legend_text, legend.title = legend_title, 
        axis.text.y = y_label, axis.text.x = x_label, axis.title = axis_titles, 
        plot.title = plot_title, strip.background = element_blank(), strip.text = element_blank(), panel.border = element_blank())
boxplot_chemics_lig

boxplot_chemics <- grid.arrange(boxplot_chemics_c, boxplot_chemics_n, boxplot_chemics_ph, boxplot_chemics_lig, ncol = 2, top = textGrob("Supplementary Figure 1", gp = gpar(cex = 0.7, fontface = "bold", col = "black")))
ggsave(boxplot_chemics, filename = "boxplot_chemics2.pdf", dpi = 300, height = 210, width = 297, units = "mm")

