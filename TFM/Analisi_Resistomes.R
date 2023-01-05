#### Resistome analysis ####
# Load data & packages
library(readxl)
library(ggplot2)
library(tibble)
library(car)
library(dplyr)
library(tidyr)
library(reshape2)
library(vegan)

LOGFILE_DNA_extractions <- read_excel("C:/Users/apico/ICRA/SCOREwater ICRA - Documents/#Metagenomics/LOGFILE_DNA_extractions.xlsx", 
                                      sheet = "Sequencing_statistics")
consum_sidiap <- read_excel("C:/Users/apico/ICRA/SCOREwater ICRA - Documents/#Metagenomics/consum_sidiap.xlsx")

## Sequencing statistics ####

# View(LOGFILE_DNA_extractions)
str(LOGFILE_DNA_extractions)
LOGFILE_DNA_extractions$Library_QC <- as.factor(LOGFILE_DNA_extractions$Library_QC)
LOGFILE_DNA_extractions$Loc <- as.factor(LOGFILE_DNA_extractions$Loc)
LOGFILE_DNA_extractions$Campaign <- as.factor(LOGFILE_DNA_extractions$Campaign)
LOGFILE_DNA_extractions$QC <- as.factor(LOGFILE_DNA_extractions$QC)
LOGFILE_DNA_extractions$SQ_set <- as.factor(LOGFILE_DNA_extractions$SQ_set)

summary(LOGFILE_DNA_extractions)
summary(LOGFILE_DNA_extractions$Library_QC)
LOGFILE_DNA_extractions<- na.omit(LOGFILE_DNA_extractions)


#### Grafic barres GC, Q30 i reads

ggplot(data=LOGFILE_DNA_extractions, aes(x=reorder(ID, GC), y=GC, fill=Loc)) +
  geom_bar(stat="identity")+
  facet_wrap(~Library_QC)+
  coord_flip()+
  geom_text(aes(label=QC), hjust=-0.3, size=2.8, col="black")+
  scale_fill_discrete(name="Barri", labels=c("Besos", "Carmel", "Poblenou", "Sant Gervasi"))+
  scale_y_continuous(expand=expansion(mult=c(0, .1)))+
  scale_x_discrete(name = "ID")+
  theme_classic()

ggplot(data=LOGFILE_DNA_extractions, aes(x=reorder(ID, Q30), y=Q30, fill=Loc)) +
  geom_bar(stat="identity")+
  facet_wrap(~Library_QC)+
  coord_flip()+
  geom_text(aes(label=QC), hjust=-0.3, size=2.8, col="black")+
  scale_fill_discrete(name="Barri", labels=c("Besos", "Carmel", "Poblenou", "Sant Gervasi"))+
  scale_y_continuous(expand=expansion(mult=c(0, .1)))+
  scale_x_discrete(name = "ID")+
  theme_classic()

ggplot(data=LOGFILE_DNA_extractions, aes(x=reorder(ID, Reads), y=Reads, fill=Loc)) +
  geom_bar(stat="identity")+
  facet_wrap(~Library_QC)+
  coord_flip()+
  geom_text(aes(label=QC), hjust=-0.3, size=2.8, col="black")+
  scale_fill_discrete(name="Barri", labels=c("Besos", "Carmel", "Poblenou", "Sant Gervasi"))+
  scale_y_continuous(expand=expansion(mult=c(0, .1)))+
  scale_x_discrete(name = "ID")+
  theme_classic()

## Diferències de reads i qualitat entre llocs

# boxplot(Reads ~ Loc, LOGFILE_DNA_extractions)
ARL <- aov(Reads ~ Loc, data=LOGFILE_DNA_extractions)
summary(ARL) # No hi ha diferències
# boxplot(Total_bp ~ Loc, LOGFILE_DNA_extractions)
ABL <- aov(Total_bp ~ Loc, data=LOGFILE_DNA_extractions)
summary(ABL) # No hi ha diferències
# boxplot(GC ~ Loc, LOGFILE_DNA_extractions)
AGCL <- aov(GC ~ Loc, data=LOGFILE_DNA_extractions)
summary(AGCL) # Hi ha diferències
TukeyHSD(AGCL) # SGVS == CRML, la resta són diferents

## Gràfic diferències GC per barri 
# ggplot(LOGFILE_DNA_extractions,aes(Loc,GC, fill=Loc))+ 
#   geom_boxplot()+
#   geom_line(data=tibble(x=c(1,2), y=c(60,60)),
#             aes(x=x, y=y), inherit.aes=FALSE)+
#   geom_text(data=tibble(x=1.5, y=60.2),
#             aes(x=x, y=y, label="*"), inherit.aes=FALSE)+
#   geom_line(data=tibble(x=c(1,3), y=c(61.5,61.5)),
#             aes(x=x, y=y), inherit.aes=FALSE)+
#   geom_text(data=tibble(x=2, y=61.7),
#             aes(x=x, y=y, label="*"), inherit.aes=FALSE)+
#   geom_line(data=tibble(x=c(1,4), y=c(63,63)),
#             aes(x=x, y=y), inherit.aes=FALSE)+
#   geom_text(data=tibble(x=2.5, y=63.2),
#             aes(x=x, y=y, label="*"), inherit.aes=FALSE)+
#   geom_line(data=tibble(x=c(2,3), y=c(52,52)),
#             aes(x=x, y=y), inherit.aes=FALSE)+
#   geom_text(data=tibble(x=2.5, y=52.2),
#             aes(x=x, y=y, label="*"), inherit.aes=FALSE)+
#   geom_line(data=tibble(x=c(3,4), y=c(54,54)),
#             aes(x=x, y=y), inherit.aes=FALSE)+
#   geom_text(data=tibble(x=3.5, y=54.2),
#             aes(x=x, y=y, label="*"), inherit.aes=FALSE)+
#   labs(x="Location", y="GC %", title="A")+
#   scale_fill_discrete(name="Location", labels=c("Besos", "Carmel", "Poblenou", "Sant Gervasi"))

shapiro.test(x=residuals(object=AGCL)) # Es compleix normalitat
leveneTest(AGCL) # Es compleix homosadesticitat

# boxplot(Q30 ~ Loc, LOGFILE_DNA_extractions)
AQL <- aov(Q30 ~ Loc, data=LOGFILE_DNA_extractions)
summary(AQL) # No hi ha diferències

## Diferències de reads i qualitat entre campanyes

# boxplot(Reads ~ Campaign, LOGFILE_DNA_extractions)
summary(aov(Reads ~ Campaign, LOGFILE_DNA_extractions)) # No hi ha diferències
# boxplot(Total_bp ~ Campaign, LOGFILE_DNA_extractions)
summary(aov(Total_bp ~ Campaign, LOGFILE_DNA_extractions)) # No hi ha diferències
# boxplot(GC ~ Campaign, LOGFILE_DNA_extractions)
summary(aov(GC ~ Campaign, LOGFILE_DNA_extractions)) # No hi ha diferències
# boxplot(Q30 ~ Campaign, LOGFILE_DNA_extractions)
summary(aov(Q30 ~ Campaign, LOGFILE_DNA_extractions)) # No hi ha diferències

## Diferències de reads i qualitat segons si han passat el control de qualitat
# boxplot(Reads ~ QC, LOGFILE_DNA_extractions)
ARQ <- aov(Reads ~ QC, data=LOGFILE_DNA_extractions)
summary(ARQ) # No hi ha diferències
# boxplot(Total_bp ~ QC, LOGFILE_DNA_extractions)
ABQ <- aov(Total_bp ~ QC, data=LOGFILE_DNA_extractions)
summary(ABQ) # No hi ha diferències
# boxplot(GC ~ QC, LOGFILE_DNA_extractions)
AGCQ <- aov(GC ~ QC, data=LOGFILE_DNA_extractions)
summary(AGCQ) # Les que han passat el control de qualitat tenen significativament més GC

## Gràfic diferències de GC en funció del control de qualitat
# ggplot(LOGFILE_DNA_extractions,aes(QC,GC, fill=QC))+ 
#   geom_boxplot()+
#   geom_line(data=tibble(x=c(1,2), y=c(60,60)),
#             aes(x=x, y=y), inherit.aes=FALSE)+
#   geom_text(data=tibble(x=1.5, y=60.2),
#             aes(x=x, y=y, label="*"), inherit.aes=FALSE)+
#   labs(x="Quality Control", y="GC %", title="B")+
#   scale_fill_discrete(name="Quality Control", labels=c("No", "Yes"))

shapiro.test(x=residuals(object=AGCQ))
leveneTest(AGCQ) #No homogeneitat, repetim el test sense assumir variances iguals
oneway.test(GC ~ QC, data=LOGFILE_DNA_extractions)


# boxplot(Q30 ~ QC, LOGFILE_DNA_extractions)
AQQ <- aov(Q30 ~ QC, data=LOGFILE_DNA_extractions)
summary(AQQ) # Les que no han passat el primer control tenen significativament més qualitat

## Gràfic diferències Q30 segons si han passat el primer control de qualitat
# ggplot(LOGFILE_DNA_extractions,aes(QC,Q30, fill=QC))+ 
#   geom_boxplot()+
#   geom_line(data=tibble(x=c(1,2), y=c(95,95)),
#             aes(x=x, y=y), inherit.aes=FALSE)+
#   geom_text(data=tibble(x=1.5, y=95.2),
#             aes(x=x, y=y, label="*"), inherit.aes=FALSE)+
#   labs(x="Quality Control", y="Q30 %", title="A")+
#   scale_fill_discrete(name="Quality Control", labels=c("No", "Yes"))

shapiro.test(x=residuals(object=AQQ)) # No normalitat repetim el test fent-ne un de no paramètric
leveneTest(AQQ)
kruskal.test(Q30 ~ QC, data=LOGFILE_DNA_extractions)

## Diferències de reads i qualitat segons el tipus de llibreria

# boxplot(Reads ~ Library_QC, LOGFILE_DNA_extractions)
ARL <- aov(Reads ~ Library_QC, data=LOGFILE_DNA_extractions)
summary(ARL) # No hi ha diferències
# boxplot(Total_bp ~ Library_QC, LOGFILE_DNA_extractions)
ABL <- aov(Total_bp ~ Library_QC, data=LOGFILE_DNA_extractions)
summary(ABL) # No hi ha diferències
# boxplot(GC ~ Library_QC, LOGFILE_DNA_extractions)
AGCL <- aov(GC ~ Library_QC, data=LOGFILE_DNA_extractions)
summary(AGCL) # Llibreria TS té significativament més GC

## Gràfic GC en funció de la llibreria 
# ggplot(LOGFILE_DNA_extractions,aes(Library_QC,GC, fill=Library_QC))+ 
#   geom_boxplot()+
#   geom_line(data=tibble(x=c(1,2), y=c(60,60)),
#             aes(x=x, y=y), inherit.aes=FALSE)+
#   geom_text(data=tibble(x=1.5, y=60.2),
#             aes(x=x, y=y, label="*"), inherit.aes=FALSE)+
#   labs(x="Library", y="GC", title="C")+
#   scale_fill_discrete(name="Library", labels=c("TS", "TSF"))

shapiro.test(x=residuals(object=AGCL)) # No hi ha normalitat, repetim test no paramètric
leveneTest(AGCL)
kruskal.test(GC ~ Library_QC, data=LOGFILE_DNA_extractions)

# boxplot(Q30 ~ Library_QC, LOGFILE_DNA_extractions)
AQL <- aov(Q30 ~ Library_QC, data=LOGFILE_DNA_extractions)
summary(AQL) # Llibreria TSF té significativament més qualitat.

## Gràfic diferències en Q30 segons llibreria
# ggplot(LOGFILE_DNA_extractions,aes(Library_QC,Q30, fill=Library_QC))+ 
#   geom_boxplot()+
#   geom_line(data=tibble(x=c(1,2), y=c(95,95)),
#             aes(x=x, y=y), inherit.aes=FALSE)+
#   geom_text(data=tibble(x=1.5, y=95.2),
#             aes(x=x, y=y, label="*"), inherit.aes=FALSE)+
#   labs(x="Library", y="Q30 %", title="B")+
#   scale_fill_discrete(name="Library", labels=c("TS", "TSF"))

shapiro.test(x=residuals(object=AQL)) # No normalitat repetim test no paramètric
leveneTest(AQL) # No homosadesticitat
kruskal.test(Q30 ~ Library_QC, data=LOGFILE_DNA_extractions)


## All sequences analysis ####

## General overview

# names(LOGFILE_DNA_extractions)

# boxplot(Bacteria ~ Loc, LOGFILE_DNA_extractions)
# boxplot(Bacteria ~ Campaign, LOGFILE_DNA_extractions)
# boxplot(Bacteria ~ Library_QC, LOGFILE_DNA_extractions)
# boxplot(Bacteria ~ SQ_set, LOGFILE_DNA_extractions)
# 
# boxplot(resistance_genes ~ Loc, LOGFILE_DNA_extractions)
# boxplot(resistance_genes ~ Campaign, LOGFILE_DNA_extractions)
# boxplot(resistance_genes ~ Library_QC, LOGFILE_DNA_extractions)
# boxplot(resistance_genes ~ SQ_set, LOGFILE_DNA_extractions)
# 
# boxplot(norm_resistance_genes ~ Loc, LOGFILE_DNA_extractions)
RGL <- aov(norm_resistance_genes ~ Loc, data=LOGFILE_DNA_extractions)
summary(RGL)

shapiro.test(x=residuals(object=RGL)) 
leveneTest(RGL) # No homosadesticitat, repetim el test sense assumir igualtat de variances
oneway.test(norm_resistance_genes ~ Loc, data=LOGFILE_DNA_extractions)
pairwise.t.test(LOGFILE_DNA_extractions$norm_resistance_genes, LOGFILE_DNA_extractions$Loc,
                p.adjust.method = "BH", pool.sd = FALSE) # SGVS significativament major que PBN. La resta no presenten diferències

ggplot(LOGFILE_DNA_extractions,aes(Loc,norm_resistance_genes, fill=Loc))+ 
  geom_boxplot()+
  geom_line(data=tibble(x=c(3,4), y=c(0.51, 0.51)),
            aes(x=x, y=y), inherit.aes=FALSE)+
  geom_text(data=tibble(x=3.5, y=0.512, Campaign="Mar"),
            aes(x=x, y=y, label="*"), inherit.aes=FALSE)+
  labs(x="Barri", y="ARG Totals/16s")+
  scale_fill_discrete(name="Barri", labels=c("Besos", "Carmel", "Poblenou", "Sant Gervasi"))+
  theme_classic()

## Per campanyes
Mar <- subset(LOGFILE_DNA_extractions, Campaign=="Mar")
Nov <- subset(LOGFILE_DNA_extractions, Campaign=="Nov")

AMar <- aov(norm_resistance_genes ~ Loc, data=Mar)
summary(AMar)
intervalsAMar <- TukeyHSD(AMar)
shapiro.test(x=residuals(object=AMar)) # Sant Gervasi diferent de poblenou
leveneTest(AMar) # No homosadesticitat
oneway.test(norm_resistance_genes ~ Loc, data=Mar)
pairwise.t.test(Mar$norm_resistance_genes, Mar$Loc, p.adjust.method="BH", pool.sd=FALSE)

ANov <- aov(norm_resistance_genes ~ Loc, data=Nov)
summary(ANov)
shapiro.test(x=residuals(object=ANov))
leveneTest(ANov) # No homosadesticitat
oneway.test(norm_resistance_genes ~ Loc, data=Nov)

ggplot(LOGFILE_DNA_extractions,aes(Loc,norm_resistance_genes, fill=Loc))+ 
  geom_boxplot()+
  facet_grid(~factor(Campaign, levels=c("Nov", "Mar")))+
  geom_line(data=tibble(x=c(3,4), y=c(0.51, 0.51), Campaign="Mar"),
            aes(x=x, y=y), inherit.aes=FALSE)+
  geom_text(data=tibble(x=3.5, y=0.512, Campaign="Mar"),
            aes(x=x, y=y, label="*"), inherit.aes=FALSE)+
  labs(x="Barri", y="ARG Totals/16s")+
  scale_fill_discrete(name="Barri", labels=c("Besos", "Carmel", "Poblenou", "Sant Gervasi"))+
  theme_classic()


## Classes de ARG ##### 

# Famílies que ens interessen (cal sumar diferent el grup others)

resistome <- data.frame(LOGFILE_DNA_extractions$Sample,
                        LOGFILE_DNA_extractions$Loc,
                        LOGFILE_DNA_extractions$Campaign,
                        LOGFILE_DNA_extractions$Day,
                        LOGFILE_DNA_extractions$Bacteria,
                        LOGFILE_DNA_extractions$Other,
                        LOGFILE_DNA_extractions$vancomycin,
                        LOGFILE_DNA_extractions$fosmidomycin,
                        LOGFILE_DNA_extractions$rifampin,
                        LOGFILE_DNA_extractions$fosfomycin,
                        LOGFILE_DNA_extractions$glycopeptide,
                        LOGFILE_DNA_extractions$trimethoprim,
                        LOGFILE_DNA_extractions$polymyxin,
                        LOGFILE_DNA_extractions$aminocoumarin,
                        LOGFILE_DNA_extractions$phenicol,
                        LOGFILE_DNA_extractions$mupirocin,
                        LOGFILE_DNA_extractions$sulfonamide,
                        LOGFILE_DNA_extractions$peptide,
                        LOGFILE_DNA_extractions$fluoroquinolone,
                        LOGFILE_DNA_extractions$aminoglycoside,
                        LOGFILE_DNA_extractions$beta_lactam,
                        LOGFILE_DNA_extractions$tetracycline,
                        LOGFILE_DNA_extractions$multidrug,
                        LOGFILE_DNA_extractions$MLS)
                        

names(resistome) <- c("ID", 
                      "Loc",
                      "Campaign",
                      "Day",
                      "Bacteria",
                      "Other",
                      "vancomycin",
                      "fosmidomycin",
                      "rifampin",
                      "fosfomycin",
                      "glycopeptide",
                      "trimethoprim",
                      "polymyxin",
                      "aminocoumarin",
                      "phenicol",
                      "mupirocin",
                      "sulfonamide",
                      "peptide",
                      "fluoroquinolone",
                      "aminoglycoside",
                      "beta_lactam",
                      "tetracycline",
                      "multidrug",
                      "MLS")

resistome2 <- resistome[,-1]
resistome2 <- resistome2[,-2:-3]
resistome2 <- na.omit(resistome2)

resistome2[,3:21] <- resistome2[,3:21]/resistome2[,2]
resistome2 <- resistome2[,-2]

resistome2 <- melt(resistome2)

ggplot(resistome2, aes(Loc, variable, fill= log10(value))) + 
  geom_tile()+
  scale_fill_distiller(palette = "RdYlBu", name="Log10(ARG reads/16s reads)")+
  scale_x_discrete(name="Barris", expand=c(0,0))+
  scale_y_discrete(name="Classes ARG")+
  theme(axis.text.x= element_text(size=8, angle=90))+
  labs(x="Location", y="ARG class")+
  theme_classic()

## SIDIAP ####

# Creating variable nddd/hab
consum_sidiap$nddd_hab <- consum_sidiap$nddd/consum_sidiap$pob_SIDIAP_mensual

# Filtering antibiotics
antibiotics_consum <- subset(consum_sidiap, Therapeutic_Group == "Antibiotics")

# Set zona as a factor for ploting
#str(antibiotics_consum)
antibiotics_consum$Zona <- as.factor(antibiotics_consum$Zona)


# Total antibiotic consumption per year (each month as a sample)
total_antibiotics <- subset(antibiotics_consum, Compound == "Qualsevol_antibiotic")
# boxplot(nddd_hab ~ Zona, data=total_antibiotics, ylab="nddd/(inhabitant · month)")

anova <- aov(nddd_hab ~ Zona, data=total_antibiotics)
summary(anova)
shapiro.test(anova$residuals)
leveneTest(anova)
TukeyHSD(anova) #PBN significativament més baix

colors <- c("#7CAE00", "#00BFC4", "#C77CFF")
ggplot(data=total_antibiotics, aes(x=Zona, y=nddd_hab, fill=Zona))+
  geom_boxplot()+
  xlab("Barri")+
  ylab("nddd/(inhabitant · month)")+
  scale_x_discrete(name="Barri", labels=c("CRML", "PBN","SGVS"))+
  scale_fill_manual(name="Barri", values=colors)+
  geom_line(data=tibble(x=c(1,2), y=c(0.51,0.51)),
            aes(x=x, y=y), inherit.aes=FALSE)+
  geom_text(data=tibble(x=1.5, y=0.515),
            aes(x=x, y=y, label="*"), inherit.aes=FALSE)+
  geom_line(data=tibble(x=c(2,3), y=c(0.52,0.52)),
            aes(x=x, y=y), inherit.aes=FALSE)+
  geom_text(data=tibble(x=2.5, y=0.525),
            aes(x=x, y=y, label="*"), inherit.aes=FALSE)+
  theme_classic()

## DIVERSITY ####

mat <- resistome[,6:14]/resistome[,5]

sppr <- diversity(mat)

diversity <- data.frame(sppr, resistome$Loc)

# boxplot(diversity$sppr ~ diversity$resistome.Loc)
# str(diversity)

divaov <- aov(diversity$sppr ~ diversity$resistome.Loc)
summary(aov(diversity$sppr ~ diversity$resistome.Loc))

## Índex de Shannon-Weaver segons el barri
# ggplot(diversity,aes(x=resistome.Loc, y=sppr, fill=resistome.Loc))+ 
#   geom_boxplot()+
#   geom_line(data=tibble(x=c(1,2), y=c(1.57, 1.57)),
#             aes(x=x, y=y), inherit.aes=FALSE)+
#   geom_text(data=tibble(x=1.5, y=1.58),
#             aes(x=x, y=y, label="*"), inherit.aes=FALSE)+
#   geom_line(data=tibble(x=c(1,3), y=c(1.6, 1.6)),
#             aes(x=x, y=y), inherit.aes=FALSE)+
#   geom_text(data=tibble(x=2, y=1.61),
#             aes(x=x, y=y, label="*"), inherit.aes=FALSE)+
#   geom_line(data=tibble(x=c(2,4), y=c(1.66, 1.66)),
#             aes(x=x, y=y), inherit.aes=FALSE)+
#   geom_text(data=tibble(x=3, y=1.67),
#             aes(x=x, y=y, label="*"), inherit.aes=FALSE)+
#   geom_line(data=tibble(x=c(1,4), y=c(1.69, 1.69)),
#             aes(x=x, y=y), inherit.aes=FALSE)+
#   geom_text(data=tibble(x=2.5, y=1.70),
#             aes(x=x, y=y, label="*"), inherit.aes=FALSE)+
#   geom_line(data=tibble(x=c(2,3), y=c(1.63, 1.63)),
#             aes(x=x, y=y), inherit.aes=FALSE)+
#   geom_text(data=tibble(x=2.5, y=1.64),
#             aes(x=x, y=y, label="*"), inherit.aes=FALSE)+
#   labs(x="Location", y="Shannon-Weaver Index")+
#   scale_fill_discrete(name="Location", labels=c("Besos", "Carmel", "Poblenou", "Sant Gervasi"))

shapiro.test(x=residuals(object=divaov)) 
leveneTest(divaov) # No hi ha homosadesticitat, repetim el test sense assumir variances iguals
oneway.test(sppr ~ resistome.Loc, data=diversity)
pairwise.t.test(diversity$sppr, diversity$resistome.Loc,
                p.adjust.method = "BH", pool.sd = FALSE)

