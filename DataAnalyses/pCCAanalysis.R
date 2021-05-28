#reset
rm(list=ls())

#read data objects from other script
source("preparingPhyloseq_VAMB.R")
#source("preparingPhyloseq.R")

### PCoA ###
#bact_dist <- vegdist(species_abundOnly,  method = "bray")
#library("ape")
#bact_pcoa <- pcoa(bact_dist)
#first50_PCs <- bact_pcoa$vectors[,1:50]

#first axis explain 15% of the variation:
#barplot(bact_pcoa$values$Relative_eig[1:10])




### CCA & pCCA ###

# Readying the tables #
#X is community data matrix (CRISPR abundances)
x_crispr <- data.frame(crispr_abundOnly)

#Y is constraining matrix (clinical variables)
#getting rid of duplicated samples for same patient (same ABCno)
y_clinical <- get_variable(newCRISPRs)[!duplicated(get_variable(newCRISPRs)$ABCno),]
#and ordering by ABCno like the other tables
y_clinical <- y_clinical[order(y_clinical$ABCno),]

#Z is the first 50 axes from the bray-curtis PCoA ordination
#z_bact <- first50_PCs
#or if the want to test bacteria abundances without PCoA
z_bact <- species_abundOnly
#or if genus only
#z_bact <- genus_abundOnly

# fish oil #
#Let's try making for fish oil consumption
y_clinical %>% group_by(Supplement_oil) %>% summarise(N=n())
#We see there's 1 NA
#making an index that will leave out the sample(s) with the NA in
index <- !is.na(y_clinical$Supplement_oil)
x_crispr_fish <- x_crispr[index,]
y_clinical_fish <- y_clinical[index,]
z_bact_fish <- z_bact[index,]

#running the pCCA on fish oil
CRISPR_fish_pcca <- cca(x_crispr_fish ~ Supplement_oil + Condition(z_bact_fish), data = y_clinical_fish)
CRISPR_fish_pcca
#with VAMB on species level: bacteria explains 73.0%, fish oil 0.07 %
#with VAMB on genus level: bacteria explains 39.02%, fish oil 0.20%
#alias(CRISPR_fish_pcca, names=TRUE)

#plotting
plot(CRISPR_fish_pcca)

#without PCoA, species level:
df <- data.frame(
  Group = c("Bacteria (73.02%)", "Fish oil (0.07%)", "Unexplained (26.91%)"),
  value = c(73.02, 0.07, 26.91)
                 )

bp_fish<- ggplot(df, aes(x="", y=value, fill=Group))+
  geom_bar(width = 1, stat = "identity")

bp_fish + 
  coord_polar("y", start=0) +
  theme(axis.ticks = element_blank(),
        axis.title = element_blank()) +
  labs(title="Variation explained by bacteria and fish oil in pCCA analysis", 
       subtitle="Species level (ANOVA test P-value: 0.901)",
       caption = "Figure B")

ggsave("plots/special/VAMB_pCCA_fishoil_species.png", width = 6.5, height =5)

#without PCoA, genus level
df <- data.frame(
  Group = c("Bacteria (39.02%)", "Fish oil (0.20%)", "Unexplained (60.78%)"),
  value = c(39.02, 0.20, 60.78)
)

bp_fish<- ggplot(df, aes(x="", y=value, fill=Group))+
  geom_bar(width = 1, stat = "identity")

bp_fish + 
  coord_polar("y", start=0) +
  theme(axis.ticks = element_blank(),
        axis.title = element_blank()) +
  labs(title="Variation explained by bacteria and fish oil in pCCA analysis", 
       subtitle="Genus level (ANOVA test P-value: 0.088)",
       caption = "Figure A") 

ggsave("plots/special/VAMB_pCCA_fishoil_genus.png", width = 6.5, height =5)



# Vitamin D #
#Let's try making for vitamin D consumption
y_clinical %>% group_by(Supplement_vitamin) %>% summarise(N=n())
#We see there's 104 NAs
#making an index that will leave out the sample(s) with the NA in
index <- !is.na(y_clinical$Supplement_vitamin)
x_crispr_vit <- x_crispr[index,]
y_clinical_vit <- y_clinical[index,]
z_bact_vit <- z_bact[index,]

#running the pCCA on vitamin D
CRISPR_vit_pcca <- cca(x_crispr_vit ~ Supplement_vitamin + Condition(z_bact_vit), data = y_clinical_vit)
CRISPR_vit_pcca
#alias(CRISPR_vit_pcca, names=TRUE)
#VAMB: bacteria explains 83.43%, vitamin D 0.018%



# Delivery #
#Let's try making for delivery method
y_clinical %>% group_by(Delivery) %>% summarise(N=n())
#We see there's 0 NAs!
#no need to make indexes then

#running the pCCA on delivery
CRISPR_del_pcca <- cca(x_crispr ~ Delivery + Condition(z_bact), data = y_clinical)
CRISPR_del_pcca
#alias(CRISPR_del_pcca, names=TRUE)
#VAMB: bacteria 72.94%, delivery 0.08



# Animal days #
#Let's try making for amount of days spent with animals. (less than or equal to 7 days spent with animal = 0)
y_clinical %>% group_by(Furred_animal_days) %>% summarise(N=n())
#We see no NAs
#so no need to make index

#running the pCCA on days with furred animals
CRISPR_ani_pcca <- cca(x_crispr ~ Furred_animal_days + Condition(z_bact), data = y_clinical)
CRISPR_ani_pcca
#alias(CRISPR_ani_pcca, names=TRUE)
#VAMB: bacteria is 72.94%, furred animal days 0.12%


## Asthma ##
# Cross #
y_clinical %>% group_by(j45_5yr_cross) %>% summarise(N=n())
#We see there's 24 NAs
#making an index that will leave out the sample(s) with the NA in
index <- !is.na(y_clinical$j45_5yr_cross)
x_crispr_cross <- x_crispr[index,]
y_clinical_cross <- y_clinical[index,]
z_bact_cross <- z_bact[index,]

#running the pCCA on asthma
CRISPR_cross_pcca <- cca(x_crispr_cross ~ j45_5yr_cross + Condition(z_bact_cross), data = y_clinical_cross)
CRISPR_cross_pcca
#VAMB: bacteria explains 75.13%, asthma cross 0.08%
#alias(CRISPR_cross_pcca, names=TRUE)

#plotting
df <- data.frame(
  Group = c("Bacteria (75.12%)", "Cross - asthma at age 5 (0.08%)", "Unexplained (24.80%)"),
  value = c(75.12, 0.08, 24.80))

bp_cross<- ggplot(df, aes(x="", y=value, fill=Group))+
  geom_bar(width = 1, stat = "identity")

bp_cross + 
  coord_polar("y", start=0) + 
  theme(axis.ticks = element_blank(),
        axis.title = element_blank()) +
  labs(title="Variation explained by bacteria and asthma (at age 5) in pCCA analysis", 
       subtitle="Species level (ANOVA test P-value: 0.89)")
ggsave("plots/special/VAMB_pCCA_asthmacross.png", width = 7, height =5)


# Ever #
y_clinical %>% group_by(j45_5yr_ever) %>% summarise(N=n())
#We see there's 24 NAs
#making an index that will leave out the sample(s) with the NA in
index <- !is.na(y_clinical$j45_5yr_ever)
x_crispr_ever <- x_crispr[index,]
y_clinical_ever <- y_clinical[index,]
z_bact_ever <- z_bact[index,]

#running the pCCA on asthma
CRISPR_ever_pcca <- cca(x_crispr_ever ~ j45_5yr_ever + Condition(z_bact_ever), data = y_clinical_ever)
CRISPR_ever_pcca
#alias(CRISPR_ever_pcca, names=TRUE)
#VAMB: bacteria is 75.13%, ever is 0.1%

# Asthmatic mother #
y_clinical %>% group_by(Asthmatic_mother) %>% summarise(N=n())
#We see there's 2 NAs
#making a index that will leave out the sample(s) with the NA in
index <- !is.na(y_clinical$Asthmatic_mother)
x_crispr_astMother <- x_crispr[index,]
y_clinical_astMother <- y_clinical[index,]
z_bact_astMother <- z_bact[index,]

#running the pCCA on house location
CRISPR_astMother_pcca <- cca(x_crispr_astMother ~ Asthmatic_mother + Condition(z_bact_astMother), data = y_clinical_astMother)
CRISPR_astMother_pcca
#alias(CRISPR_astMother_pcca, names=TRUE)
#VAMB: bacteria is 73.07%, asthmatic mother is 0.08%


# PAM clusters #
#Let's try making for pamclusters
y_clinical %>% group_by(pamcluster) %>% summarise(N=n())
#We see there's 38 NA
#making a index that will leave out the sample(s) with the NA in
index <- !is.na(y_clinical$pamcluster)
x_crispr_pam <- x_crispr[index,]
y_clinical_pam <- y_clinical[index,]
z_bact_pam <- z_bact[index,]

#running the pCCA on pam clusters
CRISPR_pam_pcca <- cca(x_crispr_pam ~ pamcluster + Condition(z_bact_pam), data = y_clinical_pam)
CRISPR_pam_pcca
#alias(CRISPR_pam_pcca, names=TRUE)
#VAMB: bacteria is 76.42%, pam cluster is 0.25%

#plotting
plot(CRISPR_pam_pcca)

#without PCoA:
df <- data.frame(
  group = c("Bacteria (76.42%)", "PAM cluster (0.26%)", "Unexplained (23.32%)"),
  value = c(76.42, 0.26, 23.32)
)
bp_pam <- ggplot(df, 
                 aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity")

bp_pam + 
  coord_polar("y", start=0) + 
  theme(axis.ticks = element_blank(),
        axis.title = element_blank()) +
  labs(title="Variation explained by bacteria and PAM clusters in pCCA analysis", 
       subtitle="Species level (ANOVA test P-value: 0.052)")
ggsave("plots/special/VAMB_pCCA_pam_species.png", width = 7, height =5)



# Rural / urban #
y_clinical %>% group_by(RuralUrbanStatus) %>% summarise(N=n())
#We see there's 1 NA
#making a index that will leave out the sample(s) with the NA in
index <- !is.na(y_clinical$RuralUrbanStatus)
x_crispr_urb <- x_crispr[index,]
y_clinical_urb <- y_clinical[index,]
z_bact_urb <- z_bact[index,]

#running the pCCA on house location
CRISPR_urb_pcca <- cca(x_crispr_urb ~ RuralUrbanStatus + Condition(z_bact_urb), data = y_clinical_urb)
CRISPR_urb_pcca
#alias(CRISPR_urb_pcca, names=TRUE)
#VAMB: bacteria is 73.13%, rural urban is 0.18%


## Antibiotics ##
# Antibiotics given to child at birth #
y_clinical %>% group_by(Antibiotics_birth_child) %>% summarise(N=n())
#We see there's 3 NAs
#making a index that will leave out the sample(s) with the NA in
index <- !is.na(y_clinical$Antibiotics_birth_child)
x_crispr_antibChild <- x_crispr[index,]
y_clinical_antibChild <- y_clinical[index,]
z_bact_antibChild <- z_bact[index,]

#running the pCCA 
CRISPR_antibChild_pcca <- cca(x_crispr_antibChild ~ Antibiotics_birth_child + Condition(z_bact_antibChild), data = y_clinical_antibChild)
CRISPR_antibChild_pcca
#alias(CRISPR_antibChild_pcca, names=TRUE)
#VAMB: bacteria is 73.12%, antibiotics at birth is 0.08%

# Antibiotics given to mother during labor (removing those who had sectio performed since these are always given antibiotics)
y_clinical %>% group_by(Antibiotics_birth_mother) %>% summarise(N=n())
#We see there's 3 NAs
#making a index that will leave out the sample(s) with the NA in
index <- !is.na(y_clinical$Antibiotics_birth_mother) & y_clinical$Delivery==0 
x_crispr_antibMother <- x_crispr[index,]
y_clinical_antibMother <- y_clinical[index,]
z_bact_antibMother <- z_bact[index,]

#running the pCCA 
CRISPR_antibMother_pcca <- cca(x_crispr_antibMother ~ Antibiotics_birth_mother + Condition(z_bact_antibMother), data = y_clinical_antibMother)
CRISPR_antibMother_pcca
#alias(CRISPR_antibMother_pcca, names=TRUE)
#VAMB: bacteria is 88.83%, antibiotics to mother is 0.11%

# Antibiotics given to child at any point up until 1 yr #
y_clinical %>% group_by(Antibiotics_ever_1yr_child) %>% summarise(N=n())
#We see there's 6 NAs
#making a index that will leave out the sample(s) with the NA in
index <- !is.na(y_clinical$Antibiotics_ever_1yr_child)
x_crispr_antibChildEver <- x_crispr[index,]
y_clinical_antibChildEver <- y_clinical[index,]
z_bact_antibChildEver <- z_bact[index,]

#running the pCCA 
CRISPR_antibChildEver_pcca <- cca(x_crispr_antibChildEver ~ Antibiotics_ever_1yr_child + Condition(z_bact_antibChildEver), data = y_clinical_antibChildEver)
CRISPR_antibChildEver_pcca
#alias(CRISPR_antibChildEver_pcca, names=TRUE)
#VAMB: bacteria is 73.50%, antibiotics to child ever is 0.08%

df <- data.frame(
  Group = c("Bacteria (73.50%)", "Child given AB 1 year up \nuntil sampling (0.08%)", "Unexplained (26.42%)"),
  value = c(73.50, 0.08, 26.42))

bp_antib1yr <- ggplot(df, aes(x="", y=value, fill=Group))+
  geom_bar(width = 1, stat = "identity")

bp_antib1yr + 
  coord_polar("y", start=0) +
  theme(axis.ticks = element_blank(),
        axis.title = element_blank()) +
  labs(title="Variation explained by bacteria and antibiotics given 1 year \nup until sampling in pCCA analysis", 
       subtitle="Species level (ANOVA test P-value: 0.833)",
       caption="Figure A")

ggsave("plots/special/VAMB_pCCA_antib1yr.png", width = 7, height =5)

# Antibiotics given to child 1 wk up to sampling #
y_clinical %>% group_by(ab_child_1wk_before) %>% summarise(N=n())
#We see there's 6 NAs
#making a index that will leave out the sample(s) with the NA in
index <- !is.na(y_clinical$ab_child_1wk_before)
x_crispr_antibChild1wk <- x_crispr[index,]
y_clinical_antibChild1wk <- y_clinical[index,]
z_bact_antibChild1wk <- z_bact[index,]

#running the pCCA 
CRISPR_antibChild1wk_pcca <- cca(x_crispr_antibChild1wk ~ ab_child_1wk_before + Condition(z_bact_antibChild1wk), data = y_clinical_antibChild1wk)
CRISPR_antibChild1wk_pcca
#alias(CRISPR_antibChild1wk_pcca, names=TRUE)
#VAMB: bacteria is 73.50%, antibiotics to child 1 wk before sampling is 0.13%
#VAMB genus level 37.39%, antibiotics 0.25%

df <- data.frame(
  Group = c("Bacteria (73.50%)", "Child given AB 1 week up \nuntil sampling (0.14%)", "Unexplained (26.36%)"),
  value = c(73.50, 0.14, 26.36))

bp_antib1wk <- ggplot(df, aes(x="", y=value, fill=Group))+
  geom_bar(width = 1, stat = "identity")

bp_antib1wk + 
  coord_polar("y", start=0) +
  theme(axis.ticks = element_blank(),
        axis.title = element_blank()) +
  labs(title="Variation explained by bacteria and antibiotics given 1 week \nup until sampling in pCCA analysis", 
       subtitle="Species level (ANOVA test P-value: 0.391)",
       caption="Figure D")

ggsave("plots/special/VAMB_pCCA_antib1wk.png", width = 7, height =5)

# Antibiotics given to child 4 wk up to sampling #
y_clinical %>% group_by(ab_child_4wk_before) %>% summarise(N=n())
#We see there's 6 NAs
#making a index that will leave out the sample(s) with the NA in
index <- !is.na(y_clinical$ab_child_4wk_before)
x_crispr_antibChild4wk <- x_crispr[index,]
y_clinical_antibChild4wk <- y_clinical[index,]
z_bact_antibChild4wk <- z_bact[index,]

#running the pCCA 
CRISPR_antibChild4wk_pcca <- cca(x_crispr_antibChild4wk ~ ab_child_4wk_before + Condition(z_bact_antibChild4wk), data = y_clinical_antibChild4wk)
CRISPR_antibChild4wk_pcca
#alias(CRISPR_antibChild4wk_pcca, names=TRUE)
#VAMB: bacteria is 73.50%, antibiotics to child 4 wks before sampling is 0.17%

df <- data.frame(
  Group = c("Bacteria (73.50%)", "Child given AB 4 weeks up \nuntil sampling (0.18%)", "Unexplained (26.32%)"),
  value = c(73.50, 0.18, 26.32))

bp_antib4wk <- ggplot(df, aes(x="", y=value, fill=Group))+
  geom_bar(width = 1, stat = "identity")

bp_antib4wk + 
  coord_polar("y", start=0) +
  theme(axis.ticks = element_blank(),
        axis.title = element_blank()) +
  labs(title="Variation explained by bacteria and antibiotics given 4 weeks \nup until sampling in pCCA analysis", 
       subtitle="Species level (ANOVA test P-value: 0.204)",
       caption="Figure C")
       
ggsave("plots/special/VAMB_pCCA_antib4wk.png", width = 7, height =5)

# Antibiotics given to child 5 weeks up to sampling #
y_clinical %>% group_by(ab_child_5wk_before) %>% summarise(N=n())
#We see there's 6 NAs
#making a index that will leave out the sample(s) with the NA in
index <- !is.na(y_clinical$ab_child_5wk_before)
x_crispr_antibChild5wk <- x_crispr[index,]
y_clinical_antibChild5wk <- y_clinical[index,]
z_bact_antibChild5wk <- z_bact[index,]

#running the pCCA 
CRISPR_antibChild5wk_pcca <- cca(x_crispr_antibChild5wk ~ ab_child_5wk_before + Condition(z_bact_antibChild5wk), data = y_clinical_antibChild5wk)
CRISPR_antibChild5wk_pcca
#alias(CRISPR_antibChild4wk_pcca, names=TRUE)
#VAMB: bacteria is 73.50%, antibiotics to child 5 wks before sampling is 0.15%

df <- data.frame(
  Group = c("Bacteria (73.50%)", "Child given AB 5 weeks up \nuntil sampling (0.15%)", "Unexplained (26.35%)"),
  value = c(73.50, 0.15, 26.35))

bp_antib5wk <- ggplot(df, aes(x="", y=value, fill=Group))+
  geom_bar(width = 1, stat = "identity")

bp_antib5wk + 
  coord_polar("y", start=0) +
  theme(axis.ticks = element_blank(),
        axis.title = element_blank()) +
  labs(title="Variation explained by bacteria and antibiotics given 5 weeks \nup until sampling in pCCA analysis", 
       subtitle="Species level (ANOVA test P-value: 0.304)",
       caption="Figure B")

ggsave("plots/special/VAMB_pCCA_antib5wk.png", width = 7, height =5)


#inspecting the pcca object
#alias(CRISPR_pam_pcca, names=TRUE)
#CRISPR_pam_pcca$anova
#vif.cca(CRISPR_pam_pcca)
#pstat <- permustats(anova(CRISPR_pam_pcca))
#summary(pstat)
#densityplot(pstat)
#pam_pcca <- plot(CRISPR_pam_pcca)



### ANOVA analysis to get P-values ###

# Set seed so we don't get different P-values every time
set.seed(123)

## Supplements ##
# Fish oil #
anova(CRISPR_fish_pcca, permutations = how(nperm=999),
      by = NULL, model = c("reduced", "direct", "full"), 
      parallel = getOption("mc.cores"), strata = NULL,
      cutoff = 1, scope = NULL)
#VAMB: P-value: species: 0.901 (not significant)
#VAMB: GENUS: 0.088 (not significant)

# Vitamin D #
set.seed(123)
anova(CRISPR_vit_pcca, permutations = how(nperm=999),
      by = NULL, model = c("reduced", "direct", "full"), 
      parallel = getOption("mc.cores"), strata = NULL,
      cutoff = 1, scope = NULL)
#VAMB: P-value: 0.316 (not significant)

# Delivery #
set.seed(123)
anova(CRISPR_del_pcca, permutations = how(nperm=999),
      by = NULL, model = c("reduced", "direct", "full"), 
      parallel = getOption("mc.cores"), strata = NULL,
      cutoff = 1, scope = NULL)
#VAMB: P-value: 0.829 (not significant)

# Animal days #
set.seed(123)
anova(CRISPR_ani_pcca, permutations = how(nperm=999),
      by = NULL, model = c("reduced", "direct", "full"), 
      parallel = getOption("mc.cores"), strata = NULL,
      cutoff = 1, scope = NULL)
#VAMB: P-value: 0.496 (not significant)

## Asthma ##
# Cross #
set.seed(123)
anova(CRISPR_cross_pcca, permutations = how(nperm=999),
      by = NULL, model = c("reduced", "direct", "full"), 
      parallel = getOption("mc.cores"), strata = NULL,
      cutoff = 1, scope = NULL)
#VAMB: P-value: 0.89 (not significant)

# Ever #
set.seed(123)
anova(CRISPR_ever_pcca, permutations = how(nperm=999),
      by = NULL, model = c("reduced", "direct", "full"), 
      parallel = getOption("mc.cores"), strata = NULL,
      cutoff = 1, scope = NULL)
#VAMB: P-value: 0.732 (not significant)

# Asthmatic mother #
set.seed(123)
anova(CRISPR_astMother_pcca, permutations = how(nperm=999),
      by = NULL, model = c("reduced", "direct", "full"), 
      parallel = getOption("mc.cores"), strata = NULL,
      cutoff = 1, scope = NULL)
#VAMB: P-value: 0.866 (not significant)

# Pam clusters #
set.seed(123)
anova(CRISPR_pam_pcca, permutations = how(nperm=999),
      by = NULL, model = c("reduced", "direct", "full"), 
      parallel = getOption("mc.cores"), strata = NULL,
      cutoff = 1, scope = NULL)
#VAMB: P-value: 0.052 (almost significant)


# Rural urban #
set.seed(123)
anova(CRISPR_urb_pcca, permutations = how(nperm=999),
      by = NULL, model = c("reduced", "direct", "full"), 
      parallel = getOption("mc.cores"), strata = NULL,
      cutoff = 1, scope = NULL)
#VAMB: P-value: 0.163 (not significant)

## Antibiotics ##
# Antibiotics given to child at birth #
set.seed(123)
anova(CRISPR_antibChild_pcca, permutations = how(nperm=999),
      by = NULL, model = c("reduced", "direct", "full"), 
      parallel = getOption("mc.cores"), strata = NULL,
      cutoff = 1, scope = NULL)
#VAMB: P-value: 0.854 (not significant)

# Antibiotics given to mom during labor #
set.seed(123)
anova(CRISPR_antibMother_pcca, permutations = how(nperm=999),
      by = NULL, model = c("reduced", "direct", "full"), 
      parallel = getOption("mc.cores"), strata = NULL,
      cutoff = 1, scope = NULL)
#VAMB: P-value: 0.771 (not significant)

# Antibiotics given to child at any time up until 1 yr #
set.seed(123)
anova(CRISPR_antibChildEver_pcca, permutations = how(nperm=999),
      by = NULL, model = c("reduced", "direct", "full"), 
      parallel = getOption("mc.cores"), strata = NULL,
      cutoff = 1, scope = NULL)
#VAMB: P-value: 0.833 (not significant)

# Antibiotics given to child 1wk up until sampling #
set.seed(123)
anova(CRISPR_antibChild1wk_pcca, permutations = how(nperm=999),
      by = NULL, model = c("reduced", "direct", "full"), 
      parallel = getOption("mc.cores"), strata = NULL,
      cutoff = 1, scope = NULL)
#VAMB: P-value: 0.391 (not significant)

# Antibiotics given to child 4wk up until sampling #
set.seed(123)
anova(CRISPR_antibChild4wk_pcca, permutations = how(nperm=999),
      by = NULL, model = c("reduced", "direct", "full"), 
      parallel = getOption("mc.cores"), strata = NULL,
      cutoff = 1, scope = NULL)
#VAMB: P-value: 0.204 (not significant)

# Antibiotics given to child 5wk up until sampling #
set.seed(123)
anova(CRISPR_antibChild5wk_pcca, permutations = how(nperm=999),
      by = NULL, model = c("reduced", "direct", "full"), 
      parallel = getOption("mc.cores"), strata = NULL,
      cutoff = 1, scope = NULL)
#VAMB: P-value: 0.304 (not significant)
