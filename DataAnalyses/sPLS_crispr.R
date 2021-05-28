rm(list=ls())

library(devtools)
library(phyloseq)
library(gtools)
library(ggplot2)
library(doMC)
library(caret)
library(pROC)
#install.packages("mixOmics")
#BiocManager::install("mixOmics")
library(mixOmics)
#install_github("jonathanth/mixOmicsCaret")
library(mixOmicsCaret)
library(tidyverse)

# Load in data ####

load('crispr/CRISPRsInClusters.RData') 
load('VAMB/210415_VAMB_phyloseq_v02.RData') 


# Pre-processing of CRISPR and Bacteria matrixes ####
# Changing sample names to match what's in the bacterial phyloseq object
clusterCRISPRs <- clusterCRISPRs %>% 
  mutate(rowname = str_replace_all(rowname, "-", ".")) %>% 
  mutate(rowname = str_replace_all(rowname, "S", "X")) %>% 
  #replace NAs with 0 
  mutate_if(is.integer, ~replace(., is.na(.),0)) %>% 
  #extract samplename into column
  column_to_rownames()

# Changing cluster_XX to clusterXX
names(clusterCRISPRs) <- str_replace(names(clusterCRISPRs),"_","")

# Extract bacterial abundances
bact.abundance <- otu_table(final.physeq) %>% data.frame()

# Get same colnames & rownames 
sort.bact.abundance <- bact.abundance[ , order(names(bact.abundance))]
sort.clusterCRISPRs <- clusterCRISPRs[ , order(names(clusterCRISPRs))]
colnames(sort.clusterCRISPRs)==colnames(sort.bact.abundance)

# Remove bacteria only found in 20% of the kids
keep_index <- colMeans(bact.abundance > 0) > 0.2
sort.bact.abundance <- sort.bact.abundance[, keep_index]
sort.clusterCRISPRs <- sort.clusterCRISPRs[, keep_index]
colnames(sort.clusterCRISPRs)==colnames(sort.bact.abundance)

# Order samples to be in same order
sort.bact.abundance <- sort.bact.abundance[order(rownames(sort.bact.abundance)) , ]
sort.clusterCRISPRs <- sort.clusterCRISPRs[order(rownames(sort.clusterCRISPRs)) , ]
rownames(sort.clusterCRISPRs)==rownames(sort.bact.abundance)


# From tibble to data frame
sort.bact.abundance <- rownames_to_column(sort.bact.abundance, var = "sample")
sort.clusterCRISPRs <- rownames_to_column(sort.clusterCRISPRs, var = "sample")

# minimum as pseudocount, minus by mean and divide by standard deviation
bacteria_transf <- sort.bact.abundance %>% 
  gather(var,val,-sample) %>%
  group_by(var) %>% 
  mutate(mn = min(val[val>0]),
         val3 = log(val + mn), 
         val3 = val3 - mean(val3), 
         val3 = val3/sd(val3)) %>% 
  select(-val,-mn) %>% 
  spread(var,val3)

# Create concatenated data
X <- bacteria_transf %>% 
  gather(var,val,-sample) %>%
  left_join(sort.clusterCRISPRs %>% 
              gather(var,crpval,-sample), by = c('sample','var')
            )

# CRISPR positive
X1 <- X %>% 
  mutate(y = val*crpval) %>% 
  select(sample,var,y) %>% 
  spread(var,y)

# CRISPR negative
X0 <- X %>% 
  mutate(y = val*as.numeric(crpval==0)) %>% 
  select(sample,var,y) %>% 
  spread(var,y)



#Add outcomes from CRISPR phyloseq #### 
load('newCRISPRs.RData')
all_clinical <- get_variable(newCRISPRs) %>% data.frame
all_clinical <- all_clinical %>% 
  select(-ABCno) %>% 
  rename(sample = SampleName) %>% 
  mutate(sample = str_replace_all(sample, "-", ".")) %>% 
  mutate(sample = str_replace_all(sample, "18097D", "X18097D"))



# Run baseline (pamcluster) #### 
pamclusters_table <- all_clinical %>% 
  select(sample, pamcluster)

# Join outcome with sample data
x <- bacteria_transf %>% 
  left_join(pamclusters_table, by = 'sample') %>% 
  mutate(class = pamcluster) %>% 
  filter(!is.na(class)) %>% 
  select(-pamcluster, -sample)

registerDoMC(cores=3)

# make training and test set
set.seed(123)
inTraining <- createDataPartition(x$class, p = .75, list = FALSE)
training <- x[ inTraining,]
testing  <- x[-inTraining,]

#10-fold cross validation, 11 times
repCV10 <- trainControl(method = "repeatedcv", 
                        number = 10, 
                        repeats = 11, 
                        returnResamp = "all", 
                        savePredictions = "all", 
                        allowParallel = T, 
                        verboseIter = F)

set.seed(123)
sim_PLS <- train(as.numeric(class == "1") ~ ., data = training,
                 method = get_mixOmics_spls(),
                 preProc = c("center", "scale"),
                 tuneGrid = expand.grid(ncomp = 1:5, 
                                        keepX = c(1, 2, 4, 8, 16, 25, 50, 125, 270, 380, 411),
                                        keepY = 1),
                 trControl = repCV10,
                 fixX = c())

sim_PLS
#ncomp = 2, keepX = 125 and keepY = 1


sim_PLS$pred %>%
  separate(Resample, c("Fold", "Rep")) %>%
  group_by(ncomp, keepX, Rep) %>%
  summarize(auc = pROC::auc(obs ~ pred, direction = '<') %>% as.numeric) %>%
  ggplot(aes(x = factor(keepX), y = auc)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = factor(ncomp)), alpha = 0.2, show.legend = FALSE) +
  geom_hline(yintercept = 0.5) +
  facet_grid(. ~ ncomp) +
  scale_color_brewer(palette = "Set1", name = NULL) +
  theme_bw() + theme(strip.background = element_blank()) +
  xlab("Number of variables") +
  ggtitle("sPLS models for PAM clusters - Bacteria only")+
  labs(caption = "Figure A")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 50, hjust = 1))
ggsave("plots/special/AUC_PAM_baseline.png", width = 13, height =5)

bestSet <- sim_PLS %>% get_best_predictions %>% group_by(Rep) 
head(bestSet)

bestSet %>% summarize(auc = pROC::auc(obs, pred) %>% as.numeric )%>% arrange(auc)

savedPreds <- sim_PLS %>% get_best_predictions %>% filter(Rep == "Rep09")
head(savedPreds)

# boxplot training set
cvauc <- auc(obs ~ pred, data = savedPreds)
qplot(ifelse(obs == 1, "Mature composition", "Immature composition"), pred, data = savedPreds, geom = c("boxplot", "jitter")) + 
  xlab("Training class") +
  ylab("Training predictions") +
  ggtitle(paste0("CV AUC = ", round(cvauc, 5)))

# boxplot test set
testPreds <- predict(sim_PLS, testing)
cvauc <- auc(predictor = testPreds, testing$class)
qplot(factor(testing$class), testPreds, geom = c("boxplot", "jitter")) + 
  xlab("Testing class (PAM cluster)") +
  ylab("Testing predictions") +
  labs(title="PAM clusters - Bacteria only", 
         subtitle=paste0("Test AUC = ", round(cvauc, 5)),
         caption="Figure A")
ggsave("plots/special/sPLS_testset_PAM_baseline.png", width = 7, height =5)



# Run sPLS for pamcluster #### 

# Join outcome with crispr data
x <- X0 %>% left_join(X1,by = 'sample', suffix = c(".without", ".with")) %>% 
  left_join(pamclusters_table, by = 'sample') %>% 
  mutate(class = pamcluster) %>% 
  filter(!is.na(class)) %>% 
  select(-pamcluster, -sample)

# Remove rare crispr types
keep_index <- colMeans(x[,-ncol(x)] > 0) > 0.05
x <- cbind(x[,-ncol(x)][, keep_index], x[,ncol(x)])


registerDoMC(cores=3)

# make training and test set
set.seed(123)
inTraining <- createDataPartition(x$class, p = .75, list = FALSE)
training <- x[ inTraining,]
testing  <- x[-inTraining,]

repCV10 <- trainControl(method = "repeatedcv", 
                        number = 10, 
                        repeats = 11, 
                        returnResamp = "all", 
                        savePredictions = "all", 
                        allowParallel = T, 
                        verboseIter = F)

set.seed(123)
sim_PLS <- train(as.numeric(class == "1") ~ ., data = training,
                 method = get_mixOmics_spls(),
                 preProc = c("center", "scale"),
                 tuneGrid = expand.grid(ncomp = 1:5, 
                                        keepX = c(1, 2, 4, 8, 16, 25, 50, 125, 270, 380, 438),
                                        keepY = 1),
                 trControl = repCV10,
                 fixX = c())

sim_PLS
#ncomp = 2, keepX = 125 and keepY = 1.

sim_PLS$pred %>%
  separate(Resample, c("Fold", "Rep")) %>%
  group_by(ncomp, keepX, Rep) %>%
  summarize(auc = pROC::auc(obs ~ pred, direction = '<') %>% as.numeric) %>%
  ggplot(aes(x = factor(keepX), y = auc)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = factor(ncomp)), alpha = 0.2, show.legend = FALSE) +
  geom_hline(yintercept = 0.5) +
  facet_grid(. ~ ncomp) +
  scale_color_brewer(palette = "Set1", name = NULL) +
  theme_bw() + theme(strip.background = element_blank()) +
  xlab("Number of variables") +
  ggtitle("sPLS models for PAM clusters - Concatenated table") + 
  labs(caption = "Figure B")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 50, hjust = 1))
ggsave("plots/special/AUC_PAM.png", width = 13, height =5)

bestSet <- sim_PLS %>% get_best_predictions %>% group_by(Rep) 
head(bestSet)

bestSet %>% summarize(auc = pROC::auc(obs, pred) %>% as.numeric )%>% arrange(auc)

savedPreds <- sim_PLS %>% get_best_predictions %>% filter(Rep == "Rep01")
head(savedPreds)

# boxplot training set
cvauc <- auc(obs ~ pred, data = savedPreds)
qplot(ifelse(obs == 1, "Present", "Not present"), pred, data = savedPreds, geom = c("boxplot", "jitter")) + 
  xlab("Training class") +
  ylab("Training predictions") +
  ggtitle(paste0("CV AUC = ", round(cvauc, 5)))

# boxplot test set
testPreds <- predict(sim_PLS, testing)
cvauc <- auc(predictor = testPreds, testing$class)
qplot(factor(testing$class), testPreds, geom = c("boxplot", "jitter")) + 
  xlab("Testing class (PAM cluster)") +
  ylab("Testing predictions") +
  labs(title="PAM clusters - Concatenated table", 
       subtitle=paste0("Test AUC = ", round(cvauc, 5)),
       caption="Figure B")

ggsave("plots/special/sPLS_testset_PAM.png", width = 7, height =5)

# do log regression on the training set and then the test set
glm(obs ~ pred, data = savedPreds, family = binomial) %>% summary
glm(testing$class ~ testPreds, family = binomial) %>% summary

# get test variables
sim_PLS %>% get_loadings %>% head # Returns a data.frame for easy plotting

sim_PLS %>% get_loadings("CV", rep = "Rep01") %>% 
  ggplot(aes(var, loading, ymin = loading - sd, ymax = loading + sd)) + 
  facet_wrap(~ comp, scales = "free_x") + 
  geom_errorbar() + 
  geom_bar(stat = "identity") + 
  ggtitle("CV, median rep (Rep1)")

sim_PLS %>% get_loadings("finalModel") %>% 
  ggplot(aes(var, loading, ymin = loading - sd, ymax = loading + sd)) + 
  facet_wrap(~ comp, scales = "free_x") + 
  geom_errorbar() + 
  geom_bar(stat = "identity") + 
  ggtitle("Full model") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 50, hjust = 1))


# Run baseline (AB) #### 

#making table for clinical variable only
AB4wk_table <- all_clinical %>% 
  select(sample, ab_child_4wk_before)

# Join outcome with sample data
x <- bacteria_transf %>% 
  left_join(AB4wk_table, by = 'sample') %>% 
  mutate(class = ab_child_4wk_before) %>% 
  filter(!is.na(class)) %>% 
  select(-ab_child_4wk_before, -sample)

registerDoMC(cores=3)

# make training and test set
set.seed(123)
inTraining <- createDataPartition(x$class, p = .75, list = FALSE)
training <- x[ inTraining,]
testing  <- x[-inTraining,]

repCV10 <- trainControl(method = "repeatedcv", 
                        number = 10, 
                        repeats = 11, 
                        returnResamp = "all", 
                        savePredictions = "all", 
                        allowParallel = T, 
                        verboseIter = F)

set.seed(123)
sim_PLS <- train(as.numeric(class == "1") ~ ., data = training,
                 method = get_mixOmics_spls(),
                 preProc = c("center", "scale"),
                 tuneGrid = expand.grid(ncomp = 1:5, 
                                        keepX = c(1, 2, 4, 8, 16, 25, 50, 125, 270, 380, 411),
                                        keepY = 1),
                 trControl = repCV10,
                 fixX = c())

sim_PLS
#ncomp = 1, keepX = 8 and keepY = 1

sim_PLS$pred %>%
  separate(Resample, c("Fold", "Rep")) %>%
  group_by(ncomp, keepX, Rep) %>%
  summarize(auc = pROC::auc(obs ~ pred, direction = '<') %>% as.numeric) %>%
  ggplot(aes(x = factor(keepX), y = auc)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = factor(ncomp)), alpha = 0.2, show.legend = FALSE) +
  geom_hline(yintercept = 0.5) +
  facet_grid(. ~ ncomp) +
  scale_color_brewer(palette = "Set1", name = NULL) +
  theme_bw() + theme(strip.background = element_blank()) +
  xlab("Number of variables") +
  ggtitle("sPLS models for antibiotics given 4 weeks before sampling - Bacteria only")+
  labs(caption = "Figure A")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 50, hjust = 1))
ggsave("plots/special/AUC_AB4wks_baseline.png", width = 13, height =5)

bestSet <- sim_PLS %>% get_best_predictions %>% group_by(Rep) 
head(bestSet)

bestSet %>% summarize(auc = pROC::auc(obs, pred) %>% as.numeric )%>% arrange(auc)

savedPreds <- sim_PLS %>% get_best_predictions %>% filter(Rep == "Rep03")
head(savedPreds)

# boxplot training set
cvauc <- auc(obs ~ pred, data = savedPreds)
qplot(ifelse(obs == 1, "Present", "Not present"), pred, data = savedPreds, geom = c("boxplot", "jitter")) + 
  xlab("Training class") +
  ylab("Training predictions") +
  ggtitle(paste0("CV AUC = ", round(cvauc, 5)))

# boxplot test set
testPreds <- predict(sim_PLS, testing)
cvauc <- auc(predictor = testPreds, testing$class)
qplot(factor(testing$class), testPreds, geom = c("boxplot", "jitter")) + 
  xlab("Testing class (AB 4 weeks before sampling)") +
  ylab("Testing predictions") +
  labs(title="Antibiotics 4 weeks before sampling - Bacteria only", 
       subtitle=paste0("Test AUC = ", round(cvauc, 5)),
       caption="Figure A")
ggsave("plots/special/sPLS_testset_AB4wks_baseline.png", width = 7, height =5)


# Run sPLS for antibiotics ####

# Join outcome with crispr data
x <- X0 %>% left_join(X1,by = 'sample', suffix = c(".without", ".with")) %>% 
  left_join(AB4wk_table, by = 'sample') %>% 
  mutate(class = ab_child_4wk_before) %>% 
  filter(!is.na(class)) %>% 
  select(-ab_child_4wk_before, -sample)

# Remove rare crispr types
keep_index <- colMeans(x[,-ncol(x)] > 0) > 0.05
x <- cbind(x[,-ncol(x)][, keep_index], x[,ncol(x)])


registerDoMC(cores=3)

# make training and test set
set.seed(123)
inTraining <- createDataPartition(x$class, p = .75, list = FALSE)
training <- x[ inTraining,]
testing  <- x[-inTraining,]

repCV10 <- trainControl(method = "repeatedcv", 
                        number = 10, 
                        repeats = 11, 
                        returnResamp = "all", 
                        savePredictions = "all", 
                        allowParallel = T, 
                        verboseIter = F)

set.seed(123)
sim_PLS <- train(as.numeric(class == "1") ~ ., data = training,
                 method = get_mixOmics_spls(),
                 preProc = c("center", "scale"),
                 tuneGrid = expand.grid(ncomp = 1:5, 
                                        keepX = c(1, 2, 4, 8, 16, 25, 50, 125, 270, 350, 437),
                                        keepY = 1),
                 trControl = repCV10,
                 fixX = c())

sim_PLS
#ncomp = 1, keepX = 8 and keepY = 1

ggplot(sim_PLS, metric = "RMSE")
ggplot(sim_PLS, metric = "RMSE") + scale_x_log10()
ggplot(sim_PLS, metric = "Rsquared")
ggplot(sim_PLS, metric = "Rsquared") + scale_x_log10()

sim_PLS$pred %>%
  separate(Resample, c("Fold", "Rep")) %>%
  group_by(ncomp, keepX, Rep) %>%
  summarize(auc = pROC::auc(obs ~ pred, direction = '<') %>% as.numeric) %>%
  ggplot(aes(x = factor(keepX), y = auc)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = factor(ncomp)), alpha = 0.2, show.legend = FALSE) +
  geom_hline(yintercept = 0.5) +
  facet_grid(. ~ ncomp) +
  scale_color_brewer(palette = "Set1", name = NULL) +
  theme_bw() + theme(strip.background = element_blank()) +
  xlab("Number of variables") +
  ggtitle("sPLS models for antibiotics given 4 weeks before sampling - Concatenated table") + 
  labs(caption = "Figure B")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 50, hjust = 1))
ggsave("plots/special/AUC_AB4wks.png", width = 13, height =5)


bestSet <- sim_PLS %>% get_best_predictions %>% group_by(Rep) 
head(bestSet)

bestSet %>% summarize(auc = pROC::auc(obs, pred) %>% as.numeric )%>% arrange(auc)

savedPreds <- sim_PLS %>% get_best_predictions %>% filter(Rep == "Rep03")
head(savedPreds)

# boxplot training set
cvauc <- auc(obs ~ pred, data = savedPreds)
qplot(ifelse(obs == 1, "Present", "Not present"), pred, data = savedPreds, geom = c("boxplot", "jitter")) + 
  xlab("Training class") +
  ylab("Training predictions") +
  ggtitle(paste0("CV AUC = ", round(cvauc, 5)))

# boxplot test set
testPreds <- predict(sim_PLS, testing)
cvauc <- auc(predictor = testPreds, testing$class)
qplot(factor(testing$class), testPreds, geom = c("boxplot", "jitter")) + 
  xlab("Testing class (AB 4 weeks before sampling)") +
  ylab("Testing predictions") +
  labs(title="Antibiotics 4 weeks before sampling - Concatenated table", 
       subtitle=paste0("Test AUC = ", round(cvauc, 5)),
       caption="Figure B")

ggsave("plots/special/sPLS_testset_AB4wks.png", width = 7, height =5)


# do log regression on the training set and then the test set
glm(obs ~ pred, data = savedPreds, family = binomial) %>% summary
glm(testing$class ~ testPreds, family = binomial) %>% summary

# get test variables
sim_PLS %>% get_loadings %>% head # Returns a data.frame for easy plotting

sim_PLS %>% get_loadings("CV", rep = "Rep03") %>% 
  ggplot(aes(var, loading, ymin = loading - sd, ymax = loading + sd)) + 
  facet_wrap(~ comp, scales = "free_x") + 
  geom_errorbar() + 
  geom_bar(stat = "identity") + 
  ggtitle("CV, median rep (Rep3)")
ggsave("plots/special/sPLS_loadings_AB4wks.png", width = 10, height =4)

sim_PLS %>% get_loadings("finalModel") %>% 
  ggplot(aes(var, loading, ymin = loading - sd, ymax = loading + sd)) + 
  facet_wrap(~ comp, scales = "free_x") + 
  geom_errorbar() + 
  geom_bar(stat = "identity") + 
  ggtitle("Full model") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1))



# Run baseline (Asthma) #### 

#making table for clinical variable only
asthmaCross_table <- all_clinical %>% 
  select(sample, j45_5yr_cross)

# Join outcome with sample data
x <- bacteria_transf %>% 
  left_join(asthmaCross_table, by = 'sample') %>% 
  mutate(class = j45_5yr_cross) %>% 
  filter(!is.na(class)) %>% 
  select(-j45_5yr_cross, -sample)

registerDoMC(cores=3)

# make training and test set
set.seed(123)
inTraining <- createDataPartition(x$class, p = .75, list = FALSE)
training <- x[ inTraining,]
testing  <- x[-inTraining,]

repCV10 <- trainControl(method = "repeatedcv", 
                        number = 10, 
                        repeats = 11, 
                        returnResamp = "all", 
                        savePredictions = "all", 
                        allowParallel = T, 
                        verboseIter = F)

set.seed(123)
sim_PLS <- train(as.numeric(class == "1") ~ ., data = training,
                 method = get_mixOmics_spls(),
                 preProc = c("center", "scale"),
                 tuneGrid = expand.grid(ncomp = 1:5, 
                                        keepX = c(1, 2, 4, 8, 16, 25, 50, 125, 270, 380, 411),
                                        keepY = 1),
                 trControl = repCV10,
                 fixX = c())

sim_PLS
#ncomp = 1, keepX = 2 and keepY = 1

sim_PLS$pred %>%
  separate(Resample, c("Fold", "Rep")) %>%
  group_by(ncomp, keepX, Rep) %>%
  summarize(auc = pROC::auc(obs ~ pred, direction = '<') %>% as.numeric) %>%
  ggplot(aes(x = factor(keepX), y = auc)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = factor(ncomp)), alpha = 0.2, show.legend = FALSE) +
  geom_hline(yintercept = 0.5) +
  facet_grid(. ~ ncomp) +
  scale_color_brewer(palette = "Set1", name = NULL) +
  theme_bw() + theme(strip.background = element_blank()) +
  xlab("Number of variables") +
  ggtitle("sPLS models for asthma diagnosed at age 5 - Bacteria only")+
  labs(caption = "Figure A")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 50, hjust = 1))

ggsave("plots/special/AUC_asthmaCross_baseline.png", width = 13, height =5)

bestSet <- sim_PLS %>% get_best_predictions %>% group_by(Rep) 
head(bestSet)

bestSet %>% summarize(auc = pROC::auc(obs, pred, direction = '<') %>% as.numeric )%>% arrange(auc)

savedPreds <- sim_PLS %>% get_best_predictions %>% filter(Rep == "Rep07")
head(savedPreds)

# boxplot training set
cvauc <- auc(obs ~ pred, data = savedPreds, direction = '<')
qplot(ifelse(obs == 1, "Present", "Not present"), pred, data = savedPreds, geom = c("boxplot", "jitter")) + 
  xlab("Training class") +
  ylab("Training predictions") +
  ggtitle(paste0("CV AUC = ", round(cvauc, 5)))

# boxplot test set
testPreds <- predict(sim_PLS, testing)
cvauc <- auc(predictor = testPreds, testing$class)
qplot(factor(testing$class), testPreds, geom = c("boxplot", "jitter")) + 
  xlab("Testing class (PAM cluster)") +
  ylab("Testing predictions") +
  labs(title="Asthma diagnosed at age 5 - Bacteria only", 
       subtitle=paste0("Test AUC = ", round(cvauc, 5)),
       caption="Figure A")
ggsave("plots/special/sPLS_testset_AsthmaCross_baseline.png", width = 7, height =5)


# Run sPLS for asthma ####

# Join outcome with crispr data
x <- X0 %>% left_join(X1,by = 'sample', suffix = c(".without", ".with")) %>% 
  left_join(asthmaCross_table, by = 'sample') %>% 
  mutate(class = j45_5yr_cross) %>% 
  filter(!is.na(class)) %>% 
  select(-j45_5yr_cross, -sample)

# Remove rare crispr types
keep_index <- colMeans(x[,-ncol(x)] > 0) > 0.05
x <- cbind(x[,-ncol(x)][, keep_index], x[,ncol(x)])


registerDoMC(cores=3)

# make training and test set
set.seed(123)
inTraining <- createDataPartition(x$class, p = .75, list = FALSE)
training <- x[ inTraining,]
testing  <- x[-inTraining,]

repCV10 <- trainControl(method = "repeatedcv", 
                        number = 10, 
                        repeats = 11, 
                        returnResamp = "all", 
                        savePredictions = "all", 
                        allowParallel = T, 
                        verboseIter = F)

set.seed(123)
sim_PLS <- train(as.numeric(class == "1") ~ ., data = training,
                 method = get_mixOmics_spls(),
                 preProc = c("center", "scale"),
                 tuneGrid = expand.grid(ncomp = 1:5, 
                                        keepX = c(1, 2, 4, 8, 16, 25, 50, 125, 270, 350, 438),
                                        keepY = 1),
                 trControl = repCV10,
                 fixX = c())

sim_PLS
#ncomp = 1, keepX = 2 and keepY = 1

ggplot(sim_PLS, metric = "RMSE")
ggplot(sim_PLS, metric = "RMSE") + scale_x_log10()
ggplot(sim_PLS, metric = "Rsquared")
ggplot(sim_PLS, metric = "Rsquared") + scale_x_log10()

sim_PLS$pred %>%
  separate(Resample, c("Fold", "Rep")) %>%
  group_by(ncomp, keepX, Rep) %>%
  summarize(auc = pROC::auc(obs ~ pred, direction = '<') %>% as.numeric) %>%
  ggplot(aes(x = factor(keepX), y = auc)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = factor(ncomp)), alpha = 0.2, show.legend = FALSE) +
  geom_hline(yintercept = 0.5) +
  facet_grid(. ~ ncomp) +
  scale_color_brewer(palette = "Set1", name = NULL) +
  theme_bw() + theme(strip.background = element_blank()) +
  xlab("Number of variables") +
  ggtitle("sPLS models for asthma diagnosed at age 5 - Concatenated table") + 
  labs(caption = "Figure B")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 50, hjust = 1))

ggsave("plots/special/AUC_AsthmaCross_crispr.png", width = 13, height =5)

bestSet <- sim_PLS %>% get_best_predictions %>% group_by(Rep) 
head(bestSet)

bestSet %>% summarize(auc = pROC::auc(obs, pred, direction = '<') %>% as.numeric )%>% arrange(auc)

savedPreds <- sim_PLS %>% get_best_predictions %>% filter(Rep == "Rep07")
head(savedPreds)

# boxplot training set
cvauc <- auc(obs ~ pred, data = savedPreds)
qplot(ifelse(obs == 1, "Present", "Not present"), pred, data = savedPreds, geom = c("boxplot", "jitter")) + 
  xlab("Training class") +
  ylab("Training predictions") +
  ggtitle(paste0("CV AUC = ", round(cvauc, 5)))

# boxplot test set
testPreds <- predict(sim_PLS, testing)
cvauc <- auc(predictor = testPreds, testing$class)
qplot(factor(testing$class), testPreds, geom = c("boxplot", "jitter")) + 
  xlab("Testing class (Asthma)") +
  ylab("Testing predictions") +
  labs(title="Asthma diagnosed at age 5 - Concatenated table", 
       subtitle=paste0("Test AUC = ", round(cvauc, 5)),
       caption="Figure B")
ggsave("plots/special/sPLS_testset_AsthmaCross.png", width = 7, height =5)

qplot(factor(testing$class), testPreds, geom = c("boxplot", "jitter")) + 
  xlab("Testing class (Asthma)") +
  ylab("Testing predictions") +
  labs(title="Asthma diagnosed at age 5 - Concatenated table", 
       subtitle=paste0("Test AUC = ", round(cvauc, 5)))+
  ylim(c(0.05,0.15))
ggsave("plots/special/sPLS_testset_AsthmaCrossZOOM.png", width = 7, height =5)


# do log regression on the training set and then the test set
glm(obs ~ pred, data = savedPreds, family = binomial) %>% summary
glm(testing$class ~ testPreds, family = binomial) %>% summary

# get test variables
sim_PLS %>% get_loadings %>% head # Returns a data.frame for easy plotting

sim_PLS %>% get_loadings("CV", rep = "Rep07") %>% 
  ggplot(aes(var, loading, ymin = loading - sd, ymax = loading + sd)) + 
  facet_wrap(~ comp, scales = "free_x") + 
  geom_errorbar() + 
  geom_bar(stat = "identity") + 
  ggtitle("CV, median repetition (Rep7)")
ggsave("plots/special/sPLS_loadings_asthma.png", width = 5, height =4)

sim_PLS %>% get_loadings("finalModel") %>% 
  ggplot(aes(var, loading, ymin = loading - sd, ymax = loading + sd)) + 
  facet_wrap(~ comp, scales = "free_x") + 
  geom_errorbar() + 
  geom_bar(stat = "identity") + 
  ggtitle("Full model - Prediction of asthma at age 5") +
  theme(axis.text.x = element_text(angle = 50, hjust = 1))
ggsave("plots/special/sPLS_loadings_asthma_full.png", width = 5, height =4)

# Run sPLS for asthma (randomized sanity check) ####
# Join outcome with crispr data

x <- X0 %>% 
  left_join(X1,by = 'sample', suffix = c(".without", ".with")) %>% 
  left_join(asthmaCross_table, by = 'sample') %>% 
  mutate(class = j45_5yr_cross) %>% 
  filter(!is.na(class)) %>% 
  select(-j45_5yr_cross, -sample)

# Shuffle data around to break any significant relations
set.seed(123)
x$class <- sample(x$class, replace=FALSE)



# Remove rare crispr types
keep_index <- colMeans(x[,-ncol(x)] > 0) > 0.05
x <- cbind(x[,-ncol(x)][, keep_index], x[,ncol(x)])


registerDoMC(cores=3)

# make training and test set
set.seed(123)
inTraining <- createDataPartition(x$class, p = .75, list = FALSE)
training <- x[ inTraining,]
testing  <- x[-inTraining,]

repCV10 <- trainControl(method = "repeatedcv", 
                        number = 10, 
                        repeats = 11, 
                        returnResamp = "all", 
                        savePredictions = "all", 
                        allowParallel = T, 
                        verboseIter = F)

set.seed(123)
sim_PLS <- train(as.numeric(class == "1") ~ ., data = training,
                 method = get_mixOmics_spls(),
                 preProc = c("center", "scale"),
                 tuneGrid = expand.grid(ncomp = 1:5, 
                                        keepX = c(1, 2, 4, 8, 16, 25, 50, 125, 270, 350, 438),
                                        keepY = 1),
                 trControl = repCV10,
                 fixX = c())

sim_PLS
#ncomp = 1, keepX = 2 and keepY = 1

sim_PLS$pred %>%
  separate(Resample, c("Fold", "Rep")) %>%
  group_by(ncomp, keepX, Rep) %>%
  summarize(auc = pROC::auc(obs ~ pred, direction = '<') %>% as.numeric) %>%
  ggplot(aes(x = factor(keepX), y = auc)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = factor(ncomp)), alpha = 0.2, show.legend = FALSE) +
  geom_hline(yintercept = 0.5) +
  facet_grid(. ~ ncomp) +
  scale_color_brewer(palette = "Set1", name = NULL) +
  theme_bw() + theme(strip.background = element_blank()) +
  xlab("Number of variables") +
  ggtitle("sPLS models for asthma diagnosed at age 5 - Concatenated table") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 50, hjust = 1))

#ggsave("plots/special/AUC_AsthmaCross_crispr.png", width = 13, height =5)

bestSet <- sim_PLS %>% get_best_predictions %>% group_by(Rep) 
head(bestSet)

bestSet %>% summarize(auc = pROC::auc(obs, pred, direction = '<') %>% as.numeric )%>% arrange(auc)

savedPreds <- sim_PLS %>% get_best_predictions %>% filter(Rep == "Rep06")
head(savedPreds)

# boxplot training set
cvauc <- auc(obs ~ pred, data = savedPreds)
qplot(ifelse(obs == 1, "Present", "Not present"), pred, data = savedPreds, geom = c("boxplot", "jitter")) + 
  xlab("Training class") +
  ylab("Training predictions") +
  ggtitle(paste0("CV AUC = ", round(cvauc, 5)))

# boxplot test set
testPreds <- predict(sim_PLS, testing)
cvauc <- auc(predictor = testPreds, testing$class)
qplot(factor(testing$class), testPreds, geom = c("boxplot", "jitter")) + 
  xlab("Testing class (Asthma)") +
  ylab("Testing predictions") +
  labs(title="Asthma diagnosed at age 5 (randomized) - Concatenated table", 
       subtitle=paste0("Test AUC = ", round(cvauc, 5)))
ggsave("plots/special/sPLS_testset_AsthmaCross_random.png", width = 7, height =5)

# Run baseline (delivery) #### 

#making table for clinical variable only
delivery_table <- all_clinical %>% 
  select(sample, Delivery)

# Join outcome with sample data
x <- bacteria_transf %>% 
  left_join(delivery_table, by = 'sample') %>% 
  mutate(class = Delivery) %>% 
  filter(!is.na(class)) %>% 
  select(-Delivery, -sample)

registerDoMC(cores=3)

# make training and test set
set.seed(123)
inTraining <- createDataPartition(x$class, p = .75, list = FALSE)
training <- x[ inTraining,]
testing  <- x[-inTraining,]

repCV10 <- trainControl(method = "repeatedcv", 
                        number = 10, 
                        repeats = 11, 
                        returnResamp = "all", 
                        savePredictions = "all", 
                        allowParallel = T, 
                        verboseIter = F)

set.seed(123)
sim_PLS <- train(as.numeric(class == "1") ~ ., data = training,
                 method = get_mixOmics_spls(),
                 preProc = c("center", "scale"),
                 tuneGrid = expand.grid(ncomp = 1:5, 
                                        keepX = c(1, 2, 4, 8, 16, 25, 50, 125, 270, 380, 411),
                                        keepY = 1),
                 trControl = repCV10,
                 fixX = c())

sim_PLS
#ncomp = 1, keepX = 25 and keepY = 1

sim_PLS$pred %>%
  separate(Resample, c("Fold", "Rep")) %>%
  group_by(ncomp, keepX, Rep) %>%
  summarize(auc = pROC::auc(obs ~ pred, direction = '<') %>% as.numeric) %>%
  ggplot(aes(x = factor(keepX), y = auc)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = factor(ncomp)), alpha = 0.2, show.legend = FALSE) +
  geom_hline(yintercept = 0.5) +
  facet_grid(. ~ ncomp) +
  scale_color_brewer(palette = "Set1", name = NULL) +
  theme_bw() + theme(strip.background = element_blank()) +
  xlab("Number of variables") +
  ggtitle("Delivery method - Bacteria only") + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave("plots/special/AUC_delivery_baseline.pdf", width = 13, height =5)

bestSet <- sim_PLS %>% get_best_predictions %>% group_by(Rep) 
head(bestSet)

bestSet %>% summarize(auc = pROC::auc(obs, pred, direction = '<') %>% as.numeric )%>% arrange(auc)

savedPreds <- sim_PLS %>% get_best_predictions %>% filter(Rep == "Rep04")
head(savedPreds)

# boxplot training set
cvauc <- auc(obs ~ pred, data = savedPreds, direction = '<')
qplot(ifelse(obs == 1, "Present", "Not present"), pred, data = savedPreds, geom = c("boxplot", "jitter")) + 
  xlab("Training class") +
  ylab("Training predictions") +
  ggtitle(paste0("CV AUC = ", round(cvauc, 5)))

# boxplot test set
testPreds <- predict(sim_PLS, testing)
cvauc <- auc(predictor = testPreds, testing$class, direction = '<')
qplot(factor(testing$class), testPreds, geom = c("boxplot", "jitter")) + 
  xlab("Testing class") +
  ylab("Testing predictions") +
  ggtitle(paste0("Test AUC = ", round(cvauc, 5)))
ggsave("plots/special/sPLS_testset_delivery_baseline.pdf", width = 7, height =5)


# Run sPLS for delivery ####

# Join outcome with crispr data
x <- X0 %>% left_join(X1,by = 'sample', suffix = c(".without", ".with")) %>% 
  left_join(delivery_table, by = 'sample') %>% 
  mutate(class = Delivery) %>% 
  filter(!is.na(class)) %>% 
  select(-Delivery, -sample)

# Remove rare crispr types
keep_index <- colMeans(x[,-ncol(x)] > 0) > 0.05
x <- cbind(x[,-ncol(x)][, keep_index], x[,ncol(x)])


registerDoMC(cores=3)

# make training and test set
set.seed(123)
inTraining <- createDataPartition(x$class, p = .75, list = FALSE)
training <- x[ inTraining,]
testing  <- x[-inTraining,]

repCV10 <- trainControl(method = "repeatedcv", 
                        number = 10, 
                        repeats = 11, 
                        returnResamp = "all", 
                        savePredictions = "all", 
                        allowParallel = T, 
                        verboseIter = F)

set.seed(123)
sim_PLS <- train(as.numeric(class == "1") ~ ., data = training,
                 method = get_mixOmics_spls(),
                 preProc = c("center", "scale"),
                 tuneGrid = expand.grid(ncomp = 1:5, 
                                        keepX = c(1, 2, 4, 8, 16, 25, 50, 125, 270, 350, 437),
                                        keepY = 1),
                 trControl = repCV10,
                 fixX = c())

sim_PLS
#ncomp = 1, keepX = 50 and keepY = 1

ggplot(sim_PLS, metric = "RMSE")
ggplot(sim_PLS, metric = "RMSE") + scale_x_log10()
ggplot(sim_PLS, metric = "Rsquared")
ggplot(sim_PLS, metric = "Rsquared") + scale_x_log10()

sim_PLS$pred %>%
  separate(Resample, c("Fold", "Rep")) %>%
  group_by(ncomp, keepX, Rep) %>%
  summarize(auc = pROC::auc(obs ~ pred, direction = '<') %>% as.numeric) %>%
  ggplot(aes(x = factor(keepX), y = auc)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = factor(ncomp)), alpha = 0.2, show.legend = FALSE) +
  geom_hline(yintercept = 0.5) +
  facet_grid(. ~ ncomp) +
  scale_color_brewer(palette = "Set1", name = NULL) +
  theme_bw() + theme(strip.background = element_blank()) +
  xlab("Number of variables") +
  ggtitle("Delivery method") + 
  theme(plot.title = element_text(hjust = 0.5))

ggsave("plots/special/AUC_delivery_crispr.pdf", width = 13, height =5)

bestSet <- sim_PLS %>% get_best_predictions %>% group_by(Rep) 
head(bestSet)

bestSet %>% summarize(auc = pROC::auc(obs, pred) %>% as.numeric )%>% arrange(auc)

savedPreds <- sim_PLS %>% get_best_predictions %>% filter(Rep == "Rep04")
head(savedPreds)

# boxplot training set
cvauc <- auc(obs ~ pred, data = savedPreds)
qplot(ifelse(obs == 1, "Present", "Not present"), pred, data = savedPreds, geom = c("boxplot", "jitter")) + 
  xlab("Training class") +
  ylab("Training predictions") +
  ggtitle(paste0("CV AUC = ", round(cvauc, 5)))

# boxplot test set
testPreds <- predict(sim_PLS, testing)
cvauc <- auc(predictor = testPreds, testing$class, direction ='<')
qplot(factor(testing$class), testPreds, geom = c("boxplot", "jitter")) + 
  xlab("Testing class") +
  ylab("Testing predictions") +
  ggtitle(paste0("Test AUC = ", round(cvauc, 5)))
ggsave("plots/special/sPLS_testset_delivery.pdf", width = 7, height =5)


# Run baseline (fish oil) #### 

#making table for clinical variable only
fishoil_table <- all_clinical %>% 
  select(sample, Supplement_oil)

# Join outcome with sample data
x <- bacteria_transf %>% 
  left_join(fishoil_table, by = 'sample') %>% 
  mutate(class = Supplement_oil) %>% 
  filter(!is.na(class)) %>% 
  select(-Supplement_oil, -sample)

registerDoMC(cores=3)

# make training and test set
set.seed(123)
inTraining <- createDataPartition(x$class, p = .75, list = FALSE)
training <- x[ inTraining,]
testing  <- x[-inTraining,]

repCV10 <- trainControl(method = "repeatedcv", 
                        number = 10, 
                        repeats = 11, 
                        returnResamp = "all", 
                        savePredictions = "all", 
                        allowParallel = T, 
                        verboseIter = F)

set.seed(123)
sim_PLS <- train(as.numeric(class == "1") ~ ., data = training,
                 method = get_mixOmics_spls(),
                 preProc = c("center", "scale"),
                 tuneGrid = expand.grid(ncomp = 1:5, 
                                        keepX = c(1, 2, 4, 8, 16, 25, 50, 125, 270, 380, 411),
                                        keepY = 1),
                 trControl = repCV10,
                 fixX = c())

sim_PLS
#ncomp = 1, keepX = 1 and keepY = 1

sim_PLS$pred %>%
  separate(Resample, c("Fold", "Rep")) %>%
  group_by(ncomp, keepX, Rep) %>%
  summarize(auc = pROC::auc(obs ~ pred, direction = '<') %>% as.numeric) %>%
  ggplot(aes(x = factor(keepX), y = auc)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = factor(ncomp)), alpha = 0.2, show.legend = FALSE) +
  geom_hline(yintercept = 0.5) +
  facet_grid(. ~ ncomp) +
  scale_color_brewer(palette = "Set1", name = NULL) +
  theme_bw() + theme(strip.background = element_blank()) +
  xlab("Number of variables") +
  ggtitle("sPLS models for fish oil - Bacteria only") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 50, hjust = 1))
ggsave("plots/special/AUC_fishoil_baseline.png", width = 13, height =5)

bestSet <- sim_PLS %>% get_best_predictions %>% group_by(Rep) 
head(bestSet)

bestSet %>% summarize(auc = pROC::auc(obs, pred, direction = '<') %>% as.numeric )%>% arrange(auc)

savedPreds <- sim_PLS %>% get_best_predictions %>% filter(Rep == "Rep07")
head(savedPreds)

# boxplot training set
cvauc <- auc(obs ~ pred, data = savedPreds, direction = '<')
qplot(ifelse(obs == 1, "Fish oil given", "No fish oil given"), pred, data = savedPreds, geom = c("boxplot", "jitter")) + 
  xlab("Training class") +
  ylab("Training predictions") +
  labs(title="Fish oil - bacteria only", subtitle=paste0("CV AUC = ", round(cvauc, 5)))

# boxplot test set
testPreds <- predict(sim_PLS, testing)
cvauc <- auc(predictor = testPreds, testing$class, direction = '<')
qplot(factor(testing$class), testPreds, geom = c("boxplot", "jitter")) + 
  xlab("Testing class (Mother given fish oil or not)") +
  ylab("Testing predictions") +
  labs(title="Fish oil - Bacteria only", 
       subtitle=paste0("Test AUC = ", round(cvauc, 5)),
       caption="Figure A")
ggsave("plots/special/sPLS_testset_fishoil_baseline.png", width = 7, height =5)


# Run sPLS for fish oil ####

# Join outcome with crispr data
x <- X0 %>% left_join(X1,by = 'sample', suffix = c(".without", ".with")) %>% 
  left_join(fishoil_table, by = 'sample') %>% 
  mutate(class = Supplement_oil) %>% 
  filter(!is.na(class)) %>% 
  select(-Supplement_oil, -sample)

# Remove rare crispr types
keep_index <- colMeans(x[,-ncol(x)] > 0) > 0.05
x <- cbind(x[,-ncol(x)][, keep_index], x[,ncol(x)])


registerDoMC(cores=3)

# make training and test set
set.seed(123)
inTraining <- createDataPartition(x$class, p = .75, list = FALSE)
training <- x[ inTraining,]
testing  <- x[-inTraining,]

repCV10 <- trainControl(method = "repeatedcv", 
                        number = 10, 
                        repeats = 11, 
                        returnResamp = "all", 
                        savePredictions = "all", 
                        allowParallel = T, 
                        verboseIter = F)

set.seed(123)
sim_PLS <- train(as.numeric(class == "1") ~ ., data = training,
                 method = get_mixOmics_spls(),
                 preProc = c("center", "scale"),
                 tuneGrid = expand.grid(ncomp = 1:5, 
                                        keepX = c(1, 2, 4, 8, 16, 25, 50, 125, 270, 350, 437),
                                        keepY = 1),
                 trControl = repCV10,
                 fixX = c())

sim_PLS
#ncomp = 1, keepX = 1 and keepY = 1

sim_PLS$pred %>%
  separate(Resample, c("Fold", "Rep")) %>%
  group_by(ncomp, keepX, Rep) %>%
  summarize(auc = pROC::auc(obs ~ pred, direction = '<') %>% as.numeric) %>%
  ggplot(aes(x = factor(keepX), y = auc)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = factor(ncomp)), alpha = 0.2, show.legend = FALSE) +
  geom_hline(yintercept = 0.5) +
  facet_grid(. ~ ncomp) +
  scale_color_brewer(palette = "Set1", name = NULL) +
  theme_bw() + theme(strip.background = element_blank()) +
  xlab("Number of variables") +
  ggtitle("sPLS models for fish oil") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 50, hjust = 1))

ggsave("plots/special/AUC_fishoil_crispr.png", width = 13, height =5)

bestSet <- sim_PLS %>% get_best_predictions %>% group_by(Rep) 
head(bestSet)

bestSet %>% summarize(auc = pROC::auc(obs, pred, direction = '<') %>% as.numeric )%>% arrange(auc)

savedPreds <- sim_PLS %>% get_best_predictions %>% filter(Rep == "Rep05")
head(savedPreds)

# boxplot training set
cvauc <- auc(obs ~ pred, data = savedPreds)
qplot(ifelse(obs == 1, "Present", "Not present"), pred, data = savedPreds, geom = c("boxplot", "jitter")) + 
  xlab("Training class") +
  ylab("Training predictions") +
  ggtitle(paste0("CV AUC = ", round(cvauc, 5)))

# boxplot test set
testPreds <- predict(sim_PLS, testing)
cvauc <- auc(predictor = testPreds, testing$class, direction ='<')
qplot(factor(testing$class), testPreds, geom = c("boxplot", "jitter")) + 
  xlab("Testing class (Mother given fish oil or not)") +
  ylab("Testing predictions") +
  labs(title="Fish oil - Concatenated table", 
       subtitle=paste0("Test AUC = ", round(cvauc, 5)),
       caption="Figure B")
ggsave("plots/special/sPLS_testset_fishoil.png", width = 7, height =5)
