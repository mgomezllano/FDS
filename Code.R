
pacman::p_load(tidyverse, metafor, ggplot2, plyr, latex2exp, dplyr, orchaRd, 
               brms, rotl, ape, phytools, readxl, MuMIn, esc, kittyR, gridExtra,
               ggtree, tidybayes)

setwd("~/Library/CloudStorage/Dropbox/FreqDep_Selection/Code")

# Get the data

pretty <- theme(axis.text = element_text(colour = 'black'),
                axis.line = element_line(colour = "black"),
                panel.background = element_blank(),
                legend.title = element_blank(),
                text = element_text(size = 13))

Da <- read.csv('Data_FDS.csv')
Da$biol_interaction <- factor(Da$biol_interaction)
Da$Sex <- factor(Da$Sex)
Da$Method <- factor(Da$Method)
Da$Fitness_class <- factor(Da$Fitness_class)
Da$Density <- factor(Da$Density)
length(unique(Da$Paper_ID))
nrow(Da)

Da$obs <- seq(1, nrow(Da))

# Create variable with effect sizes (Gs) that do not overlap 0 (Sig) and those who do not (NS)

minGs <- Da$Gs-Da$Var_Gs
maxGs <- Da$Gs+Da$Var_Gs

posneg <- maxGs*minGs

magn <- rep(NA, nrow(Da))

for(i in 1:nrow(Da)){
  if(posneg[i] < 0){
    magn[i] <- "NS" 
  } else {
    magn[i] <- "Sig"
  }
}

Da$magn <- magn

neg <- Da[which(Da$Gs<0),]
pos <- Da[which(Da$Gs>0),]

# ddply(neg, .(magn), summarise, n = length(Gs))
# ddply(pos, .(magn), summarise, n = length(Gs))
# 
# 242/nrow(Da) # perc negative
# 181/nrow(Da) # perc positive
# (nrow(Da)-(242+181))/nrow(Da) # Perc non-sig

####################################################################
########################  DATA BIASES ##############################
####################################################################

# Check the data for biases.

# This model shows there is no temporal bias (older papers showing stronger effects)
res.Year <- rma.mv(Gs, V = Var_Gs,
                   mods = ~ Pub_Year, 
                   random = list(~1|Paper_ID,
                                 ~1|obs,
                                 ~1|species), 
                   data = Da)

summary(res.Year)

# Effect sizes through time
ftime <- ggplot(data = Da, aes(x = Pub_Year, y = Gs)) + 
  geom_point(aes(size = 1/sqrt(Var_Gs)), col = "grey40", alpha = 0.4) +
  geom_smooth(method = 'lm', col = 'black') +
  xlab("Publication year") + ylab("Hedges' g") + theme_bw() +
  labs(size = TeX("Precision $\\left(\\frac{1}{\\SE}\\right)$")) +
  theme(legend.position = c(.6, 0.9), 
        legend.direction = 'horizontal',
        text = element_text(size = 12)) +
  annotate('text', label = 'B', y = 15, x = 1967, size = 7)


res.SE <- rma.mv(Gs, V = Var_Gs,
                   mods = ~ sqrt(Var_Gs), 
                   random = list(~1|Paper_ID,
                                 ~1|obs,
                                 ~1|species), 
                   data = Da)

summary(res.SE)

fsmall <- ggplot(data = Da, aes(x = sqrt(Var_Gs), y = Gs)) +
  geom_point(aes(size = 1/sqrt(Var_Gs)), col = "grey40", alpha = 0.4) +
  geom_smooth(method = 'lm', col = 'black') +
  xlab("Standard Error") + ylab("Hedges' g") + theme_bw() +
  labs(size = TeX("Precision $\\left(\\frac{1}{\\SE}\\right)$")) +
  theme(legend.position = c(.6, 0.9),
        legend.direction = 'horizontal',
        text = element_text(size = 11)) 

pdf("fS2.pdf", height = 6, width = 11)
grid.arrange(ftime, fsmall, nrow = 1)
dev.off()

####################################################################
######################  OVERALL EFFECTS ############################
####################################################################

# First I ran the general model with only the intercept. Selection is 
# mostly negative but not statistically different from 0. The model including 
# the phylogeny fits the data a bit better. 

################### OVERALL

res <- rma.mv(Gs, V = Var_Gs, 
              random = list(~1|Paper_ID, 
                            ~1|obs,
                            ~1|species),
              test = 't', 
              data = Da)

summary(res)

################### PHYLOGENETIC

# # prune tree
Da <- Da %>%
  mutate(species2 = species)

tree <- read.tree("tree")

tree$tip.label = gsub("_", " ", tree$tip.label)

Da2 = Da[Da$species %in% tree$tip.label,] # get only those with tip
tip = tree$tip.label[!tree$tip.label %in% Da2$species] # drop tips without data
tree = drop.tip(tree, tip)
phylogeny <- vcv(tree, corr = TRUE)
nrow(Da2)

resPhylo <- rma.mv(Gs, V = Var_Gs,
                   random = list(~1|Paper_ID,
                                 ~1|obs,
                                 ~1|species, # no-phylogeny
                                 ~1|species2), # phylogeny
                   test = 't', method = 'REML',
                   R = list(species2 = phylogeny), data = Da2)

summary(resPhylo)

AICc(res, resPhylo) # Phylogenetic model has lower AIC

i2_ml(resPhylo)

# Absolute mean effects

formula <- Gs|se(sqrt(Var_Gs)) ~ 1 + 
  (1|Paper_ID) + (1|obs) + (1|species) + (1|gr(species2, cov = phylogeny))
intercept <- fit_model(formula)

# Folded normal distribution
sa <- tidybayes::tidy_draws(intercept)
mu <- sa$b_Intercept
sed <- sd(sa$b_Intercept)
postfnorm <- stats::dnorm(mu, 0, sed)*2*(sed^2) + mu*(2*stats::pnorm(mu, 0, sed) -1)
mean_mpg <- mean(postfnorm)
se_mpg <- sd(postfnorm) / sqrt(length(postfnorm))

ci <- median_qi(postfnorm, .width = 0.95)
mod_table <- data.frame(mean_mpg, ci[2:3])
colnames(mod_table) <- c("mean", "lower", "upper")

pres <- 1/sqrt(Da2$Var_Gs)
Das <- data.frame(Da2$Gs, pres, Da2$magn)
Das$inter <- 'Intercept'

### FIGURE 1

f1C <- ggplot() +
  ggbeeswarm::geom_quasirandom(data = Das, aes(y = abs(Da2.Gs), inter, shape = Da2.magn, size = pres),
                               col = "grey40", alpha = 0.5, width = 0.5) +
  scale_shape_manual('Significance', values = c(1, 19)) +
  guides(shape = guide_legend(override.aes = list(size = 3)))  +
  geom_point(data = mod_table, aes(x = 1, y = mean), col = 'black', size = 4) +
  geom_linerange(data = mod_table, aes(x = 1, ymin = lower, ymax = upper), linewidth = 1.5) +
  scale_y_continuous(limits = c(0, 3.5)) +
  theme_bw() + coord_flip() +
  theme(axis.text = element_text(size = 10, colour = "black", hjust = 0.5)) +
  xlab('') + ylab("Hedge's g") +
  scale_size_continuous("Precision (1/SE)") +
  geom_hline(yintercept = 0, linetype = 'dashed', col = 'grey40') +
  theme(legend.position = c(0.6, 0.1),
        legend.direction = 'horizontal',
        legend.background = element_blank(),
        legend.spacing.y = unit(-0.2, "cm"),
        legend.spacing.x = unit(0.05, 'cm'),
        axis.text.y = element_text(angle = 90, hjust = 0.5, size = 12)) +
  annotate('text', label = 'B', y = 0.3, x = 1.55, size = 7) +
  annotate('text', label = 'K = 745', y = 2.8, x = 1.35, size = 4)

mod_res2 <- mod_results(resPhylo, mod = "1", at = NULL,  group = "Paper_ID")
mod_res2

mod_table2 <- mod_res2$mod_table
mod_table2[1] <- "Intercept"

pres <- 1/sqrt(Da2$Var_Gs)

Das2 <- data.frame(Da2$Gs, pres, Da2$magn)
Das2$inter <- 'Intercept'

f1B <- ggplot() +
  ggbeeswarm::geom_quasirandom(data = Das2, aes(y = Da2.Gs, inter, shape = Da2.magn, size = pres),
                               col = "grey40", alpha = 0.5, width = 0.5) +
  scale_shape_manual('Significance', values = c(1, 19)) +
  guides(shape = guide_legend(override.aes = list(size = 3)))  +
  geom_point(data = mod_table2, aes(x = name, y = estimate), col = 'black', size = 4) +
  geom_linerange(data = mod_table2, aes(x = name, ymin = lowerCL, ymax = upperCL), linewidth = 1.5) +
  scale_y_continuous(limits = c(-2.5, 2.5)) +
  theme_bw() + coord_flip() +
  theme(axis.text = element_text(size = 10, colour = "black", hjust = 0.5)) +
  xlab('') + ylab("Hedge's g") +
  scale_size_continuous("Precision (1/SE)") +
  geom_hline(yintercept = 0, linetype = 'dashed', col = 'grey40') +
  theme(legend.position = c(0.7, 0.1),
        legend.direction = 'horizontal',
        legend.background = element_blank(),
        legend.spacing.y = unit(-0.2, "cm"),
        legend.spacing.x = unit(0.05, 'cm'),
        axis.text.y = element_text(angle = 90, hjust = 0.5, size = 12)) +
  annotate('text', label = 'C', y = -2.2, x = 1.5, size = 7) +
  annotate('text', label = 'K = 745', y = 1.5, x = 1.3, size = 4)

p1 <- grid.arrange(f1B, f1C, ncol = 1)

pdf("f1a.pdf", width = 5, height = 9)
grid.arrange(f1B, f1C, ncol = 1)
dev.off()

###
Kingdom <- Da2$Kingdom
tip.label <- Da2$species2
d <- data.frame(tip.label, Kingdom)

f1A <- ggtree(tree) %<+% d + 
  geom_tiplab(aes(label=label, color=Kingdom), offset = 0, size = 3) + 
  xlim(NA, 1.5) + theme(legend.position = 'none') +
  scale_colour_manual(values = c('#2c7bb6', '#d73027',
                                 "#fec44f", "#54278f",
                                 '#008837')) +
  annotate('text', label = 'A', y = 64, x = 0, size = 7)
  
pdf("f1.pdf", width = 9, height = 9)
grid.arrange(f1A, p1, ncol = 2)
dev.off()

####################################################################
########################  MODERATORS ###############################
####################################################################

resPhyloM <- rma.mv(Gs, V = Var_Gs,
                   mods = ~ Method - 1,
                   random = list(~1|Paper_ID,
                                 ~1|obs,
                                 ~1|species, # non-phylogenetic
                                 ~1|species2), # phylogenetic
                   test = 't', method = 'REML',
                   R = list(species2 = phylogeny), data = Da2)

summary(resPhyloM)

# Figure methods
mod_meth <- mod_results(resPhyloM, mod = "Method", at = NULL,  group = "Paper_ID")
mod_meth <- mod_meth$mod_table
Das2 <- data.frame(Da2$Gs, Da2$Method, pres, Da2$magn)

fmeth <- ggplot() +
  ggbeeswarm::geom_quasirandom(data = Das2, aes(y = Da2.Gs, x = Da2.Method, col = Da2.Method,
                                                shape = Da2.magn, size = pres),
                               alpha = 0.5, width = 0.5) +
  scale_shape_manual('Significance', values = c(1, 19)) +
  guides(shape = guide_legend(override.aes = list(size = 3)))  +
  geom_point(data = mod_meth, aes(x = name, y = estimate), col = 'black', size = 4) +
  geom_linerange(data = mod_meth, aes(x = name, ymin = lowerCL, ymax = upperCL), linewidth = 1.5) +
  scale_y_continuous(limits = c(-2.5, 2.5)) +
  theme_bw() + coord_flip() +
  theme(axis.text = element_text(size = 10, colour = "black", hjust = 0.5)) +
  xlab('') + ylab("Hedge's g") +
  scale_size_continuous("Precision (1/SE)") +
  geom_hline(yintercept = 0, linetype = 'dashed', col = 'grey40') +
  theme(legend.position = 'none') 

 ################### SEX

resPhyloS <- rma.mv(Gs, V = Var_Gs,
                    mods = ~Sex-1,
                    random = list(~1|Paper_ID,
                                  ~1|obs,
                                  ~1|species, # non-phylogenetic
                                  ~1|species2), # phylogenetic
                    test = 't', method = 'REML',
                    R = list(species2 = phylogeny), data = Da2)

summary(resPhyloS)

# Figure sex
mod_sex <- mod_results(resPhyloS, mod = "Sex", at = NULL,  group = "Paper_ID")
mod_sex <- mod_sex$mod_table
Das2 <- data.frame(Da2$Gs, Da2$Sex, pres, Da2$magn)

fsex <- ggplot() +
  ggbeeswarm::geom_quasirandom(data = Das2, aes(y = Da2.Gs, x = Da2.Sex, col = Da2.Sex,
                                                shape = Da2.magn, size = pres),
                               alpha = 0.5, width = 0.5) +
  scale_shape_manual('Significance', values = c(1, 19)) +
  guides(shape = guide_legend(override.aes = list(size = 3)))  +
  geom_point(data = mod_sex, aes(x = name, y = estimate), col = 'black', size = 4) +
  geom_linerange(data = mod_sex, aes(x = name, ymin = lowerCL, ymax = upperCL), linewidth = 1.5) +
  scale_y_continuous(limits = c(-2.5, 2.5)) +
  theme_bw() + coord_flip() +
  theme(axis.text = element_text(size = 10, colour = "black", hjust = 0.5)) +
  xlab('') + ylab("Hedge's g") +
  scale_size_continuous("Precision (1/SE)") +
  geom_hline(yintercept = 0, linetype = 'dashed', col = 'grey40') +
  theme(legend.position = 'none') 

pdf("f2.pdf", width = 10, height = 7)
grid.arrange(fmeth, fsex, ncol = 2)
dev.off()


################### TRAIT

resPhyloT <- rma.mv(Gs, V = Var_Gs,
                    mods = ~ Trait_class-1,
                    random = list(~1|Paper_ID,
                                  ~1|obs,
                                  ~1|species, # non-phylogenetic
                                  ~1|species2), # phylogenetic
                    test = 't', method = 'REML',
                    R = list(species2 = phylogeny), data = Da2)

summary(resPhyloT)

# Figure trait
mod_trait <- mod_results(resPhyloT, mod = "Trait_class", at = NULL,  group = "Paper_ID")
mod_trait <- mod_trait$mod_table
Das2 <- data.frame(Da2$Gs, Da2$Trait_class, pres, Da2$magn)

ftrait <- ggplot() +
  ggbeeswarm::geom_quasirandom(data = Das2, aes(y = Da2.Gs, x = Da2.Trait_class, col = Da2.Trait_class,
                                                shape = Da2.magn, size = pres),
                               alpha = 0.5, width = 0.5) +
  scale_shape_manual('Significance', values = c(1, 19)) +
  guides(shape = guide_legend(override.aes = list(size = 3)))  +
  geom_point(data = mod_trait, aes(x = name, y = estimate), col = 'black', size = 4) +
  geom_linerange(data = mod_trait, aes(x = name, ymin = lowerCL, ymax = upperCL), linewidth = 1.5) +
  scale_y_continuous(limits = c(-3.5, 3.5)) +
  theme_bw() + coord_flip() +
  theme(axis.text = element_text(size = 10, colour = "black", hjust = 0.5)) +
  xlab('') + ylab("Hedge's g") +
  scale_size_continuous("Precision (1/SE)") +
  geom_hline(yintercept = 0, linetype = 'dashed', col = 'grey40') +
  theme(legend.position = 'none') 

################### FITNESS

resPhyloF <- rma.mv(Gs, V = Var_Gs,
                    mods = ~ Fitness_class-1,
                    random = list(~1|Paper_ID,
                                  ~1|obs,
                                  ~1|species, # non-phylogenetic
                                  ~1|species2), # phylogenetic
                    test = 't', method = 'REML',
                    R = list(species2 = phylogeny), data = Da2)

summary(resPhyloF)

# Figure fitness
mod_fit <- mod_results(resPhyloF, mod = "Fitness_class", at = NULL,  group = "Paper_ID")
mod_fit <- mod_fit$mod_table
Das2 <- data.frame(Da2$Gs, Da2$Fitness_class, pres, Da2$magn)

ffit <- ggplot() +
  ggbeeswarm::geom_quasirandom(data = Das2, aes(y = Da2.Gs, x = Da2.Fitness_class, col = Da2.Fitness_class,
                                                shape = Da2.magn, size = pres),
                               alpha = 0.5, width = 0.5) +
  scale_shape_manual('Significance', values = c(1, 19)) +
  guides(shape = guide_legend(override.aes = list(size = 3)))  +
  geom_point(data = mod_fit, aes(x = name, y = estimate), col = 'black', size = 4) +
  geom_linerange(data = mod_fit, aes(x = name, ymin = lowerCL, ymax = upperCL), linewidth = 1.5) +
  scale_y_continuous(limits = c(-2.5, 2.5)) +
  theme_bw() + coord_flip() +
  theme(axis.text = element_text(size = 10, colour = "black", hjust = 0.5)) +
  xlab('') + ylab("Hedge's g") +
  scale_size_continuous("Precision (1/SE)") +
  geom_hline(yintercept = 0, linetype = 'dashed', col = 'grey40') +
  theme(legend.position = 'none') 

################### BIOTIC INTERACTION

resPhyloI <- rma.mv(Gs, V = Var_Gs,
                    mods = ~ biol_interaction-1,
                    random = list(~1|Paper_ID,
                                  ~1|obs,
                                  ~1|species, # non-phylogenetic
                                  ~1|species2), # phylogenetic
                    test = 't', method = 'REML',
                    R = list(species2 = phylogeny), data = Da2)

summary(resPhyloI)

# Figure biological interaction
mod_int <- mod_results(resPhyloI, mod = "biol_interaction", at = NULL,  group = "Paper_ID")
mod_int <- mod_int$mod_table
Das2 <- data.frame(Da2$Gs, Da2$biol_interaction, pres, Da2$magn)

fint <- ggplot() +
  ggbeeswarm::geom_quasirandom(data = Das2, aes(y = Da2.Gs, x = Da2.biol_interaction, col = Da2.biol_interaction,
                                                shape = Da2.magn, size = pres),
                               alpha = 0.5, width = 0.5) +
  scale_shape_manual('Significance', values = c(1, 19)) +
  guides(shape = guide_legend(override.aes = list(size = 3)))  +
  geom_point(data = mod_int, aes(x = name, y = estimate), col = 'black', size = 4) +
  geom_linerange(data = mod_int, aes(x = name, ymin = lowerCL, ymax = upperCL), linewidth = 1.5) +
  scale_y_continuous(limits = c(-2.5, 2.5)) +
  theme_bw() + coord_flip() +
  theme(axis.text = element_text(size = 10, colour = "black", hjust = 0.5)) +
  xlab('') + ylab("Hedge's g") +
  scale_size_continuous("Precision (1/SE)") +
  geom_hline(yintercept = 0, linetype = 'dashed', col = 'grey40') +
  theme(legend.position = 'none') 

pdf("f3.pdf", height = 8, width = 13)
grid.arrange(ftrait, ffit, fint, nrow = 1, widths = c(0.9, 0.9, 1.1))
dev.off()

#########################################################
##################### SUPPLEMENTARY #####################
#########################################################

################### OUTLIERS

# This test is to see how robust was the model to outliers. 
# Did a cook distance analysis to detect highly influencial effect sizes, 
# remove those and ran the model without them, the results are qualitatively the same, 
# so the model is robust and we don't have issues with a couple of studies driving all 
# the patterns. 

cd.model1.1 <- cooks.distance(resPhylo)

plot(cd.model1.1, xlab = "Effect", ylab = "Cook's Distance")
abline(h = 4/resPhylo$k.all, col = "red", lty = 4)
text(Da2[which(cd.model1.1 > 4/resPhylo$k.all), "obs"],
     x = c(1:resPhylo$k.all)[which(cd.model1.1 > 4/resPhylo$k.all)]+5,
     y = cd.model1.1[as.vector(which(cd.model1.1 > 4/resPhylo$k.all))])

vals <- Da2[as.vector(which(cd.model1.1 > 4/resPhylo$k.all)), "obs"]


# results highly robust to outliers removed. Effect size strength increases
model1.1_refit <- rma.mv(Gs, V = Var_Gs,
                         random = list(~1|Paper_ID,
                                       ~1|obs,
                                       ~1|species, # non-phylogenetic
                                       ~1|species2), # phylogenetic
                         test = "t", data = Da2[!Da2$obs %in% vals,])

summary(model1.1_refit)

################### KINGDOM

# res.King <- rma.mv(Gs, V = Var_Gs,
#                   mods = ~ Kingdom-1,
#                   random = list(~1|Paper_ID, 
#                                 ~1|obs, 
#                                 ~1|species),
#                   data = Da)
# 
# summary(res.King)
# 
res.KingP <- rma.mv(Gs, V = Var_Gs,
                     mods = ~Kingdom-1,
                     random = list(~1|Paper_ID,
                                   ~1|obs,
                                   ~1|species,
                                   ~1|species2),
                     test = 't',
                     R = list(species2 = phylogeny), data = Da2)

summary(res.KingP)

# # Figure class
mod_king <- mod_results(res.KingP, mod = "Kingdom", at = NULL,  group = "Paper_ID")
mod_king <- mod_king$mod_table
Das2 <- data.frame(Da2$Gs, Da2$Kingdom, pres, Da2$magn)

fking <- ggplot() +
  ggbeeswarm::geom_quasirandom(data = Das2, aes(y = Da2.Gs, x = Da2.Kingdom, col = Da2.Kingdom,
                                                shape = Da2.magn, size = pres),
                               alpha = 0.5, width = 0.5) +
  scale_shape_manual('Significance', values = c(1, 19)) +
  guides(shape = guide_legend(override.aes = list(size = 3)))  +
  geom_point(data = mod_king, aes(x = name, y = estimate), col = 'black', size = 4) +
  geom_linerange(data = mod_king, aes(x = name, ymin = lowerCL, ymax = upperCL), linewidth = 1.5) +
  scale_y_continuous(limits = c(-6, 2.5)) +
  theme_bw() + coord_flip() +
  theme(axis.text = element_text(size = 10, colour = "black", hjust = 0.5)) +
  xlab('') + ylab("Hedge's g") +
  scale_size_continuous("Precision (1/SE)") +
  geom_hline(yintercept = 0, linetype = 'dashed', col = 'grey40') +
  theme(legend.position = 'none') 

pdf("fS3.pdf", height = 10, width = 8)
fking
dev.off()

# Absolute effect size
formula <- Gs|se(sqrt(Var_Gs)) ~ Kingdom-1 + 
  (1|Paper_ID) + (1|obs) + (1|species) + (1|gr(species2, cov = phylogeny))
king <- fit_model(formula)

sa <- tidybayes::tidy_draws(king)
res1 <- fold(sa, nlevs = length(unique(Da2$Kingdom))-1, colu = length(unique(Da2$Kingdom)))
res1$groups <- sort(unique(Da2$Kingdom))
res1
################### SEX * METHOD

Da2$Sex <- relevel(Da2$Sex, ref = 'Both')

resPhyloMS <- rma.mv(Gs, V = Var_Gs,
                     mods = ~ Sex * Method-1,
                     random = list(~1|Paper_ID,
                                   ~1|obs,
                                   ~1|species, # non-phylogenetic
                                   ~1|species2), # phylogenetic
                     test = 't', method = 'REML',
                     R = list(species2 = phylogeny), data = Da2)

summary(resPhyloMS)

############# ABSOLUTE EFFECT SIZES

# METHOD
formula <- Gs|se(sqrt(Var_Gs)) ~ Method-1 + 
  (1|Paper_ID) + (1|obs) + (1|species) + (1|gr(species2, cov = phylogeny))
method <- fit_model(formula)

sa <- tidybayes::tidy_draws(method)
res1 <- fold(sa, nlevs = length(unique(Da2$Method))-1, colu = length(unique(Da2$Fitness_class)))
res1$groups <- sort(unique(Da2$Method))
res1

# SEX
formula <- Gs|se(sqrt(Var_Gs)) ~ Sex-1 + 
  (1|Paper_ID) + (1|obs) + (1|species) + (1|gr(species2, cov = phylogeny))
sex <- fit_model(formula)

sa <- tidybayes::tidy_draws(sex)
res2 <- fold(sa, nlevs = length(unique(Da2$Sex))-1, colu = length(unique(Da2$Fitness_class)))
res2$groups <- sort(unique(Da2$Sex))
res2

# TRAIT
formula <- Gs|se(sqrt(Var_Gs)) ~ Trait_class-1 + 
  (1|Paper_ID) + (1|obs) + (1|species) + (1|gr(species2, cov = phylogeny))
trait <- fit_model(formula)

sa <- tidybayes::tidy_draws(trait)
res3 <- fold(sa, nlevs = length(unique(Da2$Trait_class))-1)
res3$groups <- sort(unique(Da2$Trait_class))
res3

# FITNESS
formula <- Gs|se(sqrt(Var_Gs)) ~ Fitness_class-1 + 
  (1|Paper_ID) + (1|obs) + (1|species) + (1|gr(species2, cov = phylogeny))
fitness <- fit_model(formula)

sa <- tidybayes::tidy_draws(fitness)
res4 <- fold(sa, nlevs = length(unique(Da2$Fitness_class))-1, colu = length(unique(Da2$Fitness_class)))
res4$groups <- sort(unique(Da2$Fitness_class))
res4

# INTERACTION
formula <- Gs|se(sqrt(Var_Gs)) ~ biol_interaction-1 + 
  (1|Paper_ID) + (1|obs) + (1|species) + (1|gr(species2, cov = phylogeny))
inter <- fit_model(formula)

sa <- tidybayes::tidy_draws(inter)
res5 <- fold(sa, nlevs = length(unique(Da2$biol_interaction))-1, colu = length(unique(Da2$Fitness_class)))
res5$groups <- sort(unique(Da2$biol_interaction))
res5

################################
## Non-phylo

# Method
resM <- rma.mv(Gs, V = Var_Gs,
                    mods = ~ Method - 1,
                    random = list(~1|Paper_ID,
                                  ~1|obs,
                                  ~1|species), # non-phylogenetic
                    test = 't', method = 'REML',
                    data = Da2)

summary(resM)

# Sex
resS <- rma.mv(Gs, V = Var_Gs,
                    mods = ~Sex-1,
                    random = list(~1|Paper_ID,
                                  ~1|obs,
                                  ~1|species), # non-phylogenetic
                    test = 't', method = 'REML',
                    data = Da2)

summary(resS)

# Trait
resT <- rma.mv(Gs, V = Var_Gs,
                    mods = ~ Trait_class-1,
                    random = list(~1|Paper_ID,
                                  ~1|obs,
                                  ~1|species), # non-phylogenetic
                    test = 't', method = 'REML',
                    data = Da2)

summary(resT)

# Fitness
resF <- rma.mv(Gs, V = Var_Gs,
                    mods = ~ Fitness_class-1,
                    random = list(~1|Paper_ID,
                                  ~1|obs,
                                  ~1|species), # non-phylogenetic
                    test = 't', method = 'REML',
                    data = Da2)

summary(resF)

# Interaction
resI <- rma.mv(Gs, V = Var_Gs,
                    mods = ~ biol_interaction-1,
                    random = list(~1|Paper_ID,
                                  ~1|obs,
                                  ~1|species), # non-phylogenetic
                    test = 't', method = 'REML',
                    data = Da2)

summary(resI)
