
pacman::p_load(tidyverse, metafor, ggplot2, plyr, latex2exp, dplyr, orchaRd, 
               brms, rotl, ape, phytools, readxl, MuMIn, esc, kittyR, gridExtra,
               ggtree)

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
                   data = Da2)

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

### FIGURE 1

mod_res <- mod_results(res, mod = "1", at = NULL,  group = "Paper_ID")
mod_res

f1B <- orchard_plot(mod_res, mod="1", xlab = "Hedge's g", condition.lab = "Intercept",
                   k.pos = 2, g = F,
                   trunk.size = 10, branch.size = 1.5, twig.size = 0.3) +
  scale_y_continuous(limits = c(-2.5, 2.5)) +
  scale_fill_manual(values="grey40") +
  scale_colour_manual(values="grey40") +
  xlab('Intercept') +
  annotate('text', label = 'B', y = -2.5, x = 1.5, size = 7)

mod_res2 <- mod_results(resPhylo, mod = "1", at = NULL,  group = "Paper_ID")
mod_res2

f1C <- orchard_plot(mod_res2, mod="1", xlab = "Hedge's g",
                   k.pos = 2, g = F,
                   trunk.size = 10, branch.size = 1.5, twig.size = 0.3) +
  scale_y_continuous(limits = c(-2.5, 2.5)) +
  scale_fill_manual(values="grey40") +
  scale_colour_manual(values="grey40") + xlab('Intercept') +
  annotate('text', label = 'C', y = -2.5, x = 1.5, size = 7)

p1 <- grid.arrange(f1B, f1C, ncol = 1)

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

################### METHODOLOGY

# res.Method <- rma.mv(Gs, V = Var_Gs,
#                      mods = ~ Method - 1,
#                      random = list(~1|Paper_ID, 
#                                    ~1|obs, 
#                                    ~1|species),
#                      data = Da)
# summary(res.Method)

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
fmeth <- orchard_plot(resPhyloM, mod = "Method", group = "Method", 
             xlab = "Hedges' g", g = F, angle = 0,
             trunk.size = 10, branch.size = 1.5, twig.size = 0.3,  
             k.pos = c(2, 0.5, 1)) +
  scale_y_continuous(limits = c(-2.5, 2.5)) +
  annotate('text', label = 'A', y = -2.3, x = 4.3, size = 7)


################### MINIMUM FREQUENCY

# Minimum frequency
# res.Min <- rma.mv(Gs, V = Var_Gs,
#                   mods = ~ min_freq,
#                   random = list(~1|Paper_ID, 
#                                 ~1|obs, 
#                                 ~1|species),
#                   data = Da)
# 
# summary(res.Min)
# 
resPhyloMin <- rma.mv(Gs, V = Var_Gs,
                    mods = ~ min_freq,
                    random = list(~1|Paper_ID,
                                  ~1|obs,
                                  ~1|species, # non-phylogenetic
                                  ~1|species2), # phylogenetic
                    test = 't', method = 'REML',
                    R = list(species2 = phylogeny), data = Da2)

summary(resPhyloMin)

fmin <- ggplot(data = Da2, aes(x = min_freq, y = Gs)) +
  geom_point(aes(size = 1/sqrt(Var_Gs)), col = "#54278f", alpha = 0.4) +
  geom_abline(intercept = resPhyloMin$b[1], slope = resPhyloMin$b[2]) +
  geom_abline(intercept = resPhyloMin$b[1] + resPhyloMin$se[1],
              slope = resPhyloMin$b[2] + resPhyloMin$se[2], lty = 2) +
  geom_abline(intercept = resPhyloMin$b[1] - resPhyloMin$se[1],
              slope = resPhyloMin$b[2] - resPhyloMin$se[2], lty = 2) +
  xlab("Minimum frequency") + ylab("Hedges' g") + theme_bw() +
  labs(size = TeX("Precision $\\left(\\frac{1}{\\SE}\\right)$")) +
  theme(legend.position = c(.6, 0.9),
        legend.direction = 'horizontal',
        text = element_text(size = 11)) +
  annotate('text', label = 'B', y = 14.5, x = 0.01, size = 7)


pdf("f2.pdf", width = 10, height = 7)
grid.arrange(fmeth, fmin, ncol = 2)
dev.off()


 ################### SEX

# res.Sex <- rma.mv(Gs, V = Var_Gs,
#                   mods = ~Sex-1,
#                   random = list(~1|Paper_ID, 
#                                 ~1|obs, 
#                                 ~1|species),
#                   data = Da)
# 
# summary(res.Sex)

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
fsex <- orchard_plot(resPhyloS, mod = "Sex", group = "Sex", 
             xlab = "Hedges' g", g = F, twig.size = 0.3,  angle = 0,
             trunk.size = 10, branch.size = 1.5, k.pos = c(2, 0.5, 1)) +
  scale_y_continuous(limits = c(-2.5, 2.5)) +
  annotate('text', label = 'A', y = -2.3, x = 5.3, size = 7)


################### TRAIT

# res.Trait <- rma.mv(Gs, V = Var_Gs,
#                     mods = ~ Trait_class - 1,
#                     random = list(~1|Paper_ID, 
#                                   ~1|obs, 
#                                   ~1|species),
#                     data = Da)
# 
# summary(res.Trait)

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
ftrait <- orchard_plot(resPhyloT, mod = "Trait_class", group = "Trait_class", 
             xlab = "Hedges' g", g = F, angle = 0,
             trunk.size = 10, branch.size = 1.5, k.pos = c(2, 0.5, 1),
             twig.size = 0.3) +
  scale_y_continuous(limits = c(-3.5, 2.5)) +
  annotate('text', label = 'B', y = -3.3, x = 7.25, size = 7)

pdf("f3.pdf", height = 8, width = 10)
grid.arrange(fsex, ftrait, nrow = 1)
dev.off()


################### FITNESS

# res.Fitness <- rma.mv(Gs, V = Var_Gs,
#                       mods = ~Fitness_class - 1,
#                       random = list(~1|Paper_ID, 
#                                     ~1|obs, 
#                                     ~1|species),
#                       data = Da)
# 
# summary(res.Fitness)

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
ffit <- orchard_plot(resPhyloF, mod = "Fitness_class", group = "Fitness_class", 
             xlab = "Hedges' g", g = F, twig.size = 0.3,  angle = 0,
             trunk.size = 10, branch.size = 1.5, k.pos = c(2, 0.5, 1)) +
  scale_y_continuous(limits = c(-2.5, 2.5)) +
  annotate('text', label = 'A', y = -2.35, x = 3.45, size = 7)

################### BIOTIC INTERACTION

# res.Interaction <- rma.mv(Gs, V = Var_Gs,
#                           mods = ~ biol_interaction - 1,
#                           random = list(~1|Paper_ID, 
#                                         ~1|obs, 
#                                         ~1|species),
#                           data = Da)
# 
# summary(res.Interaction)

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
fint <- orchard_plot(resPhyloI, mod = "biol_interaction", group = "biol_interaction", 
             xlab = "Hedges' g", g = F, twig.size = 0.3,  angle = 0,
             trunk.size = 10, branch.size = 1.5, k.pos = c(2, 0.5, 1)) +
  scale_y_continuous(limits = c(-2.5, 2.5)) +
  annotate('text', label = 'B', y = -2.3, x = 6.3, size = 7)

pdf("f4.pdf", height = 8, width = 10)
grid.arrange(ffit, fint, nrow = 1, widths = c(0.9, 1.1))
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
fking <- orchard_plot(res.KingP, mod = "Kingdom", group = "Kingdom",
                       xlab = "Hedges' g", g = F, twig.size = 0.3,  angle = 0,
                       trunk.size = 10, branch.size = 1.5, k.pos = c(-5, 0.5, 1)) +
  scale_y_continuous(limits = c(-6, 2.5)) 


pdf("fS3.pdf", height = 10, width = 8)
fking
dev.off()

################### FREQUENCY * METHOD

resPhyloMM <- rma.mv(Gs, V = Var_Gs,
                     mods = ~ min_freq * Method-1,
                     random = list(~1|Paper_ID,
                                   ~1|obs,
                                   ~1|species, # non-phylogenetic
                                   ~1|species2), # phylogenetic
                     test = 't', method = 'REML',
                     R = list(species2 = phylogeny), data = Da2)

summary(resPhyloMM)

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






