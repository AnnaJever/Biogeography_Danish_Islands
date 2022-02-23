### Plant dispersal characteristics shape the relationship of diversity with area and isolation

### Walentowitz et al.

### last changes applied on 17. January 2022


### Data preparation ###

# packages
require(XLConnect)
require(rms)
require(tidyverse)
require(gridExtra)
require(vegan)
require(ggeffects)
require(jtools)

# define functions
star <- function(x){a<-""
                  if(x<0.05){a<-"*"}	
                  if(x<0.01){a<-"**"}
                  if(x<0.001){a<-"***"}	
                  a}

lmp <- function(modelobject){if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
                             f <- summary(modelobject)$fstatistic
                             p <- pf(f[1],f[2],f[3],lower.tail=F)
                             attributes(p) <- NULL
                             return(p)}

# get data
data    <- loadWorkbook("Dispersal_Danish_Islands_data.xlsx")
env     <- readWorksheet(data, sheet = "env", header = TRUE)
spec    <- readWorksheet(data, sheet = "spec", header = TRUE)
trait   <- readWorksheet(data, sheet = "trait", header = TRUE)

# data structure

## table: env
  # nr = island specific ID
  # name = island name
  # mainland = name of nearest mainland
  # distance = distance of island to the mainland
  # area = island area
  # population = number of inhabitants
  # pop_density = population density (inhabitants per ha)
  # inhabited = inhabited (1) and uninhabited (0) islands (factor)

## table: species
  # names = species names
  # other column names: correspond to island id (1 indicates species presence and 0 absence)

## table: trait
  # X = project specific species ID
  # name = species name
  # following columns: dispersal syndromes and seed mass of species
  # last column: spermatophyte  = is spermatophyte (1) or is pteridophyte (0)


# Check data

  # traits
  str(trait) # check data type
  sum(!trait$name==spec$name) # should be zero, is zero
  trait$seed.mass <-as.numeric(trait$seed.mass) # make numeric
  trait$spermatophyte <- as.numeric(trait$spermatophyte)

  # environmental data
  str(env)
  env$inhabited <- as.factor(env$inhabited)
  
  # species
  str(spec)
  nrow(spec)

### Complementary calculations ###

# Percentage small islands
  nrow(env[env$area < 10,])/nrow(env)*100 #39%
  
# species number per island
  env$nspec.new   <-colSums(spec[2:ncol(spec)]>0)

# calculate community seed mass per island
  env$mean.seed.weight<-NA
  for(i in colnames(spec)[2:ncol(spec)]){
                                      env[paste("X",env$nr,sep="")==i,
                                      "mean.seed.weight"] <- mean(trait[spec[,i]>0,
                                      "seed.mass"],na.rm=T)
                                      }

# calculate number of dispersal syndrome per islands
  head(trait)
  trait$autochor    <-as.numeric(trait$autochor)
  trait$meteorochor <-as.numeric(trait$meteorochor)
  trait$nautochor   <-as.numeric(trait$nautochor)
  trait$zoochor     <-as.numeric(trait$zoochor)
  trait[,c("autochor","meteorochor", "nautochor", "zoochor")]<-trait[,c("autochor","meteorochor", "nautochor", "zoochor")]/rowSums(trait[,c("autochor","meteorochor", "nautochor", "zoochor")],na.rm=T)

# sums of species per island and dispersal syndrome
  head(trait)
  env$n.autochor    <-NA
  env$n.meteorochor <-NA
  env$n.nautochor   <-NA
  env$n.zoochor     <-NA

  for(i in colnames(spec)[2:ncol(spec)]){
    env[paste("X",env$nr,sep="")==i,"n.autochor"]   <- sum(trait[spec[,i]>0,"autochor"],na.rm=T)
    env[paste("X",env$nr,sep="")==i,"n.meteorochor"]<- sum(trait[spec[,i]>0,"meteorochor"],na.rm=T)
    env[paste("X",env$nr,sep="")==i,"n.nautochor"]  <- sum(trait[spec[,i]>0,"nautochor"],na.rm=T)
    env[paste("X",env$nr,sep="")==i,"n.zoochor"]    <- sum(trait[spec[,i]>0,"zoochor"],na.rm=T)
  }

# percentage dispersal syndrome of floral communities per island
  env$perc.autochor     <-NA
  env$perc.meteorochor  <-NA
  env$perc.nautochor    <-NA
  env$perc.zoochor      <-NA
  env[,c("perc.autochor","perc.meteorochor", "perc.nautochor", "perc.zoochor")]<-env[,c("n.autochor","n.meteorochor", "n.nautochor", "n.zoochor")]/rowSums(env[,c("n.autochor","n.meteorochor", "n.nautochor", "n.zoochor")])*100


### Log10 transformation of data ###
hist(env$nspec.new)
hist(env$area)
hist(env$distance)
hist(env$pop_density[env$inhabited == "1"])

hist(log10(env$nspec.new))
hist(log10(env$area))
hist(log10(env$distance))
hist(log10(env$pop_density[env$inhabited == "1"]))

env$nspec.new_log   <- log10(env$nspec.new)
env$area_log        <- log10(env$area)
env$distance_log    <- log10(env$distance)
env$pop_density_log <- NA
for(i in 1:nrow(env)){
  if(env$inhabited[i] == 1){
    env$pop_density_log[i] <- log10(env$pop_density[i])   
  }
}





### Analysis ###

### Generalized linear models explaining species numbers on islands ###

# Check multicollinearity model GLM_inhab: Correlation area, isolation, island inhabitation
  cor(env[,c("area_log", "distance_log")])  # R = 0.06
  t.test(env$area_log ~ env$inhabited)      # significant relationship
  t.test(env$distance_log ~ env$inhabited)  # no significant relationship

# Check multicollinearity model GLM_pop_dens: Correlation area, isolation, human density
  cor(env[which(!is.na(env$pop_density_log)),c("area_log", "pop_density_log")])  # R = 0.20

# multiple linear regression analysis explaining species richness per island with area + isolation
# selection of a glm with poisson family and poisson distribution as species number is count data

# GLM_inhab
glm_multi_full <- glm(nspec.new ~ area_log + distance_log + inhabited, data = env,
                 family = poisson(link = "log"))

  summary(glm_multi_full)
  coef(glm_multi_full)
  plot(glm_multi_full)

# backwards variable selection
  backwards <- step(glm_multi_full)
  summary(backwards)

# glm with final variable selection
  glm_multi <- glm(nspec.new ~ area_log + inhabited, data = env,
                      family = poisson(link = "log"))
  summary(glm_multi)
  coef(glm_multi)
  plot(glm_multi)

# check goodness of fit (chi-square test based on residual deviance and degrees of freedom)
  1 - pchisq(summary(glm_multi)$deviance, summary(glm_multi)$df.residual) # (has to be < 0.05)
  with(summary(glm_multi), 1 - deviance/null.deviance)# 0.85 # R²

# check residuals for normality and heteroscedasticity
  hist(glm_multi$residuals) # histogram of residuals
  cor(x = glm_multi$residuals, y = env$area_log, method = "pearson") # cor. of residuals with area
# no strong correlations of residuals with (continuous) explanatory variables

### Complementary glm with subset, including population density
env_sub <- env[env$inhabited == "1",] # create subset with inhabited islands
cor(env_sub[,c("area_log", "pop_density_log")])

# GLM_pop_dens
  glm_multi_com <- glm(nspec.new ~ area_log + pop_density_log , data = env_sub,
                 family = poisson(link = "log"))
  summary(glm_multi_com)
  coef(glm_multi_com)
  plot(glm_multi_com)

# variable selection
  backwards_com <- step(glm_multi_com)
  summary(backwards_com) # stays the same as initial model



  
  
  
### Figure 2 ###
  
# Plot results of multivariate glm
# figure 2a
ggeffect(glm_multi, terms = c("area_log", "inhabited")) %>%
  plot()

# figure 2b
plot_summs(glm_multi, scale = T, ci_level = 0.95, colors = "black")

# figure 2c
ggeffect(glm_multi_com,terms = c("area_log", "pop_density_log")) %>%
  plot()

# figure 2d
plot_summs(glm_multi_com, scale = T, ci_level = 0.95, colors = "black")






### Additional SAR calculation with Arrhenius function ###

# make results of this study comparable to other studies

nspec <- env$nspec.new
area  <- as.numeric(env$area)

# The Arrhenius model
marr <- nls(nspec ~ SSarrhenius(area, k, z))
marr

# plot SAR
plot(nspec ~ area,
     xlab = "Island Area (ha)", ylab = "Number of Species",
     ylim = c(1, max(nspec)),
     col = "black",
     pch = 16)

xtmp <- seq(min(area), max(area), len=151)
lines(xtmp, predict(marr, newdata=data.frame(area = xtmp)), lwd=2, col = "dodgerblue4")





### Figure 3 ###

# Plot spec richness and dispersal syndromes

# build linear regressions
m1          <-lm(log10(n.zoochor) ~ area_log ,data= env)
m2          <-lm(log10(n.meteorochor) ~ area_log, data= env)
m3          <-lm(log10(n.nautochor) ~ area_log, data= env)
m4          <-lm(log10(n.autochor) ~ area_log, data= env)

# plot (figure 3)
plot(log10(n.zoochor) ~ area_log,
     data= env,xlab="Area (ha, log)",
     ylab="Species richness (log)"
     ,col="goldenrod", pch = 19)

points(x = log10(env$area), y = log10(env$n.meteorochor),
       col = "darkslategrey", pch = 19)
points(x = log10(env$area), y = log10(env$n.nautochor),
       col = "darkcyan", pch = 19)
points(x = log10(env$area), y = log10(env$n.autochor),
       col = "coral2", pch = 19)

abline(m1, col = "goldenrod")
abline(m2, col = "darkslategrey")
abline(m3, col = "darkcyan")
abline(m4, col = "coral2")

mtext(text=bquote(Zoochory: italic(R)^2 == .(round(summary(m1)$adj.r, 2)) * .(star(lmp(m1)))),
      line=-1.5,cex=0.8,
      col="goldenrod",
      side = 3, at=-1.5, adj = 0)
mtext(text=bquote(Anemochory: italic(R)^2 == .(round(summary(m2)$adj.r, 2)) * .(star(lmp(m2)))),
      line=-2.5,cex=0.8,
      col="darkslategrey",
      side = 3, at=-1.5, adj = 0)
mtext(text=bquote(Hydrochory: italic(R)^2 == .(round(summary(m3)$adj.r, 2)) * .(star(lmp(m3)))),
      line=-3.5,cex=0.8,
      col="darkcyan",
      side = 3, at=-1.5, adj = 0)
mtext(text=bquote(Autochory: italic(R)^2 == .(round(summary(m4)$adj.r, 2)) * .(star(lmp(m4)))),
      line=-4.5,cex=0.8,
      col="coral2",
      side = 3, at=-1.5, adj = 0)





### Figure 4 ###

# % dispersal syndrome in relationship with area, isolation and human density

# models and plot

# a) area - seed weight (H1)
fit_a           <- glm(round(mean.seed.weight)~log10(area),data=env,family=poisson())
summary(fit_a)
m_a_pseudoR2    <- lrm(round(env$mean.seed.weight)~log10(env$area))

plot_a          <- ggplot(env, aes(log10(area), mean.seed.weight)) +
                    geom_point() +
                    geom_smooth(method = "glm", se = F, 
                                method.args = list(family = "poisson"),
                                colour = "black")+
                    theme(panel.background = element_rect(fill = "NA", colour = "black"))+
                    annotate("text", x = 1, y = 1000,
                             label= paste0("Pseud R² = ", round(m_a_pseudoR2$stats[["R2"]][1],
                                                                digits = 2), " ",
                                           star(summary(fit_a)$coefficients[8])))+
                    labs(x = "Area (ha, log)", y = "Mean seed mass (mg)")

# b) area - zoochory
fit_b           <- glm(round(perc.zoochor)~log10(area),data=env,family=poisson())
summary(fit_b)
m_b_pseudoR2    <- lrm(round(env$perc.zoochor)~log10(env$area))

plot_b          <- ggplot(env, aes(log10(area), perc.zoochor)) +
                    geom_point() +
                    geom_smooth(method = "glm", se = F, 
                                method.args = list(family = "poisson"),
                                colour = "black")+
                    theme(panel.background = element_rect(fill = "NA", colour = "black"))+
                    annotate("text", x = 1, y = 40,
                             label= paste0("Pseud R² = ", round(m_b_pseudoR2$stats[["R2"]][1],
                                                                digits = 2), " ",
                                           star(summary(fit_b)$coefficients[8])))+
                    labs(x = "Area (ha, log)", y = "Zoochory (%)")

# c) area - hydrochory
fit_c           <- glm(round(perc.nautochor)~log10(area),data=env,family=poisson())
summary(fit_c)
m_c_pseudoR2    <- lrm(round(env$perc.nautochor)~log10(env$area))

plot_c          <- ggplot(env, aes(log10(area), perc.nautochor)) +
                    geom_point() +
                    geom_smooth(method = "glm", se = F, 
                                method.args = list(family = "poisson"),
                                colour = "black")+
                    theme(panel.background = element_rect(fill = "NA", colour = "black"))+
                    annotate("text", x = 1, y = 52,
                             label= paste0("Pseud R² = ", round(m_c_pseudoR2$stats[["R2"]][1],
                                                                digits = 2), " ",
                                           star(summary(fit_c)$coefficients[8])))+
                    labs(x = "Area (ha, log)", y = "Hydrochory (%)")

# d) area - anemochory
fit_d           <- glm(round(perc.meteorochor)~log10(area),data=env,family=poisson())
summary(fit_d)

m_d_pseudoR2    <- lrm(round(env$perc.meteorochor)~log10(env$area))

plot_d          <- ggplot(env, aes(log10(area), perc.meteorochor)) +
  geom_point() +
  geom_smooth(method = "glm", se = F, 
              method.args = list(family = "poisson"),
              colour = "black")+
  theme(panel.background = element_rect(fill = "NA", colour = "black"))+
  annotate("text", x = 1, y = 15,
           label= paste0("Pseud R² = ", round(m_d_pseudoR2$stats[["R2"]][1],
                                              digits = 2), " ",
                         star(summary(fit_d)$coefficients[8])))+
  labs(x = "Area (ha, log)", y = "Anemochory (%)")

# e) area - autochory
fit_e           <- glm(round(perc.autochor)~log10(area),data=env,family=poisson())
summary(fit_e)
m_e_pseudoR2    <- lrm(round(env$perc.autochor)~log10(env$area))

plot_e          <- ggplot(env, aes(log10(area), perc.autochor)) +
                    geom_point() +
                    theme(panel.background = element_rect(fill = "NA", colour = "black"))+
                    labs(x = "Area (ha, log)", y = "Autochory (%)")

# f) isolation - seed weight (H1)
fit_f           <- glm(round(mean.seed.weight)~log10(distance),data=env,family=poisson())
summary(fit_f)
m_f_pseudoR2    <- lrm(round(env$mean.seed.weight)~log10(env$distance))

plot_f          <- ggplot(env, aes(log10(distance), mean.seed.weight)) +
                    geom_point() +
                    geom_smooth(method = "glm", se = F,
                                method.args = list(family = "poisson"),
                                colour = "black")+
                    theme(panel.background = element_rect(fill = "NA", colour = "black"))+
                    annotate("text", x = 2.3, y = 1000,
                             label= paste0("Pseud R² = ", round(m_f_pseudoR2$stats[["R2"]][1],
                                                                digits = 3), " ",
                                           star(summary(fit_f)$coefficients[8])))+
                    labs(x = "Distance (m, log)", y = "Mean seed mass (mg)")

# g) isolation - zoochory
fit_g           <- glm(round(perc.zoochor)~log10(distance),data=env,family=poisson())
summary(fit_g)

plot_g          <- ggplot(env, aes(log10(distance), perc.zoochor)) +
                    geom_point()+
                    theme(panel.background = element_rect(fill = "NA", colour = "black"))+
                    labs(x = "Distance (m, log))", y = "Zoochory (%)")

# h) isolation - hydrochory
fit_h           <- glm(round(perc.nautochor)~log10(distance),data=env,family=poisson())
summary(fit_h)

plot_h          <- ggplot(env, aes(log10(distance), perc.nautochor)) +
                    geom_point() +
                    theme(panel.background = element_rect(fill = "NA", colour = "black"))+
                    labs(x = "Distance (m, log)", y = "Hydrochory (%)")

# i) isolation - anemochory
fit_i           <- glm(round(perc.meteorochor)~log10(distance),data=env,family=poisson())
summary(fit_i)

plot_i          <- ggplot(env, aes(log10(distance), perc.meteorochor)) +
                    geom_point() +
                    theme(panel.background = element_rect(fill = "NA", colour = "black"))+
                    labs(x = "Distance (m, log)", y = "Anemochory (%)")

# j) isolation - autochory
fit_j           <- glm(round(perc.autochor)~log10(distance),data=env,family=poisson())
summary(fit_j)

plot_j          <- ggplot(env, aes(log10(distance), perc.autochor)) +
                    geom_point() +
                    theme(panel.background = element_rect(fill = "NA", colour = "black"))+
                    labs(x = "Distance (m, log)", y = "Autochory (%)")

# k) population density - seed weight

env_sub <- env[which(!is.na(env$pop_density_log)),]

fit_k           <- glm(round(mean.seed.weight)~log10(pop_density),data=env_sub,family=poisson())
summary(fit_k)

m_k_pseudoR2    <- lrm(round(env_sub$mean.seed.weight)~log10(env_sub$pop_density))

plot_k          <- ggplot(env_sub, aes(log10(pop_density), mean.seed.weight)) +
                    geom_point() +
                    geom_smooth(method = "glm", se = F,
                                method.args = list(family = "poisson"),
                                colour = "black")+
                    theme(panel.background = element_rect(fill = "NA", colour = "black"))+
                    annotate("text", x = -1, y = 1000,
                             label= paste0("Pseud R² = ", round(m_k_pseudoR2$stats[["R2"]][1],
                                                                digits = 3), " ",
                                           star(summary(fit_k)$coefficients[8])))+
                    labs(x = "Population density", y = "Mean seed mass (mg)")

# l) population density - zoochory

fit_l           <- glm(round(perc.zoochor)~log10(pop_density),data=env_sub,family=poisson())
summary(fit_l)

plot_l          <- ggplot(env_sub, aes(log10(pop_density), perc.zoochor)) +
                    geom_point() +
                    theme(panel.background = element_rect(fill = "NA", colour = "black"))+
                    labs(x = "Population density", y = "Zoochory (%)")

# m) population density - hydrochory

fit_m           <- glm(round(perc.nautochor)~log10(pop_density),data=env_sub,family=poisson())
summary(fit_m)

plot_m          <- ggplot(env_sub, aes(log10(pop_density), perc.nautochor)) +
                    geom_point() +
                    theme(panel.background = element_rect(fill = "NA", colour = "black"))+
                    labs(x = "Population density", y = "Hydrochory (%)")

# n) population density - anemochory

fit_n           <- glm(round(perc.meteorochor)~log10(pop_density),data=env_sub,family=poisson())
summary(fit_n)

plot_n          <- ggplot(env_sub, aes(log10(pop_density), perc.meteorochor)) +
                    geom_point() +
                    theme(panel.background = element_rect(fill = "NA", colour = "black"))+
                    labs(x = "Population density", y = "Anemochory (%)")

# o) population density - autochory

fit_o           <- glm(round(perc.autochor)~log10(pop_density),data=env_sub,family=poisson())
summary(fit_o)

plot_o          <- ggplot(env_sub, aes(log10(pop_density), perc.autochor)) +
                    geom_point() +
                    theme(panel.background = element_rect(fill = "NA", colour = "black"))+
                    labs(x = "Population density", y = "Autochory (%)")

# Multipanel plot
x11()
grid.arrange(plot_a, plot_b, plot_c, plot_d, plot_e,
             plot_f, plot_g, plot_h, plot_i, plot_j,
             plot_k, plot_l, plot_m, plot_n, plot_o,
             ncol=5)
dev.off()





### Supplementary analysis ###

### Repeat main analysis including distance only with spermatophytes, excluding pteridophytes

# Get data and reduce to spermatophytes
  env2      <- env
  trait2    <- trait[trait$spermatophyte == "1",]
  spec2     <- spec[which(spec$name %in% trait2$name),]

### Complementary calculations ###

# species number per island
  env2$nspec.new   <-colSums(spec2[2:ncol(spec2)]>0)

# calculate community seed mass per island
  env2$mean.seed.weight<-NA
  for(i in colnames(spec2)[2:ncol(spec2)]){
    env2[paste("X",env2$nr,sep="")==i,
        "mean.seed.weight"] <- mean(trait2[spec2[,i]>0,
                                          "seed.mass"],na.rm=T)
  }

# calculate number of dispersal syndrome per islands
  trait2$autochor    <-as.numeric(trait2$autochor)
  trait2$meteorochor <-as.numeric(trait2$meteorochor)
  trait2$nautochor   <-as.numeric(trait2$nautochor)
  trait2$zoochor     <-as.numeric(trait2$zoochor)
  trait2[,c("autochor","meteorochor", "nautochor", "zoochor")]<-trait2[,c("autochor","meteorochor", "nautochor", "zoochor")]/rowSums(trait2[,c("autochor","meteorochor", "nautochor", "zoochor")],na.rm=T)

# sums of species per island and dispersal syndrome
  for(i in colnames(spec2)[2:ncol(spec2)]){
    env2[paste("X",env2$nr,sep="")==i,"n.autochor"]   <- sum(trait2[spec2[,i]>0,"autochor"],na.rm=T)
    env2[paste("X",env2$nr,sep="")==i,"n.meteorochor"]<- sum(trait2[spec2[,i]>0,"meteorochor"],na.rm=T)
    env2[paste("X",env2$nr,sep="")==i,"n.nautochor"]  <- sum(trait2[spec2[,i]>0,"nautochor"],na.rm=T)
    env2[paste("X",env2$nr,sep="")==i,"n.zoochor"]    <- sum(trait2[spec2[,i]>0,"zoochor"],na.rm=T)
  }

# percentage dispersal syndrome of floral communities per island
  env2$perc.autochor     <-NA
  env2$perc.meteorochor  <-NA
  env2$perc.nautochor    <-NA
  env2$perc.zoochor      <-NA
  env2[,c("perc.autochor","perc.meteorochor", "perc.nautochor", "perc.zoochor")]<-env2[,c("n.autochor","n.meteorochor", "n.nautochor", "n.zoochor")]/rowSums(env2[,c("n.autochor","n.meteorochor", "n.nautochor", "n.zoochor")])*100

### Log10 transformation of data
  env2$nspec.new_log   <- log10(env2$nspec.new)
  env2$area_log        <- log10(env2$area)
  env2$distance_log    <- log10(env2$distance)


  
  
### GLMs ###
# multiple linear regression analysis explaining species richness per island with area + isolation + inhabited/uninhabited
# selection of a glm with poisson family and poisson distribution as species number is count data

glm_multi_full2 <- glm(nspec.new ~ area_log + inhabited, data = env2,
                 family = poisson(link = "log"))

  summary(glm_multi_full2)
  coef(glm_multi_full2)
  plot(glm_multi_full2)

  backwards2 <- step(glm_multi_full2)
  summary(backwards2)

glm_multi2 <- glm(nspec.new ~ area_log + inhabited, data = env2,
                       family = poisson(link = "log"))

summary(glm_multi2)
plot(glm_multi_full2)

glm_com2 <- glm(nspec.new ~ area_log + pop_density_log, data = env2,
                  family = poisson(link = "log"))

summary(glm_com2)
plot(glm_com2)

### % dispersal syndrome in relationship with isolation ###

# f) isolation - seed weight (H1)
fit_f2 <- glm(round(mean.seed.weight)~log10(distance),data=env2,family=poisson())
summary(fit_f2)

# g) isolation - zoochory
fit_g2 <- glm(round(perc.zoochor)~log10(distance),data=env2,family=poisson())
summary(fit_g2)

# h) isolation - hydrochory
fit_h2 <- glm(round(perc.nautochor)~log10(distance),data=env2,family=poisson())
summary(fit_h2)

# i) isolation - anemochory
fit_i2 <- glm(round(perc.meteorochor)~log10(distance),data=env2,family=poisson())
summary(fit_i2)

# j) isolation - autochory
fit_j2 <- glm(round(perc.autochor)~log10(distance),data=env2,family=poisson())
summary(fit_j2)

### end ###